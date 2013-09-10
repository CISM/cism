!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_front_end.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module cism_front_end_module

contains

subroutine cism_front_end()
  !*FD The CISM front-end is used to connect both the standalone driver
  !*FD (cism_driver) or the CISM interface to CESM (cism_cesm_interface),
  !*FD to the internal and external dycore interface programs.  These are
  !* cism_internal_dycore_interface and cism_external_dycore_interface.

  use parallel

  use glimmer_global
  use glide
  use glissade
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glide_nc_custom, only: glide_nc_fillall
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glide_io, only: glide_io_writeall

  use cism_internal_dycore_interface_module
  use cism_external_dycore_interface_module
  

!!  use glimmer_horiz_bcs, only : horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
!!                                horiz_bcs_unstag_scalar, horiz_bcs_stag_scalar

  use glide_stop, only: glide_finalise
  use glide_diagnostics

  implicit none


  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=dp) :: time                   ! model time in years
  real(kind=dp) :: t1,t2
  integer :: clock,clock_rate,ret

  integer :: tstep_count
    print *,'Entering CISM Front End'

  !TODO - call this only for parallel runs?
!  call parallel_initialise     

!  call glimmer_GetCommandline()

  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))
  
  ! setup paths
  call filenames_init(commandline_configname)

  ! read configuration
  call ConfigRead(commandline_configname,config)

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise profiling
  call profile_init(model%profile,'glide.profile')

  call t_startf('simple glide')

  ! initialise GLIDE
    call t_startf('glide initialization')

  call glide_config(model,config)

  ! This call is needed only if running the EISMINT test cases
  call simple_initialise(climate,config)

  if (model%options%whichdycore == DYCORE_GLIDE) then
     call glide_initialise(model)
  else       ! glam/glissade dycore	
     call glissade_initialise(model)
  endif

  call CheckSections(config)
 
 ! fill dimension variables on output files
  call glide_nc_fillall(model)

  time = model%numerics%tstart
  tstep_count = 0
  model%numerics%time = time    ! MJH added 1/10/13 - the initial diagnostic glissade solve won't know 
                                !                     the correct time on a restart unless we set it here.

     ! These calls are needed only for the EISMINT test cases, and they are not needed
     ! for initialization provided they are called at the start of each timestep.
!### ! call simple_massbalance(climate,model,time)
!### ! call simple_surftemp(climate,model,time)

  call spinup_lithot(model)
  call t_stopf('glide initialization')
 
  ! run an internal or external dycore, depending on setting external_dycore_type

print *,"external_dycore_type: ",model%options%external_dycore_type

  if (model%options%external_dycore_type .EQ. 0) then      ! NO_EXTERNAL_DYCORE) then
    call cism_internal_dycore_interface(model%options%whichdycore,model,climate,time,tstep_count)
    print *,'Exited Internal Dycore'
  else
    print *,'Using External Dycore' 
    call cism_external_dycore_interface(model%options%external_dycore_type,model)
  endif

  call t_stopf('simple glide')

  ! finalise GLIDE
  call glide_finalise(model)

  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)

#if (! defined CCSMCOUPLED && ! defined CESMTIMERS)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
#endif

  call close_log

  !TODO - call this only for parallel runs?
   call parallel_finalise
end subroutine cism_front_end

end module cism_front_end_module
