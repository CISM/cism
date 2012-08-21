!CLEANUP - simple_glide.F90
! Added calls to glissade driver

!TODO - Change name to simple_driver?

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_glide.f90 - part of the Glimmer-CISM ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!TODO - Want a program that can call either glide or glissade init, run, and finalize subroutines.
!       Call it simple_driver?

program simple_glide
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases

  use parallel

  use glimmer_global, only:rk
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glimmer_horiz_bcs, only : horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
                                horiz_bcs_unstag_scalar

  use glissade

!TODO - Change to 'glimmer_diagnostics'
  use glide_diagnostics

  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret
  character*3 suffix

  integer :: tstep_count

!TODO - parallel (glissade) only
  call parallel_initialise     

  call glimmer_GetCommandline()

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

!TODO - Initialize either glide or glissade.
  ! initialise GLIDE
    call t_startf('glide initialization')

!TODO - glide_config -> glimmer_config?

  call glide_config(model,config)

!TODO - Is this call always needed?
  call simple_initialise(climate,config)

  if (model%options%whichdycore == DYCORE_GLIDE) then
     call glide_initialise(model)
  else
     call glissade_initialise(model)
  endif

  call CheckSections(config)

!TODO - Change to glimmer_nc_fillall?
  ! fill dimension variables
  call glide_nc_fillall(model)

  time = model%numerics%tstart

!TODO - Are these calls always needed?
  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)

  call spinup_lithot(model)
    call t_stopf('glide initialization')

  tstep_count = 0

  suffix = '_t1'

  do while(time <= model%numerics%tend)

     if (tstep_count > 0) suffix = '   '

!TODO - Change to glimmer_tstep?
     call t_startf('glide_tstep'//suffix)

!TODO - Add subroutine glide_step that calls glide_tstep_p1/p2/p3?

!TODO - What does t_adj_detailf do?

     if (model%options%whichdycore == DYCORE_GLIDE) then

       call t_startf('glide_tstep_p1'//suffix)
       if (tstep_count == 0) call t_adj_detailf(+10)
         call glide_tstep_p1(model,time)
       if (tstep_count == 0) call t_adj_detailf(-10)
       call t_stopf('glide_tstep_p1'//suffix)

       call t_startf('glide_tstep_p2'//suffix)
       if (tstep_count == 0) call t_adj_detailf(+10)
         call glide_tstep_p2(model)
       if (tstep_count == 0) call t_adj_detailf(-10)
       call t_stopf('glide_tstep_p2'//suffix)

       call t_startf('glide_tstep_p3'//suffix)
       if (tstep_count == 0) call t_adj_detailf(+10)
         call glide_tstep_p3(model)
       if (tstep_count == 0) call t_adj_detailf(-10)
       call t_stopf('glide_tstep_p3'//suffix)

     else   ! glissade dycore

!TODO - Make this a single subroutine call

       call t_startf('glissade_tstep'//suffix)
       if (tstep_count == 0) call t_adj_detailf(+10)
         call glissade_tstep(model,time)
       if (tstep_count == 0) call t_adj_detailf(-10)
       call t_stopf('glissade_tstep'//suffix)

     endif   ! glide or glissade dycore

     ! override masking stuff for now  !TODO - What does this mean?

     tstep_count = tstep_count + 1

!TODO - Change to glimmer_write_diag?
     ! write ice sheet diagnostics
     if (mod(tstep_count, model%numerics%ndiag)==0)  then
        call glide_write_diag(model, time, model%numerics%idiag, &
                                           model%numerics%jdiag)
     endif

!TODO - Remove this comment?
     ! Redistribute calls here to spread the data back out.


     if (model%options%whichdycore == DYCORE_GLIDE) then

     ! Perform parallel operations for restart files
       call t_startf('glide_tstep_postp3'//suffix)
       if (tstep_count == 1) call t_adj_detailf(+10)
         call glide_tstep_postp3(model)
       if (tstep_count == 1) call t_adj_detailf(-10)
       call t_stopf('glide_tstep_postp3'//suffix)

     else   ! glissade dycore

!TODO - I think we can safely remove all these parallel halo calls.
!       Before doing so, make sure we have the required calls in glissade.F90.

       call t_startf('simple_glide_halo_upd'//suffix)
       call parallel_halo(model%stress%efvs)
       call parallel_halo(model%velocity%uvel)
       call parallel_halo(model%velocity%vvel)
       call parallel_halo(model%velocity%uflx)
       call parallel_halo(model%velocity%vflx)
       call parallel_halo(model%velocity%velnorm)
       call parallel_halo(model%geometry%thck)
       call parallel_halo(model%geomderv%stagthck)
       call parallel_halo(model%climate%acab)
       call parallel_halo(model%geomderv%dusrfdew)
       call parallel_halo(model%geomderv%dusrfdns)
       call parallel_halo(model%geomderv%dthckdew)
       call parallel_halo(model%geomderv%dthckdns)
       call parallel_halo(model%stress%tau%xx)
       call parallel_halo(model%stress%tau%yy)
       call parallel_halo(model%stress%tau%xy)
       call parallel_halo(model%stress%tau%scalar)
       call parallel_halo(model%stress%tau%xz)
       call parallel_halo(model%stress%tau%yz)
       call parallel_halo(model%geometry%topg)
       call parallel_halo(model%geometry%thkmask)
       call parallel_halo(model%geometry%marine_bc_normal)
       call parallel_halo(model%velocity%surfvel)
       call parallel_halo(model%ground%gline_flux)
       call parallel_halo(model%velocity%ubas)
       call parallel_halo(model%velocity%vbas)
       call parallel_halo(model%isos%relx)
       call parallel_halo(model%temper%flwa)
       call parallel_halo(model%climate%calving)
       call parallel_halo(model%climate%backstress)
       call parallel_halo(model%geometry%usrf)
       call parallel_halo(model%climate%backstressmap)
       call parallel_halo(model%stress%tau_x)
       call parallel_halo(model%stress%tau_y)
       call parallel_halo(model%geometry%lsrf)
       call parallel_halo(model%temper%temp)

       call horiz_bcs_stag_vector_ew(model%velocity%uvel)
       call horiz_bcs_stag_vector_ns(model%velocity%vvel)
       call horiz_bcs_unstag_scalar(model%geometry%thck)
       call horiz_bcs_unstag_scalar(model%geometry%topg)
       call horiz_bcs_unstag_scalar(model%geometry%thkmask)
       call horiz_bcs_unstag_scalar(model%temper%flwa)
       call horiz_bcs_unstag_scalar(model%temper%temp)
       call t_stopf('simple_glide_halo_upd'//suffix)

     ! Perform parallel operations for restart files
       call t_startf('glissade_post_tstep'//suffix)
       if (tstep_count == 1) call t_adj_detailf(+10)
        call glissade_post_tstep(model)
       if (tstep_count == 1) call t_adj_detailf(-10)
       call t_stopf('glissade_post_tstep'//suffix)

     endif   ! dycore = glide or glissade

     time = time + model%numerics%tinc

     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)

     call t_stopf('glide_tstep'//suffix)
  end do

  call t_stopf('simple glide')

!TODO - Write a glissade_finalise routine.
  ! finalise GLIDE
  call glide_finalise(model)

  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
#if (! defined CCSMCOUPLED && ! defined CESMTIMERS)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
#endif

  call close_log

!TODO - parallel only
  call parallel_finalise

end program simple_glide
