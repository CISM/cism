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

program simple_bisicles
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
  use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
                               horiz_bcs_stag_scalar

  use glide_diagnostics

  use glimmer_to_dycore

  implicit none

#ifdef GPTL
#include <gptl.inc>
#endif
#ifdef PAPI
#include <f90papi.h>
#endif

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

  real(kind=sp) cur_time, time_inc

  integer :: tstep_count

  ! for external dycore:
  integer*4 dycore_model_index
  integer argc
print *,"Calling parallel_initialize "
  call parallel_initialise
print *,"Returned."

  ! start gptl
#ifdef GPTL
  ret = gptlsetoption (gptlprint_method,gptlfull_tree)
  ret = gptlsetoption (PAPI_FP_OPS, 1)
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('total')
#endif

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

  ! initialise GLIDE
  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart
  time_inc = model%numerics%tinc

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)
print *,"Calling gtd_init_dycore_interface "
  call gtd_init_dycore_interface()
print *,"Calling gtd_init_dycore "
  call gtd_init_dycore(model,dycore_model_index)
print *,"Returned."
print *,"Calling gtd_run_dycore "
  call gtd_run_dycore(dycore_model_index,cur_time,time_inc)
print *,"Returned."

  tstep_count = 0
print *,"Entering while loop"
  do while(time.le.model%numerics%tend)
     call glide_tstep_p1(model,time)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     ! override masking stuff for now

     tstep_count = tstep_count + 1

     ! Redistribute calls here to spread the data back out.
     call parallel_halo(model%stress%efvs)
     call horiz_bcs_unstag_scalar(model%stress%efvs)
     call parallel_halo(model%velocity%uvel)
     call horiz_bcs_stag_vector_ew(model%velocity%uvel)
     call parallel_halo(model%velocity%vvel)
     call horiz_bcs_stag_vector_ns(model%velocity%vvel)
     call parallel_halo(model%velocity%uflx)
     call horiz_bcs_stag_vector_ew(model%velocity%uflx)
     call parallel_halo(model%velocity%vflx)
     call horiz_bcs_stag_vector_ns(model%velocity%vflx)
     call parallel_halo(model%velocity%velnorm)
     call horiz_bcs_stag_scalar(model%velocity%velnorm)
     call parallel_halo(model%geometry%thck)
     call horiz_bcs_unstag_scalar(model%geometry%thck)
     call parallel_halo(model%geomderv%stagthck)
     call horiz_bcs_stag_scalar(model%geomderv%stagthck)
     call parallel_halo(model%climate%acab)
     call horiz_bcs_unstag_scalar(model%climate%acab)
     call parallel_halo(model%geomderv%dusrfdew)
     call horiz_bcs_stag_vector_ew(model%geomderv%dusrfdew)
     call parallel_halo(model%geomderv%dusrfdns)
     call horiz_bcs_stag_vector_ns(model%geomderv%dusrfdns)
     call parallel_halo(model%geomderv%dthckdew)
     call horiz_bcs_stag_vector_ew(model%geomderv%dthckdew)
     call parallel_halo(model%geomderv%dthckdns)
     call horiz_bcs_stag_vector_ns(model%geomderv%dthckdns)
     call parallel_halo(model%stress%tau%xx)
     call horiz_bcs_unstag_scalar(model%stress%tau%xx)
     call parallel_halo(model%stress%tau%yy)
     call horiz_bcs_unstag_scalar(model%stress%tau%yy)
     call parallel_halo(model%stress%tau%xy)
     call horiz_bcs_unstag_scalar(model%stress%tau%xy)
     call parallel_halo(model%stress%tau%scalar)
     call horiz_bcs_unstag_scalar(model%stress%tau%scalar)
     call parallel_halo(model%stress%tau%xz)
     call horiz_bcs_unstag_scalar(model%stress%tau%xz)
     call parallel_halo(model%stress%tau%yz)
     call horiz_bcs_unstag_scalar(model%stress%tau%yz)
     call parallel_halo(model%geometry%topg)
     call horiz_bcs_unstag_scalar(model%geometry%topg)
     call parallel_halo(model%geometry%thkmask)
     call horiz_bcs_unstag_scalar(model%geometry%thkmask)
     call parallel_halo(model%geometry%marine_bc_normal)
     call horiz_bcs_unstag_scalar(model%geometry%marine_bc_normal)
     call parallel_halo(model%velocity%surfvel)
     call horiz_bcs_stag_scalar(model%velocity%surfvel)
     call parallel_halo(model%ground%gline_flux)
     call horiz_bcs_unstag_scalar(model%ground%gline_flux)
     call parallel_halo(model%velocity%ubas)
     call horiz_bcs_stag_vector_ew(model%velocity%ubas)
     call parallel_halo(model%velocity%vbas)
     call horiz_bcs_stag_vector_ns(model%velocity%vbas)
     call parallel_halo(model%isos%relx)
     call horiz_bcs_unstag_scalar(model%isos%relx)
     call parallel_halo(model%temper%flwa)
     call horiz_bcs_unstag_scalar(model%temper%flwa)
     call parallel_halo(model%climate%calving)
     call horiz_bcs_unstag_scalar(model%climate%calving)
     call parallel_halo(model%climate%backstress)
     call horiz_bcs_unstag_scalar(model%climate%backstress)
     call parallel_halo(model%geometry%usrf)
     call horiz_bcs_unstag_scalar(model%geometry%usrf)
     call parallel_halo(model%climate%backstressmap)
     call horiz_bcs_unstag_scalar(model%climate%backstressmap)
     call parallel_halo(model%stress%tau_x)
     call horiz_bcs_stag_vector_ew(model%stress%tau_x)
     call parallel_halo(model%stress%tau_y)
     call horiz_bcs_stag_vector_ew(model%stress%tau_y)
     call parallel_halo(model%geometry%lsrf)
     call horiz_bcs_unstag_scalar(model%geometry%lsrf)

     if (model%options%whichtemp == TEMP_REMAP_ADV) then
        call parallel_halo(model%temper%temp)
        call horiz_bcs_unstag_scalar(model%temper%temp)
     else
        call parallel_halo_temperature(model%temper%temp)
        call horiz_bcs_unstag_scalar(model%temper%temp)
     endif

     ! Perform parallel operations for restart files
     call glide_tstep_postp3(model)

     time = time + model%numerics%tinc
     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)     
  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log
  call parallel_finalise

  ! stop gptl
#ifdef GPTL
  ret = gptlstop ('total')
  ret = gptlpr (0)
  ret = gptlfinalize ()
#endif
  
end program simple_bisicles
