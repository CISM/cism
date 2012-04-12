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

  use glide_diagnostics

  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

  integer :: tstep_count

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

#if (defined CCSMCOUPLED || defined CESMTIMERS)
  ! initialise profiling
  call glide_prof_init(model)
#endif

  call t_startf('simple glide')

  ! initialise GLIDE
  call t_startf('glide initialization')
  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)
  call t_stopf('glide initialization')

  tstep_count = 0

  do while(time.le.model%numerics%tend)
    call t_startf('glide_tstep_p1')
     call glide_tstep_p1(model,time)
    call t_stopf('glide_tstep_p1')

    call t_startf('glide_tstep_p2')
     call glide_tstep_p2(model)
    call t_stopf('glide_tstep_p2')

    call t_startf('glide_tstep_p3')
     call glide_tstep_p3(model)
    call t_stopf('glide_tstep_p3')

     ! override masking stuff for now
     tstep_count = tstep_count + 1

     ! Redistribute calls here to spread the data back out.
    call t_startf('simple_glide_halo_upd')
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

     if (model%options%whichtemp == TEMP_REMAP_ADV) then
        call parallel_halo(model%temper%temp)
     else
        call parallel_halo_temperature(model%temper%temp)
     endif
    call t_stopf('simple_glide_halo_upd')

     ! Perform parallel operations for restart files
    call t_startf('glide_tstep_postp3')
     call glide_tstep_postp3(model)
    call t_stopf('glide_tstep_postp3')

     time = time + model%numerics%tinc

    call t_startf('simple_massbalance')
     call simple_massbalance(climate,model,time)
    call t_stopf('simple_massbalance')

    call t_startf('simple_surftemp')
     call simple_surftemp(climate,model,time)     
    call t_stopf('simple_surftemp')
  end do

  call t_stopf('simple glide')

  ! finalise GLIDE
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log
  call parallel_finalise

end program simple_glide
