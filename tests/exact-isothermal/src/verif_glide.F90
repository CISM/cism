! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  verif_glide.F90 - part of the Glimmer-CISM ice model     + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2006, 2009
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

! testing steady state

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program verifglide

  ! load various modules
  use glimmer_global, only:rk ! precision of the model
  use glide                   ! main glide module
  use glimmer_log             ! module for logging messages
  use glimmer_config          ! module for handling configuration files
  use verif
  use verif_io
  use glimmer_commandline
  use glimmer_writestats
  implicit none

  ! some variables
  type(glide_global_type) :: model    
  type(ConfigSection), pointer :: config

  type(verif_type) :: veri

  real(kind=rk) time ! current time
  real(kind=dp) t1,t2
  integer clock,clock_rate
  
  call glimmer_GetCommandline()
  
  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))
  
  ! read configuration
  call ConfigRead(commandline_configname,config)
  call glide_config(model,config)
  call verif_config(config,veri)
  call verif_printconfig(veri)
  call CheckSections(config)

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise test setup
  call verif_init(model, veri)
  
  ! initialise GLIDE
  call glide_initialise(model)
  ! fill dimension variables
  ! create verif variables
  call verif_io_createall(model)
  call glide_nc_fillall(model)
  ! get current time from start time
  time = get_tstart(model)

  ! initial conditions
  call verif_update(model, veri, time)
  call verif_initthk(model, veri)

  ! loop over times
  do while(time.le.model%numerics%tend)
     call verif_update(model, veri, time)

     ! calculate temperature and velocity distribution
     call glide_tstep_p1(model,time)
     call verif_io_writeall(veri,model)
     ! write to netCDF file, move ice
     call glide_tstep_p2(model)
     ! calculate isostatic adjustment
     call glide_tstep_p3(model)
     ! increment time counter
     time = time + get_tinc(model)
  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log

end program verifglide
