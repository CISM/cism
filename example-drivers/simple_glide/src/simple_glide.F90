!TODO - Change name to cism_driver?

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
!       Call it cism_driver?

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
!!  use glimmer_horiz_bcs, only : horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
!!                                horiz_bcs_unstag_scalar, horiz_bcs_stag_scalar

  use glissade

!TODO - Change to 'cism_diagnostics'?
  use glide_diagnostics

  implicit none

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
!WHLTSTEP - Changed time to dp
!  real(kind=rk) time
  real(kind=dp) :: time                   ! model time in years
  real(kind=dp) :: t1,t2
  integer :: clock,clock_rate,ret

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
! MJH says: simple_initialise is only needed when using simple_massbalance and simple_surftemp.  
! If those are removed, this should be removed too.
  call simple_initialise(climate,config)

  if (model%options%whichdycore == DYCORE_GLIDE) then
     call glide_initialise(model)
  else       ! glam/glissade dycore	
     call glissade_initialise(model)
  endif

  call CheckSections(config)

!TODO - Change to glimmer_nc_fillall?
  ! fill dimension variables on output files
  call glide_nc_fillall(model)

  time = model%numerics%tstart
  tstep_count = 0
  model%numerics%time = time    ! MJH added 1/10/13 - the initial diagnostic glissade solve won't know the correct time on a restart unless we set it here.

!### !TODO - Are these calls always needed?
!### !TODO MJH These aren't needed for the initial state - only during time-stepping.  
     !         So delete the calls here and only call them in the time-stepping routine
!### !  call simple_massbalance(climate,model,time)
!### !  call simple_surftemp(climate,model,time)

  call spinup_lithot(model)
    call t_stopf('glide initialization')

!###   !tstep_count = 0 !TODO MJH delete this line - moved up above

  !MJH Created this block here to fill out initial state without needing to enter time stepping loop.  This allows
  ! a run with tend=tstart to be run without time-stepping at all.  It requires solving all diagnostic (i.e. not
  ! time depdendent) variables (most important of which is velocity) for the initial state and then writing the 
  !initial state as time 0 (or more accurately, as time=tstart).  Also, halo updates need to occur after the diagnostic
  ! variables are calculated.

  ! ------------- Calculate initial state and output it -----------------
  if (model%options%whichdycore == DYCORE_GLIDE) then
     call t_startf('glide_initial_diag_var_solve')
      ! disable further profiling in normal usage
      call t_adj_detailf(+10)

     call glide_init_state_diagnostic(model)

      ! restore profiling to normal settings
      call t_adj_detailf(-10)
     call t_stopf('glide_initial_diag_var_solve')
  else   ! glissade dycore
     call t_startf('glissade_initial_diag_var_solve')
      ! disable further profiling in normal usage
      call t_adj_detailf(+10)

     ! solve the remaining diagnostic variables for the initial state
     call glissade_diagnostic_variable_solve(model)  !velocity, usrf, etc.
     ! TODO HALO UPDATES?  Hopefully these are done in the subroutine

      ! restore profiling to normal settings
      call t_adj_detailf(-10)
     call t_stopf('glissade_initial_diag_var_solve')
  end if

!WHL - Write initial diagnostic output to log file
  call glide_write_diag(model, time, model%numerics%idiag_global,  &
                                     model%numerics%jdiag_global)

  ! --- Output the initial state -------------
  ! TODO MJH Copied this below from glissade_post_tstep().  May want to make a subroutine that just has this 
  !block in it.  It could be called glimmer_write_output and be in simple_glide if it can be used by both glide
  ! and glissade.  Or else separate routines at the glissade/glide module level.
  !TODO - Change to glimmer_io_writeall?
  !TODO - the write operation in post_step is inside an if-construct that checks an optional 'nowrite' logical variable. 
  ! However the call to that subroutine at the end of this module does not supply the optional variable.  Therefore I am 
  !leaving out that if-construct here.  If simple_glide actually does support a nowrite option, then a check for it would
  ! need to occur here!
  call t_startf('glide_io_writeall')
!WHLTSTEP - Changed time to dp
!  call glide_io_writeall(model,model, time=REAL(time,4))  
                                                          
  call glide_io_writeall(model, model, time=time)          ! MJH The optional time argument needs to be supplied 
                                                           !     since we have not yet set model%numerics%time
  call t_stopf('glide_io_writeall')
!TODO Do we want to write out litho output here too?  I think so, but I'm not sure what it is exactly.
  if (model%options%gthf > 0) then
     call t_startf('glide_lithot_io_writeall')
!WHLTSTEP - Changed time to dp
!     call glide_lithot_io_writeall(model,model, time=REAL(time,4)) 
                                                                    
     call glide_lithot_io_writeall(model, model, time=time)          ! MJH The optional time argument needs to be supplied
                                                                     !     since we have not yet set model%numerics%time

     call t_stopf('glide_lithot_io_writeall')
  endif 

!==============================
! if-block here is for forcing glide code to match old glide code.  It is left here 
! just for testing purposes, but otherwise should not be used.  In the old code,
! velocity fields were all 0 on the first temperature solve.  In the new 
! organization of time-stepping there is a diagnostic solve for the initial state
! prior to time-stepping, which results in velocity being non-zero on the first 
! temperature solve.  This changes the result due to dissipation and sliding heating.
  !if (model%options%whichdycore == DYCORE_GLIDE) then
  ! zero out velocity fields for GLIDE dycore if you want to get same answers as old versions of the code.  
  ! Zeroing all of these might not be necessary.
  !  model%velocity%uvel = 0.0
  !  model%velocity%vvel = 0.0
  !  model%velocity%velnorm = 0.0
  !  model%velocity%wvel = 0.0
  !  model%velocity%wgrd = 0.0
  !  model%velocity%surfvel = 0.0
  !  model%velocity%uflx = 0.0
  !  model%velocity%vflx = 0.0
  !  model%velocity%diffu = 0.0
  !  model%velocity%ubas = 0.0
  !  model%velocity%vbas = 0.0
  !  model%velocity%btrc = 0.0
  !endif
!==============================

  ! ------------- Begin time step loop -----------------
!### !  do while(time <= model%numerics%tend)
  do while(time < model%numerics%tend)

     ! --- First assign forcings ----
     ! Because this is Forward Euler, the forcings should be from the previous time step (e.g. H1 = f(H0, V0, SMB0))
     ! TODO Write generic forcing subroutines that could call simple_massbalance/surftemp or some other forcing module.  
     ! simple_massbalance/surftemp are only used for the EISMINT experiments.  
     ! If they are called without EISMINT options in the config file, they will do nothing.
     ! TODO May want to move them to the glissade/glide time steppers.  
     ! If so be careful about using the right time level (i.e. the old one).
     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)


     ! --- Increment time before performing time step operations ---
     ! We are solving variables at the new time level using values from the previous time level.
     ! TODO Can we just use model%numerics%time and model%numerics%timecounter?  
     !      That would be less confusing than having time-keeping done separately here and in glissade.
     time = time + model%numerics%tinc
     tstep_count = tstep_count + 1

     call t_startf('glide_tstep')

!TODO - Add subroutine glide_step that calls glide_tstep_p1/p2/p3?

     if (model%options%whichdycore == DYCORE_GLIDE) then

       call t_startf('glide_tstep_p1')
       call glide_tstep_p1(model,time)
       call t_stopf('glide_tstep_p1')

       call t_startf('glide_tstep_p2')
       call glide_tstep_p2(model)
       call t_stopf('glide_tstep_p2')

       call t_startf('glide_tstep_p3')
       call glide_tstep_p3(model)
       call t_stopf('glide_tstep_p3')

     else   ! glam/glissade dycore

!TODO - Make this a single subroutine call

       call t_startf('glissade_tstep')
       call glissade_tstep(model,time)
       call t_stopf('glissade_tstep')

     endif   ! glide v. glam/glissade dycore

     ! override masking stuff for now  !TODO - What does this mean?

!TODO - Change to cism_write_diag?
!TODO Can this be moved to end of time loop so that we can merge the two if-statements into a single construct?
     ! write ice sheet diagnostics
     if (mod(tstep_count, model%numerics%ndiag)==0)  then
        call glide_write_diag(model, time, model%numerics%idiag_global,  &
                                           model%numerics%jdiag_global)
     endif

!TODO - Remove this comment?
     ! Redistribute calls here to spread the data back out.

     if (model%options%whichdycore == DYCORE_GLIDE) then

     ! Perform parallel operations for restart files
       call t_startf('glide_tstep_postp3')
       call glide_tstep_postp3(model)
       call t_stopf('glide_tstep_postp3')

     else   ! glam/glissade dycore

!TODO - I think we can safely remove all these parallel halo calls.
!       Before doing so, make sure we have the required calls in glissade.F90.

!WHL - Replaced some parallel_halo calls with staggered_parallel_halo as appropriate       
       call t_startf('simple_glide_halo_upd')

!WHL - Commented out parallel_halo calls here.
!      Whatever halo updates are needed should be done within the appropriate dycore.
!!       call parallel_halo(model%stress%efvs)
!!       call horiz_bcs_unstag_scalar(model%stress%efvs)
!!       call staggered_parallel_halo(model%velocity%uvel)
!!       call horiz_bcs_stag_vector_ew(model%velocity%uvel)
!!       call staggered_parallel_halo(model%velocity%vvel)
!!       call horiz_bcs_stag_vector_ns(model%velocity%vvel)
!!       call staggered_parallel_halo(model%velocity%uflx)
!!       call horiz_bcs_stag_vector_ew(model%velocity%uflx)
!!       call staggered_parallel_halo(model%velocity%vflx)
!!       call horiz_bcs_stag_vector_ns(model%velocity%vflx)
!!       call staggered_parallel_halo(model%velocity%velnorm)
!!       call horiz_bcs_stag_scalar(model%velocity%velnorm)
!!       call parallel_halo(model%geometry%thck)
!!       call horiz_bcs_unstag_scalar(model%geometry%thck)
!!       call staggered_parallel_halo(model%geomderv%stagthck)
!!       call horiz_bcs_stag_scalar(model%geomderv%stagthck)
!!       call parallel_halo(model%climate%acab)
!!       call horiz_bcs_unstag_scalar(model%climate%acab)
!!       call staggered_parallel_halo(model%geomderv%dusrfdew)
!!       call horiz_bcs_stag_vector_ew(model%geomderv%dusrfdew)
!!       call staggered_parallel_halo(model%geomderv%dusrfdns)
!!       call horiz_bcs_stag_vector_ns(model%geomderv%dusrfdns)
!!       call staggered_parallel_halo(model%geomderv%dthckdew)
!!       call horiz_bcs_stag_vector_ew(model%geomderv%dthckdew)
!!       call staggered_parallel_halo(model%geomderv%dthckdns)
!!       call horiz_bcs_stag_vector_ns(model%geomderv%dthckdns)
!!       call parallel_halo(model%stress%tau%xx)
!!       call horiz_bcs_unstag_scalar(model%stress%tau%xx)
!!       call parallel_halo(model%stress%tau%yy)
!!       call horiz_bcs_unstag_scalar(model%stress%tau%yy)
!!       call parallel_halo(model%stress%tau%xy)
!!       call horiz_bcs_unstag_scalar(model%stress%tau%xy)
!!       call parallel_halo(model%stress%tau%scalar)
!!       call horiz_bcs_unstag_scalar(model%stress%tau%scalar)
!!       call parallel_halo(model%stress%tau%xz)
!!       call horiz_bcs_unstag_scalar(model%stress%tau%xz)
!!       call parallel_halo(model%stress%tau%yz)
!!       call horiz_bcs_unstag_scalar(model%stress%tau%yz)
!!       call parallel_halo(model%geometry%topg)
!!       call horiz_bcs_unstag_scalar(model%geometry%topg)
!!       call parallel_halo(model%geometry%thkmask)
!!       call horiz_bcs_unstag_scalar(model%geometry%thkmask)
!!       call parallel_halo(model%geometry%marine_bc_normal)
!!       call horiz_bcs_unstag_scalar(model%geometry%marine_bc_normal)
!!       call staggered_parallel_halo(model%velocity%surfvel)
!!       call horiz_bcs_stag_scalar(model%velocity%surfvel)
!!       call parallel_halo(model%ground%gline_flux)
!!       call horiz_bcs_unstag_scalar(model%ground%gline_flux)
!!       call staggered_parallel_halo(model%velocity%ubas)
!!       call horiz_bcs_stag_vector_ew(model%velocity%ubas)
!!       call staggered_parallel_halo(model%velocity%vbas)
!!       call horiz_bcs_stag_vector_ns(model%velocity%vbas)
!!       call parallel_halo(model%isos%relx)
!!       call horiz_bcs_unstag_scalar(model%isos%relx)
!!       call parallel_halo(model%temper%flwa)
!!       call horiz_bcs_unstag_scalar(model%temper%flwa)
!!       call parallel_halo(model%climate%calving)
!!       call horiz_bcs_unstag_scalar(model%climate%calving)
!!       call parallel_halo(model%climate%backstress)
!!       call horiz_bcs_unstag_scalar(model%climate%backstress)
!!       call parallel_halo(model%geometry%usrf)
!!       call horiz_bcs_unstag_scalar(model%geometry%usrf)
!!       call parallel_halo(model%climate%backstressmap)
!!       call horiz_bcs_unstag_scalar(model%climate%backstressmap)
!!       call parallel_halo(model%stress%tau_x)
!!       call horiz_bcs_stag_vector_ew(model%stress%tau_x)
!!       call parallel_halo(model%stress%tau_y)
!!       call horiz_bcs_stag_vector_ew(model%stress%tau_y)
!!       call parallel_halo(model%geometry%lsrf)
!!       call horiz_bcs_unstag_scalar(model%geometry%lsrf)
!!       call parallel_halo(model%temper%temp)
!!       call horiz_bcs_unstag_scalar(model%temper%temp)

       call t_stopf('simple_glide_halo_upd')

     ! Perform parallel operations for restart files
     ! TODO post_tstep has nothing left in it but the output writing.  
     !      I think it should be eliminated and the output writing occur here in simple_glide.
       call t_startf('glissade_post_tstep')
       call glissade_post_tstep(model)
       call t_stopf('glissade_post_tstep')

     endif   ! glide v. glam/glissade dycore

     call t_stopf('glide_tstep')
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
