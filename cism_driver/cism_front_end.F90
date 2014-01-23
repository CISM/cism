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

module cism_front_end
  !*FD The CISM front-end is used to connect both the standalone driver
  !*FD (cism_driver) or the CISM interface to CESM (cism_cesm_interface),
  !*FD to the internal and external dycore interface programs.  These are
  !* cism_internal_dycore_interface and cism_external_dycore_interface.

contains

subroutine cism_init_dycore(model)

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

  use cism_external_dycore_interface

!  use glimmer_to_dycore

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

  integer*4 external_dycore_model_index

  integer :: wd
  logical :: do_glide_init

  integer :: tstep_count

    print *,'Entering cism_init_dycore'


  !TODO - call this only for parallel runs?
  ! call parallel_initialise     

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
 
  ! initialise GLIDE
    call t_startf('glide initialization')

  call glide_config(model,config)

  ! This call is needed only if running the EISMINT test cases
  call simple_initialise(climate,config)

  wd = model%options%whichdycore 
!  do_glide_init = (wd == DYCORE_GLIDE) .OR. (wd == DYCORE_BISICLES) .OR. (wd == DYCORE_ALBANYFELIX)
  do_glide_init = (wd == DYCORE_GLIDE)

  if (do_glide_init) then
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

  if ((model%options%whichdycore == DYCORE_BISICLES) .OR. (model%options%whichdycore == DYCORE_ALBANYFELIX)) then
    call cism_init_external_dycore(model%options%external_dycore_type,model)
  endif

  if (model%options%whichdycore .ne. DYCORE_BISICLES) then
  !MJH Created this block here to fill out initial state without needing to enter time stepping loop.  This allows
  ! a run with tend=tstart to be run without time-stepping at all.  It requires solving all diagnostic (i.e. not
  ! time depdendent) variables (most important of which is velocity) for the initial state and then writing the 
  ! initial state as time 0 (or more accurately, as time=tstart).  Also, halo updates need to occur after the 
  ! diagnostic variables are calculated.

  ! ------------- Calculate initial state and output it -----------------

    select case (model%options%whichdycore)
      case (DYCORE_GLIDE,DYCORE_ALBANYFELIX)

        call t_startf('glide_initial_diag_var_solve')
        ! disable further profiling in normal usage
        call t_adj_detailf(+10)

        ! Don't call glide_init_state_diagnostic when running old glide
        ! Instead, start with zero velocity
        if (.not. oldglide) then
          call glide_init_state_diagnostic(model)
        endif

        ! restore profiling to normal settings
        call t_adj_detailf(-10)
        call t_stopf('glide_initial_diag_var_solve')

      case (DYCORE_GLAM)
        call t_startf('glissade_initial_diag_var_solve')
        ! disable further profiling in normal usage
        call t_adj_detailf(+10)

        ! solve the remaining diagnostic variables for the initial state
        call glissade_diagnostic_variable_solve(model)  !velocity, usrf, etc.

        ! restore profiling to normal settings
        call t_adj_detailf(-10)
        call t_stopf('glissade_initial_diag_var_solve')

        case default

    end select

    ! Write initial diagnostic output to log file

    call glide_write_diagnostics(model,        time,       &
                                 tstep_count = tstep_count)

    ! --- Output the initial state -------------

    call t_startf('glide_io_writeall')                                                          
    call glide_io_writeall(model, model, time=time)          ! MJH The optional time argument needs to be supplied 
                                                             !     since we have not yet set model%numerics%time
                                                             !WHL - model%numerics%time is now set above
    call t_stopf('glide_io_writeall')
  end if ! whichdycore .ne. DYCORE_BISICLES

end subroutine cism_init_dycore


subroutine cism_run_dycore(model)

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

  use cism_external_dycore_interface
  

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

  integer*4 :: external_dycore_model_index

!  external_dycore_model_index = this_rank + 1
  external_dycore_model_index = 1

  time = model%numerics%tstart
  tstep_count = 0


  ! ------------- Begin time step loop -----------------
 
  ! run an internal or external dycore, depending on setting external_dycore_type
  do while(time < model%numerics%tend)

    ! TODO: need to fix program to get initialized climate variable
    ! NOTE: these only do something when an EISMINT case is run
    ! call simple_massbalance(climate,model,time)
    ! call simple_surftemp(climate,model,time)
 
    if (model%options%whichdycore /= DYCORE_BISICLES) then
      time = time + model%numerics%tinc
      tstep_count = tstep_count + 1
    endif
! print *,"external_dycore_type: ",model%options%external_dycore_type

    !if (model%options%external_dycore_type .EQ. 0) then      ! NO_EXTERNAL_DYCORE) then
    !  if (model%options%whichdycore == DYCORE_GLIDE) then
    select case (model%options%whichdycore)
      case (DYCORE_GLIDE)
        call t_startf('glide_tstep')

        call t_startf('glide_tstep_p1')
        call glide_tstep_p1(model,time)
        call t_stopf('glide_tstep_p1')

        call t_startf('glide_tstep_p2')
        call glide_tstep_p2(model)
        call t_stopf('glide_tstep_p2')

        call t_startf('glide_tstep_p3')
        call glide_tstep_p3(model)
        call t_stopf('glide_tstep_p3')

        call t_stopf('glide_tstep')

      case (DYCORE_GLAM)
        ! glam/glissade dycore

        call t_startf('glissade_tstep')
        call glissade_tstep(model,time)
        call t_stopf('glissade_tstep')

      case (DYCORE_BISICLES,DYCORE_ALBANYFELIX)
        print *,'Using External Dycore'
        ! The time variable gets incremented within this call:
        call cism_run_external_dycore(model%options%external_dycore_model_index, &
                                      time,model%numerics%tinc)
        ! time = time + model%numerics%tinc
      case default
    end select
    !endif


    ! write ice sheet diagnostics to log file at desired interval (model%numerics%dt_diag)

    call glide_write_diagnostics(model,        time,       &
                                  tstep_count = tstep_count)

    ! Write to output netCDF files at desired intervals

    call t_startf('glide_io_writeall')
    call glide_io_writeall(model,model)
    call t_stopf('glide_io_writeall')
  
  end do   ! time < model%numerics%tend
end subroutine cism_run_dycore

subroutine cism_finalize_dycore(model)

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

  use cism_external_dycore_interface
  

!!  use glimmer_horiz_bcs, only : horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
!!                                horiz_bcs_unstag_scalar, horiz_bcs_stag_scalar

  use glide_stop, only: glide_finalise
  use glide_diagnostics

  implicit none

  type(glide_global_type) :: model        ! model instance
  integer :: clock,clock_rate,ret


  call t_stopf('simple glide')

  ! finalise GLIDE
  call glide_finalise(model)

  call system_clock(clock,clock_rate)
!  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)

#if (! defined CCSMCOUPLED && ! defined CESMTIMERS)
!  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
#endif

  call close_log

  !TODO - call this only for parallel runs?
  ! call parallel_finalise
end subroutine cism_finalize_dycore

end module cism_front_end
