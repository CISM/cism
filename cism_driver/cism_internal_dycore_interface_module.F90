!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_internal_dycore_interface.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module cism_internal_dycore_interface_module

contains


subroutine cism_internal_dycore_interface(whichdycore,model,climate,time,tstep_count)


  use parallel
  use glimmer_global
  use glide
  use glissade
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init
  use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_vector_ew, &
                                horiz_bcs_stag_vector_ns, horiz_bcs_stag_scalar

  use glide_diagnostics

  implicit none

    integer, intent(in) :: whichdycore        ! internale dycore selector
    type(glide_global_type) :: model
    type(simple_climate) :: climate         ! climate
    real(kind=dp), intent(inout) :: time      ! model time in years
    integer :: tstep_count


  real(kind=dp) cur_time, time_inc

  print *,'In internal_dycore_interface'
  print *,'time, tend = ',time,model%numerics%tend

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

     ! Don't call glide_init_state_diagnostic when running old glide
     ! Instead, start with zero velocity

    if (.not. oldglide) then
      call glide_init_state_diagnostic(model)
    endif

    ! restore profiling to normal settings
    call t_adj_detailf(-10)
    call t_stopf('glide_initial_diag_var_solve')

  else   ! glam/glissade dycore
     call t_startf('glissade_initial_diag_var_solve')
      ! disable further profiling in normal usage
      call t_adj_detailf(+10)

     ! solve the remaining diagnostic variables for the initial state
     call glissade_diagnostic_variable_solve(model)  !velocity, usrf, etc.

      ! restore profiling to normal settings
      call t_adj_detailf(-10)
     call t_stopf('glissade_initial_diag_var_solve')
  end if

  ! Write initial diagnostic output to log file

  call glide_write_diagnostics(model,        time,       &
                               tstep_count = tstep_count)

  ! --- Output the initial state -------------

  call t_startf('glide_io_writeall')                                                          
  call glide_io_writeall(model, model, time=time)          ! MJH The optional time argument needs to be supplied 
                                                           !     since we have not yet set model%numerics%time
                                                           !WHL - model%numerics%time is now set above
  call t_stopf('glide_io_writeall')


  ! -------- Begin time step loop ----------------

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

     ! TODO Should we just use model%numerics%time and model%numerics%timecounter?  

     time = time + model%numerics%tinc
     tstep_count = tstep_count + 1

     if (model%options%whichdycore == DYCORE_GLIDE) then

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

     else   ! glam/glissade dycore

       call t_startf('glissade_tstep')
       call glissade_tstep(model,time)
       call t_stopf('glissade_tstep')

     endif   ! glide v. glam/glissade dycore

     ! write ice sheet diagnostics to log file at desired interval (model%numerics%dt_diag)

     call glide_write_diagnostics(model,        time,       &
                                  tstep_count = tstep_count)

     ! Write to output netCDF files at desired intervals

     call t_startf('glide_io_writeall')
     call glide_io_writeall(model,model)
     call t_stopf('glide_io_writeall')

  end do   ! time < model%numerics%tend

end subroutine cism_internal_dycore_interface

end module cism_internal_dycore_interface_module
