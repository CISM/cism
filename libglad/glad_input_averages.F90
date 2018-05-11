!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_input_averages.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glad_input_averages

  !> This module defines a type and related operations for working with inputs from the
  !> GCM. Its main purpose is to produce temporal averages of these inputs.

  ! Note that this module has some functionality in common with glad_mbal_coupling, but
  ! they are used at different stages in the time loops.

  ! The accumulate_averages routine does NOT need to be called every time glad is called,
  ! if time-averaging is occurring at a higher level (e.g., in the coupler). However, the
  ! expectation is that the averaging frequency divides evenly into the mass balance time
  ! step frequency. In addition, each input gets equal weighting, so we are assuming all
  ! inputs apply over the same time interval.

  ! If accumulate_averages is called multiple times within a mass balance time step, then
  ! you cannot restart in the middle of a mass balance time step. (More precisely: You
  ! cannot restart after accumulate_averages has been called one or more times, but
  ! before the mass balance time step occurs.) You can check if it's an okay time to
  ! write a restart file by calling averages_okay_to_restart.

  use glimmer_global, only : dp
  use glimmer_paramets, only: GLC_DEBUG, stdout
  use glimmer_log
  use parallel, only : main_task
  
  implicit none
  private
  
  type, public :: glad_input_averages_type
     private

     integer :: av_start_time = 0  ! Value of time from the last occasion averaging was restarted (hours)
     integer :: av_steps      = 0  ! Number of times glimmer has been called in current round of averaging
     
     real(dp),pointer,dimension(:,:) :: tot_qsmb => null()  ! running total surface mass balance (kg m-2 s-1)
     real(dp),pointer,dimension(:,:) :: tot_tsfc => null()  ! running total surface temperature (deg C)
  
  end type glad_input_averages_type

  public :: initialize_glad_input_averages
  public :: get_av_start_time
  public :: accumulate_averages
  public :: calculate_averages
  public :: reset_glad_input_averages
  public :: averages_okay_to_restart

contains

  subroutine initialize_glad_input_averages(glad_inputs, ewn, nsn, next_av_start)
    ! Initialize a glad_inputs instance
    type(glad_input_averages_type), intent(inout) :: glad_inputs

    ! dimensions of local grid
    integer, intent(in) :: ewn
    integer, intent(in) :: nsn

    ! Starting time of next averaging period (hours)
    integer, intent(in) :: next_av_start

    allocate(glad_inputs%tot_qsmb(ewn,nsn));  glad_inputs%tot_qsmb = 0.d0
    allocate(glad_inputs%tot_tsfc(ewn,nsn));  glad_inputs%tot_tsfc = 0.d0

    glad_inputs%av_start_time = next_av_start
  end subroutine initialize_glad_input_averages

  integer function get_av_start_time(glad_inputs)
    ! Get value of time from the last occasion averaging was restarted (hours)
    type(glad_input_averages_type), intent(in) :: glad_inputs

    get_av_start_time = glad_inputs%av_start_time
  end function get_av_start_time
    
  subroutine accumulate_averages(glad_inputs, qsmb, tsfc, time)
    ! Accumulate averages based on one set of inputs.
    ! Should be called every time we have new inputs from the climate model.
    type(glad_input_averages_type), intent(inout) :: glad_inputs
    real(dp),dimension(:,:),intent(in)  :: qsmb     ! flux of glacier ice (kg/m^2/s)
    real(dp),dimension(:,:),intent(in)  :: tsfc     ! surface ground temperature (C)
    integer, intent(in) :: time  ! Current model time
    
    glad_inputs%tot_qsmb(:,:) = glad_inputs%tot_qsmb(:,:) + qsmb(:,:)
    glad_inputs%tot_tsfc(:,:) = glad_inputs%tot_tsfc(:,:) + tsfc(:,:)
    glad_inputs%av_steps = glad_inputs%av_steps + 1
  end subroutine accumulate_averages

  subroutine calculate_averages(glad_inputs, qsmb, tsfc)
    ! Calculate averages over the averaging period
    type(glad_input_averages_type), intent(in) :: glad_inputs
    real(dp), dimension(:,:), intent(out) :: qsmb  ! average surface mass balance (kg m-2 s-1)
    real(dp), dimension(:,:), intent(out) :: tsfc  ! average surface temperature (deg C)
    
    qsmb(:,:) = glad_inputs%tot_qsmb(:,:) / real(glad_inputs%av_steps,dp)
    tsfc(:,:) = glad_inputs%tot_tsfc(:,:) / real(glad_inputs%av_steps,dp)
  end subroutine calculate_averages

  subroutine reset_glad_input_averages(glad_inputs, next_av_start)
    ! Resets this glad_inputs instance
    ! Should be called at the end of an averaging period, in order to prepare for the
    ! next averaging period
    type(glad_input_averages_type), intent(inout) :: glad_inputs
    integer, intent(in) :: next_av_start  ! start time for next averaging period (hours)

    glad_inputs%tot_qsmb(:,:) = 0.d0
    glad_inputs%tot_tsfc(:,:) = 0.d0
    glad_inputs%av_steps      = 0
    glad_inputs%av_start_time = next_av_start
  end subroutine reset_glad_input_averages

  pure logical function averages_okay_to_restart(glad_inputs)
    ! Returns true if this is an okay time to write a restart file based on these
    ! glad_inputs, false if not.
    !
    ! It is not an okay time to write a restart file if we have some accumulated averages,
    ! since we currently do not write these partial averages to restart files.
    type(glad_input_averages_type), intent(in) :: glad_inputs

    averages_okay_to_restart = (glad_inputs%av_steps == 0)
  end function averages_okay_to_restart

end module glad_input_averages
