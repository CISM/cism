
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                            +
! +  glimmer_paramets.f90 - part of the Glimmer-CISM ice model + 
! +                                                            +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
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

!> model scaling constants
module glimmer_paramets

  use glimmer_global, only : sp, dp
  use glimmer_physcon, only : scyr

  implicit none; save

!lipscomb - TO DO - redundant output units  
!           It is redundant to define both stdout (which is public) and 
!            glimmer_unit (which is private to glimmer_log.F90).
!           However, it is sometimes convenient to write to stdout in Glimmer
!            without calling write_log.  
!           May want to delete this later (and declare stdout in glc_constants 
!            for CESM runs).

  integer :: stdout = 6

! some parameters for debugging and diagnostics
   integer, parameter ::   &
      itest = 133, jtest = 84,  &          ! in Greenland (FV2), lat 67.3 N, lon 330 E
                  jjtest = 97 - jtest,  &  ! reversed for N to S indexing (FV2, ny = 96)
      itest_local = 60, jtest_local = 54   ! Greenland 20 deg grid, initial usrf = 491 m

   integer, parameter :: idiag = 30, jdiag = 50  ! point for diagnostic output

! logical flag to turn on special DEBUG output (related to test points), false by default
   logical :: GLC_DEBUG = .false.

! scaling parameters
! TO DO - Remove these at some point?

  real(dp), parameter :: thk0 = 2000.0d0          ! m
  real(dp), parameter :: len0 = 200.0d3        ! m
  real(dp), parameter :: vel0 = 500.0 / scyr    ! m yr^{-1} converted to S.I. units
  !real(dp), parameter :: vis0 = 5.70d-18 / scyr  ! yr^{-1} Pa^{-3} converted to S.I. units
  real(dp), parameter :: vis0 = 1d-16 / scyr 
  real(dp), parameter :: acc0 = thk0 * vel0 / len0  ! m s^{-1} 
  ! ** for zero order model real(dp), parameter :: tim0 = thk0 / acc0      ! s
  real(dp), parameter :: tim0 = len0 / vel0      ! s
  real(dp) :: tau0                        ! Pa note cannot define here as f90 wont allow
                                          ! parameters with noninteger powers in - look
                                          ! in initial in blah.f90 (not sure this applies now...)

  real(sp), parameter :: conv = tim0 / scyr

end module glimmer_paramets
