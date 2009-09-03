
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_paramets.f90 - part of the GLIMMER ice model     + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module paramets

  use glimmer_global, only : sp, dp
  use physcon, only : scyr

  implicit none; save

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

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_PARAMETS
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_PARAMETS
!MH!#endif
!MH!
!MH!#ifdef RESTARTS
!MH!contains
!MH!
!MH!#define RST_PARAMETS
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_PARAMETS
!MH!#endif

end module paramets
