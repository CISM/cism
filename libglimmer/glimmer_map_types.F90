! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_map_types.f90 - part of the GLIMMER ice model    + 
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

module glimmer_map_types

  !*FD New combined map-projection handling code for GLIMMER:
  !*FD This module contains derived types.

  use glimmer_global, only: rk

  implicit none

  type glimmap_proj
     logical :: found = .false.
     type(proj_laea),  pointer :: laea  => NULL() !*FD Pointer to Lambert azimuthal equal area type
     type(proj_aea),   pointer :: aea   => NULL() !*FD Pointer to Albers equal area conic type
     type(proj_lcc),   pointer :: lcc   => NULL() !*FD Pointer to Lambert conic conformal type
     type(proj_stere), pointer :: stere => NULL() !*FD Pointer to Stereographic type
  end type glimmap_proj

  !-------------------------------------------------------------

  type proj_laea
     !*FD Lambert Azimuthal Equal Area
     real(rk) :: longitude_of_central_meridian
     real(rk) :: latitude_of_projection_origin
     real(rk) :: false_easting
     real(rk) :: false_northing
     real(rk) :: sinp    !*FD Sine of latitude_of_projection_origin
     real(rk) :: cosp    !*FD Cosine of latitude_of_projection_origin
     integer :: pole !*FD Set to 1 for N pole, -1 for S pole, 0 otherwise
  end type proj_laea

  !-------------------------------------------------------------

  type proj_aea
     !*FD Albers Equal-Area Conic
     real(rk),dimension(2) :: standard_parallel
     real(rk) :: longitude_of_central_meridian
     real(rk) :: latitude_of_projection_origin
     real(rk) :: false_easting
     real(rk) :: false_northing
     real(rk) :: rho0 !*FD Convenience constant
     real(rk) :: rho0_R !*FD Convenience constant (is rho0/EQ_RAD)
     real(rk) :: c    !*FD Convenience constant
     real(rk) :: n    !*FD Convenience constant
     real(rk) :: i_n  !*FD Convenience constant (inverse of n)
  end type proj_aea

  !-------------------------------------------------------------

  type proj_lcc
     !*FD Lambert Conic Conformal
     real(rk),dimension(2) :: standard_parallel
     real(rk) :: longitude_of_central_meridian
     real(rk) :: latitude_of_projection_origin
     real(rk) :: false_easting
     real(rk) :: false_northing
     real(rk) :: rho0 !*FD Convenience constant
     real(rk) :: f    !*FD Convenience constant
     real(rk) :: n    !*FD Convenience constant
     real(rk) :: i_n  !*FD Convenience constant (inverse of n)
  end type proj_lcc

  !-------------------------------------------------------------

  type proj_stere
     !*FD Stereographic projection derived type 
     real(rk) :: longitude_of_central_meridian          
     real(rk) :: latitude_of_projection_origin          
     real(rk) :: scale_factor_at_proj_origin = 0. 
     real(rk) :: standard_parallel = 0.                 
     real(rk) :: false_easting          
     real(rk) :: false_northing
     integer :: pole  !*FD Set to 1 for N pole, -1 for S pole, 0 otherwise
     logical :: equatorial !*FD Set true if equatorial aspect
     real(rk) :: k0 !*FD scale factor or std par converted to scale factor
     real(rk) :: ik0 !*FD inverse of k0
     real(rk) :: sinp,cosp !*FD sin and cos of latitude_of_projection_origin
  end type proj_stere

  ! Global mapping parameters ----------------------------------

  real(rk),parameter :: pi         = 3.141592654    !*FD The value of $\pi$.
  real(rk),parameter :: M_PI_4     = pi/4           !*FD The value of $\pi/4$.
  real(rk),parameter :: M_PI_2     = pi/2           !*FD The value of $\pi/2$.
  real(rk),parameter :: D2R        = pi/180.0       !*FD Degrees-to-radians conversion factor.
  real(rk),parameter :: R2D        = 180.0/pi       !*FD Radians-to-degrees conversion factor.
  real(rk),parameter :: EQ_RAD     = 6.37e6         !*FD Radius of the earth (m)
  real(rk),parameter :: i_EQ_RAD   = 1.0_rk/EQ_RAD  !*FD Inverse radius of the earth (m^-1)
  real(rk),parameter :: CONV_LIMIT = 1.0e-8         !*FD Convergence limit (a small number).

  integer, parameter :: GMAP_LAEA=1
  integer, parameter :: GMAP_AEA=2
  integer, parameter :: GMAP_LCC=3
  integer, parameter :: GMAP_STERE=4

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_MAP_TYPES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_MAP_TYPES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_MAP_TYPES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_MAP_TYPES
!MH!#endif

  function glimmap_allocated(proj)

    !*FD return true if structure contains a known projection

    implicit none
    type(glimmap_proj) :: proj
    logical glimmap_allocated

    glimmap_allocated = proj%found
  end function glimmap_allocated

  subroutine glimmap_diag(proj)

    !*FD This is incomplete diagnostics code to output full
    !*FD content of projection type. Only does
	!*FD Stereographic projections so far.

    use glimmer_log

    type(glimmap_proj) :: proj

    if (associated(proj%stere)) then
       call glimmap_diag_stere(proj%stere)
    else
       call write_log('Stereographic projection not found')
    end if

  end subroutine glimmap_diag

  subroutine glimmap_diag_stere(params)

    use glimmer_log

    type(proj_stere) :: params
    character(80) :: message
    
    call write_log('***** Stereographic *****')
    write(message,*)'longitude_of_central_meridian:', params%longitude_of_central_meridian
    call write_log(message)
    write(message,*)'latitude_of_projection_origin:', params%latitude_of_projection_origin
    call write_log(message)
    write(message,*)'scale_factor_at_proj_origin:', params%scale_factor_at_proj_origin
    call write_log(message)
    write(message,*)'standard_parallel:', params%standard_parallel
    call write_log(message)
    write(message,*)'false_easting:', params%false_easting
    call write_log(message)
    write(message,*)'false_northing:', params%false_northing
    call write_log(message)
    write(message,*)'pole:', params%pole
    call write_log(message)
    write(message,*)'equatorial:', params%equatorial
    call write_log(message)
    write(message,*)'k0:', params%k0
    call write_log(message)
    write(message,*)'ik0:', params%ik0
    call write_log(message)
    write(message,*)'sinp:', params%sinp
    call write_log(message)
    write(message,*)'cosp:', params%cosp
    call write_log(message)

  end subroutine glimmap_diag_stere

end module glimmer_map_types
