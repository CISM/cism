! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_scales.f90 - part of the GLIMMER ice model       + 
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

!> this module holds scales for various fields
module glimmer_scales

  use glimmer_global, only : dp

  real(dp) :: scale2d_f1, scale2d_f2, scale2d_f3, scale2d_f4, scale2d_f5, scale2d_f6, scale2d_f7, scale2d_f8, scale2d_f9
  real(dp) :: scale3d_f1, scale3d_f2, scale3d_f3, scale3d_f4, scale3d_f5, scale3d_f6, scale3d_f7, scale3d_f8

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_SCALES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_SCALES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_SCALES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_SCALES
!MH!#endif

  !> calculate scale factors (can't have non-integer powers)
  subroutine glimmer_init_scales
    use physcon, only : scyr, gn
    use glimmer_paramets, only : thk0, tim0, vel0, vis0, len0, tau0, acc0
    implicit none

    scale2d_f1 = scyr * thk0 / tim0
    scale2d_f2 = scyr * vel0 * thk0
    scale2d_f3 = vel0 / (vis0 * len0)
    scale2d_f4 = vel0 * scyr * len0
    scale2d_f5 = scyr * vel0
    scale2d_f6 = scyr * vel0 * len0 / (thk0**2)
    scale2d_f7 = tau0
    scale2d_f8 = tau0 * len0 / (scyr * vel0)
    scale2d_f9 = scyr * acc0

    scale3d_f1 = scyr * vel0
    scale3d_f2 = vis0 * (vel0/len0)**(gn - 1)
    scale3d_f3 = scyr * thk0
    scale3d_f4 = vel0/(vis0*len0)
    scale3d_f5 = 1.0d0/scale3d_f2**(1.0/gn)
    scale3d_f6 = scale3d_f4**(1.0/gn)
    scale3d_f7 = scyr * thk0/tim0
    scale3d_f8 = vis0*scyr
  end subroutine glimmer_init_scales
end module glimmer_scales

