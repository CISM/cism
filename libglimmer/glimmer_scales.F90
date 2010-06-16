! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_scales.f90 - part of the Glimmer-CISM ice model  + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

!> this module holds scales for various fields
module glimmer_scales

  use glimmer_global, only : dp

  real(dp) :: scale2d_f1, scale2d_f2, scale2d_f3, scale2d_f4, scale2d_f5, scale2d_f6, scale2d_f7, scale2d_f8, scale2d_f9
  real(dp) :: scale3d_f1, scale3d_f2, scale3d_f3, scale3d_f4, scale3d_f5, scale3d_f6, scale3d_f7, scale3d_f8

contains

  !> calculate scale factors (can't have non-integer powers)
  subroutine glimmer_init_scales
    use glimmer_physcon, only : scyr, gn
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

