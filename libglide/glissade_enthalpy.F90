!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   glissade_enthalpy.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)
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
!
! This module computes temperature diffusion, strain heating, and local
!  melting and refreezing in each ice column using enthalpy as the
!  primary state variable.
!
! It is under construction; just a stub for now
!

#include "glide_mask.inc"

module glissade_enthalpy

    use glimmer_global, only : dp
    use glide_types
    use glimmer_log

    implicit none

    private
    public :: glissade_init_enthalpy, glissade_enthalpy_driver

!--------------------------------------------------------------------------------

contains

!--------------------------------------------------------------------------------

  subroutine glissade_init_enthalpy (model)

    ! initialization for the enthalpy scheme

    type(glide_global_type),intent(inout) :: model      ! ice model parameters

    call write_log('Glissade enthalpy scheme is still under construction', GM_FATAL)

  end subroutine glissade_init_enthalpy

!--------------------------------------------------------------------------------

  subroutine glissade_enthalpy_driver (model)

    ! driver for the enthalpy scheme

    type(glide_global_type),intent(inout) :: model      ! ice model parameters

    call write_log('Glissade enthalpy scheme is still under construction', GM_FATAL)

  end subroutine glissade_enthalpy_driver

!--------------------------------------------------------------------------------

end module glissade_enthalpy

!--------------------------------------------------------------------------------
