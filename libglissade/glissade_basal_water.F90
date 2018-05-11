!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_water.F90 - part of the Community Ice Sheet Model (CISM)  
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

!TODO - Support Jesse's water-routing code (or something similar) in parallel?
!       Currently supported only for serial Glide runs, in module glide_bwater.F90

module glissade_basal_water

   use glimmer_global, only: dp
   use glide_types

   implicit none

   private
   public :: glissade_basal_water_init, glissade_calcbwat

contains

  subroutine glissade_basal_water_init(model)

    use glimmer_paramets, only: thk0

    type(glide_global_type) :: model

    select case (model%options%which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate

    case(HO_BWAT_CONSTANT)

       ! Set a constant water thickness where ice is present
       where (model%geometry%thck > model%numerics%thklim)
          model%temper%bwat(:,:) = model%basal_physics%const_bwat / thk0
       elsewhere
          model%temper%bwat(:,:) = 0.0d0
       endwhere

    case default

       ! currently nothing to do for other options

    end select

  end subroutine glissade_basal_water_init


  subroutine glissade_calcbwat(which_ho_bwat,    &
                               basal_physics,    &
                               dt,               &
                               thck,             &
                               thklim,           &
                               bmlt,             &
                               bwat)

    ! Driver for updating basal water
    ! Note: This subroutine assumes SI units.
    ! Currently, only a few simple options are supported.

    use glimmer_physcon, only: rhow, scyr
    use glide_types

    !WHL - debug
    use parallel

    integer, intent(in) :: &
         which_ho_bwat     !> basal water options

    type(glide_basal_physics), intent(in)   :: basal_physics      ! basal physics object

    real(dp), intent(in) :: &
         dt,             & !> time step (s) 
         thklim            !> threshold for dynamically active ice (m)

    real(dp), dimension(:,:), intent(in) ::  &
         thck,           & !> ice thickness (m)
         bmlt              !> basal melt rate (m/s of ice)

    real(dp), dimension(:,:), intent(inout) ::  &
         bwat              !> basal water depth (m of water)

    ! local variables

    integer :: nx, ny    ! horizontal grid dimensions
    integer :: i, j

    real(dp) ::  &
         dbwat_dt        ! rate of change of bwat (m/s of water)

    select case (which_ho_bwat)

    ! HO_BWAT_NONE:         basal water depth = 0
    ! HO_BWAT_CONSTANT:     basal water depth = prescribed constant
    ! HO_BWAT_LOCAL_TILL:   local basal till model with prescribed drainage rate

    case(HO_BWAT_NONE)

       bwat(:,:) = 0.0d0

    case(HO_BWAT_CONSTANT)

       ! Use a constant water thickness where ice is present, to force Tbed = Tpmp
       where (thck > thklim)
          bwat(:,:) = basal_physics%const_bwat
       elsewhere
          bwat(:,:) = 0.0d0
       endwhere

     case(HO_BWAT_LOCAL_TILL)

        nx = size(bwat,1)
        ny = size(bwat,2)

        do j = 1, ny
           do i = 1, nx

              ! compute new bwat, given source (bmlt) and sink (drainage)
              ! Note: bmlt > 0 for ice melting. Freeze-on will reduce bwat.
              dbwat_dt = bmlt(i,j)*rhoi/rhow - basal_physics%c_drainage/scyr  ! convert c_drainage from m/yr to m/s
              bwat(i,j) = bwat(i,j) + dbwat_dt*dt

              ! limit to the range [0, bwat_till_max]
              bwat(i,j) = min(bwat(i,j), basal_physics%bwat_till_max)
              bwat(i,j) = max(bwat(i,j), 0.0d0)

           enddo
        enddo

    end select

  end subroutine glissade_calcbwat

end module glissade_basal_water
