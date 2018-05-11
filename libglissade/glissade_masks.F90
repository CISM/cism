!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_masks.F90 - part of the Community Ice Sheet Model (CISM)  
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
!
! This module contains routines for computing various masks used by the Glissade 
! velocity solver.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_masks

    use glimmer_global, only: dp
    use glimmer_log
    use glimmer_physcon, only: rhoi, rhoo
    use glide_types
    use parallel

    implicit none

    ! All subroutines in this module are public

  contains

!****************************************************************************

  subroutine glissade_get_masks(nx,          ny,          &
                                thck,        topg,        &
                                eus,         thklim,      &
                                ice_mask,                 &
                                floating_mask,            &
                                ocean_mask,               &
                                land_mask,                &
                                grounding_line_mask,      &
                                active_ice_mask,          &
                                which_ho_calving_front,   &
                                calving_front_mask,       &
                                thck_calving_front,       &
                                marine_cliff_mask)
                                  
    !----------------------------------------------------------------
    ! Compute various masks for the Glissade dycore.
    !
    ! There are different ways to handle masks.
    ! The approach in Glide was to define an integer cell_mask array,
    !  in which each bit encodes information about whether a cell is ice-covered
    !  or ice-free, floating or grounded, etc.
    ! The approach here is to compute a separate 2D integer array for
    !  each kind of mask. This uses more memory but also is more transparent.
    !
    ! The basic masks used for Glissade dynamic calculations are as follows:
    !
    ! (1) ice_mask = 1 where ice is present (thck > thklim), else = 0
    ! (2) floating_mask = 1 if ice is present (thck > thklim) and floating, else = 0
    ! (3) ocean_mask = 1 if the topography is below sea level (topg < eus) and thk <= thklim, else = 0
    ! (4) land_mask = 1 if the topography is at or above sea level (topg >= eus), else = 0
    ! (5) grounding_line_mask = 1 if a cell is adjacent to the grounding line, else = 0
    ! (5) active_ice_mask = 1 for dynamically active cells, else = 0
    !     With the subgrid calving front scheme, cells that lie on the calving front and have
    !     thck < thck_calving_front are inactive. Otherwise, all cells with ice_mask = 1 are active.
    ! (6) calving_front_mask = 1 for floating cells that border at least one ocean cell, else = 0
    ! (7) marine_cliff_mask = 1 for grounded marine-based cells that border at least one ocean or
    !      inactive calving_front cell, else = 0
    !
    ! where thck = ice thickness
    !       thklim = threshold thickness for ice to be dynamically active
    !       topg = bed topography
    !       eus = eustatic sea level (= 0 by default)
    !       rhoi = ice density
    !       rhoo = ocean density
    !       thck_calving_front = effective thickness defined by adjacent cells not on the calving front
    !
    ! Notes:
    ! (1) thck, thklim, topg and eus can either have units of length (e.g., meters)
    !     or be dimensionless, as long as they are defined consistently.
    ! (2) ice_mask is always computed; the other masks are optional.
    ! (3) Thermal calculations may have a different threshold, thklim_temp
    !     (where generally thklim_temp < thklim).
    !     This mask can be computed by replacing thklim with thklim_temp in the subroutine call.
    ! (4) For some calculations it may be useful to call this subroutine with thklim = 0
    !     so as to identify all cells with nonzero ice thickness, not just dynamically active cells.
    !----------------------------------------------------------------
    
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
         nx,  ny                ! number of grid cells in each direction

    ! Default dimensions are meters, but this subroutine will work for
    ! any units as long as thck, topg, eus and thklim have the same units.

    real(dp), dimension(nx,ny), intent(in) ::  &
         thck,                & ! ice thickness (m)
         topg                   ! elevation of topography (m)

    real(dp), intent(in) :: &
         eus,                 & ! eustatic sea level (m), = 0. by default
         thklim                 ! minimum ice thickness for active cells (m)

    integer, dimension(nx,ny), intent(out) ::  &
         ice_mask               ! = 1 if thck > thklim, else = 0  

    integer, dimension(nx,ny), intent(out), optional ::  &
         floating_mask,       & ! = 1 if thck > thklim and ice is floating, else = 0
         ocean_mask,          & ! = 1 if topg is below sea level and thk <= thklim, else = 0
         land_mask,           & ! = 1 if topg is at or above sea level, else = 0
         grounding_line_mask, & ! = 1 if a cell is adjacent to the grounding line, else = 0
         active_ice_mask,     & ! = 1 if dynamically active, else = 0
         calving_front_mask,  & ! = 1 if ice is floating and borders at least one ocean cell, else = 0
         marine_cliff_mask      ! = 1 if ice is grounded and marine_based and borders at least one ocean
                                !  or inactive calving_front cell, else = 0  

    integer, intent(in), optional :: &
         which_ho_calving_front ! subgrid calving front option

    real(dp), dimension(nx,ny), intent(out), optional :: &
         thck_calving_front     ! effective ice thickness at the calving front

    !TODO - Make eps10 a model parameter?
    real(dp), parameter :: eps10 = 1.0d-10  ! small number

    !----------------------------------------------------------------
    ! Local arguments
    !----------------------------------------------------------------

    integer :: i, j, ii, jj

    !TODO - Make either of these optional output arguments?
    integer, dimension(nx,ny) :: &
         floating_interior_mask, &! = 1 if ice is floating and borders no ocean cells; else = 0
         grounded_mask            ! = 1 if ice is present and grounded, else = 0

    !----------------------------------------------------------------
    ! Compute masks in cells
    !----------------------------------------------------------------

    ice_mask(:,:) = 0

    do j = 1, ny
       do i = 1, nx

          if (thck(i,j) > thklim) then
             ice_mask(i,j) = 1
          else
             ice_mask(i,j) = 0
          endif

          if (present(ocean_mask)) then
             if (topg(i,j) < eus .and. ice_mask(i,j) == 0) then
                ocean_mask(i,j) = 1
             else
                ocean_mask(i,j) = 0
             endif
          endif

          if (present(floating_mask)) then
             if (topg(i,j) - eus < (-rhoi/rhoo)*thck(i,j) .and. ice_mask(i,j) == 1) then
                floating_mask(i,j) = 1
             else
                floating_mask(i,j) = 0
             endif
          endif

          if (present(land_mask)) then
             if (topg(i,j) >= eus) then
                land_mask(i,j) = 1
             else
                land_mask(i,j) = 0
             endif
          endif

          ! Note: active_ice_mask will be overwritten if the subgrid calving front scheme is used
          if (present(active_ice_mask)) then
             active_ice_mask(i,j) = ice_mask(i,j)
          endif

       enddo  ! i
    enddo  ! j

    ! halo updates
    ! Note: These are not strictly needed because the above loops include halo cells.
    !       However, they are included in case the user calls this subroutine without
    !        first updating thck in halo cells.

    call parallel_halo(ice_mask)
    if (present(floating_mask)) call parallel_halo(floating_mask)
    if (present(active_ice_mask)) call parallel_halo(active_ice_mask)

    ! Identify grounded cells; this mask is used in some calculations below
    if (present(floating_mask)) then
       where (ice_mask == 1 .and. floating_mask == 0)
          grounded_mask = 1
       elsewhere
          grounded_mask = 0
       endwhere
    endif

    ! Optionally, compute grounding line mask using grounded_mask, floating_mask and ocean_mask

    if (present(grounding_line_mask)) then

       if (.not.present(floating_mask) .or. .not.present(ocean_mask)) then
          call write_log('Need floating_mask and ocean_mask to compute grounding_line_mask', GM_FATAL)
       endif

       grounding_line_mask(:,:) = 0

       do j = 2, ny-1
          do i = 2, nx-1

             if (grounded_mask(i,j) == 1) then
                ! check whether one or more neighbors is a floating or ocean cell
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (floating_mask(ii,jj) == 1 .or. ocean_mask(ii,jj) == 1) then
                         grounding_line_mask(i,j) = 1
                      endif
                   enddo
                enddo
             elseif (floating_mask(i,j) == 1) then
                ! check whether one or more neighbors is a grounded cell
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (grounded_mask(ii,jj) == 1) then
                         grounding_line_mask(i,j) = 1
                      endif
                   enddo
                enddo
             endif   ! grounded_mask or floating_mask

          enddo   ! i
       enddo   ! j

       call parallel_halo(grounding_line_mask)

    endif   ! present(grounding_line_mask)

    ! Note: Halo calls are not included for the ocean and land masks.
    !       Halo values will still be correct, provided that topg is correct in halo cells.
    !       The reason not to include these calls is that for outflow global BCs,
    !        we may not want to set ocean_mask and land_mask = 0 in halo cells (as would be
    !        done automatically for outflow BCs). Instead, we want to compute ocean_mask
    !        and land_mask in halo cells based on topg (which for outflow BCs is extrapolated
    !        to halo cells from adjacent physical cells). 
    !       In particular, setting ocean_mask = 1 in the global halo ensures that calving_front
    !        cells are treated correctly just inside the global halo.

!!    if (present(ocean_mask)) call parallel_halo(ocean_mask)
!!    if (present(land_mask)) call parallel_halo(land_mask)

    ! Optionally, compute the calving_front mask and effective calving_front thickness

    if (present(calving_front_mask)) then

       if (.not.present(which_ho_calving_front)) then
          call write_log('Need which_ho_calving_front to compute calving_front_mask', GM_FATAL)
       endif

       if (which_ho_calving_front == HO_CALVING_FRONT_SUBGRID) then

          if (.not.present(floating_mask) .or. .not.present(ocean_mask)) then
             call write_log('Need floating_mask and ocean_mask to compute calving_front_mask', GM_FATAL)
          endif

          calving_front_mask(:,:) = 0
          floating_interior_mask(:,:) = 0

          ! Identify calving front cells (floating cells that border ice-free ocean)
          ! and floating interior cells (floating cells not at the calving front).
          do j = 2, ny-1
             do i = 2, nx-1
                if (floating_mask(i,j) == 1) then
                   if (ocean_mask(i-1,j) == 1 .or. ocean_mask(i+1,j) == 1 .or. &
                        ocean_mask(i,j-1) == 1 .or. ocean_mask(i,j+1) == 1) then
                      calving_front_mask(i,j) = 1
                   else
                      floating_interior_mask(i,j) = 1
                   endif
                endif
             enddo
          enddo

          call parallel_halo(calving_front_mask)
          call parallel_halo(floating_interior_mask)

          if (present(thck_calving_front)) then

             ! Compute an effective thickness in calving-front cells.
             ! This is set to the minimum nonzero thickness in a marine-based neighbor (either floating
             ! or grounded) that is not at the calving front.
             thck_calving_front(:,:) = 0.0d0

             do j = 2, ny-1
             do i = 2, nx-1
                if (calving_front_mask(i,j) == 1) then

                   ! loop over edge neighbors (floating interior or marine-based grounded)
                   !  and choose the minimum nonzero thickness
                   do jj = j-1, j+1
                      do ii = i-1, i+1
                         if ((ii == i .or. jj == j) .and. (ii /= i .or. jj /= j)) then  ! edge neighbors
                            if ( (grounded_mask(ii,jj) == 1 .and. topg(ii,jj) < eus) .or. &
                                 floating_interior_mask(ii,jj) == 1) then
                               if (thck_calving_front(i,j) > 0.0d0) then
                                  thck_calving_front(i,j) = min(thck_calving_front(i,j), thck(ii,jj))
                               else
                                  thck_calving_front(i,j) = thck(ii,jj)
                               endif
                            endif
                         endif
                      enddo   ! ii
                   enddo   !! jj

                   ! loop over corner neighbors if necessary
                   ! This will set thck_calving_front for cells at corners of the calving front,
                   !  with no interior edge neighbors
                   if (thck_calving_front(i,j) == 0.0d0) then
                      do jj = j-1, j+1
                         do ii = i-1, i+1
                            if (abs(ii-i) == 1 .and. abs(jj-j) == 1) then ! corner neighbors
                               if ( (grounded_mask(ii,jj) == 1 .and. topg(ii,jj) < eus) .or. &
                                    floating_interior_mask(ii,jj) == 1) then
                                  if (thck_calving_front(i,j) > 0.0d0) then
                                     thck_calving_front(i,j) = min(thck_calving_front(i,j), thck(ii,jj))
                                  else
                                     thck_calving_front(i,j) = thck(ii,jj)
                                  endif
                               endif
                            endif
                         enddo   ! ii
                      enddo   ! jj
                   endif   ! thck_calving_front = 0

                endif   ! calving front cell
             enddo   ! i
             enddo   ! j

             call parallel_halo(thck_calving_front)

          endif   ! present(thck_calving_front)

          ! Optionally, update the active_ice_mask so that calving_front cells with thck < thck_calving_front are inactive,
          ! but those with thck >= thck_calving_front are active.

          if (present(active_ice_mask)) then

             if (.not.present(thck_calving_front)) then
                call write_log  &
                     ('Must pass thck_calving_front to compute active_ice_mask with calving-front option', GM_FATAL)
             endif

             ! reset active_ice_mask
             active_ice_mask(:,:) = 0

             ! Mark ice-filled cells as active.
             ! Calving-front cells, however, are inactive, unless they have thck >= thck_calving front or
             !  are adjacent to grounded cells.

             do j = 2, ny-1
                do i = 2, nx-1
                   if (ice_mask(i,j) == 1) then
                      if (calving_front_mask(i,j) == 0) then
                         active_ice_mask(i,j) = 1
                      elseif (calving_front_mask(i,j) == 1) then
                         !WHL - There is a possible rounding issue here, if two adjacent cells are both being restored
                         !      (via inversion for bmlt_float) to the same thickness. For this reason, let the
                         !      cell be active if thck is very close to thck_calving front, but slightly less.

                         if (thck_calving_front(i,j) > 0.0d0 .and. &
                             thck(i,j)*(1.0d0 + eps10) >= thck_calving_front(i,j)) then
                            active_ice_mask(i,j) = 1
                         elseif (grounded_mask(i-1,j) == 1 .or. grounded_mask(i+1,j) == 1 .or. &
                                 grounded_mask(i,j-1) == 1 .or. grounded_mask(i,j+1) == 1 .or. &
                                 grounded_mask(i-1,j+1) == 1 .or. grounded_mask(i+1,j+1) == 1 .or. &
                                 grounded_mask(i-1,j-1) == 1 .or. grounded_mask(i+1,j-1) == 1) then
                            active_ice_mask(i,j) = 1
                         endif
                      endif
                   endif  ! ice_mask
                enddo
             enddo

             call parallel_halo(active_ice_mask)

          endif   ! active_ice_mask is present
 
       else   ! no subgrid calving front

          calving_front_mask(:,:) = 0
          if (present(thck_calving_front)) thck_calving_front(:,:) = 0.0d0

          ! Note: active_ice_mask, if present, was set above and need not be reset

       endif  ! which_ho_calving_front

    endif   ! calving_front_mask and which_ho_calving_front are present

    ! Optionally, compute the marine_cliff mask

    if (present(marine_cliff_mask)) then

       ! Make sure other required masks are present
       if (.not.present(floating_mask)   .or. .not.present(land_mask) .or. .not.present(active_ice_mask) ) then
          call write_log &
               ('Need floating, land and active_ice masks to compute marine_cliff_mask', GM_FATAL)
       endif

       marine_cliff_mask(:,:) = 0

       do j = 2, ny-1
          do i = 2, nx-1
             if (ice_mask(i,j) == 1 .and. land_mask(i,j) == 0 .and. floating_mask(i,j) == 0) then ! grounded marine-based ice
                if ( (land_mask(i-1,j) == 0 .and. active_ice_mask(i-1,j) == 0) .or. &  ! adjacent to inactive CF or ocean 
                     (land_mask(i+1,j) == 0 .and. active_ice_mask(i+1,j) == 0) .or. &
                     (land_mask(i,j-1) == 0 .and. active_ice_mask(i,j-1) == 0) .or. &
                     (land_mask(i,j+1) == 0 .and. active_ice_mask(i,j+1) == 0) ) then
                   marine_cliff_mask(i,j) = 1
                endif   ! marine cliff cell
             endif  ! grounded marine-based ice
          enddo  ! i
       enddo   ! j

       call parallel_halo(marine_cliff_mask)

    endif

  end subroutine glissade_get_masks

!****************************************************************************

  end module glissade_masks

!****************************************************************************

