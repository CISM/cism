!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grounding_line.F90 - part of the Community Ice Sheet Model (CISM)  
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
! This module contains routines for computing grounding-line fields and diagnostics
! for the Glissade solver.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  module glissade_grounding_line

    use glimmer_global, only: dp
    use glimmer_physcon, only: rhoi, rhoo
    use glide_types  ! grounding line options
    use glimmer_log
    use parallel

    implicit none

    ! All subroutines in this module are public

  contains

!****************************************************************************

  subroutine glissade_grounded_fraction(nx,            ny,                      &
                                        ice_mask,                               &
                                        floating_mask, land_mask,               &
                                        whichground,   f_ground)

    !----------------------------------------------------------------
    ! Compute fraction of ice that is grounded.
    ! This fraction is computed at vertices based on the thickness and
    !  topography of the four neighboring cell centers.
    !
    ! There are three options for computing the grounded fraction, based on the value of whichground:
    ! (0) HO_GROUND_NO_GLP: f_ground = 1 for vertices with grounded and/or land-based neighbor cells
    !                       f_ground = 0 for vertices with floating neighbors only
    ! (1) HO_GROUND_GLP: 0 <= f_ground <= 1 based on grounding-line parameterization
    !        A flotation function is interpolated over the bounding box of each vertex
    !        and analytically integrated to compute the grounded and floating fractions.
    ! (2) HO_GROUND_ALL: f_ground = 1 for all vertices with ice-covered neighbor cells
    !
    ! NOTE: Option 1 = HO_GROUND_GLP is not suppported for this release.
    !
    !----------------------------------------------------------------
    
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny                ! number of grid cells in each direction

    integer, dimension(nx,ny), intent(in) ::   &
       ice_mask,            & ! = 1 if ice is present (thck > thklim), else = 0
       floating_mask,       & ! = 1 if ice is present (thck > thklim) and floating, else = 0
       land_mask              ! = 1 if topg is at or above sea level

    ! see comments above for more information about these options
    integer, intent(in) ::     &
       whichground            ! option for computing f_ground

    real(dp), dimension(nx-1,ny-1), intent(out) ::  &
       f_ground               ! grounded ice fraction at vertex, 0 <= f_ground <= 1

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
           
    integer :: i, j

    integer, dimension(nx-1,ny-1) ::   &
       vmask                     ! = 1 for vertices neighboring at least one cell where ice is present, else = 0

    logical, dimension(nx,ny) :: &
       cground             ! true if a cell is land and/or has grounded ice, else = false

    !----------------------------------------------------------------
    ! Compute ice mask at vertices (= 1 if any surrounding cells have ice or are land)
    !----------------------------------------------------------------

    do j = 1, ny-1
       do i = 1, nx-1
          if (ice_mask(i,j+1)==1  .or. ice_mask(i+1,j+1)==1  .or.   &
              ice_mask(i,j)  ==1  .or. ice_mask(i+1,j)  ==1  .or.   &
              land_mask(i,j+1)==1 .or. land_mask(i+1,j+1)==1 .or.   &
              land_mask(i,j)  ==1 .or. land_mask(i+1,j+1)==1) then
             vmask(i,j) = 1
          else
             vmask(i,j) = 0
          endif
       enddo
    enddo

    ! initialize f_ground
    f_ground(:,:) = 0.0d0

    ! Compute f_ground according to the value of whichground

    select case(whichground)

    case(HO_GROUND_NO_GLP)   ! default: no grounding-line parameterization
                             ! f_ground = 1 at a vertex if any neighbor cell is land or has grounded ice

       ! compute a mask that is true for cells that are land and/or have grounded ice
       do j = 1, ny
          do i = 1, nx
             if ((ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .or. land_mask(i,j) == 1) then
                cground(i,j) = .true.
             else
                cground(i,j) = .false.
             endif
          enddo
       enddo

       ! vertices are grounded if any neighbor cell is land and/or has grounded ice, else are floating

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j) == 1) then
                if (cground(i,j+1) .or. cground(i+1,j+1) .or. cground(i,j) .or. cground(i+1,j)) then
                   f_ground(i,j) = 1.d0
                else
                   f_ground(i,j) = 0.d0
                endif
             endif
           enddo
        enddo

    case(HO_GROUND_ALL)

       ! all vertices with ice-covered or land-based neighbors are assumed grounded, regardless of thck and topg

       do j = 1, ny-1
          do i = 1, nx-1
             if (vmask(i,j) == 1) then
                f_ground(i,j) = 1.d0
             endif
          enddo
       enddo

    case(HO_GROUND_GLP)      ! grounding-line parameterization

       call write_log('The GLP option for which_ho_ground is not supported for this release', GM_FATAL)

    end select

  end subroutine glissade_grounded_fraction

!=======================================================================

  subroutine glissade_grounding_line_flux(nx,                       ny,            &
                                          dx,                       dy,            &
                                          sigma,                                   &
                                          thck,                                    &
                                          uvel,                     vvel,          &
                                          ice_mask,                 floating_mask, &
                                          ocean_mask,                              &
                                          gl_flux_east,             gl_flux_north, &
                                          gl_flux                                   )

    ! Computes northward and eastward land ice fluxes at grounding lines,
    !  and a cell-based grounding-line flux field.
    ! Note: Since the GL thicknesses are approximated, the GL fluxes will not exactly 
    !        match the fluxes computed by the transport scheme.
    !       Also, the GL fluxes do not include thinning/calving of grounded marine cliffs.

    use parallel, only: nhalo
    use glimmer_paramets, only: thk0, vel0, len0

    implicit none

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::                     &
        nx, ny                                   !> horizontal grid dimensions

    real(dp), intent(in) ::                    &
        dx, dy                                   !> horizontal grid spacing

    real(dp), dimension(:), intent(in) ::      &
        sigma                                    !> vertical sigma coordinate

    real(dp), dimension(nx,ny), intent(in) ::  &
        thck                                     !> ice thickness

    real(dp), dimension(:,:,:), intent(in) ::  &
        uvel, vvel                               !> ice velocity in x and y directions

    integer, dimension(nx,ny), intent(in) ::  &
        ice_mask,                              & !> = 1 where ice is present, else = 0
        floating_mask,                         & !> = 1 where ice is present and floating, else = 0
        ocean_mask                               !> = 1 for ice-free ocean, else = 0

    ! Note: gl_flux_east and gl_flux_north are directional 
    !       (positive for eastward/northward, negative for westward/southward)
    !       gl_flux is a cell-based quantity based on flux magnitudes on each edge
    !       (so gl_flux >= 0)

    real(dp), dimension(:,:), intent(out) ::   &
        gl_flux_east,                          & !> grounding line flux on east edges
        gl_flux_north,                         & !> grounding line flux on north edges
        gl_flux                                  !> grounding line flux per grid cell


    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer  :: i,j,k                                     !> local cell indices
    integer  :: upn                                       !> vertical grid dimension
    real(dp), dimension(:), allocatable :: uavg, vavg     !> local horizontal velocity averages
    real(dp) :: thck_gl                                   !> GL thickness derived from topg_gl

    upn = size(sigma)

    allocate(uavg(upn), vavg(upn))

    ! Initialize
    gl_flux_east(:,:)  = 0.d0
    gl_flux_north(:,:) = 0.d0
    gl_flux(:,:)       = 0.d0

    ! Compute grounding line fluxes on east and north edges.
    ! Look for edges with a grounded cell on one side and a floating cell on the other.

    do j = nhalo+1, ny-nhalo
        do i = nhalo+1, nx-nhalo

            ! check east edge
           if ( (   (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .and.   &  ! (i,j) grounded
                (ocean_mask(i+1,j) == 1 .or.  floating_mask(i+1,j) == 1) )     &  ! (i+1,j) floating or ocean
                                        .or.                                   &
                ( (ice_mask(i+1,j) == 1 .and. floating_mask(i+1,j) == 0) .and. &  ! (i+1,j) grounded
                  (ocean_mask(i,j) == 1  .or. floating_mask(i,j) == 1) ) ) then   ! (i,j) floating or ocean

                uavg(:) = (uvel(:,i,j) + uvel(:,i,j-1)) / 2.d0
                if (ice_mask(i,j) == 1 .and. ice_mask(i+1,j) == 1) then
                   ! set GL thickness to the average thickness of the two cells
                   thck_gl = (thck(i,j) + thck(i+1,j)) / 2.d0
                else
                   ! set GL thickness to the thickness of the ice-filled cell
                   thck_gl = max(thck(i,j), thck(i+1,j))
                endif

                do k = 1, upn-1
                    gl_flux_east(i,j) = gl_flux_east(i,j) &
                                        + thck_gl * (sigma(k+1) - sigma(k)) * (uavg(k) + uavg(k+1))/2.d0
                enddo
            endif

            ! check north edge
           if ( (   (ice_mask(i,j) == 1 .and. floating_mask(i,j) == 0) .and.   &  ! (i,j) grounded
                (ocean_mask(i,j+1) == 1 .or.  floating_mask(i,j+1) == 1) )     &  ! (i,j+1) floating or ocean
                                        .or.                                   &
                ( (ice_mask(i,j+1) == 1 .and. floating_mask(i,j+1) == 0) .and. &  ! (i,j+1) grounded
                  (ocean_mask(i,j) == 1  .or. floating_mask(i,j) == 1) ) ) then   ! (i,j) floating or ocean

                vavg(:) = (vvel(:,i-1,j) + vvel(:,i,j)) / 2.d0
                if (ice_mask(i,j) == 1 .and. ice_mask(i,j+1) == 1) then
                   ! set GL thickness to the average thickness of the two cells
                   thck_gl = (thck(i,j) + thck(i,j+1)) / 2.d0
                else
                   ! set GL thickness to the thickness of the ice-filled cell
                   thck_gl = max(thck(i,j), thck(i,j+1))
                endif

                do k = 1, upn-1
                    gl_flux_north(i,j) = gl_flux_north(i,j) &
                                        + thck_gl * (sigma(k+1) - sigma(k)) * (vavg(k) + vavg(k+1))/2.d0
                enddo
             endif

        enddo   ! i
    enddo   ! j

    ! Compute mass flux through grounding line in each cell.
    ! Only a grounded cell can lose mass. We need to check the direction of the fluxes.

    do j = nhalo+1,ny-nhalo
        do i = nhalo+1,nx-nhalo

            ! Check the sign for east-west flow and assign the flux accordingly
            if (gl_flux_east(i,j) < 0.d0) then
                ! The ice is flowing westward and the flux belongs to the right adjacent cell
                gl_flux(i+1,j) = gl_flux(i+1,j) - gl_flux_east(i,j)
            else
                ! The ice is flowing eastward and the flux belongs to this cell
                gl_flux(i,j) = gl_flux(i,j) + gl_flux_east(i,j)
            endif

            ! Check the sign for north-south flow and assign the flux accordingly
            if (gl_flux_north(i,j) < 0.d0) then
                ! The ice is flowing southward and the flux belongs to the top adjacent cell
                gl_flux(i,j+1) = gl_flux(i,j+1) - gl_flux_north(i,j)
            else
                ! The ice is flowing northward and the flux belongs to this cell
                gl_flux(i,j) = gl_flux(i,j) + gl_flux_north(i,j)
            endif

        enddo   ! i
    enddo   ! j

    ! Convert from model units to kg/m/s
    gl_flux_east  = gl_flux_east  * rhoi*thk0*vel0
    gl_flux_north = gl_flux_north * rhoi*thk0*vel0
    gl_flux       = gl_flux       * rhoi*thk0*vel0

    deallocate(uavg, vavg)

  end subroutine glissade_grounding_line_flux

!****************************************************************************

end module glissade_grounding_line

!****************************************************************************

