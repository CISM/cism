! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_mask.f90 - part of the GLIMMER ice model           + 
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
! GLIMMER is hosted on NeSCForge:
!
! http://forge.nesc.ac.uk/projects/glimmer/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glide_mask
    use glimmer_global, only : dp, sp, NaN
  !*FD masking ice thicknesses

!EIB! from glimmer-cism2/glimmer-cism-lanl
!EIB! needed by glide_lithot1d -> haven't merged yet
  integer, parameter :: glide_mask_ocean          = -2
  integer, parameter :: glide_mask_land           = -1
  integer, parameter :: glide_mask_boundary       = 0
  integer, parameter :: glide_mask_thin_ice       = 1
  integer, parameter :: glide_mask_interior       = 2
  integer, parameter :: glide_mask_shelf          = 4
  integer, parameter :: glide_mask_stream         = 8
  integer, parameter :: glide_mask_grounding_line = 16
  integer, parameter :: glide_mask_stream_margin  = 32
  integer, parameter :: glide_mask_land_margin    = 64
  integer, parameter :: glide_mask_shelf_front    = 128
  integer, parameter :: glide_mask_marine_edge    = 256

contains
  subroutine glide_set_mask(numerics, thck, topg, ewn, nsn, eus, mask, iarea, ivol)
    use glide_types
    use glimmer_physcon, only : rhoi, rhoo
    implicit none
    type(glide_numerics), intent(in) :: numerics !Numerical parameters structure
    real(dp), dimension(:,:), intent(in) :: thck !Ice thickness
    real(dp), dimension(:,:), intent(in) :: topg !Bedrock topography (not lower surface!)
    integer, intent(in) :: ewn, nsn !Grid size
    real(sp), intent(in) :: eus !Sea level
    integer, dimension(:,:), intent(inout) :: mask !Output mask
    real(dp), intent(inout), optional :: ivol, iarea !Area and volume of ice

    ! local variables
    integer ew,ns
    real(dp), parameter :: con = - rhoi / rhoo

    !Create an array to "fake" the boundaries of the mask so that boundary
    !finding can work even on the boundaries of the real mask.
    integer, dimension(0:ewn+1,0:nsn+1) :: maskWithBounds;

    mask = 0
    if (present(iarea)) then
        iarea = 0.
    end if

    if (present(ivol)) then
        ivol = 0.
    end if

    !Identify points with any ice
    where (thck > 0)
        mask = ior(mask, GLIDE_MASK_HAS_ICE)
    endwhere

    !Identify points where the ice is below the ice dynamics limit
    where (thck > 0 .and. thck < numerics%thklim)
        mask = ior(mask, GLIDE_MASK_THIN_ICE)
    endwhere

    !Identify points where the ice is floating or where there is open ocean
    where (topg - eus < con * thck)
        mask = ior(mask, GLIDE_MASK_OCEAN)
    elsewhere
        mask = ior(mask, GLIDE_MASK_LAND)
    endwhere

    if (present(iarea) .and. present(ivol)) then
        call get_area_vol(thck, numerics%dew, numerics%dns, iarea, ivol)
    end if

    maskWithBounds = 0
    maskWithBounds(1:ewn, 1:nsn) = MASK
    maskWithBounds(0,1:nsn) = mask(1,:)
    maskWithBounds(1:ewn,0) = mask(:,1)
    maskWithBounds(ewn+1,1:nsn) = mask(ewn,:)
    maskWithBounds(1:ewn,nsn+1) = mask(:,nsn)
    maskWithBounds(0,0) = mask(1,1)
    maskWithBounds(ewn+1,nsn+1) = mask(ewn,nsn)
    maskWithBounds(0,nsn+1) = mask(1,nsn)
    maskWithBounds(ewn+1,0) = mask(ewn,1)
    ! finding boundaries
    do ns=1,nsn
       do ew = 1,ewn
          !Find the grounding line
          if (GLIDE_IS_GROUND(MASK(ew,ns))) then
             if (GLIDE_IS_FLOAT(maskWithBounds(ew-1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew+1,ns)) .or. &
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns-1)) .or. & 
                  GLIDE_IS_FLOAT(maskWithBounds(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_GROUNDING_LINE)
             end if
          end if

          ! Ice margin
          ! *tb* A point is now masked even if it touches the ocean on one corner.
          if ( GLIDE_HAS_ICE(mask(ew, ns)) .and. &
              (GLIDE_NO_ICE(maskWithBounds(ew-1,ns))   .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew,ns-1))   .or. GLIDE_NO_ICE(maskWithBounds(ew,ns+1))   .or. &
               GLIDE_NO_ICE(maskWithBounds(ew-1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew-1,ns+1)) .or. &
               GLIDE_NO_ICE(maskWithBounds(ew+1,ns-1)) .or. GLIDE_NO_ICE(maskWithBounds(ew+1,ns+1)))) then
             MASK(ew,ns) = ior(MASK(ew,ns),GLIDE_MASK_MARGIN)
          end if

         !Mark domain boundaries
         if (ns == 1 .or. ns == nsn .or. ew == 1 .or. ew == ewn) then
            mask(ew, ns) = ior(mask(ew, ns), GLIDE_MASK_COMP_DOMAIN_BND)
         end if
       end do
    end do
  end subroutine glide_set_mask

!EIB! from glimmer-cism2/glimmer-cism-lanl
!EIB! needed by glide_lithot1d -> haven't merged yet
  logical elemental function is_ocean(mask)
    !*FD returns .true. if node is ocean
    implicit none
    integer, intent(in) :: mask 

    is_ocean = mask.eq.glide_mask_ocean
  end function is_ocean

  logical elemental function is_land(mask)
    !*FD returns .true. if node is land
    implicit none
    integer, intent(in) :: mask 

    is_land = mask.eq.glide_mask_land
  end function is_land

  logical elemental function has_ice(mask)
    !*FD returns .true. if node contains ice
    implicit none
    integer, intent(in) :: mask 

    has_ice = mask .gt. 0
  end function has_ice

  logical elemental function is_thin(mask)
    !*FD returns .true. if node is below dynamic limit
    implicit none
    integer, intent(in) :: mask

    is_thin = (iand(mask,glide_mask_thin_ice) .gt. 0 .and. mask.gt.0)
  end function is_thin

  logical elemental function is_float(mask)
    !*FD returns .true. if node is floating
    implicit none
    integer, intent(in) :: mask

    is_float = (iand(mask,glide_mask_shelf) .gt. 0 .and. mask.gt.0)
  end function is_float

  logical elemental function is_ground(mask)
    !*FD returns .true. if node is grounded ice
    implicit none
    integer, intent(in) :: mask

    is_ground = (iand(mask,glide_mask_interior) .gt. 0 .and. mask.gt.0)
  end function is_ground

  logical elemental function is_calving(mask)
    !*FD return .true. if node is at the shelf front
    implicit none
    integer, intent(in) :: mask

    is_calving = (iand(mask,glide_mask_shelf_front) .gt. 0 .and. mask.gt.0)
  end function is_calving

  logical elemental function is_marine_ice_edge(mask)
    !*FD return .true. if node is at edge of ice and topgraphy is
    !*FD below sea-level
    implicit none
    integer, intent(in) :: mask

    is_marine_ice_edge = (iand(mask,glide_mask_marine_edge) .gt. 0 .and. mask.gt.0)
  end function is_marine_ice_edge
  
  logical elemental function is_grounding_line(mask)
    !*FD returns .true. if node is grounding line
    implicit none
    integer, intent(in) :: mask 

    is_grounding_line = mask.eq.glide_mask_grounding_line
  end function is_grounding_line

!EIB! back to new stuff

  subroutine augment_kinbc_mask(mask, kinbcmask)
    !*FD Augments the Glide mask with the location of kinematic (dirichlet) boundary
    !*FD conditions.  These locations cannot be determined by the model a priori, and
    !*FD must be specified through a field in a NetCDF file.
    integer, dimension(:,:), target :: mask
    integer, dimension(:,:) :: kinbcmask

    integer, dimension(:,:), pointer :: maskp

    !Because the kinematic boundary conditions are specified on the staggered grid,
    !there may be a size mismatch here depending on whether we are computing a mask
    !for the staggered grid.
    if (size(mask, 1) /= size(kinbcmask, 1)) then
        maskp => mask(1:size(mask,1) - 1, 1:size(mask,2) - 1)
    else
        maskp => mask
    end if

    where (kinbcmask /= 0)
        maskp = ior(maskp, GLIDE_MASK_DIRICHLET_BC)
    endwhere
  end subroutine

  subroutine get_area_vol(thck, dew, dns, iarea, ivol)
    real(dp), dimension(:,:) :: thck
    real(dp) :: dew, dns
    real(dp) :: iarea, ivol

    integer :: i,j

    do i = 1, size(thck,1)
        do j = 2, size(thck,2)
            if (thck(i,j) > 0) then
                iarea = iarea + 1
                ivol = ivol + thck(i,j)
            end if
        end do
    end do

    iarea = iarea  * dew * dns
    ivol = ivol * dew * dns
    
  end subroutine get_area_vol
 
  subroutine calc_iareaf_iareag(dew, dns, iarea, mask, iareaf, iareag)
    
    implicit none
    real(dp), intent(in) :: dew, dns
    real(dp), intent(in) :: iarea
    real(dp), intent(out) :: iareaf, iareag
    integer, dimension(:,:), intent(in) :: mask 
    integer :: i,j
 
    iareaf = 0.0
    iareag = 0.0 

    do i = 1, size(mask,1)
        do j = 2, size(mask,2)
            if (GLIDE_IS_FLOAT(mask(i,j))) then
            iareaf = iareaf + 1
            else if(GLIDE_IS_GROUND_OR_GNDLINE(mask(i,j))) then
            iareag = iareag + 1
            end if
        end do
    end do

    iareaf = iareaf  * dew * dns
    iareag = iareag  * dew * dns
    
  end subroutine calc_iareaf_iareag

    subroutine glide_marine_margin_normal(thck, mask, marine_bc_normal)
        use glimmer_physcon, only:pi
        implicit none
        !*FD This subroutine derives from the given mask the normal to an ice shelf
        !*FD each point on the marine margin.
        real(dp), dimension(:,:), intent(in) :: thck
        integer, dimension(:,:), intent(in) :: mask
        real(dp), dimension(:,:), intent(out) :: marine_bc_normal

        integer :: i, j, dx, dy, k

        real(dp), dimension(size(thck,1), size(thck,2)) :: direction_x, direction_y
        
        real(dp), dimension(-1:1, -1:1) :: angle_lookup

                !direction_y =    -1       0       1        !direction_x = 
        angle_lookup(-1, :) = (/ 3*pi/4,   pi/2,   pi/4 /)  !-1
        angle_lookup( 0, :) = (/   pi,     0D0,  2*pi   /)  ! 0
        angle_lookup( 1, :) = (/ 5*pi/4, 3*pi/2, 7*pi/4 /)  ! 1
        call upwind_from_mask(mask, direction_x, direction_y)

        !Set up a thickness variable with "ghost cells" so that we don't go out
        !of bounds with the vectorized operation below
        !thckWithBounds(1:size(thck,1), 1:size(thck,2)) = thck
        !thckWithBounds(:,0) = thckWithBounds(:,1)
        !thckWithBounds(0,:) = thckWithBounds(1,:)
        !thckWithBounds(size(thck,1)+1,:) = thckWithBounds(size(thck,1),:)
        !thckWithBounds(:,size(thck,2)+1) = thckWithBounds(:,size(thck,2))
        do i = 1, size(mask, 1)
            do j = 1, size(mask, 2)
                if (GLIDE_IS_CALVING(mask(i,j))) then
                    dx = int(direction_x(i,j))
                    dy = int(direction_y(i,j))
                    if (dx == 0 .and. dy == 0) then
                        write(*,*)"A shelf front point has been identified at:"
                        write(*,*)"x = ",i
                        write(*,*)"y = ",j
                        write(*,*)"But neither x nor y derivatives have been marked as upwinded."
                        write(*,*)"This should never happen, if this error appears it is a bug"
                        write(*,*)"and should be reported."
                        write(*,*)"The mask around this point follows:"
                        write(*,*)"--------------------------"

                        !Write a header row with a * in the column corresponding to the center
                        do k = -4, 4
                            if (k==0) then
                                write(*,"(A)",advance="no")"           *"
                            else if (i+k > 0 .and. i+k <= size(mask,1)) then
                                write(*,"(A)",advance="no")"            "
                            end if
                        end do 
                        write(*,*)

                        do k=4, -4, -1
                            if (j+k > 0 .and. j+k <= size(mask, 2)) then
                                if (k == 0) then
                                    write(*,*) "*", mask(max(1,i-4):min(size(mask,1),i+4),j+k)
                                else
                                    write(*,*) " ", mask(max(1,i-4):min(size(mask,1),i+4),j+k)
                                end  if
                            end if
                        end do
                        write(*,*)"--------------------------"
                        write(*,*)"Have a nice day!"
                        !stop
                    end if
                    marine_bc_normal(i,j) = angle_lookup(dx, dy) 
                    !marine_bc_normal(i,j) = calc_normal_45deg(thckWithBounds(i-1:i+1,j-1:j+1))
                else
                    marine_bc_normal(i,j) = NaN
                end if
            end do
        end do
        
    end subroutine

    function calc_normal_45deg(thck3x3)
        use glimmer_physcon, only: pi
        
        !*FD Computes the angle of the normal vector, in radians, for the given
        !*FD 3x3 segment of ice geometry.
        !*FD The normal is given in increments of 45 degrees (no nicer
        !*FD interpolation is currently done)
        !*FD This is based on the Payne and Price GLAM code, if/when this is
        !*FD integrated into CISM it should probably be refactored to use this.
        real(dp), dimension(3,3) :: thck3x3

        real(dp) :: calc_normal_45deg
         
        real (kind = dp), dimension(3,3) :: mask, maskcorners
        real (kind = dp), dimension(3,3) :: thckmask
        real (kind = dp), dimension(3) :: testvect
        real (kind = dp) :: phi, deg2rad

        deg2rad = pi / 180.0d0
        loc_latbc = 0; phi = 0
        mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
        mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
        mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
        maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
        maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
        maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

        ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
        ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
        ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
        ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
        ! and 'crosshorizmainbc_lat'

        ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
        ! boundary shape

        where( thck3x3 .ne. 0.0d0 )
            thckmask = 0.0_dp
        elsewhere( thck3x3 .eq. 0.0d0 )
            thckmask = 1.0d0
        endwhere

        testvect = sum( thckmask * mask, 1 )

        !if( up .eq. 3 )then ! temporary code for debugging
        !  do i = 3,1,-1
        !  print *, 'thck = ', thck(:,i)
        !  end do
        !  print *, ' '
        !
        !  do i = 3,1,-1
        !      print *, 'thckmask = ', thckmask(:,i)
        !  end do
        !  print *, ' '
        !
        !  print *, 'testvect =  ', testvect
        !  print *, ' '
        !end if

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) .eq. 1.0d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect .eq. 360.0d0 ) )then
                if( sum( testvect ) .eq. 450.0d0 )then
                    phi = 45.0d0
                elseif( sum( testvect ) .eq. 630.0d0 )then
                    phi = 315.0d0
                else
                    phi = 0.0d0
                end if
            elseif( all( testvect .ne. 360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect .ne. 0.0d0 )
            end if
        end if

        calc_normal_45deg = deg2rad * phi
        
        !Tim's Note: This appears to actually compute 0 at 6 O'clock according
        !to Glimmer's coordinate system.  90 deg. is still 3 O'clock.
        !I'm going to correct for this here rather than dig through the code
        !above
        !(TODO: correct it in the code above!)
        calc_normal_45deg = pi - calc_normal_45deg 
        if (calc_normal_45deg < 0) calc_normal_45deg = calc_normal_45deg + 2*pi

    end function

    !Fills a field of differencing directions suitable to give a field
    !derivative routine.  Uses centered differencing everywhere except for the
    !marine ice margin, where upwinding and downwinding is used to avoid
    !differencing across the boundary.
    subroutine upwind_from_mask(mask, direction_x, direction_y)
        integer, dimension(:,:), intent(in) :: mask
        double precision, dimension(:,:), intent(out) :: direction_x, direction_y

        integer :: i,j

        direction_x = 0
        direction_y = 0

        !Detect locations of the marine margin
        do i = 1, size(mask,1)
            do j = 1, size(mask,2)
                if (GLIDE_IS_CALVING(mask(i,j))) then
                    !Detect whether we need to upwind or downwind in the Y
                    !direction
                    if (i > 1) then
                        if (.not. GLIDE_HAS_ICE(mask(i-1,j))) then
                            direction_x(i,j) = 1
                        end if
                    end if

                    if (i < size(mask, 1)) then
                        if (.not. GLIDE_HAS_ICE(mask(i+1,j))) then
                            direction_x(i,j) = -1
                        end if
                    end if

                    !Detect whether we need to upwind or downwind in the X
                    !direction
                    if (j > 1) then
                        if (.not. GLIDE_HAS_ICE(mask(i,j-1))) then
                            direction_y(i,j) = 1
                        end if
                    end if
                    
                    if (j < size(mask, 2)) then
                        if (.not. GLIDE_HAS_ICE(mask(i,j+1))) then
                            direction_y(i,j) = -1
                        end if
                    end if

                    !If we are at a point that is "interior" to two other boundary points, 
                    !such as the lower right of:
                    !o b i
                    !b b i
                    !(o = ocean, b = boundary, i = interior), then we will not detect the need
                    !to upwind or downwind.  However, we still should for consistency with other
                    !mask points (in some cases, not doing so can lead to a singular calculation
                    !at the marine ice front)
                    !
                    !We can think of this operation as avoiding calving points where there is 
                    !a non-calving point to upwind into.
                    !
                    !TODO: We need a better way to detect interior points.  Right now I am just using
                    !points that are floating, and that works, but this doesn't work for two reasons:
                    !1. Boundary points are also floating
                    !2. Could fail for a very thin ice shelf
                    if (int(direction_x(i,j)) == 0 .and. int(direction_y(i,j)) == 0 .and. &
                        i > 1 .and. j > 1 .and. i < size(mask, 1) .and. j < size(mask, 2)) then
                        if (.not. GLIDE_HAS_ICE(mask(i-1, j-1))) then
                            direction_x(i,j) = 1
                            direction_y(i,j) = 1
                        else if (.not. GLIDE_HAS_ICE(mask(i-1, j+1))) then
                            direction_x(i,j) = 1
                            direction_y(i,j) = -1
                        else if (.not. GLIDE_HAS_ICE(mask(i+1, j-1))) then
                            direction_x(i,j) = -1
                            direction_y(i,j) = 1
                        else if (.not. GLIDE_HAS_ICE(mask(i+1, j+1))) then
                            direction_x(i,j) = -1
                            direction_y(i,j) = -1
                        end if
                    end if
                end if
            end do
        end do
    end subroutine

end module glide_mask
