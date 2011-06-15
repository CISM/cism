!Helper module containing routines to move between staggered and
!unstaggered grids
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_nan.inc"
#include "glide_mask.inc"

module glide_grids
    use glimmer_global, only : dp, NaN
    implicit none

contains
 
  !> A special staggering algorithm that is meant to conserve mass when operating on thickness fields
  !! This incorporates Ann LeBroque's nunatak fix and the calving front fix.
  !!
  !! \param ipvr The input thickness field (ewn x nsn)
  !! \param opvr The output (staggered) thickness field (ewn - 1 x nsn - 1)
  !! \param ewn
  !! \param nsn
  !! \param usrf   Surface elevation field, non-staggered.
  !! \param thklim Minimum thickness to enable ice dynamics
  !! \param mask   Geometry mask field (used for determining the location of shelf fronts)
  !<
  subroutine stagthickness(ipvr,opvr,ewn,nsn,usrf,thklim,mask)
    use glimmer_paramets, only: thk0
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    
    real(dp), intent(in), dimension(:,:) :: usrf
    real(dp) :: thklim
    integer, intent(in), dimension(:,:) :: mask
    
    integer :: ewn,nsn,ew,ns,n
    real(dp) :: tot

        do ns = 1,nsn-1
            do ew = 1,ewn-1

                !If any of our staggering points are shelf front, ignore zeros when staggering
                !if (any(GLIDE_IS_CALVING(mask(ew:ew+1, ns:ns+1)))) then

                !Use the "only nonzero thickness" staggering criterion for ALL marginal ice. For
                ! reasons that are not entirely clear, this corrects an error whereby the land ice 
                ! margin is defined incorrectly as existing one grid cell too far inland from where 
                ! it should be.  
                if (any(GLIDE_HAS_ICE(mask(ew:ew+1,ns:ns+1)))) then
                    n = 0
                    tot = 0
    
                    if (abs(ipvr(ew,ns)) > 1e-10) then
                        tot = tot + ipvr(ew,ns)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew+1,ns)) > 1e-10) then
                        tot = tot + ipvr(ew+1,ns)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew,ns+1)) > 1e-10) then
                        tot = tot + ipvr(ew,ns+1)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew+1,ns+1)) > 1e-10) then
                        tot = tot + ipvr(ew+1,ns+1)
                        n   = n   + 1
                    end if
                    if (n > 0) then
                        opvr(ew,ns) = tot/n
                    else
                        opvr(ew,ns) = 0
                    end if
                !The following cases relate to Anne LeBroque's fix for nunataks
                !ew,ns cell is ice free:
                else if (ipvr(ew,ns) <= thklim/thk0 .and. &
                   ((usrf(ew,ns) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >= thklim/thk0) &
                    .or. (usrf(ew,ns) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0

                !ew+1,ns cell is ice free:
                else if (ipvr(ew+1,ns) <= thklim/thk0 .and. &
                    ((usrf(ew+1,ns) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
                    .or. (usrf(ew+1,ns) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0
    
                !ew,ns+1 cell is ice free:
                else if (ipvr(ew,ns+1) <= thklim/thk0 .and. &
                    ((usrf(ew,ns+1) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim/thk0) &
                    .or. (usrf(ew,ns+1) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim/thk0))) then
                        opvr(ew,ns) = 0.0
    
                !ew+1,ns+1 cell is ice free:
                else if (ipvr(ew+1,ns+1) <= thklim/thk0 .and. &
                    ((usrf(ew+1,ns+1) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >=thklim/thk0) &
                    .or. (usrf(ew+1,ns+1) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >=thklim/thk0))) then
                        opvr(ew,ns) = 0.0
                
!                !Standard Staggering   !! Not needed if only-nonzero-thickness staggering scheme is used
!                else
!                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
!                               ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
                end if
  
        end do
    end do

  end subroutine stagthickness

 
  !> Moves a variable from the ice grid to the velocity grid by averaging onto the centroids
  !! \param ipvr Input variable (on the ice grid)
  !! \param opvr Output variable (on the velocity grid)
  !! \param ewn
  !! \param nsn
  !<
  subroutine stagvarb(ipvr,opvr,ewn,nsn)
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn

        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0
  end subroutine stagvarb

  subroutine stagvarb_3d(ipvr, opvr, ewn, nsn, upn)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer :: k

    do k = 1, upn
        call stagvarb(ipvr(k,:,:), opvr(k,:,:), ewn, nsn)
    end do
  end subroutine stagvarb_3d

  subroutine stagvarb_mask(ipvr,opvr,ewn,nsn,geometry_mask)
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:)  :: ipvr
    
    integer, intent(in) :: ewn,nsn
    integer, intent(in), dimension(:,:) :: geometry_mask
    integer :: ew,ns,n
    real(dp) :: tot

        opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                                 ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0
        do ns = 1,nsn-1
            do ew = 1,ewn-1

                !If any of our staggering points are shelf front, ignore zeros when staggering
                if (any(GLIDE_NO_ICE(geometry_mask(ew:ew+1, ns:ns+1)))) then
                    n = 0
                    tot = 0
    
                    if (GLIDE_HAS_ICE(geometry_mask(ew,ns))) then
                        tot = tot + ipvr(ew,ns)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew+1,ns))) then
                        tot = tot + ipvr(ew+1,ns)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew,ns+1))) then
                        tot = tot + ipvr(ew,ns+1)
                        n   = n   + 1
                    end if
                    if (GLIDE_HAS_ICE(geometry_mask(ew+1,ns+1))) then
                        tot = tot + ipvr(ew+1,ns+1)
                        n   = n   + 1
                    end if
                    if (n > 0) then
                        opvr(ew,ns) = tot/n
                    else
                        opvr(ew,ns) = 0
                    end if
                
                !Standard Staggering
                else
                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
                               ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0
                end if
  
        end do
    end do
  end subroutine stagvarb_mask

  subroutine stagvarb_3d_mask(ipvr, opvr, ewn, nsn, upn, geometry_mask)
    real(dp), intent(in), dimension(:,:,:) :: ipvr
    real(dp), intent(out), dimension(:,:,:) :: opvr
    integer, intent(in) :: ewn, nsn, upn
    integer, intent(in), dimension(:,:) :: geometry_mask
    integer :: k

    do k = 1, upn
        call stagvarb_mask(ipvr(k,:,:), opvr(k,:,:), ewn, nsn, geometry_mask)
    end do
  end subroutine stagvarb_3d_mask




    !*FD Copies a staggered grid onto a nonstaggered grid.  This verion
    !*FD assumes periodic boundary conditions.
    subroutine unstagger_field_2d(f_stag, f, periodic_x, periodic_y)
        real(dp), dimension(:,:), intent(in) :: f_stag
        real(dp), dimension(:,:), intent(out) :: f
        logical, intent(in) :: periodic_x, periodic_y

        real(dp), dimension(4) :: pts

        real(dp) :: s,n

        integer :: i,j, k,i1, i2, j1, j2, ni, nj
        
        ni = size(f, 1)
        nj = size(f, 2)

        do i = 1, size(f, 1)
            do j = 1, size(f, 2)
                s = 0
                n = 0
                
                i1 = i-1
                i2 = i
                
                !If we're unstaggering with periodic boundaries, we cross over to the
                !other side of the domain when we "de-average".  Otherwise, we just ignore
                !the point that's off the domain.
                if (i1 == 0) then
                    if (periodic_y) then
                        i1 = ni - 1
                    else
                        i1 = 1
                    end if
                end if
    
                if (i2 == ni) then
                    if (periodic_y) then
                        i2 = 1
                    else
                        i2 = ni - 1
                    end if
                end if
    
                j1 = j-1
                j2 = j
    
                if (j1 == 0) then
                    if (periodic_y) then
                        j1 = nj - 1
                    else
                        j1 = 1
                    end if
                end if
    
                if (j2 == nj) then
                    if (periodic_x) then
                        j2 = 1
                    else
                        j2 = nj - 1
                    end if
                end if
                
                !Place the points into an array, loop over them, and average
                !all the points that AREN'T NaN.
                pts = (/f_stag(i1, j1), f_stag(i2, j1), f_stag(i1, j2), f_stag(i2, j2)/)
            
                do k=1,4
                    if (.not. (IS_NAN(pts(k)))) then
                        s = s + pts(k)
                        n = n + 1
                    end if
                end do
                
                if (n /= 0) then
                    f(i,j) = s/n
                else
                    f(i,j) = NaN
                end if

                !If the upper left of this location is not a number, then the staggered point must be
                !not a number.  This is to prevent dirichlet boundary conditions from being duplicated.
                if (IS_NAN(f_stag(i2,j2))) then
                    f(i,j) = NaN
                end if
            end do
        end do
    
    end subroutine unstagger_field_2d

    subroutine unstagger_field_3d(f, f_stag, periodic_x, periodic_y)
        real(dp), dimension(:,:,:) :: f, f_stag
        logical, intent(in) :: periodic_x, periodic_y

        integer :: i

        do i = 1,size(f,1)
            call unstagger_field_2d(f(i,:,:), f_stag(i,:,:), periodic_x, periodic_y)
        end do
        
    end subroutine unstagger_field_3d


    subroutine periodic_boundaries(m, apply_to_x, apply_to_y, nlayers_arg)
        !*FD Applies periodic boundary conditions to a 2D array
        real(dp), dimension(:,:), intent(inout) :: m
        integer :: maxx, maxy
        logical :: apply_to_x, apply_to_y
        integer, optional :: nlayers_arg

        integer :: nlayers 

        if (present(nlayers_arg)) then
            nlayers = nlayers_arg
        else
            nlayers = 1
        end if

        maxx = size(m, 1)
        maxy = size(m, 2)
       
        if (apply_to_x) then
            m( 1 : nlayers, : ) = m( maxx-nlayers*2 + 1 : maxx - nlayers, :)
            m( maxx-nlayers+1 : maxx, : ) = m(nlayers + 1 : nlayers*2, : )
        end if

        if (apply_to_y) then
            m( :, 1 : nlayers ) = m( :, maxy-nlayers*2 + 1 : maxy - nlayers )
            m( :, maxy-nlayers+1 : maxy ) = m( :, nlayers + 1 : nlayers*2 )
        end if

        
        !If both directions are periodic, treat the corners specially.
        if(apply_to_x .and. apply_to_y) then
            m(1:nlayers, 1:nlayers) = m(maxx-nlayers*2+1:maxx-nlayers, maxy-nlayers*2+1:maxy-nlayers)
            m(nlayers+1:2*nlayers, nlayers+1:2*nlayers) = m(maxx-nlayers+1:maxx, maxy-nlayers+1:maxy)
            
            m(1:nlayers, maxy-nlayers+1:maxy) = m(maxx-nlayers*2+1:maxx-nlayers, nlayers+1:2*nlayers)
            m(nlayers+1:2*nlayers, maxy-nlayers*2+1:maxy-nlayers) = m(maxx-nlayers+1:maxx, 1:nlayers)
        end if
    end subroutine periodic_boundaries
    
    subroutine periodic_boundaries_3d(m, apply_to_x, apply_to_y, nlayers_arg)
        !*FD Applies periodic boundary conditions to a 3D array
        real(dp), dimension(:,:,:), intent(inout) :: m
        logical :: apply_to_x, apply_to_y
        integer, optional :: nlayers_arg
    
        integer :: i
        
        do i = 1, size(m,1)
            call periodic_boundaries(m(i,:,:), apply_to_x, apply_to_y, nlayers_arg)
        end do
    end subroutine periodic_boundaries_3d
 
end module glide_grids
