!Module for the thickness mask, pulled out of glide_setup to avoid circular dependencies now
!that the mask function is needed in glide_thck
module glide_thckmask
   implicit none
contains
  subroutine glide_maskthck(thck,massb,include_adjacent,thklim,dom,pointno,totpts,empty)
    
    !*FD Calculates the contents of the mask array.

    use glimmer_global, only : dp, sp 

    !-------------------------------------------------------------------------
    ! Subroutine arguments
    !-------------------------------------------------------------------------

    real(dp),dimension(:,:),intent(in)  :: thck      !*FD Ice thickness
    real(sp),dimension(:,:),intent(in)  :: massb      !*FD Mass balance
    logical,                intent(in)  :: include_adjacent 
    !*FD If true, points with no ice but that are adjacent to points with ice 
    !*FD are included in the mask
    real(dp),               intent(in)  :: thklim
    !*FD Ice dynamics thickness limit, set to 0 to include all ice.
    integer, dimension(:),  intent(out) :: dom        
    integer, dimension(:,:),intent(out) :: pointno    !*FD Output mask
    integer,                intent(out) :: totpts     !*FD Total number of points
    logical,                intent(out) :: empty      !*FD Set if no mask points set.

    !-------------------------------------------------------------------------
    ! Internal variables
    !-------------------------------------------------------------------------

    integer,dimension(size(thck,2),2) :: band
    logical,dimension(size(thck,2))   :: full
    integer :: covtot 
    integer :: ew,ns,ewn,nsn

    !-------------------------------------------------------------------------

    ewn=size(thck,1) ; nsn=size(thck,2)

    pointno = 0
    covtot  = 0 

    !-------------------------------------------------------------------------

    do ns = 1,nsn

      full(ns) = .false.

      do ew = 1,ewn
        if ( thckcrit(thck(max(1,ew-1):min(ewn,ew+1),max(1,ns-1):min(nsn,ns+1)),massb(ew,ns)) ) then

          covtot = covtot + 1
          pointno(ew,ns) = covtot 

          if ( .not. full(ns) ) then
            band(ns,1) = ew
            full(ns)   = .true.
          else
            band(ns,2) = ew
          end if
               
        end if
      end do
    end do
  
    totpts = covtot
                                             
    dom(1:2) = (/ewn,1/); empty = .true.

    do ns = 1,nsn
           
      if (full(ns)) then

        if (empty) then
          empty  = .false.
          dom(3) = ns
        end if
        dom(4) = ns
        dom(1) = min0(dom(1),band(ns,1))
        dom(2) = max0(dom(2),band(ns,2))
      end if
    end do

  contains
    logical function thckcrit(ca,cb)

      implicit none

      real(dp),dimension(:,:),intent(in) :: ca 
      real(sp),               intent(in) :: cb

      if (.not. include_adjacent) then
        ! Include only points with ice in the mask

        if ( ca(2,2) > thklim .or. cb > thklim) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      else

        ! If the thickness in the region under consideration
        ! or the mass balance is positive, thckcrit is .true.
        ! This means that the mask includes points that have no
        ! ice but are adjacent to points that do have ice
        if ( any((ca(:,:) > thklim)) .or. cb > thklim ) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      end if

    end function thckcrit

  end subroutine glide_maskthck


end module glide_thckmask
