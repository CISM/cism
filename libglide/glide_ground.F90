!TODO - Consider a glissade_ground (or glissade_calving?) module that contains
!        a parallel-friendly version of glide_marinlim.
!       Other subroutines in this module not currently used by parallel code.  
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_ground.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"
module glide_ground

  use glide_types
  use glimmer_global
  use parallel

  implicit none

contains
!-------------------------------------------------------------------------------  
!WHL - Removed this subroutine, which was used only by old case(5)

!!!  subroutine glide_initialise_backstress(thck,backstressmap,backstress,sigmabin,sigmabout)
!!!  end subroutine glide_initialise_backstress

!-------------------------------------------------------------------------

!WHL - Removed old case (5) as recommended by Jesse Johnson
!      Then changed old case(7) to new case(5) to avoid having a gap in the case numbering.
!TODO - Remove any other cases?
!       At some point, we should add a proper calving law for parallel code.

  subroutine glide_marinlim(which,                        &
                            thck,       relx,             &    
                            topg,       mask,             &
                            mlimit,     calving_fraction, &    
                            eus,        calving_field,    &
                            ground,                       &
                            dew,        dns,              &
                            nsn,        ewn)

    ! Remove non-grounded ice according to one of several alternative methods

    use glimmer_global, only : dp, sp

    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

!TODO: real(dp) for calving_field, eus?
!TODO: Change mask to thkmask?

    integer,                intent(in)    :: which   !*FD Calving method option
    real(dp),dimension(:,:),intent(inout) :: thck    !*FD Ice thickness
    real(dp),dimension(:,:),intent(in)    :: relx    !*FD Relaxed topography
    real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography
    integer, dimension(:,:), intent(in)   :: mask    !*FD grid type mask
    real(dp), intent(in)                  :: mlimit  !*FD Lower limit on topography elevation for
                                                     !*FD ice to be present. 
    real(dp), intent(in) :: calving_fraction         !*FD fraction of ice lost when calving; used with 
                                                     !*FD $\mathtt{which}=3$.
    real, intent(in) :: eus                          !*FD eustatic sea level
    real(sp),dimension(:,:),intent(out) :: calving_field ! thickness lost due to calving
    real(dp), intent(in) :: dew,dns
    integer, intent(in) ::  nsn,ewn

    type(glide_grnd), intent(inout) :: ground        !*FD ground instance

    integer :: ew,ns
    !---------------------------------------------------------------------
   
!HALO - Moved parallel_halo updates from this subroutine to glissade.F90

    calving_field(:,:) = 0.d0   ! using dp for constants in case calving_field changed to dp

    select case (which)

!WHL - debug - temporary case(0) to exercise this subroutine in the dome problem
    case(0) ! set thickness to zero at the ice margin, land or ocean

!SFP - commented out for now so that this option remains consistent w/ prev. versions of the code
! where case(0) was "do nothing".
!
!      print*, ' '
!      print*, 'Calving case 0: Remove ice at margin'
!      where (GLIDE_IS_MARGIN(mask))
!        calving_field = thck
!        thck = 0.0d0
!      end where
        
    case(1) ! Set thickness to zero if ice is floating

!LOOP TODO - Change to do loops over scalar cells?

      where (GLIDE_IS_FLOAT(mask))
        calving_field = thck
        thck = 0.0d0
      end where

    case(2) ! Set thickness to zero if relaxed bedrock is below a given level

!LOOP TODO - Change to do loops over scalar cells?

      where (relx <= mlimit+eus)
         calving_field = thck
         thck = 0.0d0
      end where

    case(3) ! remove fraction of ice when floating

!LOOP TODO - Why is the outer row of cells skipped here?

      do ns = 2,size(thck,2)-1
        do ew = 2,size(thck,1)-1
          if (GLIDE_IS_CALVING(mask(ew,ns))) then
            calving_field(ew,ns) = (1.d0-calving_fraction)*thck(ew,ns)
            thck(ew,ns) =  calving_fraction*thck(ew,ns)
            !mask(ew,ns) = ior(mask(ew,ns), GLIDE_MASK_OCEAN)
          end if
        end do
      end do

      ! if uncomment above mask update, then call parallel_halo(mask)

!TODO - Cases 2 and 4 are very similar; can we combine them?

    case(4) ! Set thickness to zero at marine edge if present bedrock is below a given level

!LOOP TODO - Change to do loops over scalar cells?

      where (GLIDE_IS_MARINE_ICE_EDGE(mask) .and. topg < mlimit+eus)
        calving_field = thck
        thck = 0.0d0
      end where

!WHL - Removed old case (5) based on recommendation from Jesse Johnson
!      Then changed old case(7) to new case(5) to avoid having a gap in the case numbering.

!TODO - not sure we want to support this case

    !Huybrechts grounding line scheme for Greenland initialization

    case(5)   ! used to be case(7)

      if(eus > -80.d0) then
        where (relx <= 2.d0*eus)
          calving_field = thck
          thck = 0.0d0
        end where
      elseif (eus <= -80.d0) then
        where (relx <= (2.d0*eus - 0.25d0*(eus + 80.d0)**2.d0))
          calving_field = thck
          thck = 0.0d0
        end where
      end if


!TODO - not sure we want to support this case
    case(6)

      ! not serial as far as I can tell as well; for parallelization, issues
      ! arise from components of ground being updated, and corresponding halos
      ! also need to be updated? Waiting until serial fixes are implemented

      call not_parallel(__FILE__, __LINE__) ! not serial as far as I can tell as well
      call update_ground_line(ground, topg, thck, eus, dew, dns, ewn, nsn, mask)

      where (GLIDE_IS_FLOAT(mask))
        calving_field = thck
        thck = 0.0d0
      end where
    
    end select

  end subroutine glide_marinlim

!HALO - This routine may not be supported.  I doubt it's accurate for HO flow.
!       If we do leave in this routine (or something similar), the loop should be 
!        over locally owned velocity points.

  !simple subroutine to calculate the flux at the grounding line

  subroutine calc_gline_flux(stagthk, velnorm, mask, gline_flux, ubas, vbas, dew)

    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
    implicit none

    !JEFF removing pointer attribute integer, dimension(:,:),pointer       :: mask    !*FD grid type mask
    integer, dimension(:,:)       :: mask                ! grid type mask
    real(dp),dimension(:,:),intent(in) :: stagthk        ! Ice thickness (scaled)
    real(dp),dimension(:,:,:), intent(in) :: velnorm     ! horizontal ice speed
    real(dp),dimension(:,:), intent(inout) :: gline_flux ! Grounding Line flux
    real(dp),dimension(:,:), intent(in) :: ubas          ! basal velocity in u-dir
    real(dp),dimension(:,:), intent(in) :: vbas          ! basal velocity in v-dir
    real(dp),intent(in)                 :: dew           !gridspacing  
    integer :: ewn, nsn

    !TODO: get the grounding line flux on the velo grid - right now it seems
    !to be using both the ice grid and the velo grid.
    ewn = size(gline_flux, 1)
    nsn = size(gline_flux, 2)
       
!LOOP - Loop over velocity points?

    where (GLIDE_IS_GROUNDING_LINE(mask))
         gline_flux = stagthk * ((4.0/5.0)* velnorm(1,:,:) + &
         (ubas**2.0 + vbas**2.0)**(1.0/2.0))  * dew  
    end where

!HALO - Pretty sure this is not needed.  gline_flux is just a diagnostic.
    call parallel_halo(gline_flux)
    call horiz_bcs_unstag_scalar(gline_flux)

  end subroutine calc_gline_flux

!TODO - The rest of this module is associated with case 6.
!       Not sure if it should be supported.

!-------------------------------------------------------------------------

!TODO - Is this function needed? Currently not called from anywhere I can find.

  !This function returns the correct grounding line using the data given 
  ! the mask reference point.  dir is specifying 'ew' or 'ns', but can be 
  ! left null if there's only one option.

  real function get_ground_line(ground,ew1,ns1,ew2,ns2)

     use glide_types
     implicit none
     type(glide_grnd) :: ground       !*FD glide ground instance
     integer ns1,ew1,ns2,ew2,slot_ns,slot_ew !grounding line in ns/ew direction
     real appr_ground !grounding line
     
     if (ns1 == ns2) then
         slot_ew = min(ew1,ew2)
         appr_ground = ground%gl_ns(slot_ew,ns1)
     else if (ew1 == ew2) then
         slot_ns = min(ns1,ns2)
         appr_ground = ground%gl_ew(ew1,slot_ns)
     end if
     get_ground_line = appr_ground
     return

  end function get_ground_line
    
!-------------------------------------------------------------------------

!TODO - Is this function needed? Only called from update_ground_line (case 6).

  subroutine set_ground_line(ground,ew1,ns1,ew2,ns2,value)

     use glide_types
     implicit none

     type(glide_grnd) :: ground        !*FD model instance
     integer, intent(in) :: ns1 !grounding line in ns direction
     integer, intent(in) :: ew1 !grounding line in ew direction
     integer, intent(in) :: ns2 !grounding line in ns direction
     integer, intent(in) :: ew2 !grounding line in ew direction
     real(sp), intent(in) :: value !grounding line in ew direction
     integer slot_ew, slot_ns !integers to compute the min
     
     if (ns1 == ns2) then
         slot_ew = min(ew1,ew2)
         ground%gl_ew(slot_ew,ns1) = value
     else if (ew1 == ew2) then
         slot_ns = min(ns1,ns2)
         ground%gl_ns(ew1,slot_ns) = value
     end if
  end subroutine set_ground_line

!-------------------------------------------------------------------------
!TODO - Is this needed? Only called from update_ground_line (case 6).

  !does the pattyn interpolation for the grounding line

  real function lin_reg_xg(topg, thck, eus, dew, dns, ew, ns, j1ew, j1ns)

     use glide_types
     use glimmer_physcon, only : rhoi, rhoo
     real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
     real(dp),dimension(:,:),intent(in)    :: thck    !*FD Present thickness (scaled)
     real, intent(in) :: eus                       !*FD eustatic sea level
     real(dp), intent(in) ::  dew, dns
     integer, intent(in) :: ns !grounding line in ns direction
     integer, intent(in) :: ew !grounding line in ew direction
     integer, intent(in) :: j1ns !ice shelf in ns direction
     integer, intent(in) :: j1ew !ice shelf line in ew direction
     real(sp) ::  xg                     !grounding line
     real(sp) ::  dx                      !distance between gridpts
     real(sp) ::  xj                        !grounding line
     real(sp) :: fj                        !f at grid pnt j
     real(sp) :: fj_1                      !f evaluated at j (+/-) 1
     real(sp) :: df                        !delta f of fj,jf_1
     
     if (ew == j1ew) then
        dx = dns 
        xj = ns*dx
     else
        dx = dew
        xj = ew*dx
     end if
     !set the pattyn f function - assuming ocean water 
     fj = (eus - topg(ew,ns))*rhoo/(rhoi*thck(ew,ns))
     if (thck(j1ew,j1ns) > 0.d0) then
         fj_1 = (eus - topg(j1ew,j1ns))*rhoo/(rhoi*thck(j1ew,j1ns))
         df = (fj_1 - fj)/dx
         xg = (1 - fj + df*xj)/df
     else
         xg = xj
     end if
     
     lin_reg_xg = xg
     return 
  end function lin_reg_xg

!-------------------------------------------------------------------------

!TODO - Is this needed? Only called from case 6.

  !Loops through the mask and does the interpolation for all the grounding lines

  subroutine update_ground_line(ground, topg, thck, eus, dew, dns, ewn, nsn, mask)

     implicit none
     type(glide_grnd) :: ground        !*FD ground instance
     real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
     real(dp),dimension(:,:),intent(in)    :: thck    !*FD Present thickness (scaled)
     real, intent(in) :: eus                       !*FD eustatic sea level
     real(dp), intent(in) ::  dew, dns
     integer, intent(in) ::  ewn, nsn
     !JEFF remove pointer attribute integer, dimension(:,:),pointer :: mask    !*FD grid type mask
     integer, dimension(:,:) :: mask    !*FD grid type mask
     integer ew,ns,jns,jew,j1ns,j1ew
     real(sp) :: xg                        !grounding line
     !this is assuming the grounding line is the last grounded pt on the mask
     !reset grounding line data to zero
     ground%gl_ew = 0.0
     ground%gl_ns = 0.0
     do ns = 1,nsn
        do ew = 1,ewn
            if (GLIDE_IS_GROUNDING_LINE(mask(ew,ns))) then
                !the grounding line always rounds down so it is grounded.
                !southern grounding line
                if (GLIDE_IS_OCEAN(mask(ew,ns - 1)) &
                        .or. (GLIDE_IS_FLOAT(mask(ew,ns - 1)))) then
                    xg = lin_reg_xg(topg,thck,eus,dew,dns,ew,ns,ew,ns-1)
                    call set_ground_line(ground,ew,ns,ew,ns-1,xg)
                !northern grounding line
                else if (GLIDE_IS_OCEAN(mask(ew,ns + 1)) &
                        .or. (GLIDE_IS_FLOAT(mask(ew,ns + 1)))) then
                    xg = lin_reg_xg(topg,thck,eus,dew,dns,ew,ns,ew,ns+1)
                    call set_ground_line(ground,ew,ns,ew,ns+1,xg) 
                end if 
                
                !western grounding line
                if (GLIDE_IS_OCEAN(mask(ew - 1,ns)) &
                        .or. GLIDE_IS_FLOAT(mask(ew - 1,ns))) then
                    xg = lin_reg_xg(topg,thck,eus,dew,dns,ew,ns,ew - 1,ns)
                    call set_ground_line(ground,ew,ns,ew-1,ns,xg)
                !eastern grounding line
                else if (GLIDE_IS_OCEAN(mask(ew + 1,ns)) &
                        .or. GLIDE_IS_FLOAT(mask(ew + 1,ns))) then
                    xg = lin_reg_xg(topg,thck,eus,dew,dns,ew,ns,ew + 1,ns)
                    call set_ground_line(ground,ew,ns,ew + 1,ns,xg)
                end if
            end if 
        end do
     end do

  end subroutine update_ground_line

  !-------------------------------------------------------------------------

!TODO - Is this needed?  Currently not called.

  real function get_ground_thck(ground,topg,usrf,dew,dns,ew1,ns1,ew2,ns2)

     use glide_types
     implicit none
     type(glide_grnd) :: ground        !*FD ground instance
     real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
     real(dp),dimension(:,:),intent(in)    :: usrf    !*FD surface height
     real(dp), intent(in) ::  dew, dns
     integer ns1,ew1,ns2,ew2,min_ns,min_ew,max_ns,max_ew !grounding line in ns/ew direction
     real(sp) ::  xg                        !grounding line
     real(sp) ::  tg                        !topographic height at grounding line
     real(sp) ::  ig                        !ice height at grounding line
     real(sp) ::  hg                        !thickness at the grounding line
     real(sp) ::  x1                        !pts for linear interpolation
     real(sp) ::  x0
     real(sp) ::  y1                        
     real(sp) ::  y0
     !using lin. interpolation to find top at grounding line
     if (ns1 == ns2) then
         min_ew = min(ew1,ew2)
         max_ew = max(ew1,ew2)
         min_ns = ns1
         max_ns = ns1
         x0 = min_ew*dew !model%numerics%dew
         x1 = max_ew*dew
     else if (ew1 == ew2) then
         min_ns = min(ns1,ns2)
         max_ns = max(ns1,ns2)
         min_ew = ew1
         max_ew = ew1
         x0 = min_ns*dns !model%numerics%dns
         x1 = max_ns*dns
     end if
     !get grounding line
     xg = ground%gl_ns(min_ew,min_ns)
     !find top height at xg
     y0 = topg(min_ew,min_ns) !model%geometry%topg
     y1 = topg(max_ew,max_ns)
     tg = y0 + (xg - x0)*((y1 - y0)/(x1 - x0))
     !find ice at xg
     y0 = usrf(min_ew,min_ns) !model%geometry%usrf
     y1 = usrf(max_ew,max_ns)
     ig = y0 + (xg - x0)*((y1 - y0)/(x1 - x0))
     !thickness
     hg = ig - tg
     get_ground_thck = hg
     return
  end function get_ground_thck
!---------------------------------------------------------------------------
end module glide_ground
