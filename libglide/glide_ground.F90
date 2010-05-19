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
  implicit none
contains
  
!-------------------------------------------------------------------------------  
  subroutine glide_initialise_backstress(thck,backstressmap,backstress,sigmabin,sigmabout)
  
     implicit none
     
     real(dp), dimension(:,:), intent(in) :: thck !*FD Ice thickness
     logical, dimension(:,:), intent(inout) :: backstressmap !*FD Backstress map
     real(sp), dimension(:,:), intent(inout):: backstress !*FD Backstress
     real(sp) :: sigmabin,sigmabout
     backstress = 0.0
     backstressmap = .False. 
     where(thck > 0.0) 
         backstressmap = .True. 
     end where
     where (thck > 0.0)
         backstress = sigmabin
     elsewhere
         backstress = sigmabout
     end where
  end subroutine glide_initialise_backstress
!-------------------------------------------------------------------------

  subroutine glide_marinlim(which,thck,relx,topg,flwa,levels,mask,mlimit,calving_fraction,eus,ablation_field,backstress, & 
                 tempanmly,dew,dns,backstressmap,stressout,stressin,ground,nsn,ewn,usrf)


    !*FD Removes non-grounded ice, according to one of two altenative
    !*FD criteria, and sets upper surface of non-ice-covered points 
    !*FD equal to the topographic height, or sea-level, whichever is higher.

    use glimmer_global, only : dp, sp
    use glimmer_physcon, only : rhoi, rhoo, grav, gn
    use glide_vertint, only : vertint_output2d
    use glimmer_paramets, only: thk0
    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer,                intent(in)    :: which   !*FD Option to choose ice-removal method
                                                     !*FD \begin{description}
                                                     !*FD \item[0] Set thickness to zero if 
                                                     !*FD relaxed bedrock is below a given level.
                                                     !*FD \item[1] Set thickness to zero if
                                                     !*FD ice is floating.
                                                     !*FD \end{description}
    real(dp),dimension(:,:),intent(inout) :: thck    !*FD Ice thickness (scaled)
    real(dp),dimension(:,:),intent(in)    :: relx    !*FD Relaxed topography (scaled)
    real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
    real(dp),dimension(:,:,:),intent(in)  :: flwa    !*FD Vertically averaged ice hardness
    real(dp),dimension(:)    ,intent(in)  :: levels    !*FD Vertically averaged ice hardness
    integer, dimension(:,:),pointer       :: mask    !*FD grid type mask
    real(dp)                              :: mlimit  !*FD Lower limit on topography elevation for
                                                     !*FD ice to be present (scaled). Used with 
                                                     !*FD $\mathtt{which}=0$.
    real(dp), intent(in) :: calving_fraction         !*FD fraction of ice lost when calving Used with 
                                                     !*FD $\mathtt{which}=3$.
    real, intent(inout) :: eus                       !*FD eustatic sea level
    real(sp),dimension(:,:),intent(inout) :: ablation_field !*FD this is passed as climate%calving

    real(dp), parameter :: con = - rhoi / rhoo
    real(dp), parameter :: sigmaxx = 0.5 * rhoi * grav * (1.0 - rhoi / rhoo)
    real(dp), parameter :: theta = 0.125
    real(dp), dimension(2,2) :: A
    real(sp) :: pi =  3.141592654 

    real(sp), dimension(:,:), intent(inout) :: backstress
    real(sp) :: tempanmly
    real(dp), intent(in) ::  dew,dns
    integer, intent(in) ::  nsn,ewn
    real(dp),dimension(:,:),intent(in)    :: usrf    !*FD usrf
    logical, dimension(:,:), intent(in)   :: backstressmap !*FD map of the
                                                           !*FD backstresses for the initial map
    integer ew,ns
    real(sp) :: stressout,stressin
    type(glide_grnd) :: ground        !*FD ground instance
    !---------------------------------------------------------------------
   

    ablation_field=0.0

    select case (which)
        
    case(1) ! Set thickness to zero if ice is floating
      where (GLIDE_IS_FLOAT(mask))
        ablation_field=thck
        thck = 0.0d0
      end where

    case(2) ! Set thickness to zero if relaxed bedrock is below a 
       ! given level
       where (relx <= mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where
    
    case(3) ! remove fraction of ice when floating
       do ns = 2,size(thck,2)-1
          do ew = 2,size(thck,1)-1
             if (GLIDE_IS_CALVING(mask(ew,ns))) then
                ablation_field(ew,ns)=(1.0-calving_fraction)*thck(ew,ns)
                thck(ew,ns) =  calving_fraction*thck(ew,ns)
                !mask(ew,ns) = ior(mask(ew,ns), GLIDE_MASK_OCEAN)
             end if
          end do
       end do

    case(4) ! Set thickness to zero at marine edge if present
            ! bedrock is below a given level
       where (GLIDE_IS_MARINE_ICE_EDGE(mask).and.topg<mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where
    case(5) ! Relation based on computing the horizontal stretching
            ! of the unconfined ice shelf (\dot \epsilon_{xx}) and multiplying by H.
            ! 
        do ns = 2, size(backstress,2)-1
           do ew = 2, size(backstress,1)-1
                if(.not. backstressmap(ew,ns)) then
                   !should be > -1.0 if using log10
                   if (tempanmly > 0.0) then
                      backstress(ew,ns) = stressout
                   else
                      !backstress(ew,ns) = stressout + (1-stressout)*log10(-tempanmly)
                      !( 1-exp(tempanmly))
                      backstress(ew,ns) = stressout + (1-stressout)*abs(tempanmly/9.2)
                      ! backstress(ew,ns) =sigmabout + (1-sigmabout)*atan(-tempanmly)/(pi/2)
                     
                   end if
                else
                   !should be > -1.0 if using log10
                   if (tempanmly > 0.0) then
                      backstress(ew,ns) = stressin
                   else
                     backstress(ew,ns) = stressin + (1-stressin)*abs(tempanmly/9.2)
                      !backstress(ew,ns) =stressin + (1-stressin)*log10(-tempanmly)

                     
                     !backstress(ew,ns) =sigmabin + (1-sigmabin)*atan(-tempanmly)/(pi/2)
                     
                   end if
                end if
           end do
        end do 
        
        !do ns = 2, size(backstress,2)-1
        !   do ew = 2, size(backstress,1)-1
        !where (backstressmap) 
        !    if (tempanmly > 0.0) then
        !      backstress(ew,ns) = sigmabin
        !    else
            !backstress = sigmabin + (1-sigmabin)*abs(tempanmly/9.2)
        !      backstress(ew,ns) =sigmabin + (1-sigmabin)*( 1-exp(tempanmly))
        !    end if 
        !   end do
       ! end do
        !end where 
        
        where(backstress > 1.0)
          backstress = 1.0
        end where
           
           do ns = 2,size(thck,2)-1
              do ew = 2,size(thck,1)-1
                 if (GLIDE_IS_GROUNDING_LINE(mask(ew,ns))) then
                    call vertint_output2d(flwa(:,ew-1:ew,ns-1:ns),A, levels * thck(ew,ns))
                    ablation_field(ew,ns)= ((dew*dns)/(50.0d3)**2)* theta * A(2,2) * (sigmaxx * &
                    thck(ew,ns)  * (1 - backstress(ew,ns))) ** gn
                    if ((thck(ew,ns) - ablation_field(ew,ns)) >= 0.0) then
                      thck(ew,ns) = thck(ew,ns) - ablation_field(ew,ns) 
                    else 
                      ablation_field(ew,ns) = thck(ew,ns)
                      thck(ew,ns) = 0.0d0
                      
                    end if
                end if
              end do
           end do
          
          !where (GLIDE_IS_FLOAT(mask).and.relx<mlimit+eus)
          !    ablation_field=thck
          !    thck = 0.0d0
          ! end where
          
          !remove all ice that is outside of a single grid layer of floating ice
          !adjacent to the grounding line
          do ns = 2,size(thck,2)-1
             do ew = 2,size(thck,1)-1
               if (GLIDE_IS_FLOAT(mask(ew,ns)) .and. .not. backstressmap(ew,ns))then 
          
                 if ((.not. GLIDE_IS_GROUNDING_LINE(mask(ew-1,ns))) &
                  .and. (.not. GLIDE_IS_GROUNDING_LINE(mask(ew+1,ns))) .and. &
                  (.not. GLIDE_IS_GROUNDING_LINE(mask(ew,ns-1))) .and. &
                  (.not. GLIDE_IS_GROUNDING_LINE(mask(ew,ns+1)))) then
                       ablation_field(ew,ns) = thck(ew,ns)
                       
                       thck(ew,ns) = 0.0
                 end if
               end if
             end do
           end do
    case(6)
       call update_ground_line(ground, topg, thck, eus, dew, dns, ewn, nsn, mask)
       where (GLIDE_IS_FLOAT(mask))
        ablation_field=thck
        thck = 0.0d0
       end where
    
    !Huybrechts grounding line scheme for Greenland initialization
    case(7)
       if(eus > -80.0) then
       where(relx <= 2.0*eus)
      
          ablation_field=thck
          thck = 0.0d0
       end where
       elseif(eus <= -80.0) then
       where ( relx <= (2.0*eus - 0.25*(eus + 80.0)**2.0))
          ablation_field = thck
          thck = 0.0d0
       end where
       end if


    end select
  end subroutine glide_marinlim

  !simple subroutine to calculate the flux at the grounding line
  subroutine calc_gline_flux(stagthk, surfvel, mask, gline_flux, ubas, vbas, &
  dew)
    implicit none
    integer, dimension(:,:),pointer       :: mask    !*FD grid type mask
    real(dp),dimension(:,:),intent(in) :: stagthk    !*FD Ice thickness (scaled)
    real(dp),dimension(:,:,:), intent(in) :: surfvel !*FD Surface velocity
    real(dp),dimension(:,:), intent(inout) :: gline_flux !*FD Grounding Line flux
    real(dp),dimension(:,:), intent(in) :: ubas !*FD basal velocity in u-dir
    real(dp),dimension(:,:), intent(in) :: vbas !*FD basal velocity in v-dir
    real(dp),intent(in)                 :: dew !*FD gridspacing  
    integer :: ewn, nsn

    !TODO: get the grounding line flux on the velo grid - right now it seems
    !to be using both the ice grid and the velo grid.
    ewn = size(gline_flux, 1)
    nsn = size(gline_flux, 2)

    
       
    where (GLIDE_IS_GROUNDING_LINE(mask))
         gline_flux = stagthk * ((4.0/5.0)* surfvel(1,:,:) + &
         (ubas**2.0 + vbas**2.0)**(1.0/2.0))  * dew  
    end where
  end subroutine calc_gline_flux

!-------------------------------------------------------------------------
  !This function returns the correct grounding line using the data given 
  ! the mask reference point.  dir is specifying 'ew' or 'ns', but can be 
  ! left null if there's only one option.
  real function get_ground_line(ground,ew1,ns1,ew2,ns2)
     use glide_types
     implicit none
     type(glide_grnd) :: ground       !*FD glide ground instance
     integer ns1,ew1,ns2,ew2,slot_ns,slot_ew !grounding line in ns/ew direction
     real appr_ground !grounding line
     
     if (ns1 .eq. ns2) then
         slot_ew = min(ew1,ew2)
         appr_ground = ground%gl_ns(slot_ew,ns1)
     else if (ew1 .eq. ew2) then
         slot_ns = min(ns1,ns2)
         appr_ground = ground%gl_ew(ew1,slot_ns)
     end if
     get_ground_line = appr_ground
     return
  end function get_ground_line
    
!-------------------------------------------------------------------------
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
     
     if (ns1 .eq. ns2) then
         slot_ew = min(ew1,ew2)
         ground%gl_ew(slot_ew,ns1) = value
     else if (ew1 .eq. ew2) then
         slot_ns = min(ns1,ns2)
         ground%gl_ns(ew1,slot_ns) = value
     end if
  end subroutine set_ground_line
!-------------------------------------------------------------------------
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
     
     if (ew .eq. j1ew) then
        dx = dns 
        xj = ns*dx
     else
        dx = dew
        xj = ew*dx
     end if
     !set the pattyn f function - assuming ocean water 
     fj = (eus - topg(ew,ns))*rhoo/(rhoi*thck(ew,ns))
     if (thck(j1ew,j1ns) .gt. 0.0) then
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
  !Loops through the mask and does the interpolation for all the grounding lines
  subroutine update_ground_line(ground, topg, thck, eus, dew, dns, ewn, nsn, mask)
     implicit none
     type(glide_grnd) :: ground        !*FD ground instance
     real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
     real(dp),dimension(:,:),intent(in)    :: thck    !*FD Present thickness (scaled)
     real, intent(in) :: eus                       !*FD eustatic sea level
     real(dp), intent(in) ::  dew, dns
     integer, intent(in) ::  ewn, nsn
     integer, dimension(:,:),pointer :: mask    !*FD grid type mask
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
     if (ns1 .eq. ns2) then
         min_ew = min(ew1,ew2)
         max_ew = max(ew1,ew2)
         min_ns = ns1
         max_ns = ns1
         x0 = min_ew*dew !model%numerics%dew
         x1 = max_ew*dew
     else if (ew1 .eq. ew2) then
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
