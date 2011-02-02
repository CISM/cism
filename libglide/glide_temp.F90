! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_temp.f90 - part of the Glimmer-CISM ice model    + 
! +  glide_temp.f90 - part of the GLIMMER ice model           +
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010
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

! JCC - Turn some of these changes into if-else statements so the
! "velocity_hom%uvel" and "velocity_hom%vvel" are only used when higher-order
! physics are on.
!
!whl - Dec. 2010: Moved several subroutines (finddisp, corrpmpt, calcpmpt, calcpmptb,
!       calcflwa, patebudd) to a new module, glide_temp_utils.
!      Those subroutines can then be called from an alternate temperature driver
!       without using glide_temp.
!      Also removed subroutines swappnt, swapbndt, and temperature_smoothing,
!       which were not being used.
! *sfp* this version contains hacks to replace SIA calc horiz vel fields 
! with their higher-order counterparts. That is, all "velocity%uvel" and
! "velocity%vvel" were replaced with "velocity_hom%uvel" and "velocity_hom%vvel".
! Note that dissip and basal melt calcs still need to be updated to account for
! higher-order stress and velocity fields.

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

! some macros used to disable parts of the temperature equation
! vertical diffusion
#ifdef NO_VERTICAL_DIFFUSION
#define VERT_DIFF 0.
#else
#define VERT_DIFF 1.
#endif

! horizontal advection
#ifdef NO_HORIZONTAL_ADVECTION
#define HORIZ_ADV 0.
#else
#define HORIZ_ADV 1.
#endif

! vertical advection
#ifdef NO_VERTICAL_ADVECTION
#define VERT_ADV 0.
#else
#define VERT_ADV 1.
#endif

! strain heating
#ifdef NO_STRAIN_HEAT
#define STRAIN_HEAT 0.
#else
#define STRAIN_HEAT 1.
#endif

module glide_temp

  use glide_types

contains

!------------------------------------------------------------------------------------

  subroutine glide_init_temp(model)

    !*FD initialise temperature module
    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow
    use glimmer_paramets, only : tim0, thk0, acc0, len0, vis0, vel0, tau0_glam
    use glimmer_global, only : dp 
    use glimmer_log
    use glide_bwater, only : find_dt_wat
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    integer, parameter :: p1 = gn + 1  
    integer up
    real(dp) :: estimate

    if (VERT_DIFF.eq.0.) call write_log('Vertical diffusion is switched off')
    if (HORIZ_ADV.eq.0.) call write_log('Horizontal advection is switched off')
    if (VERT_ADV.eq.0.) call write_log('Vertical advection is switched off')
    if (STRAIN_HEAT.eq.0.) call write_log('Strain heating is switched off')

    ! horizontal advection stuff
    allocate(model%tempwk%hadv_u(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%hadv_v(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%initadvt(model%general%upn,model%general%ewn,model%general%nsn))

    allocate(model%tempwk%inittemp(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%dissip(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%compheat(model%general%upn,model%general%ewn,model%general%nsn))
    model%tempwk%compheat = 0.0d0
    allocate(model%tempwk%dups(model%general%upn,3))

    allocate(model%tempwk%c1(model%general%upn))

    allocate(model%tempwk%dupa(model%general%upn),model%tempwk%dupb(model%general%upn))
    allocate(model%tempwk%dupc(model%general%upn))

    allocate(model%tempwk%smth(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%wphi(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatu(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatv(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%fluxew(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%fluxns(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bint(model%general%ewn-1,model%general%nsn-1))

    model%tempwk%advconst(1) = HORIZ_ADV*model%numerics%dttem / (16.0d0 * model%numerics%dew)
    model%tempwk%advconst(2) = HORIZ_ADV*model%numerics%dttem / (16.0d0 * model%numerics%dns)

    model%tempwk%dups = 0.0d0

    do up = 2, model%general%upn-1
       model%tempwk%dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) * &
            (model%numerics%sigma(up)   - model%numerics%sigma(up-1)))
       model%tempwk%dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) *  &
            (model%numerics%sigma(up+1) - model%numerics%sigma(up)))
       model%tempwk%dups(up,3) = 1.d0/(model%numerics%sigma(up+1)  - model%numerics%sigma(up-1))
    end do

    model%tempwk%zbed = 1.0d0 / thk0
    model%tempwk%dupn = model%numerics%sigma(model%general%upn) - model%numerics%sigma(model%general%upn-1)
    model%tempwk%wmax = 5.0d0 * tim0 / (scyr * thk0)

    model%tempwk%cons = (/ 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2), &
         model%numerics%dttem / 2.0d0, &
         VERT_DIFF*2.0d0 * tim0 * model%numerics%dttem / (thk0 * rhoi * shci), &
         VERT_ADV*tim0 * acc0 * model%numerics%dttem / coni, &
         ( tau0_glam * vel0 / len0 ) / ( rhoi * shci ) * ( model%numerics%dttem * tim0 ) /)  
         !*sfp* added last term to vector above for use in HO & SSA dissip. cacl

    model%tempwk%c1 = STRAIN_HEAT *(model%numerics%sigma * rhoi * grav * thk0**2 / len0)**p1 * &
         2.0d0 * vis0 * model%numerics%dttem * tim0 / (16.0d0 * rhoi * shci)

    model%tempwk%dupc = (/ (model%numerics%sigma(2) - model%numerics%sigma(1)) / 2.0d0, &
         ((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) / 2.0d0, &
         up=2,model%general%upn-1), (model%numerics%sigma(model%general%upn) - &
         model%numerics%sigma(model%general%upn-1)) / 2.0d0  /)
    model%tempwk%dupa = (/ 0.0d0, 0.0d0, &
         ((model%numerics%sigma(up) - model%numerics%sigma(up-1)) / &
         ((model%numerics%sigma(up-2) - model%numerics%sigma(up-1)) * &
         (model%numerics%sigma(up-2) - model%numerics%sigma(up))), &
         up=3,model%general%upn)  /)
    model%tempwk%dupb = (/ 0.0d0, 0.0d0, &
         ((model%numerics%sigma(up) - model%numerics%sigma(up-2)) / &
         ((model%numerics%sigma(up-1) - model%numerics%sigma(up-2)) * &
         (model%numerics%sigma(up-1) - model%numerics%sigma(up))), &
         up=3,model%general%upn)  /)
    
    model%tempwk%f = (/ tim0 * coni / (thk0**2 * lhci * rhoi), &
         tim0 / (thk0 * lhci * rhoi), &
         tim0 * thk0 * rhoi * shci /  (thk0 * tim0 * model%numerics%dttem * lhci * rhoi), &
         tim0 * thk0**2 * vel0 * grav * rhoi / (4.0d0 * thk0 * len0 * rhoi * lhci), &
         tim0 * vel0 * tau0_glam / (4.0d0 * thk0 * rhoi * lhci) /)      
         !*sfp* added the last term in the vect above for HO and SSA dissip. calc. 

    ! setting up some factors for sliding contrib to basal heat flux
    model%tempwk%slide_f = (/ VERT_DIFF * grav * thk0 * model%numerics%dttem/ shci, & ! vert diffusion
         VERT_ADV * rhoi*grav*acc0*thk0*thk0*model%numerics%dttem/coni /)             ! vert advection

    select case(model%options%whichbwat)
       case(0)
          model%paramets%hydtim = tim0 / (model%paramets%hydtim * scyr)
          estimate = 0.2d0 / model%paramets%hydtim
          !EIB! following not in lanl glide_temp
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 
          
          model%tempwk%c = (/ model%tempwk%dt_wat, 1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
               1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 
       case(1)
          !EIB! not in lanl
          !model%tempwk%watvel = model%paramets%hydtim * tim0 / (scyr * len0)
          !estimate = (0.2d0 * model%tempwk%watvel) / min(model%numerics%dew,model%numerics%dns)
          !call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 
          
          !print *, model%numerics%dttem*tim0/scyr, model%tempwk%dt_wat*tim0/scyr, model%tempwk%nwat

          !model%tempwk%c = (/ rhow * grav, rhoi * grav, 2.0d0 * model%numerics%dew, 2.0d0 * model%numerics%dns, &
          !     0.25d0 * model%tempwk%dt_wat / model%numerics%dew, 0.25d0 * model%tempwk%dt_wat / model%numerics%dns, &
          !     0.5d0 * model%tempwk%dt_wat / model%numerics%dew, 0.5d0 * model%tempwk%dt_wat / model%numerics%dns /)
          
       end select
       ! JCC - was in LANL but not in parallel
       !      model%temper%temp = -10.0

  end subroutine glide_init_temp

!****************************************************    

  subroutine glide_temp_driver(model,whichtemp)

    !*FD Calculates the ice temperature, according to one
    !*FD of several alternative methods.

    use glimmer_utils, only: tridiag
    use glimmer_global, only : dp
    use glimmer_paramets, only : thk0
    use glide_velo
    use glide_thck
    use glide_grids
    use glide_bwater
    use glide_temp_utils

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.
    integer,                intent(in)    :: whichtemp   !*FD Flag to choose method.

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp),dimension(size(model%numerics%sigma)) :: subd, diag, supd, rhsd
    real(dp),dimension(size(model%numerics%sigma)) :: prevtemp, iteradvt, diagadvt
    real(dp) :: tempresid

    integer :: iter
    integer :: ew,ns

    real(dp),parameter :: tempthres = 0.001d0, floatlim = 10.0d0 / thk0
    integer, parameter :: mxit = 100
    integer, parameter :: ewbc = 1, nsbc = 1 

    real(dp), dimension(size(model%numerics%sigma)) :: weff

    !------------------------------------------------------------------------------------
    ! ewbc/nsbc set the type of boundary condition aplied at the end of
    ! the domain. a value of 0 implies zero gradient.
    !------------------------------------------------------------------------------------

    select case(whichtemp)

    case(0) ! Set column to surface air temperature -------------------------------------

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns)))
          end do
       end do

    case(1) ! Do full temperature solution ---------------------------------------------

!lipscomb - restart mod - These routines are now called at the end of tstep_p3, so that wgrd
!                         can be written to the hotstart file and used for restart.
                         
       ! Calculate time-derivatives of thickness and upper surface elevation ------------

       call timeders(model%thckwk,   &
            model%geometry%thck,     &
            model%geomderv%dthckdtm, &
            model%geometry%mask,     &
            model%numerics%time,     &
            1)

       call timeders(model%thckwk,   &
            model%geometry%usrf,     &
            model%geomderv%dusrfdtm, &
            model%geometry%mask,     &
            model%numerics%time,     &
            2)
       !JCC - Don't use the ho velo fields unless we're using ho physics
       if (model%options%which_ho_diagnostic .EQ. 0 ) then
          ! Calculate the vertical velocity of the grid ------------------------------------
    
          call gridwvel(model%numerics%sigma,  &
                model%numerics%thklim, &
                model%velocity%uvel,   &
                model%velocity%vvel,   &
                model%geomderv,        &
                model%geometry%thck,   &
                model%velocity%wgrd)
    
          ! Calculate the actual vertical velocity; method depends on whichwvel ------------
    
          select case(model%options%whichwvel)
          case(0) 
    
              ! Usual vertical integration
    
              call wvelintg(model%velocity%uvel,                        &
                  model%velocity%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%velocity%wgrd(model%general%upn,:,:), &
                  model%geometry%thck,                        &
                  model%temper%bmlt,                          &
                  model%velocity%wvel)
    
          case(1)
    
              ! Vertical integration constrained so kinematic upper BC obeyed.
    
              call wvelintg(model%velocity%uvel,                        &
                  model%velocity%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%velocity%wgrd(model%general%upn,:,:), &
                  model%geometry%thck,                        &
                  model%temper%  bmlt,                        &
                  model%velocity%wvel)
    
              call chckwvel(model%numerics,                             &
                  model%geomderv,                             &
                  model%velocity%uvel(1,:,:),                 &
                  model%velocity%vvel(1,:,:),                 &
                  model%velocity%wvel,                        &
                  model%geometry%thck,                        &
                  model%climate% acab)
          end select
       else ! using ho physics
          ! Calculate the vertical velocity of the grid ------------------------------------
    
          call gridwvel(model%numerics%sigma,  &
                model%numerics%thklim, &
                model%velocity_hom%uvel,   &
                model%velocity_hom%vvel,   &
                model%geomderv,        &
                model%geometry%thck,   &
                model%velocity%wgrd)
    
          ! Calculate the actual vertical velocity; method depends on whichwvel ------------
    
          select case(model%options%whichwvel)
          case(0) 
    
              ! Usual vertical integration
    
              call wvelintg(model%velocity_hom%uvel,                        &
                  model%velocity_hom%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%velocity%wgrd(model%general%upn,:,:), &
                  model%geometry%thck,                        &
                  model%temper%bmlt,                          &
                  model%velocity%wvel)
    
          case(1)
    
              ! Vertical integration constrained so kinematic upper BC obeyed.
    
              call wvelintg(model%velocity_hom%uvel,                        &
                  model%velocity_hom%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%velocity%wgrd(model%general%upn,:,:), &
                  model%geometry%thck,                        &
                  model%temper%  bmlt,                        &
                  model%velocity%wvel)
    
              call chckwvel(model%numerics,                             &
                  model%geomderv,                             &
                  model%velocity_hom%uvel(1,:,:),                 &
                  model%velocity_hom%vvel(1,:,:),                 &
                  model%velocity%wvel,                        &
                  model%geometry%thck,                        &
                  model%climate% acab)
          end select
       end if ! model%options%which_ho_diagnostic .EQ. 0

       ! Calculate the vertical velocity of the grid ------------------------------------

       call gridwvel(model%numerics%sigma,  &
            model%numerics%thklim, &
            model%velocity_hom%uvel,   &
            model%velocity_hom%vvel,   &
            model%geomderv,        &
            model%geometry%thck,   &
            model%velocity%wgrd)

       ! Calculate the actual vertical velocity; method depends on whichwvel ------------

       select case(model%options%whichwvel)

       case(0) 

          ! Usual vertical integration

          call wvelintg(model%velocity_hom%uvel,           &
               model%velocity_hom%vvel,                    &
               model%geomderv,                             &
               model%numerics,                             &
               model%velowk,                               &
               model%velocity%wgrd(model%general%upn,:,:), &
               model%geometry%thck,                        &
               model%temper%bmlt,                          &
               model%velocity%wvel)

       case(1)

          ! Vertical integration constrained so kinematic upper BC obeyed.

          call wvelintg(model%velocity_hom%uvel,           &
               model%velocity_hom%vvel,                    &
               model%geomderv,                             &
               model%numerics,                             &
               model%velowk,                               &
               model%velocity%wgrd(model%general%upn,:,:), &
               model%geometry%thck,                        &
               model%temper%  bmlt,                        &
               model%velocity%wvel)

          call chckwvel(model%numerics,                    &
               model%geomderv,                             &
               model%velocity_hom%uvel(1,:,:),             &
               model%velocity_hom%vvel(1,:,:),             &
               model%velocity%wvel,                        &
               model%geometry%thck,                        &
               model%climate% acab)

       end select   ! whichwvel

       ! apply periodic ew BC
       if (model%options%periodic_ew) then
          call wvel_ew(model)
       end if

       model%tempwk%inittemp = 0.0d0
       model%tempwk%initadvt = 0.0d0
       !*MH model%tempwk%dissip   = 0.0d0  is also set to zero in finddisp
       ! ----------------------------------------------------------------------------------

       call finddisp(model,          &
            model%geometry%thck,     &
            model%options%which_disp,&
            model%velocity_hom%efvs, &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%temper%flwa)
    
       !JCC - Don't use ho velocity fields unless we're using the ho model
       if (model%options%which_ho_diagnostic .EQ. 0 ) then
          ! translate velo field
          do ns = 2,model%general%nsn-1
              do ew = 2,model%general%ewn-1
                model%tempwk%hadv_u(:,ew,ns) = model%tempwk%advconst(1) * ( model%velocity%uvel(:,ew-1,ns-1) &
                    + model%velocity%uvel(:,ew-1,ns) + model%velocity%uvel(:,ew,ns-1) + model%velocity%uvel(:,ew,ns) )
                model%tempwk%hadv_v(:,ew,ns) = model%tempwk%advconst(2) * ( model%velocity%vvel(:,ew-1,ns-1) &
                    + model%velocity%vvel(:,ew-1,ns) + model%velocity%vvel(:,ew,ns-1) + model%velocity%vvel(:,ew,ns) )
              end do
          end do
      else ! Using ho physics
          ! translate velo field
          do ns = 2,model%general%nsn-1
              do ew = 2,model%general%ewn-1
                model%tempwk%hadv_u(:,ew,ns) = model%tempwk%advconst(1) * ( model%velocity_hom%uvel(:,ew-1,ns-1) &
                    + model%velocity_hom%uvel(:,ew-1,ns) + model%velocity_hom%uvel(:,ew,ns-1) + model%velocity_hom%uvel(:,ew,ns) )
                model%tempwk%hadv_v(:,ew,ns) = model%tempwk%advconst(2) * ( model%velocity_hom%vvel(:,ew-1,ns-1) &
                    + model%velocity_hom%vvel(:,ew-1,ns) + model%velocity_hom%vvel(:,ew,ns-1) + model%velocity_hom%vvel(:,ew,ns) )
              end do
          end do
      end if ! model%options%which_ho_diagnostic .EQ. 0

       call hadvall(model, &
            model%temper%temp, &
            model%geometry%thck)

       ! zeroth iteration
       iter = 0
       tempresid = 0.0d0

       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1
             if(model%geometry%thck(ew,ns)>model%numerics%thklim) then

                weff = model%velocity%wvel(:,ew,ns) - model%velocity%wgrd(:,ew,ns)
                if (maxval(abs(weff)) > model%tempwk%wmax) then
                   weff = 0.0d0
                end if

                call hadvpnt(iteradvt,                       &
                     diagadvt,                               &
                     model%temper%temp(:,ew-2:ew+2,ns),      &
                     model%temper%temp(:,ew,ns-2:ns+2),      &
                     model%tempwk%hadv_u(:,ew,ns), &
                     model%tempwk%hadv_v(:,ew,ns))
               
                call findvtri(model,ew,ns,subd,diag,supd,diagadvt, &
                     weff, &
                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                call findvtri_init(model,ew,ns,subd,diag,supd,weff,model%temper%temp(:,ew,ns), &
                     model%geometry%thck(ew,ns),GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                call findvtri_rhs(model,ew,ns,model%climate%artm(ew,ns),iteradvt,rhsd, &
                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                prevtemp(:) = model%temper%temp(:,ew,ns)

                call tridiag(subd(1:model%general%upn), &
                             diag(1:model%general%upn), &
                             supd(1:model%general%upn), &
                             model%temper%temp(1:model%general%upn,ew,ns), &
                             rhsd(1:model%general%upn))

                call corrpmpt(model%temper%temp(:,ew,ns),     &
                              model%geometry%thck(ew,ns),     &
                              model%temper%bwat(ew,ns),       &
                              model%numerics%sigma,           &
                              model%general%upn)

                tempresid = max(tempresid,maxval(abs(model%temper%temp(:,ew,ns)-prevtemp(:))))

             endif
          end do
       end do

       do while (tempresid.gt.tempthres .and. iter.le.mxit)
          tempresid = 0.0d0

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if(model%geometry%thck(ew,ns)>model%numerics%thklim) then

                   weff = model%velocity%wvel(:,ew,ns) - model%velocity%wgrd(:,ew,ns)
                   if (maxval(abs(weff)) > model%tempwk%wmax) then
                      weff = 0.0d0
                   end if

                   call hadvpnt(iteradvt,                       &
                        diagadvt,                               &
                        model%temper%temp(:,ew-2:ew+2,ns),      &
                        model%temper%temp(:,ew,ns-2:ns+2),      &
                        model%tempwk%hadv_u(:,ew,ns), &
                        model%tempwk%hadv_v(:,ew,ns))

                   call findvtri(model,ew,ns,subd,diag,supd,diagadvt, &
                        weff, &
                        GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                   call findvtri_rhs(model,ew,ns,model%climate%artm(ew,ns),iteradvt,rhsd, &
                        GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                   prevtemp(:) = model%temper%temp(:,ew,ns)

                   call tridiag(subd(1:model%general%upn), &
                        diag(1:model%general%upn), &
                        supd(1:model%general%upn), &
                        model%temper%temp(1:model%general%upn,ew,ns), &
                        rhsd(1:model%general%upn))

                   call corrpmpt(model%temper%temp(:,ew,ns),     &
                                 model%geometry%thck(ew,ns),     &
                                 model%temper%bwat(ew,ns),       &
                                 model%numerics%sigma,           &
                                 model%general%upn)


                   tempresid = max(tempresid, maxval(abs(model%temper%temp(:,ew,ns)-prevtemp(:))))

                endif
             end do
          end do

          iter = iter + 1
       end do

       model%temper%niter = max(model%temper%niter, iter )
       
       ! set temperature of thin ice to the air temperature and set ice free nodes to zero
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
                model%temper%temp(:,ew,ns) = min(0.0d0,dble(model%climate%artm(ew,ns)))
             else if (GLIDE_NO_ICE(model%geometry%thkmask(ew,ns))) then
                !JCC - one of these
                !model%temper%temp(:,ew,ns) = min(0.0d0,dble(model%climate%artm(ew,ns)))
                model%temper%temp(:,ew,ns) = 0.0d0
             end if
          end do
       end do

       ! apply periodic ew BC
       if (model%options%periodic_ew) then
          model%temper%temp(:,0,:) = model%temper%temp(:,model%general%ewn-2,:)
          model%temper%temp(:,1,:) = model%temper%temp(:,model%general%ewn-1,:)
          model%temper%temp(:,model%general%ewn,:) = model%temper%temp(:,2,:)
          model%temper%temp(:,model%general%ewn+1,:) = model%temper%temp(:,3,:)
       end if

       ! Calculate basal melt rate --------------------------------------------------

       call calcbmlt(model, &
            model%options%which_bmelt, &
            model%temper%temp, &
            model%geometry%thck, &
            model%geomderv%stagthck, &
            model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns, &
            model%velocity%ubas, &
            model%velocity%vbas, &
            model%temper%bmlt, &
            GLIDE_IS_FLOAT(model%geometry%thkmask))

       ! Calculate basal water depth ------------------------------------------------

       !JCC - Using lanl's bwater calculation routines.
       call calcbwat(model, &
            model%options%whichbwat, &
            model%temper%bmlt, &
            model%temper%bwat, &
            model%temper%bwatflx, &
            model%geometry%thck, &
            model%geometry%topg, &
            model%temper%temp(model%general%upn,:,:), &
            GLIDE_IS_FLOAT(model%geometry%thkmask), &
            model%tempwk%wphi)

       !EIB! rest of case not present in lanl
       ! Transform basal temperature and pressure melting point onto velocity grid -
!    case(2) ! *sfp* stealing this un-used option ... 

!        call stagvarb(model%temper%temp(model%general%upn,1:model%general%ewn,1:model%general%nsn), &
!             model%temper%stagbtemp ,&
!             model%general%  ewn, &
!             model%general%  nsn)
!        
!        call calcbpmp(model,model%geometry%thck,model%temper%bpmp)
! 
!        call stagvarb(model%temper%bpmp, &
!             model%temper%stagbpmp ,&
!             model%general%  ewn, &
!             model%general%  nsn)

!    case(2) ! Do something else, unspecified ---------------------------------------
!
!       do ns = 1,model%general%nsn
!          do ew = 1,model%general%ewn
!             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns))) * (1.0d0 - model%numerics%sigma)
!             call corrpmpt(model%temper%temp(:,ew,ns),model%geometry%thck(ew,ns),model%temper%bwat(ew,ns),&
!                  model%numerics%sigma,model%general%upn)
!          end do
!       end do


    case(2) ! *sfp* stealing this un-used option ... 

        ! DO NOTHING. That is, hold T const. at initially assigned value

    end select   ! whichtemp

    ! Calculate Glenn's A --------------------------------------------------------

    call calcflwa(model%numerics%sigma,        &
                  model%numerics%thklim,       &
                  model%temper%flwa,           &
                  model%temper%temp(:,1:model%general%ewn,1:model%general%nsn), &
                  model%geometry%thck,         &
                  model%paramets%flow_factor,  &
                  model%paramets%default_flwa, &
                  model%options%whichflwa) 

    ! Output some information ----------------------------------------------------

#ifdef DEBUG
    print *, "* temp ", model%numerics%time, iter, model%temper%niter, &
         real(model%temper%temp(model%general%upn,model%general%ewn/2+1,model%general%nsn/2+1))
#endif

  end subroutine glide_temp_driver

  !-------------------------------------------------------------------------

  subroutine hadvpnt(iteradvt,diagadvt,tempx,tempy,u,v)

    use glimmer_global, only : dp

    implicit none

    real(dp), dimension(:), intent(in) :: u,v
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    if (u(1) > 0.0d0) then
       iteradvt = u * (- 4.0d0*tempx(:,2) + tempx(:,1))
       diagadvt = u * 3.0d0
    else if (u(1) < 0.0d0) then
       iteradvt = u * (4.0d0*tempx(:,4) - tempx(:,5))
       diagadvt = - u * 3.0d0
    end if

    if (v(1) > 0.0d0) then
       iteradvt = iteradvt + v * (- 4.0d0*tempy(:,2) + tempy(:,1))
       diagadvt = diagadvt + v * 3.0d0
    else if (v(1) < 0.0d0) then
       iteradvt = iteradvt + v * (4.0d0*tempy(:,4) - tempy(:,5))
       diagadvt = diagadvt - v * 3.0d0
    end if

  end subroutine hadvpnt

  !-------------------------------------------------------------------------

  subroutine fohadvpnt(tempwk,iteradvt,diagadvt,tempx,tempy,uvel,vvel)

    use glimmer_global, only : dp
    use glimmer_utils, only: hsum

    implicit none

    type(glide_tempwk) :: tempwk
    real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt

    real(dp), dimension(size(iteradvt)) :: u, v

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    u = tempwk%advconst(1) * hsum(uvel(:,:,:))
    v = tempwk%advconst(2) * hsum(vvel(:,:,:))

    if (u(1) > 0.0d0) then
       iteradvt = - u * 2.0d0 * tempx(:,1)
       diagadvt = 2.0d0 * u 
    else if (u(1) < 0.0d0) then
       iteradvt = u * 2.0d0 * tempx(:,3)
       diagadvt = - 2.0d0 * u 
    end if

    if (v(1) > 0.0d0) then
       iteradvt = iteradvt - v * 2.0d0 * tempy(:,1) 
       diagadvt = diagadvt + 2.0d0 * v 
    else if (v(1) < 0.0d0) then
       iteradvt = iteradvt + v * 2.0d0 * tempy(:,3)
       diagadvt = diagadvt - 2.0d0 * v 
    end if

  end subroutine fohadvpnt

  !-------------------------------------------------------------------------

  subroutine hadvall(model,temp,thck)

    use glimmer_global, only : dp 

    implicit none

    type(glide_global_type) :: model
    real(dp), dimension(:,0:,0:), intent(in) :: temp
    real(dp), dimension(:,:), intent(in) :: thck

    real(dp), dimension(size(temp,dim=1)) :: diagadvt

    integer :: ew,ns

    model%tempwk%initadvt = 0.0d0

    do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then

             call hadvpnt(model%tempwk%initadvt(:,ew,ns), &
                  diagadvt,                       &
                  temp(:,ew-2:ew+2,ns),           &
                  temp(:,ew,ns-2:ns+2),           &
                  model%tempwk%hadv_u(:,ew,ns), &
                  model%tempwk%hadv_v(:,ew,ns))
          end if
       end do
    end do

  end subroutine hadvall

  !-------------------------------------------------------------------------

  subroutine findvtri(model,ew,ns,subd,diag,supd,diagadvt,weff,float)

    use glimmer_global, only : dp

    implicit none

    type(glide_global_type) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(in) :: weff,  diagadvt
    real(dp), dimension(:), intent(out) :: subd, diag, supd
    logical, intent(in) :: float

    real(dp) :: fact(3)

! These constants are precomputed:
! model%tempwk%cons(1) = 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)
! model%tempwk%cons(2) = model%numerics%dttem / 2.0d0 

    fact(1) = VERT_DIFF*model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2
    fact(2) = VERT_ADV*model%tempwk%cons(2) / model%geometry%thck(ew,ns)    
    
    subd(2:model%general%upn-1) = fact(2) * weff(2:model%general%upn-1) * &
         model%tempwk%dups(2:model%general%upn-1,3)

    supd(2:model%general%upn-1) = - subd(2:model%general%upn-1) - fact(1) * &
         model%tempwk%dups(2:model%general%upn-1,2)

    subd(2:model%general%upn-1) = subd(2:model%general%upn-1) - fact(1) * &
         model%tempwk%dups(2:model%general%upn-1,1)

    diag(2:model%general%upn-1) = 1.0d0 - subd(2:model%general%upn-1) &
         - supd(2:model%general%upn-1) &
         + diagadvt(2:model%general%upn-1)

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0

    ! now do the basal boundary
    ! for grounded ice, a heat flux is applied
    ! for floating ice, temperature held constant

     if (float) then

       supd(model%general%upn) = 0.0d0
       subd(model%general%upn) = 0.0d0
       diag(model%general%upn) = 1.0d0

     else 
 
        supd(model%general%upn) = 0.0d0 
        subd(model%general%upn) = -0.5*fact(1)/(model%tempwk%dupn**2)
        diag(model%general%upn) = 1.0d0 - subd(model%general%upn) + diagadvt(model%general%upn)
 
     end if

  end subroutine findvtri

  !-------------------------------------------------------------------------

  subroutine findvtri_init(model,ew,ns,subd,diag,supd,weff,temp,thck,float)
    !*FD called during first iteration to set inittemp
    use glimmer_global, only : dp
    implicit none
    type(glide_global_type) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(in) :: temp,diag,subd,supd,weff
    real(dp), intent(in) :: thck
    logical, intent(in) :: float    

    ! local variables
    real(dp) :: slterm
    integer ewp,nsp
    integer slide_count

    model%tempwk%inittemp(2:model%general%upn-1,ew,ns) = temp(2:model%general%upn-1) * &
         (2.0d0 - diag(2:model%general%upn-1)) &
         - temp(1:model%general%upn-2) * subd(2:model%general%upn-1) &
         - temp(3:model%general%upn) * supd(2:model%general%upn-1) & 
         - model%tempwk%initadvt(2:model%general%upn-1,ew,ns) &
         + model%tempwk%dissip(2:model%general%upn-1,ew,ns)
    
    if (float) then
       model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) 
       !EIB old!model%tempwk%inittemp(model%general%upn,ew,ns) = pmpt(thck)
    else 
       ! sliding contribution to basal heat flux
       slterm = 0.
       slide_count = 0

       !whl - BUG! - The following expression for taub*ubas is valid only for the SIA
       !             Need a different expression for HO dynamics

       ! only include sliding contrib if temperature node is surrounded by sliding velo nodes
       do nsp = ns-1,ns
          do ewp = ew-1,ew
             if (abs(model%velocity%ubas(ewp,nsp)).gt.0.000001 .or. abs(model%velocity%vbas(ewp,nsp)).gt.0.000001) then
                slide_count = slide_count + 1
                slterm = slterm + (&
                     model%geomderv%dusrfdew(ewp,nsp) * model%velocity%ubas(ewp,nsp) + &
                     model%geomderv%dusrfdns(ewp,nsp) * model%velocity%vbas(ewp,nsp))
             end if
          end do
       end do
       if (slide_count.ge.4) then
          slterm = 0.25*slterm
       else
          slterm = 0.
       end if
       model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) * &
            (2.0d0 - diag(model%general%upn)) &
            - temp(model%general%upn-1) * subd(model%general%upn) &
            - 0.5*model%tempwk%cons(3) * model%temper%bheatflx(ew,ns) / (thck * model%tempwk%dupn) & ! geothermal heat flux (diff)
            - model%tempwk%slide_f(1)*slterm/ model%tempwk%dupn &                                    ! sliding heat flux    (diff)
            - model%tempwk%cons(4) * model%temper%bheatflx(ew,ns) * weff(model%general%upn) &        ! geothermal heat flux (adv)
            - model%tempwk%slide_f(2)*thck*slterm* weff(model%general%upn) &                         ! sliding heat flux    (adv)
            - model%tempwk%initadvt(model%general%upn,ew,ns)  &
            + model%tempwk%dissip(model%general%upn,ew,ns)
    end if

  end subroutine findvtri_init

  !-----------------------------------------------------------------------

  subroutine findvtri_rhs(model,ew,ns,artm,iteradvt,rhsd,float)

    !*FD RHS of temperature tri-diag system
    use glimmer_global, only : dp, sp 
    implicit none
    type(glide_global_type) :: model
    integer, intent(in) :: ew, ns
    real(sp), intent(in) :: artm 
    real(dp), dimension(:), intent(in) :: iteradvt
    real(dp), dimension(:), intent(out) :: rhsd
    logical, intent(in) :: float    

    ! upper boundary condition
    rhsd(1) = artm
    if (float) then
       rhsd(model%general%upn) = model%tempwk%inittemp(model%general%upn,ew,ns)    
    else
       rhsd(model%general%upn) = model%tempwk%inittemp(model%general%upn,ew,ns) - iteradvt(model%general%upn)
    end if
    rhsd(2:model%general%upn-1) = model%tempwk%inittemp(2:model%general%upn-1,ew,ns) - iteradvt(2:model%general%upn-1)

  end subroutine findvtri_rhs

  !-----------------------------------------------------------------------

  subroutine calcbmlt(model,whichbmelt,temp,thck,stagthck,dusrfdew,dusrfdns,ubas,vbas,bmlt,floater)

    use glimmer_global, only : dp 
    use glide_temp_utils, only: calcpmpt

    implicit none 

    type(glide_global_type) :: model
    real(dp), dimension(:,0:,0:), intent(in) :: temp
    real(dp), dimension(:,:), intent(in) :: thck,  stagthck, dusrfdew, dusrfdns, ubas, vbas  
    real(dp), dimension(:,:), intent(inout) :: bmlt
    logical, dimension(:,:), intent(in) :: floater
    integer, intent(in) :: whichbmelt

    real(dp), dimension(size(model%numerics%sigma)) :: pmptemp
    real(dp) :: slterm, newmlt
 
    integer :: ewp, nsp, up, ew, ns

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim .and. .not. floater(ew,ns)) then

             call calcpmpt(pmptemp,thck(ew,ns),model%numerics%sigma)

             if (abs(temp(model%general%upn,ew,ns)-pmptemp(model%general%upn)) .lt. 0.001) then

                slterm = 0.0d0

                select case( whichbmelt )    !*sfp* added for calculating differently based on model physics

                case( SIA_BMELT )                   ! 0-order SIA approx. --> Tau_d = Tau_b                                     

                    do nsp = ns-1,ns
                        do ewp = ew-1,ew
                            slterm = slterm - stagthck(ewp,nsp) * &
                            (dusrfdew(ewp,nsp) * ubas(ewp,nsp) + dusrfdns(ewp,nsp) * vbas(ewp,nsp))
                        end do
                    end do

                    !*sfp* NOTE that multiplication by this term has been moved up from below
                    slterm = model%tempwk%f(4) * slterm 


                case( FIRSTORDER_BMELT )                   ! 1st-order SIA approx. (HO model)
                    do nsp = ns-1,ns 
                        do ewp = ew-1,ew
                        !*sfp* Note that vel and stress arrays have diff vert dims (upn for vel, upn-1 for stress)
                        ! so that for now, basal vel at upn is multiplied by basal stress at upn-1 to get frictional
                        ! heating term. This may not be entirely correct ... 
                             slterm = slterm - &
                                
                                !! NEW version: uses consistent basal tractions                
                                ( -model%velocity_hom%btraction(1,ewp,nsp) * &
                                   model%velocity_hom%uvel(model%general%upn,ewp,nsp)  &
                                  -model%velocity_hom%btraction(2,ewp,nsp) * &
                                   model%velocity_hom%vvel(model%general%upn,ewp,nsp) )
                                
                                !!!! OLD version: uses HO basal shear stress calc. from FD                
                                !( model%velocity_hom%tau%xz(model%general%upn-1,ewp,nsp) * &
                                !  model%velocity_hom%uvel(model%general%upn,ewp,nsp) + &
                                !  model%velocity_hom%tau%yz(model%general%upn-1,ewp,nsp) * &
                                !  model%velocity_hom%vvel(model%general%upn,ewp,nsp) )

                        end do
                    end do

                    slterm = model%tempwk%f(5) * slterm

!                case( SSA_BMELT )                  ! 1st-order, depth-integrated approx. (SSA) 
!                                            ! NOTE: need to pass 2d basal shear stress arrays from SSA model
!                   do nsp = ns-1,ns 
!                         do ewp = ew-1,ew
!                             slterm = slterm - &
!                                 ( taubxs(upn,ewp,nsp) * ubas(ewp,nsp) + taubys(upn,ewp,nsp) * vbas(ewp,nsp))
!                         end do
!                   end do
!
!                   slterm = model%tempwk%f(5) * slterm

                end select

                bmlt(ew,ns) = 0.0d0

                !*sfp* changed this so that 'slterm' is multiplied by f(4) const. above ONLY for the 0-order SIA case,
                ! since for the HO and SSA cases a diff. const. needs to be used

                ! OLD version
!                newmlt = model%tempwk%f(4) * slterm - model%tempwk%f(2)*model%temper%bheatflx(ew,ns) + model%tempwk%f(3) * &
!                     model%tempwk%dupc(model%general%upn) * &
!                     thck(ew,ns) * model%tempwk%dissip(model%general%upn,ew,ns)

                ! NEW version (sfp)
                newmlt = slterm - model%tempwk%f(2)*model%temper%bheatflx(ew,ns)   &
                        + model%tempwk%f(3) * model%tempwk%dupc(model%general%upn) * &
                          thck(ew,ns) * model%tempwk%dissip(model%general%upn,ew,ns)

                up = model%general%upn - 1

                do while (abs(temp(up,ew,ns)-pmptemp(up)) .lt. 0.001 .and. up .ge. 3)
                   bmlt(ew,ns) = bmlt(ew,ns) + newmlt
                   newmlt = model%tempwk%f(3) * model%tempwk%dupc(up) * thck(ew,ns) * model%tempwk%dissip(up,ew,ns)
                   up = up - 1
                end do

                up = up + 1

                if (up == model%general%upn) then
                   bmlt(ew,ns) = newmlt - &
                        model%tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * model%tempwk%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * model%tempwk%dupb(up) ) / thck(ew,ns) 
                else
                   bmlt(ew,ns) = bmlt(ew,ns) + max(0.0d0, newmlt - &
                        model%tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * model%tempwk%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * model%tempwk%dupb(up) ) / thck(ew,ns)) 
                end if

             else

                bmlt(ew,ns) = 0.0d0

             end if

          !EIB! else if (model%options%use_plume == 1) then

          ! do nothing because the plume model will have written the bmlt field
          else

              bmlt(ew,ns) = 0.0d0

          end if
       end do
    end do

    ! apply periodic BC

    if (model%options%periodic_ew) then
       do ns = 2,model%general%nsn-1
          bmlt(1,ns) = bmlt(model%general%ewn-1,ns)
          bmlt(model%general%ewn,ns) = bmlt(2,ns)
       end do
    end if
  end subroutine calcbmlt

!-------------------------------------------------------------------

end module glide_temp

