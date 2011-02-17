! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glissade_temp.f90 - part of the GLIMMER ice model        +
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

! This module computes temperature diffusion and strain heating
!  in a local column without doing horizontal or vertical advection.
! Temperature advection is done separately, e.g. using the incremental 
!  remapping transport scheme.
! It is assumed here that temperature values are staggered in the
!  vertical compared to the velocity.  That is, the temperature lives
!  at the midpoint of each layer instead of at layer interfaces.
! As in the unstaggered case, the temperature is also defined at the 
!  upper and lower surfaces with appropriate boundary conditions.   

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade_temp

  use glide_types

contains

!****************************************************    

  subroutine glissade_init_temp (model)

    ! initialization subroutine for the case whichtemp = TEMP_REMAP_ADV,
    !      with a vertically staggered temperature that is advected elsewhere
    !      (i.e., not here or in glide_temp)

    !*FD initialise temperature module
    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow
    use glimmer_paramets, only : tim0, thk0, len0, vis0, vel0, tau0_glam
    use glimmer_global, only : dp 
    use glimmer_log

    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    integer, parameter :: p1 = gn + 1  
    integer up
    real(dp) :: estimate

!whl - Note vertical dimensions here.  Dissipation is computed for each of (upn-1) layers.
!      Temperature is defined at midpoint of each layer, plus upper and lower surfaces.
    allocate(model%tempwk%dups(model%general%upn+1,2))
    allocate(model%tempwk%inittemp(model%general%upn+1,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%dissip  (model%general%upn-1,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%compheat(model%general%upn-1,model%general%ewn,model%general%nsn))
    model%tempwk%compheat = 0.0d0

    allocate(model%tempwk%c1(model%general%upn-1))   ! upn-1 for staggered grid

    allocate(model%tempwk%dupa(model%general%upn))
    allocate(model%tempwk%dupb(model%general%upn))
    allocate(model%tempwk%dupc(model%general%upn))

    allocate(model%tempwk%smth(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%wphi(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatu(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatv(model%general%ewn,model%general%nsn))

    model%tempwk%dups = 0.0d0

!whl - Note that the 'dups' grid coefficients are not the same as for unstaggered temperatures.

    up = 1
    model%tempwk%dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                    (model%numerics%stagsigma(up) - model%numerics%sigma(up)) )
    do up = 2, model%general%upn-1
       model%tempwk%dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                       (model%numerics%stagsigma(up) - model%numerics%stagsigma(up-1)) )
    enddo

    do up = 1, model%general%upn-2
       model%tempwk%dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                       (model%numerics%stagsigma(up+1) - model%numerics%stagsigma(up)) )
    end do
    up = model%general%upn-1
    model%tempwk%dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                    (model%numerics%sigma(up+1) - model%numerics%stagsigma(up)) )

    model%tempwk%zbed = 1.0d0 / thk0
    model%tempwk%dupn = model%numerics%sigma(model%general%upn) - model%numerics%sigma(model%general%upn-1)

!whl - We need only two of these (cons(1) and cons(5)) for the staggered case with no advection.
!!    model%tempwk%cons = (/ 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2), &
!!         model%numerics%dttem / 2.0d0, &
!!         VERT_DIFF*2.0d0 * tim0 * model%numerics%dttem / (thk0 * rhoi * shci), &
!!         VERT_ADV*tim0 * acc0 * model%numerics%dttem / coni, &
!!         0.0d0, &   !whl - no vertical advection
!!         ( tau0_glam * vel0 / len0 ) / ( rhoi * shci ) * ( model%numerics%dttem * tim0 ) /)  
         !*sfp* added last term to vector above for use in HO & SSA dissip. cacl

!whl - The factor of 2 in the numerator in the original code can be traced to a missing factor of 0.5
!      in the denominator of the dups coefficients.  On the vertically staggered grid, there is no
!      factor of 0.5 in the dups coefficients, so there is no factor of 2 here.
!whl - The factor of 2 in the denominator is a Crank-Nicolson averaging factor.
!!    model%tempwk%cons(1) = 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)
    model%tempwk%cons(1) = tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)

    model%tempwk%cons(5) = (tau0_glam * vel0 / len0 ) / (rhoi * shci) * (model%numerics%dttem * tim0)

!whl - Note: stagsigma here instead of sigma
    model%tempwk%c1(1:model%general%upn-1) =   &
                             (model%numerics%stagsigma(1:model%general%upn-1)   &
                             * rhoi * grav * thk0**2 / len0)**p1 * 2.0d0 * vis0              &
                             * model%numerics%dttem * tim0 / (16.0d0 * rhoi * shci)

    model%tempwk%f = (/ tim0 * coni / (thk0**2 * lhci * rhoi), &
                        tim0 / (thk0 * lhci * rhoi), &
                        tim0 * thk0 * rhoi * shci /  (thk0 * tim0 * model%numerics%dttem * lhci * rhoi), &
                        tim0 * thk0**2 * vel0 * grav * rhoi / (4.0d0 * thk0 * len0 * rhoi * lhci), &
                        tim0 * vel0 * tau0_glam / (4.0d0 * thk0 * rhoi * lhci) /)      
                        !*sfp* added the last term in the vect above for HO and SSA dissip. calc. 

    select case(model%options%whichbwat)
       case(0)
          model%paramets%hydtim = tim0 / (model%paramets%hydtim * scyr)
          estimate = 0.2d0 / model%paramets%hydtim
          
          model%tempwk%c = (/ model%tempwk%dt_wat,   &
                              1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
                              1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
                              0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 
       case(1)
          
    end select
  
  end subroutine glissade_init_temp

!****************************************************    

  subroutine glissade_temp_driver(model)

    ! Calculates the ice temperature 

    use glimmer_utils,  only : tridiag
    use glimmer_global, only : dp
    use glimmer_paramets, only : thk0, tim0
    use glimmer_physcon, only: shci, coni, rhoi
    use glide_mask
    use glide_bwater
    use glide_temp_utils
    use glimmer_log

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model       !*FD Ice model parameters.

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    integer :: ew, ns, up, upn
    character(len=100) :: message

    ! These arrays have the same size as the vertical dimension of temperature:
    !  upn+1 on the staggered grid.
    real(dp), dimension(size(model%temper%temp,1)) :: subd, diag, supd, rhsd

    ! These have the same dimensions as staggered temperature
    real(dp),dimension(0:model%general%upn) :: Tstagsigma, prevtemp_stag

    ! for energy conservation check
    real(dp) :: einit, efinal, delta_e, dTtop, dTbot

    upn = model%general%upn

    ! local column calculation
    ! No horizontal or vertical advection; vertical diffusion and strain heating only.
    ! Temperatures are vertically staggered relative to velocities.  That is, the temperature is
    !  defined at the midpoint of each layer (and at the top and bottom surfaces).

    ! Set Tstagsigma (= stagsigma except that it has values at the top and bottom surfaces).

    Tstagsigma(0) = 0.d0
    Tstagsigma(1:model%general%upn-1) = model%numerics%stagsigma(1:model%general%upn-1)
    Tstagsigma(model%general%upn) = 1.d0

    model%tempwk%inittemp = 0.0d0

    ! Calculate interior heat dissipation -------------------------------------

    call finddisp( model,                   &
                   model%geometry%thck,     &
                   model%options%which_disp,&
                   model%velocity_hom%efvs, &
                   model%geomderv%stagthck, &
                   model%geomderv%dusrfdew, &
                   model%geomderv%dusrfdns, &
                   model%temper%flwa)

    ! Calculate heating from basal friction -----------------------------------

    call calcbfric( model,                        &
                    model%options%which_bmelt,      &
                    model%geometry%thck,          &
                    model%velocity_hom%btraction, &
                    model%geomderv%dusrfdew,      &
                    model%geomderv%dusrfdns,      &
                    model%velocity%ubas,          &
                    model%velocity%vbas,          &
                    GLIDE_IS_FLOAT(model%geometry%thkmask) )

    ! Note: No iteration is needed here since we are doing a local tridiagonal solve without advection.

    do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          if(model%geometry%thck(ew,ns) > model%numerics%thklim) then

             ! compute initial internal energy in column (for energy conservation check)
             einit = 0.0d0
             do up = 1, upn-1
                einit = einit + model%temper%temp(up,ew,ns) *  &
                               (model%numerics%sigma(up+1) -   &
                                model%numerics%sigma(up) )
             enddo
             einit = einit * rhoi * shci * model%geometry%thck(ew,ns)*thk0

             ! compute matrix elements

             call glissade_findvtri( model, ew,   ns,             &
                                     subd,  diag, supd, rhsd,     &            
                                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)), &
                                     model%options%which_bmelt )
                
             prevtemp_stag(:) = model%temper%temp(:,ew,ns)

             ! solve the tridiagonal system

             ! Note: Temperature is indexed from 0 to upn, with indices 1 to upn-1 colocated
             !  with stagsigma values of the same index.
             ! However, the matrix elements are indexed 1 to upn+1, with the first row
             !  corresponding to the surface temperature, temp(0,:,:).

             call tridiag(subd(1:model%general%upn+1), &
                          diag(1:model%general%upn+1), &
                          supd(1:model%general%upn+1), &
                          model%temper%temp(0:model%general%upn,ew,ns), &
                          rhsd(1:model%general%upn+1))

             ! Check that the net input of energy to the column is equal to the difference
             !  between the initial and final internal energy.
             !whl - to do - Make this check optional and/or move it to a subroutine.

             ! first compute the final internal energy

             efinal = 0.0d0
             do up = 1, upn-1
                efinal = efinal + model%temper%temp(up,ew,ns) *  &
                                 (model%numerics%sigma(up+1) - model%numerics%sigma(up))
             enddo
             efinal = efinal * rhoi*shci * model%geometry%thck(ew,ns)*thk0

             ! compute net heat flux to the column

             ! conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down

             dTtop = 0.5d0 * ( model%temper%temp(1,ew,ns) - model%temper%temp(0,ew,ns) &
                             +     prevtemp_stag(1)       -     prevtemp_stag(0) )
             model%temper%ucondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                           * dTtop / (Tstagsigma(1) - Tstagsigma(0))


             dTbot = 0.5d0 * ( model%temper%temp(upn,ew,ns) - model%temper%temp(upn-1,ew,ns) &
                            +     prevtemp_stag(upn)       -   prevtemp_stag(upn-1) )
             model%temper%lcondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                           * dTbot / (Tstagsigma(upn) - Tstagsigma(upn-1))

             ! total dissipation in column

             model%temper%dissipcol = 0.0d0
             do up = 1, upn-1
                model%temper%dissipcol = model%temper%dissipcol + &
                                         model%tempwk%dissip(up,ew,ns)  &
                                        * (model%numerics%sigma(up+1) - model%numerics%sigma(up))  
             enddo 
             model%temper%dissipcol = model%temper%dissipcol     &
                                     * thk0*model%geometry%thck(ew,ns)*rhoi*shci / (tim0*model%numerics%dttem)  

             ! Verify that the net input of energy into the column is equal to the change in
             ! internal energy.  

             delta_e = (model%temper%ucondflx(ew,ns) - model%temper%lcondflx(ew,ns)  &
                      + model%temper%dissipcol(ew,ns)) * tim0*model%numerics%dttem

!whl - would need a different threshold if running in single precision
             if ( abs((efinal-einit-delta_e)/(tim0*model%numerics%dttem)) > 1.0e-8 ) then
                write(message,*) 'WARNING: Energy conservation error, ew, ns =', ew, ns
                call write_log(message)
! Can uncomment the following for diagnostics
!                write(50,*) 'Interior fluxes:'
!                write(50,*) 'ftop (pos up)=', -model%temper%ucondflx(ew,ns) 
!                write(50,*) 'fbot (pos up)=', -model%temper%lcondflx(ew,ns)
!                write(50,*) 'fdissip =',       model%temper%dissipcol(ew,ns)
!                write(50,*) 'Net flux =', delta_e/(tim0*model%numerics%dttem)
!                write(50,*) ' '
!                write(50,*) 'delta_e =', delta_e
!                write(50,*) 'einit =',  einit
!                write(50,*) 'efinal =', efinal
!                write(50,*) 'einit + delta_e =', einit + delta_e
!                write(50,*) ' '
!                write(50,*) 'Energy imbalance =', efinal - einit - delta_e
!                write(50,*) ' '
!                write(50,*) 'Basal fluxes:'
!                write(50,*) 'ffric =', model%temper%bfricflx(ew,ns)
!                write(50,*) 'fgeo =', -model%temper%bheatflx(ew,ns)
!                write(50,*) 'flux for bottom melting =', model%temper%bfricflx(ew,ns)   &
!                                                             - model%temper%bheatflx(ew,ns)   &
!                                                             + model%temper%lcondflx(ew,ns)
             endif

!whl - No call here to corrpmpt.  Temperatures above pmpt are set to pmpt 
!      in calcbmlt_remapadv (conserving energy).

          endif  ! thck > thklim
       end do    ! ew
    end do       ! ns

    ! set temperature of thin ice to the air temperature and set ice-free nodes to zero
    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn

          if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
             model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
          else if (model%geometry%thkmask(ew,ns) < 0) then
             model%temper%temp(:,ew,ns) = min(0.0d0, dble(model%climate%artm(ew,ns)))
          !else if (model%geometry%thkmask(ew,ns) < -1) then
          !   model%temper%temp(:,ew,ns) = 0.0d0
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

    ! Calculate basal melt rate
    ! Temperature above the pressure melting point are reset to Tpmp,
    !  with excess heat contributing to melting.

    call glissade_calcbmlt( model,                     &
                            model%options%which_bmelt, &
                            model%temper%temp,         &
                            Tstagsigma,                &
                            model%geometry%thck,       &
                            model%geomderv%stagthck,   &
                            model%geomderv%dusrfdew,   &
                            model%geomderv%dusrfdns,   &
                            model%velocity%ubas,       &
                            model%velocity%vbas,       &
                            model%temper%bmlt,         &
                            GLIDE_IS_FLOAT(model%geometry%thkmask))

    !whl - to do - Should ice thickness be reduced as a result of basal melting?

    ! Calculate basal water depth ------------------------------------------------

    call calcbwat( model,                     &
                   model%options%whichbwat,   &
                   model%temper%bmlt,         &
                   model%temper%bwat,         &
                   model%temper%bwatflx,      &
                   model%geometry%thck,       &
                   model%geometry%topg,       &
                   model%temper%temp(model%general%upn,:,:), &
                   GLIDE_IS_FLOAT(model%geometry%thkmask),   &
                   model%tempwk%wphi)

    ! Calculate Glenn's A --------------------------------------------------------

    call calcflwa(model%numerics%stagsigma,    &
                  model%numerics%thklim,       &
                  model%temper%flwa,           &
                  model%temper%temp(1:model%general%upn-1,:,:),  &
                  model%geometry%thck,         &
                  model%paramets%flow_factor,  &
                  model%paramets%default_flwa, &
                  model%options%whichflwa) 

  end subroutine glissade_temp_driver

  !-------------------------------------------------------------------------

  subroutine glissade_findvtri (model, ew,   ns,          &
                                subd,  diag, supd, rhsd,  &
                                float, whichbmlt)

    ! compute matrix elements for the tridiagonal solve

    use glimmer_global, only : dp
    use glimmer_paramets, only : thk0
    use glimmer_physcon,  only : rhoi, grav, coni
    use glide_temp_utils, only: calcpmptb

    implicit none

!whl - Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
!             whereas temperature is indexed from 0 to upn.
!            The first row of the matrix is the equation for temp(0,ew,ns),
!             the second row is the equation for temp(1,ew,ns), and so on.

    type(glide_global_type), intent(inout) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    logical, intent(in) :: float
    integer, intent(in) :: whichbmlt   ! SIA or higher-order

    ! local variables

    real(dp) :: pmptempb  ! pressure melting temp at bed
    real(dp) :: fact

    ! set surface temperature

    model%temper%temp(0,ew,ns) = dble(model%climate%artm(ew,ns))

    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements

    ! upper boundary: set to surface air temperature

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = model%temper%temp(0,ew,ns)

    ! ice interior. layers 1:upn-1  (matrix elements 2:upn)

    ! model%tempwk%cons(1) = 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)

    fact = model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2
    subd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,1)
    supd(2:model%general%upn) = -fact * model%tempwk%dups(1:model%general%upn-1,2)
    diag(2:model%general%upn) = 1.0d0 - subd(2:model%general%upn)     &
                                      - supd(2:model%general%upn)

    model%tempwk%inittemp(1:model%general%upn-1,ew,ns) =   &
           model%temper%temp(1:model%general%upn-1,ew,ns) * (2.0d0 - diag(2:model%general%upn)) &
         - model%temper%temp(0:model%general%upn-2,ew,ns) * subd(2:model%general%upn) &
         - model%temper%temp(2:model%general%upn,  ew,ns) * supd(2:model%general%upn) & 
         + model%tempwk%dissip(1:model%general%upn-1,ew,ns)
    
    rhsd(2:model%general%upn) = model%tempwk%inittemp(1:model%general%upn-1,ew,ns)

    ! basal boundary:
    ! for grounded ice, a heat flux is applied
    ! for floating ice, the basal temperature is held constant

!whl - This lower BC is different from the one in standard glide_temp.
!      If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
!       That is, net heat flux at lower boundary must equal zero.
!      If T(upn) >= Tpmp, then set T(upn) = Tpmp

    if (float) then

       supd(model%general%upn+1) = 0.0d0
       subd(model%general%upn+1) = 0.0d0
       diag(model%general%upn+1) = 1.0d0

       model%tempwk%inittemp(model%general%upn,ew,ns) = model%temper%temp(model%general%upn,ew,ns) 
       rhsd(model%general%upn+1) = model%temper%temp(model%general%upn,ew,ns)

    else    ! grounded ice

       call calcpmptb(pmptempb, model%geometry%thck(ew,ns))

       if (abs(model%temper%temp(model%general%upn,ew,ns) - pmptempb) < 0.001d0) then  ! melting

          ! hold basal temperature at pressure melting point

          supd(model%general%upn+1) = 0.0d0
          subd(model%general%upn+1) = 0.0d0
          diag(model%general%upn+1) = 1.0d0

          model%tempwk%inittemp(model%general%upn,ew,ns) = pmptempb
          rhsd(model%general%upn+1) = pmptempb

       else   ! frozen at bed
              ! maintain balance of heat sources and sinks
              ! (conductive flux, geothermal flux, and basal friction)
              ! Note: Heat fluxes are positive down, so slterm <= 0 and bheatflx <= 0.

          ! matrix elements corresponding to dT/dsigma
          ! 0.5 is a Crank-Nicolson factor

          subd(model%general%upn+1) = -0.5d0 / (1.0d0 - model%numerics%stagsigma(model%general%upn-1))
          supd(model%general%upn+1) =  0.0d0 
          diag(model%general%upn+1) = -subd(model%general%upn+1)

          ! Note: The heat source due to basal sliding (bfricflx) is computed in subroutine calcbfric.
          ! Also note that bheatflx is generally <= 0, since defined as positive down.

          model%tempwk%inittemp(model%general%upn,ew,ns) =    &
                - model%temper%temp(model%general%upn-1,ew,ns) * subd(model%general%upn+1)  &
                - model%temper%temp(model%general%upn,  ew,ns) * diag(model%general%upn+1)  &
                - model%geometry%thck(ew,ns)*thk0/coni * model%temper%bheatflx(ew,ns) & ! geothermal (H/k)*G
                + model%geometry%thck(ew,ns)*thk0/coni * model%temper%bfricflx(ew,ns)   ! sliding (H/k)*taub*ub

          rhsd(model%general%upn+1) = model%tempwk%inittemp(model%general%upn,ew,ns)

       endif   ! melting or frozen

    end if     ! floating or grounded

  end subroutine glissade_findvtri

  !-----------------------------------------------------------------------

  subroutine calcbfric (model,    whichbmlt,   &
                        thck,     btraction,   &
                        dusrfdew, dusrfdns,    &
                        ubas,     vbas,        &
                        float)

    ! compute frictional heat source due to sliding at the bed
    ! whl - to do - move to glide_temp_utils?

    use glimmer_global,   only: dp 
    use glimmer_physcon,  only: rhoi, grav
    use glimmer_paramets, only: thk0

    implicit none 

    type(glide_global_type) :: model
    integer, intent(in) :: whichbmlt
    real(dp), dimension(:,:), intent(in) :: thck, dusrfdew, dusrfdns
    real(dp), dimension(:,:), intent(in) :: ubas, vbas
    real(dp), dimension(:,:,:), intent(in) :: btraction
    logical, dimension(:,:), intent(in) :: float

    real(dp) :: slterm       ! sliding friction
 
    integer :: ewp, nsp, ew, ns
    integer :: slide_count   ! number of neighbor cells with nonzero sliding

       ! compute heat source due to basal friction
       ! Note: slterm and bfricflx are defined to be >= 0

       do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          slterm = 0.d0
          slide_count = 0

          select case( whichbmlt)

          case( SIA_BMELT)      ! taub*ub = -rhoi * g * H * (grad(S) * ubas) 

             if (thck(ew,ns) > model%numerics%thklim .and. .not. float(ew,ns)) then

                do nsp = ns-1,ns
                do ewp = ew-1,ew
                   if (abs(model%velocity%ubas(ewp,nsp)) > 1.0d-6 .or.  &
                       abs(model%velocity%vbas(ewp,nsp)) > 1.0d-6) then
                      slide_count = slide_count + 1
                      slterm = slterm - (&
                               model%geomderv%dusrfdew(ewp,nsp) * model%velocity%ubas(ewp,nsp) - &
                               model%geomderv%dusrfdns(ewp,nsp) * model%velocity%vbas(ewp,nsp))
                   end if
                end do
                end do

                slterm = slterm * rhoi * grav * model%geometry%thck(ew,ns)*thk0

             endif  ! thk > thklim, not floating

          case( FIRSTORDER_BMELT) 

             !whl - copied Steve Price's formulation from calcbmlt
             ! btraction is computed in glam_strs2.F90

             if (thck(ew,ns) > model%numerics%thklim .and. .not. float(ew,ns)) then
                do nsp = ns-1,ns
                do ewp = ew-1,ew
                   if (abs(model%velocity%ubas(ewp,nsp)) > 1.0d-6 .or.   &
                       abs(model%velocity%vbas(ewp,nsp)) > 1.0d-6) then
                      slide_count = slide_count + 1
                      slterm = slterm + model%velocity_hom%btraction(1,ewp,nsp) * &
                                        model%velocity_hom%uvel(model%general%upn,ewp,nsp) &
                                      + model%velocity_hom%btraction(2,ewp,nsp) * &
                                        model%velocity_hom%vvel(model%general%upn,ewp,nsp) 
                   end if
                end do
                end do

             endif  ! thk > thklim, not floating

          end select  ! whichbmlt

          ! include sliding contrib only if temperature node is surrounded by sliding velo nodes
          !whl - to do - This may result in non-conservation of energy.
          !              Why not include all nonzero terms? 

          if (slide_count >= 4) then
             slterm = 0.25d0 * slterm
          else
             slterm = 0.0d0
          end if

          model%temper%bfricflx(ew,ns) = slterm

       enddo    ! ns
       enddo    ! ew

  end subroutine calcbfric

  !-----------------------------------------------------------------------------------

  subroutine glissade_calcbmlt( model,    whichbmelt,    &
                                temp,     stagsigma,     &
                                thck,     stagthck,      &
                                dusrfdew, dusrfdns,      &
                                ubas,     vbas,          &
                                bmlt,     floater)

    ! Compute the amount of basal melting.
    ! Any temperatures above the pressure melting point are reset here to the
    !  pmp temperature, with excess energy applied toward basal melting.

    use glimmer_global, only : dp 
    use glimmer_physcon, only: shci, rhoi, lhci
    use glide_temp_utils, only: calcpmpt, calcpmptb

    implicit none 

    type(glide_global_type) :: model

    real(dp), dimension(0:,:,:), intent(inout) :: temp
    real(dp), dimension(0:),     intent(in) :: stagsigma
    real(dp), dimension(:,:),    intent(in) :: thck,  stagthck, dusrfdew, dusrfdns, ubas, vbas  
    real(dp), dimension(:,:),    intent(out):: bmlt    ! scaled melt rate (m/s * tim0/thk0)
    logical,  dimension(:,:),    intent(in) :: floater
    integer,  intent(in) ::      whichbmelt

    real(dp), dimension(model%general%upn) :: pmptemp   ! pressure melting point temperature
    real(dp) :: bflx    ! heat flux available for basal melting (W/m^2)
    real(dp) :: hmlt    ! scaled depth of internal melting (m/thk0)
    integer :: up, ew, ns

    bmlt(:,:) = 0.0d0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1

          if (thck(ew,ns) > model%numerics%thklim .and. .not. floater(ew,ns)) then

             ! Basal friction term is computed above in subroutine calcbfric

             ! Compute basal melting
             ! f(2) = tim0 / (thk0 * lhci * rhoi)

             bflx = model%temper%bfricflx(ew,ns) + model%temper%lcondflx(ew,ns) - model%temper%bheatflx(ew,ns)
             bmlt(ew,ns) = bflx * model%tempwk%f(2)

            ! Add internal melting associated with temp > pmptemp

            ! whl - to do - adjust layer thickness?
            ! If internal melting is rare, it should be OK to remove ice from the lowest layer only.
            ! But for temperate ice we should do something more realistic.

             call calcpmpt(pmptemp(:), thck(ew,ns), stagsigma(1:model%general%upn) )

             do up = 1, model%general%upn-1
                 if (temp(up,ew,ns) > pmptemp(up)) then
                    hmlt = (shci * thck(ew,ns) * (temp(up,ew,ns) - pmptemp(up))) / (rhoi * lhci) 
                    bmlt(ew,ns) = bmlt(ew,ns) + hmlt / model%numerics%dttem 
                    temp(up,ew,ns) = pmptemp(up)
                 endif
             enddo

             ! Reset basal temp to pmptemp, if necessary

             up = model%general%upn
             call calcpmptb(pmptemp(up), thck(ew,ns))
             temp(up,ew,ns) = min (temp(up,ew,ns), pmptemp(up))

          endif   ! thk > thklim

       enddo
    enddo

    ! apply periodic BC

    if (model%options%periodic_ew) then
       do ns = 2,model%general%nsn-1
          bmlt(1,ns) = bmlt(model%general%ewn-1,ns)
          bmlt(model%general%ewn,ns) = bmlt(2,ns)
       end do
    end if

  end subroutine glissade_calcbmlt

!-------------------------------------------------------------------
 
end module glissade_temp

