!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_therm.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
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

  ! This module combines two previous modules: glissade_temp.F90 and glissade_enthalpy.F90.
  ! It was created in Nov. 2014 by William Lipscomb (LANL).
  ! It gives the same answers as the two previous modules, but with modified data structures
  !  and streamlined organization.
  ! 
  ! The module computes temperature diffusion and strain heating in a local column 
  !  without doing horizontal or vertical advection.
  ! Temperature advection is done separately, e.g. using the incremental 
  !  remapping transport scheme.
  ! It is assumed here that temperature (and enthalpy) values are staggered in the
  !  vertical compared to the velocity.  That is, the temperature lives
  !  at the midpoint of each layer instead of at layer interfaces.
  ! Temperature and enthalpy are also defined at the upper and lower surfaces with 
  !  appropriate boundary conditions.   
  !
  ! The glissade_temp.F90 module was based on glide_temp.F90, with modifications for Glissade.
  ! The glissade_enthalpy.F90 module contained modifications for a new enthalpy option;
  !  it was originally written by Brian Macpherson (CU) under the supervision of Hari Rajaram.

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade_therm

    use glimmer_global, only : dp 
    use glide_types
    use glimmer_log
    use parallel, only: this_rank

    implicit none

    private
    public :: glissade_init_therm, glissade_therm_driver, glissade_therm_calcflwa, glissade_therm_calcbpmp, &
              glissade_enth2temp, glissade_temp2enth

    !WHL - debug
    integer :: itest, jtest, rtest

    ! time stepping scheme

    !NOTE:  For the dome test case, the Crank-Nicolson scheme can give unstable 
    !        temperature fluctuations for thin ice immediately after the ice 
    !        becomes thick enough for the temperature calculation.
    !       The fully implicit scheme has been stable for all cases (but is only
    !        first-order accurate in time). 

    logical, parameter::   &
         crank_nicolson = .false.  ! if true, use Crank-Nicolson time-stepping
                                   ! if false, use fully implicit

    ! max and min allowed temperatures
    ! Temperatures sometimes go below -100 for cases where Crank-Nicholson is unstable
    real(dp), parameter ::   &
       maxtemp_threshold = 1.d11,   &
       mintemp_threshold = -100.d0

    real(dp), dimension(:,:), allocatable :: dups   ! vertical grid quantities

  contains

!****************************************************    

  subroutine glissade_init_therm (model)

    ! initialization subroutine for the case that temperature lives on the
    ! vertically staggered grid (i.e., at layer centers)

    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow, trpt
    use glimmer_paramets, only : tim0, thk0, len0, vis0, vel0, tau0
    use parallel, only: lhalo, uhalo
!!    use glissade_enthalpy, only: glissade_init_enthalpy

    type(glide_global_type),intent(inout) :: model       !> Ice model parameters.

    integer, parameter :: p1 = gn + 1  
    integer up, ns, ew, upn, nsn, ewn

    upn = model%general%upn
    ewn = model%general%ewn
    nsn = model%general%nsn

    !TODO - Should these allocations be done in glide_allocarr?

    ! Note vertical dimensions here.  Dissipation is computed for each of (upn-1) layers.
    ! Temperature is defined at midpoint of each layer, plus upper and lower surfaces.
    allocate(dups(upn+1,2))
    dups(:,:) = 0.0d0

    !Note: The 'dups' grid coefficients are not the same as for unstaggered temperatures.
    up = 1
    dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                    (model%numerics%stagsigma(up) - model%numerics%sigma(up)) )
    do up = 2, model%general%upn-1
       dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                       (model%numerics%stagsigma(up) - model%numerics%stagsigma(up-1)) )
    enddo

    do up = 1, model%general%upn-2
       dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                       (model%numerics%stagsigma(up+1) - model%numerics%stagsigma(up)) )
    end do
    up = model%general%upn-1
    dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up)) * &
                                    (model%numerics%sigma(up+1) - model%numerics%stagsigma(up)) )

      !==== Initialize ice temperature.============
      !This block of code is similar to that in glide_init_temp

    ! Five possibilities:
    ! (1) Set ice temperature to 0 C everywhere in column (TEMP_INIT_ZERO)
    ! (2) Set ice temperature to surface air temperature everywhere in column (TEMP_INIT_ARTM)
    ! (3) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
    !     at the bed (TEMP_INIT_LINEAR). 
    !     A parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
    ! (4) Read ice temperature from an initial input file.
    ! (5) Read ice temperature from a restart file.
    !
    ! If restarting, we always do (5).
    ! If not restarting and the temperature field is present in the input file, we do (4).
    ! If (4) or (5), then the temperature field should already have been read from a file,
    !  and the rest of this subroutine will do nothing.
    ! Otherwise, the initial temperature is controlled by model%options%temp_init,
    !  which can be read from the config file.
    !
    !TODO - When reading temperature from restart or input file, make sure that halo values are correct.

    if (model%options%is_restart == RESTART_TRUE) then

       ! Temperature has already been initialized from a restart file. 
       ! (Temperature is always a restart variable.)

       call write_log('Initializing ice temperature from the restart file')

    elseif ( minval(model%temper%temp(1:model%general%upn, &
                    1+lhalo:model%general%ewn-lhalo, 1+uhalo:model%general%nsn-uhalo)) > &
                    (-1.0d0 * trpt) ) then    ! trpt = 273.15 K
                                              ! Default initial temps in glide_types are -999

       if ( (maxval(model%temper%temp(model%general%upn+1, &
                    1+lhalo:model%general%ewn-lhalo, 1+uhalo:model%general%nsn-uhalo)) == &
                    -999.0d0 ) .and.    &
                    (model%options%whichdycore /= DYCORE_BISICLES) )then
          ! Throw a fatal error if we think the user has supplied temp instead of tempstag
          ! (We don't want to implicitly shift the vertical layers from one coordinate system
          !  to another without the user knowing.)
          ! This case will look like good data in all the layers except the top layer.
          ! MJH: Letting BISICLES run with this situation as per Dan Martin's request.
          call write_log("The variable 'temp' has been read from an input file, but it only is appropriate " &
             // "for the GLIDE dycore.  Use the 'tempstag' variable with higher-order dycores instead.", GM_FATAL)
       else
          ! Temperature has already been initialized from an input file.
          ! (We know this because the default initial temps of -999 have been overwritten.)

          call write_log('Initializing ice temperature from an input file')
       endif

    else   ! not reading temperature from restart or input file
           ! initialize it here basee on model%options%temp_init

       ! First set T = 0 C everywhere

       model%temper%temp(:,:,:) = 0.0d0                                              
                                                    
       if (model%options%temp_init == TEMP_INIT_ZERO) then

          call write_log('Initializing ice temperature to 0 deg C')

          ! No call is needed to glissade_init_temp_column because the
          ! ice temperature has been set to zero above

       elseif (model%options%temp_init == TEMP_INIT_ARTM) then

          ! Initialize ice column temperature to min(artm, 0 C).

          !Note: Old glide sets temp = artm everywhere without regard to whether ice exists in a column.

          call write_log('Initializing ice temperature to the surface air temperature')

          do ns = 1, model%general%nsn
             do ew = 1, model%general%ewn

                call glissade_init_temp_column(model%options%temp_init,         &
                                               model%numerics%stagsigma(:),     &
                                               dble(model%climate%artm(ew,ns)), &
                                               model%geometry%thck(ew,ns),      &
                                               model%temper%temp(:,ew,ns) )
             end do
          end do

       elseif (model%options%temp_init == TEMP_INIT_LINEAR) then

          ! Initialize ice column temperature with a linear profile:
          ! T = artm at the surface, and T <= Tpmp at the bed.
          ! Loop over physical cells where artm is defined (not temperature halo cells)

          call write_log('Initializing ice temperature to a linear profile in each column')

          do ns = 1, model%general%nsn
             do ew = 1, model%general%ewn

                call glissade_init_temp_column(model%options%temp_init,         &
                                               model%numerics%stagsigma(:),     &
                                               dble(model%climate%artm(ew,ns)), &
                                               model%geometry%thck(ew,ns),      &
                                               model%temper%temp(:,ew,ns) )

             end do
          end do

       endif ! model%options%temp_init

    endif    ! restart file, input file, or other options

  end subroutine glissade_init_therm

!=======================================================================

  subroutine glissade_init_temp_column(temp_init,                 &
                                       stagsigma,   artm,         &
                                       thck,        temp)

  ! Initialize temperatures in a column based on the value of temp_init.
  ! Three possibilities:
  ! (1) Set ice temperature in column to 0 C (TEMP_INIT_ZERO)
  ! (2) Set ice temperature in column to surface air temperature (TEMP_INIT_ARTM)
  ! (3) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
  !     at the bed (TEMP_INIT_LINEAR). 
  !     A local parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
  !
  ! This subroutine is functionally equivalent to glide_init_temp_column.
  ! The only difference is that temperature is staggered in the vertical
  !  (i.e., located at layer midpoints as well as the top and bottom surfaces).

  ! In/out arguments
 
  integer, intent(in) :: temp_init          ! option for temperature initialization

  real(dp), dimension(:), intent(in)     :: stagsigma  ! staggered vertical coordinate
                                                       ! includes layer midpoints, but not top and bottom surfaces
  real(dp), intent(in)                   :: artm   ! surface air temperature (deg C)
                                                   ! Note: artm should be passed in as double precision
  real(dp), intent(in)                   :: thck   ! ice thickness
  real(dp), dimension(0:), intent(inout) :: temp   ! ice column temperature (deg C)
                                                   ! Note first index of zero
                                                   
  ! Local variables and parameters

  real(dp) :: pmptb                              ! pressure melting point temp at the bed
  real(dp), dimension(size(stagsigma)) :: pmpt   ! pressure melting point temp thru the column
  integer :: upn                                 ! number of vertical levels (deduced from temp array)

  !TODO - Define pmpt_offset elsewhere?
  real(dp), parameter :: pmpt_offset = 2.d0  ! offset of initial Tbed from pressure melting point temperature (deg C)
                                             ! Note: pmtp_offset is positive for T < Tpmp

  upn = size(temp) - 1     ! temperature array has dimension (0:model%general%upn)

  ! Set the temperature in the column

  select case(temp_init)

  case(TEMP_INIT_ZERO)     ! set T = 0 C

     temp(:) = 0.d0

  case(TEMP_INIT_ARTM)     ! initialize ice-covered areas to the min of artm and 0 C
                           ! set ice-free areas to T = 0 C

     if (thck > 0.0d0) then
        temp(:) = min(0.0d0, artm)
     else
        temp(:) = 0.d0
     endif

  case(TEMP_INIT_LINEAR)

     ! Tsfc = artm, Tbed = Tpmp = pmpt_offset, linear profile in between 

     temp(0) = artm

     call glissade_calcpmpt_bed (pmptb, thck)
     temp(upn) = pmptb - pmpt_offset

     temp(1:upn-1) = temp(0) + (temp(upn) - temp(0))*stagsigma(:)
                               
     ! Make sure T <= Tpmp - pmpt_offset in column interior
     ! TODO: Change condition to T <= Tpmp?

     call glissade_calcpmpt(pmpt(:), thck, stagsigma(:))
     temp(1:upn-1) = min(temp(1:upn-1), pmpt(1:upn-1) - pmpt_offset)

  end select

  end subroutine glissade_init_temp_column

!=======================================================================

  subroutine glissade_therm_driver(model, whichtemp)

    ! Calculates the new ice temperature 

    use glimmer_utils,  only : tridiag
    use glimmer_paramets, only : thk0, tim0
    use glimmer_physcon, only: shci, coni, rhoi
    use glide_mask
!!    use glissade_enthalpy
    use glissade_grid_operators, only: glissade_stagger
    use glissade_masks, only: glissade_get_masks

    !TODO - Modify glissade_temp_driver to compute over locally owned cells only?
    !       This would make the module a bit cheaper but would require halo updates at the end.

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model       ! Ice model parameters
    integer,                intent(in)    :: whichtemp   ! Flag to choose method

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    integer :: ew, ns, up
    integer :: ewn, nsn, upn
    character(len=100) :: message

    ! These arrays have the same size as the vertical dimension of temperature:
    !  upn+1 on the staggered grid.
    real(dp), dimension(size(model%temper%temp,1)) :: subd, diag, supd, rhsd

    ! These have the same dimensions as staggered temperature
    real(dp),dimension(0:model%general%upn) :: prevtemp_stag, enthalpy

    ! for energy conservation check
    real(dp) :: einit, efinal, delta_e, dTtop, dTbot, denth_top, denth_bot

    real(dp) :: maxtemp, mintemp   ! max and min temps in column

    real(dp), dimension(0:model%general%upn) :: pmptemp   ! pressure melting pt temperature

    real(dp), dimension(1:model%general%upn) :: alpha_enth   ! diffusivity at interfaces (m2/s) for enthalpy solver
                                                             ! = coni / (rhoi*shci) for cold ice

    integer, dimension(model%general%ewn,model%general%nsn) ::  &
         ice_mask,      &! = 1 where ice velocity is computed (thck > thklim), else = 0
         ice_mask_temp, &! = 1 where ice temperature is computed (thck > thklim_temp), else = 0
         floating_mask   ! = 1 where ice is floating, else = 0

    logical, parameter:: verbose_temp = .false.
    logical :: verbose_column
    integer :: k

    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    select case(whichtemp)

    case(TEMP_SURFACE_AIR_TEMP)  ! Set column to surface air temperature ------------------

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns)))
          end do
       end do

    case(TEMP_STEADY)   ! do nothing

    case(TEMP_PROGNOSTIC, TEMP_ENTHALPY)  ! Local column calculation (with advection done elsewhere)

       print*, 'Temp driver, whichtemp =', whichtemp

       ! Compute masks: ice_mask = 1 where thck > thklim, floating_mask = 1 where ice is floating

       call glissade_get_masks(ewn,                  nsn,                    &
                               model%geometry%thck,  model%geometry%topg,    &
                               model%climate%eus,    model%numerics%thklim,  &
                               ice_mask,             floating_mask)

       ! Compute ice mask for temperature: ice_mask_temp = 1 where thck > thklim_temp

       call glissade_get_masks(ewn,                  nsn,                         &
                               model%geometry%thck,  model%geometry%topg,         &
                               model%climate%eus,    model%numerics%thklim_temp,  &
                               ice_mask_temp)

       ! No horizontal or vertical advection; vertical diffusion and strain heating only.
       ! Temperatures are vertically staggered relative to velocities.  
       ! That is, the temperature is defined at the midpoint of each layer 
       ! (and at the top and bottom surfaces).
       
       ! Calculate interior heat dissipation -------------------------------------

       model%temper%dissip(:,:,:) = 0.d0

       if (model%options%which_ho_disp == HO_DISP_SIA) then

          call glissade_interior_dissipation_sia(model%numerics%dttem,       &
                                                 ewn,          nsn,          &
                                                 upn,          model%numerics%stagsigma(:), &
                                                 ice_mask,                   &
                                                 model%geomderv%stagthck,    &
                                                 model%temper%flwa,          &
                                                 model%geomderv%dusrfdew,    &
                                                 model%geomderv%dusrfdns,   &
                                                 model%temper%dissip)

       else    ! first-order dissipation

          call glissade_interior_dissipation_first_order(model%numerics%dttem,       &
                                                         ewn,          nsn,          &
                                                         upn,                        &
                                                         ice_mask,                   &
                                                         model%stress%tau%scalar,    &
                                                         model%stress%efvs,          &
                                                         model%temper%dissip)
          
       endif

       ! Calculate heating from basal friction
       ! (only for Glam; bfricflx is computed by the Glissade velocity solver)

       if (model%options%whichdycore == DYCORE_GLAM) then
          call glissade_basal_friction(ewn,         nsn,                &
                                       ice_mask,    floating_mask,      &
                                       model%velocity%uvel(upn,:,:),    &
                                       model%velocity%vvel(upn,:,:),    &
                                       model%velocity%btraction(:,:,:), &
                                       model%temper%bfricflx(:,:) )
       endif

       do ns = 2, nsn-1
       do ew = 2, ewn-1

          if (verbose_temp .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
             verbose_column = .true.
          else
             verbose_column = .false.
          endif

          if (ice_mask_temp(ew,ns) == 1) then

             ! Set surface temperature

             model%temper%temp(0,ew,ns) = dble(model%climate%artm(ew,ns))

             if (whichtemp == TEMP_ENTHALPY) then

                ! Given temperature and waterfrac, compute enthalpy (dimension 0:upn).
                ! Assume waterfrac = 0 at upper and lower surfaces.

                call glissade_temp2enth(model%temper%enthalpy(0:upn,ew,ns),    &
                                        model%temper%temp(0:upn,ew,ns),        &
                                        model%temper%waterfrac(1:upn-1,ew,ns), &
                                        model%geometry%thck(ew,ns),            &
                                        model%numerics%stagsigma(1:upn-1))

                if (verbose_column) then
                   print*, ' '
                   print*, 'Before prognostic enthalpy, i, j =', ew, ns
                   print*, 'thck =', model%geometry%thck(ew,ns)*thk0
                   print*, 'Temp, waterfrac, enthalpy:'
                   k = 0
                   print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
                   do k = 1, upn-1
                      print*, k, model%temper%temp(k,ew,ns), model%temper%waterfrac(k,ew,ns), model%temper%enthalpy(k,ew,ns)
                   enddo
                   k = upn
                   print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
                endif

                ! compute initial internal energy in column (for energy conservation check)
                einit = 0.0d0
                do up = 1, upn-1
                   einit = einit + model%temper%enthalpy(up,ew,ns) *  &
                                  (model%numerics%sigma(up+1) - model%numerics%sigma(up) )
                enddo
                einit = einit * model%geometry%thck(ew,ns)*thk0

                ! compute matrix elements using enthalpy gradient method

                call glissade_enthalpy_findvtri(model%numerics%dttem,            &
                                                upn,                             &
                                                model%numerics%stagsigma,        &
                                                subd,  diag, supd, rhsd,         &
                                                dups(:,:),                       &
                                                floating_mask(ew,ns),            &
                                                model%geometry%thck(ew,ns),      &
                                                model%temper%temp(:,ew,ns),      &
                                                model%temper%waterfrac(:,ew,ns), &
                                                model%temper%enthalpy(:,ew,ns),  &
                                                model%temper%dissip(:,ew,ns),    &
                                                model%temper%bheatflx(ew,ns),    &
                                                model%temper%bfricflx(ew,ns),    &
                                                alpha_enth,                      &
                                                verbose_column)

                !WHL - debug
                if (verbose_column) then
                   print*, ' '
                   print*, 'After vtri, i, j =', ew, ns
                   print*, 'k, subd, diag, supd, rhs/(rhoi*ci):'
                   do k = 1, upn+1
                      print*, k-1, subd(k), diag(k), supd(k), rhsd(k)/(rhoi*shci)
                   enddo
                endif

                ! solve the tridiagonal system
                ! Note: Enthalpy is indexed from 0 to upn, with indices 1 to upn-1 colocated
                ! with stagsigma values of the same index.
                ! However, the matrix elements are indexed 1 to upn+1, with the first row
                ! corresponding to the surface enthalpy, enthalpy(0).

                call tridiag(subd(1:upn+1),   &
                             diag(1:upn+1),   &
                             supd(1:upn+1),   &
                             enthalpy(0:upn), &
                             rhsd(1:upn+1))
               
                ! Copy the local enthalpy array into the global derived type
                model%temper%enthalpy(:,ew,ns) = enthalpy(:)

                ! BDM convert back to temperature and water content
                call glissade_enth2temp(model%temper%enthalpy(0:upn,ew,ns),    &
                                        model%temper%temp(0:upn,ew,ns),        &
                                        model%temper%waterfrac(1:upn-1,ew,ns), &
                                        model%geometry%thck(ew,ns),            &
                                        model%numerics%stagsigma(1:upn-1))

                if (verbose_column) then
                   print*, ' '
                   print*, 'After prognostic enthalpy, i, j =', ew, ns
                   print*, 'thck =', model%geometry%thck(ew,ns)*thk0
                   print*, 'Temp, waterfrac, enthalpy:'
                   k = 0
                   print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
                   do k = 1, upn-1
                      print*, k, model%temper%temp(k,ew,ns), model%temper%waterfrac(k,ew,ns), model%temper%enthalpy(k,ew,ns)
                   enddo
                   k = upn
                   print*, model%temper%temp(k,ew,ns), 0.d0, model%temper%enthalpy(k,ew,ns)
                endif
                
                ! compute the final internal energy

                efinal = 0.0d0
                do up = 1, upn-1
                   efinal = efinal + model%temper%enthalpy(up,ew,ns) *  &
                                    (model%numerics%sigma(up+1) - model%numerics%sigma(up) )
                enddo
                efinal = efinal * model%geometry%thck(ew,ns)*thk0

                ! compute net heat flux to the column
                
                ! conductive flux = (alpha/H * denth/dsigma) at upper and lower surfaces; positive down.
                ! Here alpha = coni / (rhoi*shci) for cold ice, with a smaller value for temperate ice.
                ! Assume fully implicit backward Euler time step
                
                denth_top = enthalpy(1) - enthalpy(0)
                denth_bot = enthalpy(upn) - enthalpy(upn-1)
                
                model%temper%ucondflx(ew,ns) = -alpha_enth(1) / (model%geometry%thck(ew,ns)*thk0)         &
                                              * denth_top / ( model%numerics%stagsigma(1) )

                model%temper%lcondflx(ew,ns) = -alpha_enth(upn) / (model%geometry%thck(ew,ns)*thk0)         &
                                              * denth_bot / (1.d0 - model%numerics%stagsigma(upn-1))

             else   ! whichtemp = TEMP_PROGNOSTIC

                if (verbose_column) then
                   print*, ' '
                   print*, 'Before prognostic temp, i, j =', ew, ns
                   print*, 'thck =', model%geometry%thck(ew,ns)*thk0
                   print*, 'Temp:'
                   do k = 0, upn
                      print*, k, model%temper%temp(k,ew,ns)
                   enddo
                endif
                
                ! compute initial internal energy in column (for energy conservation check)
                einit = 0.0d0
                do up = 1, upn-1
                   einit = einit + model%temper%temp(up,ew,ns) *  &
                        (model%numerics%sigma(up+1) -   &
                        model%numerics%sigma(up) )
                enddo
                einit = einit * rhoi * shci * model%geometry%thck(ew,ns)*thk0
                
                ! compute matrix elements

                call glissade_findvtri(model%numerics%dttem,            &
                                       upn,                             &
                                       model%numerics%stagsigma,        &
                                       subd,  diag, supd, rhsd,         &            
                                       floating_mask(ew,ns),            &
                                       model%geometry%thck(ew,ns),      &
                                       model%temper%temp(:,ew,ns),      &
                                       model%temper%dissip(:,ew,ns),    &
                                       model%temper%bheatflx(ew,ns),    &
                                       model%temper%bfricflx(ew,ns))
                
                if (verbose_column) then
                   print*, 'After glissade_findvtri, i, j =', ew,ns
                   print*, 'k, subd, diag, supd, rhsd:'
                   do k = 1, upn+1
                      print*, k, subd(k), diag(k), supd(k), rhsd(k)
                   enddo
                endif

                prevtemp_stag(:) = model%temper%temp(:,ew,ns)

                ! solve the tridiagonal system

                ! Note: Temperature is indexed from 0 to upn, with indices 1 to upn-1 colocated
                !  with stagsigma values of the same index.
                ! However, the matrix elements are indexed 1 to upn+1, with the first row
                !  corresponding to the surface temperature, temp(0,:,:).

                call tridiag(subd(1:upn+1), &
                             diag(1:upn+1), &
                             supd(1:upn+1), &
                             model%temper%temp(0:upn,ew,ns), &
                             rhsd(1:upn+1))

                if (verbose_column) then
                   print*, ' '
                   print*, 'After prognostic temp, i, j =', ew, ns
                   print*, 'Temp:'
                   do k = 0, upn
                      print*, k, model%temper%temp(k,ew,ns)
                   enddo
                endif
                
                ! compute the final internal energy
                
                efinal = 0.0d0
                do up = 1, upn-1
                   efinal = efinal + model%temper%temp(up,ew,ns) *  &
                        (model%numerics%sigma(up+1) - model%numerics%sigma(up))
                enddo
                efinal = efinal * rhoi*shci * model%geometry%thck(ew,ns)*thk0
                
                ! compute net heat flux to the column

                ! conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down

                if (crank_nicolson) then
                   ! average temperatures between start and end of timestep
                   dTtop = 0.5d0 * ( model%temper%temp(1,ew,ns) - model%temper%temp(0,ew,ns) &
                                   +     prevtemp_stag(1)       -     prevtemp_stag(0) )
                   dTbot = 0.5d0 * ( model%temper%temp(upn,ew,ns) - model%temper%temp(upn-1,ew,ns) &
                                   +     prevtemp_stag(upn)       -   prevtemp_stag(upn-1) )
                else    ! fully implicit
                   ! use temperatures at end of timestep
                   dTtop = model%temper%temp(1,ew,ns) - model%temper%temp(0,ew,ns)
                   dTbot = model%temper%temp(upn,ew,ns) - model%temper%temp(upn-1,ew,ns)
                endif

                model%temper%ucondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                              * dTtop / (model%numerics%stagsigma(1))

                model%temper%lcondflx(ew,ns) = (-coni / (model%geometry%thck(ew,ns)*thk0) )         &
                                              * dTbot / (1.d0 - model%numerics%stagsigma(upn-1))

             endif   ! whichtemp

             ! Compute total dissipation in column (W/m^2)

             model%temper%dissipcol(ew,ns) = 0.0d0
             do up = 1, upn-1
                model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew,ns) + &
                                                model%temper%dissip(up,ew,ns)  &
                                             * (model%numerics%sigma(up+1) - model%numerics%sigma(up))  
             enddo
             model%temper%dissipcol(ew,ns) = model%temper%dissipcol(ew, ns)     &
                                           * thk0*model%geometry%thck(ew,ns)*rhoi*shci / (tim0*model%numerics%dttem)  

             ! Verify that the net input of energy into the column is equal to the change in
             ! internal energy.  

             delta_e = (model%temper%ucondflx(ew,ns) - model%temper%lcondflx(ew,ns)  &
                      + model%temper%dissipcol(ew,ns)) * tim0*model%numerics%dttem

             if ( abs((efinal-einit-delta_e)/(tim0*model%numerics%dttem)) > 1.0d-8 ) then

                if (verbose_column) then
                   print*, 'Ice thickness:', thk0*model%geometry%thck(ew,ns)
                   print*, 'thklim_temp:', thk0*model%numerics%thklim_temp
                   print*, ' '
                   print*, 'Interior fluxes:'
                   print*, 'ftop (pos up)=', -model%temper%ucondflx(ew,ns) 
                   print*, 'fbot (pos up)=', -model%temper%lcondflx(ew,ns)
                   print*, 'fdissip =',       model%temper%dissipcol(ew,ns)
                   print*, 'Net flux =', delta_e/(tim0*model%numerics%dttem)
                   print*, ' '
                   print*, 'delta_e =', delta_e
                   print*, 'einit =',  einit
                   print*, 'efinal =', efinal
                   print*, 'einit + delta_e =', einit + delta_e
                   print*, ' '
                   print*, 'Energy imbalance =', efinal - einit - delta_e
                   print*, ' '
                   print*, 'Basal fluxes:'
                   print*, 'ffric =', model%temper%bfricflx(ew,ns)
                   print*, 'fgeo =', -model%temper%bheatflx(ew,ns)
                   print*, 'flux for bottom melting =', model%temper%bfricflx(ew,ns)   &
                                                      - model%temper%bheatflx(ew,ns)   &
                                                      + model%temper%lcondflx(ew,ns)
                endif   ! verbose_temp

                write(message,*) 'WARNING: Energy conservation error, ew, ns =', ew, ns
                call write_log(message,GM_FATAL)
             endif

             !WHL - No call here to corrpmpt.  Temperatures above pmpt are set to pmpt 
             !      in glissade_calcbmlt (conserving energy).

          endif  ! thck > thklim_temp
       end do    ! ew
       end do    ! ns

       ! Set temperature of thin ice to the air temperature and set ice-free nodes to zero

       do ns = 1, model%general%nsn
          do ew = 1, model%general%ewn

              if (model%geometry%thck(ew,ns) <= model%numerics%thklim_temp) then
                 model%temper%temp(:,ew,ns) = min(0.d0, model%climate%artm(ew,ns))
              endif

              !TODO - Maybe it should be done in the following way, so that the temperature profile for thin ice
              !       is consistent with the temp_init option, with T = 0 for ice-free cells.

             ! NOTE: Calling this subroutine will maintain a sensible temperature profile
             !        for thin ice, but in general does *not* conserve energy.
             !       To conserve energy, we need either thklim_temp = 0, or some additional
             !        energy accounting and correction.
 
!             if (model%geometry%thck(ew,ns) <= model%numerics%thklim_temp) then
!                call glissade_init_temp_column(model%options%temp_init,         &
!                                               model%numerics%stagsigma(:),     &
!                                               dble(model%climate%artm(ew,ns)), &
!                                               model%geometry%thck(ew,ns),      &
!                                               model%temper%temp(:,ew,ns) )
!             else if (model%geometry%thkmask(ew,ns) < 0) then
!                model%temper%temp(:,ew,ns) = 0.d0
!             end if

          end do
       end do

       ! Calculate basal melt rate
       ! Temperature above the pressure melting point are reset to Tpmp,
       !  with excess heat contributing to melting.

       call glissade_calcbmlt(whichtemp,                     &
                              model%numerics%dttem,          &
                              ewn,  nsn,  upn,               &
                              model%numerics%sigma,          &
                              model%numerics%stagwbndsigma,  &
                              ice_mask_temp,                 &
                              floating_mask,                 &
                              model%geometry%thck,           &
                              model%temper%temp,             &
                              model%temper%waterfrac,        &
                              model%temper%bfricflx,         &
                              model%temper%bheatflx,         &
                              model%temper%lcondflx,         &
                              model%temper%bwat,             &
                              model%temper%bmlt)

       ! Interpolate basal temperature and pressure melting point onto velocity grid

       call glissade_stagger(model%general%ewn,     model%general%nsn,  &
                             model%temper%temp(upn,:,:),  &
                             model%temper%stagbtemp(:,:),               &
                             ice_mask_temp(:,:),                        &
                             stagger_margin_in = 1)

       call glissade_stagger(model%general%ewn,     model%general%nsn,  &
                             model%temper%bpmp(:,:),                    &
                             model%temper%stagbpmp(:,:),                &
                             ice_mask_temp(:,:),                        &
                             stagger_margin_in = 1)

    end select

    ! Check for temperatures that are physically unrealistic.
    ! Thresholds are set at the top of this module.

    do ns = 1, model%general%nsn
       do ew = 1, model%general%ewn

          maxtemp = maxval(model%temper%temp(:,ew,ns))
          mintemp = minval(model%temper%temp(:,ew,ns))
          
          if (maxtemp > maxtemp_threshold) then
             write(message,*) 'maxtemp > 0: i, j, maxtemp =', ew, ns, maxtemp
             call write_log(message,GM_FATAL)
          endif
          
          if (mintemp < mintemp_threshold) then
             !uncommment these lines to get more info
!             print*, 'thck =', model%geometry%thck(ew,ns) * thk0
!             print*, 'temp:'
!             do k = 1, upn
!                print*, k, model%temper%temp(k,ew,ns)
!             enddo
             write(message,*) 'mintemp < mintemp_threshold: i, j, mintemp =', ew, ns, mintemp
             call write_log(message,GM_FATAL)
          endif
          
       enddo
    enddo

    ! Rescale dissipation term to deg C/s (instead of deg C)
    !TODO - Treat dissip above as a rate (deg C/s) instead of deg C?
    model%temper%dissip(:,:,:) =  model%temper%dissip(:,:,:) /  (model%numerics%dttem*tim0)

  end subroutine glissade_therm_driver

!=======================================================================
 
  subroutine glissade_interior_dissipation_sia(dttem,                &
                                               ewn,       nsn,       &
                                               upn,       stagsigma, &
                                               ice_mask,             &
                                               stagthck,  flwa,      &
                                               dusrfdew,  dusrfdns,  &
                                               dissip)

    ! Compute the dissipation source term associated with strain heating,
    ! based on the shallow-ice approximation.
    
    use glimmer_physcon, only : gn   ! Glen's n
    use glimmer_paramets, only: len0, tim0, thk0, vis0

    real(dp), intent(in) :: dttem          ! time step (model time units)
    integer, intent(in) :: ewn, nsn, upn   ! grid dimensions

    real(dp), dimension(upn-1), intent(in) :: stagsigma   ! staggered vertical grid for temperature

    real(dp), dimension(:,:), intent(in) :: stagthck, dusrfdew, dusrfdns

    integer, dimension(:,:), intent(in) :: ice_mask    ! = 1 where ice is present (thck > thklim), else = 0

    real(dp), dimension(:,:,:), intent(in) ::  &
         flwa         ! flow factor, Pa^(-n) yr^(-1)  

    real(dp), dimension(:,:,:), intent(out) ::  &
         dissip       ! interior heat dissipation (deg C during timestep)
    
    integer, parameter :: p1 = gn + 1  

    integer :: ew, ns
    real(dp), dimension(upn-1) :: sia_dissip_fact  ! factor in SIA dissipation calculation
    real(dp) :: geom_fact         ! geometric factor

    ! Two methods of doing this calculation: 
    ! 1. find dissipation at u-pts and then average
    ! 2. find dissipation at H-pts by averaging quantities from u-pts
    ! (2) works best for eismint divide (symmetry) but (1) may be better for full expts
    ! This subroutine uses (2).

    if (size(dissip,1) /= upn-1) then  ! staggered vertical grid
        !TODO - Write an error message and exit gracefully
    endif

    dissip(:,:,:) = 0.0d0

    ! Note: Factor of 16 is for averaging flwa
    sia_dissip_fact(1:upn-1) = (stagsigma(1:upn-1) * rhoi * grav * thk0**2 / len0)**p1 &
                             * 2.0d0 * vis0 * dttem * tim0 / (16.0d0 * rhoi * shci)

    do ns = 2, nsn-1
       do ew = 2, ewn-1
          !Note: ice_mask = 1 where thck > thcklim.  Elsewhere, dissipation is assumed to be zero.
          if (ice_mask(ew,ns) == 1) then
             geom_fact = (0.25d0*sum(stagthck(ew-1:ew,ns-1:ns)) * dsqrt((0.25d0*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                       + (0.25d0*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             dissip(:,ew,ns) = geom_fact * sia_dissip_fact *   & 
                              (flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                              2.d0*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                              4.d0*flwa(:,ew,ns))
          end if
       end do
    end do

  end subroutine glissade_interior_dissipation_sia

!=======================================================================
 
  subroutine glissade_interior_dissipation_first_order(dttem,                &
                                                       ewn,       nsn,       &
                                                       upn,                  &
                                                       ice_mask,             &
                                                       tau_eff,   efvs,      &
                                                       dissip)

    ! Compute the dissipation source term associated with strain heating.
    ! Note that the dissipation is computed in the same way on either a staggered or an
    !  unstaggered vertical grid.  
    ! Note also that dissip and flwa must have the same vertical dimension 
    !  (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
    
    use glimmer_paramets, only: tau0, len0, tim0, vel0

    real(dp), intent(in) :: dttem          ! time step (model time units)
    integer, intent(in) :: ewn, nsn, upn   ! grid dimensions
    integer, dimension(:,:), intent(in) :: ice_mask    ! = 1 where ice is present (thck > thklim), else = 0

    real(dp), dimension(:,:,:), intent(in) ::  &
         tau_eff,    & ! effective stress, Pa
         efvs          ! effective viscosity, Pa yr

    real(dp), dimension(:,:,:), intent(out) ::  &
         dissip       ! interior heat dissipation (deg C during timestep)
    
    integer :: ew, ns, k
    real(dp) :: ho_dissip_fact    ! factor in higher-order dissipation calculation

    ! 3D, 1st-order case
    ! Note: Glissade computes efvs and tau%scalar using the strain rate terms appropriate for the approximation.
    ! E.g, the SIA quantities are computed based on (du_dz, dv_dz) only, and the SSA quantities
    !  are computed based on (du_dx, du_dy, dv_dx, dv_dy) only.
    ! So this computation should give the appropriate heating for whichapprox = HO_APPROX_SIA,
    !  HO_APPROX_SSA, HO_APPROX_L1L2 or HO_APPROX_BP.
    !
    if (size(dissip,1) /= upn-1) then  ! staggered vertical grid
        !TODO - Write an error message and exit gracefully
    endif

    dissip(:,:,:) = 0.0d0
    ho_dissip_fact = (tau0*vel0/len0)/(rhoi*shci) * (dttem*tim0)

    do ns = 1, nsn
       do ew = 1, ewn
          !Note: ice_mask = 1 where thck > thcklim.  Elsewhere, dissipation is assumed to be zero.
          if (ice_mask(ew,ns) == 1) then
             do k = 1, upn-1
                if (efvs(k,ew,ns) /= 0.d0) then
                   dissip(k,ew,ns) = (tau_eff(k,ew,ns)**2 / efvs(k,ew,ns)) * ho_dissip_fact
                endif
             enddo
          endif
       enddo
    enddo

  end subroutine glissade_interior_dissipation_first_order

!=======================================================================

  subroutine glissade_basal_friction (ewn,        nsn,             &
                                      ice_mask,   floating_mask,   &
                                      ubas,       vbas,            &
                                      btraction,  bfricflx)

    ! Compute frictional heat source due to sliding at the bed
    ! Note: This subroutine is for Glam only; it is not scientifically supported.

    use glimmer_paramets, only: vel0, vel_scale

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: ewn, nsn         ! grid dimensions
    integer, dimension(:,:), intent(in) ::   &
         ice_mask,      & ! = 1 if thck > thklim, else = 0
         floating_mask    ! = 1 if ice is floating, else = 0
    real(dp), dimension(:,:), intent(in) :: ubas, vbas   ! basal velocity
    real(dp), dimension(:,:,:), intent(in) :: btraction  ! basal traction
    real(dp), dimension(:,:), intent(out) :: bfricflx    ! basal friction heat flux (W m-2)

    !-----------------------------------------------------------------
    ! Local arguments
    !-----------------------------------------------------------------

    real(dp) :: slterm       ! sliding friction
    integer :: ew, ns, i, j
    integer :: slide_count   ! number of neighbor cells with nonzero sliding

    bfricflx(:,:) = 0.d0

    ! compute heat source due to basal friction
    ! Note: slterm and bfricflx are defined to be >= 0

    do ns = 2, nsn-1
       do ew = 2, ewn-1

          slterm = 0.d0
          slide_count = 0

          ! Note: btraction is computed in glam_strs2.F90

          !WHL - Using thklim instead of thklim_temp because ice thinner than thklim
          !      is assumed to be at rest.

          if (ice_mask(ew,ns)==1 .and. floating_mask(ew,ns)==0) then
             do j = ns-1,ns
                do i = ew-1,ew

                   !SCALING - WHL: Multiplied ubas by vel0/vel_scale so we get the same result in these two cases:
                   !           (1) With scaling:     vel0 = vel_scale = 500/scyr, and ubas is non-dimensional
                   !           (2) Without scaling:  vel0 = 1, vel_scale = 500/scyr, and ubas is in m/s.

!!!                   if (abs(ubas(i,j)) > 1.0d-6 .or.   &
!!!                       abs(vbas(i,j)) > 1.0d-6) then
                   if ( abs(ubas(i,j))*(vel0/vel_scale) > 1.0d-6 .or.   &
                        abs(vbas(i,j))*(vel0/vel_scale) > 1.0d-6 ) then
                      slide_count = slide_count + 1
                      slterm = slterm + btraction(1,i,j)*ubas(i,j) + btraction(2,i,j)*vbas(i,j) 
                   end if
                end do
             end do

          endif  ! ice_mask = 1, floating_mask = 0

          ! include sliding contrib only if temperature node is surrounded by sliding velo nodes
          !NOTE - This may result in non-conservation of energy. 

          if (slide_count == 4) then
             slterm = 0.25d0 * slterm
          else
             slterm = 0.0d0
          end if

          bfricflx(ew,ns) = slterm

       enddo    ! ns
    enddo       ! ew

  end subroutine glissade_basal_friction

!=======================================================================

  subroutine glissade_findvtri (dttem,                       &
                                upn,   stagsigma,            &
                                subd,  diag, supd, rhsd,     &
                                floating_mask,               &
                                thck,        temp,           &
                                dissip,                      &
                                bheatflx,  bfricflx)

    ! compute matrix elements for the tridiagonal solve

    use glimmer_paramets, only : thk0, tim0
    use glimmer_physcon,  only : rhoi, grav, coni

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    !             whereas temperature is indexed from 0 to upn.
    !            The first row of the matrix is the equation for temp(0,ew,ns),
    !             the second row is the equation for temp(1,ew,ns), and so on.

    real(dp), intent(in) :: dttem       ! time step
    integer, intent(in) :: upn          ! number of layer interfaces
    real(dp), dimension(upn-1), intent(in) :: stagsigma    ! sigma coordinate at temp nodes
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    integer, intent(in) :: floating_mask
    real(dp), intent(in) ::  thck       ! ice thickness
    real(dp), dimension(0:upn), intent(in) ::  temp     ! ice temperature
    real(dp), dimension(upn-1), intent(in) :: dissip     ! interior heat dissipation (deg)
    real(dp), intent(in) :: bheatflx    ! geothermal flux (W m-2), positive down
    real(dp), intent(in) :: bfricflx    ! basal friction heat flux (W m-2), >= 0

    ! local variables

    real(dp) :: pmptempb  ! pressure melting temp at bed
    real(dp) :: fact
    real(dp) :: dsigbot   ! bottom layer thicknes in sigma coords

    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements

    ! upper boundary: set to surface air temperature

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = temp(0)

    ! ice interior, layers 1:upn-1  (matrix elements 2:upn)

    if (crank_nicolson) then  ! C-N can lead to oscillations in thin ice; currently deprecated

       fact = (tim0*dttem * coni / (2.d0 * rhoi*shci * thk0**2)) / thck**2
       subd(2:upn) = -fact * dups(1:upn-1,1)
       supd(2:upn) = -fact * dups(1:upn-1,2)
       diag(2:upn) = 1.0d0 - subd(2:upn) - supd(2:upn)
       rhsd(2:upn) =  temp(1:upn-1) * (2.0d0 - diag(2:upn)) &
                    - temp(0:upn-2) * subd(2:upn) &
                    - temp(2:upn  ) * supd(2:upn) & 
                    + dissip(1:upn-1)

    else   ! fully implicit

       !TODO - Remove redundant factors of 2 (results in roundoff-level answer change)
       fact = 2.d0 * (tim0*dttem * coni / (2.d0 * rhoi*shci * thk0**2)) / thck**2
       subd(2:upn) = -fact * dups(1:upn-1,1)
       supd(2:upn) = -fact * dups(1:upn-1,2)
       diag(2:upn) = 1.0d0 - subd(2:upn) - supd(2:upn)
       rhsd(2:upn) = temp(1:upn-1) + dissip(1:upn-1)

    endif    ! crank_nicolson

    ! basal boundary:
    ! for grounded ice, a heat flux is applied
    ! for floating ice, the basal temperature is held constant

    !NOTE: This lower BC is different from the one in glide_temp.
    !      If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
    !       That is, net heat flux at lower boundary must equal zero.
    !      If T(upn) >= Tpmp, then set T(upn) = Tpmp

    if (floating_mask == 1) then

       supd(upn+1) = 0.0d0
       subd(upn+1) = 0.0d0
       diag(upn+1) = 1.0d0
       rhsd(upn+1) = temp(upn) 

    else    ! grounded ice

       call glissade_calcpmpt_bed(pmptempb, thck)

       if (abs(temp(upn) - pmptempb) < 0.001d0) then  ! melting

          ! hold basal temperature at pressure melting point

          supd(upn+1) = 0.0d0
          subd(upn+1) = 0.0d0
          diag(upn+1) = 1.0d0
          rhsd(upn+1) = pmptempb

       else   ! frozen at bed
              ! maintain balance of heat sources and sinks
              ! (conductive flux, geothermal flux, and basal friction)
              ! Note: Heat fluxes are positive down, so slterm <= 0 and bheatflx <= 0.

          ! Note: bheatflx is generally <= 0, since defined as positive down.

          ! calculate dsigma for the bottom layer between the basal boundary and the temp. point above
          dsigbot = (1.0d0 - stagsigma(upn-1))

          ! =====Backward Euler flux basal boundary condition=====
           ! MJH: If Crank-Nicolson is desired for the b.c., it is necessary to
           ! ensure that the i.c. temperature for the boundary satisfies the
           ! b.c. - otherwise oscillations will occur because the C-N b.c. only
           ! specifies the basal flux averaged over two consecutive time steps.
          subd(upn+1) = -1.0d0
          supd(upn+1) =  0.0d0 
          diag(upn+1) =  1.0d0 
          rhsd(upn+1) = (bfricflx - bheatflx) * dsigbot*thck*thk0 / coni

       endif   ! melting or frozen

    end if     ! floating or grounded

  end subroutine glissade_findvtri

!=======================================================================

  subroutine glissade_enthalpy_findvtri (dttem,                       &
                                         upn,       stagsigma,        &
                                         subd,      diag,             &
                                         supd,      rhsd,             &
                                         dups,      floating_mask,    &
                                         thck,                        &
                                         temp,      waterfrac,        &
                                         enthalpy,  dissip,           &
                                         bheatflx,  bfricflx,         &
                                         alpha_enth,                  &
                                         verbose_column_in)

    ! solve for tridiagonal entries of sparse matrix

    use glimmer_physcon,  only : rhoi, shci, lhci, rhow, coni
    use glimmer_paramets, only : thk0, tim0

    !WHL - debug
    use parallel, only: this_rank

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    ! whereas temperature/enthalpy is indexed from 0 to upn.
    ! The first row of the matrix is the equation for enthalpy(0),
    ! the last row is the equation for enthalpy(upn), and so on.

    !I/O variables
    real(dp), intent(in) :: dttem       ! time step
    integer, intent(in) :: upn          ! number of layer interfaces
    real(dp), dimension(upn-1), intent(in) :: stagsigma    ! sigma coordinate at temp/enthalpy nodes
    real(dp), dimension(:,:), intent(in) :: dups   ! vertical grid quantities
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd
    integer, intent(in) :: floating_mask
    real(dp), intent(in) :: thck        ! ice thickness
    real(dp), dimension(0:upn), intent(in) :: temp       ! temperature (deg C)
    real(dp), dimension(upn-1), intent(in) :: waterfrac  ! water fraction (unitless)
    real(dp), dimension(0:upn), intent(in) :: enthalpy   ! specific enthalpy (J/m^3)
    real(dp), dimension(upn-1), intent(in) :: dissip     ! interior heat dissipation (deg)
    real(dp), intent(in) :: bheatflx   ! geothermal flux (W m-2), positive down
    real(dp), intent(in) :: bfricflx   ! basal friction heat flux (W m-2), >= 0
    real(dp), dimension(:), intent(out) :: alpha_enth  ! half-node diffusivity (m^2/s) for enthalpy
	                                               ! located halfway between temperature points

    logical, intent(in), optional :: verbose_column_in   ! if true, print debug statements for this column
    
    ! local variables
    real(dp) :: dsigbot  ! bottom layer thicknes in sigma coords.
    real(dp) :: alphai ! cold ice diffusivity
    real(dp) :: alpha0 ! temperate ice diffusivity
    real(dp) :: fact ! coefficient in tridiag
    integer  :: up
    real(dp), dimension(0:upn) :: pmptemp    ! pmptemp all nodes (deg C)
    real(dp), dimension(0:upn) :: enth_T     ! temperature part of specific enthalpy (J/m^3)
    real(dp) :: denth    ! enthalpy difference between adjacent layers
    real(dp) :: denth_T  ! difference in temperature component of enthalpy between adjacent layers
    real(dp) :: alpha_fact  ! factor for averaging diffusivity, 0 <= fact <= 1
    logical :: verbose_column ! if true, print debug statements for this column

    logical, parameter :: &
         alpha_harmonic_avg  = .false.  ! if true, take harmonic average of alpha in adjacent layers
                                        ! if false, take arithmetic average

    if (present(verbose_column_in)) then
       verbose_column = verbose_column_in
    else
       verbose_column = .false.
    endif

    ! define diffusivities alpha_i and alpha_0
    alphai = coni / rhoi / shci
    alpha0 = alphai / 100.0d0
	
    ! find pmptemp for this column (interior nodes and boundary)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:upn-1), thck, stagsigma(1:upn-1))
    call glissade_calcpmpt_bed(pmptemp(upn), thck)

    !WHL - debug                                                                                                                       
    if (verbose_column) then
       print*, ' '
       print*, 'Starting enthalpy calc'
       print*, 'k, temp, wfrac, enthalpy/(rhoi*ci), pmpt:'
       up = 0
       print*, up, temp(up), 0.d0, enthalpy(up)/(rhoi*shci)
       do up = 1, upn-1
          print*, up, temp(up), waterfrac(up), &
               enthalpy(up)/(rhoi*shci), pmptemp(up)
       enddo
       up = upn
       print*, up, temp(up), 0.d0, enthalpy(up)/(rhoi*shci), pmptemp(up)
    endif

    !WHL - Commenting out the following and replacing it with a new way of computing alpha.
    !      The commented-out code can result in sudden large changes in alpha that
    !       lead to oscillations in the thickness, temperature and velocity fields.
    !      These oscillations have a period of ~1 yr or more, spatial scale of
    !       many grid cells, and amplitude of ~10 m in thickness, 1 deg in temperature,
    !       and 2 m/s in velocity.

    ! create a column vector of size (0:upn) of diffusivity based on 
    ! previous timestep's temp.  Boundary nodes need a value so half-node
    ! diffusivity can be calculated at interior nodes (1:upn-1)

!    do up = 0,upn
!       if (temp(up) < pmptemp(up)) then
!          alpha(up) = alphai
!       else
!          alpha(up) = alpha0
!       endif
!    end do
    
    ! Find half-node diffusivity using harmonic average between nodes.
    ! The vector will be size (1:upn) - the first value is the half-node
    ! between nodes 0 and 1, the last value is the half-node between
    ! nodes upn-1 and upn. 
    
!    do up = 1,upn
!       alpha_enth(up) = 2.d0 / ((1.d0/alpha(up-1)) + (1.d0/alpha(up)))
!    end do
       
    !--------------------------------------------------------------------
    !WHL - Trying a different approach to the diffusivity at layer interfaces.
    ! Let d(enth)/dz = the gradient of enthalpy
    ! Can write 
    !    d(enth)/dz = d(enth_T)/dz + d(enth_w)/dz,
    ! where
    !    enth_T = (1-phi_w) * rhoi*ci*T
    !    enth_w =    phi_w  * rhow*(L + ci*Tpmp)
    !
    ! Now let f = d(enth_T)/z / d(enth)/dz
    !   (f -> 0 if f is computed to be negative)
    ! For cold ice, f = 1 and alpha = alphai
    ! For temperate ice, f ~ 0 and alpha = alpha0
    ! At the interface between cold and temperate ice,
    !  f ~ 0 if the temperate ice has large phi_w, but
    !  f ~ 1 if the temperate ice has close to zero phi_w.
    ! Two ways to average:
    ! (1) arithmetic average:  alpha = f*alphai + (1-f)*alpha0
    ! (2) harmonic average:    alpha = 1 / (f/alphai + (1-f)/alpha0).
    ! Both methods have the same asymptotic values at f = 0 or 1,
    !  but the arithmetic average gives greater diffusivity for
    !  intermediate values.
    !
    ! Still to be determined which is more accurate.
    ! The harmonic average allows large temperature gradients between the 
    !  bottom layer and the next layer up; the arithmetic average gives
    !  smoother gradients.
    !--------------------------------------------------------------------
    !
    ! At each temperature point, compute the temperature part of the enthalpy.
    ! enth_T = enth for cold ice, enth_T < enth for temperate ice

    do up = 0, upn
       enth_T(up) = (1.d0 - waterfrac(up)) * rhoi*shci*temp(up)
    enddo

!WHL - debug
    if (verbose_column) then
       print*, ' '
       print*, 'k, denth_T/(rhoi*shci), denth/(rhoi*shci), alpha_fact, alpha_enth(up):'
    endif

    ! Compute factors relating the temperature gradient to the total enthalpy gradient.
    ! Use these factors to average the diffusivity between adjacent temperature points.
    do up = 1,upn
       denth   = enthalpy(up) - enthalpy(up-1)
       denth_T = enth_T(up) - enth_T(up-1)   ! = denth in cold ice, < denth in temperate ice
       if (abs(denth) > 1.d-20 * rhow*lhci) then
          alpha_fact = max(0.d0, denth_T/denth)
          alpha_fact = min(1.d0, alpha_fact)
       else
          alpha_fact = 0.d0
       endif

       if (alpha_harmonic_avg) then  ! take a harmonic average
                                     ! This gives slower cooling of temperate layers and allows
                                     !  large temperature gradients between cold and temperate layers
          alpha_enth(up) = 1.d0 / ((alpha_fact/alphai) + (1.d0-alpha_fact)/alpha0)
       else   ! take an arithmetic average
              ! This gives faster cooling of temperate layers and smaller gradients 
          alpha_enth(up) = alpha_fact*alphai + (1.d0-alpha_fact)*alpha0
       endif

!WHL - debug
       if (verbose_column) then
          print*, up, denth_T/(rhoi*shci), denth/(rhoi*shci), alpha_fact, alpha_enth(up)
       endif

    end do

    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements
    
    ! upper boundary: set to surface air temperature
    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = min(0.0d0,temp(0)) * rhoi*shci
  
    ! RJH - Multiplied fact by a factor of 2 to become backward Euler coefficients
    fact = 2.d0 * (tim0*dttem / (2.0d0 * thk0**2)) / thck**2

    ! RJH - Altered rhsd to become fully implicit (backward Euler). 
    !       This included deleting subd and supd terms and dropping enthalpy(1:upn-1) coefficient

    ! ice interior. layers 1:upn-1  (matrix elements 2:upn)
    subd(2:upn) = -fact * alpha_enth(1:upn-1) * dups(1:upn-1,1)                                
    supd(2:upn) = -fact * alpha_enth(2:upn) * dups(1:upn-1,2)                                
    diag(2:upn) = 1.0d0 - subd(2:upn) - supd(2:upn)                                
    rhsd(2:upn) = enthalpy(1:upn-1) + dissip(1:upn-1) * rhoi * shci
                              
    ! BDM I'm assuming that dissip has units of phi/rhoi/shci.
    ! For an enthalpy calc, we want just phi, hence dissip * rhoi * shci
	
    ! basal boundary:
    ! for grounded ice, a heat flux is applied
    ! for floating ice, the basal temperature is held constant

    !NOTE: This lower BC is different from the one in standard glide_temp.
    !      If T(upn) < T_pmp, then require dT/dsigma = H/k * (G + taub*ubas)
    !       That is, net heat flux at lower boundary must equal zero.
    !      If T(upn) >= Tpmp, then set T(upn) = Tpmp

    if (floating_mask == 1) then

       supd(upn+1) = 0.0d0 
       subd(upn+1) = 0.0d0
       diag(upn+1) = 1.0d0
       rhsd(upn+1) = enthalpy(upn)
    
    else    ! grounded ice

       !WHL - debug
       if (verbose_column) then
          up = upn-1
          print*, 'temp(upn-1), pmptemp(upn-1):', temp(up), pmptemp(up)
          up = upn
          print*, 'temp(upn), pmptemp(upn):', temp(up), pmptemp(up)
       endif

    !Positive-Thickness Basal Temperate Boundary Layer

    !WHL - Not sure whether this condition is right.  
    !      It implies that the enthalpy at the bed (upn) = enthalpy in layer (upn-1). 
       if (abs(temp(upn-1) - pmptemp(upn-1)) < 0.001d0) then   
       
          subd(upn+1) = -1.0d0
          supd(upn+1) =  0.0d0 
          diag(upn+1) = 1.0d0 
          rhsd(upn+1) = 0.0d0

          !WHL - debug
          if (verbose_column) then
             print*, 'basal BC: branch 1 (finite-thck BL)'
          endif

          !Zero-Thickness Basal Temperate Boundary Layer
       elseif (abs(temp(upn) -  pmptemp(upn)) < 0.001d0) then  ! melting
          
          ! hold basal temperature at pressure melting point
          supd(upn+1) = 0.0d0
          subd(upn+1) = 0.0d0
          diag(upn+1) = 1.0d0
          rhsd(upn+1) = pmptemp(upn) * rhoi * shci
          
          !WHL - debug
          if (verbose_column) then
             print*, 'basal BC: branch 2 (zero-thck BL)'
          endif
          
       else  
          
          !WHL - debug
          if (verbose_column) then
             print*, 'basal BC: branch 3 (cold ice)'
          endif
          
          ! frozen at bed
          ! maintain balance of heat sources and sinks
          ! (conductive flux, geothermal flux, and basal friction)
          ! Note: Heat fluxes are positive down, so slterm <= 0 and bheatflx <= 0.
          
          ! Note: The heat source due to basal sliding (bfricflx) is computed in subroutine calcbfric.
          ! Also note that bheatflx is generally <= 0, since defined as positive down.
          
          ! calculate dsigma for the bottom layer between the basal boundary and the temp. point above
          dsigbot = (1.0d0 - stagsigma(upn-1))                                                                  
          
          ! =====Backward Euler flux basal boundary condition=====
          ! MJH: If Crank-Nicolson is desired for the b.c., it is necessary to
          ! ensure that the i.c. temperature for the boundary satisfies the
          ! b.c. - otherwise oscillations will occur because the C-N b.c. only
          ! specifies the basal flux averaged over two consecutive time steps.
          subd(upn+1) = -1.0d0
          supd(upn+1) =  0.0d0 
          diag(upn+1) = 1.0d0 
          rhsd(upn+1) = (bfricflx - bheatflx) * dsigbot*thck*thk0 * rhoi*shci/coni
          ! BDM temp approach should work out to be dT/dsigma, so enthalpy approach
          ! should just need dT/dsigma * rhoi * shci for correct units

       endif   ! melting or frozen
       
    end if     ! floating or grounded

  end subroutine glissade_enthalpy_findvtri

!=======================================================================

  subroutine glissade_calcbmlt(whichtemp,     dttem,         &
                               ewn,           nsn,           &
                               upn,                          &
                               sigma,         stagwbndsigma, &
                               ice_mask,      floating_mask, &
                               thck,                         &
                               temp,          waterfrac,     &
                               bfricflx,      bheatflx,      &
                               lcondflx,                     &
                               bwat,          bmlt)

    ! Compute the amount of basal melting.
    ! The basal melting computed here is applied to the ice thickness
    !  by glissade_transport_driver, conserving mass and energy.
    !
    ! Any internal temperatures above the pressure melting point are reset to the
    !  pmp temperature, with excess energy applied toward basal melting.
    !  Hopefully this is rare.
    ! TODO: Moving all internal melting to the basal surface is not very realistic 
    !       and should be revisited.

    use glimmer_physcon, only: shci, rhoi, lhci
    use glimmer_paramets, only : thk0, tim0

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: whichtemp          ! temperature method (TEMP_PROGNOSTIC or TEMP_ENTHALPY)
    real(dp), intent(in) :: dttem             ! time step
    integer, intent(in) :: ewn, nsn, upn      ! grid dimensions
    real(dp), dimension(:),      intent(in) :: sigma           ! vertical sigma coordinate
    real(dp), dimension(0:),     intent(in) :: stagwbndsigma   ! staggered vertical coordinate for temperature
    real(dp), dimension(0:,:,:), intent(inout) :: temp         ! temperature (deg C)
    real(dp), dimension(:,:,:),  intent(inout) :: waterfrac    ! water fraction
    real(dp), dimension(:,:),    intent(in) :: &
         thck,                 & ! ice thickness
         bfricflx,             & ! basal frictional heating flux (W m-2), >= 0
         bheatflx,             & ! geothermal heating flux (W m-2), positive down
         lcondflx,             & ! heat conducted from ice interior to bed (W m-2), positive down
         bwat                    ! depth of basal water

    integer, dimension(:,:), intent(in) ::  &
         ice_mask,           &! = 1 where ice exists (thck > thklim_temp), else = 0
         floating_mask        ! = 1 where ice is floating, else = 0

    real(dp), dimension(:,:),    intent(out):: bmlt    ! scaled melt rate (m/s * tim0/thk0)
                                                       ! > 0 for melting, < 0 for freeze-on

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    real(dp), dimension(size(stagwbndsigma))    :: pmptemp   ! pressure melting point temperature
    real(dp) :: bflx    ! heat flux available for basal melting (W/m^2)
    integer :: up, ew, ns

    real(dp) :: layer_thck    ! layer thickness (m)
    real(dp) :: melt_energy   ! energy available for internal melting (J/m^2)
    real(dp) :: internal_melt_rate   ! internal melt rate, transferred to bed (m/s)
    real(dp) :: melt_fact     ! factor for bmlt calculation
    real(dp) :: hmlt          ! melt thickness associated with excess meltwater

    real(dp), parameter :: max_waterfrac = 0.01d0   ! maximum allowed water fraction
                                                    ! excess water drains to the bed

    bmlt(:,:) = 0.0d0
    melt_fact = tim0 / (thk0 * lhci * rhoi) 
    
    !TODO - Change to 1, nsn
    do ns = 2, nsn
       do ew = 2, ewn

          if (ice_mask(ew,ns)==1 .and. floating_mask(ew,ns)==0) then

             ! Basal friction term is computed above in subroutine glissade_calcbfric,
             !  or in the Glissade velocity solver.
             !
             ! Compute basal melting
             ! Note: bmlt > 0 for melting, < 0 for freeze-on
             !       bfricflx >= 0 by definition
             !       bheatflx is positive down, so usually bheatflx < 0 (with negative values contributing to melt)
             !       lcondflx is positive down, so lcondflx < 0 for heat is flowing from the bed toward the surface
             !
             !       This equation allows for freeze-on (bmlt < 0) if the conductive term 
             !        (lcondflx, positive down) is carrying enough heat away from the boundary.  
             !       But freeze-on requires a local water supply, bwat > 0.
             !       When bwat = 0, we reset the bed temperature to a value slightly below the melting point.

             bflx = bfricflx(ew,ns) + lcondflx(ew,ns) - bheatflx(ew,ns)  ! W/m^2
             bmlt(ew,ns) = bflx * melt_fact

             ! Add internal melting

             if (whichtemp == TEMP_ENTHALPY) then

                ! Add internal melting associated with waterfrac > max_waterfrac (1%)
                !TODO - Define max_waterfrac above

!!             call glissade_calcpmpt(pmptemp(:), thck(ew,ns),   &
!!                                    stagsigma(:) )

                !WHL - Any correction for rhoi/rhow here?
                do up = 1, upn-1
                   if (waterfrac(up,ew,ns) > max_waterfrac) then
                      ! compute melt rate associated with excess water
                      hmlt = (waterfrac(up,ew,ns) - max_waterfrac) * (thck(ew,ns)*thk0) * (sigma(up+1) - sigma(up))  ! m
                      internal_melt_rate = hmlt / (dttem*tim0)          ! m/s
                      ! transfer meltwater to the bed
                      bmlt(ew,ns) = bmlt(ew,ns) + internal_melt_rate * tim0/thk0
                      ! reset waterfrac to max value
                      waterfrac(up,ew,ns) = max_waterfrac
                   endif
                enddo

             else   ! whichtemp = TEMP_PROGNOSTIC

                ! Add internal melting associated with temp > pmptemp
                ! Note: glissade_calcpmpt does not compute pmpt at the top surface or the bed.

                !TODO - Make sure the '0' starting index is handled correctly
                call glissade_calcpmpt(pmptemp(:), thck(ew,ns),   &
                                       stagwbndsigma(:) )

                do up = 1, upn-1
                   if (temp(up,ew,ns) > pmptemp(up)) then
                      ! compute excess energy available for melting
                      layer_thck = thck(ew,ns)*thk0 * (sigma(up+1) - sigma(up))  ! m
                      melt_energy = rhoi * shci * (temp(up,ew,ns) - pmptemp(up)) * layer_thck         ! J/m^2
                      internal_melt_rate = melt_energy / (rhoi * lhci * dttem*tim0)  ! m/s
                      ! transfer internal melting to the bed
                      bmlt(ew,ns) = bmlt(ew,ns) + internal_melt_rate * tim0/thk0  ! m/s * tim0/thk0
                      ! reset T to Tpmp
                      temp(up,ew,ns) = pmptemp(up)
                   endif
                enddo

             endif    ! whichtemp

             ! Cap basal temp at pmptemp, if necessary

             call glissade_calcpmpt_bed(pmptemp(upn), thck(ew,ns))
             temp(upn,ew,ns) = min (temp(upn,ew,ns), pmptemp(upn))

             ! If freeze-on was computed above (bmlt < 0) and Tbed = Tpmp but no basal water is present, then set T(upn) < Tpmp.
             ! Note: In subroutine findvtri, we solve for Tbed (instead of holding it at Tpmp) when Tbed < 0.001.
             !       With an offset here of 0.01, we will solve for T_bed at the next timestep.
             ! Note: Energy is not exactly conserved here.

             if (bmlt(ew,ns) < 0.d0 .and. bwat(ew,ns)==0.d0 .and. temp(upn,ew,ns) >= pmptemp(upn)) then
                temp(upn,ew,ns) = pmptemp(upn) - 0.01d0
             endif

          endif   ! thk > thklim_temp

       enddo
    enddo

  end subroutine glissade_calcbmlt

!=======================================================================

  !TODO - Inline glissade_calcpmpt and glissade_calcbpmp above?

  subroutine glissade_calcpmpt(pmptemp, thck, stagsigma)

    ! Compute the pressure melting point temperature in the column
    ! (but not at the surface or bed).
    ! Note: pmptemp and stagsigma should have dimensions (1:upn-1).

    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), dimension(:), intent(out) :: pmptemp  ! pressure melting point temperature (deg C)
    real(dp), intent(in) :: thck                    ! ice thickness
    real(dp), intent(in), dimension(:) :: stagsigma ! staggered vertical coordinate
                                                    ! (defined at layer midpoints)

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp(:) = fact * thck * stagsigma(:)

  end subroutine glissade_calcpmpt

!=======================================================================

  subroutine glissade_therm_calcbpmp(ewn,   nsn,  &
                                     thck,  bpmp)

    ! Calculate the pressure melting point at the base of the ice sheet

    integer, intent(in) ::  ewn, nsn    ! grid dimensions

    real(dp), dimension(:,:), intent(in)  :: thck  ! ice thickness (dimensionless)
    real(dp), dimension(:,:), intent(out) :: bpmp  ! bed pressure melting point (deg C)

    integer :: ew,ns

    bpmp(:,:) = 0.d0

    do ns = 1, nsn
       do ew = 1, ewn
          call glissade_calcpmpt_bed(bpmp(ew,ns),thck(ew,ns))
       end do
    end do

  end subroutine glissade_therm_calcbpmp

!=======================================================================

  subroutine glissade_calcpmpt_bed(pmptemp_bed, thck)

    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), intent(out) :: pmptemp_bed ! pressure melting point temp at bed (deg C)
    real(dp), intent(in) :: thck         ! ice thickness

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp_bed = fact * thck 

  end subroutine glissade_calcpmpt_bed

!=======================================================================

  subroutine glissade_therm_calcflwa(stagsigma,   thklim,   &
                                     flwa,        temp,     &
                                     thck,        flow_factor, &
                                     default_flwa_arg,      &
                                     flag,        waterfrac)

    ! Calculates Glen's $A$ over the three-dimensional domain,
    ! using one of three possible methods.
    !
    ! The primary method is to use this equation from \emph{Paterson and Budd} [1982]:
    ! \[
    ! A(T^{*})=a \exp \left(\frac{-Q}{RT^{*}}\right)
    ! \]
    ! This is equation 9 in {\em Payne and Dongelmans}. $a$ is a constant of proportionality,
    ! $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
    ! The pressure-corrected temperature, $T^{*}$ is given by:
    ! \[
    ! T^{*}=T-T_{\mathrm{pmp}}+T_0
    ! \] 
    ! \[
    ! T_{\mathrm{pmp}}=T_0-\sigma \rho g H \Phi
    ! \]
    ! $T$ is the ice temperature, $T_{\mathrm{pmp}}$ is the pressure melting point 
    ! temperature, $T_0$ is the triple point of water, $\rho$ is the ice density, and 
    ! $\Phi$ is the (constant) rate of change of melting point temperature with pressure.

    use glimmer_physcon
    use glimmer_paramets, only : thk0, vis0

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

!   Note: The flwa, temp, and stagsigma arrays should have the same vertical dimension
!         (1:upn-1 on the staggered vertical grid).

    real(dp),dimension(:),      intent(in)    :: stagsigma ! vertical coordinate at layer midpoints
    real(dp),                   intent(in)    :: thklim    ! thickness threshold
    real(dp),dimension(:,:,:),  intent(in)    :: temp      ! 3D temperature field
    real(dp),dimension(:,:),    intent(in)    :: thck      ! ice thickness
    real(dp)                                  :: flow_factor ! fudge factor in Arrhenius relationship
    real(dp),                   intent(in)    :: default_flwa_arg ! Glen's A to use in isothermal case 
                                                                  ! Units: Pa^{-n} yr^{-1} 
    integer,                    intent(in)    :: flag      !> Flag to select the method of calculation
    real(dp),dimension(:,:,:),  intent(out)   :: flwa      !> The calculated values of $A$
    real(dp),dimension(:,:,:),  intent(in), optional :: waterfrac !> internal water content fraction, 0 to 1

    !> \begin{description}
    !> \item[0] {\em Paterson and Budd} relationship.
    !> \item[1] {\em Paterson and Budd} relationship, with temperature set to -5$^{\circ}$C.
    !> \item[2] Set to prescribed constant value.
    !> \end{description}

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp) :: default_flwa
    integer :: ew, ns, up, ewn, nsn, uflwa
    real(dp) :: tempcor

    real(dp), parameter :: fact = grav * rhoi * pmlt * thk0

    real(dp),dimension(4), parameter ::  &
       arrfact = (/ arrmlh / vis0,      &   ! Value of A when T* is above -263K
                    arrmll / vis0,      &   ! Value of A when T* is below -263K
                   -actenh / gascon,    &   ! Value of -Q/R when T* is above -263K
                   -actenl / gascon/)       ! Value of -Q/R when T* is below -263K
    
    real(dp), parameter :: const_temp = -5.0d0
    real(dp), parameter :: flwa_waterfrac_enhance_factor = 181.25d0

    !------------------------------------------------------------------------------------ 
   
    uflwa=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    ! Check that the temperature array has the desired vertical dimension

    if (size(temp,1) /= size(flwa,1)) then
       call write_log('glissade_calcflwa: temp and flwa must have the same vertical dimensions', GM_FATAL)
    endif

    ! Scale the default rate factor (default value has units Pa^{-n} yr^{-1}).
    ! Also multiply by fudge factor

    default_flwa = flow_factor * default_flwa_arg / (vis0*scyr)
    !write(*,*)"Default flwa = ",default_flwa

    select case(flag)

    case(FLWA_PATERSON_BUDD)

      ! This is the Paterson and Budd relationship
      ! BDM added waterfrac relationship for whichtemp=TEMP_ENTHALPY case

      do ns = 1,nsn
         do ew = 1,ewn

            if (thck(ew,ns) > thklim) then
            
               do up = 1, uflwa   ! uflwa = upn - 1 (values at layer midpoints)

                  ! Calculate the corrected temperature

                  tempcor = min(0.0d0, temp(up,ew,ns) + thck(ew,ns)*fact*stagsigma(up))
                  tempcor = max(-50.0d0, tempcor)

                  ! Calculate Glen's A (including flow fudge factor)

                  if (tempcor >= -10.d0) then
                     flwa(up,ew,ns) = flow_factor * arrfact(1) * exp(arrfact(3)/(tempcor + trpt))
                  else
                     flwa(up,ew,ns) = flow_factor * arrfact(2) * exp(arrfact(4)/(tempcor + trpt))
                  endif

                  ! BDM added correction for a liquid water fraction 
                  ! Using Greve and Blatter, 2009 formulation for Glen's A flow rate factor:
                  !    A = A(theta_PMP) * (1 + 181.25 * waterfrac)
		  ! RJH - commenting out waterfrac correction to explore causes of
		  ! oscillations in thk and vel for EISMINT-2 test cases
                  if (present(waterfrac)) then
                     if (waterfrac(up,ew,ns) > 0.0d0) then
                        flwa(up,ew,ns) = flwa(up,ew,ns) * (1.d0 + flwa_waterfrac_enhance_factor * waterfrac(up,ew,ns))      
                     endif
                  endif

               enddo   ! up

            else   ! thck < thklim

               flwa(:,ew,ns) = default_flwa

            end if

         end do
      end do

    case(FLWA_PATERSON_BUDD_CONST_TEMP)

      ! This is the Paterson and Budd relationship, but with the temperature held constant
      ! at -5 deg C

      do ns = 1,nsn
         do ew = 1,ewn

            if (thck(ew,ns) > thklim) then

               ! Calculate Glen's A with a fixed temperature (including flow fudge factor)

               if (const_temp >= -10.d0) then
                  flwa(:,ew,ns) = flow_factor * arrfact(1) * exp(arrfact(3)/(const_temp + trpt))
               else
                  flwa(:,ew,ns) = flow_factor * arrfact(2) * exp(arrfact(4)/(const_temp + trpt))
               endif

            else

               flwa(:,ew,ns) = default_flwa

            end if

         end do
      end do

    case(FLWA_CONST_FLWA) 

      flwa(:,:,:) = default_flwa
  
    end select

  end subroutine glissade_therm_calcflwa

!=======================================================================

  subroutine glissade_enth2temp (enthalpy, temp, waterfrac, thck, stagsigma)

    ! BDM convert from specific enthalpy to ice temperature and water content
    ! takes a vertical column of size enthalpy(dup-1,1,1) and converts to
    ! temp(1:dup-1,1,1) and waterfrac(1:dup-1,1,1)
    ! enthalpy(1:dup-1,1,1)
    ! temp(1:dup-1,1,1)
    ! waterfrac(1:dup-1,1,1)
    ! thck(1,1,1)
    ! stagsigma(1:dup-1,1,1)
	
    use glimmer_physcon, only : rhoi, shci, lhci, rhow

    ! I/O variables
    real(dp), dimension(0:), intent(inout)               :: enthalpy !enthalpy is (0:upn)
    real(dp), intent(in)                                 :: thck
    real(dp), dimension(:), intent(in)                   :: stagsigma !stagsigma is (1:upn-1)
    real(dp), dimension(0:size(enthalpy)-1), intent(out) :: temp !temp is (0:upn)
    real(dp), dimension(1:size(enthalpy)-2), intent(out) :: waterfrac !waterfrac is (1:upn-1)

    ! local variables
    real(dp), dimension(0:size(enthalpy)-1)              :: pmptemp ! (0:upn)
    real(dp), dimension(0:size(enthalpy)-1)              :: pmpenthalpy ! (0:upn)
    integer                                              :: upn !used for convenience
    integer                                              :: up
	
    upn = size(enthalpy)-1
	
    ! Find pmpenthalpy(0:upn)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:upn-1), thck, stagsigma(1:upn-1))
    call glissade_calcpmpt_bed(pmptemp(upn), thck)
	
    pmpenthalpy = pmptemp * rhoi * shci
	
    !solve for temp and waterfrac
    if(enthalpy(0) >= pmpenthalpy(0)) then ! temperate ice
       temp(0) = pmptemp(0)                ! temperate ice
       !WHL - Resetting enthalpy so that it's consistent with the new temperature
       !      This is consistent with energy conservation because the top surface
       !       is infinitesimally thin.
       enthalpy(0) = pmpenthalpy(0)
    else
       temp(0) = enthalpy(0) / (rhoi*shci) ! temp is cold
    endif
	
    do up = 1, upn-1
       if(enthalpy(up) >= pmpenthalpy(up)) then ! temperate ice
          temp(up) = pmptemp(up)                ! temp = pmptemp
          waterfrac(up) = (enthalpy(up)-pmpenthalpy(up)) /                 &
                          ((rhow-rhoi) * shci * pmptemp(up) + rhow * lhci)
       else ! cold ice

          !WHL - debug
          if (waterfrac(up) > 0.d0) then
             print*, 'Zeroing out waterfrac: k, waterfrac =', up, waterfrac(up)
          endif

          temp(up) = enthalpy(up) / (rhoi*shci) ! temp is cold
          waterfrac(up) = 0.0d0                 ! waterfrac = 0
       endif
    end do
	
    if(enthalpy(upn) >= pmpenthalpy(upn)) then  ! temperate ice
       temp(upn) = pmptemp(upn)                 ! temp = pmptemp
    else
       temp(upn) = enthalpy(upn) / (rhoi*shci)  ! temp is cold
       !WHL - Resetting enthalpy so that it's consistent with the new temperature
       !      This is consistent with energy conservation because the top surface
       !       is infinitesimally thin.
       enthalpy(upn) = pmpenthalpy(upn)
    endif

  end subroutine glissade_enth2temp

!=======================================================================

  subroutine glissade_temp2enth (enthalpy, temp, waterfrac, thck, stagsigma)

    ! BDM convert from temperature and water content and converts to specific enthalpy
    ! takes a vertical column of size temp(0:dup) and converts to enthalpy(0:dup,1,1)
    ! waterfrac is only size(1:dup-1), so will assume no waterfrac at boundaries

    use glimmer_physcon, only : rhoi, shci, lhci, rhow

    ! I/O variables
    real(dp), dimension(0:), intent(out)        :: enthalpy !enthalpy is (0:upn)
    real(dp), intent(in)                        :: thck
    real(dp), dimension(:), intent(in)          :: stagsigma !stagsigma is (1:upn-1)
    real(dp), dimension(0:), intent(in)         :: temp !temp is (0:upn)
    real(dp), dimension(:), intent(in)          :: waterfrac !waterfrac is (1:upn-1)

    ! local variables
    real(dp), dimension(0:size(temp)-1)         :: pmptemp !(0:upn)
    real(dp), dimension(0:size(temp)-1)         :: pmpenthalpy !(0:upn)
    integer                                     :: up
    integer                                     :: upn !used for convenience
		
    upn = size(temp)-1
	
    ! Find pmpenthalpy(0:dup,1,1)
    pmptemp(0) = 0.0d0
    call glissade_calcpmpt(pmptemp(1:upn-1), thck, stagsigma(1:upn-1))
    call glissade_calcpmpt_bed(pmptemp(upn), thck)
    
    !WHL - This variable is not used below
    pmpenthalpy = rhoi * shci * pmptemp
	
    ! solve for enthalpy
    ! assume waterfrac = 0 at upper and lower ice surfaces
    enthalpy(0) = temp(0) * rhoi * shci
    do up = 1, upn-1
       enthalpy(up) = ((1 - waterfrac(up)) * rhoi * shci * temp(up))          &
                      + waterfrac(up) * rhow * ((shci * pmptemp(up)) + lhci)
    end do
    enthalpy(upn) = temp(upn) * rhoi * shci
	
  end subroutine glissade_temp2enth

!=======================================================================

end module glissade_therm

!=======================================================================
