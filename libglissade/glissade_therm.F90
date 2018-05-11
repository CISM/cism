!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_therm.F90 - part of the Community Ice Sheet Model (CISM)  
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

  !-----------------------------------------------------------------------------
  ! This module combines two previous modules: glissade_temp.F90 and glissade_enthalpy.F90.
  ! It was created in Nov. 2014 by William Lipscomb (LANL).
  ! It is functionally equivalent to the two previous modules, but with modified data structures,
  !  streamlined organization, and exclusively SI units.
  ! 
  ! The module computes temperature diffusion and strain heating in a local column 
  !  without doing horizontal or vertical advection.
  ! Temperature advection is done separately, e.g. using the incremental remapping transport scheme.
  ! It is assumed here that temperature (and enthalpy) values are staggered in the
  !  vertical compared to the velocity.  That is, the temperature lives at layer midpoints 
  !  instead of layer interfaces (as in Glide).
  ! Temperature and enthalpy are also defined at the upper and lower surfaces with 
  !  appropriate boundary conditions.   
  !
  ! The glissade_temp.F90 module was based on glide_temp.F90, with modifications for Glissade.
  ! The glissade_enthalpy.F90 module contained modifications for a new enthalpy option;
  !  it was originally written by Brian Macpherson (CU) under the supervision of Hari Rajaram.
  !
  ! NOTE: SI units are used throughout this module.
  !-----------------------------------------------------------------------------

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade_therm

    use glimmer_global, only : dp 
    use glide_types
    use glimmer_log
    use parallel

    implicit none

    private
    public :: glissade_init_therm, glissade_therm_driver, glissade_flow_factor, glissade_pressure_melting_point, &
              glissade_interior_dissipation_sia, glissade_interior_dissipation_first_order,                         &
              glissade_enth2temp, glissade_temp2enth

    ! time-stepping scheme

    !NOTE:  For the dome test case, the Crank-Nicolson scheme can give unstable 
    !        temperature fluctuations for thin ice immediately after the ice 
    !        becomes thick enough for the temperature calculation.
    !       The fully implicit scheme has been stable for all cases (but is only
    !        first-order accurate in time). 
    !NOTE:  Crank-Nicolson is not supported for the enthalpy scheme.

    logical, parameter::   &
         crank_nicolson = .false.  ! if true, use Crank-Nicolson time-stepping
                                   ! if false, use fully implicit

    ! max and min allowed temperatures
    ! Temperatures sometimes go below -100 for cases where Crank-Nicholson is unstable
    real(dp), parameter ::   &
         maxtemp_threshold = 1.d11,   &
         mintemp_threshold = -100.d0

    real(dp), dimension(:,:), allocatable :: dups   ! vertical grid quantities

    ! local parameter for debugging
    logical, parameter:: verbose_therm = .false.  ! set to true for diagnostic column output

  contains

!****************************************************    

  subroutine glissade_init_therm (temp_init,  is_restart,        &
                                  ewn,        nsn,        upn,   &
                                  itest,      jtest,      rtest, &
                                  sigma,      stagsigma,         &
                                  thck,                          &
                                  artm,                          &
                                  acab,                          &
                                  bheatflx,                      &
                                  pmp_offset,                    &
                                  temp,       tempunstag)

    ! initialization subroutine for higher-order dycores, where temperature is defined at
    ! the midpoint of each layer plus the upper and lower surfaces

    use glimmer_physcon, only : trpt
    use glimmer_paramets, only : unphys_val

    ! In/out arguments

    integer, intent(in) :: &
         temp_init,       &! method for initializing the temperature
         is_restart        ! = 1 if restarting, else = 0

    integer, intent(in) :: &
         ewn, nsn, upn,      & ! grid dimensions
         itest, jtest, rtest   ! coordinates of diagnostic point

    real(dp), dimension(:), intent(in) ::   &
         sigma,           &! vertical coordinate, located at layer interfaces
         stagsigma         ! staggered vertical coordinate, located at the center of each layer

    real(dp), dimension(:,:), intent(in) ::  &
         thck,            &! ice thickness (m)
         artm,            &! surface air temperature (deg C)
         acab,            &! surface mass balance (m/s)
         bheatflx          ! basal heat flux (W/m^2), positive down

    real(dp), intent(in) :: &
         pmp_offset        ! offset of initial Tbed from pressure melting point temperature (deg C)

    real(dp), dimension(0:,:,:), intent(inout) ::  &
         temp              ! ice temperature
                           ! intent(inout) because it might have been read already from an input file,
                           !  but otherwise is set in this subroutine

    real(dp), dimension(:,:,:), intent(inout) ::  &
         tempunstag        ! ice temperature on unstaggered grid
                           ! may be used to initialize from an unstaggered glide temperature field
                           ! intent(inout) because it might have been read from an input file,
                           ! but otherwise it is set (diagnostically) in this subroutine
         
    ! Local variables

    character(len=100) :: message

    integer :: up, ns, ew

    logical :: verbose_column    ! if true, then write diagnostic info for the column

    ! Precompute some grid quantities used in the vertical temperature solve
 
    allocate(dups(upn+1,2))   !TODO - upn-1 instead?
    dups(:,:) = 0.0d0

    up = 1
    dups(up,1) = 1.d0/((sigma(up+1) - sigma(up)) * (stagsigma(up) - sigma(up)) )
                            
    do up = 2, upn-1
       dups(up,1) = 1.d0/((sigma(up+1) - sigma(up)) * (stagsigma(up) - stagsigma(up-1)) )
    enddo

    do up = 1, upn-2
       dups(up,2) = 1.d0/((sigma(up+1) - sigma(up)) * (stagsigma(up+1) - stagsigma(up)) )
    end do

    up = upn-1
    dups(up,2) = 1.d0/((sigma(up+1) - sigma(up)) * (sigma(up+1) - stagsigma(up)) )
                                    

    !==== Initialize ice temperature.============
    ! Following possibilities:
    ! (0) Set ice temperature to 0 C everywhere in column (TEMP_INIT_ZERO).
    ! (1) Set ice temperature to surface air temperature everywhere in column (TEMP_INIT_ARTM).
    ! (2) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
    !     at the bed (TEMP_INIT_LINEAR).
    !     A parameter (pmpt_offset) controls how far below Tpmp the initial bed temp is set.
    ! (3) Set up a temperature profile based on advective-diffusive balance, with T = artm
    !     at the surface and dT/dz = -F_geo/k at the bed (TEMP_INIT_ADVECTIVE_DIFFUSIVE).
    !     The temperature at each level is capped at the value computed by method (3).
    ! (4) Read ice temperature from external file (TEMP_INIT_EXTERNAL).
    !   (4a) If variable tempstag is present: Set ice temperature to tempstag from input file.
    !   (4b) If variable tempunstag is present (and tempstag is not): interpolate tempunstag
    !        to the staggered ice temperature needed by the Glissade dycore.
    !
    ! The default is (1).
    ! Methods (0-3) overwrite any temperature given in input files.
    ! Method (3) may be optimal for reducing spinup time in the interior of large ice sheets.
    ! Option (4) requires that temperature is present in the input file.

    if (is_restart == RESTART_TRUE) then

       ! Temperature has already been initialized from a restart file.
       ! (Temperature is always a restart variable.)

       call write_log('Initializing ice temperature from the restart file')

    elseif (temp_init == TEMP_INIT_EXTERNAL) then

       ! Temperature from external file

       ! Check for a possible input error.  If the user supplies a file with the 'temp' field, which has
       ! vertical dimension (1:upn), then the temperature in layers 1:upn may appear correct (though
       ! staggered incorrectly), but the temperature in layer 0 will remain at an unphysical value.
       ! Let the user know if this has happened.
       !WHL - Nov. 2014 - I verified that the code aborts here if temp (rather than tempstag) is in the restart file.

       if (minval(temp(0,:,:)) < (-1.d0*trpt) .and. minval(temp(1:upn,:,:)) > (-1.d0*trpt)) then
          call write_log('Error, temperature field has been read incorrectly. Note that the '  &
             // 'Glissade dycore must be initialized with tempstag, not temp.' &
             // 'You can rename temp to tempunstag if you want it to be interpolated ' &
             // 'to the vertically staggered tempstag (loss of detail).', GM_FATAL)
       endif


       if ( minval(temp) > (-1.0d0 * trpt)) then  ! temperature has been read from an input file
                                                      ! Note: trpt = 273.15 K

          ! Temperature has already been initialized from a restart or input file.
          ! We know this because the default initial temps of unphys_val
          ! (a large negative number) have been overwritten.

          ! Initialise tempunstag, the temperature on unstaggered vertical grid
          ! (not used for calculations)
          ! Use sigma weighted interpolation from temp to layer interfaces
          do up = 2, upn-1
            tempunstag(up,:,:) = temp(up-1,:,:) + (temp(up,:,:) - temp(up-1,:,:)) *   &
                 (sigma(up) - stagsigma(up-1)) / (stagsigma(up) - stagsigma(up-1))
          end do
          ! boundary conditions are identical on both grids, but temp starts at index 0
          tempunstag(1,:,:) = temp(0,:,:)
          tempunstag(upn,:,:) = temp(upn,:,:)

          call write_log('Initializing ice temperature from an input file')

       elseif ( maxval(tempunstag(:,:,:)) > (-1.0d0 * trpt)) then

          ! Test if any temperature in physical range
          ! If yes, we have data from restart or input for unstaggered tempunstag that we use to initialise temp

          call write_log('Initializing ice temperature by interpolating tempunstag from restart/input file')

          ! Set temp to linear interpolation from tempunstag at layer interfaces
          do up = 1, upn-1
            temp(up,:,:) = (tempunstag(up,:,:) + tempunstag(up+1,:,:)) * 0.5d0
          end do
          ! boundary conditions are identical on both grids, but temp starts at index 0
          temp(0,:,:) = tempunstag(1,:,:)
          temp(upn,:,:) = tempunstag(upn,:,:)

       else
          ! fatal error, neither tempstag nor tempunstag are present in the input files
          call write_log('Error, temp_init = 4 requires temperature variable tempstag or ' &
             // 'tempunstag specified in the restart or input files. ', GM_FATAL)
       endif


    else   ! not reading temperature from restart or input file
           ! initialize it here based on temp_init = (0-3)

       ! initialize T = 0 C everywhere
       temp(:,:,:) = 0.0d0
                   
       ! set temperature in each column based on the value of temp_init
                                 
       if (temp_init == TEMP_INIT_ZERO) then
          call write_log('Initializing ice temperature to 0 deg C')
       elseif (temp_init == TEMP_INIT_ARTM) then   ! initialize ice column temperature to min(artm, 0 C)
          call write_log('Initializing ice temperature to the surface air temperature')
       elseif (temp_init == TEMP_INIT_LINEAR) then ! initialize ice column temperature with a linear profile:
                                                   ! T = artm at the surface, and T <= Tpmp at the bed
          call write_log('Initializing ice temperature to a linear profile in each column')
          write(message,*) 'Offset from pressure melting point temperature =', pmp_offset
          call write_log(message)
       elseif (temp_init == TEMP_INIT_ADVECTIVE_DIFFUSIVE) then  ! initialize ice column temperature based on advective-diffusive balance
                                                                 ! T = artm at the surface, and dT/dz = -F_geo/k at the bed
          call write_log('Initializing ice temperature based on advective-diffusive balance in each column')
          write(message,*) 'Offset from pressure melting point temperature =', pmp_offset
          call write_log(message)
       else
          call write_log('Error: invalid temp_init option in glissade_init_therm. It is possible that the temperature' &
                      // 'was not read correctly from the input file, resulting in unphysical values', GM_FATAL)
       endif

       do ns = 1, nsn
          do ew = 1, ewn

             if (verbose_therm .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
                verbose_column = .true.
                print*, ' '
                print*, 'Initial temperature diagnostics: rank, i, j =', this_rank, ew, ns
             else
                verbose_column = .false.
             endif

             call glissade_init_temp_column(temp_init,          &
                                            stagsigma(:),       &
                                            thck(ew,ns),        &
                                            artm(ew,ns),        &
                                            acab(ew,ns),        &
                                            bheatflx(ew,ns),    &
                                            pmp_offset,         &
                                            temp(:,ew,ns),      &
                                            verbose_column)

          end do
       end do
       
    endif    ! restart file, input file, or other options

  end subroutine glissade_init_therm

!=======================================================================

  subroutine glissade_init_temp_column(temp_init,     &
                                       stagsigma,     &
                                       thck,          &
                                       artm,          &
                                       acab,          &
                                       bheatflx,      &
                                       pmp_offset,    &
                                       temp,          &
                                       verbose_column)

    use glimmer_physcon, only : pi, rhoi, shci, coni

    ! Initialize temperatures in a column based on the value of temp_init.
    ! Four possibilities:
    ! (1) Set ice temperature in column to 0 C (TEMP_INIT_ZERO).
    ! (2) Set ice temperature in column to surface air temperature (TEMP_INIT_ARTM).
    ! (3) Set up a linear temperature profile, with T = artm at the surface and T <= Tpmp
    !     at the bed (TEMP_INIT_LINEAR). 
    !     A config parameter (pmp_offset) controls how far below Tpmp the initial bed temp is set.
    ! (4) Set up a temperature profile based on advective-diffusive balance, with T = artm
    !     at the surface and dT/dz = -F_geo/k at the bed (TEMP_INIT_ADVECTIVE_DIFFUSIVE).
    !     The temperature at each level is capped at the value computed by method (3).
    !
    ! In/out arguments
 
    integer, intent(in) :: temp_init          ! option for temperature initialization

    real(dp), dimension(:), intent(in) ::  &
         stagsigma               ! staggered vertical coordinate; includes layer midpoints, but not top and bottom surfaces

    real(dp), intent(in) :: &
         artm,                &  ! surface air temperature (deg C)
         thck,                &  ! ice thickness (m)
         acab,                &  ! surface mass balance (m/s)
         bheatflx,            &  ! basal heat flux (W/m^2), positive down
         pmp_offset              ! minimum offset of initial Tbed from pressure melting point temperature (deg C)
                                 ! Note: pmp_offset is positive for Tbed < bpmp

    real(dp), dimension(0:), intent(inout) :: &
         temp                    ! ice column temperature (deg C)
                                 ! Note first index of zero

    logical, intent(in) :: verbose_column    ! if true, then write diagnostic info for the column

    ! Local variables and parameters

    real(dp) :: pmptemp_bed                           ! pressure melting point temp at the bed
    real(dp), dimension(size(stagsigma)) :: pmptemp   ! pressure melting point temp thru the column
    real(dp) :: erf_sfc, erf_k                        ! error function evaluated at surface and in layers
    real(dp) :: length_scale                          ! length scale for advective-diffusive balance
                                                      ! length_scale^2 = 2 * k * H / (rhoi * ci * b)
    real(dp) :: temp_advective_diffusive              ! temperature computed based on advective-diffusive balance
    real(dp) :: factor

    integer :: upn                                    ! number of vertical levels (deduced from temp array)
    integer :: k

    upn = size(temp) - 1     ! temperature array has dimension (0:model%general%upn)

    ! Set the temperature in the column

    if (temp_init == TEMP_INIT_ZERO) then        ! set T = 0 C
       
       temp(:) = 0.0d0

    elseif (temp_init == TEMP_INIT_ARTM) then    ! initialize ice-covered areas to the min of artm and 0 C

       temp(:) = min(0.0d0, artm)
       
    elseif (temp_init == TEMP_INIT_LINEAR) then  ! Tsfc = artm, Tbed = Tpmp - pmp_offset, linear profile in between

       if (thck > 0.0d0) then

          temp(0) = min(0.0d0, artm)

          call glissade_pressure_melting_point(thck, pmptemp_bed)
          temp(upn) = pmptemp_bed - pmp_offset

          temp(1:upn-1) = temp(0) + (temp(upn) - temp(0))*stagsigma(:)
                               
          ! Make sure T <= Tpmp - pmp_offset in column interior

          call glissade_pressure_melting_point_column(thck, stagsigma(1:upn-1), pmptemp(1:upn-1))
          temp(1:upn-1) = min(temp(1:upn-1), pmptemp(1:upn-1) - pmp_offset)

       else

          temp(:) = min(0.0d0, artm)

       endif  ! thck > 0

    elseif (temp_init == TEMP_INIT_ADVECTIVE_DIFFUSIVE) then

       ! Where acab > 0, assume a balance between advection and diffusion in the column:
       !
       ! w dT/dz = k/(rhoi*ci) * d2T/dz2
       !
       ! where w = vertical velocity = -b * (z_s - z) / H
       !       H = z_s - z_b = thck
       !       b = surface mass balance = acab
       !
       ! with boundary conditions Tsfc = artm at z = z_s
       !                          dT/dZ|b = -F_geo/k at z = z_b
       !
       ! where F_geo = geothermal heat flux = -bheatflx
       !
       ! This has an exact solution (see Cuffey & Paterson, Physics of Glaciers, 2010, Section 9.5.1):
       !
       ! T(z) = Tsfc + 1/2 * sqrt(pi) * (z*) * dT/dz|b * [erf((z-z_b)/z*) - erf(H/z*)]
       !
       ! where erf(z) = error function = (2/sqrt(pi)) * int_0^z [exp(-y^2)] dy
       !   and z* is a length scale given by z*^2 = (2 * k * H) / (rhoi * ci * b)
       !
       ! The resulting temperature profile is typical of the ice divide, where horizontal advection is small.
       !
       ! The temperature is not allowed to exceed the linear value that would be computed using method (3).
       !
       ! Where acab <= 0, we use the linear profile of method 3.
       ! If acab is not supplied in the input file, then it is initialized to zero,
       ! and the temperature profile will default to linear.
       !
       ! The purpose of this option is to reduce the spinup time by prescribing initial temperatures
       !  that are a good approximation of the steady-state values near the ice divide.
       ! If the initial temperatures are too warm in the upper layers of the ice interior
       !  (as may be the case with a linear profile), we can have thawing where the bed should be frozen.

       if (thck > 0.0d0) then

          ! Initialize the column to a linear profile as for TEMP_INIT_LINEAR
          temp(0) = min(0.0d0, artm)

          call glissade_pressure_melting_point(thck, pmptemp_bed)
          temp(upn) = pmptemp_bed - pmp_offset

          temp(1:upn-1) = temp(0) + (temp(upn) - temp(0))*stagsigma(:)

          ! Make sure T <= Tpmp - pmp_offset in column interior

          call glissade_pressure_melting_point_column(thck, stagsigma(1:upn-1), pmptemp(1:upn-1))
          temp(1:upn-1) = min(temp(1:upn-1), pmptemp(1:upn-1) - pmp_offset)

          if (acab > 0.0d0) then   ! SMB is positive; switch to the advective-diffusive balance

             length_scale = sqrt( (2.0d0*coni*thck)/(rhoi*shci*acab) )
             factor = 0.5d0 * sqrt(pi) * length_scale * (bheatflx/coni)  ! note bheatflx <= 0
 
             !TODO - Evaluate the error function using the Fortran intrinsic erf(x)?
             !       This function is supported in Fortran 2008 (and the Gnu compiler)
             !        but might not be portable to all compilers.

!             erf_sfc = erf(thck/length_scale)  ! might not be portable
             erf_sfc = error_function(thck/length_scale)   ! use hand-coded function instead

             do k = 1, upn-1
!                erf_k = erf ( thck/length_scale * (1.0d0 - stagsigma(k)))  ! might not be portable
                erf_k = error_function( thck/length_scale * (1.0d0 - stagsigma(k)))
                temp_advective_diffusive = temp(0) + factor * (erf_k - erf_sfc)
                temp(k) = min(temp(k), temp_advective_diffusive)
             enddo

             ! bottom surface
             k = upn
             erf_k = 0.0d0
             temp_advective_diffusive = temp(0) + factor * (erf_k - erf_sfc)
             temp(k) = min(temp(k), temp_advective_diffusive)

          endif   ! acab > 0

       else

          temp(:) = min(0.0d0, artm)

       endif  ! thck > 0

       if (verbose_column) then
          print*, 'thck (m), acab (m/yr), artm(deg C) =', thck, acab*scyr, artm
          print*, ' '
          print*, 'Initial temperature profile: k, temperature:'
          do k = 0, upn
             print*, k, temp(k)
          enddo
       endif

    endif  ! temp_init

  end subroutine glissade_init_temp_column

!=======================================================================

  function error_function(x)

    use glimmer_physcon, only : pi

    real(dp) :: x
    real(dp) :: error_function

    ! Approximate the error function erf(x) using an analytic formula from Wikipedia:
    ! 
    ! erf(x) = sgn(x) * sqrt {1 - exp[-x^2 * (4/pi + a*x^2) / (1 + a*x^2)] }
    !
    ! where a = 8*(pi-3) / (3*pi*(4-pi))
    !
    ! According to Wikipedia (21 March 2016):
    ! "This is designed to be very accurate in a neighborhood of 0 and a neighborhood of infinity, 
    !  and the error is less than 0.00035 for all x."
    ! I verified this statement by comparing results to the intrinsic erf(x) provided with the Gnu compiler.

    real(dp) :: a, exp_factor

    a = 8.d0 * (pi - 3.d0) / (3.d0 * pi * (4.d0 - pi))

    exp_factor = -x**2 * (4.d0/pi + a*x**2) / (1.d0 + a*x**2)

    error_function = sqrt(1.d0 - exp(exp_factor))

    if (x < 0.d0) error_function = -error_function

  end function error_function

!=======================================================================

  subroutine glissade_therm_driver(whichtemp,                         &
                                   temp_init,                         &
                                   dttem,                             &
                                   ewn,             nsn,       upn,   &
                                   itest,           jtest,     rtest, &
                                   sigma,           stagsigma,        &
                                   thklim_temp,                       &
                                   thck,            topg,             &
                                   lsrf,            eus,              &
                                   artm,                              &
                                   bheatflx,        bfricflx,         &
                                   dissip,                            &
                                   pmp_threshold,                     &
                                   bwat,                              &
                                   temp,            waterfrac,        &
                                   bpmp,                              &
                                   bmlt)

    ! Calculate the new ice temperature by one of several methods:
    ! (0) set to surface air temperature
    ! (1) standard prognostic temperature solve
    ! (2) hold temperature steady
    ! (3) prognostic solve for enthalpy (a function of temperature and waterfrac)

    ! Note: SI units are used throughout this subroutine

    use glimmer_utils,  only : tridiag
    use glimmer_physcon, only: shci, coni, rhoi, tocnfrz_sfc, dtocnfrz_dh
    use glide_mask
    use glissade_masks, only: glissade_get_masks

    !------------------------------------------------------------------------------------
    ! Input/output arguments
    !------------------------------------------------------------------------------------

    integer, intent(in) ::   &
         whichtemp             ! option for computing temperature

    integer, intent(in) ::   &
         temp_init             ! option for initializing the temperature (used for thin ice)

    integer, intent(in) ::   &
         ewn, nsn, upn,      & ! grid dimensions
         itest, jtest, rtest   ! coordinates of diagnostic point
 
    real(dp), intent(in) ::   &
         dttem,           &! time step for temperature solve (s)
         thklim_temp,     &! minimum ice thickness (m) for temperature calculation
         eus               ! eustatic sea level (m), = 0. by default  

    real(dp), dimension(:), intent(in) ::   &
         sigma,           &! vertical coordinate, located at layer interfaces
         stagsigma         ! staggered vertical coordinate, located at the center of each layer

    real(dp), dimension(:,:), intent(in) ::  &
         thck,            &! ice thickness (m)
         topg,            &! elevation of basal topography (m)
         lsrf,            &! elevation of lower ice surface (m)
         artm,            &! surface air temperature (deg C)
         bwat,            &! basal water depth (m)
         bheatflx,        &! geothermal flux (W m-2), positive down
         bfricflx          ! basal friction heat flux (W m-2), >= 0
         
    real(dp), dimension(:,:,:), intent(in) ::  &
         dissip            ! interior heat dissipation (deg/s)

    real(dp), intent(in) ::  &
         pmp_threshold     ! bed is assumed thawed where Tbed >= pmptemp - pmp_threshold (deg C)

    real(dp), dimension(0:,:,:), intent(out) ::  &
         temp              ! ice temperature (deg C)

    real(dp), dimension(:,:,:), intent(out) ::  &
         waterfrac         ! internal water fraction (unitless)

    real(dp), dimension(:,:), intent(out) ::  &
         bpmp,            &! basal pressure melting point temperature (deg C)
         bmlt              ! basal melt rate (m/s), > 0 for melting

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    character(len=150) :: message

    real(dp), dimension(0:upn,ewn,nsn) ::  &
         enthalpy          ! specific enthalpy (J m-3)

    real(dp), dimension(upn+1) :: subd, diag, supd, rhsd   ! matrix coefficients

    real(dp),dimension(0:upn) :: prevtemp   ! previous temperature in column

    real(dp) ::                &
         einit, efinal,        &! initial and final internal energy
         delta_e,              &! net energy input to ice
         dTtop, dTbot,         &! temperature differences
         denth_top, denth_bot, &! enthalpy differences
         maxtemp, mintemp,     &! max and min temps in column
         depth                  ! depth at base of ice shelf (m)

    real(dp), dimension(1:upn) :: alpha_enth   ! diffusivity at interfaces (m2/s) for enthalpy solver
                                               ! = coni / (rhoi*shci) for cold ice

    integer, dimension(ewn,nsn) ::  &
         ice_mask,      &! = 1 where ice temperature is computed (thck > thklim_temp), else = 0
         floating_mask, &! = 1 where ice is present and floating, else = 0
         ocean_mask      ! = 1 where topg is below sea level and ice is absent

    !TODO - ucondflx may be needed for coupling; make it an output argument?
    real(dp), dimension(ewn,nsn) ::  &
         ucondflx,     & ! conductive heat flux (W/m^2) at upper sfc (positive down)
         lcondflx,     & ! conductive heat flux (W/m^2) at lower sfc (positive down)
         dissipcol       ! total heat dissipation rate (W/m^2) in column (>= 0)

    integer :: ew, ns, up
    integer :: k

    logical :: lstop = .false.   ! flag for energy conservation error
    integer :: istop, jstop      ! local location of energy conservation error
    integer :: istop_global, jstop_global    ! global location of energy conservation error

    logical :: verbose_column    ! if true, then write diagnostic info for the column

    !------------------------------------------------------------------------------------
    ! Compute the new temperature profile in each column
    !------------------------------------------------------------------------------------

    ! initialize some fluxes
    lcondflx(:,:) = 0.d0
    ucondflx(:,:) = 0.d0
    dissipcol(:,:) = 0.d0
    
    ! Compute masks: ice_mask = 1 where thck > thklim_temp;
    !                floating_mask = 1 where ice is present (thck > thklim_temp) and floating;
    !                ocean_mask = 1 where topg is below sea level and thck <= thklim_temp

    call glissade_get_masks(ewn,           nsn,            &
                            thck,          topg,           &
                            eus,           thklim_temp,    &
                            ice_mask,                      &
                            floating_mask = floating_mask, &
                            ocean_mask = ocean_mask)
      
    ! Compute basal pressure melting point temperature
    ! (needed for temperature/enthalpy calculation below, and also needed later for
    !  certain basal sliding laws)

    do ns = 1, nsn
       do ew = 1, ewn
          if (ice_mask(ew,ns) == 1) then    ! thck > thklim_temp
             call glissade_pressure_melting_point(thck(ew,ns), bpmp(ew,ns))
          else
             bpmp(ew,ns) = 0.d0
          endif
       enddo
    enddo

    select case(whichtemp)

    case(TEMP_SURFACE_AIR_TEMP)  ! Set column to surface air temperature ------------------

       do ns = 1, nsn
          do ew = 1, ewn
             temp(:,ew,ns) = min(0.0d0, artm(ew,ns))
          end do
       end do

    case(TEMP_STEADY)   ! do nothing

    case(TEMP_PROGNOSTIC, TEMP_ENTHALPY)  ! Local column calculation

       ! No horizontal or vertical advection; vertical diffusion and strain heating only.
       ! Temperatures are vertically staggered relative to velocities.  
       ! That is, the temperature is defined at the midpoint of each layer 
       ! (and at the top and bottom surfaces).

       !Note: Interior heat dissipation used to be calculated here.
       !      Now it is computed at the end of the previous time step, after solving for velocity.
       !
       !      Basal frictional heating also was calculated here.
       !      Now it is computed in the velocity solver (for Glissade) or just after
       !       the velocity solution (for Glam).

       ! Set error flag
       lstop = .false.

       ! loop over cells
       do ns = 1, nsn
       do ew = 1, ewn

          if (verbose_therm .and. this_rank==rtest .and. ew==itest .and. ns==jtest) then
             verbose_column = .true.
          else
             verbose_column = .false.
          endif

          if (ice_mask(ew,ns) == 1) then   ! thck > thklim_temp

             ! Set surface temperature

             temp(0,ew,ns) = min(0.d0, artm(ew,ns))

             ! For floating ice, set the basal temperature to the freezing temperature of seawater
             ! Values based on Ocean Water Freezing Point Calculator with S = 35 PSU
             if (floating_mask(ew,ns) == 1) then
                depth = thck(ew,ns) * (rhoi/rhow)
                temp(upn,ew,ns) = tocnfrz_sfc + dtocnfrz_dh * depth
             endif

             if (whichtemp == TEMP_ENTHALPY) then

                ! Given temperature and waterfrac, compute enthalpy (dimension 0:upn)
                ! Assume waterfrac = 0 at upper and lower surfaces.

                call glissade_temp2enth(stagsigma(1:upn-1),                            &
                                        temp(0:upn,ew,ns),  waterfrac(1:upn-1,ew,ns),  &
                                        thck(ew,ns),        enthalpy(0:upn,ew,ns))

                if (verbose_column) then
                   print*, ' '
                   print*, 'Before prognostic enthalpy, i, j =', ew, ns
                   print*, 'thck =', thck(ew,ns)
                   print*, 'Temp, waterfrac, enthalpy/(rhoi*shci):'
                   k = 0
                   print*, k, temp(k,ew,ns), 0.d0, enthalpy(k,ew,ns)/(rhoi*shci)
                   do k = 1, upn-1
                      print*, k, temp(k,ew,ns), waterfrac(k,ew,ns), enthalpy(k,ew,ns)/(rhoi*shci)
                   enddo
                   k = upn
                   print*, k, temp(k,ew,ns), 0.d0, enthalpy(k,ew,ns)/(rhoi*shci)
                endif

                ! compute initial internal energy in column (for energy conservation check)
                einit = 0.0d0
                do up = 1, upn-1
                   einit = einit + enthalpy(up,ew,ns) * (sigma(up+1) - sigma(up))
                enddo
                einit = einit * thck(ew,ns)

                ! compute matrix elements using enthalpy gradient method

                call glissade_enthalpy_matrix_elements(dttem,                     &
                                                       upn,         stagsigma,    &
                                                       subd,        diag,         &
                                                       supd,        rhsd,         &
                                                       dups(:,:),                 &
                                                       floating_mask(ew,ns),      &
                                                       thck(ew,ns),               &
                                                       bpmp(ew,ns),               &
                                                       temp(:,ew,ns),             &  !TODO - 0:upn?
                                                       waterfrac(:,ew,ns),        &
                                                       enthalpy(0:upn,ew,ns),     &
                                                       dissip(:,ew,ns),           &
                                                       bheatflx(ew,ns),           &
                                                       bfricflx(ew,ns),           &
                                                       bwat(ew,ns),               &
                                                       pmp_threshold,             &
                                                       alpha_enth,                &
                                                       verbose_column)
                
                !WHL - debug
                if (verbose_column) then
                   print*, ' '
                   print*, 'After matrix elements, i, j =', ew, ns
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
                             enthalpy(0:upn,ew,ns), &
                             rhsd(1:upn+1))
               
                ! Compute conductive fluxes = (alpha/H * denth/dsigma) at upper and lower surfaces; positive down.
                ! Here alpha = coni / (rhoi*shci) for cold ice, with a smaller value for temperate ice.
                ! Assume fully implicit backward Euler time step.
                ! Note: These fluxes should be computed before calling glissade_enth2temp (which might reset the bed enthalpy).

                denth_top = enthalpy(1,ew,ns) - enthalpy(0,ew,ns)
                denth_bot = enthalpy(upn,ew,ns) - enthalpy(upn-1,ew,ns)
                
                ucondflx(ew,ns) = -alpha_enth(1)  /thck(ew,ns) * denth_top/( stagsigma(1))
                lcondflx(ew,ns) = -alpha_enth(upn)/thck(ew,ns) * denth_bot/(1.d0 - stagsigma(upn-1))
                                              
                ! convert enthalpy back to temperature and water content
                call glissade_enth2temp(stagsigma(1:upn-1),                          &
                                        thck(ew,ns),       enthalpy(0:upn,ew,ns),    &
                                        temp(0:upn,ew,ns), waterfrac(1:upn-1,ew,ns))

                if (verbose_column) then
                   print*, ' '
                   print*, 'After prognostic enthalpy, i, j =', ew, ns
                   print*, 'thck =', thck(ew,ns)
                   print*, 'Temp, waterfrac, enthalpy/(rhoi*shci):'
                   k = 0
                   print*, k, temp(k,ew,ns), 0.d0, enthalpy(k,ew,ns)/(rhoi*shci)
                   do k = 1, upn-1
                      print*, k, temp(k,ew,ns), waterfrac(k,ew,ns), enthalpy(k,ew,ns)/(rhoi*shci)
                   enddo
                   k = upn
                   print*, k, temp(k,ew,ns), 0.d0, enthalpy(k,ew,ns)/(rhoi*shci)
                endif
                
                ! compute the final internal energy

                efinal = 0.0d0
                do up = 1, upn-1
                   efinal = efinal + enthalpy(up,ew,ns) * (sigma(up+1) - sigma(up) )
                enddo
                efinal = efinal * thck(ew,ns)

             else   ! whichtemp = TEMP_PROGNOSTIC

                if (verbose_column) then
                   print*, ' '
                   print*, 'Before prognostic temp, i, j =', ew, ns
                   print*, 'thck =', thck(ew,ns)
                   print*, 'Temp:'
                   do k = 0, upn
                      print*, k, temp(k,ew,ns)
                   enddo
                endif
                
                ! compute initial internal energy in column (for energy conservation check)
                einit = 0.0d0
                do up = 1, upn-1
                   einit = einit + temp(up,ew,ns) * (sigma(up+1) - sigma(up))
                enddo
                einit = einit * rhoi * shci * thck(ew,ns)
                
                ! compute matrix elements
                !TODO - Pass dups?
                call glissade_temperature_matrix_elements(dttem,                 &
                                                          upn,     stagsigma,    &
                                                          subd,    diag,         &
                                                          supd,    rhsd,         &            
                                                          floating_mask(ew,ns),  &
                                                          thck(ew,ns),           &
                                                          bpmp(ew,ns),           &
                                                          temp(:,ew,ns),         &
                                                          dissip(:,ew,ns),       &
                                                          bheatflx(ew,ns),       &
                                                          bfricflx(ew,ns),       &
                                                          bwat(ew,ns),           &
                                                          pmp_threshold)
                
                if (verbose_column) then
                   print*, 'After matrix elements, i, j =', ew,ns
                   print*, 'k, subd, diag, supd, rhsd:'
                   do k = 1, upn+1
                      print*, k, subd(k), diag(k), supd(k), rhsd(k)
                   enddo
                endif

                prevtemp(:) = temp(:,ew,ns)

                ! solve the tridiagonal system

                ! Note: Temperature is indexed from 0 to upn, with indices 1 to upn-1 colocated
                !  with stagsigma values of the same index.
                ! However, the matrix elements are indexed 1 to upn+1, with the first row
                !  corresponding to the surface temperature, temp(0,:,:).

                call tridiag(subd(1:upn+1), &
                             diag(1:upn+1), &
                             supd(1:upn+1), &
                             temp(0:upn,ew,ns), &
                             rhsd(1:upn+1))

                ! conductive flux = (k/H * dT/dsigma) at upper and lower surfaces; positive down

                if (crank_nicolson) then
                   ! average temperatures between start and end of timestep
                   dTtop = 0.5d0 * (temp(1,ew,ns) - temp(0,ew,ns) + prevtemp(1) - prevtemp(0))
                   dTbot = 0.5d0 * (temp(upn,ew,ns) - temp(upn-1,ew,ns) + prevtemp(upn) - prevtemp(upn-1))
                else    ! fully implicit
                   ! use temperatures at end of timestep
                   dTtop = temp(1,ew,ns) - temp(0,ew,ns)
                   dTbot = temp(upn,ew,ns) - temp(upn-1,ew,ns)
                endif

                ucondflx(ew,ns) = (-coni/thck(ew,ns) ) * dTtop / (stagsigma(1))
                lcondflx(ew,ns) = (-coni/thck(ew,ns) ) * dTbot / (1.d0 - stagsigma(upn-1))

                if (verbose_column) then
                   print*, ' '
                   print*, 'After prognostic temp, i, j =', ew, ns
                   print*, 'Temp:'
                   do k = 0, upn
                      print*, k, temp(k,ew,ns)
                   enddo
                endif
                
                ! compute the final internal energy
                
                efinal = 0.0d0
                do up = 1, upn-1
                   efinal = efinal + temp(up,ew,ns) * (sigma(up+1) - sigma(up))
                enddo
                efinal = efinal * rhoi*shci * thck(ew,ns)
                
             endif   ! whichtemp

             ! Compute total dissipation rate in column (W/m^2)

             dissipcol(ew,ns) = 0.0d0
             do up = 1, upn-1
                dissipcol(ew,ns) = dissipcol(ew,ns) + dissip(up,ew,ns) * (sigma(up+1) - sigma(up))
             enddo
             dissipcol(ew,ns) = dissipcol(ew,ns) * thck(ew,ns)*rhoi*shci

             ! Verify that the net input of energy into the column is equal to the change in
             ! internal energy.  

             delta_e = (ucondflx(ew,ns) - lcondflx(ew,ns) + dissipcol(ew,ns)) * dttem

             if (abs((efinal-einit-delta_e)/dttem) > 1.0d-8) then

                if (verbose_column) then
                   print*, 'Ice thickness:', thck(ew,ns)
                   print*, 'thklim_temp:', thklim_temp
                   print*, ' '
                   print*, 'Interior fluxes:'
                   print*, 'ftop (pos up)=', -ucondflx(ew,ns) 
                   print*, 'fbot (pos up)=', -lcondflx(ew,ns)
                   print*, 'fdissip =',       dissipcol(ew,ns)
                   print*, 'Net flux =', delta_e/dttem
                   print*, ' '
                   print*, 'delta_e =', delta_e
                   print*, 'einit =',  einit
                   print*, 'efinal =', efinal
                   print*, 'einit + delta_e =', einit + delta_e
                   print*, ' '
                   print*, 'Energy imbalance =', efinal - einit - delta_e
                   print*, ' '
                   print*, 'Basal fluxes:'
                   print*, 'bfricflx =', bfricflx(ew,ns)
                   print*, 'bheatflx =', -bheatflx(ew,ns)
                   print*, 'flux for bottom melting =', bfricflx(ew,ns) - bheatflx(ew,ns) + lcondflx(ew,ns)
                endif   ! verbose_column

                lstop = .true.
                istop = ew
                jstop = ns
            
             endif  ! energy conservation error

          endif  ! thck > thklim_temp
       end do    ! ew
       end do    ! ns

       if (lstop) then
          call parallel_globalindex(istop, jstop, istop_global, jstop_global)
          call broadcast(istop_global, proc=this_rank)
          call broadcast(istop_global, proc=this_rank)
          print*, 'ERROR: Energy not conserved in glissade_therm, rank, i, j =', this_rank, istop, jstop
          print*, 'Global i, j:', istop_global, jstop_global
          write(message,*) 'ERROR: Energy not conserved in glissade_therm' //  &
               ' (could be caused by a CFL violation), global i, j =', istop_global, jstop_global
          call write_log(message,GM_FATAL)
       endif

       ! Set the temperature of thin ice columns where the temperature is not computed above.
       ! NOTE: This procedure will maintain a sensible temperature profile for thin ice,
       !        but in general does *not* conserve energy.
       !       To conserve energy, we need either thklim_temp = 0, or some additional
       !        energy accounting and correction.

       do ns = 1, nsn
          do ew = 1, ewn

             if (thck(ew,ns) <= thklim_temp) then

                temp(:,ew,ns) = min(0.0d0, artm(ew,ns))

             endif

          end do
       end do

    end select   ! whichtemp

    ! Calculate the basal melt rate for grounded ice.
    ! Basal melt for floating ice is handled in a separate module.
    ! Note: This calculation includes internal melting for all ice, either grounded or floating.
    !       For the prognostic temperature scheme, temperatures above the pressure melting point 
    !        are reset to Tpmp, with excess heat contributing to basal melt.
    !       For the enthalpy scheme, internal meltwater in excess of the prescribed maximum
    !        fraction (0.01 by default) is drained to the bed.

    call glissade_basal_melting_ground(whichtemp,                         &
                                       dttem,                             &
                                       ewn,              nsn,             &
                                       upn,                               &
                                       sigma,            stagsigma,       &
                                       ice_mask,         floating_mask,   &
                                       thck,             bpmp,            &
                                       temp,             waterfrac,       &
                                       enthalpy,                          &
                                       bfricflx,         bheatflx,        &
                                       lcondflx,         bwat,            &
                                       pmp_threshold,                     &
                                       bmlt)

    ! Check for temperatures that are physically unrealistic.
    ! Thresholds are set at the top of this module.

    do ns = 1, nsn
       do ew = 1, ewn

          maxtemp = maxval(temp(:,ew,ns))
          mintemp = minval(temp(:,ew,ns))
          
          if (maxtemp > maxtemp_threshold) then
             write(message,*) 'maxtemp > 0: i, j, maxtemp =', ew, ns, maxtemp
             call write_log(message,GM_FATAL)
          endif
          
          if (mintemp < mintemp_threshold) then
             !uncommment these lines to get more info
!             print*, 'thck =', thck(ew,ns)
!             print*, 'temp:'
!             do k = 1, upn
!                print*, k, temp(k,ew,ns)
!             enddo
             write(message,*) 'mintemp < mintemp_threshold: i, j, mintemp =', ew, ns, mintemp
             call write_log(message,GM_FATAL)
          endif
          
       enddo
    enddo

  end subroutine glissade_therm_driver

!=======================================================================

  subroutine glissade_temperature_matrix_elements(dttem,                        &
                                                  upn,          stagsigma,      &
                                                  subd,         diag,           &
                                                  supd,         rhsd,           &
                                                  floating_mask,                &
                                                  thck,         bpmp,           &
                                                  temp,         dissip,         &
                                                  bheatflx,     bfricflx,       &
                                                  bwat,         pmp_threshold)

    ! compute matrix elements for the tridiagonal solve

    use glimmer_physcon,  only : rhoi, grav, coni

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    !        whereas temperature is indexed from 0 to upn.
    !       The first row of the matrix is the equation for temp(0,ew,ns),
    !        the second row is the equation for temp(1,ew,ns), and so on.

    real(dp), intent(in) :: dttem       ! time step (s)
    integer, intent(in) :: upn          ! number of layer interfaces
    real(dp), dimension(upn-1), intent(in) :: stagsigma    ! sigma coordinate at temp nodes

    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd

    integer, intent(in) :: floating_mask
    real(dp), intent(in) ::  thck       ! ice thickness (m)
    real(dp), intent(in) ::  bpmp       ! basal pressure melting point temperature (deg C)

    real(dp), dimension(0:upn), intent(in) :: temp       ! ice temperature (deg C)
    real(dp), dimension(upn-1), intent(in) :: dissip     ! interior heat dissipation (deg/s)
    real(dp), intent(in) :: bheatflx    ! geothermal flux (W m-2), positive down
    real(dp), intent(in) :: bfricflx    ! basal friction heat flux (W m-2), >= 0
    real(dp), intent(in) :: bwat        ! basal water thickness (m)

    real(dp), intent(in) :: &
         pmp_threshold       ! bed is assumed thawed where Tbed >= pmptemp - pmp_threshold (deg C)

    ! local variables

    real(dp) :: fact
    real(dp) :: dsigbot      ! bottom layer thicknes in sigma coords

    ! Compute subdiagonal, diagonal, and superdiagonal matrix elements

    ! upper boundary: set to surface temperature

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = temp(0)

    ! ice interior, layers 1:upn-1  (matrix elements 2:upn)

    if (crank_nicolson) then  ! C-N can lead to oscillations in thin ice; currently deprecated

       fact = dttem * coni / (2.d0 * rhoi*shci) / thck**2
       subd(2:upn) = -fact * dups(1:upn-1,1)
       supd(2:upn) = -fact * dups(1:upn-1,2)
       diag(2:upn) = 1.0d0 - subd(2:upn) - supd(2:upn)
       rhsd(2:upn) =  temp(1:upn-1) * (2.0d0 - diag(2:upn)) &
                    - temp(0:upn-2) * subd(2:upn) &
                    - temp(2:upn  ) * supd(2:upn) & 
                    + dissip(1:upn-1)

    else   ! fully implicit

       fact = dttem * coni / (rhoi*shci) / thck**2
       subd(2:upn) = -fact * dups(1:upn-1,1)
       supd(2:upn) = -fact * dups(1:upn-1,2)
       diag(2:upn) = 1.0d0 - subd(2:upn) - supd(2:upn)
       rhsd(2:upn) = temp(1:upn-1) + dissip(1:upn-1)*dttem

    endif    ! crank_nicolson

    ! basal boundary:
    ! For floating ice, the basal temperature is held constant.
    ! For grounded ice, a heat flux is applied. The bed temperature is held at bpmp if it is already
    !  at or near bpmp or if basal water is present; else the bed temperature is computed based on
    !  a balance of fluxes.

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

       if (temp(upn) >= bpmp - pmp_threshold .or. bwat > 0.0d0) then

          ! hold basal temperature at pressure melting point

          supd(upn+1) = 0.0d0
          subd(upn+1) = 0.0d0
          diag(upn+1) = 1.0d0
          rhsd(upn+1) = bpmp

       else   ! frozen at bed
              ! maintain balance of heat sources and sinks
              ! (conductive flux, geothermal flux, and basal friction)

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
          rhsd(upn+1) = (bfricflx - bheatflx) * dsigbot*thck / coni

       endif   ! melting or frozen

    end if     ! floating or grounded

  end subroutine glissade_temperature_matrix_elements

!=======================================================================

  subroutine glissade_enthalpy_matrix_elements(dttem,                       &
                                               upn,       stagsigma,        &
                                               subd,      diag,             &
                                               supd,      rhsd,             &
                                               dups,      floating_mask,    &
                                               thck,      bpmp,             &
                                               temp,      waterfrac,        &
                                               enthalpy,  dissip,           &
                                               bheatflx,  bfricflx,         &
                                               bwat,      pmp_threshold,    &
                                               alpha_enth,                  &
                                               verbose_column_in)

    ! solve for tridiagonal entries of sparse matrix

    use glimmer_physcon,  only : rhoi, shci, lhci, rhow, coni

    ! Note: Matrix elements (subd, supd, diag, rhsd) are indexed from 1 to upn+1,
    ! whereas temperature/enthalpy is indexed from 0 to upn.
    ! The first row of the matrix is the equation for enthalpy(0),
    ! the last row is the equation for enthalpy(upn), and so on.

    ! input/output variables
    real(dp), intent(in) :: dttem       ! time step (s)
    integer, intent(in) :: upn          ! number of layer interfaces
    real(dp), dimension(upn-1), intent(in) :: stagsigma    ! sigma coordinate at temp/enthalpy nodes

    real(dp), dimension(:,:), intent(in) :: dups   ! vertical grid quantities
    real(dp), dimension(:), intent(out) :: subd, diag, supd, rhsd

    integer, intent(in) :: floating_mask
    real(dp), intent(in) :: thck        ! ice thickness (m)
    real(dp), intent(in) :: bpmp        ! basal pressure melting point temperature (deg C)

    real(dp), dimension(0:upn), intent(in) :: temp       ! temperature (deg C)
    real(dp), dimension(upn-1), intent(in) :: waterfrac  ! water fraction (unitless)
    real(dp), dimension(0:upn), intent(in) :: enthalpy   ! specific enthalpy (J/m^3)
    real(dp), dimension(upn-1), intent(in) :: dissip     ! interior heat dissipation (deg/s)

    real(dp), intent(in) :: bheatflx   ! geothermal flux (W m-2), positive down
    real(dp), intent(in) :: bfricflx   ! basal friction heat flux (W m-2), >= 0
    real(dp), intent(in) :: bwat       ! basal water thickness (m)

    real(dp), intent(in) :: &
         pmp_threshold           ! bed is assumed thawed where Tbed >= pmptemp - pmp_threshold (deg C)

    real(dp), dimension(:), intent(out) :: alpha_enth  ! half-node diffusivity (m^2/s) for enthalpy
	                                               ! located halfway between temperature points

    logical, intent(in), optional :: verbose_column_in   ! if true, print debug statements for this column
    
    ! local variables
    real(dp) :: dsigbot  ! bottom layer thicknes in sigma coords.
    real(dp) :: alphai ! cold ice diffusivity
    real(dp) :: alpha0 ! temperate ice diffusivity
    real(dp) :: fact ! coefficient in tridiag
    integer  :: up
    real(dp), dimension(1:upn-1) :: pmptemp    ! pressure melting point temp in interior (deg C)
    real(dp), dimension(0:upn) :: enth_T       ! temperature part of specific enthalpy (J/m^3)
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
	
    ! find pmptemp for this column
    ! Note: The basal pmp temperature (bpmp) is an input argument
    call glissade_pressure_melting_point_column(thck, stagsigma(1:upn-1), pmptemp(1:upn-1))

    !WHL - debug                                                                                                                       
    if (verbose_column) then
       print*, ' '
       print*, 'Computing enthalpy matrix elements'
       print*, 'k, temp, wfrac, enthalpy/(rhoi*ci), pmpt:'
       up = 0
       print*, up, temp(up), 0.d0, enthalpy(up)/(rhoi*shci)
       do up = 1, upn-1
          print*, up, temp(up), waterfrac(up), &
               enthalpy(up)/(rhoi*shci), pmptemp(up)
       enddo
       up = upn
       print*, up, temp(up), 0.d0, enthalpy(up)/(rhoi*shci), bpmp
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
    ! Assume backward Euler time stepping
    
    ! upper boundary: set to surface air temperature*rhoi*shci
    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0
    rhsd(1) = min(0.0d0,temp(0)) * rhoi*shci
  
    ! ice interior, layers 1:upn-1  (matrix elements 2:upn)

    fact = dttem / thck**2

    subd(2:upn) = -fact * alpha_enth(1:upn-1) * dups(1:upn-1,1)                                
    supd(2:upn) = -fact * alpha_enth(2:upn) * dups(1:upn-1,2)                                
    diag(2:upn) = 1.0d0 - subd(2:upn) - supd(2:upn)                                
    rhsd(2:upn) = enthalpy(1:upn-1) + dissip(1:upn-1)*dttem * rhoi * shci
                              
    ! BDM I'm assuming that dissip has units of phi/rhoi/shci.
    ! For an enthalpy calc, we want just phi, hence dissip * rhoi * shci
	
    ! basal boundary:
    ! For floating ice, the basal temperature is held constant.
    ! For grounded ice, a heat flux is applied. The bed temperature is held at bpmp if it is already
    !  at or near bpmp or if basal water is present; else the bed temperature is computed based on 
    !  a balance of fluxes.

    !NOTE: This lower BC is different from the one in glide_temp.
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
          print*, 'temp(upn), pmptemp(upn):', temp(up), bpmp
       endif

    ! Positive-Thickness Basal Temperate Boundary Layer

    !WHL - Not sure whether this condition is ideal.
    !      It implies that the enthalpy at the bed (upn) = enthalpy in layer (upn-1). 
       if (temp(upn-1) >=  pmptemp(upn-1) - pmp_threshold .or. bwat > 0.0d0) then
       
          subd(upn+1) = -1.0d0
          supd(upn+1) =  0.0d0 
          diag(upn+1) = 1.0d0 
          rhsd(upn+1) = 0.0d0

          !WHL - debug
          if (verbose_column) then
             print*, 'basal BC: branch 1 (finite-thck BL)'
          endif

       !Zero-Thickness Basal Temperate Boundary Layer
       elseif (temp(upn) >= bpmp - pmp_threshold) then  ! melting
          
          ! hold basal temperature at pressure melting point
          supd(upn+1) = 0.0d0
          subd(upn+1) = 0.0d0
          diag(upn+1) = 1.0d0
          rhsd(upn+1) = bpmp * rhoi * shci
          
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
          rhsd(upn+1) = (bfricflx - bheatflx) * dsigbot*thck * rhoi*shci/coni
          ! BDM temp approach should work out to be dT/dsigma, so enthalpy approach
          ! should just need dT/dsigma * rhoi * shci for correct units

       endif   ! melting or frozen
       
    end if     ! floating or grounded

  end subroutine glissade_enthalpy_matrix_elements

!=======================================================================

  subroutine glissade_basal_melting_ground(whichtemp,                         &
                                           dttem,                             &
                                           ewn,              nsn,             &
                                           upn,                               &
                                           sigma,            stagsigma,       &
                                           ice_mask,         floating_mask,   &
                                           thck,             bpmp,            &
                                           temp,             waterfrac,       &
                                           enthalpy,                          &
                                           bfricflx,         bheatflx,        &
                                           lcondflx,         bwat,            &
                                           pmp_threshold,                     &
                                           bmlt_ground)

    ! Compute the rate of basal melting for grounded ice.
    !
    ! For the standard prognostic temperature scheme, any internal temperatures 
    !  above the pressure melting point are reset to Tpmp.  Excess energy 
    !  is applied toward melting with immediate drainage to the bed. 
    ! For the enthalpy scheme, any meltwater in excess of the maximum allowed
    !  meltwater fraction (0.01 by default) is drained to the bed.

    use glimmer_physcon, only: shci, rhoi, lhci

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: whichtemp                  ! temperature method (TEMP_PROGNOSTIC or TEMP_ENTHALPY)

    real(dp), intent(in) :: dttem                     ! time step (s)

    integer, intent(in) :: ewn, nsn, upn              ! grid dimensions

    real(dp), dimension(upn),    intent(in) :: sigma           ! vertical sigma coordinate
    real(dp), dimension(upn-1),  intent(in) :: stagsigma       ! staggered vertical coordinate for temperature

    real(dp), dimension(0:,:,:), intent(inout) :: temp         ! temperature (deg C)
    real(dp), dimension(:,:,:),  intent(inout) :: waterfrac    ! water fraction
    real(dp), dimension(0:,:,:), intent(in) :: enthalpy        ! enthalpy

    real(dp), dimension(:,:),    intent(in) :: &
         thck,                 & ! ice thickness (m)
         bpmp,                 & ! basal pressure melting point temperature (deg C)
         bfricflx,             & ! basal frictional heating flux (W m-2), >= 0
         bheatflx,             & ! geothermal heating flux (W m-2), positive down
         lcondflx,             & ! heat conducted from ice interior to bed (W m-2), positive down
         bwat                    ! depth of basal water (m)

    integer, dimension(:,:), intent(in) ::  &
         ice_mask,             & ! = 1 where ice is present (thck > thklim_temp), else = 0
         floating_mask           ! = 1 where ice is present and floating, else = 0

    real(dp), intent(in) :: &
         pmp_threshold           ! bed is assumed thawed where Tbed >= pmptemp - pmp_threshold (deg C)

    real(dp), dimension(:,:), intent(out):: &
         bmlt_ground             ! basal melt rate for grounded ice (m/s)

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    real(dp), dimension(ewn,nsn) :: &
         bflx_mlt             ! heat flux available for basal melting (W/m^2)

    integer :: up, ew, ns
    real(dp), dimension(upn-1)  :: pmptemp   ! pressure melting point temp in ice interior
    real(dp) :: layer_thck    ! layer thickness (m)
    real(dp) :: melt_energy   ! energy available for internal melting (J/m^2)
    real(dp) :: internal_melt_rate   ! internal melt rate, transferred to bed (m/s)
    real(dp) :: melt_fact     ! factor for bmlt calculation
    real(dp) :: hmlt          ! melt thickness associated with excess meltwater (m)

    real(dp), parameter :: max_waterfrac = 0.01d0   ! maximum allowed water fraction
                                                    ! excess water drains to the bed
    real(dp), parameter :: eps11 = 1.d-11     ! small number

    bmlt_ground(:,:) = 0.0d0

    ! Compute the heat flux available to melt grounded ice
    ! The basal friction term is computed above in subroutine glissade_calcbfric,
    !  or in the Glissade velocity solver.
    ! Note: bflx_mlt > 0 for melting, < 0 for freeze-on
    !       bfricflx >= 0 by definition
    !       bheatflx is positive down, so usually bheatflx < 0 (with negative values contributing to melt)
    !       lcondflx is positive down, so lcondflx < 0 for heat flowing from the bed toward the surface
    !
    !       This equation allows for freeze-on (bmlt_ground < 0) if the conductive term
    !        (lcondflx, positive down) is carrying enough heat away from the boundary.  
    !       But freeze-on requires a local water supply, bwat > 0.
    !       When bwat = 0, we reset the bed temperature to a value slightly below the melting point.

    bflx_mlt(:,:) = bfricflx(:,:) + lcondflx(:,:) - bheatflx(:,:)  ! W/m^2
    
    ! bflx_mlt might be slightly different from zero because of rounding errors; if so, then zero out
    where (abs(bflx_mlt) < eps11)
       bflx_mlt = 0.d0
    endwhere

    ! Compute the basal melt rate for grounded ice
    ! Note: If the temperature/enthalpy is not prognosed, then there is no basal melting for grounded ice.

    if (whichtemp == TEMP_ENTHALPY) then

       do ns = 1, nsn
          do ew = 1, ewn

             !TODO - For the enthalpy scheme, deal with the rare case that the bottom layer melts completely
             !       and overlying layers with a different enthalpy also melt.

             if (ice_mask(ew,ns) == 1 .and. floating_mask(ew,ns) == 0) then   ! ice is present and grounded
                bmlt_ground(ew,ns) = bflx_mlt(ew,ns) / (lhci*rhoi - enthalpy(upn,ew,ns))  !TODO - Use enthalpy in layer upn-1? 
             endif

             ! Add internal melting associated with waterfrac > max_waterfrac (1%)
             ! Note: It is possible to have internal melting for floating ice.
             !       Rather than have a separate calculation in subroutine glissade_basal_melting_float,
             !        internal melting for all ice (both grounded and floating) is computed here. 

             if (ice_mask(ew,ns) == 1) then  ! ice is present

                !TODO - Any correction for rhoi/rhow here?  Or melting ice that is already partly melted?
                do up = 1, upn-1
                   if (waterfrac(up,ew,ns) > max_waterfrac) then

                      ! compute melt rate associated with excess water
                      hmlt = (waterfrac(up,ew,ns) - max_waterfrac) * thck(ew,ns) * (sigma(up+1) - sigma(up))  ! m
                      internal_melt_rate = hmlt / dttem          ! m/s

                      ! reset waterfrac to max value
                      waterfrac(up,ew,ns) = max_waterfrac

                      ! transfer meltwater to the bed
                      bmlt_ground(ew,ns) = bmlt_ground(ew,ns) + internal_melt_rate      ! m/s

                   endif  ! waterfrac > max_waterfrac
                enddo   ! up

             endif   ! ice is present

          enddo   ! ew
       enddo   ! ns

    elseif (whichtemp == TEMP_PROGNOSTIC) then

       do ns = 1, nsn
          do ew = 1, ewn

             if (ice_mask(ew,ns) == 1 .and. floating_mask(ew,ns) == 0) then   ! ice is present and grounded
                bmlt_ground(ew,ns) = bflx_mlt(ew,ns) / (lhci*rhoi)   ! m/s
             endif

             ! Add internal melting associated with T > Tpmp
             ! Note: It is possible to have internal melting for floating ice.
             !       Rather than have a separate calculation in subroutine glissade_basal_melting_float,
             !        internal melting for all ice (both grounded and floating) is computed here. 

             if (ice_mask(ew,ns) == 1) then  ! ice is present

                call glissade_pressure_melting_point_column(thck(ew,ns), stagsigma(:), pmptemp(:))

                do up = 1, upn-1
                   if (temp(up,ew,ns) > pmptemp(up)) then

                      ! compute melt rate associated with T > Tpmp
                      layer_thck = thck(ew,ns) * (sigma(up+1) - sigma(up))  ! m
                      melt_energy = rhoi * shci * (temp(up,ew,ns) - pmptemp(up)) * layer_thck         ! J/m^2
                      internal_melt_rate = melt_energy / (rhoi * lhci * dttem)  ! m/s

                      ! reset T to Tpmp
                      temp(up,ew,ns) = pmptemp(up)

                      ! transfer internal melting to the bed
                      bmlt_ground(ew,ns) = bmlt_ground(ew,ns) + internal_melt_rate  ! m/s

                   endif   ! temp > pmptemp
                enddo   ! up

             endif   ! ice is present

          enddo   ! ew
       enddo   ! ns

    endif   ! whichtemp

    ! Cap basal temperature at pressure melting point

    do ns = 1, nsn
       do ew = 1, ewn

          if (ice_mask(ew,ns) == 1 .and. floating_mask(ew,ns) == 0) then  ! ice is present and grounded

             temp(upn,ew,ns) = min (temp(upn,ew,ns), bpmp(ew,ns))

             ! If freeze-on was computed above (bmlt_ground < 0) and Tbed = pmptemp but no basal water is present, then set T(upn) < pmptemp.
             ! Specifically, set Tbed to the temperature of the layer nearest the bed.
             ! In most cases, Tbed will then be computed instead of prescribed during the next time step. 
             ! Note: Energy conservation is not violated here, because no energy is associated with
             !       the infinitesimally thin layer at the bed.

             if (bmlt_ground(ew,ns) < 0.d0 .and. bwat(ew,ns)==0.d0 .and. temp(upn,ew,ns) >= bpmp(ew,ns)) then
                temp(upn,ew,ns) = temp(upn-1,ew,ns)
                bmlt_ground(ew,ns) = 0.d0   ! Set freeze-on to zero since no water is present
             endif

          endif   ! ice is present and grounded

       enddo   ! ew
    enddo   ! ns

  end subroutine glissade_basal_melting_ground

!=======================================================================

  subroutine glissade_pressure_melting_point(depth, pmptemp)

    ! Compute the pressure melting point temperature at a given depth

    use glimmer_physcon, only : rhoi, grav, pmlt 

    real(dp), intent(in) :: depth      ! depth in column (model thickness units)
    real(dp), intent(out) :: pmptemp   ! pressure melting point temp (deg C)

    pmptemp = - rhoi * grav * pmlt * depth

  end subroutine glissade_pressure_melting_point

!=======================================================================

  subroutine glissade_pressure_melting_point_column(thck, stagsigma, pmptemp)

    ! Compute the pressure melting point temperature in a column
    ! Note: pmptemp and stagsigma should have the same dimension

    use glimmer_physcon, only : rhoi, grav, pmlt 

    real(dp), intent(in) :: thck                    ! ice thickness (m)
    real(dp), dimension(:), intent(in) :: stagsigma ! staggered vertical sigma coordinate
                                                    ! (defined at layer midpoints)
    real(dp), dimension(:), intent(out) :: pmptemp  ! pressure melting point temperature (deg C)

    pmptemp(:) = - rhoi * grav * pmlt * thck * stagsigma(:)

  end subroutine glissade_pressure_melting_point_column

!=======================================================================

  subroutine glissade_enth2temp (stagsigma,           &
                                 thck,    enthalpy,   &
                                 temp,    waterfrac)

    ! Convert from specific enthalpy to ice temperature and water content
	
    use glimmer_physcon, only : rhoi, shci, lhci, rhow

    ! I/O variables
    real(dp), dimension(:), intent(in)                      :: stagsigma ! (1:upn-1)
    real(dp), intent(in)                                    :: thck
    real(dp), dimension(0:size(stagsigma)+1), intent(inout) :: enthalpy  ! (0:upn)
    real(dp), dimension(0:size(stagsigma)+1), intent(out)   :: temp      ! (0:upn)
    real(dp), dimension(size(stagsigma)), intent(out)       :: waterfrac ! (1:upn-1)

    ! local variables
    real(dp), dimension(size(stagsigma))                 :: pmptemp     ! (1:upn-1)
    real(dp)                                             :: pmptemp_bed
    real(dp), dimension(0:size(stagsigma)+1)             :: pmpenthalpy ! (0:upn)
    integer :: up, upn
	
    upn = size(stagsigma) + 1
	
    ! find pmpenthalpy(0:upn)
    call glissade_pressure_melting_point_column(thck, stagsigma(1:upn-1), pmptemp(1:upn-1))
    call glissade_pressure_melting_point(thck, pmptemp_bed)
    pmpenthalpy(0) = 0.d0
    pmpenthalpy(1:upn-1) = pmptemp(1:upn-1) * rhoi*shci
    pmpenthalpy(upn) = pmptemp_bed * rhoi*shci

    ! solve for temp and waterfrac
    if (enthalpy(0) >= pmpenthalpy(0)) then ! temperate ice
       temp(0) = 0.d0                      ! temperate ice
       ! Reset enthalpy to be consistent with the surface temperature.
       ! This is consistent with energy conservation because the top surface
       !  is infinitesimally thin.
       enthalpy(0) = pmpenthalpy(0)
    else
       temp(0) = enthalpy(0) / (rhoi*shci) ! cold ice
    endif
	
    do up = 1, upn-1
       if (enthalpy(up) >= pmpenthalpy(up)) then ! temperate ice
          temp(up) = pmptemp(up)
          waterfrac(up) = (enthalpy(up)-pmpenthalpy(up)) /                 &
                          ((rhow-rhoi) * shci * pmptemp(up) + rhow * lhci)
       else ! cold ice
          temp(up) = enthalpy(up) / (rhoi*shci)
          waterfrac(up) = 0.0d0
       endif
    end do
	
    if (enthalpy(upn) >= pmpenthalpy(upn)) then  ! temperate ice
       temp(upn) = pmptemp_bed
       ! Reset enthalpy to be consistent with the bed temperature.
       ! This is consistent with energy conservation because the basal surface
       !  is infinitesimally thin.
       enthalpy(upn) = pmpenthalpy(upn)
    else
       temp(upn) = enthalpy(upn) / (rhoi*shci)   ! cold ice
    endif

  end subroutine glissade_enth2temp

!=======================================================================

  subroutine glissade_temp2enth (stagsigma,          &
                                 temp, waterfrac,    &
                                 thck, enthalpy)

    ! Convert from temperature and water fraction to specific enthalpy

    use glimmer_physcon, only : rhoi, shci, lhci, rhow

    ! I/O variables
    real(dp), dimension(:), intent(in)                    :: stagsigma ! (1:upn-1)
    real(dp), dimension(0:size(stagsigma)+1), intent(in)  :: temp      ! (0:upn)
    real(dp), dimension(1:size(stagsigma)), intent(in)    :: waterfrac ! (1:upn-1)
    real(dp), intent(in)                                  :: thck
    real(dp), dimension(0:size(stagsigma)+1), intent(out) :: enthalpy  ! (0:upn)

    ! local variables
    real(dp), dimension(size(stagsigma))        :: pmptemp  !(1:upn-1)
    integer :: up, upn
		
    upn = size(stagsigma) + 1
	
    ! find pmptemp in column
    call glissade_pressure_melting_point_column (thck, stagsigma(1:upn-1), pmptemp(1:upn-1))
    
    ! solve for enthalpy
    ! assume waterfrac = 0 at upper and lower ice surfaces
    enthalpy(0) = temp(0) * rhoi * shci
    do up = 1, upn-1
       enthalpy(up) = ((1.d0 - waterfrac(up)) * rhoi * shci * temp(up))          &
                      + waterfrac(up) * rhow * ((shci * pmptemp(up)) + lhci)
    end do
    enthalpy(upn) = temp(upn) * rhoi * shci
	
  end subroutine glissade_temp2enth

! The remaining subroutines are called from glissade.F90
!=======================================================================
 
  subroutine glissade_interior_dissipation_sia(ewn,       nsn,       &
                                               upn,       stagsigma, &
                                               ice_mask,             &
                                               stagthck,  flwa,      &
                                               dusrfdew,  dusrfdns,  &
                                               dissip)

    ! Compute the dissipation source term associated with strain heating,
    ! based on the shallow-ice approximation.
    
    use glimmer_physcon, only : gn   ! Glen's n

    integer, intent(in) :: ewn, nsn, upn   ! grid dimensions

    real(dp), dimension(upn-1), intent(in) :: stagsigma   ! staggered vertical grid for temperature

    real(dp), dimension(:,:), intent(in) :: stagthck, dusrfdew, dusrfdns

    integer, dimension(:,:), intent(in) :: ice_mask    ! = 1 where ice is present (thck > thklim), else = 0

    real(dp), dimension(:,:,:), intent(in) ::  &
         flwa         ! flow factor, Pa^(-n) yr^(-1)  

    real(dp), dimension(:,:,:), intent(out) ::  &
         dissip       ! interior heat dissipation (deg/s)
    
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
       call write_log('Error, glissade SIA dissipation: dissip has the wrong vertical dimension', GM_FATAL)
    endif

    dissip(:,:,:) = 0.0d0

    ! Note: Factor of 16 is for averaging flwa
    sia_dissip_fact(1:upn-1) = (stagsigma(1:upn-1) * rhoi * grav)**p1 * 2.0d0 / (16.0d0 * rhoi * shci)

    do ns = 2, nsn-1
       do ew = 2, ewn-1
          !Note: ice_mask = 1 where thck > thcklim.  Elsewhere, dissipation is assumed to be zero.
          if (ice_mask(ew,ns) == 1) then
             geom_fact = (0.25d0*sum(stagthck(ew-1:ew,ns-1:ns)) * sqrt((0.25d0*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                                                                     + (0.25d0*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             dissip(:,ew,ns) = geom_fact * sia_dissip_fact(:) *   & 
                              (flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                              2.d0*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                              4.d0*flwa(:,ew,ns))
          end if
       end do
    end do

  end subroutine glissade_interior_dissipation_sia

!=======================================================================

  subroutine glissade_interior_dissipation_first_order(ewn,       nsn,       &
                                                       upn,                  &
                                                       ice_mask,             &
                                                       tau_eff,   efvs,      &
                                                       dissip)

    ! Compute the dissipation source term associated with strain heating.
    ! Note that the dissipation is computed in the same way on either a staggered or an
    !  unstaggered vertical grid.  
    ! Note also that dissip and flwa must have the same vertical dimension 
    !  (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
    
    integer, intent(in) :: ewn, nsn, upn   ! grid dimensions
    integer, dimension(:,:), intent(in) :: ice_mask    ! = 1 where ice is present (thck > thklim), else = 0

    real(dp), dimension(:,:,:), intent(in) ::  &
         tau_eff,    & ! effective stress, Pa
         efvs          ! effective viscosity, Pa s

    real(dp), dimension(:,:,:), intent(out) ::  &
         dissip       ! interior heat dissipation (deg/s)
    
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
       call write_log('Error, glissade 1st order dissipation: dissip has the wrong vertical dimension',GM_FATAL)
    endif

    dissip(:,:,:) = 0.0d0
    ho_dissip_fact = 1.d0 / (rhoi*shci)

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

  !TODO - For damage-based calving, try multiplying flwa by a damage factor, (1 - damage)

  subroutine glissade_flow_factor(whichflwa,               whichtemp,  &
                                  stagsigma,                           &
                                  thck,                                &
                                  temp,                                &
                                  flwa,                                &
                                  default_flwa,                        &
                                  flow_enhancement_factor,             &
                                  flow_enhancement_factor_float,       &
                                  floating_mask,                       &
                                  waterfrac)

    ! Calculate Glen's $A$ over the 3D domain, using one of three possible methods.
    !
    ! The primary method is to use this equation from \emph{Paterson and Budd} [1982]:
    ! \[
    ! A(T^{*})=a \exp \left(\frac{-Q}{RT^{*}}\right)
    ! \]
    ! This is equation 9 in {\em Payne and Dongelmans}. $a$ is a constant of proportionality,
    ! $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
    ! The pressure-corrected temperature, $T^{*}$ is given by:
    ! \[
    ! T^{*} = T - T_{\mathrm{pmp}} + T_0
    ! \] 
    ! \[
    ! T_{\mathrm{pmp}} = T_0- \sigma \rho g H \Phi
    ! \]
    ! $T$ is the ice temperature, $T_{\mathrm{pmp}}$ is the pressure melting point 
    ! temperature, $T_0$ is the triple point of water, $\rho$ is the ice density, and 
    ! $\Phi$ is the (constant) rate of change of melting point temperature with pressure.

    use glimmer_physcon, only: scyr, arrmlh, arrmll, actenh, actenl, gascon, trpt
    use glimmer_paramets, only: vis0

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

!   Note: The flwa, temp, and stagsigma arrays should have vertical dimension 1:upn-1.
!         The temperatures at the upper surface (k=1) and bed (k=upn) are not included in the input array.
!   Note: A flow factor is computed for all cells, including cells of zero thickness.
!         Provided zero-thickness cells have a sensible default temperature (e.g., artm),
!          the resulting flwa should be physically reasonable if used in any calculations.

    integer,                    intent(in)    :: whichflwa !> which method of calculating A
    integer,                    intent(in)    :: whichtemp !> which method of calculating temperature;
                                                           !> include waterfrac in calculation if using enthalpy method
    real(dp),dimension(:),      intent(in)    :: stagsigma !> vertical coordinate at layer midpoints
    real(dp),dimension(:,:),    intent(in)    :: thck      !> ice thickness (m)
    real(dp),dimension(:,:,:),  intent(in)    :: temp      !> 3D temperature field (deg C)
    real(dp),dimension(:,:,:),  intent(inout) :: flwa      !> output $A$, in units of Pa^{-n} s^{-1}, allow input for data option
    real(dp), intent(in)                      :: default_flwa  !> Glen's A to use in isothermal case, Pa^{-n} s^{-1} 
    real(dp), intent(in), optional            :: flow_enhancement_factor       !> flow enhancement factor in Arrhenius relationship
    real(dp), intent(in), optional            :: flow_enhancement_factor_float !> flow enhancement factor for floating ice
    integer, dimension(:,:),   intent(in), optional :: floating_mask !> = 1 where ice is present and floating, else = 0
    real(dp),dimension(:,:,:), intent(in), optional :: waterfrac     !> internal water content fraction, 0 to 1

    !> \begin{description}
    !> \item[0] Set to prescribed constant value.
    !> \item[1] {\em Paterson and Budd} relationship, with temperature set to -5$^{\circ}$C.
    !> \item[2] {\em Paterson and Budd} relationship.
    !> \end{description}

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    integer :: ew, ns, up, ewn, nsn, nlayers
    real(dp), dimension(size(stagsigma)) :: pmptemp   ! pressure melting point temperature

    real(dp), dimension(:,:), allocatable :: &
         enhancement_factor      ! flow enhancement factor in Arrhenius relationship

    real(dp) :: tempcor                 ! temperature relative to pressure melting point

    real(dp),dimension(4), parameter ::  &
       arrfact = (/ arrmlh,             &   ! Value of a when T* is above -263K, Pa^{-n} s^{-1}
                    arrmll,             &   ! Value of a when T* is below -263K, Pa^{-n} s^{-1}
                   -actenh / gascon,    &   ! Value of -Q/R when T* is above -263K
                   -actenl / gascon/)       ! Value of -Q/R when T* is below -263K
    
    real(dp), parameter :: const_temp = -5.0d0   ! deg C
    real(dp), parameter :: flwa_waterfrac_enhance_factor = 181.25d0

    !------------------------------------------------------------------------------------
   
    nlayers = size(flwa,1)   ! upn - 1
    ewn = size(flwa,2)
    nsn = size(flwa,3)

    allocate(enhancement_factor(ewn,nsn))

    if (present(flow_enhancement_factor)) then
       if (present(flow_enhancement_factor_float) .and. present(floating_mask)) then
          do ns = 1, nsn
             do ew = 1, ewn
                if (floating_mask(ew,ns) == 1) then
                   enhancement_factor(ew,ns) = flow_enhancement_factor_float
                else
                   enhancement_factor(ew,ns) = flow_enhancement_factor
                endif
             enddo
          enddo
       else  ! no separate factor for floating ice
          enhancement_factor(:,:) = flow_enhancement_factor
       endif
    else
       enhancement_factor(:,:) = 1.d0
    endif

    ! Check that the temperature array has the desired vertical dimension

    if (size(temp,1) /= size(flwa,1)) then
       call write_log('glissade_flow_factor: temp and flwa must have the same vertical dimensions', GM_FATAL)
    endif

    ! Multiply the default rate factor by the enhancement factor if applicable
    ! Note: Here, default_flwa is assumed to have units of Pa^{-n} s^{-1},
    !       whereas model%paramets%default_flwa has units of Pa^{-n} yr^{-1}.

    ! initialize
    if (whichflwa /= FLWA_INPUT) then
       do ns = 1, nsn
          do ew = 1, ewn
             flwa(:,ew,ns) = enhancement_factor(ew,ns) * default_flwa
          enddo
       enddo
    endif

    select case(whichflwa)

    case(FLWA_PATERSON_BUDD)

      ! This is the Paterson and Budd relationship

      do ns = 1,nsn
         do ew = 1,ewn
            
            call glissade_pressure_melting_point_column (thck(ew,ns), stagsigma, pmptemp)

            do up = 1, nlayers   ! nlayers = upn - 1

               ! Calculate the corrected temperature
               tempcor = min(0.0d0, temp(up,ew,ns) - pmptemp(up))   ! pmptemp < 0
               tempcor = max(-50.0d0, tempcor)

               ! Calculate Glen's A (including flow enhancement factor)

               if (tempcor >= -10.d0) then
                  flwa(up,ew,ns) = enhancement_factor(ew,ns) * arrfact(1) * exp(arrfact(3)/(tempcor + trpt))
               else
                  flwa(up,ew,ns) = enhancement_factor(ew,ns) * arrfact(2) * exp(arrfact(4)/(tempcor + trpt))
               endif

               ! BDM added correction for a liquid water fraction
               ! Using Greve and Blatter (2009) formulation for Glen's A flow rate factor:
               !    A = A(theta_PMP) * (1 + 181.25 * waterfrac)
               if (whichtemp == TEMP_ENTHALPY .and. present(waterfrac)) then
                  if (waterfrac(up,ew,ns) > 0.0d0) then
                     flwa(up,ew,ns) = flwa(up,ew,ns) * (1.d0 + flwa_waterfrac_enhance_factor * waterfrac(up,ew,ns))
                  endif
               endif

            enddo   ! up

         end do     ! ew
      end do        ! ns

    case(FLWA_PATERSON_BUDD_CONST_TEMP)

      ! This is the Paterson and Budd relationship, but with the temperature held to a constant parameter.

      do ns = 1,nsn
         do ew = 1,ewn

            ! Calculate Glen's A with a fixed temperature (including flow enhancement factor)

            if (const_temp >= -10.d0) then
               flwa(:,ew,ns) = enhancement_factor(ew,ns) * arrfact(1) * exp(arrfact(3)/(const_temp + trpt))
            else
               flwa(:,ew,ns) = enhancement_factor(ew,ns) * arrfact(2) * exp(arrfact(4)/(const_temp + trpt))
            endif

         end do
      end do

    case(FLWA_CONST_FLWA)

       ! do nothing (flwa is initialized to default_flwa above)

    case(FLWA_INPUT)
      ! do nothing - use flwa from input or forcing file
      print *, 'FLWA', minval(flwa), maxval(flwa)

    end select

    ! This logic assumes that the input flwa is already in dimensionless model units.
    ! TODO: Make a different assumption about input units?
    if (whichflwa /= FLWA_INPUT) then
       ! Change flwa to model units (glissade_flow_factor assumes SI units of Pa{-n} s^{-1})
       flwa(:,:,:) = flwa(:,:,:) / vis0
    endif

  end subroutine glissade_flow_factor

!=======================================================================

end module glissade_therm

!=======================================================================
