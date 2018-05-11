!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_traction.F90 - part of the Community Ice Sheet Model (CISM)  
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

#include "glide_mask.inc"
#include "config.inc"

  module glissade_basal_traction

  !-----------------------------------------------------------------------------
  ! Compute or prescribe the basal traction coefficient 'beta' as required by
  ! the higher-order velocity solver.
  ! 
  ! Note that beta is assumed to be a positive constant.  In earlier versions of
  ! the code it was called 'betasquared'.
  !
  ! The units are Pa/(m/yr) if we assume a linear sliding law of the form
  !    taub_x = -beta * u, taub_y = -beta * v
  !
  ! However, the units are Pa if beta is treated as a till yield stress.
  !
  ! See glide_types.F90 for the meaning of the various options.
  !
  !-----------------------------------------------------------------------------

  use glimmer_paramets, only : dp
  use glimmer_physcon,  only : scyr
  use glimmer_paramets, only : vel0, tau0
  use glimmer_log
  use glide_types
  use parallel,         only : staggered_parallel_halo  
  use glissade_grid_operators

  implicit none

  private
  public :: calcbeta, calc_effective_pressure

!***********************************************************************

contains

!***********************************************************************

  subroutine calcbeta (whichbabc,                    &
                       dew,           dns,           &
                       ewn,           nsn,           &
                       thisvel,       othervel,      &
                       basal_physics,                &
                       flwa_basal,    thck,          &
                       mask,          beta_external, &
                       beta,                         &
                       topg,          eus,           &
                       ice_mask,                     &
                       floating_mask,                &
                       land_mask,                    &
                       f_ground,                     &
                       itest, jtest,  rtest)

  ! subroutine to calculate map of beta sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").
   
  ! NOTE: Previously, the input arguments were assumed to be dimensionless
  ! and were rescaled in this routine.  Now the input arguments are
  ! assumed to have the units given below.
     
  use glimmer_paramets, only: len0
  use glimmer_physcon, only: gn, pi
  use parallel, only: nhalo, this_rank
  use parallel, only: parallel_globalindex, global_ewn, global_nsn
  use parallel, only: distributed_scatter_var, parallel_halo, main_task

  implicit none

  ! Input/output arguments

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real(dp), intent(in)                    :: dew, dns           ! m
  real(dp), intent(in), dimension(:,:)    :: thisvel, othervel  ! basal velocity components (m/yr)
  type(glide_basal_physics), intent(in)   :: basal_physics      ! basal physics object
  real(dp), intent(in), dimension(:,:)    :: flwa_basal         ! flwa for the basal ice layer (Pa^{-3} yr^{-1})
  real(dp), intent(in), dimension(:,:)    :: thck               ! ice thickness (m)
  integer,  intent(in), dimension(:,:)    :: mask               ! staggered grid mask
  real(dp), intent(in), dimension(:,:)    :: beta_external      ! fixed beta read from external file (Pa yr/m)
  real(dp), intent(inout), dimension(:,:) :: beta               ! basal traction coefficient (Pa yr/m)
                                                                ! Note: This is beta_internal in glissade

  ! Note: Adding fields for parallel ISHOM-C test case
  real(dp), dimension(:,:), allocatable :: beta_global          ! beta on the global grid
  real(dp), dimension(:,:), allocatable :: beta_extend          ! beta extended to the ice grid (dimensions ewn, nsn)

  ! Note: optional arguments topg and eus are used for pseudo-plastic sliding law
  !TODO - Make these argument non-optional? Can do this after removing the call to calcbeta from Glam.
  real(dp), intent(in), dimension(:,:), optional :: topg         ! bed topography (m)
  real(dp), intent(in), optional :: eus                          ! eustatic sea level (m) relative to z = 0

  integer, intent(in), dimension(:,:), optional :: &
       ice_mask,        & ! = 1 where ice is present (thck > thklim), else = 0
       floating_mask,   & ! = 1 where ice is present and floating, else = 0
       land_mask          ! = 1 where topg > eus

  real(dp), intent(in), dimension(:,:), optional :: f_ground     ! grounded ice fraction, 0 <= f_ground <= 1
  integer, intent(in), optional :: itest, jtest, rtest           ! coordinates of diagnostic point

  ! Local variables

  real(dp) :: smallnum = 1.0d-2  ! m/yr

  real(dp) :: Ldomain   ! size of full domain
  real(dp) :: omega     ! frequency of beta field
  real(dp) :: dx, dy
  integer :: ew, ns

  real(dp), dimension(size(beta,1), size(beta,2)) :: speed      ! ice speed, sqrt(uvel^2 + vvel^2), m/yr

  ! variables for power law
  real(dp) :: powerlaw_p, powerlaw_q

  ! variables for Coulomb friction law
  real(dp) :: coulomb_c   ! Coulomb law friction coefficient (unitless)
  real(dp) :: powerlaw_c  ! power law friction coefficient (Pa m^{-1/3} yr^{1/3})
  real(dp) :: lambda_max  ! wavelength of bedrock bumps at subgrid scale (m)
  real(dp) :: m_max       ! maximum bed obstacle slope (unitless)
  real(dp) :: m           ! exponent m in power law
  integer, dimension(size(thck,1), size(thck,2)) :: &
       ice_or_land_mask, &! = 1 where ice_mask = 1 or land_mask = 1, else = 0       
       imask              ! = 1 where thck > 0, else = 1

  real(dp), dimension(size(beta,1), size(beta,2)) ::  &
       big_lambda,                 &  ! bedrock characteristics
       flwa_basal_stag                ! basal flwa interpolated to the staggered grid (Pa^{-n} yr^{-1})

  ! variables for Tsai et al. parameterization
  real(dp) :: taub_powerlaw  ! basal shear stress given by a power law as in Tsai et al. (2015)
  real(dp) :: taub_coulomb   ! basal shear stress given by Coulomb friction as in Tsai et al. (2015)

  ! variables for pseudo-plastic law
  real(dp) :: q              ! exponent for pseudo-plastic law (unitless)
                             ! q = 1 for linear, q = 0 for plastic, 0 < q < 1 for power law
  real(dp) :: u0             ! threshold velocity for pseudo-plastic law (m/yr)
  real(dp) :: phi            ! phi for pseudoplastic law (degress, 0 < phi < 90)
  real(dp) :: tanphi         ! tan(phi) for pseudo-plastic law (unitless)
  real(dp) :: bed            ! bed elevation, topg - eus (m)
  real(dp) :: phimin, phimax ! min and max values of phi for pseudo-plastic law (degrees)
  real(dp) :: bedmin, bedmax ! bed elevations (m) below which phi = phimin and above which phi = phimax
  real(dp) :: tau_c          ! yield stress for pseudo-plastic law (unitless)
  real(dp) :: numerator, denominator

  character(len=300) :: message

  integer :: iglobal, jglobal

  logical, parameter :: verbose_beta = .false.

  ! Compute the ice speed: used in power laws where beta = beta(u).
  ! Enforce a minimum speed to prevent beta from become very large when velocity is small.
  speed(:,:) = dsqrt(thisvel(:,:)**2 + othervel(:,:)**2 + smallnum**2)

  ! If beta_powerlaw_umax is set to a nonzero value, then limit the speed to this value.
  ! Note: The actual ice speed can be greater than umax.  This is just a way of shutting off the feedback
  !        between beta and ice speed (beta down as speed up) when the ice speed is large.
  !       It helps make the model more stable.
  if (basal_physics%beta_powerlaw_umax > 0.0d0) then
     speed(:,:) = min(speed(:,:), basal_physics%beta_powerlaw_umax)
  endif

  ! Compute beta based on whichbabc

  select case(whichbabc)

    case(HO_BABC_BETA_CONSTANT)   ! spatially uniform beta value; useful for debugging and test cases

          beta(:,:) = basal_physics%ho_beta_const  ! Pa yr/m

    case(HO_BABC_BETA_BPMP)  ! large value for frozen bed, lower value for bed at pressure melting point

          where(basal_physics%bpmp_mask == 1)      ! bed is at pressure melting point
             beta(:,:) = basal_physics%ho_beta_small    ! Pa yr/m
          elsewhere 
             beta(:,:) = basal_physics%ho_beta_large    ! Pa yr/m
          endwhere

    case(HO_BABC_PSEUDO_PLASTIC)

       ! Pseudo-plastic sliding law from PISM:
       !
       ! (tau_bx,tau_by) = -tau_c * (u,v) / (u_0^q * |u|^(1-q))
       ! where the yield stress tau_c = tan(phi) * N
       ! N = effective pressure, computed in subroutine calc_effective_pressure
       ! q, u0 and phi are user-configurable parameters:
       !    q = exponent (q = 1 for linear sliding, q = 0 for a plastic bed, 0 < q < 1 for power-law behavior), default = 1/3
       !    u0 = threshold velocity (the velocity at which tau_b = tau_c), default = 100 m/yr
       !    0 < tan(phi) < 1
       ! As in PISM, phi is allowed to vary with bed elevation
       ! See Aschwanden et al. (2013), The Cryosphere, 7, 1083-1093, Supplement; see also the PISM Users Guide.

       phimin = basal_physics%pseudo_plastic_phimin
       phimax = basal_physics%pseudo_plastic_phimax
       bedmin = basal_physics%pseudo_plastic_bedmin
       bedmax = basal_physics%pseudo_plastic_bedmax

       q = basal_physics%pseudo_plastic_q
       u0 = basal_physics%pseudo_plastic_u0

       if (present(topg) .and. present(eus)) then

          do ns = 1, nsn-1
             do ew = 1, ewn-1

                ! compute tan(phi) based on bed elevation
                bed = topg(ew,ns) - eus
                if (bed <= bedmin) then
                   phi = phimin
                elseif (bed >= bedmax) then
                   phi = phimax
                else   ! bed elevation is between bedmin and bedmax
                   phi = phimin + ((bed - bedmin)/(bedmax - bedmin)) * (phimax - phimin)
                endif
                tanphi = tan(phi * pi/180.d0)

                ! compute beta based on tan(phi), N and u
                tau_c = tanphi * basal_physics%effecpress_stag(ew,ns) 
                beta(ew,ns) = tau_c / (u0**q * speed(ew,ns)**(1.0d0 - q))

                !WHL - debug
                if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
                   if (this_rank == rtest .and. ew == itest .and. ns == jtest) then
                      write(6,*) 'i, j, bed, phi, tanphi, tau_c, speed, beta:', &
                           ew, ns, bed, phi, tanphi, tau_c, speed(ew,ns), beta(ew,ns)
                   endif
                endif

             enddo
          enddo

       else

          call write_log('Pseudo-plastic sliding law requires topg and eus as input arguments', GM_FATAL)

       endif

    case(HO_BABC_YIELD_PICARD)  ! take input value for till yield stress and force beta to be implemented such
                                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option could provide the till yield stress as calculate from the basal processes submodel.
      !!!       Currently, to enable sliding over plastic till, simply specify the value of "beta" as 
      !!!       if it were the till yield stress (in units of Pascals).
      
      beta(:,:) = basal_physics%mintauf(:,:)*tau0 &                                      ! plastic yield stress (converted to Pa)
                         / dsqrt( thisvel(:,:)**2 + othervel(:,:)**2 + (smallnum)**2 )   ! velocity components (m/yr)

      !!! since beta is updated here, communicate that info to halos
      call staggered_parallel_halo(beta)

    case(HO_BABC_BETA_LARGE)      ! frozen (u=v=0) ice-bed interface

       !Note: This option is redundant in that it could be implemented using HO_BETA_CONST
       !      But keeping it for historical reasons since many config files use it

       beta(:,:) = basal_physics%ho_beta_large      ! Pa yr/m  (= 1.0d10 by default)

    case(HO_BABC_ISHOMC)          ! prescribe according to ISMIP-HOM test C

       !Note: Ideally, beta would be read in from an external netCDF file.
       !      However, this is not possible given that the global velocity grid is smaller
       !       than the ice grid and hence not able to fit the full beta field.
       !      The following code sets beta on the full grid as prescribed by Pattyn et al. (2008).

       ! Allocate a global array on the main task only.
       ! On other tasks, allocate a size 0 array, since distributed_scatter_var wants to deallocate on all tasks.
       if (main_task) then
          allocate(beta_global(global_ewn, global_nsn))
       else
          allocate(beta_global(0,0))
       endif

       ! Prescribe beta as in Pattyn et al., The Cryosphere, 2008
       ! Note: These beta values live at vertices, not cell centers.
       !       They need a global array of size (ewn,nsn) to hold values on the global boundary.
       if (main_task) then

          Ldomain = global_ewn * dew   ! size of full domain (must be square)
          omega = 2.d0*pi / Ldomain

          beta_global(:,:) = 0.0d0
          do ns = 1, global_nsn
             do ew = 1, global_ewn
                dx = dew * ew
                dy = dns * ns
                beta_global(ew,ns) = 1000.d0 + 1000.d0 * sin(omega*dx) * sin(omega*dy)
             enddo
          enddo

       endif

       ! Scatter the global beta values back to local arrays
       ! Note: beta_extend has dimensions (ewn,nsn), so it can receive scattered data from beta_global.
       allocate(beta_extend(ewn, nsn))
       beta_extend(:,:) = 0.d0
       call distributed_scatter_var(beta_extend, beta_global)

       ! distributed_scatter_var does not update the halo, so do an update here
       call parallel_halo(beta_extend)

       ! Copy beta_extend to beta on the local processor.
       ! This is done since beta lives on the velocity grid and has dimensions (ewn-1,nsn-1).
       beta(:,:) = 0.d0
       do ns = 1, nsn-1
          do ew = 1, ewn-1
             beta(ew,ns) = beta_extend(ew, ns)
          enddo
       enddo

      ! beta_extend is no longer needed (beta_global is deallocated in distributed_scatter_var)
      deallocate(beta_extend)

    case(HO_BABC_BETA_EXTERNAL)   ! use beta value from external file

       ! set beta to the prescribed external value
       ! Note: This assumes that beta_external has units of Pa yr/m on input.
       beta(:,:) = beta_external(:,:)

       ! beta is initialized to a negative value; we can use that fact to check whether
       ! it has been read correctly from the file
       if (maxval(beta) < 0.d0) then
          call write_log('ERROR: Trying to use HO_BABC_EXTERNAL_BETA, but all beta values are < 0,')
          call write_log('which implies that beta could not be read from the input file.')
          call write_log('Make sure that beta is in the cism input file,')
          call write_log('or change which_ho_babc to a different option.')
          call write_log('Invalid value for beta. See log file for details.', GM_FATAL)
       end if

    case(HO_BABC_POWERLAW)   ! a simple power law
       !   Assume taub = C * ub^(1/m)
       ! implying beta = C * ub^(1/m - 1) 
       ! m should be a positive exponent

       ! Set beta assuming a spatially uniform value of powerlaw_c
       beta(:,:) = basal_physics%powerlaw_c * speed(:,:)**(1.0d0/basal_physics%powerlaw_m - 1.0d0)

    case(HO_BABC_POWERLAW_EFFECPRESS)   ! a power law that uses effective pressure
       !TODO - Remove POWERLAW_EFFECPRESS option? Rarely if ever used.
       ! See Cuffey & Paterson, Physics of Glaciers, 4th Ed. (2010), p. 240, eq. 7.17
       ! This is based on Weertman's classic sliding relation (1957) augmented by the bed-separation index described by Bindschadler (1983)
       !   ub = k taub^p N^-q
       ! rearranging for taub gives:
       !   taub = k^(-1/p) ub^(1/p) N^(q/p)

       ! p and q should be _positive_ exponents. If p/=1, this is nonlinear in velocity.
       ! Cuffey & Paterson recommend p=3 and q=1, and k dependent on thermal & mechanical properties of ice and inversely on bed roughness.   
       !TODO - Change powerlaw_p to powerlaw_m, and make powerlaw_q a config parameter

       powerlaw_p = 3.0d0
       powerlaw_q = 1.0d0

       beta(:,:) = basal_physics%friction_powerlaw_k**(-1.0d0/powerlaw_p) &
            * basal_physics%effecpress_stag(:,:)**(powerlaw_q/powerlaw_p) &
            * speed(:,:)**(1.0d0/powerlaw_p - 1.0d0)

    case(HO_BABC_COULOMB_FRICTION)

      ! Basal stress representation using Coulomb friction law
      ! Coulomb sliding law: Schoof 2005 PRS, eqn. 6.2  (see also Pimentel, Flowers & Schoof 2010 JGR)

       ! Set up parameters needed for the friction law
       m_max = basal_physics%coulomb_bump_max_slope       ! maximum bed obstacle slope(unitless)
       lambda_max = basal_physics%coulomb_bump_wavelength ! wavelength of bedrock bumps (m)
       coulomb_c = basal_physics%coulomb_c                ! basal shear stress factor (Pa (m^-1 y)^1/3)

       ! Need flwa of the basal layer on the staggered grid
       !TODO - Pass in ice_mask instead of computing imask here?
       !       (Small difference: ice_mask = 1 where thck > thklim rather than thck > 0)
       where (thck > 0.d0)
          imask = 1
       elsewhere
          imask = 0
       end where
       call glissade_stagger(ewn,         nsn,               &
                             flwa_basal,  flwa_basal_stag,   &
                             imask,       stagger_margin_in = 1)
       ! TODO Not sure if a halo update is needed on flwa_basal_stag!  I don't think so if nhalo>=2.

       ! Compute biglambda = wavelength of bedrock bumps [m] * flwa [Pa^-n yr^-1] / max bed obstacle slope [dimensionless]
       big_lambda(:,:) = (lambda_max / m_max) * flwa_basal_stag(:,:)

       ! Note: For MISMIP3D, coulomb_c is multiplied by a spatial factor (C_space_factor) which is
       !       read in during initialization. This factor is typically between 0 and 1.
       !       If this factor is not present in the input file, it is set to 1 everywhere.

       ! Compute beta
       ! gn = Glen's n from physcon module
       beta(:,:) = coulomb_c * basal_physics%C_space_factor_stag(:,:) * &
            basal_physics%effecpress_stag(:,:) * speed(:,:)**(1.0d0/gn - 1.0d0) * &
            (speed(:,:) + basal_physics%effecpress_stag(:,:)**gn * big_lambda)**(-1.0d0/gn)

       ! Limit for numerical stability
       !TODO - Is limiting needed?
       where (beta > 1.0d8)
          beta = 1.0d8
       end where

    case(HO_BABC_COULOMB_POWERLAW_SCHOOF)

       ! Use a constant value of basal flwa.
       ! This allows several Coulomb parameters (lambda_max, m_max and flwa_basal)
       !  to be combined into a single parameter powerlaw_c, as in the Tsai power law below.
       !
       ! The equation for tau_b = beta * u_b is
       ! 
       !                    powerlaw_c * coulomb_c * N
       ! tau_b = ---------------------------------------------- u_b^{1/m}
       !         [powerlaw_c^m * u_b + (coulomb_c * N)^m]^{1/m}
       !
       ! where m = powerlaw_m
       !
       ! This is the second modified basal traction law in MISMIP+. See Eq. 11 of Asay-Davis et al. (2016).
       ! Note: powerlaw_c corresponds to beta^2 in their notation, and coulomb_c corresponds to alpha^2.

       ! use constant powerlaw_c and coulomb_c
       powerlaw_c = basal_physics%powerlaw_c
       coulomb_c = basal_physics%coulomb_c
       m = basal_physics%powerlaw_m

       do ns = 1, nsn-1
          do ew = 1, ewn-1

             numerator = powerlaw_c * coulomb_c * basal_physics%effecpress_stag(ew,ns)
             denominator = ( powerlaw_c**m * speed(ew,ns) +  &
                            (coulomb_c * basal_physics%effecpress_stag(ew,ns))**m )**(1.d0/m)
             beta(ew,ns) = (numerator/denominator) * speed(ew,ns)**(1.d0/m - 1.d0)
          enddo
       enddo

       ! Limit for numerical stability
       !TODO - Is limiting needed?
       where (beta > 1.0d8)
          beta = 1.0d8
       end where

       !WHL - debug - Write values along a flowline
!      write(6,*) ' '
!      write(6,*) 'Apply Coulomb friction: i, j, speed, N_stag, beta, taub:'
!      ns = jtest
!      do ew = itest, itest+15
!         write(6,*) ew, ns, speed(ew,ns), basal_physics%effecpress_stag(ew,ns), beta(ew,ns), beta(ew,ns)*speed(ew,ns)
!      enddo

    case(HO_BABC_COULOMB_POWERLAW_TSAI)

      ! Basal stress representation based on Tsai et al. (2015)
      ! The basal stress is the minimum of two values:
      ! (1) power law:          tau_b = powerlaw_c * |u_b|^(1/powerlaw_m)
      ! (2) Coulomb friction:   tau_b = coulomb_c * N
      !                             N = effective pressure = rhoi*g*(H - H_f)
      !                           H_f = flotation thickness = (rhow/rhoi)*(eus-topg)
      ! This value of N is obtained by setting basal_water = BWATER_OCEAN_PENETRATION = 4 
      !  with p_ocean_penetration = 1.0 in the config file.
      ! The other parameters (powerlaw_c, powerlaw_m and coulomb_c) can also be set in the config file.

       !WHL - debug - write out basal stresses
!       write(6,*) ' '
!       write(6,*) 'powerlaw_c, powerlaw_m, Coulomb_c =', &
!           basal_physics%powerlaw_c, basal_physics%powerlaw_m, basal_physics%coulomb_c
!       write(6,*) 'Apply Tsai parameterization: i, j, speed, beta, taub, taub_powerlaw, taub_coulomb, effecpress:'

       do ns = 1, nsn-1
          do ew = 1, ewn-1
             
             taub_powerlaw = basal_physics%powerlaw_c * speed(ew,ns)**(1.d0/basal_physics%powerlaw_m)
             taub_coulomb  = basal_physics%coulomb_c * basal_physics%effecpress_stag(ew,ns)

             if (taub_coulomb <= taub_powerlaw) then   ! apply Coulomb stress, which is smaller
                beta(ew,ns) = taub_coulomb / speed(ew,ns)
             else  ! apply power-law stress
                beta(ew,ns) = taub_powerlaw / speed(ew,ns)
             endif

!             !WHL - debug - Write values along a flowline
!             if (ns == jtest .and. ew >= itest .and. ew <= itest+15) then
!                write(6,*) ew, ns, speed(ew,ns), beta(ew,ns), speed(ew,ns)*beta(ew,ns), &
!                     taub_powerlaw, taub_coulomb, basal_physics%effecpress_stag(ew,ns)
!             endif

          enddo   ! ew
       enddo   ! ns

    case(HO_BABC_SIMPLE)    ! simple pattern; also useful for debugging and test cases
                            ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      beta(:,:) = 1.d4        ! Pa yr/m

      !TODO - Change this loop to work in parallel (set beta on the global grid and scatter to local)
      do ns = 5, nsn-5
         do ew = 1, ewn-1
            beta(ew,ns) = 100.d0      ! Pa yr/m
         end do
      end do

    case default
       ! do nothing

   end select

   ! If f_ground is passed in (as for Glissade), then multiply beta by f_ground (0 <= f_ground <= 1) 
   !  to reduce the basal traction in regions that are partially or totally floating.
   ! Note: With a GLP, f_ground will have values between 0 and 1 at vertices adjacent to the GL.
   !       Without a GLP, f_ground = 0 or 1 everywhere based on a flotation criterion.
   !       By convention, f_ground = 0 where no ice is present.
   !
   ! If f_ground in not passed in (as for Glam), then check for areas where the ice is floating
   !  and make sure beta in these regions is 0. 
   !TODO - Replace GLIDE_IS_FLOAT with floating_mask

   if (present(f_ground)) then   ! Multiply beta by grounded ice fraction

      beta(:,:) = beta(:,:) * f_ground(:,:)

   else    ! set beta = 0 where Glide mask says the ice is floating

      do ns = 1, nsn-1
         do ew = 1, ewn-1 
            if (GLIDE_IS_FLOAT(mask(ew,ns))) then
               beta(ew,ns) = 0.d0
            endif
         end do
      end do

   endif   ! present(f_ground)

   ! For beta close to 0 beneath grounded ice, it is possible to generate unrealistically fast flow.
   ! This could happen, for example, when reading beta from an external file.
   ! To prevent this, set beta to a minimum value beneath grounded ice.
   ! The default value of beta_grounded_min = 0.0, in which case this loop has no effect.
   ! However, beta_grounded_min can be set to a nonzero value in the config file.
  
   if (present(f_ground)) then
      do ns = 1, nsn-1
         do ew = 1, ewn-1
            if (f_ground(ew,ns) > 0.d0 .and. beta(ew,ns) < basal_physics%beta_grounded_min) then
               beta(ew,ns) = basal_physics%beta_grounded_min
            endif
         enddo
      enddo
   endif   ! present(f_ground)

   ! Bug check: Make sure beta >= 0
   ! This check will find negative values as well as NaNs
   do ns = 1, nsn-1
      do ew = 1, ewn-1 
         if (beta(ew,ns) >= 0.d0) then
            ! do nothing
         else
            call parallel_globalindex(ew, ns, iglobal, jglobal)
            write(message,*) 'Invalid beta value in calcbeta: this_rank, i, j, iglobal, jglobal, beta, f_ground:', &
                 this_rank, ew, ns, iglobal, jglobal, beta(ew,ns), f_ground(ew,ns)
            call write_log(trim(message), GM_FATAL)
         endif
      end do
   end do
   
   ! halo update
   call staggered_parallel_halo(beta)

                !WHL - debug
                if (verbose_beta .and. present(rtest) .and. present(itest) .and. present(jtest)) then
                   if (this_rank == rtest) then
                      ew = itest; ns = jtest
!!                      ew = itest-1; ns = jtest
                      write(6,*) 'End of calcbeta, r, i, j, speed, beta:', &
                           rtest, ew, ns, speed(ew,ns), beta(ew,ns)
                   endif
                endif

  end subroutine calcbeta

!***********************************************************************

  subroutine calc_effective_pressure (which_effecpress,             &
                                      ewn,           nsn,           &
                                      basal_physics,                &
                                      ice_mask,      floating_mask, &
                                      thck,          topg,          &
                                      eus,                          &
                                      delta_bpmp,                   &
                                      bmlt,          bwat)

    ! Calculate the effective pressure at the bed

    use glimmer_physcon, only: rhoi, grav, rhoo

    use parallel

    implicit none

    ! Input/output arguments

    integer, intent(in) :: &
         which_effecpress    ! input option for effective pressure

    integer, intent(in) :: &
         ewn, nsn            ! grid dimensions

    type(glide_basal_physics), intent(inout) :: &
         basal_physics       ! basal physics object
                             ! includes effecpress, effecpress_stag and various parameters

    integer, dimension(:,:), intent(in) :: &
         ice_mask,         & ! = 1 where ice is present (thk > thklim), else = 0
         floating_mask       ! = 1 where ice is present and floating, else = 0
 
    real(dp), dimension(:,:), intent(in) ::  &
         thck,             & ! ice thickness (m)
         topg                ! bed topography (m)

    real(dp), intent(in) ::  &
         eus                 ! eustatic sea level (m) relative to z = 0

    !NOTE: If used, the following 2D fields (delta_bpmp, bmlt, bwat, thck and topg) need to be correct in halos.

    real(dp), dimension(:,:), intent(in), optional ::  &
         delta_bpmp          ! Tpmp - T at the bed (deg C)
                             ! used for HO_EFFECPRESS_BPMP option

    real(dp), dimension(:,:), intent(in), optional ::  &
         bmlt                ! basal melt rate at the bed (m/yr)
                             ! used for HO_EFFECPRESS_BMLT option

    real(dp), dimension(:,:), intent(in), optional ::  &
         bwat                ! basal water thickness at the bed (m)
                             ! used for HO_EFFECPRESS_BWAT option

    ! Local variables

    real(dp) :: &
         bpmp_factor,     &  ! factor between 0 and 1, used in linear ramp based on bpmp
         bmlt_factor,     &  ! factor between 0 and 1, used in linear ramp based on bmlt
         relative_bwat       ! ratio bwat/bwat_till_max, limited to range [0,1]

    real(dp), dimension(ewn,nsn) ::  &
         overburden          ! overburden pressure, rhoi*g*H

    real(dp) :: ocean_p           ! exponent in effective pressure parameterization, 0 <= ocean_p <= 1
    real(dp) :: f_pattyn          ! rhoo*(eus-topg)/(rhoi*thck)
                                  ! = 1 at grounding line, < 1 for grounded ice, > 1 for floating ice
    real(dp) :: f_pattyn_capped   ! f_pattyn capped to lie in range [0,1]

    integer :: i, j

    ! Initialize the effective pressure N to the overburden pressure, rhoi*g*H

    overburden(:,:) = rhoi*grav*thck(:,:)

    select case(which_effecpress)

    case(HO_EFFECPRESS_OVERBURDEN)

       basal_physics%effecpress(:,:) = overburden(:,:)

       ! Note: Here we assume (unrealistically) that N = rhoi*g*H even for floating ice.

    case(HO_EFFECPRESS_BPMP)

       if (present(delta_bpmp)) then

          ! Reduce N where the basal temperature is near the pressure melting point,
          !  as defined by delta_bpmp = bpmp - Tbed.
          ! bpmp_factor = 0 where the bed is thawed (delta_bpmp <= 0)
          ! bpmp_factor = 1 where the bed is frozen (delta_bpmp >= effecpress_bpmp_threshold)
          ! 0 < bpmp_factor < 1 where 0 < delta_bpmp < bpmp_threshold 

          do j = 1, nsn
             do i = 1, ewn

                bpmp_factor = max(0.0d0, min(1.0d0, delta_bpmp(i,j)/basal_physics%effecpress_bpmp_threshold))
                basal_physics%effecpress(i,j) = overburden(i,j) * &
                     (basal_physics%effecpress_delta + bpmp_factor * (1.0d0 - basal_physics%effecpress_delta))

                ! set to zero for floating ice
                if (floating_mask(i,j) == 1) basal_physics%effecpress(i,j) = 0.0d0

             enddo
          enddo

       endif   ! present(delta_bpmp)

    case(HO_EFFECPRESS_BMLT)

       if (present(bmlt)) then

          ! Reduce N where there is melting at the bed.
          ! The effective pressure ramps down from full overburden for bmlt = 0
          !  to a small value for bmlt >= effecpress_bmlt_threshold.
          ! Both bmlt and effecpress_bmlt_threshold have units of m/yr.
          ! bmlt_factor = 0 where there is no basal melting (bmlt = 0)
          ! bmlt_factor = 1 where there is large basal melting (bmlt >= effecpress_bmlt_threshold)
          ! 0 < bmlt_factor < 1 where 0 < bmlt < bmlt_threshold 

          do j = 1, nsn
             do i = 1, ewn

                bmlt_factor = max(0.0d0, min(1.0d0, bmlt(i,j)/basal_physics%effecpress_bmlt_threshold))
                basal_physics%effecpress(i,j) = overburden(i,j) * &
                     (basal_physics%effecpress_delta + (1.0d0 - bmlt_factor) * (1.0d0 - basal_physics%effecpress_delta))

                ! set to zero for floating ice
                if (floating_mask(i,j) == 1) basal_physics%effecpress(i,j) = 0.0d0

             enddo
          enddo

       endif   ! present(bmlt)

    case(HO_EFFECPRESS_BWAT)

       ! Initialize for the case where bwat isn't present, and also for points with bwat == 0

       basal_physics%effecpress(:,:) = overburden(:,:)

       if (present(bwat)) then

          ! Reduce N where basal water is present.
          ! The effective pressure decreases from overburden P_0 for bwat = 0 to a small value for bwat = bwat_till_max.
          ! Note: Instead of using a linear ramp for the variation between overburden and the small value
          !       (as for the BPMP and BMLT options above), we use the published formulation of Bueler & van Pelt (2015).
          !       This formulation has N = P_0 for bwat up to ~0.6*bwat_till_max; then N decreases as bwat => bwat_till_max.
          !       See Fig. 1b of Bueler & van Pelt (2015).

          do j = 1, nsn
             do i = 1, ewn

                if (bwat(i,j) > 0.0d0) then

                   relative_bwat = max(0.0d0, min(bwat(i,j)/basal_physics%bwat_till_max, 1.0d0))

                   ! Eq. 23 from Bueler & van Pelt (2015)
                   basal_physics%effecpress(i,j) = basal_physics%N_0  &
                        * (basal_physics%effecpress_delta * overburden(i,j) / basal_physics%N_0)**relative_bwat  &
                        * 10.d0**((basal_physics%e_0/basal_physics%C_c) * (1.0d0 - relative_bwat))

                   ! The following line (if uncommented) would implement Eq. 5 of Aschwanden et al. (2016).
                   ! Results are similar to Bueler & van Pelt, but the dropoff in N from P_0 to delta*P_0 begins
                   !  with a larger value of bwat (~0.7*bwat_till_max instead of 0.6*bwat_till_max).

!!                 basal_physics%effecpress(i,j) = basal_physics%effecpress_delta * overburden(i,j)  &
!!                      * 10.d0**((basal_physics%e_0/basal_physics%C_c) * (1.0d0 - relative_bwat))

                   !WHL - Uncomment to try a linear ramp in place of the Bueler & van Pelt relationship.
                   !      This might lead to smoother variations in N with spatial variation in bwat.
!!                 basal_physics%effecpress(i,j) = overburden(i,j) * &
!!                      (basal_physics%effecpress_delta + (1.0d0 - relative_bwat) * (1.0d0 - basal_physics%effecpress_delta))

                   ! limit so as not to exceed overburden
                   basal_physics%effecpress(i,j) = min(basal_physics%effecpress(i,j), overburden(i,j))
                end if
             enddo
          enddo

       endif   ! present(bwat)

       where (floating_mask == 1)
          ! set to zero for floating ice
          basal_physics%effecpress = 0.0d0
       end where

    case(HO_EFFECPRESS_OCEAN_PENETRATION)

       ! Reduce N for ice grounded below sea level based on connectivity of subglacial water to the ocean
       ! p = 1 => full connectivity
       ! 0 < p < 1 => partial connectivity
       ! p = 0 => no connectivity; p_w = 0

       ocean_p = basal_physics%p_ocean_penetration

       if (ocean_p > 0.0d0) then

          do j = 1, nsn
             do i = 1, ewn
                if (thck(i,j) > 0.0d0) then
                   f_pattyn = rhoo*(eus-topg(i,j)) / (rhoi*thck(i,j))    ! > 1 for floating, < 1 for grounded
                   f_pattyn_capped = max( min(f_pattyn,1.0d0), 0.0d0)    ! capped to lie in the range [0,1]
                   basal_physics%effecpress(i,j) = overburden(i,j) * (1.0d0 - f_pattyn_capped)**ocean_p
                else
                   basal_physics%effecpress(i,j) = 0.0d0
                endif
             enddo
          enddo

       else

          basal_physics%effecpress(:,:) = overburden(:,:)

       endif

       ! Note(WHL): If ocean_p = 0, then we have N = rhoi*grav*H for floating ice (f_pattyn_capped = 1).
       !            Equivalently, we are defining 0^0 = 1 for purposes of the Leguy et al. effective pressure parameterization.
       !            This is OK for the Schoof basal friction law provided that the resulting beta is multiplied by f_ground,
       !             where f_ground is the fraction of floating ice at a vertex, with f_ground = 0 if all four neighbor cells are floating.
       !            If we were to set N = 0 where f_pattyn_capped = 1 (i.e., defining 0^0 = 0), then we would have a
       !             sudden sharp increase in N_stag (the effective pressure at the vertex) when f_pattyn_capped at a cell center
       !             falls from 1 to a value slightly below 1.  This sudden increase would occur despite the use of a GLP.

    end select

    ! Cap the effective pressure at 0x and 1x overburden pressure to avoid strange values going to the friction laws.
    ! This capping may not be necessary, but is included as a precaution.

    where (basal_physics%effecpress < 0.0d0)
       basal_physics%effecpress = 0.0d0
    elsewhere (basal_physics%effecpress > overburden)
       basal_physics%effecpress = overburden
    endwhere

    ! Interpolate the effective pressure to the staggered grid.
    ! stagger_margin_in = 0: Interpolate using values in all cells, including ice-free cells
    ! (to give a smooth transition in N_stag as a cell switches from ice-free to ice-covered)
    !TODO - Does ice_mask need to be passed in? Modify glissade_stagger so it can be called without a mask.

    call glissade_stagger(ewn,                       nsn,                             &
                          basal_physics%effecpress,  basal_physics%effecpress_stag,   &
                          ice_mask,                  stagger_margin_in = 0)

 end subroutine calc_effective_pressure

!=======================================================================

end module glissade_basal_traction

!=======================================================================

