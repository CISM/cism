!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade.F90 - part of the Community Ice Sheet Model (CISM)  
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

! WJS (1-30-12): The following (turning optimization off) is needed as a workaround for an
! xlf compiler bug, at least in IBM XL Fortran for AIX, V12.1 on bluefire
#ifdef CPRIBM
@PROCESS OPT(0)
#endif

!CLEANUP - glissade.F90
!
! NOTE: MJH Lines that start with !### are ones I have identified to be deleted.
!
! This is a new module, originally copied from glide.F90 (William Lipscomb, June 2012)
! Removed SIA-specific code, leaving only the HO code with remapping transport
! Whenever possible, parallel_halo updates should go in this module rather
!  than at lower levels.
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glissade.f90 - part of the CISM ice model        + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade

  ! Driver for Glissade (parallel, higher-order) dynamical core

  use glimmer_global, only: dp
  use glimmer_log
  use glide_types
  use glide_io
  use glide_lithot
  use glimmer_config

  implicit none

  integer :: &
     ntracer = 1   ! number of tracers to transport (just temperature for now)
                   ! BDM change ntracer from parameter so it's not a constant
                   !TODO - Declare elsewhere?

  integer, private, parameter :: dummyunit=99

  !WHL - debug
  logical, parameter :: test_transport = .false.   ! if true, call test_transport subroutine
  real(dp), parameter :: thk_init = 500.d0         ! initial thickness (m) for test_transport

contains

!=======================================================================

! Note: There is no glissade_config subroutine; glide_config works for all dycores.

!=======================================================================

  subroutine glissade_initialise(model)

    ! initialise Glissade model instance

!TODO - Are all of these needed?
    use parallel
    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio
    use glide_velo, only: init_velo  !TODO - Remove this
    use glissade_temp, only: glissade_init_temp
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground
    use glide_thck, only : glide_calclsrf
    use glam_strs2, only : glam_velo_init
    use glimmer_coordinates, only: coordsystem_new
    use glide_grid_operators, only: stagvarb
    use glissade_velo_higher, only: glissade_velo_higher_init
    use glide_diagnostics, only: glide_init_diag
    use felix_dycore_interface, only: felix_velo_init
    use glide_bwater

!WHL - debug
    use glimmer_paramets, only: tau0, vel0, thk0
    use glimmer_physcon, only: scyr

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    !TODO - build glimmer_vers file or put this character elsewhere?
    character(len=100), external :: glimmer_version_char

!WHL - debug
    logical, parameter :: test_halo = .false.   ! if true, call test_halo subroutine

    integer :: i, j, nx, ny

!WHL - for artificial adjustment to ismip-hom surface elevation
    logical, parameter :: ismip_hom_adjust_usrf = .false.
    real(dp) :: usrf_ref

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters (some conversions to SI units)
    call glide_scale_params(model)

    ! set up coordinate systems
    ! time to change to the parallel values of ewn and nsn

    !WHL - added choice between periodic and open global boundary conditions
    if (model%general%global_bc == GLOBAL_BC_OUTFLOW) then
       call distributed_grid(model%general%ewn,model%general%nsn,outflow_bc_in=.true.)
    else
       call distributed_grid(model%general%ewn,model%general%nsn)
    endif

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0,               &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    !TODO - Change 2. to 2.d0
    model%general%velo_grid = coordsystem_new(model%numerics%dew/2., model%numerics%dns/2., &
                                              model%numerics%dew,    model%numerics%dns,    &
                                              model%general%ewn-1,   model%general%nsn-1)


    ! allocate arrays
    call glide_allocarr(model)

    !TODO - May be able to eliminate the bed softness parameter 
    !       and set btrc to model%velowo%btrac_const in glide_velo

    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

    ! set uniform basal heat flux (positive down)
    model%temper%bheatflx = model%paramets%geot

    ! compute sigma levels or load from external file
    ! (if not already read from config file)
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)

    ! and read first time slice
    call glide_io_readall(model,model)

    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)

    select case(model%options%whichrelaxed)
    case(RELAXED_TOPO_INPUT)   ! Supplied topography is relaxed
       model%isostasy%relx = model%geometry%topg
    case(RELAXED_TOPO_COMPUTE) ! Supplied topography is in equilibrium
                               !TODO - Test this case
       call not_parallel(__FILE__,__LINE__)
       call isos_relaxed(model)
    end select

    ! open all output files
    call openall_out(model)

    ! create glide variables
    call glide_io_createall(model, model)

    ! If a 2D bheatflx field is present in the input file, it will have been written 
    !  to model%temper%bheatflx.  For the case model%options%gthf = 0, we want to use
    !  a uniform heat flux instead.
    ! If no bheatflx field is present in the input file, then we default to the 
    !  prescribed uniform value, model%paramets%geot.

    if (model%options%gthf == GTHF_UNIFORM) then

       ! Check to see if this flux was present in the input file
       ! (by checking whether the flux is nonuniform over the domain)
       if (abs(maxval(model%temper%bheatflx) - minval(model%temper%bheatflx)) > 1.d-6) then  
          call write_log('Setting uniform prescribed geothermal flux')
          call write_log('(Set gthf = 1 to read geothermal flux field from input file)')
       endif

       ! set uniform basal heat flux (positive down)
       model%temper%bheatflx = model%paramets%geot

    endif

    ! initialise glissade components

    !TODO - Most of what's done in init_velo is needed for SIA only
    !       Can remove this call provided velowk is not used elsewhere (e.g., to call wvelintg)
    call init_velo(model)

    !BDM - Just call glissade_init_temp, which will call glissade_init_enthalpy
    !      if necessary
    call glissade_init_temp(model) 

    ! Initialize basal hydrology model, if enabled
    call bwater_init(model)

    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call init_lithot(model)
    end if


    ! Dycore-specific velocity solver initialization
    select case (model%options%whichdycore)
    case ( DYCORE_GLAM )   ! glam finite-difference

       call glam_velo_init(model%general%ewn,    model%general%nsn,  &
                           model%general%upn,                        &
                           model%numerics%dew,   model%numerics%dns, &
                           model%numerics%sigma)

    case ( DYCORE_GLISSADE )   ! glissade finite-element

       call glissade_velo_higher_init

    case ( DYCORE_ALBANYFELIX)

       call felix_velo_init(model)

    end select

    ! If unstagbeta (i.e., beta on the scalar ice grid) was read from an input file,
    ! then interpolate it to beta on the staggered grid.
    ! NOTE: unstagbeta is initialized to -999.d0, so its maxval will be > 0 only if
    !       the field is read in.
    ! We make an exception for ISHOM case C, in which case for greater accuracy we want to
    !  set beta in subroutine calcbeta instead of interpolating from unstagbeta.

    if (maxval(model%velocity%unstagbeta) > 0.d0 .and.   &
               model%options%which_ho_babc /= HO_BABC_ISHOMC) then  ! interpolate to staggered grid
       call write_log('Interpolating beta from unstaggered (unstagbeta) to staggered grid (beta)')
       if (maxval(model%velocity%beta) > 0.0d0 ) then
         call write_log('Warning: the input "beta" field will be overwritten with values interpolated from the input "unstagbeta" field!')
       endif

       call parallel_halo(model%velocity%unstagbeta)    ! fill in halo values
       call stagvarb(model%velocity%unstagbeta,  &      ! interpolate
                     model%velocity%beta,        &
                     model%general%ewn,          &
                     model%general%nsn)

       !WHL - debug
       print*, ' '
       print*, 'Interpolating beta from unstaggered to staggered grid'
       print*, ' '
!       print*, 'unstagbeta, Pa/(m/yr):'
!       do j = model%general%nsn, 1, -1
!          do i = 1, model%general%ewn
!             write(6,'(f8.0)',advance='no') model%velocity%unstagbeta(i,j) * tau0/vel0/scyr
!          enddo
!          print*, ' '
!       enddo

!       print*, ' '
!       print*, 'beta, Pa/(m/yr):'
!       do j = model%general%nsn-1, 1, -1
!          do i = 1, model%general%ewn-1
!             write(6,'(f8.0)',advance='no') model%velocity%beta(i,j) * tau0/vel0/scyr
!          enddo
!          print*, ' '
!       enddo

    endif  ! unstagbeta > 0

!WHL - This option is disabled for now.
    ! *mb* added; initialization of basal proc. module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%ntem)
!!    end if      

    ! calculate mask
       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)

    !TODO- Not sure why this needs to be called here.
    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    ! and calculate lower and upper ice surface

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

!WHL - The following is a hack to ensure that usrf is uniform (to double precision) for 
!      a given value of i for the ISMIP-HOM tests.  We take one value in each column as 
!      the benchmark and set all other values in that column to the same value.
!      Then we correct the thickness to ensure that lsrf + thck = usrf to double precision.
!      
!      A better way would be to read double-precision data from the input files.

    if (ismip_hom_adjust_usrf) then
       do i = 1, model%general%ewn
          usrf_ref = model%geometry%usrf(i,nhalo+1)
          do j = 1, model%general%nsn
             model%geometry%usrf(i,:) = usrf_ref  ! same usrf for all j in this column
          enddo
       enddo
       model%geometry%thck(:,:) = model%geometry%usrf(:,:) - model%geometry%lsrf(:,:) 
    endif

    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)

    ! initialise model diagnostics                                                                                                                 \
                                                                                                                                                  
    call glide_init_diag(model)

     !WHL - debug
    if (test_halo) then
       call glissade_test_halo (model)
       call parallel_finalise
    endif
     
    if (test_transport) then
       where (model%geometry%thck > model%numerics%thklim)
          model%geometry%thck = thk_init/thk0
       elsewhere
          model%geometry%thck = 0.d0
       endwhere
    endif


    ! Initial solve of calcbwat
    ! TODO sort out if this should go here or in diagnostic solve routine, and make sure consistent with glide.
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

  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model, time, no_write)

    ! Perform time-step of an ice model instance with glissade dycore
    !TODO - Reorganize to put isostasy and calving at start of step?

    use parallel

    use glimmer_paramets, only: tim0, len0, vel0, thk0
    use glimmer_physcon, only: scyr
    use glissade_temp, only: glissade_temp_driver
    use glide_mask, only: glide_set_mask, calc_iareaf_iareag
    use glide_ground, only: glide_marinlim
    use glide_grid_operators
    use isostasy
    use glissade_enthalpy
    use glissade_transport, only: glissade_transport_driver, glissade_check_cfl
    use glissade_grid_operators
    use glide_thck, only: glide_calclsrf
    use glide_bwater

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp),  intent(in)   :: time         ! Current time in years

    !TODO - Remove this argument; it is not used 
    logical, optional, intent(in) :: no_write

    logical nw

    ! --- Local Variables ---

    integer :: sc  ! subcycling index

    ! temporary thck array in SI units (m)
    !TODO - SCALING - Remove if and when scaling is removed
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_unscaled

    ! temporary acab array in SI units (m/s)
    real(dp), dimension(model%general%ewn,model%general%nsn) :: acab_unscaled

    ! temporary variables needed to reset geometry for the EVOL_NO_THICKNESS option
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_old
    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) :: stagthck_old

    ! temporary bmlt array
    real(dp), dimension(model%general%ewn,model%general%nsn) :: &
       bmlt_continuity  ! = bmlt if basal mass balance is included in continuity equation
                        ! else = 0

    logical :: do_upwind_transport  ! Logical for whether transport code should do upwind transport or incremental remapping

    !WHL - debug
    integer :: i, j

    ! ========================

    ! Update internal clock
    model%numerics%time = time  
    model%numerics%timecounter = model%numerics%timecounter + 1
    model%temper%newtemps = .false.

!WHL - debug     
    if (test_transport) then
       call glissade_test_transport (model)
       return
    endif

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    !TODO Not sure if this is in the right place.  G1=f(G0,T0) and T1=g(G0,T0)  
    !     If we update G1 now, then we will be doing T1=g(G1,T0).
    if (model%options%gthf == GTHF_COMPUTE) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! Vertical diffusion and strain heating only; no advection
    ! ------------------------------------------------------------------------ 

    ! Note: These times have units of years.

    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

!HALO TODO - Modify glissade_temp_driver to compute over locally owned cells only?

      call t_startf('glissade_temp_driver')
      call glissade_temp_driver(model, model%options%whichtemp)
      call t_stopf('glissade_temp_driver')

       model%temper%newtemps = .true.

      ! Update hydrology, if needed
      call calcbwat( model,                                    &
                     model%options%whichbwat,                  &
                     model%temper%bmlt,                        &
                     model%temper%bwat,                        &
                     model%temper%bwatflx,                     &
                     model%geometry%thck,                      &
                     model%geometry%topg,                      &
                     model%temper%temp(model%general%upn,:,:), &
                     GLIDE_IS_FLOAT(model%geometry%thkmask),   &
                     model%tempwk%wphi)

    end if

    ! ------------------------------------------------------------------------ 
    ! Halo updates
    ! ------------------------------------------------------------------------ 

    call parallel_halo(model%temper%bwat)    !HALO: not sure if this is needed

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
    ! MJH: This now uses velocity from the previous time step, which is appropriate for a Forward Euler time-stepping scheme
    ! WHL: We used to have EVOL_NO_THICKNESS = -1 as a Glide option, used to hold the ice surface elevation fixed during CESM runs.  
    !      This option has been replaced by a Glint option, evolve_ice.
    !      We now have EVOL_NO_THICKESS = 5 as a glam/glissade option.  It is used to hold the ice surface elevation fixed
    !       while allowing temperature to evolve, which can be useful for model spinup.  This option might need more testing.
    !TODO - Is select case needed here?  There is only one case statement.

    select case(model%options%whichevol)

       case(EVOL_INC_REMAP, EVOL_UPWIND, EVOL_NO_THICKNESS) 

       if (model%options%whichevol == EVOL_UPWIND) then
          do_upwind_transport = .true.
       else
          do_upwind_transport = .false.
       endif

       ! Use incremental remapping scheme for advecting ice thickness ---
       ! (and temperature too, if whichtemp = TEMP_PROGNOSTIC)
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       ! (but not required) to use IR to advect temperature when thickness evolution is turned off.

       ! TODO  MJH If we really want to support no evolution, then we may want to implement it so that IR does not occur 
       !       at all - right now a run can fail because of a CFL violation in IR even if evolution is turned off.  Do we want
       !       to support temperature evolution without thickness evolution?  If so, then the current implementation may be preferred approach.

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! store old thickness
          thck_old = model%geometry%thck
          stagthck_old = model%geomderv%stagthck
       endif

      call t_startf('inc_remap_driver')

       if (main_task) then
          print *, ' '
          print *, 'Compute dH/dt'
       endif

         !WHL - Testing a new subroutine that updates all the key scalars (thck, temp, etc.) at once
         !TODO - Do we need updates of lsrf, usrf, or topg?

!!       call parallel_halo_scalars(model%geometry%thck,   &
!!                                  model%temper%temp)


!       Velocity update may be needed if velo was not updated in halo at the end of the previous diagnostic solve
!        (just to be on the safe side).

        ! Halo updates for velocities, thickness and tracers
      call t_startf('new_remap_halo_upds')

       call staggered_parallel_halo(model%velocity%uvel)
       call staggered_parallel_halo(model%velocity%vvel)

       call parallel_halo(model%geometry%thck)

       if (model%options%whichtemp == TEMP_PROGNOSTIC) then
          call parallel_halo(model%temper%temp)
       endif

       if (model%options%whichtemp == TEMP_ENTHALPY) then
          call parallel_halo(model%temper%temp)
          call parallel_halo(model%temper%waterfrac)
       endif

      call t_stopf('new_remap_halo_upds')

      call t_startf('glissade_transport_driver')


       !TODO  It would be less confusing to just store the subcycling dt in a local/module variable - 
       !       really only needs to be calculated once on init

       model%numerics%dt = model%numerics%dt / model%numerics%subcyc

       if (model%options%basal_mbal == BASAL_MBAL_CONTINUITY) then    ! include bmlt in continuity equation
          bmlt_continuity(:,:) = model%temper%bmlt(:,:) * thk0/tim0   ! convert to m/s
       else                                                           ! do not include bmlt in continuity equation
          bmlt_continuity(:,:) = 0.d0
       endif

      ! --- First determine CFL limits ---
      ! Note we are using the subcycled dt here (if subcycling is on).
      ! (see note above about the EVOL_NO_THICKNESS option and how it is affected by a CFL violation)
      ! stagthck, dusrfdew/ns and u/vvel need to be from the previous time step (and are at this point)
      call glissade_check_cfl(model%general%ewn, model%general%nsn, model%general%upn-1,                   &
                           model%numerics%dew * len0, model%numerics%dns * len0, model%numerics%sigma,     &
                           model%geomderv%stagthck * thk0, model%geomderv%dusrfdew*thk0/len0, model%geomderv%dusrfdns*thk0/len0, &
                           model%velocity%uvel * scyr * vel0, model%velocity%vvel * scyr * vel0,           &
                           model%numerics%dt * tim0 / scyr,                                                &
                           model%numerics%adv_cfl_dt, model%numerics%diff_cfl_dt )

       ! Call the transport driver.
       ! Note: This subroutine assumes SI units:
       !       * dt (s)
       !       * dew, dns, thck (m)
       !       * uvel, vvel, acab, blmt (m/s)
       !       Since thck has intent(inout), we create and pass a temporary array with units of m.

       !TODO: Test this driver for 1st order upwind transport.

       do sc = 1 , model%numerics%subcyc
          if (model%numerics%subcyc > 1) write(*,*) 'Subcycling transport: Cycle ',sc

          ! temporary in/out arrays in SI units (m)                               
          thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0
          acab_unscaled(:,:) = model%climate%acab(:,:) * thk0/tim0

          if (model%options%whichtemp == TEMP_PROGNOSTIC) then  ! Use IR to transport thickness, temperature
                                                                ! (and other tracers, if present)
                                                                ! Note: We are passing arrays in SI units.

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            ntracer,                                              &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            acab_unscaled(:,:),                                   &
                                            bmlt_continuity(:,:),                                 &
                                            model%temper%temp(:,:,:),                             &
                                            upwind_transport_in = do_upwind_transport )

             ! convert thck and acab back to scaled units
             model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0
             model%climate%acab(:,:) = acab_unscaled(:,:) / (thk0/tim0)

          elseif (model%options%whichtemp == TEMP_ENTHALPY) then  ! Use IR to transport thickness, temperature,
                                                                ! and waterfrac.  Also set ntracer = 3
                                                                ! Note: We are passing arrays in SI units.

             ! BDM Set ntracer = 3 since temp and waterfrac need to be passed (and maybe ice age).
             !     If no ice age is present, a dummy array will be passed for tracer = 2
             ntracer = 3

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      & 
                                            ntracer,                                              &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            acab_unscaled(:,:),                                   &
                                            bmlt_continuity(:,:),                                 &
                                            model%temper%temp(:,:,:),                             & 
                                            model%temper%waterfrac(:,:,:),                        &
                                            upwind_transport_in = do_upwind_transport )

          else  ! Use IR to transport thickness only
                ! Note: In glissade_transport_driver, the ice thickness is transported layer by layer,
                !       which is inefficient if no tracers are being transported.  (It would be more
                !       efficient to transport thickness in one layer only, using a vertically
                !       averaged velocity.)  Not sure if this option will be used in practice.

             !TODO - Will we continue to support this option?

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            ntracer,                                              &
                                            model%velocity%uvel(:,:,:) * vel0,                    & 
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            thck_unscaled(:,:),                                   &
                                            acab_unscaled(:,:),                                   &
                                            bmlt_continuity(:,:) ,                                &
                                            upwind_transport_in = do_upwind_transport )

          endif  ! whichtemp

          ! convert thck and acab back to scaled units
          model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0
          model%climate%acab(:,:) = acab_unscaled(:,:) / (thk0/tim0)
         
!WHL - debug
!    print*, ' '
!    print*, 'After glissade_transport_driver:'
!    print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
!    print*, 'max, min acab (m/yr) =', maxval(model%climate%acab)*scale_acab, minval(model%climate%acab)*scale_acab
!    print*, 'max, min artm =', maxval(model%climate%artm), minval(model%climate%artm)
!    print*, 'thklim =', model%numerics%thklim * thk0
!    print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
!    print*, ' '
!    print*, 'thck:'
!    do j = model%general%nsn, 1, -1
!       write(6,'(23f7.2)') model%geometry%thck(7:29,j)*thk0
!    enddo

          ! Update halos of modified fields

         call t_startf('after_remap_haloupds')

         call parallel_halo(model%geometry%thck)
         call parallel_halo(model%temper%temp)

          ! Halo updates of other tracers, if present, would need to go here
         if (model%options%whichtemp == TEMP_ENTHALPY) then
            call parallel_halo(model%temper%waterfrac)
         endif

         call t_stopf('after_remap_haloupds')

       enddo     ! subcycling

!TODO: Don't divide and multiply this variable
       model%numerics%dt = model%numerics%dt * model%numerics%subcyc

      call t_stopf('glissade_transport_driver')

      call t_stopf('inc_remap_driver')

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck = thck_old
          model%geomderv%stagthck = stagthck_old
       endif

    end select

    call parallel_halo(model%geometry%topg)

    ! ------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    ! ------------------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

!       TODO: Make sure a call to calc_flwa (based on post-remap temperature) 
!             is not needed for glide_marinlim.  Should be based on post-remap temperature.

    ! --- Calculate updated mask because marinlim calculation needs a mask.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask)

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

!HALO - Look at marinlim more carefully and see which fields need halo updates before it is called.
!       It appears that marinlim only needs the halo of thkmask for case 5 (which was removed).  
!       If that case is removed, a thkmask halo update does not need to occur here.

    ! TODO: glide_set_mask includes a halo update of model%geometry%thkmask; move it here?
    call parallel_halo(model%geometry%thkmask) 
    call parallel_halo(model%isostasy%relx)

    call glide_marinlim(model%options%whichmarn, &
         model%geometry%thck,  &
         model%isostasy%relx,      &
         model%geometry%topg,  &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving,  &
         model%ground, &
         model%numerics%dew,    &
         model%numerics%dns, &
         model%general%nsn, &
         model%general%ewn)
         ! model%geometry%usrf) not used in routine

!HALO TODO: Test the effect of these updates with a nonzero calving field

    ! halo updates

    call parallel_halo(model%geometry%thck)    ! Updated halo values of thck are needed below in calc_lsrf

!TODO - Remove this call to glide_set_mask?
!      This subroutine is called at the beginning of glissade_velo_driver,
!       so a call here is not needed for the velo diagnostic solve.
!      The question is whether it is needed for the isostasy.
!      And the isostasy may itself be in the wrong place.

    ! --- marinlim adjusts thickness for calved ice.  Therefore the mask needs to be recalculated.
    ! --- This time we want to calculate the optional arguments iarea and ivol because thickness 
    ! --- will not change further during this time step.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

!HALO TODO - glide_set_mask includes a halo update of model%geometry%thkmask at end of call
!       That update should be moved here if needed later (but may not be needed).

    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call

    ! --- Calculate area of ice that is floating and grounded.
    !TODO This subroutine does not use iarea - remove from the call/subroutine.
    !TODO May want to only calculate iarea, iareaf, iareag in glide_write_diag() and remove those calculations here.  

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%iarea,  model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    !TODO - Need a global sum here? (currently done inside calc_iareaf_iareag)

!TODO These isostasy calls may be in the wrong place.  
! Consider for a forward Euler time step:
! With a relaxing mantle model, topg is a prognostic (time-evolving) variable (I think):
! topg1 = f(topg0, thk0, ...) 
! However, for a fluid mantle where the adjustment is instantaneous, topg is a diagnostic variable 
!(comparable to calculating floatation height of ice in the ocean):
! topg1 = f(thk1)
! In either case, the topg update should be separate from the thickness evolution (because thk1 = f(thk0, vel0=g(topg0,...)).
! However, if the isostasy calculation needs topg0, the icewaterload call should be made BEFORE thck is updated.  
! If the isostasy calculation needs topg1, the icewaterload call should be made AFTER thck is updated.  
! Also, we need think carefully about when marinlim, usrf, lsrf, derivatives should be calculated relative to the topg update via isostasy.

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       if (model%numerics%time >= model%isostasy%next_calc) then
          model%isostasy%next_calc = model%isostasy%next_calc + model%isostasy%period
          call isos_icewaterload(model)
          model%isostasy%new_load = .true.
       end if
    end if
   
      ! calculate isostatic adjustment and upper and lower ice surface

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    !TODO - Test the local isostasy schemes in the parallel model.
    !       The elastic lithosphere scheme is not expected to work in parallel.

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call isos_compute(model)
    end if

    ! ------------------------------------------------------------------------
    ! Calculate diagnostic variables, including velocity
    ! ------------------------------------------------------------------------

    call glissade_diagnostic_variable_solve(model)

    !TODO - Any halo updates needed at the end?  

  end subroutine glissade_tstep

!=======================================================================
!MJH added this diagnostic solve subroutine so it can be called from init.  

  subroutine glissade_diagnostic_variable_solve(model) 

     ! Solve diagnostic (not time-dependent) variables.  This is needed at the end of each time step once the 
     !  prognostic variables (thickness, tracers) have been updated.  
     ! It is also needed to fill out the initial state from the fields that have been read in.

    use parallel

    use glimmer_paramets, only: tim0, len0, vel0
    use glimmer_physcon, only: scyr
    use glide_thck, only: glide_calclsrf
    use glissade_temp, only: glissade_calcflwa
    use glam_velo, only: glam_velo_driver
    use glissade_velo, only: glissade_velo_driver
    use glide_velo, only: wvelintg

    use glam_grid_operators, only: glam_geometry_derivs, stagthickness
    use felix_dycore_interface, only: felix_velo_driver

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

!WHL - debug
    integer :: i, j
    logical, parameter :: verbose_glissade = .false.

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 1. First part of diagnostic solve: 
    !    Now that advection is done, update geometry- and temperature-related 
    !    diagnostic fields that are needed for the velocity solve.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    !WHL - Moved this calculation here from glissade_temp, since flwa is a diagnostic variable.

    ! Calculate Glen's A --------------------------------------------------------
    !
    ! Note: because flwa is not a restart variable in glissade, no check is included 
    !       here for whether to calculate it on initial time (as is done in glide).
    ! 
    ! Note: We are passing in only vertical elements (1:upn-1) of the temp array,
    !       so that it has the same vertical dimensions as flwa.

    ! BDM - adding a call for whichtemp = TEMP_ENTHALPY that includes waterfrac as input
    ! TODO - May be OK to use a single call (with waterfrac) for either TEMP option.

    if (model%options%whichtemp == TEMP_ENTHALPY) then

       call glissade_calcflwa(model%numerics%stagsigma,    &
                              model%numerics%thklim,       &
                              model%temper%flwa,           &
                              model%temper%temp(1:model%general%upn-1,:,:),  &
                              model%geometry%thck,         &
                              model%paramets%flow_factor,  &
                              model%paramets%default_flwa, &
                              model%options%whichflwa,      &
                              model%temper%waterfrac(:,:,:))
    else

       call glissade_calcflwa(model%numerics%stagsigma,    &
                              model%numerics%thklim,       &
                              model%temper%flwa,           &
                              model%temper%temp(1:model%general%upn-1,:,:),  &
                              model%geometry%thck,         &
                              model%paramets%flow_factor,  &
                              model%paramets%default_flwa, &
                              model%options%whichflwa)
    endif

    ! Halo update for flwa
    call parallel_halo(model%temper%flwa)

    ! ------------------------------------------------------------------------
    ! Halo updates for ice topography and thickness
    !
    !WHL - Note the optional argument periodic_offset_ew for topg.
    !      This is for ismip-hom experiments. A positive EW offset means that 
    !       the topography in west halo cells will be raised, and the topography 
    !       in east halo cells will be lowered.  This ensures that the topography
    !       and upper surface elevation are continuous between halo cells
    !       and locally owned cells at the edge of the global domain.
    !      In other cases (anything but ismip-hom), periodic_offset_ew = periodic_offset_ns = 0, 
    !       and this argument will have no effect.
    ! ------------------------------------------------------------------------

    call parallel_halo(model%geometry%thck)
    call parallel_halo(model%geometry%topg, periodic_offset_ew = model%numerics%periodic_offset_ew)

    ! ------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    ! ------------------------------------------------------------------------

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    ! ------------------------------------------------------------------------
    ! Update some geometry derivatives
    ! ------------------------------------------------------------------------
    !TODO - The fields computed by glam_geometry_derivs are not required by
    !       the glissade solver (which computes them internally).  
    !       However, some of the fields (stagthck, dusrfdew and dusrfdns) 
    !       are needed during the next timestep by glissade_temp
    !       if we're doing shallow-ice dissipation.  If dissipation is always
    !       higher-order, then I don't think this call is needed here.
    !       (The glam_velo driver includes its own call to glam_geometry_derivs.) 

    call glam_geometry_derivs(model)


    !WHL - Moved glam-specific geometry calculations to glam_velo_driver in glam_velo.F90.

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 2. Second part of diagnostic solve: 
    !    Now that geometry- and temperature-related diagnostic fields are updated, 
    !    solve velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Do not solve velocity for initial time on a restart because that breaks an exact restart.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model % numerics % time == model % numerics % tstart) ) then
  
       call write_log('Using uvel, vvel from restart file at initial time')

    else

       ! If this is not a restart or we are not at the initial time, then proceed normally.

       if ( (model % numerics % time == model % numerics % tstart) .and. &
         ( (maxval(abs(model%velocity%uvel))/=0.0d0) .or. & 
           (maxval(abs(model%velocity%vvel))/=0.0d0) ) ) then
          ! If velocity was input and this is NOT a restart, then use the input field as the first guess at the initial time.
          ! This happens automatically, but let the user know.
          ! Using this value versus not will only change the answer within the tolerance of the nonlinear solve.  
          ! If a user already has a good guess from a previous run, they may wish to start things off with it to speed the initial solution.
          call write_log('Using uvel, vvel from input file as initial guess at initial time.')
          call write_log('If this is not desired, please remove those fields from the input file.')
       endif

       if (main_task) then
          print *, ' '
          print *, 'Compute ice velocities, time =', model%numerics%time
       endif

       !! extrapolate value of mintauf into halos to enforce periodic lateral bcs (only if field covers entire domain)
       if (model%options%which_ho_babc == HO_BABC_YIELD_PICARD) then
          call staggered_parallel_halo_extrapolate(model%basalproc%mintauf)
       endif

       !WHL - Broke up into separate calls to glam_velo_driver and glissade_velo_driver
       !      Previously had a single call to glissade_velo_driver
       ! MJH - modified Bill's organization to add felix, as well.

       ! Call the appropriate velocity solver
       select case (model%options%whichdycore)
       case ( DYCORE_GLAM )   ! glam finite-difference

          call t_startf('glam_velo_driver')
          call glam_velo_driver(model)
          call t_stopf('glam_velo_driver')

       case ( DYCORE_GLISSADE )   ! glissade finite-element

          call t_startf('glissade_velo_driver')
          call glissade_velo_driver(model)
          call t_stopf('glissade_velo_driver')

       case ( DYCORE_ALBANYFELIX)

          call t_startf('felix_velo_driver')
          call felix_velo_driver(model)
          call t_stopf('felix_velo_driver')

       end select
 
!WHL - debug
       if (main_task .and. verbose_glissade) then

          print*, ' '
          print*, 'After glissade velocity solve: uvel, k = 1:'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.2)',advance='no') model%velocity%uvel(1,i,j) * (vel0*scyr)
             enddo 
             print*, ' '
          enddo

          print*, ' '
          print*, 'After glissade velocity solve: vvel, k = 1:'
          do i = 1, model%general%ewn-1
             write(6,'(i8)',advance='no') i
          enddo
          print*, ' '
          do j = model%general%nsn-1, 1, -1
             write(6,'(i4)',advance='no') j
             do i = 1, model%general%ewn-1
                write(6,'(f8.2)',advance='no') model%velocity%vvel(1,i,j) * (vel0*scyr)
             enddo 
             print*, ' '
          enddo

       endif  ! main_task

    endif     ! is_restart

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 3. Third part of diagnostic solve: 
    ! Now that velocity is solved, calculate any diagnostic fields that are
    ! a function of velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! compute the velocity norm (for diagnostic output)

    model%velocity%velnorm = sqrt(model%velocity%uvel**2 + model%velocity%vvel**2)

    ! WHL - Copy uvel and vvel to arrays uvel_icegrid and vvel_icegrid.
    !       These arrays have horizontal dimensions (nx,ny) instead of (nx-1,ny-1).
    !       Thus they are better suited for I/O if we have periodic BC,
    !        where the velocity field we are solving for has global dimensions (nx,ny).
    !       Since uvel and vvel are not defined for i = nx or for j = ny, the
    !        uvel_icegrid and vvel_icegrid arrays will have values of zero at these points.
    !       But these are halo points, so when we write netCDF I/O it shouldn't matter;
    !        we should have the correct values at physical points.
    
    model%velocity%uvel_icegrid(:,:,:) = 0.d0
    model%velocity%vvel_icegrid(:,:,:) = 0.d0

    do j = 1, model%general%nsn-1
       do i = 1, model%general%ewn-1
          model%velocity%uvel_icegrid(:,i,j) = model%velocity%uvel(:,i,j)
          model%velocity%vvel_icegrid(:,i,j) = model%velocity%vvel(:,i,j)             
       enddo
    enddo
 
    ! Calculate wvel, assuming grid velocity is 0.
    ! This calculated wvel relative to ice sheet base, rather than a fixed reference location
    ! Note, this current implementation for wvel only supports whichwvel=VERTINT_STANDARD
    call wvelintg(model%velocity%uvel,                            &
                      model%velocity%vvel,                        &
                      model%geomderv,                             &
                      model%numerics,                             &
                      model%velowk,                               &
                      model%geometry%thck * 0.0d0,                &  ! Just need a 2d array of all 0's for wgrd
                      model%geometry%thck,                        &
                      model%temper%bmlt,                          &
                      model%velocity%wvel_ho)
    ! Note: halos may be wrong for wvel_ho, but since it is currently only used as an output diagnostic variable, that is fine.

!TODO - Don't think we need to update ubas, vbas, or velnorm,
!       because these can be derived directly from the 3D uvel and vvel arrays

    call staggered_parallel_halo(model%velocity%velnorm)
    call staggered_parallel_halo(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)
    call parallel_halo(model%stress%efvs)

     !WHL - Removed this call because tau is now computed more accurately in glissade_velo_higher.F90.
!!    call glide_calcstrsstr( model )       !*sfp* added for populating stress tensor w/ HO fields

  end subroutine glissade_diagnostic_variable_solve

!=======================================================================

  !TODO: Not currently called; may not be worth having a separate subroutine for this

  subroutine parallel_halo_scalars(thck,     temp,   &
                                   lsrf,     usrf,   &
                                   topg,     tracers)

    ! Do parallel halo updates for the main scalar state variables

    use parallel

    real(dp), intent(inout), dimension(:,:) :: thck   
    real(dp), intent(inout), dimension(:,:,:) :: temp

    real(dp), intent(inout), dimension(:,:), optional ::  &
       lsrf,       & ! lower ice surface
       usrf,       & ! upper ice surface
       topg          ! basal topography

    real(dp), intent(inout), dimension(:,:,:,:), optional :: tracers

    integer :: nt   ! tracer index

    call parallel_halo(thck)
    call parallel_halo(temp)

    ! optional updates for geometry variables

    if (present(lsrf)) call parallel_halo(lsrf)

    if (present(usrf)) then
       if (present(lsrf)) then   ! compute usrf = lsrf + thck; no halo call needed
          usrf(:,:) = lsrf(:,:) + thck(:,:)
       else
          call parallel_halo(usrf)
       endif
    endif

    if (present(topg)) call parallel_halo(topg)

    ! optional update for 3D tracers (e.g., ice age)

    if (present(tracers)) then

       do nt = 1, ntracer
          call parallel_halo(tracers(:,:,:,nt))
       enddo

    endif

  end subroutine parallel_halo_scalars

!=======================================================================

    subroutine glissade_test_halo(model)

    use parallel

    ! various tests of parallel halo updates

    type(glide_global_type), intent(inout) :: model      ! model instance

    integer, dimension (:,:), allocatable    ::  pgID       ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDr4     ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDr8     ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDr8_3d  ! unique global ID for parallel runs  

    integer, dimension (:,:), allocatable    ::  pgIDstagi  ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable   ::  pgIDstagr  ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDstagr3 ! unique global ID for parallel runs  

    logical, dimension(:,:), allocatable    :: logvar
    integer, dimension(:,:), allocatable    :: intvar
    real, dimension(:,:), allocatable       :: r4var
    real(dp), dimension(:,:), allocatable   :: r8var
    real(dp), dimension(:,:,:), allocatable :: r8var_3d

    integer :: i, j, k, n
    integer :: nx, ny, nz

    integer, parameter :: rdiag = 0         ! rank for diagnostic prints 

!    real(dp) :: global_row, global_col, global_ID

    print*, ' '
    print*, 'In glissade_test_halo, this_rank =', this_rank

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn

    allocate(logvar(nx,ny))
    allocate(intvar(nx,ny))
    allocate(r4var(nx,ny))
    allocate(r8var(nx,ny))
    allocate(r8var_3d(nz,nx,ny))

    allocate(pgID(nx,ny))
    allocate(pgIDr4(nx,ny))
    allocate(pgIDr8(nx,ny))
    allocate(pgIDr8_3d(nz,nx,ny))
    allocate(pgIDstagi(nx-1,ny-1))
    allocate(pgIDstagr(nx-1,ny-1))
    allocate(pgIDstagr3(nz,nx-1,ny-1))

    if (main_task) then
       print*, ' '
       print*, 'nx, ny, nz =', nx, ny, nz
       print*, 'uhalo, lhalo =', uhalo, lhalo
       print*, 'global_ewn, global_nsn =', global_ewn, global_nsn
       print*, ' '
    endif

    print*, 'this_rank, global_row/col offset =', this_rank, global_row_offset, global_col_offset

    ! Test the 5 standard parallel_halo routines for scalars: logical_2d, integer_2d, real4_2d, real8_2d, real8_3d

    ! logical 2D field

    logvar(:,:) = .false.

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       logvar(i,j) = .true.
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'logvar: this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35l3)') logvar(:,j)
       enddo
    endif

    call parallel_halo(logvar)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35l3)') logvar(:,j)
       enddo
    endif

!WHL - Skip the next few sections since the global IDs aren't correct in parallel.
!TODO - Fix this up or remove?
    go to 100

    ! integer 2D field

    intvar(:,:) = 1

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       intvar(i,j) = (nx-2*nhalo)*(j-nhalo-1) + i-nhalo
    enddo
    enddo 

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'intvar: this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35i5)') intvar(:,j)
       enddo
    endif

    call parallel_halo(intvar)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update:'
       do j = ny, 1, -1
          write(6,'(35i5)') intvar(:,j)
       enddo
    endif

    ! real 2D
    
    r4var(:,:) = 1.

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       r4var(i,j) = (nx-2.*real(nhalo))*(real(j)-real(nhalo)-1.) + real(i)-real(nhalo)
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'r4var: this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') r4var(:,j)
       enddo
    endif

    call parallel_halo(r4var)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update:'
       do j = ny, 1, -1
          write(6,'(35f6.0)') r4var(:,j)
       enddo
    endif

    ! double 2D
    
    r8var(:,:) = 1.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       r8var(i,j) = (real(nx,dp)-2.d0*real(nhalo,dp))*(real(j,dp)-real(nhalo,dp)-1.d0) + real(i,dp)-real(nhalo,dp)
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'r8var: this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') r8var(:,j)
       enddo
    endif

    call parallel_halo(r8var)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update:'
       do j = ny, 1, -1
          write(6,'(35f6.0)') r8var(:,j)
       enddo
    endif

    ! double 3D

    r8var_3d(:,:,:) = 1.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          r8var_3d(k,i,j) = (real(nx,dp)-2.d0*real(nhalo,dp))*(real(j,dp)-real(nhalo,dp)-1.d0) + real(i,dp)-real(nhalo,dp)
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'r8var_3d: this_rank, k =', this_rank, k
       do j = ny, 1, -1
          write(6,'(35f6.0)') r8var_3d(k,:,j)
       enddo
    endif

    call parallel_halo(r8var_3d)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update:'
       do j = ny, 1, -1
          write(6,'(35f6.0)') r8var_3d(k,:,j)
       enddo
    endif

100 continue

    ! The next part of the code concerns parallel global IDs

    ! Compute parallel global ID for each grid cell

    pgID(:,:) = 0   ! integer

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgID(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (integer), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35i5)') pgID(:,j)
       enddo
    endif

    call parallel_halo(pgID)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35i5)') pgID(:,j)
       enddo
    endif

    ! real 2D
    
    pgIDr4(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr4(i,j) = real(parallel_globalID_scalar(i,j,nz))
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (r4 2D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') pgIDr4(:,j)
       enddo
    endif

    call parallel_halo(pgIDr4)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') pgIDr4(:,j)
       enddo
    endif

    ! double 2D
    
    pgIDr8(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr8(i,j) = real(parallel_globalID_scalar(i,j,nz), dp)
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (r8 2D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') pgIDr8(:,j)
       enddo
    endif

    call parallel_halo(pgIDr8)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') pgIDr8(:,j)
       enddo
    endif

    ! double 3D

    pgIDr8_3d(:,:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDr8_3d(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (real 3D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') pgIDr8_3d(k,:,j)
       enddo
    endif

    call parallel_halo(pgIDr8_3d)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(35f6.0)') pgIDr8_3d(k,:,j)
       enddo
    endif

    ! Repeat for staggered variables

    ! First for an integer 2D field

    pgIDstagi(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagi(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    ! Print
    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (integer), this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(34i5)') pgIDstagi(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagi)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(34i5)') pgIDstagi(:,j)
       enddo
    endif

    ! Then for a real 2D field

    pgIDstagr(:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDstagr(i,j) = real(parallel_globalID_scalar(i,j,nz),dp)    ! function in parallel_mpi.F90
    enddo
    enddo

    ! Print
    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (real 2D), this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(34f6.0)') pgIDstagr(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(34f6.0)') pgIDstagr(:,j)
       enddo
    endif

    ! Then for a real 3D field

    pgIDstagr3(:,:,:) = 0.d0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDstagr3(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Staggered parallel global ID (real 3D), k, this_rank =', k, this_rank
       do j = ny-1, 1, -1
          write(6,'(34f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr3)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(34f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    deallocate(logvar)
    deallocate(intvar)
    deallocate(r4var)
    deallocate(r8var)
    deallocate(r8var_3d)
    deallocate(pgID)
    deallocate(pgIDr4)
    deallocate(pgIDr8)
    deallocate(pgIDr8_3d)
    deallocate(pgIDstagi)
    deallocate(pgIDstagr)
    deallocate(pgIDstagr3)

    end subroutine glissade_test_halo

!=======================================================================

    subroutine glissade_test_transport(model)

    use parallel
    use glissade_transport, only: glissade_transport_driver
    use glimmer_paramets, only: len0, thk0, tim0
    use glimmer_physcon, only: pi, scyr

    use glide_diagnostics

    !-------------------------------------------------------------------
    ! Test transport of a cylinder or block of ice once around the domain and
    ! back to the starting point, assuming uniform motion in a straight line.
    !
    ! Instructions for use:
    ! (1) At the top of this module, set test_transport = .true. and choose a
    !     value for thk_init.
    ! (2) Set the config file to run with the glissade dycore for one timestep, 
    !     with CF output written every timestep. 
    !     Note: Whatever the initial ice geometry in the config file
    !     (e.g., the dome test case), the ice extent will be preserved, but
    !     the thickness will be set everywhere to thk_init, giving a steep
    !     front at the ice margin that is challenging for transport schemes.
    ! (3) Comment out the call to glissade_diagnostic_variable_solve in
    !     cism_driver or simple_glide.  (It probable won't hurt to
    !     have it turned on, but will just use extra cycles.)
    !
    ! During the run, the following will happen:
    ! (1) During glissade_initialise, the ice thickness will be set to
    !     thk_init (everywhere there is ice).
    ! (2) During the first call to glissade_tstep, this subroutine 
    !     (glissade_test_transport) will be called.  The ice will be transported 
    !     around to its initial starting point (assuming periodic BCs).  
    ! (3) Then the model will return from glissade_tstep before doing any 
    !     other computations. Output will be written to netCDF and the code 
    !     will complete.
    !
    ! There should be two time slices in the output netCDF file.
    ! Compare the ice geometry at these two slices to see how diffusive 
    !  and/or dispersive the transport scheme is. A perfect scheme
    !  would leave the geometry unchanged.  
    !-------------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model      ! model instance

    real(dp), dimension(:,:,:), allocatable :: uvel, vvel   ! uniform velocity field (m/yr)

    integer :: i, j, k, n
    integer :: nx, ny, nz
    real(dp) :: dx, dy

    integer, parameter :: rdiag = 0         ! rank for diagnostic prints 

!    real(dp), parameter :: umag = 100.     ! uniform speed (m/yr)
    real(dp), parameter :: umag = 1000.     ! uniform speed (m/yr)

    ! Set angle of motion
    !WHL - Tested all three of these angles (eastward, northward, and northeastward)
!    real(dp), parameter :: theta = 0.d0      ! eastward
!     real(dp), parameter :: theta = pi/4.d0   ! northeastward
!    real(dp), parameter :: theta = pi/2.d0   ! northward
    real(dp), parameter :: theta = pi         ! westward

    real(dp), parameter :: thk = 500.d0

    real(dp) :: dt          ! time step in yr

    integer :: ntstep       ! run for this number of timesteps
    real(dp) :: theta_c    ! angle dividing paths through EW walls from paths thru NS walls
    real(dp) :: lenx       ! length of shortest path thru EW walls
    real(dp) :: leny       ! length of shortest path thru NS walls
    real(dp) :: len_path   ! length of path back to starting point
    real(dp) :: adv_cfl    ! advective CFL number

    logical :: do_upwind_transport  ! if true, do upwind transport

    ! Initialize

    dx = model%numerics%dew * len0
    dy = model%numerics%dns * len0

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn

    allocate(uvel(nz,nx-1,ny-1))
    allocate(vvel(nz,nx-1,ny-1))

    ! Find the length of the path around the domain and back to the starting point

    lenx = global_ewn * dx
    leny = global_nsn * dy
    theta_c = atan(leny/lenx)   ! 0 <= theta_c <= pi/2

    if ( (theta >= -theta_c   .and. theta <= theta_c) .or.   &
         (theta >= pi-theta_c .and. theta <= pi+theta_c) ) then
       ! path will cross east/west wall
       len_path = lenx / abs(cos(theta))
    else
       ! path will cross north/south wall
       len_path = leny / abs(sin(theta))
    endif

    ! Choose the number of time steps such that the ice will travel
    ! less than one grid cell per time step

    ntstep = len_path/dx + 10   ! 10 is arbitrary

    ! Choose the time step such that the ice will go around the domain
    ! exactly once in the chosen number of time steps.

    dt = len_path / (umag*ntstep)

    ! CFL check, just to be sure

    adv_cfl = max (dt*umag*cos(theta)/dx, dt*umag*sin(theta)/dy)
    
    if (adv_cfl >= 1.d0) then
       print*, 'dt is too big for advective CFL; increase ntstep to', ntstep * adv_cfl
       stop
    endif

    ! Print some diagnostics

    if (main_task) then
       print*, ' '
       print*, 'In glissade_test_transport'
       print*, 'nx, ny, nz =', nx, ny, nz
       print*, 'len_path =', len_path
       print*, 'umag (m/yr) =', umag
       print*, 'dt (yr) =', dt
       print*, 'ntstep =', ntstep
       print*, 'theta (deg) =', theta * 180.d0/pi
    endif

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'Initial thck'
       do j = ny, 1, -1
          write(6,'(19f7.2)') model%geometry%thck(1:19,j) * thk0
       enddo
!          write(6,*) ' '
!          write(6,*) 'New layer 1 temp, n =', n
!          do j = ny, 1, -1
!             write(6,'(19f7.2)') model%temper%temp(1,1:19,j)
!          enddo
    endif

    ! Set uniform ice speed everywhere

    do j = 1, ny-1
    do i = 1, nx-1
       do k = 1, nz
          uvel(k,i,j) = umag * cos(theta)
          vvel(k,i,j) = umag * sin(theta)
       enddo
    enddo
    enddo

    ! Determine which transport scheme

    if (model%options%whichevol == EVOL_UPWIND) then
       do_upwind_transport = .true.
    else
       do_upwind_transport = .false.
    endif

    ! Transport the ice around the domain

    do n = 1, ntstep

       ! call transport scheme
       ! Note: glissade_transport_driver expects dt in seconds, uvel/vvel in m/s

       call glissade_transport_driver(dt*scyr,                                              &
                                      dx,                        dy,                        &
                                      nx,                        ny,                        &
                                      nz-1,                      model%numerics%sigma,      &
                                      ntracer,                                              &
                                      uvel(:,:,:)/scyr,          vvel(:,:,:)/scyr,          &
                                      model%geometry%thck(:,:),                             &
                                      model%climate%acab(:,:),                              &
                                      model%temper%bmlt(:,:),                               &
                                      model%temper%temp(:,:,:),                             &
                                      upwind_transport_in = do_upwind_transport)

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'New thck, n =', n
          do j = ny, 1, -1
             write(6,'(19f7.2)') model%geometry%thck(1:19,j) * thk0
          enddo
!          write(6,*) ' '
!          write(6,*) 'New layer 1 temp, n =', n
!          do j = ny, 1, -1
!             write(6,'(19f7.2)') model%temper%temp(1,1:19,j)
!          enddo
       endif

    enddo  ! ntstep

    if (main_task) print*, 'Done in glissade_test_parallel'

    deallocate(uvel)
    deallocate(vvel)

    end subroutine glissade_test_transport

!=======================================================================


end module glissade
