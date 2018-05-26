!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade.F90 - part of the Community Ice Sheet Model (CISM)  
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
  use glissade_test, only: glissade_test_halo, glissade_test_transport
  use glide_thck, only: glide_calclsrf  ! TODO - Make this a glissade subroutine, or inline

  implicit none

  integer, private, parameter :: dummyunit=99

  logical, parameter :: verbose_glissade = .false.

  ! Change any of the following logical parameters to true to carry out simple tests
  logical, parameter :: test_transport = .false.   ! if true, call test_transport subroutine
  real(dp), parameter :: thk_init = 500.d0         ! initial thickness (m) for test_transport
  logical, parameter :: test_halo = .false.        ! if true, call test_halo subroutine

contains

!=======================================================================

! Note: There is no glissade_config subroutine; glide_config works for all dycores.

!=======================================================================

  subroutine glissade_initialise(model, evolve_ice)

    ! initialise Glissade model instance

    use parallel
    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio
    use glide_velo, only: init_velo  !TODO - Remove call to init_velo?
    use glissade_therm, only: glissade_init_therm
    use glissade_transport, only: glissade_overwrite_acab_mask
    use glissade_basal_water, only: glissade_basal_water_init
    use glissade_masks, only: glissade_get_masks
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glam_strs2, only : glam_velo_init
    use glimmer_coordinates, only: coordsystem_new
    use glissade_grid_operators, only: glissade_stagger
    use glissade_velo_higher, only: glissade_velo_higher_init
    use glide_diagnostics, only: glide_init_diag
    use glissade_calving, only: glissade_calving_mask_init, glissade_calve_ice
    use glimmer_paramets, only: thk0, len0, tim0
    use felix_dycore_interface, only: felix_velo_init

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    logical, intent(in), optional :: evolve_ice       ! whether ice evolution is turned on (if not present, assumed true)

    !TODO - Is glimmer_version_char still needed?
    character(len=100), external :: glimmer_version_char

    character(len=100) :: message

    real(dp) :: var_maxval   ! max value of a given variable; = 0 if not yet read in
    integer :: i, j, k
    logical :: l_evolve_ice  ! local version of evolve_ice
    integer, dimension(:,:), allocatable :: &
         ice_mask,          & ! = 1 where ice is present, else = 0
         floating_mask,     & ! = 1 where ice is present and floating, else = 0
         ocean_mask,        & ! = 1 if topg is below sea level and ice is absent, else = 0
         land_mask            ! = 1 if topg is at or above sea level, else = 0

    integer :: itest, jtest, rtest

    if (present(evolve_ice)) then
       l_evolve_ice = evolve_ice
    else
       l_evolve_ice = .true.
    end if

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters
    call glide_scale_params(model)

    ! set up coordinate systems, and change to the parallel values of ewn and nsn

    ! With outflow BCs, scalars in the halos are set to zero.
    if (model%general%global_bc == GLOBAL_BC_OUTFLOW) then
       call distributed_grid(model%general%ewn, model%general%nsn, outflow_bc_in=.true.)
    elseif (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
       ! Note: In this case, halo updates are treated the same as for periodic BC
       !       In addition, we set masks that are used to zero out the outflow velocities in the Glissade velocity solvers
       !       The no-penetration BC is not supported for the Glam solver
       call distributed_grid(model%general%ewn, model%general%nsn)
    else
       call distributed_grid(model%general%ewn, model%general%nsn)
    endif

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0,               &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.d0, model%numerics%dns/2.d0, &
                                              model%numerics%dew,      model%numerics%dns,      &
                                              model%general%ewn-1,     model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

    ! set masks at global boundary for no-penetration boundary conditions
    ! this subroutine includes a halo update
    if (model%general%global_bc == GLOBAL_BC_NO_PENETRATION) then
       call staggered_no_penetration_mask(model%velocity%umask_no_penetration, model%velocity%vmask_no_penetration)
    endif

    ! set uniform basal heat flux (positive down)
    model%temper%bheatflx = model%paramets%geot

    ! compute sigma levels or load from external file
    ! (if not already read from config file)
    call glide_load_sigma(model,dummyunit)

    ! initialize the time step counter
    ! For restart, tstep_count will be overwritten from the restart file.
    ! Alternatively, we could initialize tstep_count as follows:
    !    model%numerics%tstep_count = nint(model%numerics%time/model%numerics%tinc)
    ! But reading from a restart file should prevent roundoff issues.
    model%numerics%tstep_count = 0

    ! open all input files
    call openall_in(model)

    ! read first time slice
    call glide_io_readall(model,model)

    ! Optionally, compute area scale factors for stereographic map projection.
    ! This should be done after reading the input file, in case the input file contains mapping info.
    ! Note: Not yet enabled for other map projections.
    ! Note: Area factors currently are used only when computing diagnostics for the log file.
    !       If area factors are computed based on projection info, then the ice area and volume
    !         computed in CISM's dycore are corrected for area distortions, giving a better estimate
    !        of the true ice area and volume.
    !       However, applying scale factors will give a mass conservation error (total dmass_dt > 0)
    !        in the diagnostics, because horizontal transport does not account for area factors.
    !        Transport conserves mass only under the assumption of rectangular grid cells.
    ! TODO - Tested only for Greenland (N. Hem.; projection origin offset from N. Pole). Test for other grids.

    if (associated(model%projection%stere)) then

       call glimmap_stere_area_factor(model%projection%stere,  &
                                      model%general%ewn,       &
                                      model%general%nsn,       &
                                      model%numerics%dew*len0, &
                                      model%numerics%dns*len0)
    endif

    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! If SMB input units are mm/yr w.e., then convert to units of acab.
    ! Note: In this case the input field should be called 'smb', not 'acab'.

    if (model%options%smb_input == SMB_INPUT_MMYR_WE) then

       ! make sure a nonzero SMB was read in
       var_maxval = maxval(abs(model%climate%smb))
       var_maxval = parallel_reduce_max(var_maxval)
       if (var_maxval < 1.0d-11) then
          write(message,*) 'Error: Failed to read in a nonzero SMB field with smb_input =', SMB_INPUT_MMYR_WE
          call write_log(trim(message), GM_FATAL)
       endif

       ! Convert units from mm/yr w.e. to m/yr ice
       model%climate%acab(:,:) = model%climate%smb(:,:) * (rhow/rhoi) / 1000.d0

       ! Convert acab from m/yr ice to model units
       model%climate%acab(:,:) = model%climate%acab(:,:) / scale_acab

    else
       ! assume acab was read in with units of m/yr ice; do nothing
    endif

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first

    call init_isostasy(model)

    select case(model%isostasy%whichrelaxed)

    case(RELAXED_TOPO_INPUT)   ! supplied input topography is relaxed

       model%isostasy%relx = model%geometry%topg

    case(RELAXED_TOPO_COMPUTE) ! supplied topography is in equilibrium
                               !TODO - Test the case RELAXED_TOPO_COMPUTE

       call isos_relaxed(model)

    end select

    ! open all output files
    call openall_out(model)

    ! create glide variables
    call glide_io_createall(model, model)

    ! Compute the cell areas of the grid
    model%geometry%cell_area = model%numerics%dew*model%numerics%dns

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

    ! Make sure the basal heat flux follows the positive-down sign convention
    if (maxval(model%temper%bheatflx) > 0.0d0) then
       write(message,*) 'Error, Input basal heat flux has positive values: '
       call write_log(trim(message))
       write(message,*) 'this_rank, maxval =', this_rank, maxval(model%temper%bheatflx)
       call write_log(trim(message))
       write(message,*) 'Basal heat flux is defined as positive down, so should be <= 0 on input'
       call write_log(trim(message), GM_FATAL)
    endif

    ! initialize model diagnostics
    call glide_init_diag(model)

!!    if (this_rank == model%numerics%rdiag_local) then
!!       i = model%numerics%idiag_local
!!       j = model%numerics%jdiag_local
!!       print*, ' '
!!       print*, 'this_rank, i, j, new acab:', this_rank, i, j, model%climate%acab(i,j) * scyr*thk0/tim0
!!    endif

    ! initialize glissade components

    ! Set some variables in halo cells
    ! Note: We need thck and artm in halo cells so that temperature will be initialized correctly
    !        (if not read from input file).
    !       We do an update here for temp in case temp is read from an input file.
    !       If temp is computed below in glissade_init_therm (based on the value of options%temp_init),
    !        then the halos will receive the correct values.
    call parallel_halo(model%geometry%thck)
    call parallel_halo(model%climate%artm)
    call parallel_halo(model%temper%temp)
    call parallel_halo(model%temper%tempunstag)

    ! calculate the lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    ! save starting ice thickness for diagnostics
    model%geometry%thck_old(:,:) = model%geometry%thck(:,:)

    ! Note: For outflow BCs, most fields (thck, usrf, temp, etc.) are set to zero in the global halo,
    !        to create ice-free conditions. However, we might not want to set topg = 0 in the global halo,
    !        because then the global halo will be interpreted as ice-free land, whereas we may prefer to
    !        treat it as ice-free ocean. For this reason, topg is extrapolated from adjacent cells.
    ! Note: For periodic BCs, there is an optional argument periodic_offset_ew for topg.
    !       This is for ismip-hom experiments. A positive EW offset means that
    !        the topography in west halo cells will be raised, and the topography
    !        in east halo cells will be lowered.  This ensures that the topography
    !        and upper surface elevation are continuous between halo cells
    !        and locally owned cells at the edge of the global domain.
    !       In other cases (anything but ismip-hom), periodic_offset_ew = periodic_offset_ns = 0,
    !        and this argument will have no effect.

    if (model%general%global_bc == GLOBAL_BC_OUTFLOW) then
       call parallel_halo_extrapolate(model%geometry%topg)
    else  ! other global BCs, including periodic
       call parallel_halo(model%geometry%topg, periodic_offset_ew = model%numerics%periodic_offset_ew)
    endif

    if (model%options%whichtemp == TEMP_ENTHALPY) call parallel_halo(model%temper%waterfrac)

    ! halo update for kinbcmask (= 1 where uvel and vvel are prescribed, elsewhere = 0)
    ! Note: Instead of assuming that kinbcmask is periodic, we extrapolate it into the global halo
    !       (and also into the north and east rows of the global domain, which are not included 
    !       on the global staggered grid).
    call staggered_parallel_halo_extrapolate (model%velocity%kinbcmask)  ! = 1 for Dirichlet BCs

    !TODO - Remove call to init_velo in glissade_initialise?
    !       Most of what's done in init_velo is needed for SIA only, but still need velowk for call to wvelintg
    call init_velo(model)

    ! Initialize basal hydrology, if needed
    call glissade_basal_water_init(model)

    ! Initialize the temperature profile in each column
    call glissade_init_therm(model%options%temp_init,    model%options%is_restart,  &
                             model%general%ewn,          model%general%nsn,         &
                             model%general%upn,                                     &
                             model%numerics%idiag_local, model%numerics%jdiag_local,&
                             model%numerics%rdiag_local,                            &
                             model%numerics%sigma,       model%numerics%stagsigma,  &
                             model%geometry%thck*thk0,                              & ! m
                             model%climate%artm,                                    & ! deg C
                             model%climate%acab*thk0/tim0,                          & ! m/s
                             model%temper%bheatflx,                                 & ! W/m^2, positive down
                             model%temper%pmp_offset,                               & ! deg C
                             model%temper%temp,                                     & ! deg C
                             model%temper%tempunstag)                                 ! deg C

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

    !TODO - Add halo updates of state variables here?

    ! If unstagbeta (i.e., beta on the scalar ice grid) was read from an input file,
    !  then interpolate it to beta on the staggered grid.
    ! NOTE: unstagbeta is initialized to unphys_val (a large negative number),
    !       so its maxval will be > 0 only if the field is read in.
    ! We can make an exception for ISHOM case C; for greater accuracy we set beta in 
    !  subroutine calcbeta instead of interpolating from unstagbeta (one processor only).

    if (maxval(model%velocity%unstagbeta) > 0.d0 .and.   &
               model%options%which_ho_babc /= HO_BABC_ISHOMC) then  ! interpolate to staggered grid
       call write_log('Interpolating beta from unstaggered (unstagbeta) to staggered grid (beta)')
       if (maxval(model%velocity%beta) > 0.0d0 ) then
          call write_log('Warning: the input "beta" field will be overwritten with values interpolated from the input &
               &"unstagbeta" field!')
       endif

       ! do a halo update and interpolate to the staggered grid
       ! Note: stagger_margin_in = 0 => use all values in the staggering, including where ice is absent
       call parallel_halo(model%velocity%unstagbeta)
       call glissade_stagger(model%general%ewn,          model%general%nsn,      &
                             model%velocity%unstagbeta,  model%velocity%beta,    &
                             stagger_margin_in = 0)

    endif  ! unstagbeta > 0

    ! Note: If the beta field has been read from an external file, then model%velocity%beta should not be modified.
    !       If beta is weighted by f_ground or otherwise modified, the modified field is model%velocity%beta_internal.

    ! The MISMIP 3D test case requires reading in a spatial factor that multiplies Coulomb_C.
    ! This factor is read in on the unstaggered grid, then interpolated to the staggered grid.
    ! If this factor is not present in the input file, then set it to 1 everywhere.

    if (maxval(model%basal_physics%C_space_factor) > tiny(0.d0)) then  ! C_space_factor was read in

       if (main_task) print*, 'stagger Coulomb spatial factor C'
       ! do a halo update and interpolate to the staggered grid
       ! Note: stagger_margin_in = 0 => use all values in the staggering, including where ice is absent
       call parallel_halo(model%basal_physics%C_space_factor)
       call glissade_stagger(model%general%ewn,                  model%general%nsn,                       &
                             model%basal_physics%C_space_factor, model%basal_physics%C_space_factor_stag, &
                             stagger_margin_in = 0)

    else  ! C_space_factor was not read in; set to 1 everywhere

       model%basal_physics%C_space_factor(:,:) = 1.d0
       model%basal_physics%C_space_factor_stag(:,:) = 1.d0

    endif

    ! Note: The basal process option is currently disabled.
    ! initialize basal process module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%dttem)
!!    end if      

    ! If acab is to be overwritten for some cells, then set overwrite_acab_mask = 1 for these cells.
    ! We can overwrite the input acab with a fixed value (typically negative) where
    ! (1) the input acab = 0 at initialization, or
    ! (2) the input thck <= overwrite_acab_minthck at initialization
    ! Note: This option is designed for standalone runs, and should be used only with caution for coupled runs.
    !       On restart, overwrite_acab_mask is read from the restart file.

    if (model%climate%overwrite_acab_value /= 0 .and. model%options%is_restart == RESTART_FALSE) then

!!       print*, 'Setting acab = overwrite value (m/yr):', model%climate%overwrite_acab_value * scyr*thk0/tim0

       call glissade_overwrite_acab_mask(model%options%overwrite_acab,          &
                                         model%climate%acab,                    &
                                         model%geometry%thck,                   &
                                         model%climate%overwrite_acab_minthck,  &
                                         model%climate%overwrite_acab_mask)

    endif

    ! calculate mask
    ! Note: This call includes a halo update for thkmask
    !TODO - Replace with glissade masks?
    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    ! compute halo for relaxed topography
    !TODO - Is this halo update necessary?
    call parallel_halo(model%isostasy%relx)

    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)

    ! optional unit tests

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

    ! Set debug diagnostics
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    ! initial calving, if desired
    ! Note: Do this only for a cold start with evolving ice, not for a restart
    if (l_evolve_ice .and. &
         model%options%calving_init == CALVING_INIT_ON .and. &
         model%options%is_restart == RESTART_FALSE) then

       ! ------------------------------------------------------------------------
       ! Note: The initial calving solve is treated differently from the runtime calving solve.
       !       In particular, calving-front culling is done at initialization only.
       !       Culling may also be used to remove a row of thin cells (~1 m)
       !        at the calving front, as present in some input data sets.
       !       But we do not want to remove calving_front cells every timestep.
       ! ------------------------------------------------------------------------

       call glissade_calving_solve(model, .true.)   ! init_calving = .true.

       ! The mask needs to be recalculated after calving.
       ! Note: glide_set_mask includes a halo update for thkmask.
       !TODO - Delete this call? Glissade dycore does not use thkmask.
       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)

    endif  ! initial calving

    ! Initialize the no-advance calving_mask, if desired
    ! Note: This is done after initial calving, which may include iceberg removal or calving-front culling.
    !       The calving front that exists after initial culling is the one that is held fixed during the simulation.
    ! Note: calving_front_x and calving_front_y already have units of m, so do not require multiplying by len0.

    if (model%options%whichcalving == CALVING_GRID_MASK .and. model%options%is_restart == RESTART_FALSE) then

       call glissade_calving_mask_init(&
                                       model%numerics%dew*len0,       model%numerics%dns*len0,        &
                                       model%geometry%thck*thk0,      model%geometry%topg*thk0,       &
                                       model%climate%eus*thk0,        model%numerics%thklim*thk0,     &
                                       model%calving%calving_front_x, model%calving%calving_front_y,  &
                                       model%calving%calving_mask)
    endif

    ! Note: The DIVA solver needs a halo update for effective viscosity.
    !       This is done at the end of glissade_diagnostic_variable_solve, which in most cases is sufficient.
    !       However, if we are (1) reading efvs from an input file and (2) solving for velocity before
    !        the halo update, then we need a halo update here too, to avoid symmetry issues in the velocity solver.
    !       I ran into this issue when running MISMIP+, which does cold starts (restart = 0) from files containing efvs.
    !       An update is done here regardless of code options, just to be on the safe side.
    call parallel_halo(model%stress%efvs)

    ! recalculate the lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0, model%geometry%thck + model%geometry%lsrf)

    !WHL - debug
    if (main_task) print*, 'Done in glissade_initialise'

  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model, time)

    ! Perform time-step of an ice model instance with the Glissade dycore

    use parallel

    use glimmer_paramets, only: tim0, len0, thk0
    use glimmer_physcon, only: scyr
    use glide_mask, only: glide_set_mask, calc_iareaf_iareag

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp), intent(in) :: time         ! current time in years

    ! --- Local variables ---

    integer :: i, j
    integer :: itest, jtest, rtest

    ! ========================

    ! Set debug diagnostics
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    ! Update internal clock
    model%numerics%time = time  
    model%numerics%tstep_count = model%numerics%tstep_count + 1

    if (main_task .and. verbose_glissade) print*, 'glissade_tstep, test_count =', model%numerics%tstep_count

    ! optional transport test
    ! code execution will end when this is done
    if (test_transport) then
       call glissade_test_transport (model)
       return
    endif

    ! save old ice thickness for diagnostics
    ! also used to reset thickness for the no-evolution option
    model%geometry%thck_old(:,:) = model%geometry%thck(:,:)

    ! ------------------------------------------------------------------------
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------
    !
    ! Note: This call used to be near the end of the glissade time step, between
    !       calving and the velocity solve. But this can be problematic, because
    !       a cell identified as grounded for calving purposes can become floating
    !       as a result of isostatic adjustment, or vice versa.
    !       It is better to compute isostasy just after the velocity solve,
    !       at the start of the next time step.
    !
    ! Matt Hoffman writes:
    ! Is this isostasy call in the right place?
    ! Consider for a forward Euler time step:
    ! With a relaxing mantle model, topg is a prognostic (time-evolving) variable:
    !      topg1 = f(topg0, thk0, ...)
    ! However, for a fluid mantle where the adjustment is instantaneous, topg is a diagnostic variable
    !(comparable to calculating floatation height of ice in the ocean):
    !      topg1 = f(thk1)
    ! In either case, the topg update should be separate from the thickness evolution (because thk1 = f(thk0, vel0=g(topg0,...)).
    ! However, if the isostasy calculation needs topg0, the icewaterload call should be made BEFORE thck is updated.
    ! If the isostasy calculation needs topg1, the icewaterload call should be made AFTER thck is updated.
    ! Also, we should think about when marinlim, usrf, lsrf, derivatives should be calculated relative to the topg update via isostasy.
    !
    ! WHL writes (May 2017):
    ! When isostasy is turned on, it is usually run with a relaxing mantle.
    ! With the call moved to the start of the time step, both the icewaterload call (if needed) and
    !  the relaxation are done before the ice thickness update. So we have
    !       topg1 = f(topg0, thk0, ...)
    !  followed by
    !       thk1  = f(thk0, vel0=g(topg0,...)
    ! I think this is what is desired.
    ! ------------------------------------------------------------------------

    call glissade_isostasy_solve(model)

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
    ! Do the vertical thermal solve if it is time to do so.
    ! Vertical diffusion and strain heating only; no temperature advection.
    ! Note: model%numerics%tinc and model%numerics%time have units of years.
    !       dttem has scaled units, so multiply by tim0/scyr to convert to years.
    ! ------------------------------------------------------------------------ 

    if ( model%numerics%tinc > mod(model%numerics%time, model%numerics%dttem*tim0/scyr)) then

       if (model%options%which_ho_thermal_timestep == HO_THERMAL_BEFORE_TRANSPORT) then

          ! vertical thermal solve before transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0)   ! convert dt from model units to s

       elseif (model%options%which_ho_thermal_timestep == HO_THERMAL_SPLIT_TIMESTEP) then

          ! vertical thermal solve split into two parts, before and after transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0/2.0d0)

       endif

    end if

    ! ------------------------------------------------------------------------ 
    ! Compute the basal melt rate beneath floating ice.
    ! (The basal melt rate beneath grounded ice is part of the thermal solve.)
    ! ------------------------------------------------------------------------ 

    call glissade_bmlt_float_solve(model)

    ! Add bmlt_float to bmlt_ground to determine the total bmlt

    ! Note: bmlt = bmlt_ground + bmlt_float may not be equal to the applied melt rate in a given cell,
    !       if ice is thin or absent in the cell.
    ! Note: bmlt_ground is computed in glissade_therm_driver.
    !       If glissade_therm_driver is called twice per time step, then bmlt_ground from
    !        the second time is not included in the transport solve, but is diagnostic only.
    !       That is, the transport scheme assumes that the bmlt_ground rate computed during the
    !        first call is applied during the entire time step.
    !       This might lead to small violations of energy conservation.
    !       TODO: Separate the bmlt_ground computation from the temperature computation?

    model%basal_melt%bmlt(:,:) = model%basal_melt%bmlt_ground(:,:) + model%basal_melt%bmlt_float(:,:)

    ! ------------------------------------------------------------------------ 
    ! Calculate ice thickness and tracer evolution under horizontal transport.
    ! The surface and basal mass balances are also applied here.
    ! ------------------------------------------------------------------------ 

    call glissade_thickness_tracer_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate iceberg calving
    ! ------------------------------------------------------------------------ 

    call glissade_calving_solve(model, .false.)   ! init_calving = .false.

    ! ------------------------------------------------------------------------
    ! Clean up variables in ice-free columns.
    ! This subroutine should be called after transport and calving, which may
    !  have set thck = 0 in some cells without zeroing out basal water and tracers.
    ! ------------------------------------------------------------------------ 

    call glissade_cleanup_icefree_cells(model)

    ! ------------------------------------------------------------------------ 
    ! Increment the ice age.
    ! If a cell becomes ice-free, the age is reset to zero.
    ! Note: Internally, the age has the same units as dt, but on output it will be converted to years.
    ! ------------------------------------------------------------------------ 
    
    if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
       do j = 1, model%general%nsn 
          do i = 1, model%general%ewn 
             if (model%geometry%thck(i,j) > 0.0d0) then
                model%geometry%ice_age(:,i,j) = model%geometry%ice_age(:,i,j) + model%numerics%dt
             else
                model%geometry%ice_age(:,i,j) = 0.0d0
             endif
          enddo
       enddo
    endif

    ! glissade_calve_ice adjusts thickness for calved ice.  Therefore the mask needs to be recalculated.
    ! Note: glide_set_mask includes a halo update of thkmask
    ! This time we want to calculate the optional arguments iarea and ivol because thickness 
    ! will not change further during this time step.
    !TODO - Remove this call to glide_set_mask?
    !       This subroutine is called at the beginning of glissade_velo_driver,
    !        so a call here is not needed for the velo diagnostic solve.
    !       The question is whether it is needed for the isostasy.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    ! --- Calculate global area of ice that is floating and grounded.
    !TODO  May want to calculate iareaf and iareag in glide_write_diag and remove those calculations here.  

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%thkmask,                        &
                            model%geometry%iareaf, model%geometry%iareag)

    ! ------------------------------------------------------------------------
    ! Do the vertical thermal solve if it is time to do so.
    ! Note: A thermal solve should be done here (using option HO_THERMAL_AFTER_TRANSPORT 
    !       or HO_THERMAL_SPLIT_TIMESTEP) if it is desired to update the bed temperature 
    !       and pmp temperature after transport and before the velocity solve.
    ! ------------------------------------------------------------------------

    if ( model%numerics%tinc > mod(model%numerics%time, model%numerics%dttem*tim0/scyr)) then

       if (model%options%which_ho_thermal_timestep == HO_THERMAL_AFTER_TRANSPORT) then

          ! vertical thermal solve after transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0)   ! convert dt from model units to s

       elseif (model%options%which_ho_thermal_timestep == HO_THERMAL_SPLIT_TIMESTEP) then

          ! vertical thermal solve split into two parts, before and after transport
          call glissade_thermal_solve(model,  &
                                      model%numerics%dttem*tim0/2.0d0)

       endif

    end if  ! take a temperature time step

    ! ------------------------------------------------------------------------
    ! Calculate diagnostic variables, including ice velocity
    ! ------------------------------------------------------------------------

    call glissade_diagnostic_variable_solve(model)

    !TODO - Any halo updates needed at the end of glissade_tstep?

  end subroutine glissade_tstep

!=======================================================================

  subroutine glissade_bmlt_float_solve(model)

    ! Solve for basal melting beneath floating ice.

    use glimmer_paramets, only: tim0, thk0, len0
    use glissade_bmlt_float, only: glissade_basal_melting_float
    use glissade_transport, only: glissade_add_mbal_anomaly
    use glissade_masks, only: glissade_get_masks

    use parallel

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

    integer, dimension(model%general%ewn, model%general%nsn) ::   &
       ice_mask,          & ! = 1 if ice is present (thck > 0, else = 0
       floating_mask        ! = 1 if ice is present (thck > 0) and floating

    real(dp) :: previous_time

    integer :: i, j

    integer :: itest, jtest, rtest

    ! Set debug diagnostics
    rtest = model%numerics%rdiag_local
    itest = model%numerics%idiag_local
    jtest = model%numerics%jdiag_local

    ! ------------------------------------------------------------------------
    ! Compute the basal melt rate beneath floating ice.
    ! Note: model%basal_melt is a derived type with various fields and parameters
    ! ------------------------------------------------------------------------

    !WHL - Put other simple options in this subroutine instead of glissade_basal_melting_float subroutine?

    if (main_task .and. verbose_glissade) print*, 'Call glissade_basal_melting_float'

    ! Compute masks:
    ! - ice_mask = 1 where thck > 0
    ! - floating_mask = 1 where thck > 0 and ice is floating;
    ! - ocean_mask = 1 where topg is below sea level and ice is absent
    !Note: The '0.0d0' argument is thklim. Here, any ice with thck > 0 gets ice_mask = 1.

    call glissade_get_masks(model%general%ewn,   model%general%nsn,     &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   0.0d0,                 &  ! thklim = 0
                            ice_mask,                                   &
                            floating_mask = floating_mask)

    ! Compute bmlt_float depending on the whichbmlt_float option

    if (model%options%whichbmlt_float == BMLT_FLOAT_NONE) then

       model%basal_melt%bmlt_float(:,:) = 0.0d0

    elseif (model%options%whichbmlt_float == BMLT_FLOAT_EXTERNAL) then

       ! Apply the external melt rate

       model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float_external(:,:)

       ! Optionally, multiply bmlt_float by a scalar adjustment factor
       if (model%basal_melt%bmlt_float_factor /= 1.0d0) then
          model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float(:,:) * model%basal_melt%bmlt_float_factor
       endif

    else  ! other options include BMLT_FLOAT_CONSTANT, BMLT_FLOAT_MISMIP, BMLT_FLOAT_DEPTH AND BMLT_FLOAT_MISOMIP
          !TODO - Call separate subroutines for each of these three options?

       !WHL - May want to comment out temporarily, if doing basal melting in the diagnostic solve for testing
       call glissade_basal_melting_float(model%options%whichbmlt_float,                                &
                                         model%general%ewn,          model%general%nsn,                &
                                         model%numerics%dew*len0,    model%numerics%dns*len0,          &
                                         model%numerics%idiag_local, model%numerics%jdiag_local,       &
                                         model%numerics%rdiag_local,                                   &
                                         model%general%x1,                                             & ! m
                                         model%geometry%thck*thk0,                                     & ! m
                                         model%geometry%lsrf*thk0,                                     & ! m
                                         model%geometry%topg*thk0,                                     & ! m
                                         model%climate%eus*thk0,                                       & ! m
                                         model%basal_melt)                                               ! bmlt_float in m/s

       ! Convert bmlt_float from SI units (m/s) to scaled model units
       model%basal_melt%bmlt_float(:,:) = model%basal_melt%bmlt_float(:,:) * tim0/thk0

    endif  ! whichbmlt_float

    ! If desired, add a bmlt_anomaly field.
    ! This is done for the initMIP Greenland and Antarctic experimennts.

    if (model%options%enable_bmlt_anomaly) then

       ! Compute the previous time
       ! Note: When being ramped up, the anomaly is not incremented until after the final time step of the year.
       !       This is the reason for passing the previous time to the subroutine.
       previous_time = model%numerics%time - model%numerics%dt * tim0/scyr

       ! Add the bmlt_float anomaly where ice is present and floating
       call glissade_add_mbal_anomaly(model%basal_melt%bmlt_float,              &   ! scaled model units
                                      model%basal_melt%bmlt_float_anomaly,      &   ! scaled model units 
                                      model%basal_melt%bmlt_anomaly_timescale,  &   ! yr
                                      previous_time)                                ! yr

    endif

    ! Limit the melting to cells where ice is present and floating.
    ! Note: For this to work correctly, basal melting must be applied before the floating mask changes
    !       (i.e., before horizontal transport).

    where (floating_mask == 0)
       model%basal_melt%bmlt_float = 0.0d0
    endwhere

    if (model%options%which_ho_ground_bmlt == HO_GROUND_BMLT_GLP) then

       ! Set bmlt_float = 0 in any grid cells that are grounded, or are adjacent to grounded cells.
       ! This prevents grounding-line retreat driven by spurious thinning of grounded ice.

       do j = 2, model%general%nsn-1
          do i = 2, model%general%ewn-1
             if (floating_mask(i,j) == 1) then
                ! check for grounded edge neighbors
                ! Note: There might be grounded corner neighbors, but we assume the GL passes
                !       through the cell only if it has a grounded edge neighbor.
                if ( (ice_mask(i-1,j)==1 .and. floating_mask(i-1,j)==0) .or. &
                     (ice_mask(i+1,j)==1 .and. floating_mask(i+1,j)==0) .or. &
                     (ice_mask(i,j-1)==1 .and. floating_mask(i,j-1)==0) .or. &
                     (ice_mask(i,j+1)==1 .and. floating_mask(i,j+1)==0) ) then
                   model%basal_melt%bmlt_float(i,j) = 0.0d0
                endif
             endif
          enddo
       enddo

    endif

  end subroutine glissade_bmlt_float_solve

!=======================================================================

  subroutine glissade_thermal_solve(model, dt)

    ! Do the vertical thermal solve.
    ! First call a driver subroutine for vertical temperature or enthalpy evolution,
    ! and then update the basal water.
    use parallel

    use glimmer_paramets, only: tim0, thk0, len0
    use glissade_therm, only: glissade_therm_driver
    use glissade_basal_water, only: glissade_calcbwat

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    real(dp), intent(in) :: dt   ! time step (s)
    
    ! unscaled model parameters (SI units)
    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
       bmlt_ground_unscaled,   & ! basal melt rate for grounded ice (m/s)
       bwat_unscaled             ! basal water thickness (m)

    integer :: up

    call t_startf('glissade_thermal_solve')

    if (main_task .and. verbose_glissade) print*, 'Call glissade_therm_driver'

    ! Note: glissade_therm_driver uses SI units
    !       Output arguments are temp, waterfrac, bpmp and bmlt_ground
    call glissade_therm_driver (model%options%whichtemp,                                      &
                                model%options%temp_init,                                      &
                                dt,                                                           & ! s
                                model%general%ewn,          model%general%nsn,                &
                                model%general%upn,                                            &
                                model%numerics%idiag_local, model%numerics%jdiag_local,       &
                                model%numerics%rdiag_local,                                   &
                                model%numerics%sigma,       model%numerics%stagsigma,         &
                                model%numerics%thklim_temp*thk0,                              & ! m
                                model%geometry%thck*thk0,                                     & ! m
                                model%geometry%topg*thk0,                                     & ! m
                                model%geometry%lsrf*thk0,                                     & ! m
                                model%climate%eus*thk0,                                       & ! m
                                model%climate%artm,                                           & ! deg C    
                                model%temper%bheatflx,      model%temper%bfricflx,            & ! W/m2
                                model%temper%dissip,                                          & ! deg/s
                                model%temper%pmp_threshold,                                   & ! deg C
                                model%temper%bwat*thk0,                                       & ! m
                                model%temper%temp,                                            & ! deg C
                                model%temper%waterfrac,                                       & ! unitless
                                model%temper%bpmp,                                            & ! deg C
                                bmlt_ground_unscaled)                                           ! m/s
                                     
    ! Update basal hydrology, if needed
    ! Note: glissade_calcbwat uses SI units

    if (main_task .and. verbose_glissade) print*, 'Call glissade_calcbwat'

    ! convert bwat to SI units for input to glissade_calcbwat
    bwat_unscaled(:,:) = model%temper%bwat(:,:) * thk0

    call glissade_calcbwat(model%options%which_ho_bwat,      &
                           model%basal_physics,              &
                           dt,                               &  ! s
                           model%geometry%thck*thk0,         &  ! m
                           model%numerics%thklim_temp*thk0,  &  ! m
                           bmlt_ground_unscaled,             &  ! m/s
                           bwat_unscaled)                       ! m

    ! convert bmlt and bwat from SI units (m/s and m) to scaled model units
    model%basal_melt%bmlt_ground(:,:) = bmlt_ground_unscaled(:,:) * tim0/thk0
    model%temper%bwat(:,:) = bwat_unscaled(:,:) / thk0

    ! Update tempunstag as sigma weighted interpolation from temp to layer interfaces
    do up = 2, model%general%upn-1
      model%temper%tempunstag(up,:,:) = model%temper%temp(up-1,:,:) +       &
           (model%temper%temp(up,:,:) - model%temper%temp(up-1,:,:)) *      &
           (model%numerics%sigma(up) - model%numerics%stagsigma(up-1)) /    &
           (model%numerics%stagsigma(up) - model%numerics%stagsigma(up-1))
    end do
    ! boundary conditions are identical on both grids, but temp starts at index 0
    model%temper%tempunstag(1,:,:) = model%temper%temp(0,:,:)
    model%temper%tempunstag(model%general%upn,:,:) = model%temper%temp(model%general%upn,:,:)

    !------------------------------------------------------------------------ 
    ! Halo updates
    !------------------------------------------------------------------------ 
    
    ! Note: bwat is needed in halos to compute effective pressure if which_ho_effecpress = HO_EFFECPRESS_BWAT
    call parallel_halo(model%temper%bwat)

    call t_stopf('glissade_thermal_solve')
    
  end subroutine glissade_thermal_solve

!=======================================================================

  subroutine glissade_thickness_tracer_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate ice thickness and tracer evolution, including horizontal transport and surface and basal mass balance.
    ! MJH: This subroutine uses velocity from the previous time step, which is appropriate for a Forward Euler time-stepping scheme.
    ! WHL: We used to have EVOL_NO_THICKNESS = -1 as a Glide option, used to hold the ice surface elevation fixed during CESM runs. 
    !      This option has been replaced by a Glint/Glad option, evolve_ice.
    !      We now have EVOL_NO_THICKESS = 5 as a glissade option.  It is used to hold the ice surface elevation fixed
    !       while allowing temperature to evolve, which can be useful for model spinup.  This option might need more testing.
    ! ------------------------------------------------------------------------ 

    use parallel

    use glimmer_paramets, only: tim0, thk0, vel0, len0
    use glimmer_physcon, only: scyr
    use glimmer_scales, only: scale_acab
    use glissade_therm, only: glissade_temp2enth, glissade_enth2temp
    use glissade_transport, only: glissade_mass_balance_driver, &
                                  glissade_transport_driver, &
                                  glissade_check_cfl,  &
                                  glissade_transport_setup_tracers, &
                                  glissade_transport_finish_tracers, &
                                  glissade_overwrite_acab,  &
                                  glissade_add_mbal_anomaly
    use glissade_masks, only: glissade_get_masks

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    integer :: sc  ! subcycling index

    ! temporary arrays in SI units
    real(dp), dimension(model%general%ewn,model%general%nsn) ::   &
       thck_unscaled,     & ! ice thickness (m)
       topg_unscaled,     & ! bedrock topography (m)
       thck_new_unscaled, & ! expected new ice thickness, after mass balance (m)
       acab_unscaled,     & ! surface mass balance (m/s)
       bmlt_unscaled        ! = bmlt (m/s) if basal mass balance is included in continuity equation, else = 0

    ! masks
    integer, dimension(model%general%ewn, model%general%nsn) ::   &
       ice_mask,             & ! = 1 if thck > 0, else = 0
       floating_mask,        & ! = 1 where ice is present and floating, else = 0
       ocean_mask,           & ! = 1 if topg is below sea level and thck = 0, else = 0
       land_mask,            & ! = 1 if topg is at or above sea level, else = 0
       grounding_line_mask,  & ! = 1 if a cell is adjacent to the grounding line, else = 0
       calving_front_mask      ! = 1 where ice is floating and borders an ocean cell, else = 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
       thck_flotation,       & ! thickness at which ice is exactly floating
       thck_calving_front,   & ! effective thickness of ice at the calving front
       effective_areafrac      ! effective fractional area of ice at the calving front

    real(dp) :: previous_time       ! time (yr) at the start of this time step
                                    ! (The input time is the time at the end of the step.)

    real(dp) :: advective_cfl       ! advective CFL number
                                    ! If advective_cfl > 1, the model is unstable without subcycling
    real(dp) :: dt_transport        ! time step (s) for transport; = model%numerics%dt*tim0 by default

    integer :: nsubcyc              ! number of times to subcycle advection

    logical :: do_upwind_transport  ! logical for whether transport code should do upwind transport or incremental remapping
                                    ! set to true for EVOL_UPWIND, else = false

    integer :: ntracers             ! number of tracers to be transported

    integer :: i, j, k
    integer :: ewn, nsn, upn
    integer :: itest, jtest, rtest

    rtest = -999
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

    select case(model%options%whichevol)

    case(EVOL_INC_REMAP, EVOL_UPWIND, EVOL_NO_THICKNESS) 

       if (model%options%whichevol == EVOL_UPWIND) then
          do_upwind_transport = .true.
       else
          do_upwind_transport = .false.
       endif

       !-------------------------------------------------------------------------
       ! First transport ice thickness and temperature, given the horizontal velocity (u,v).
       ! Then apply the surface and basal mass balance in each grid cell.
       ! Note: The main reason to do both transport and mass balance in one subroutine is
       !  that both operations require tracer array setup and cleanup (e.g., copying the 
       !  various tracer fields into generic tracer arrays). With this arrangement,
       !  the tracer operations need to be done only once.
       !-------------------------------------------------------------------------
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       !      (but not required) to use IR to advect temperature when thickness evolution is turned off.
       ! TODO  MJH If we really want to support no evolution, then we may want to implement it so that IR does not occur 
       !       at all - right now a run can fail because of a CFL violation in IR even if evolution is turned off.  Do we want
       !       to support temperature evolution without thickness evolution?  If so, then the current implementation may be preferred approach.

       call t_startf('inc_remap_driver')

       if (verbose_glissade .and. main_task) print *, 'Compute dH/dt'

       call t_startf('glissade_transport_driver')

       ! ------------------------------------------------------------------------
       ! Compute some masks prior to horizontal transport.
       ! ------------------------------------------------------------------------

       call glissade_get_masks(ewn,                      nsn,              &
                               model%geometry%thck*thk0,                   &   ! m
                               model%geometry%topg*thk0,                   &   ! m
                               model%climate%eus*thk0,                     &   ! m
                               0.0d0,                                      &   ! thklim = 0
                               ice_mask,                                   &
                               floating_mask = floating_mask,              &
                               ocean_mask = ocean_mask,                    &
                               land_mask = land_mask,                      &
                               grounding_line_mask = grounding_line_mask)

       ! For the enthalpy option, derive enthalpy from temperature and waterfrac.
       ! Must transport enthalpy rather than temperature/waterfrac to conserve energy.
       !TODO: Move this code to glissade_transport_setup_tracers?

       if (model%options%whichtemp == TEMP_ENTHALPY) then  ! Use IR to transport enthalpy
          ! Note: glissade_temp2enth expects SI units
          do j = 1, nsn
             do i = 1, ewn
                call glissade_temp2enth (model%numerics%stagsigma(1:upn-1),        &
                                         model%temper%temp(0:upn,i,j),     model%temper%waterfrac(1:upn-1,i,j),   &
                                         model%geometry%thck(i,j)*thk0,    model%temper%enthalpy(0:upn,i,j))
             enddo
          enddo
       endif    ! TEMP_ENTHALPY

       ! copy tracers (temp/enthalpy, etc.) into model%geometry%tracers
       call glissade_transport_setup_tracers (model)

       ! temporary in/out arrays in SI units
       thck_unscaled(:,:) = model%geometry%thck(:,:) * thk0
       topg_unscaled(:,:) = model%geometry%topg(:,:) * thk0

       ! pre-transport halo updates for thickness and tracers
       call parallel_halo(thck_unscaled)
       call parallel_halo(topg_unscaled)
       call parallel_halo_tracers(model%geometry%tracers)
       call parallel_halo_tracers(model%geometry%tracers_usrf)
       call parallel_halo_tracers(model%geometry%tracers_lsrf)

       ! pre-transport halo updates for velocity
       ! Velocity update might be needed if velo has not been updated in the halo since the previous diagnostic solve.
       !  (just to be on the safe side).

       call staggered_parallel_halo(model%velocity%uvel)
       call staggered_parallel_halo(model%velocity%vvel)

       ! --- Determine CFL limits ---
       ! Note: We are using the subcycled dt here (if subcycling is on).
       !  (See note above about the EVOL_NO_THICKNESS option and how it is affected by a CFL violation)
       !  stagthck, dusrfdew/ns and u/vvel need to be from the previous time step (and are at this point)
       ! Note: If using adaptive subcycling (with adaptive_cfl_threshold > 0), then dt_transport should
       !       be equal to dt (which is the case by default).
       !TODO - Remove the dt_transport option and simply rely on adaptive subcycling as needed?

       call glissade_check_cfl(ewn,                       nsn,                       upn-1,                    &
                               model%numerics%dew * len0, model%numerics%dns * len0, model%numerics%sigma,     &
                               model%geomderv%stagthck * thk0,                                                 &
                               model%geomderv%dusrfdew*thk0/len0, model%geomderv%dusrfdns*thk0/len0,           &
                               model%velocity%uvel * scyr * vel0, model%velocity%vvel * scyr * vel0,           &
                               model%numerics%dt_transport * tim0 / scyr,                                      &
                               model%numerics%adv_cfl_dt,         model%numerics%diff_cfl_dt)

       ! Set the transport timestep.
       ! The timestep is model%numerics%dt by default, but optionally can be reduced for subcycling

       if (model%numerics%adaptive_cfl_threshold > 0.0d0) then

          !WHL - debug
!          if (main_task) then
!             print*, 'Check advective CFL threshold'
!             print*, 'model dt (yr) =', model%numerics%dt * tim0/scyr
!             print*, 'adv_cfl_dt    =', model%numerics%adv_cfl_dt
!          endif

          advective_cfl = model%numerics%dt*(tim0/scyr) / model%numerics%adv_cfl_dt
          if (advective_cfl > model%numerics%adaptive_cfl_threshold) then

             ! compute the number of subcycles
             ! If advective_cfl > advective_cfl_threshold, then nsubcyc >= 2
             ! The larger the ratio, the larger the value of nsubcyc
             nsubcyc = ceiling(advective_cfl / model%numerics%adaptive_cfl_threshold)

             if (main_task) then
                print*, 'WARNING: adv_cfl_dt exceeds threshold; CFL =', advective_cfl
                print*, 'Ratio =', advective_cfl / model%numerics%adaptive_cfl_threshold
                print*, 'nsubcyc =', nsubcyc
             endif
          else
             nsubcyc = 1
          endif
          dt_transport = model%numerics%dt * tim0 / real(nsubcyc,dp)   ! convert to s

       else  ! no adaptive subcycling
          nsubcyc = model%numerics%subcyc
          dt_transport = model%numerics%dt_transport * tim0  ! convert to s
       endif

       !-------------------------------------------------------------------------
       ! Compute horizontal transport of mass and tracers, subcycling as needed.
       !-------------------------------------------------------------------------

       do sc = 1, nsubcyc

          if (nsubcyc > 1 .and. main_task) write(*,*) 'Subcycling transport: Cycle', sc

          ! Call the transport driver subroutine.
          ! (includes a halo update for thickness: thck_unscaled in this case)
          !
          ! Note: This subroutine assumes SI units:
          !       * dt (s)
          !       * dew, dns, thck (m)
          !       * uvel, vvel (m/s)
          !       Since thck has intent(inout), we create and pass a temporary array (thck_unscaled) with units of m.
          ! Note: tracers_ursf and tracers_lsrf are not transported, but they provide upper and lower BCs
          !       for vertical remapping. They are intent(in).

          call glissade_transport_driver(dt_transport,                                         &  ! s
                                         model%numerics%dew * len0, model%numerics%dns * len0, &
                                         ewn,          nsn,         upn-1,                     &
                                         model%numerics%sigma,                                 &
                                         model%velocity%uvel(:,:,:) * vel0,                    &  ! m/s
                                         model%velocity%vvel(:,:,:) * vel0,                    &  ! m/s
                                         thck_unscaled(:,:),                                   &  ! m
                                         model%geometry%ntracers,                              &
                                         model%geometry%tracers(:,:,:,:),                      &
                                         model%geometry%tracers_usrf(:,:,:),                   &
                                         model%geometry%tracers_lsrf(:,:,:),                   &
                                         model%options%which_ho_vertical_remap,                &
                                         upwind_transport_in = do_upwind_transport)

          ! halo updates for thickness and tracers
          call parallel_halo(thck_unscaled)
          call parallel_halo_tracers(model%geometry%tracers)

       enddo     ! subcycling of transport

       !-------------------------------------------------------------------------
       ! Prepare the surface and basal mass balance terms.
       ! Note: The basal mass balance has been computed in subroutine glissade_bmlt_float_solve.
       !-------------------------------------------------------------------------

       ! Compute a corrected acab field that includes any prescribed anomalies.
       ! Typically, acab_corrected = acab, but sometimes (e.g., for initMIP) it includes a time-dependent anomaly.
       ! Note that acab itself does not change in time.

       ! initialize
       model%climate%acab_corrected(:,:) = model%climate%acab(:,:)

       ! Optionally, multiply acab by a scalar adjustment factor
       if (model%climate%acab_factor /= 1.0d0) then
          model%climate%acab_corrected(:,:) = model%climate%acab_corrected(:,:) * model%climate%acab_factor
       endif

       if (model%options%enable_acab_anomaly) then

          ! Note: When being ramped up, the anomaly is not incremented until after the final time step of the year.
          !       This is the reason for passing the previous time to the subroutine.
          previous_time = model%numerics%time - model%numerics%dt * tim0/scyr

          call glissade_add_mbal_anomaly(model%climate%acab_corrected,          &   ! scaled model units
                                         model%climate%acab_anomaly,            &   ! scaled model units
                                         model%climate%acab_anomaly_timescale,  &   ! yr
                                         previous_time)

          !WHL - debug
!!          if (this_rank==rtest) then
!!             i = model%numerics%idiag
!!             j = model%numerics%jdiag
!!             print*, 'i, j, total anomaly (m/yr), previous_time, new acab (m/yr):', &
!!                      i, j, model%climate%acab_anomaly(i,j)*thk0*scyr/tim0, previous_time, model%climate%acab_corrected(i,j)
!!          endif

       endif

       ! Optionally, overwrite acab_corrected where overwrite_acab_mask = 1.

       if (model%options%overwrite_acab /= 0) then

          call glissade_overwrite_acab(model%climate%overwrite_acab_mask,  &
                                       model%climate%overwrite_acab_value, &
                                       model%climate%acab_corrected)
       endif

       ! Convert acab_corrected to a temporary array in SI units (m/s)
       acab_unscaled(:,:) = model%climate%acab_corrected(:,:) * thk0/tim0


       ! Convert bmlt in SI units (m/s)
       ! Note: bmlt is the sum of bmlt_ground (computed in glissade_thermal_solve) and bmlt_float
       !       (computed in glissade_bmlt_float_solve).
       ! Note: bmlt can be turned off by setting options%basal_mbal = BASAL_MBAL_NO_CONTINUITY

       if (model%options%basal_mbal == BASAL_MBAL_CONTINUITY) then    ! include bmlt in continuity equation
          bmlt_unscaled(:,:) = model%basal_melt%bmlt(:,:) * thk0/tim0
       else                                                           ! do not include bmlt in continuity equation
          bmlt_unscaled(:,:) = 0.0d0
       endif

       ! ------------------------------------------------------------------------
       ! Get masks used for the mass balance calculation.
       ! Pass thklim = 0 to identify cells with thck > 0 (not thck > thklim).
       ! Use ocean_mask to identify ocean cells where positive acab should not be applied.
       ! Use thck_calving_front to compute a fractional area for calving_front cells.
       ! ------------------------------------------------------------------------

       call glissade_get_masks(ewn,                      nsn,              &
                               thck_unscaled,                              &   ! m
                               topg_unscaled,                              &   ! m
                               model%climate%eus*thk0,                     &   ! m
                               0.0d0,                                      &   ! thklim = 0
                               ice_mask,                                   &
                               floating_mask = floating_mask,              &
                               land_mask = land_mask,                      &
                               ocean_mask = ocean_mask,                    &
                               which_ho_calving_front = model%options%which_ho_calving_front, &
                               calving_front_mask = calving_front_mask,    &
                               thck_calving_front = thck_calving_front)

       ! Compute the effective fractional area of calving_front cells.

       do j = 1, nsn
          do i = 1, ewn
             if (calving_front_mask(i,j) == 1 .and. thck_calving_front(i,j) > 0.0d0) then
                effective_areafrac(i,j) = thck_unscaled(i,j) / thck_calving_front(i,j)
                effective_areafrac(i,j) = min(effective_areafrac(i,j), 1.0d0)
             elseif (ocean_mask(i,j) == 1) then
                effective_areafrac(i,j) = 0.0d0  ! acab and bmlt not applied to ice-free ocean cells
             else  ! non-CF ice-covered cells and/or land cells
                effective_areafrac(i,j) = 1.0d0
             endif
          enddo
       enddo

       ! Initialize the applied acab and bmlt.
       ! Note: These are smaller in magnitude than the input acab and bmlt for cells where either
       !       (1) the full column melts, and energy remains for melting, or
       !       (2) a positive mass balance is ignored, because a cell is ice-free ocean

       model%climate%acab_applied = 0.0d0
       model%basal_melt%bmlt_applied = 0.0d0

       ! ------------------------------------------------------------------------
       ! Apply the surface mass balance (acab) and basal mass balance (bmlt).
       ! Note: This subroutine assumes SI units:
       !       * dt (s)
       !       * dew, dns, thck (m)
       !       * acab, bmlt (m/s)
       ! ------------------------------------------------------------------------

       call glissade_mass_balance_driver(model%numerics%dt * tim0,                             &
                                         model%numerics%dew * len0, model%numerics%dns * len0, &
                                         ewn,         nsn,          upn-1,                     &
                                         model%numerics%sigma,                                 &
                                         thck_unscaled(:,:),                                   &  ! m
                                         acab_unscaled(:,:),                                   &  ! m/s
                                         bmlt_unscaled(:,:),                                   &  ! m/s
                                         model%climate%acab_applied(:,:),                      &  ! m/s
                                         model%basal_melt%bmlt_applied(:,:),                   &  ! m/s
                                         ocean_mask(:,:),                                      &
                                         effective_areafrac(:,:),                              &
                                         model%geometry%ntracers,                              &
                                         model%geometry%tracers(:,:,:,:),                      &
                                         model%geometry%tracers_usrf(:,:,:),                   &
                                         model%geometry%tracers_lsrf(:,:,:),                   &
                                         model%options%which_ho_vertical_remap)

       !WHL - debug
       call parallel_halo(thck_unscaled)

       !-------------------------------------------------------------------------
       ! Cleanup
       !-------------------------------------------------------------------------

       ! copy tracers (temp/enthalpy, etc.) from model%geometry%tracers back to standard arrays
       call glissade_transport_finish_tracers(model)

       ! convert applied mass balance from m/s back to scaled model units
       model%climate%acab_applied(:,:) = model%climate%acab_applied(:,:)/thk0 * tim0
       model%basal_melt%bmlt_applied(:,:) = model%basal_melt%bmlt_applied(:,:)/thk0 * tim0

       ! convert thck back to scaled units
       ! (acab_unscaled is intent(in) above, so no need to scale it back)
       model%geometry%thck(:,:) = thck_unscaled(:,:) / thk0

       ! For the enthalpy option, convert enthalpy back to temperature/waterfrac.

       if (model%options%whichtemp == TEMP_ENTHALPY) then

          ! Note: glissade_enth2temp expects SI units
          do j = 1, nsn
             do i = 1, ewn
                call glissade_enth2temp(model%numerics%stagsigma(1:upn-1),                                    &
                                        model%geometry%thck(i,j)*thk0,    model%temper%enthalpy(0:upn,i,j),   &
                                        model%temper%temp(0:upn,i,j),     model%temper%waterfrac(1:upn-1,i,j))
             enddo
          enddo
         
       endif    ! TEMP_ENTHALPY

       if (this_rank==rtest .and. verbose_glissade) then
          print*, ' '
          print*, 'After glissade_transport_driver:'
          print*, 'max, min thck (m)=', maxval(model%geometry%thck)*thk0, minval(model%geometry%thck)*thk0
          print*, 'max, min acab (m/yr) =', &
                  maxval(model%climate%acab_corrected)*scale_acab, &
                  minval(model%climate%acab_corrected)*scale_acab
          print*, 'max, min artm =', maxval(model%climate%artm), minval(model%climate%artm)
          print*, 'thklim =', model%numerics%thklim * thk0
          print*, 'max, min temp =', maxval(model%temper%temp), minval(model%temper%temp)
          print*, ' '
          print*, 'thck:'
          write(6,'(a6)',advance='no') '      '
          do i = itest-5, itest+5
             write(6,'(i14)',advance='no') i
          enddo
          write(6,*) ' '
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-5, itest+5
                write(6,'(f14.7)',advance='no') model%geometry%thck(i,j) * thk0
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          k = upn
          print*, 'temp, k =', k
          write(6,'(a6)',advance='no') '      '
          do i = itest-5, itest+5
             write(6,'(i14)',advance='no') i
          enddo
          write(6,*) ' '
          do j = jtest+2, jtest-2, -1
             write(6,'(i6)',advance='no') j
             do i = itest-5, itest+5
                write(6,'(f14.7)',advance='no') model%temper%temp(k,i,j)
             enddo
             write(6,*) ' '
          enddo
       endif

       call t_stopf('glissade_transport_driver')

       call t_stopf('inc_remap_driver')

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck(:,:) = model%geometry%thck_old(:,:)
       endif

    end select

    !------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    !TODO - Not sure this update is needed here.  It is done at the start
    !       of the diagnostic solve, but may not be needed for calving.
    !------------------------------------------------------------------------
    
    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       &
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

  end subroutine glissade_thickness_tracer_solve

!=======================================================================

  subroutine glissade_calving_solve(model, init_calving)

    ! ------------------------------------------------------------------------ 
    ! Calculate iceberg calving
    ! ------------------------------------------------------------------------ 

    use parallel

    use glimmer_paramets, only: thk0, tim0, len0
    use glissade_calving, only: glissade_calve_ice
    use glide_mask, only: glide_set_mask

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    logical, intent(in) :: init_calving  ! true when this subroutine is called at initialization

    ! --- Local variables ---

    real(dp), dimension(model%general%ewn, model%general%nsn) :: &
         thck_unscaled              ! model%geometry%thck converted to m

    logical :: cull_calving_front   ! true iff init_calving = T and options%cull_calving_front = T

    integer :: i, j

    !TODO - Make sure no additional halo updates are needed before glissade_calve_ice

    ! Determine whether the calving front should be culled
    if (init_calving .and. model%options%cull_calving_front) then
       cull_calving_front = .true.
    else
       cull_calving_front = .false.
    endif

    ! ------------------------------------------------------------------------ 
    ! Calve ice, based on the value of whichcalving 
    !       Pass in thck, topg, etc. with units of meters.
    ! ------------------------------------------------------------------------ 

    thck_unscaled(:,:) = model%geometry%thck(:,:)*thk0

    call glissade_calve_ice(model%options%whichcalving,        &
                            model%options%calving_domain,      &
                            model%options%which_ho_calving_front, &
                            model%options%remove_icebergs,     &
                            model%options%limit_marine_cliffs, &
                            cull_calving_front,                &
                            model%calving,                     &        ! calving object; includes calving_thck (m)
                            model%numerics%idiag_local,        &
                            model%numerics%jdiag_local,        &
                            model%numerics%rdiag_local,        &
                            model%numerics%dt*tim0,            &        ! s
                            model%numerics%dew*len0,           &        ! m
                            model%numerics%dns*len0,           &        ! m
                            model%numerics%sigma,              &
                            model%numerics%thklim*thk0,        &        ! m
                            thck_unscaled,                     &        ! m
                            model%isostasy%relx*thk0,          &        ! m
                            model%geometry%topg*thk0,          &        ! m
                            model%climate%eus*thk0)                     ! m
    
    ! Convert geometry%thck and calving%calving_thck to scaled model units
    model%geometry%thck(:,:) = thck_unscaled(:,:)/thk0
    model%calving%calving_thck(:,:) = model%calving%calving_thck(:,:)/thk0

    !TODO: Are any other halo updates needed after calving?
    ! halo updates
    call parallel_halo(model%geometry%thck)    ! Updated halo values of thck are needed below in calclsrf

    ! Eliminate ice from cells where a no-advance mask prohibits it.
    ! Add this ice to the calving field for mass conservation diagnostics.
    ! Note: The calving_mask option accomplishes the same thing. However, if we want to use a different
    !       calving law while also limiting calving-front advance, we can use the no_advance mask.
    do j = 1, model%general%nsn
       do i = 1, model%general%ewn
          if (model%climate%no_advance_mask(i,j) == 1 .and. model%geometry%thck(i,j) > 0.0d0) then
             model%calving%calving_thck(i,j) = model%calving%calving_thck(i,j) + model%geometry%thck(i,j)
             model%geometry%thck(i,j) = 0.d0
             ! Note: Tracers (temp, etc.) are cleaned up later by call to glissade_cleanup_icefree_cells
          endif
       enddo
    enddo

    ! update the upper and lower surfaces

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       &
                        model%climate%eus,   model%geometry%lsrf)
    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

  end subroutine glissade_calving_solve

!=======================================================================

  subroutine glissade_isostasy_solve(model)

    ! ------------------------------------------------------------------------ 
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------ 

    use parallel
    use isostasy

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! --- Local variables ---

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! Note: Suppose the update period is 100 years, and the time step is 1 year.
    !       Then the update will be done on the first time step of the simulation,
    !        (model%numerics%tstep_count = 1) and again on step 101, 201, etc.
    !       The update will not be done before writing output at t = 100, when
    !        model%numerics%tstep_count = 100.
    !       Thus the output file will contain the load that was applied during the
    !        preceding years, not the new load.
    !       In older code versions, the new load would have been computed on step 100.
    ! ------------------------------------------------------------------------

    if (model%options%isostasy == ISOSTASY_COMPUTE) then

       if (model%isostasy%nlith > 0) then
          if (mod(model%numerics%tstep_count-1, model%isostasy%nlith) == 0) then
             if (main_task) then
                print*, 'Update lithospheric load: tstep_count, nlith =', &
                     model%numerics%tstep_count, model%isostasy%nlith
             endif
             call isos_icewaterload(model)
             model%isostasy%new_load = .true.
          end if
       endif  ! nlith > 0

    end if
   
    ! ------------------------------------------------------------------------ 
    ! Calculate isostatic adjustment
    ! ------------------------------------------------------------------------ 

    if (model%options%isostasy == ISOSTASY_COMPUTE) then
       call isos_compute(model)

       ! update topography in halo cells
       ! Note: For outflow BCs, most fields (thck, usrf, temp, etc.) are set to zero in the global halo,
       !        to create ice-free conditions. However, we might not want to set topg = 0 in the global halo,
       !        because then the global halo will be interpreted as ice-free land, whereas we may prefer to
       !        treat it as ice-free ocean. For this reason, topg is extrapolated from adjacent cells.
       ! Note: The topg halo update at initialization has an optional argument periodic_ew,
       !        which is needed for ismip-hom. I doubt ismip-hom will be run with active isostasy,
       !        but the argument is included to be on the safe side.

       if (model%general%global_bc == GLOBAL_BC_OUTFLOW) then
          call parallel_halo_extrapolate(model%geometry%topg)
       else  ! other global BCs, including periodic
          call parallel_halo(model%geometry%topg, periodic_offset_ew = model%numerics%periodic_offset_ew)
       endif
    end if

  end subroutine glissade_isostasy_solve

!=======================================================================

  subroutine glissade_diagnostic_variable_solve(model) 

     ! Solve diagnostic (not time-dependent) variables, in particular the ice velocity.
     ! This is needed at the end of each time step once the prognostic variables (thickness, tracers) have been updated.  
     ! It is also needed to fill out the initial state from the fields that have been read in.

    use parallel

    use glimmer_paramets, only: tim0, len0, vel0, thk0, vis0, tau0, evs0
    use glimmer_physcon, only: scyr
    use glimmer_scales, only: scale_acab
    use glide_thck, only: glide_calclsrf
    use glam_velo, only: glam_velo_driver, glam_basal_friction
    use glissade_velo, only: glissade_velo_driver
    use glide_velo, only: wvelintg
    use glissade_masks, only: glissade_get_masks
    use glissade_grounding_line, only: glissade_grounded_fraction, glissade_grounding_line_flux
    use glissade_therm, only: glissade_interior_dissipation_sia,  &
                              glissade_interior_dissipation_first_order, &
                              glissade_flow_factor,  &
                              glissade_pressure_melting_point
    use glissade_calving, only: verbose_calving
    use glam_grid_operators, only: glam_geometry_derivs
    use felix_dycore_interface, only: felix_velo_driver

    !WHL - debug
    use glissade_bmlt_float, only: glissade_basal_melting_float

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

    integer :: i, j, k, n
    integer :: itest, jtest, rtest

    integer, dimension(model%general%ewn, model%general%nsn) :: &
         ice_mask,           & ! = 1 where thck > thklim, else = 0
         floating_mask,      & ! = 1 where ice is present and floating, else = 0
         calving_front_mask, & ! = 1 where ice is floating and borders an ocean cell, else = 0
         ocean_mask,         & ! = 1 where topg is below sea level and ice is absent
         land_mask,          & ! = 1 where topg is at or above sea level
         active_ice_mask       ! = 1 where ice is dynamically active, else = 0

    real(dp), dimension(model%general%ewn, model%general%nsn) ::  &
       thck_calving_front      ! effective thickness of ice at the calving front

    real(dp) :: &
         dsigma,                   & ! layer thickness in sigma coordinates
         tau_xx, tau_yy, tau_xy,   & ! stress tensor components
         strain_rate_xx, strain_rate_yy, strain_rate_xy  ! strain rate tensor components
 
    real(dp) :: &
         a, b, c, root,   & ! terms in quadratic formula
         lambda1, lambda2   ! eigenvalues of horizontal strain rate tensor

    integer :: ii, jj

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 1. First part of diagnostic solve: 
    !    Now that advection is done, update geometry- and temperature-related 
    !    diagnostic fields that are needed for the velocity solve.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! ------------------------------------------------------------------------
    ! Halo update for ice thickness
    ! Note: The halo update for topg is done at initialization and again (as needed)
    !       after computing isostasy.
    ! ------------------------------------------------------------------------

    call parallel_halo(model%geometry%thck)

    ! ------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    ! ------------------------------------------------------------------------
    !TODO - These are currently updated after transport. Needed for calving/isostasy, or not until here?

    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       & 
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

    ! ------------------------------------------------------------------------
    ! Update some geometry derivatives
    ! ------------------------------------------------------------------------
    !Note - The fields computed by glam_geometry_derivs are not required by
    !       the Glissade velocity solver (which computes them internally).  
    !       However, some of the fields (stagthck, dusrfdew and dusrfdns) 
    !       are needed during the next timestep by glissade_therm
    !       if we're doing shallow-ice dissipation.
    !TODO - Replace this glam_geometry_derivs call with calls to Glissade subroutines?
    !       (The glam_velo driver includes its own call to glam_geometry_derivs.) 

    call glam_geometry_derivs(model)

    ! ------------------------------------------------------------------------
    ! Update some masks that are used for subsequent calculations
    ! ------------------------------------------------------------------------

    call glissade_get_masks(model%general%ewn,   model%general%nsn,     &
                            model%geometry%thck, model%geometry%topg,   &
                            model%climate%eus,   model%numerics%thklim, &
                            ice_mask,                                   &
                            floating_mask = floating_mask,              &
                            ocean_mask = ocean_mask,                    &
                            land_mask = land_mask,                      &
                            which_ho_calving_front = model%options%which_ho_calving_front, &
                            calving_front_mask = calving_front_mask,    &
                            thck_calving_front = thck_calving_front)

    ! ------------------------------------------------------------------------
    ! Compute the fraction of grounded ice in each cell
    ! (requires that thck and topg are up to date in halo cells).
    ! This is used in the velocity solver to compute the basal stress BC.
    !
    ! See comments in subroutine glissade_grounded_fraction for details
    ! on the whichground and whichflotation_function options.
    !
    ! Computing f_ground here ensures that it is always available as a diagnostic, even if
    ! the velocity solver is not called (e.g., on the first time step of a restart).
    ! ------------------------------------------------------------------------

    call glissade_grounded_fraction(model%general%ewn,             &
                                    model%general%nsn,             &
                                    ice_mask,                      &
                                    floating_mask,                 &
                                    land_mask,                     &
                                    model%options%which_ho_ground, &
                                    model%geometry%f_ground)

    ! ------------------------------------------------------------------------ 
    ! Calculate Glen's A
    !
    ! Notes:
    ! (1) Because flwa is not a restart variable in Glissade, no check is included 
    !      here for whether to calculate it on initial time (as is done in Glide).
    ! (2) We are passing in only vertical elements (1:upn-1) of the temp array,
    !       so that it has the same vertical dimensions as flwa.
    ! (3) The flow enhancement factor is 1 by default.
    ! (4) The waterfrac field is ignored unless whichtemp = TEMP_ENTHALPY.
    ! (5) Inputs and outputs of glissade_flow_factor should have SI units.
    ! ------------------------------------------------------------------------

    call glissade_flow_factor(model%options%whichflwa,            &
                              model%options%whichtemp,            &
                              model%numerics%stagsigma,           &
                              model%geometry%thck * thk0,         &  ! scale to m
                              model%temper%temp(1:model%general%upn-1,:,:),  &
                              model%temper%flwa,                  &  ! Pa^{-n} s^{-1}
                              model%paramets%default_flwa / scyr, &  ! scale to Pa^{-n} s^{-1}
                              model%paramets%flow_enhancement_factor,       &
                              model%paramets%flow_enhancement_factor_float, &
                              floating_mask,                      &
                              model%temper%waterfrac(:,:,:))

    !TODO - flwa halo update not needed?
    ! Halo update for flwa
    call parallel_halo(model%temper%flwa)

    ! ------------------------------------------------------------------------
    ! Do some additional operations if this is the first time step.
    ! The model thickness and temperature fields will have been initialized, but the
    !  thermal and transport solvers have not been called yet.
    ! Note: Some operations must be done in this subroutine when restarting; others are skipped.
    ! ------------------------------------------------------------------------

    if (model%numerics%time == model%numerics%tstart) then

       ! Compute the pressure melting point temperature, which is needed
       ! by certain basal sliding laws.

       do j = 1, model%general%nsn
          do i = 1, model%general%ewn
             call glissade_pressure_melting_point(model%geometry%thck(i,j) * thk0, &
                                                  model%temper%bpmp(i,j))
          enddo
       enddo

       ! If the velocity fields have been read in on the extended staggered mesh,
       ! then copy them to the standard staggered mesh.
       !
       ! Note: For problems with nonzero velocity along the global boundaries (e.g., MISMIP on a periodic domain),
       !        exact restart requires that the restart velocity field lies on an extended staggered mesh with
       !        an extra row and column of velocity points along the north and east boundary of the domain.
       !       In that case, uvel_extend and vvel_extend should be written to the restart file by setting
       !        restart_extend_velo = 1 in the config file. They will then be read as input fields and at this
       !        point need to be copied into uvel and vvel.
       !       (It would have been cleaner to give uvel and vvel the same dimensions as the scalar mesh,
       !        but that design decision was made many years ago and would take a lot of work to change.)

       ! If uvel_extend and vvel_extend are input fields, then copy them into uvel and vvel.
       ! Halo updates are then needed to make sure the velocities are correct along the boundaries.
 
       if  ( (maxval(abs(model%velocity%uvel_extend)) /= 0.0d0) .or. & 
             (maxval(abs(model%velocity%vvel_extend)) /= 0.0d0) ) then
          call write_log('Using uvel_extend, vvel_extend from input or restart file at initial time')
          model%velocity%uvel(:,:,:) = model%velocity%uvel_extend(:,1:model%general%ewn-1,1:model%general%nsn-1)
          model%velocity%vvel(:,:,:) = model%velocity%vvel_extend(:,1:model%general%ewn-1,1:model%general%nsn-1)
!       elseif ( (maxval(abs(model%velocity%uvel)) /= 0.0d0) .or. & 
!                (maxval(abs(model%velocity%vvel)) /= 0.0d0) ) then
!          call write_log('Using uvel, vvel from input or restart file at initial time')
       endif

       call staggered_parallel_halo(model%velocity%uvel)
       call staggered_parallel_halo(model%velocity%vvel)

       ! The DIVA solver option requires some additional fields on the staggered mesh for exact restart.
       ! If these fields were input on the extended mesh, they need to be copied to the standard staggered mesh.
       ! Halo updates are then needed to make sure they have the correct values along the boundaries.

       if (model%options%which_ho_approx == HO_APPROX_DIVA) then

          if  ( (maxval(abs(model%velocity%uvel_2d_extend)) /= 0.0d0) .or. & 
                (maxval(abs(model%velocity%vvel_2d_extend)) /= 0.0d0) ) then
             call write_log('Using uvel_2d_extend, vvel_2d_extend from input or restart file at initial time')
             model%velocity%uvel_2d(:,:) = model%velocity%uvel_2d_extend(1:model%general%ewn-1,1:model%general%nsn-1)
             model%velocity%vvel_2d(:,:) = model%velocity%vvel_2d_extend(1:model%general%ewn-1,1:model%general%nsn-1)
!          elseif ( (maxval(abs(model%velocity%uvel_2d)) /= 0.0d0) .or. & 
!                   (maxval(abs(model%velocity%vvel_2d)) /= 0.0d0) ) then
!             call write_log('Using uvel_2d, vvel_2d from input or restart file at initial time')
          endif

          if  ( (maxval(abs(model%stress%btractx_extend)) /= 0.0d0) .or. & 
                (maxval(abs(model%stress%btracty_extend)) /= 0.0d0) ) then
             model%stress%btractx(:,:) = model%stress%btractx_extend(1:model%general%ewn-1,1:model%general%nsn-1)
             model%stress%btracty(:,:) = model%stress%btracty_extend(1:model%general%ewn-1,1:model%general%nsn-1)
          endif

          call staggered_parallel_halo(model%velocity%uvel_2d)
          call staggered_parallel_halo(model%velocity%vvel_2d)
          call staggered_parallel_halo(model%stress%btractx)
          call staggered_parallel_halo(model%stress%btracty)

       endif   ! DIVA approx
             
    endif   ! time = tstart

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 2. Second part of diagnostic solve: 
    !    Now that geometry- and temperature-related diagnostic fields are updated, 
    !    solve velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Do not solve velocity for initial time on a restart because that breaks an exact restart.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model%numerics%time == model%numerics%tstart) ) then
  
       ! Do not solve for velocity, because this would break exact restart

    else

       ! If this is not a restart or we are not at the initial time, then proceed normally

       if ( (model%numerics%time == model%numerics%tstart) .and. &
         ( (maxval(abs(model%velocity%uvel)) /= 0.0d0) .or. & 
           (maxval(abs(model%velocity%vvel)) /= 0.0d0) ) ) then
          ! If velocity was input and this is NOT a restart, then use the input field as the first guess at the initial time.
          ! This happens automatically, but let the user know.
          ! Using this value will change the answer only within the tolerance of the nonlinear solve.  
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
          call staggered_parallel_halo_extrapolate(model%basal_physics%mintauf)
       endif

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
 
       ! Compute internal heat dissipation
       ! This is used in the prognostic temperature calculation during the next time step.
       ! Note: These glissade subroutines assume SI units on input and output

       model%temper%dissip(:,:,:) = 0.d0

       if (model%options%which_ho_disp == HO_DISP_SIA) then

          call glissade_interior_dissipation_sia(model%general%ewn,              &
                                                 model%general%nsn,              &
                                                 model%general%upn,              &
                                                 model%numerics%stagsigma(:),    &
                                                 ice_mask,                       &
                                                 model%geomderv%stagthck * thk0, & ! scale to m
                                                 model%temper%flwa * vis0,       & ! scale to Pa^{-n} s^{-1}
                                                 model%geomderv%dusrfdew * thk0/len0, & ! scale to m/m
                                                 model%geomderv%dusrfdns * thk0/len0, & ! scale to m/m
                                                 model%temper%dissip)
          
       else    ! first-order dissipation                                                                                                                                                               
          call glissade_interior_dissipation_first_order(model%general%ewn,          &
                                                         model%general%nsn,          &
                                                         model%general%upn,          &
                                                         ice_mask,                   &
                                                         model%stress%tau%scalar * tau0,  &  ! scale to Pa
                                                         model%stress%efvs * evs0,   &  ! scale to Pa s
                                                         model%temper%dissip)
          
       endif    ! which_ho_disp
          
       ! If running Glam, compute the basal friction heat flux
       ! (Glissade computes this flux as part of the velocity solution.)
       
       if (model%options%whichdycore == DYCORE_GLAM) then
          call glam_basal_friction(model%general%ewn,                             &
                                   model%general%nsn,                             &
                                   ice_mask,                                      &
                                   floating_mask,                                 &
                                   model%velocity%uvel(model%general%upn,:,:),    &
                                   model%velocity%vvel(model%general%upn,:,:),    &
                                   model%velocity%btraction(:,:,:),               &
                                   model%temper%bfricflx(:,:) )
       endif
       
    endif     ! is_restart

    if (this_rank==rtest .and. verbose_glissade) then
       i = itest
       j = jtest
       print*, 'itest, jtest =', i, j
       print*, 'k, dissip (deg/yr):'
       do k = 1, model%general%upn-1
          print*, k, model%temper%dissip(k,i,j)*scyr
       enddo
       print*, 'ubas, vbas =', model%velocity%uvel(model%general%upn,i,j),  &
            model%velocity%vvel(model%general%upn,i,j)
       print*, 'btraction =',  model%velocity%btraction(:,i,j)
       print*, 'bfricflx =', model%temper%bfricflx(i,j)
       print*, ' '
       print*, 'After glissade velocity solve (or restart): uvel, k = 1:'
       write(6,'(a8)',advance='no') '          '
!!          do i = 1, model%general%ewn-1
       do i = itest-5, itest+5
          write(6,'(i12)',advance='no') i
       enddo
       print*, ' '
!!          do j = model%general%nsn-1, 1, -1
       do j = jtest+2, jtest-2, -1
          write(6,'(i8)',advance='no') j
!!             do i = 1, model%general%ewn-1
          do i = itest-5, itest+5
             write(6,'(f12.3)',advance='no') model%velocity%uvel(1,i,j) * (vel0*scyr)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'After glissade velocity solve (or restart): vvel, k = 1:'
       write(6,'(a8)',advance='no') '          '
!!          do i = 1, model%general%ewn-1
       do i = itest-5, itest+5
          write(6,'(i12)',advance='no') i
       enddo
       print*, ' '
!!          do j = model%general%nsn-1, 1, -1
       do j = jtest+2, jtest-2, -1
          write(6,'(i8)',advance='no') j
!!             do i = 1, model%general%ewn-1
          do i = itest-5, itest+5
             write(6,'(f12.3)',advance='no') model%velocity%vvel(1,i,j) * (vel0*scyr)
          enddo
          print*, ' '
       enddo
       
    endif  ! this_rank = rtest & verbose_glissade

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 3. Third part of diagnostic solve: 
    ! Now that velocity is solved, calculate any diagnostic fields that are
    ! a function of velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Calculate wvel, assuming grid velocity is 0.
    ! This is calculated relative to the ice sheet base, rather than a fixed reference location.
    ! Note: This current implementation for wvel only supports whichwvel=VERTINT_STANDARD
    ! Note: For glissade, wvel is diagnostic only.
    !       Pass in the basal melt rate for grounded ice only, as in Glide.
    call wvelintg(model%velocity%uvel,                        &
                  model%velocity%vvel,                        &
                  model%geomderv,                             &
                  model%numerics,                             &
                  model%velowk,                               &
                  model%geometry%thck * 0.0d0,                &  ! Just need a 2d array of all 0's for wgrd
                  model%geometry%thck,                        &
                  model%basal_melt%bmlt,                      &
                  model%velocity%wvel)
    ! Note: halos may be wrong for wvel, but since it is currently only used as an output diagnostic variable, that is OK.

    ! compute the velocity norm (for diagnostic output)
    model%velocity%velnorm(:,:,:) = sqrt(model%velocity%uvel(:,:,:)**2 + model%velocity%vvel(:,:,:)**2)

    ! compute the mean velocity

    model%velocity%uvel_mean(:,:) = 0.0d0
    model%velocity%vvel_mean(:,:) = 0.0d0

    k = 1    ! top surface velocity associated with top half of layer 1
    model%velocity%uvel_mean(:,:) = model%velocity%uvel_mean(:,:) &
                                  + model%numerics%stagsigma(k) * model%velocity%uvel(k,:,:)
    model%velocity%vvel_mean(:,:) = model%velocity%vvel_mean(:,:) &
                                  + model%numerics%stagsigma(k) * model%velocity%vvel(k,:,:)

    do k = 2, model%general%upn-1
       model%velocity%uvel_mean(:,:) = model%velocity%uvel_mean(:,:) &
                                     + (model%numerics%stagsigma(k) - model%numerics%stagsigma(k-1)) * model%velocity%uvel(k,:,:)
       model%velocity%vvel_mean(:,:) = model%velocity%vvel_mean(:,:) &
                                     + (model%numerics%stagsigma(k) - model%numerics%stagsigma(k-1)) * model%velocity%vvel(k,:,:)
    enddo

    k = model%general%upn  ! basal velocity associated with bottom half of layer (upn-1)
    model%velocity%uvel_mean(:,:) = model%velocity%uvel_mean(:,:) &
                                  + (1.0d0 - model%numerics%stagsigma(k-1)) * model%velocity%uvel(k,:,:)
    model%velocity%vvel_mean(:,:) = model%velocity%vvel_mean(:,:) &
                                  + (1.0d0 - model%numerics%stagsigma(k-1)) * model%velocity%vvel(k,:,:)

    ! Compute the vertically integrated stress tensor (Pa) and its eigenvalues.
    ! These are used for some calving schemes.

    if ( (model%options%is_restart == RESTART_TRUE) .and. &
         (model%numerics%time == model%numerics%tstart) ) then

       ! do nothing, since the tau eigenvalues are read from the restart file

    else  ! compute the eigenvalues given the stress just computed in the velocity solver

       model%calving%tau_eigen1(:,:) = 0.0d0
       model%calving%tau_eigen2(:,:) = 0.0d0

       do j = 1, model%general%nsn
          do i = 1, model%general%ewn

             ! compute vertically averaged stress components
             tau_xx = 0.0d0
             tau_yy = 0.0d0
             tau_xy = 0.0d0

             do k = 1, model%general%upn-1
                dsigma = model%numerics%sigma(k+1) - model%numerics%sigma(k)
                tau_xx = tau_xx + tau0 * model%stress%tau%xx(k,i,j) * dsigma
                tau_yy = tau_yy + tau0 * model%stress%tau%yy(k,i,j) * dsigma
                tau_xy = tau_xy + tau0 * model%stress%tau%xy(k,i,j) * dsigma
             enddo

             ! compute the eigenvalues of the vertically integrated stress tensor
             a = 1.0d0
             b = -(tau_xx + tau_yy)
             c = tau_xx*tau_yy - tau_xy*tau_xy
             if (b*b - 4.0d0*a*c > 0.0d0) then   ! two real eigenvalues
                root = sqrt(b*b - 4.0d0*a*c)
                lambda1 = (-b + root) / (2.0d0*a)
                lambda2 = (-b - root) / (2.0d0*a)
                if (lambda1 > lambda2) then
                   model%calving%tau_eigen1(i,j) = lambda1
                   model%calving%tau_eigen2(i,j) = lambda2
                else
                   model%calving%tau_eigen1(i,j) = lambda2
                   model%calving%tau_eigen2(i,j) = lambda1
                endif
             endif  ! b^2 - 4ac > 0

          enddo   ! i
       enddo   ! j

       ! Extrapolate tau eigenvalues to inactive CF cells where the stress tensor is not computed.
       do j = 2, model%general%nsn-1
          do i = 2, model%general%ewn-1
             if (calving_front_mask(i,j) == 1 .and. &
                  model%calving%tau_eigen1(i,j) == 0.0d0 .and. model%calving%tau_eigen2(i,j) == 0.0d0) then

                !  Look for nonzero values in an upstream cell
                do jj = j-1, j+1
                   do ii = i-1, i+1
                      if (thck_calving_front(i,j) > 0.0d0 .and. &
                           model%geometry%thck(ii,jj) == thck_calving_front(i,j)) then
                         model%calving%tau_eigen1(i,j) = model%calving%tau_eigen1(ii,jj)
                         model%calving%tau_eigen2(i,j) = model%calving%tau_eigen2(ii,jj)
                      endif
                   enddo
                enddo

             endif  ! inactive CF cell
          enddo
       enddo

    endif   ! restart

    call parallel_halo(model%calving%tau_eigen1)
    call parallel_halo(model%calving%tau_eigen2)

    !WHL - debug
    if (this_rank == rtest .and. verbose_calving) then
       print*, ' '
       print*, 'tau eigen1 (Pa), i, j, rtest =:', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(e11.3)',advance='no') model%calving%tau_eigen1(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'tau eigen2 (Pa), i, j, rtest =:', itest, jtest, rtest
       do j = jtest+3, jtest-3, -1
          write(6,'(i8)',advance='no') j
          do i = itest-3, itest+3
             write(6,'(e11.3)',advance='no') model%calving%tau_eigen2(i,j)
          enddo
          print*, ' '
       enddo
    endif  ! this_rank = rtest

    ! Compute the 3D strain rate tensor (s^{-1})
    ! Note: The stress tensor tau is derived by taking strain rates at quadrature points in the velocity solve.
    !       The strain rate tensor is simply diagnosed from the stress tensor.
    where (model%stress%efvs > 0.0d0) 
       model%velocity%strain_rate%scalar = tau0 * model%stress%tau%scalar / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%xz = tau0 * model%stress%tau%xz / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%yz = tau0 * model%stress%tau%yz / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%xx = tau0 * model%stress%tau%xx / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%yy = tau0 * model%stress%tau%yy / (2.d0 * evs0 * model%stress%efvs)
       model%velocity%strain_rate%xy = tau0 * model%stress%tau%xy / (2.d0 * evs0 * model%stress%efvs)
    elsewhere
       model%velocity%strain_rate%scalar = 0.0d0
       model%velocity%strain_rate%xz = 0.0d0
       model%velocity%strain_rate%yz = 0.0d0
       model%velocity%strain_rate%xx = 0.0d0
       model%velocity%strain_rate%yy = 0.0d0
       model%velocity%strain_rate%xy = 0.0d0
    endwhere

    ! Compute various vertical means.
    ! TODO - Write a utility subroutine for vertical averaging

    ! Compute the vertical mean effective viscosity
    model%stress%efvs_vertavg = 0.0d0
    do j = 1, model%general%nsn
       do i = 1, model%general%ewn
          do k = 1, model%general%upn-1
             model%stress%efvs_vertavg(i,j) = model%stress%efvs_vertavg(i,j)  &
                                            + model%stress%efvs(k,i,j) * (model%numerics%sigma(k+1) - model%numerics%sigma(k))
          enddo
       enddo
    enddo

    ! Compute the vertically integrated divergence of the horizontal velocity field.
    ! Note: Units of divu and strain_rate components are s^{-1}.
    model%velocity%divu(:,:) = 0.0d0

    do j = 1, model%general%nsn
       do i = 1, model%general%ewn
          do k = 1, model%general%upn-1
             dsigma = model%numerics%sigma(k+1) - model%numerics%sigma(k)
             model%velocity%divu(i,j) = model%velocity%divu(i,j) + &
                  (model%velocity%strain_rate%xx(k,i,j) + model%velocity%strain_rate%yy(k,i,j)) * dsigma
          enddo
       enddo
    enddo

    ! magnitude of basal traction
    model%stress%btract(:,:) = sqrt(model%stress%btractx(:,:)**2 + model%stress%btracty(:,:)**2)

    ! Copy uvel and vvel to arrays uvel_extend and vvel_extend.
    ! These arrays have horizontal dimensions (nx,ny) instead of (nx-1,ny-1).
    ! They are needed for exact restart if we have nonzero velocities along the
    !  north and east edges of the global domain, as in some test problems.
    
    model%velocity%uvel_extend(:,:,:) = 0.d0
    model%velocity%vvel_extend(:,:,:) = 0.d0

    do j = 1, model%general%nsn-1
       do i = 1, model%general%ewn-1
          model%velocity%uvel_extend(:,i,j) = model%velocity%uvel(:,i,j)
          model%velocity%vvel_extend(:,i,j) = model%velocity%vvel(:,i,j)             
       enddo
    enddo

    ! Copy some additional 2D arrays to the extended grid if using the DIVA solver.

    if (model%options%which_ho_approx == HO_APPROX_DIVA) then

       model%velocity%uvel_2d_extend(:,:) = 0.d0
       model%velocity%vvel_2d_extend(:,:) = 0.d0
       do j = 1, model%general%nsn-1
          do i = 1, model%general%ewn-1
             model%velocity%uvel_2d_extend(i,j) = model%velocity%uvel_2d(i,j)
             model%velocity%vvel_2d_extend(i,j) = model%velocity%vvel_2d(i,j)             
          enddo
       enddo

       model%stress%btractx_extend(:,:) = 0.d0
       model%stress%btracty_extend(:,:) = 0.d0
       do j = 1, model%general%nsn-1
          do i = 1, model%general%ewn-1
             model%stress%btractx_extend(i,j) = model%stress%btractx(i,j)
             model%stress%btracty_extend(i,j) = model%stress%btracty(i,j)   
          enddo
       enddo

    endif

    ! If beta is not passed in from an external file, then copy beta_internal to beta.
    ! Note: beta_internal, which is weighted by the grounded ice fraction, is the actual beta field
    !        in the glissade velocity calculation.  But if users specify 'beta' instead of
    !        'beta_internal' as an output field, this copy ensures that they get the output
    !         they expect.
    !       The copy would break exact restart, however, if beta is read from an external file.
    !        In that case, users must specify 'beta_internal' to see the weighted beta field.

    if (model%options%which_ho_babc /= HO_BABC_BETA_EXTERNAL) then
       model%velocity%beta(:,:) = model%velocity%beta_internal(:,:)
    endif

    ! DIVA needs a halo update for efvs, since the latest guess (in both local and halo cells)
    ! is used to start iterating for efvs in the next time step.
    call parallel_halo(model%stress%efvs)

    !TODO - I don't think we need to update ubas, vbas, or velnorm, since these are diagnostic only
    call staggered_parallel_halo(model%velocity%velnorm)
    call staggered_parallel_halo(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 4. Fourth part of diagnostic solve: 
    ! Diagnose some quantities that are not velocity-dependent, but may be desired for output
    ! ------------------------------------------------------------------------
    ! ------------------------------------------------------------------------

    ! basal ice temperature
    ! This is the same as temp(upn,:,:), the lowest-level of the prognostic temperature array.
    ! However, it is set to zero for ice-free columns (unlike temp(upn) = min(artm,0.0) for ice-free columns)
    ! TODO - Make btemp a prognostic array, and limit the 3D temp array to internal layer temperatures?
    do j = 1, model%general%nsn
       do i = 1, model%general%ewn
          if (model%geometry%thck(i,j) > 0.0d0) then
             model%temper%btemp(i,j) = model%temper%temp(model%general%upn,i,j)
          else
             model%temper%btemp(i,j) = 0.0d0
          endif
       enddo
    enddo

    ! thickness tendency dH/dt from one step to the next (m/s)
    ! Note: This diagnostic will not be correct on the first step of a restart

    do j = 1, model%general%nsn
       do i = 1, model%general%ewn
          model%geometry%dthck_dt(i,j) = (model%geometry%thck(i,j) - model%geometry%thck_old(i,j))*thk0 &
                                       / (model%numerics%dt * tim0)
       enddo
    enddo

    ! surface mass balance in units of mm/yr w.e.
    ! (model%climate%acab * scale_acab) has units of m/yr of ice
    ! Note: This is not necessary (and can destroy exact restart) if the SMB was already input in units of mm/yr
    if (model%options%smb_input /= SMB_INPUT_MMYR_WE) then
       model%climate%smb(:,:) = (model%climate%acab(:,:) * scale_acab) * (1000.d0 * rhoi/rhow)
    endif

    ! surface, basal and calving mass fluxes (kg/m^2/s)
    ! positive for mass gain, negative for mass loss
    model%geometry%sfc_mbal_flux(:,:) = rhoi * model%climate%acab_applied(:,:)*thk0/tim0
    model%geometry%basal_mbal_flux(:,:) = rhoi * (-model%basal_melt%bmlt_applied(:,:)) * thk0/tim0
    model%geometry%calving_flux(:,:) = rhoi * (-model%calving%calving_thck(:,:)*thk0) / (model%numerics%dt*tim0)

    ! calving rate (m/yr ice; positive for calving)
    model%calving%calving_rate(:,:) = (model%calving%calving_thck(:,:)*thk0) / (model%numerics%dt*tim0/scyr)

    ! set integer masks in the geometry derived type

    ! unstaggered grid
    do j = 1, model%general%nsn
       do i = 1, model%general%ewn
          if (ice_mask(i,j) == 1) then
             model%geometry%ice_mask(i,j) = 1
             if (floating_mask(i,j) == 1) then
                model%geometry%grounded_mask(i,j) = 0
                model%geometry%floating_mask(i,j) = 1
             else
                model%geometry%grounded_mask(i,j) = 1
                model%geometry%floating_mask(i,j) = 0
             endif
          else  ! ice_mask = 0
             model%geometry%ice_mask(i,j) = 0
             model%geometry%grounded_mask(i,j) = 0
             model%geometry%floating_mask(i,j) = 0
          endif
       enddo
    enddo

    ! staggered grid
    ! set ice_mask_stag = 1 at vertices with ice_mask = 1 in any neighbor cell
    do j = 1, model%general%nsn - 1
       do i = 1, model%general%ewn - 1
          if (ice_mask(i,j+1)==1 .or. ice_mask(i+1,j+1)==1 .or. &
              ice_mask(i,j)  ==1 .or. ice_mask(i+1,j)  ==1) then
             model%geometry%ice_mask_stag(i,j) = 1
          else
             model%geometry%ice_mask_stag(i,j) = 0
          endif
       enddo
    enddo

    ! Compute grounding line fluxes
    ! Note: gl_flux_east and gl_flux_north are signed fluxes computed at cell edges;
    !       gl_flux is cell-based and is found by summing magnitudes of edge fluxes.

    call glissade_grounding_line_flux(model%general%ewn,    model%general%nsn,   &
                                      model%numerics%dew,   model%numerics%dns,  &
                                      model%numerics%sigma,                      &
                                      model%geometry%thck,                       &
                                      model%velocity%uvel,  model%velocity%vvel, &
                                      ice_mask,             floating_mask,       &
                                      ocean_mask,                                &
                                      model%geometry%gl_flux_east,               &
                                      model%geometry%gl_flux_north,              &
                                      model%geometry%gl_flux                      )

    !------------------------------------------------------------------------
    ! Update the upper and lower ice surface
    ! Note that glide_calclsrf loops over all cells, including halos,
    !  so halo updates are not needed for lsrf and usrf.
    !
    !
    ! TODO(wjs, 2017-05-21) I don't think we should need to update lsrf and usrf
    ! here. However, glissade_velo_higher_solve and glissade_velo_sia_solve (called from
    ! glissade_velo_driver) multiply/divide topg (and other variables) by their scale
    ! factors on entry to / exit from the routine. This can lead to roundoff-level changes
    ! in topg and other variables.
    !
    ! If we don't update usrf here, then we can get roundoff-level changes in exact
    ! restart tests when running inside a climate model: In the straight-through run
    ! (without an intervening restart), the value of usrf sent to the coupler is the one
    ! set earlier in this routine, which doesn't incorporate these roundoff-level changes
    ! to topg. The restarted run, in contrast, reads the slightly-modified topg from the
    ! restart file and recomputes usrf in initialization; thus, the values of usrf that
    ! the coupler sees in the first year differ slightly from those in the
    ! straight-through run.
    !
    ! A cleaner solution could be to avoid applying these rescalings to the fundamental
    ! model variables in glissade_velo_higher_solve and glissade_velo_sia_solve - instead,
    ! introducing temporary variables in those routines to hold the scaled
    ! quantities. Then I think it would be safe to remove the following code that updates
    ! lsrf and usrf. Or, if we completely removed these scale factors from CISM, then
    ! again I think it would be safe to remove the following code.
    ! ------------------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg,       &
                        model%climate%eus,   model%geometry%lsrf)
    model%geometry%usrf(:,:) = max(0.d0, model%geometry%thck(:,:) + model%geometry%lsrf(:,:))

  end subroutine glissade_diagnostic_variable_solve

!=======================================================================

  subroutine glissade_cleanup_icefree_cells(model)

    ! Clean up prognostic variables in ice-free cells.
    ! This means seting most tracers to zero (or min(artm,0) for the case of temperature).

    use parallel, only: parallel_halo

    type(glide_global_type), intent(inout) :: model   ! model instance

    integer :: nx, ny
    integer :: i, j

    nx = model%general%ewn
    ny = model%general%nsn

    ! Make sure the ice thickness is updated in halo cells
    call parallel_halo(model%geometry%thck)

    ! Set prognostic variables in ice-free columns to default values (usually zero).
    do j = 1, ny
       do i = 1, nx

          if (model%geometry%thck_old(i,j) > 0.0d0 .and. model%geometry%thck(i,j) == 0.0d0) then

             ! basal water
             model%temper%bwat(i,j) = 0.0d0

             ! thermal variables
             if (model%options%whichtemp == TEMP_INIT_ZERO) then
                model%temper%temp(:,i,j) = 0.0d0
             else
                model%temper%temp(:,i,j) = min(model%climate%artm(i,j), 0.0d0)
             endif

             if (model%options%whichtemp == TEMP_ENTHALPY) then
                model%temper%waterfrac(:,i,j) = 0.0d0
             endif

             ! other tracers
             ! Note: Tracers should be added here as they are added to the model

             if (model%options%whichcalving == CALVING_DAMAGE) then
                model%calving%damage(:,i,j) = 0.0d0
             endif

             if (model%options%which_ho_ice_age == HO_ICE_AGE_COMPUTE) then
                model%geometry%ice_age(:,i,j) = 0.0d0
             endif

          endif    ! thck = 0

       enddo
    enddo

  end subroutine glissade_cleanup_icefree_cells

!=======================================================================

end module glissade

!=======================================================================
