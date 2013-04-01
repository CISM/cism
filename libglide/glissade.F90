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
! In general, all parallel_halo updates should go in this module, except
!  for the uvel and vvel updates required by the velocity solvers.
! Also, we would like to eliminate use of the 'model' derived type below this level.
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glissade.f90 - part of the Glimmer-CISM ice model        + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

module glissade

  ! Driver for Glissade (parallel, higher-order) dynamical core

  use glide_types
  use glide_io
  use glide_lithot_io
  use glide_lithot
  use glimmer_config
  use glimmer_global

!TODO - Remove scaling

  implicit none

  integer, private, parameter :: dummyunit=99

contains

!=======================================================================

  subroutine glissade_initialise(model)

    ! initialise Glissade model instance

!TODO - Are all of these needed?
    use parallel
    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio
    use glide_velo, only: init_velo  !TODO - Remove this
    use glissade_temp, only: glissade_init_temp, glissade_calcflwa
    use glimmer_log
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground
    use glam_strs2, only : glam_velo_init
    use glissade_velo_higher, only: glissade_velo_higher_init

!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
!!    use fo_upwind_advect, only : fo_upwind_advect_init
!!    use glam_Basal_Proc, only : Basal_Proc_init

!WHL - debug
    use glimmer_paramets, only: thk0

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

!lipscomb - TODO - build glimmer_vers file or put this character elsewhere?
    character(len=100), external :: glimmer_version_char

!WHL - debug
!    logical, parameter :: test_parallel = .true.   ! if true, call test_parallel subroutine
    logical, parameter :: test_parallel = .false.   ! if true, call test_parallel subroutine
    integer :: i, j, nx, ny

!WHL - for artificial adjustment to ismip-hom surface elevation
    logical, parameter :: ismip_hom_adjust_usrf = .false.
    real(dp) :: usrf_ref

    call write_log(trim(glimmer_version_char()))

!SCALING - to be removed?
    ! initialise scales
    call glimmer_init_scales

    ! scale parameters (some conversions to SI units)
    call glide_scale_params(model)

    ! set up coordinate systems
    ! time to change to the parallel values of ewn and nsn

    call distributed_grid(model%general%ewn,model%general%nsn)

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

!TODO - Can remove; PBJ only
    !Initialize boundary condition fields to be NaN everywhere 
!!    model%geometry%marine_bc_normal = NaN
    
    ! set uniform basal heat flux 
    model%temper%bheatflx = model%paramets%geot

!TODO - glimmer_load_sigma?
    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)

!TODO - glimmer_io_readall?
    ! and read first time slice
    call glide_io_readall(model,model)

    ! Write projection info to log
    call glimmap_printproj(model%projection)

    ! read lithot if required
    if (model%options%gthf > 0) then
       call glide_lithot_io_readall(model,model)
    end if

    ! handle relaxed/equilibrium topo
    ! Initialise isostasy first
    call init_isostasy(model)

    select case(model%options%whichrelaxed)
    case(1) ! Supplied topography is relaxed
       model%isos%relx = model%geometry%topg
    case(2) ! Supplied topography is in equilibrium
       call not_parallel(__FILE__,__LINE__)
       call isos_relaxed(model)
    end select

    ! open all output files
    call openall_out(model)

!TODO - glimmer_io_createall?
    ! create glide variables
    call glide_io_createall(model)

    ! initialise glissade components

!TODO - Most of what's done in init_velo is needed for SIA only
!       Can remove this call provided velowk is not used elsewhere (e.g., to call wvelintg)
    call init_velo(model)

    call glissade_init_temp(model)  ! temperature lives at layer centers

!WHL - Removed call to glide_initialise_backstress; used only by old glide_marinlim case(5)

    if (model%options%gthf > 0) then
       call not_parallel(__FILE__,__LINE__)
       call glide_lithot_io_createall(model)
       call init_lithot(model)
    end if

    if (model%options%whichdycore == DYCORE_GLAM ) then  ! glam finite-difference

        call glam_velo_init(model%general%ewn,    model%general%nsn,  &
                            model%general%upn,                        &
                            model%numerics%dew,   model%numerics%dns, &
                            model%numerics%sigma)

    elseif (model%options%whichdycore == DYCORE_GLISSADE ) then  ! glissade finite-element

!WHL - Removed scaling of dew and dns
        call glissade_velo_higher_init

    endif

!WHL - This option is disabled for now.
    ! *mb* added; initialization of basal proc. module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!!                              model%numerics%ntem)
!!    end if      

!TODO - Change to glissade_set_mask?

    ! calculate mask
       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)
!       call horiz_bcs_unstag_scalar(model%geometry%thkmask)

!TODO- Not sure why this needs to be called here.
    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

!TODO - inline calclsrf?
    ! and calculate lower and upper ice surface

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

!WHL - The following is a hack to ensure that usrf is uniform (to double precision) for 
!      a given value of i for the ISMIP-HOM tests.  We take one value in each column as 
!      the benchmark and set all other values in that column to the same value.
!      Then we correct the thickness to ensure that lsrf + thck = usrf to double precision.
!      
!      A better way would be to read double-precision data from the input files.
!

    if (ismip_hom_adjust_usrf) then
       do i = 1, model%general%ewn
          usrf_ref = model%geometry%usrf(i,nhalo+1)
          do j = 1, model%general%nsn
             model%geometry%usrf(i,:) = usrf_ref  ! same usrf for all j in this column
          enddo
       enddo
       model%geometry%thck(:,:) = model%geometry%usrf(:,:) - model%geometry%lsrf(:,:) 
    endif

!TODO - Pretty sure these are SIA only
    ! initialise thckwk variables; used in timeders subroutine
!!    model%thckwk%olds(:,:,1) = model%geometry%thck(:,:)
!!    model%thckwk%olds(:,:,2) = model%geometry%usrf(:,:)

    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)

!WHL - debug
    if (test_parallel) then
       call glissade_test_parallel (model)
       call parallel_finalise
    endif
     
  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model, time, no_write)

    ! Perform time-step of an ice model instance with glissade dycore

    use glimmer_global, only : rk
    use glide_setup

    use glissade_temp, only: glissade_temp_driver
    use glissade_transport, only: glissade_transport_driver,  &
                                  nghost_transport, ntracer_transport

    use glimmer_paramets, only: tim0, len0, vel0
    use glimmer_physcon, only: scyr

    use glide_mask, only: glide_set_mask, calc_iareaf_iareag
    use glide_ground, only: glide_marinlim
    use glide_grid_operators
    use isostasy

    use parallel
!!    use glimmer_horiz_bcs, only: horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns, &
!!                                 horiz_bcs_unstag_scalar

    use parallel
    use glide_setup
    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(dp),  intent(in)   :: time         !*FD Current time in years

!TODO - Change no_write to write, to avoid double negatives?
    logical, optional, intent(in) :: no_write

    logical nw

    ! --- Local Variables ---

    integer :: sc  ! subcycling index

    ! temporary variables needed to reset geometry for the EVOL_NO_THICKNESS option
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_old
    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) :: stagthck_old

    ! ========================

    ! Update internal clock
    model%numerics%time = time  
    model%numerics%timecounter = model%numerics%timecounter + 1
    model%temper%newtemps = .false.

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    !TODO Not sure if this is in the right place.  G1=f(G0,T0) and T1=g(G0,T0)  If we update G1 now, then we will be doing T1=g(G1,T0).
    if (model%options%gthf > 0) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! Vertical diffusion and strain heating only; no advection
    ! ------------------------------------------------------------------------ 

    ! Note: These times have units of years.

    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

!TODO - Remove model derived type from argument list?
!HALO - Will modify glissade_temp_driver to compute over locally owned cells only.

      call t_startf('glissade_temp_driver')
       call glissade_temp_driver(model, model%options%whichtemp)
      call t_stopf('glissade_temp_driver')

       model%temper%newtemps = .true.

    end if

    ! ------------------------------------------------------------------------ 
    ! Halo updates
    ! ------------------------------------------------------------------------ 

    call parallel_halo(model%temper%bwat)    !HALO: not sure if this is needed

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
    ! MJH: This now uses velocity from the previous time step, which is appropriate for a Forward Euler time-stepping scheme

    select case(model%options%whichevol)

       case(EVOL_INC_REMAP, EVOL_NO_THICKNESS) 

       ! Use incremental remapping scheme for advecting ice thickness ---
       ! (and temperature too, if whichtemp = TEMP_REMAP_ADV)
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       ! (but not required) to use IR to advect temperature when thickness evolution is turned off.

       ! TODO  MJH If we really want to support no evolution, then we may want to implement it so that IR does not occur 
       !at all - right now a run can fail because of a CFL violation in IR even if evolution is turned off.  Do we want
       ! to support temperature evolution without thickness evolution?  If so, then the current implementation may be preferred approach.

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

!WHL - Only the new remapping scheme supported here

!WHL - Testing a new subroutine that updates all the key scalars (thck, temp, etc.) at once
!TODO - Do we need updates of lsrf, usrf, or topg?

!!       call parallel_halo_scalars(model%geometry%thck,   &
!!                                  model%temper%temp)


!HALO - Here we need halo updates for thck and temp.  

!       Velocity update may be needed if velo was not updated in halo at the end of the previous diagnostic solve
!        (just to be on the safe side).

        ! Halo updates for velocities, thickness and tracers
      call t_startf('new_remap_halo_upds')

       call staggered_parallel_halo(model%velocity%uvel)
!       call horiz_bcs_stag_vector_ew(model%velocity%uvel)

       call staggered_parallel_halo(model%velocity%vvel)
!       call horiz_bcs_stag_vector_ns(model%velocity%vvel)

       call parallel_halo(model%geometry%thck)
!       call horiz_bcs_unstag_scalar(model%geometry%thck)

       if (model%options%whichtemp == TEMP_REMAP_ADV) then
          call parallel_halo(model%temper%temp)
!          call horiz_bcs_unstag_scalar(model%temper%temp)
       endif

      call t_stopf('new_remap_halo_upds')

      call t_startf('glissade_transport_driver')

!TODO  It would be less confusing to just store the subcycling dt in a local/module variable - 
!       really only needs to be calculated once on init

       model%numerics%dt = model%numerics%dt / model%numerics%subcyc

       do sc = 1 , model%numerics%subcyc
          if (model%numerics%subcyc > 1) write(*,*) 'Subcycling transport: Cycle ',sc

          if (model%options%whichtemp == TEMP_REMAP_ADV) then  ! Use IR to transport thickness, temperature
                                                                ! (and other tracers, if present)

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            nhalo,                     ntracer_transport,         &
                                            model%velocity%uvel(:,:,:) * vel0,                    &
                                            model%velocity%vvel(:,:,:) * vel0,                    &
                                            model%geometry%thck(:,:),                             &
                                            model%temper%temp(1:model%general%upn-1,:,:) )

!TODO - Will we continue to support this option?  May not be needed.

          else  ! Use IR to transport thickness only
                ! Note: In glissade_transport_driver, the ice thickness is transported layer by layer,
                !       which is inefficient if no tracers are being transported.  (It would be more
                !       efficient to transport thickness in one layer only, using a vertically
                !       averaged velocity.)  But this option probably will not be used in practice;
                !       it is left in the code just to ensure backward compatibility with an
                !       older remapping scheme for transporting thickness only.

             call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                            model%numerics%dew * len0, model%numerics%dns * len0, &
                                            model%general%ewn,         model%general%nsn,         &
                                            model%general%upn-1,       model%numerics%sigma,      &
                                            nhalo,                     1,                         &
                                            model%velocity%uvel(:,:,:) * vel0,                &
                                            model%velocity%vvel(:,:,:) * vel0,                &
                                            model%geometry%thck(:,:))

          endif  ! whichtemp

          ! Update halos of modified fields

         call t_startf('after_remap_haloupds')

          call parallel_halo(model%geometry%thck)
!          call horiz_bcs_unstag_scalar(model%geometry%thck)

          call parallel_halo(model%temper%temp)
!          call horiz_bcs_unstag_scalar(model%temper%temp)

          ! Halo updates of other tracers, if present, would need to go here

         call t_stopf('after_remap_haloupds')

       enddo     ! subcycling

!TODO: Don't divide and multiple this variable
       model%numerics%dt = model%numerics%dt * model%numerics%subcyc

      call t_stopf('glissade_transport_driver')

      call t_stopf('inc_remap_driver')

!TODO - Would this be the appropriate place to add/remove ice at the upper and lower surfaces?
!       Should probably do this before vertical remapping.
!       (Vertical remapping currently lives in glissade_transport_driver)
!       Note that vertical remapping is needed to return to standard sigma levels,
!        as assumed by both the temperature and velocity solvers.

       if (model%options%whichevol == EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck = thck_old
          model%geomderv%stagthck = stagthck_old
       endif

    end select

!HALO - What halo updates needed here?
!       We could put the various geometry halo updates here, after thickness and temperature evolution.
!       This would include thck and tracer at a minimum.
!       Also include topg if basal topography is evolving.
!       
!       TODO: Make sure a call to calc_flwa (based on post-remap temperature) 
!             is not needed for glide_marinlim.Should be based on post-remap temperature.

    call parallel_halo(model%geometry%topg)
!    call horiz_bcs_unstag_scalar(model%geometry%topg)

    ! --- Calculate updated mask because marinlim calculation needs a mask.

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask)

!WHL - debug - print thickness field before calving
!    if (main_task) then
!       print*, ' '
!       print*, 'size(thck) =', size(model%geometry%thck,1), size(model%geometry%thck,2)
!       print*, 'Thickness before calving:'
!       do j = model%general%nsn, 1, -1
!          write(6,100) model%geometry%thck(:,j)
!       enddo 
!    endif

!100 format(34f6.2)
100 format(19f6.2)

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

!HALO - Look at marinlim more carefully and see which fields need halo updates before it is called.
!       It appears that marinlim only needs the halo of thkmask for case 5 (which was removed).  
!       If that case is removed, a thkmask halo update does not need to occur here.

    ! TODO: glide_set_mask includes a halo update of model%geometry%thkmask; move it here?
    call parallel_halo(model%geometry%thkmask) 
!    call horiz_bcs_unstag_scalar(model%geometry%thkmask)

    call parallel_halo(model%isos%relx)
!    call horiz_bcs_unstag_scalar(model%isos%relx)

!WHL - Removed old case(5), allowing for removal of some arguments

!WHL - debug - temporary case(0) to test glide_marinlim for dome problem
!!    call glide_marinlim(0, &
    call glide_marinlim(model%options%whichmarn, &
         model%geometry%thck,  &
         model%isos%relx,      &
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

!WHL - debug - print thickness field after calving
!  if (main_task) then
!    print*, ' '
!    print*, 'Thickness after calving:'
!    do j = model%general%nsn, 1, -1
!       write(6,100) model%geometry%thck(:,j)
!    enddo 
!  endif

!WHL - debug - print calving field
!  if (main_task) then
!    print*, ' '
!    print*, 'Calving:'
!    do j = model%general%nsn, 1, -1
!       write(6,100) model%climate%calving(:,j)
!    enddo 
!  endif

!HALO TODO: Test the effect of these updates with a nonzero calving field

    ! halo updates

    call parallel_halo(model%geometry%thck)    ! Updated halo values of thck are needed below in calc_lsrf
!    call horiz_bcs_unstag_scalar(model%geometry%thck)   

!WHL - debug - print thickness field after halo update
!  if (main_task) then
!    print*, ' '
!    print*, 'Thickness after halo update:'
!    do j = model%general%nsn, 1, -1
!       write(6,100) model%geometry%thck(:,j)
!    enddo 
!  endif

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

!HALO - glide_set_mask includes a halo update of model%geometry%thkmask at end of call
!       That update should be moved here if needed later (but may not be needed).

    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call
    ! call horiz_bcs_unstag_scalar(model%geometry%thkmask)

    ! --- Calculate area of ice that is floating and grounded.
!TODO This subroutine does not use iarea - remove from the call/subroutine.
!TODO May want to only calculate iarea, iareaf, iareag in glide_write_diag() and remove those calculations here.  
!     It seems hazardous to make those calculations in two different places.

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%iarea,  model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

!HALO - Need a global sum here (currently done inside calc_iareaf_iareag)

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

!TODO - Are we supporting an isostasy calculation in the parallel model?
!       While this may be a low priority in the near term, we should do so eventually.

    if (model%isos%do_isos) then
       !JEFF the isos_icewaterload() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       if (model%numerics%time >= model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model)
          model%isos%new_load = .true.
       end if
    end if
   
      ! calculate isostatic adjustment and upper and lower ice surface

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    if (model%isos%do_isos) then
       !JEFF the isos_isostasy() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       call isos_isostasy(model)
    end if

    ! ------------------------------------------------------------------------
    ! Calculate diagnostic variables, including velocity
    ! ------------------------------------------------------------------------

    call glissade_diagnostic_variable_solve(model)

    !---------------------------------------------------------------------
    ! write to netCDF file
    ! ------------------------------------------------------------------------

    if (present(no_write)) then
       nw = no_write
    else
       nw = .false.
    end if
   
    if (.not. nw) then
      call t_startf('glide_io_writeall')
       call glide_io_writeall(model,model)
      call t_stopf('glide_io_writeall')

       if (model%options%gthf > 0) then
         call t_startf('glide_lithot_io_writeall')
          call glide_lithot_io_writeall(model,model)
         call t_stopf('glide_lithot_io_writeall')
       end if

    end if

!TODO - Any halo updates needed at the end?  

  end subroutine glissade_tstep

!=======================================================================
!MJH added this diagnostic solve subroutine so it can be called from init.  

  subroutine glissade_diagnostic_variable_solve(model) 

     ! Solve diagnostic (not time-dependent) variables.  This is needed at the end of each time step once the 
     !  prognostic variables (thickness, tracers) have been updated.  
     ! It is also needed to fill out the initial state from the fields that have been read in.

    use parallel
    use glide_setup

    use glissade_temp, only: glissade_calcflwa
    use glissade_velo, only: glissade_velo_driver
    use glimmer_paramets, only: tim0, len0, vel0
    use glimmer_physcon, only: scyr
    use glide_stress, only : glide_calcstrsstr
!!    use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar, horiz_bcs_stag_scalar, horiz_bcs_stag_vector_ew, horiz_bcs_stag_vector_ns
    use glam_grid_operators,  only: glam_geometry_derivs, df_field_2d_staggered

    use glide_grid_operators, only: stagvarb    !TODO - Is this needed?  Seems redundant with df_field_2d_staggered

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

    ! Local variables

!WHL - debug
    integer :: i, j

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

    call glissade_calcflwa(model%numerics%stagsigma,    &
                           model%numerics%thklim,       &
                           model%temper%flwa,           &
                           model%temper%temp(1:model%general%upn-1,:,:),  &
                           model%geometry%thck,         &
                           model%paramets%flow_factor,  &
                           model%paramets%default_flwa, &
                           model%options%whichflwa)

    ! Halo update for flwa

    call parallel_halo(model%temper%flwa)

!WHL - Moved mask update and calving from beginning of diagnostic subroutine
!      to end of main glissade_tstep subroutine.  This ensures that for a
!      simple diagnostic case, the velocities are based on the input geometry
!      rather than a geometry that evolves due to calving.

!      Note that glide_set_mask is called near the beginning of glissade_velo_driver,
!      so that call is not needed here.

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

!WHL - uncommented these halo updates so that all the geometry updates are grouped together
!TODO - Remove the halo update of thck in glide_marinlim if it's not needed.
!       The topg update might not be needed here if it's been updated above (before the call to glide_set_mask).
!       But topg definitely needs a halo update after the isostasy calculation.

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

! MJH: Next 53 lines copied from start of glissade_tstep.
!      Notes below indicate it is unclear which of these derivatives are actually needed.  
!      But geometry derivatives are diagnostic fields and should be calculated here.

    ! ------------------------------------------------------------------------ 
    ! Now that geometry is finalized for the time step, calculate derivatives 
    ! that may be needed for the velocity solve.
    ! ------------------------------------------------------------------------     

!HALO - Make sure these geometry derivs are computed everywhere they are needed
!       (all locally owned velocity points?)

!TODO - Remove this call to geometry_derivs?  It is repeated below.
!
! SFP: Note that NOT all of the vars calculated in "geometry_derivs" are calculated in geomtry_derivs by default
! I added explicit calls to the missing ones below, along w/ appropriate halo updates prior to their call. I agree
! that we should comment this out here (and elsewhere). Better to make explicit which fields are being calculated
! and where the necessary halo updates for those fields are being made.
! An additional call to this has been commented out below.

    call glam_geometry_derivs(model)
    
    !EIB! from gc2 - think this was all replaced by geometry_derivs??
!TODO - The subroutine geometry_derivs calls subroutine stagthickness to compute stagthck.
!       Similarly for dthckdew/ns and dusrfdew/ns
!       No need to call the next three subroutines as well as geometry_derivs
!       This calculation of stagthck differs from that in geometry_derivs which calls stagthickness() 
!        in the glide_grids.F90  Which do we want to use?  
!        stagthickness() seems to be noisier but there are notes in there about some issue related to margins.

    ! SFP: not sure if these are all needed here or not. Halo updates for usrf and thck are needed in order 
    ! for periodic bcs to work. Otherwise, global halos do not contain correct values and, presumably, the gradients
    ! calculated below are incorrect in and near the global halos.
    ! Calls were added here for other staggered variables (stagusrf, stagtopg, and staglsrf), first providing halo
    ! updates to the non-stag vars, then calc. their stag values. This was done because debug lines show that these
    ! stag fields did not have the correct values in their global halos. This may be ok if they are not used at all 
    ! by the dycores called here, but I added them for consistency. More testing needed to determine if they are
    ! essential or not.

    !WHL: I am commenting out these four parallel_halo updates.
    !     Halo updates for topg and thck are done above, after the call to glide_calcflwa and before
    !      the call to geometry_derivs.  These halo updates are immediately followed by calculations
    !      of lsrf and usrf.  Since usrf and lsrf are calculated for all grid cells based on thck and topg,
    !      they do not require halo updates.
    !     Steve may want to verify that moving the halo updates to before the call to geometry_derivs
    !      does not change the results.
!!    call parallel_halo (model%geometry%usrf)
!!    call parallel_halo (model%geometry%lsrf)
!!    call parallel_halo (model%geometry%topg)
!!    call parallel_halo (model%geometry%thck)

    ! SFP: for consistency, I added these calls, so that all scalars interpolated to the stag mesh
    ! first have had their global halos updated. As w/ above calls to halo updates, these may be better 
    ! placed elsewhere. The only call originally here was the one to calc stagthck.

    !TODO - Should we replace these with calls to df_field_2d_staggered?

    call stagvarb(model%geometry%usrf, model%geomderv%stagusrf,&
                  model%general%ewn,   model%general%nsn)

    call stagvarb(model%geometry%lsrf, model%geomderv%staglsrf,&
                  model%general%ewn,   model%general%nsn)

    call stagvarb(model%geometry%topg, model%geomderv%stagtopg,&
                  model%general%ewn,   model%general%nsn)

    call stagvarb(model%geometry%thck, model%geomderv%stagthck,&    ! SFP: this call was already here. Calls to calc 
                  model%general%ewn,   model%general%nsn)           ! stagusrf, staglsrf, and stagtop were added


    call df_field_2d_staggered(model%geometry%usrf, &
                               model%numerics%dew,      model%numerics%dns, &
                               model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                               model%geometry%thck,     model%numerics%thklim )

    call df_field_2d_staggered(model%geometry%thck, &
                               model%numerics%dew,      model%numerics%dns, &
                               model%geomderv%dthckdew, model%geomderv%dthckdns, &
                               model%geometry%thck,     model%numerics%thklim )

!SFP: W.r.t WHL comment below, I went the other route above - that is, did halo updates for the non-stag
!fields first, then did the subroutine calls to calc. fields on the unstag mesh. I think this makes sure
!you are not populating the stag field global halos with bad information that may have been sitting in the 
!associated non-stag field halos in the case that you forgot to update them. Maybe?

!WHL - Changed these updates from parallel_halo to staggered_parallel_halo
!TODO - Not sure these are needed.
       !Halo updates required for inputs to glide_stress?
       call staggered_parallel_halo (model%geomderv%dusrfdew)
!       call horiz_bcs_stag_vector_ew(model%geomderv%dusrfdew)

       call staggered_parallel_halo (model%geomderv%dusrfdns)
!       call horiz_bcs_stag_vector_ns(model%geomderv%dusrfdns)

       call staggered_parallel_halo (model%geomderv%dthckdew)
!       call horiz_bcs_stag_vector_ew(model%geomderv%dthckdew)

       call staggered_parallel_halo (model%geomderv%dthckdns)
!       call horiz_bcs_stag_vector_ns(model%geomderv%dthckdns)

!TODO - It should be possible to compute stagthck without a halo update,
!       provided that thck has been updated.
!SFP: This is now done above, so I commented out the call here (that is, above we first call halo updates
! on the non-stag thickness, then calc the stagthck field. I have confirmed that stagthck field calc in this way
! has correct values in global halos.
!    call staggered_parallel_halo(model%geomderv%stagthck)
!    call horiz_bcs_stag_scalar(model%geomderv%stagthck)

    ! call parallel_halo(model%geometry%thkmask) in earlier glide_set_mask call

!###     ! basal shear stress calculations
!### 
!### !HALO - If these values are needed, it should be possible to compute them without halo updates,
!### !        provided that thck and usrf have been updated in halos.
!### !       But I think they are not needed, because calc_basal_shear is SIA-only.
!###     ! call staggered_parallel_halo(model%geomderv%stagthck) prior to calc_gline_flux
!###     ! call staggered_parallel_halo(model%geomderv%dusrfdew) prior to glide_calcstrsstr
!###     ! call staggered_parallel_halo(model%geomderv%dusrfdns) prior to glide_calcstrsstr

        ! Compute the new geometry derivatives for this time step
! TODO Merge the geom derivs with above.  Are they needed for the velocity solve?  
!      I assume so, but need to make sure they happen in the right order and have halo updates if needed.

!TODO - Would be better to have just one derivative call per field.
!       Also could consider computing derivatives in glissade_velo.F90.
!       Do halo updates as needed (e.g., thck) before computing derivatives.
!       If we need a derivative on the staggered grid (for all locally owned cells),
!        then we need one layer of halo scalars before calling the derivative routine.

!SFP: For some reason, this next call IS needed. It does not affect the results of the periodic ismip-hom test case either
! way (that is, if it is active or commented out), or the dome test case. But for some reason, if it is not active, it
! messes up both shelf test cases. There must be some important derivs being calculated within this call that are NOT
! being explicitly calculated above. 

          call glam_geometry_derivs(model)

!TODO - Pretty sure this is not needed
!SFP: this calls appear to be overwritting calls made above which explicitly make sure that 
! global halos are filled first ... comment these out?
!          call geometry_derivs_unstag(model)

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 2. Second part of diagnostic solve: 
    !    Now that geometry- and temperature-related diagnostic fields are updated, 
    !    solve velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

    ! Do not solve velocity for initial time on a restart because that breaks an exact restart.
    if ( (model%options%is_restart == 1) .and. &
         (model % numerics % time == model % numerics % tstart) ) then
       call write_log('Using uvel, vvel from restart file at initial time.')
       model%velocity%is_velocity_valid = .true.  ! TODO I don't think this flag is used anywhere, but I am setting it anyway.
    else
       ! If this is not a restart or we are not at the initial time, then proceed normally.

       if ( (model % numerics % time == model % numerics % tstart) .and. &
         ( (maxval(abs(model%velocity%uvel))/=0.0d0) .or. & 
           (maxval(abs(model%velocity%vvel))/=0.0d0) ) ) then
          ! If velocity was input and this is NOT a restart, then use the input field as the first guess at the initial time.
          ! This happens automatically, but let the user know.
          ! Using this value versus not will only change the answer within the tolerance of the nonlinear solve.  
          ! If a user already has a good guess from a previous run, they may wish to start things off with it to speed the initial solution.
          call write_log('Using uvel, vvel from input file as initial guess at initial time.  If this is not desired, please remove those fields from the input file.')
       endif

       if (main_task) then
          print *, ' '
          print *, 'Compute higher-order ice velocities, time =', model%numerics%time
       endif

       call t_startf('glissade_velo_driver')

       call glissade_velo_driver(model)

       call t_stopf('glissade_velo_driver')

!TODO - Don't think we need to update ubas, vbas, or velnorm,
!       because these can be derived directly from the 3D uvel and vvel arrays

       call staggered_parallel_halo(model%velocity%velnorm)
!       call horiz_bcs_stag_scalar(model%velocity%velnorm)
       call staggered_parallel_halo(model%velocity%ubas)
!       call horiz_bcs_stag_vector_ew(model%velocity%ubas)
       call staggered_parallel_halo(model%velocity%vbas)
!       call horiz_bcs_stag_vector_ns(model%velocity%vbas)

    endif

    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 
    ! 3. Third part of diagnostic solve: 
    ! Now that velocity is solved, calculate any diagnostic fields that are
    ! a function of velocity.
    ! ------------------------------------------------------------------------ 
    ! ------------------------------------------------------------------------ 

!TODO - Not sure if we will support calc_gline_flux.  
!       It computes a diagnostic flux, but doesn't compute it accurately.
!       Comment out for now.

    !calculate the grounding line flux after the mask is correct
    !Halo updates required for inputs to calc_gline_flux?

!!    call calc_gline_flux(model%geomderv%stagthck,  model%velocity%velnorm,  &
!!                         model%geometry%thkmask,   model%ground%gline_flux, &
!!                         model%velocity%ubas,      model%velocity%vbas,     &
!!                         model%numerics%dew)
    !Includes a halo update of model%ground%gline_flux at end of call
!TODO - Halo update of gline_flux (if needed) should go here.
!       But I think this update is not needed.

       call parallel_halo(model%stress%efvs)
!       call horiz_bcs_unstag_scalar(model%stress%efvs)

       !Tau is calculated in glide_stress and initialized in glide_types.

       call glide_calcstrsstr( model )       !*sfp* added for populating stress tensor w/ HO fields

       ! Includes halo updates of 
       ! model%stress%tau%xx, model%stress%tau%yy, model%stress%tau%xy,
       ! model%stress%tau%scalar, model%stress%tau%xz, model%stress%tau%yz

!HALO - If the stress%tau halo updates are needed, they should go here (in glissade.F90)
!       But I think they are not needed.


    ! --- A calculation of wvel could go here if we want to calculate it.
    ! --- For now, such a calculation is not needed.

  end subroutine glissade_diagnostic_variable_solve

!=======================================================================

!WHL - Moved this code back to glissade_tstep to avoid having an extra
!      subroutine call.

!!  subroutine glissade_post_tstep(model, no_write)

!!  end subroutine glissade_post_tstep

!=======================================================================

  subroutine parallel_halo_scalars(thck,     temp,   &
                                   lsrf,     usrf,   &
                                   topg,             &
                                   ntracers, tracers)

    ! Do parallel halo updates for the main scalar state variables

    use parallel

    real(dp), intent(inout), dimension(:,:) :: thck   
    real(dp), intent(inout), dimension(:,:,:) :: temp

    real(dp), intent(inout), dimension(:,:), optional ::  &
       lsrf,       & ! lower ice surface
       usrf,       & ! upper ice surface
       topg          ! basal topography

    integer, intent(in), optional :: ntracers
    real(dp), intent(inout), dimension(:,:,:,:), optional :: tracers

    integer :: nt   ! tracer index

    call parallel_halo(thck)
!    call horiz_bcs_unstag_scalar(thck)

    call parallel_halo(temp)
!!   call horiz_bcs_unstag_scalar(temp)

    ! optional updates for geometry variables

    if (present(lsrf)) then
       call parallel_halo(lsrf)
!       call horiz_bcs_unstag_scalar(lsrf)
    endif

    if (present(usrf)) then
       if (present(lsrf)) then   ! compute usrf = lsrf + thck; no halo call needed
          usrf(:,:) = lsrf(:,:) + thck(:,:)
       else
          call parallel_halo(usrf)
!          call horiz_bcs_unstag_scalar(usrf)
       endif
    endif

    if (present(topg)) then
       call parallel_halo(topg)
!       call horiz_bcs_unstag_scalar(topg)
    endif

    ! optional update for 3D tracers (e.g., ice age)

    if (present(tracers)) then

       do nt = 1, ntracers
          call parallel_halo(tracers(:,:,:,nt))
!          call horiz_bcs_unstag_scalar(tracers(:,:,:,nt))
       enddo

    endif

  end subroutine parallel_halo_scalars

!=======================================================================

    subroutine glissade_test_parallel(model)

    use parallel
    use glissade_transport, only: glissade_transport_driver,  &
                                  nghost_transport, ntracer_transport
    use glimmer_paramets, only: len0

    ! various tests of parallel model

    type(glide_global_type), intent(inout) :: model      ! model instance

    integer, dimension (:,:), allocatable ::  pgID    ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable ::  pgIDr    ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDr3    ! unique global ID for parallel runs  

    integer, dimension (:,:), allocatable ::  pgIDstagi    ! unique global ID for parallel runs  
    real(dp), dimension (:,:), allocatable ::  pgIDstagr    ! unique global ID for parallel runs  
    real(dp), dimension (:,:,:), allocatable ::  pgIDstagr3    ! unique global ID for parallel runs  

    real(dp), dimension(:,:,:), allocatable :: uvel, vvel   ! uniform velocity field

    logical, dimension(:,:), allocatable :: logvar
 
    integer :: i, j, k, n
    integer :: nx, ny, nz

    integer, parameter :: rdiag = 0    ! rank for diagnostic prints 

    real(dp), parameter :: pi = 3.14159265358979

    real(dp), parameter :: dt = 1.0       ! time step in yr
    integer, parameter  :: ntstep = 10     ! run for this number of timesteps

!    real(dp), parameter :: umag = 100.    ! uniform speed (m/yr)
    real(dp), parameter :: umag = 1000.   ! uniform speed (m/yr)

    !WHL - Tested all three of these angles (eastward, northward, and northeastward)
!    real(dp), parameter :: theta = 0.d0     ! eastward
    real(dp), parameter :: theta = pi/4.d0   ! northeastward
!    real(dp), parameter :: theta = pi/2.d0  ! northward


    real(dp) :: global_row, global_col, global_ID

    print*, ' '
    print*, 'In test_parallel, this_rank =', this_rank

    nx = model%general%ewn
    ny = model%general%nsn
    nz = model%general%upn

    allocate(logvar(nx,ny))
    allocate(pgID(nx,ny))
    allocate(pgIDr(nx,ny))
    allocate(pgIDr3(nz,nx,ny))
    allocate(pgIDstagi(nx-1,ny-1))
    allocate(pgIDstagr(nx-1,ny-1))
    allocate(pgIDstagr3(nz,nx-1,ny-1))
    allocate(uvel(nz,nx-1,ny-1), vvel(nz,nx-1,ny-1))

    if (main_task) then
       print*, ' '
       print*, 'nx, ny, nz =', nx, ny, nz
       print*, 'uhalo, lhalo =', uhalo, lhalo
       print*, 'global_ewn, global_nsn =', global_ewn, global_nsn
       print*, ' '
    endif

    print*, 'this_rank, global_row/col offset =', this_rank, global_row_offset, global_col_offset

    ! logical 2D field

    logvar(:,:) = .false.

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       logvar(i,j) = .true.
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Logical field, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,*) logvar(1:34,j)
       enddo
    endif

    call parallel_halo(logvar)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,*) logvar(:,j)
       enddo
    endif

    ! integer 2D field

    ! Compute parallel global ID for each grid cell

    pgID(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgID(i,j) = parallel_globalID_scalar(i,j,nz)    ! function in parallel_mpi.F90
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (integer), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34i5)') pgID(:,j)
       enddo
    endif

    call parallel_halo(pgID)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34i5)') pgID(:,j)
       enddo
    endif

    ! real 2D
    
    pgIDr(:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       pgIDr(i,j) = real(parallel_globalID_scalar(i,j,nz), dp)
    enddo
    enddo

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (real 2D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr(:,j)
       enddo
    endif

    call parallel_halo(pgIDr)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr(:,j)
       enddo
    endif

    ! real 3D

    pgIDr3(:,:,:) = 0

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          pgIDr3(k,i,j) = real(parallel_globalID_scalar(i,j,nz),dp) + real(k,dp)    ! function in parallel_mpi.F90
       enddo
    enddo
    enddo

    k = 1

    if (this_rank == rdiag) then
       write(6,*) ' '
       print*, 'Parallel global ID (real 3D), this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr3(k,:,j)
       enddo
    endif

    call parallel_halo(pgIDr3)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After parallel_halo_update, this_rank =', this_rank
       do j = ny, 1, -1
          write(6,'(34f6.0)') pgIDr3(k,:,j)
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
          write(6,'(33i5)') pgIDstagi(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagi)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33i5)') pgIDstagi(:,j)
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
          write(6,'(33f6.0)') pgIDstagr(:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr(:,j)
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
          write(6,'(33f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    call staggered_parallel_halo(pgIDstagr3)

    if (this_rank == rdiag) then
       write(6,*) ' '
       write(6,*) 'After staggered_parallel_halo_update, this_rank =', this_rank
       do j = ny-1, 1, -1
          write(6,'(33f6.0)') pgIDstagr3(k,:,j)
       enddo
    endif

    ! Run remapping routine

    uvel(:,:,:) = 0.d0
    vvel(:,:,:) = 0.d0

    ! Set velocity in locally owned cells

    do j = 1+lhalo, ny-uhalo
    do i = 1+lhalo, nx-uhalo
       do k = 1, nz
          uvel(k,i,j) = umag * cos(theta)
          vvel(k,i,j) = umag * sin(theta)
       enddo
    enddo
    enddo

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'Before halo update:'
          write(6,*) 'uvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') uvel(1,:,j)
          enddo
          write(6,*) ' '
          write(6,*) 'vvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') vvel(1,:,j)
          enddo
       endif

    ! staggered halo update

    call staggered_parallel_halo(uvel)
    call staggered_parallel_halo(vvel)

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'After halo update:'
          write(6,*) 'uvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') uvel(1,:,j)
          enddo
          write(6,*) ' '
          write(6,*) 'vvel, this_rank =', this_rank
          do j = ny-1, 1, -1
             write(6,'(33f6.1)') vvel(1,:,j)
          enddo
       endif

    do n = 1, ntstep

       call glissade_transport_driver(dt,                                                   &
                                      model%numerics%dew * len0, model%numerics%dns * len0, &
                                      model%general%ewn,         model%general%nsn,         &
                                      model%general%upn-1,       model%numerics%sigma,      &
                                      nghost_transport,          ntracer_transport,         &
                                      uvel(:,:,:),               vvel(:,:,:),               &
                                      model%geometry%thck(:,:),                             &
                                      model%temper%temp(1:model%general%upn-1,:,:) )

       if (this_rank == rdiag) then
          write(6,*) ' '
          write(6,*) 'New thck, n =', n
          do j = ny, 1, -1
             write(6,'(19e10.2)') model%geometry%thck(1:19,j)
          enddo
          write(6,*) ' '
          write(6,*) 'New layer 1 temp, n =', n
          do j = ny, 1, -1
             write(6,'(19f10.2)') model%temper%temp(1,1:19,j)
          enddo
       endif

    enddo  ! ntstep

    if (main_task) print*, 'Done in parallel diagnostic test'

    deallocate(logvar)
    deallocate(pgID)
    deallocate(pgIDr)
    deallocate(pgIDr3)
    deallocate(pgIDstagi)
    deallocate(pgIDstagr)
    deallocate(pgIDstagr3)
    deallocate(uvel)
    deallocate(vvel)

    end subroutine glissade_test_parallel

!=======================================================================


end module glissade
