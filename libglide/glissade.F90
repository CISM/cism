!CLEANUP - glissade.F90
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

!TODO - Should some of these be moved into subroutines?
  use glide_types
  use glide_nc_custom
  use glide_io
  use glide_lithot_io
  use glide_lithot
!  use glide_profile
  use glimmer_config
  use glimmer_global

!TODO - Remove old remapping
  use glam, only : old_remapping

!TODO - Remove scaling

  implicit none

  integer, private, parameter :: dummyunit=99

contains

!=======================================================================

!TODO - Currently identical to glide_config.
!       Change to glimmer_config?

  subroutine glissade_config(model,config,fileunit)

    !*FD read glide configuration from file and print it to the log
    use glide_setup
    use isostasy
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames
    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    type(ConfigSection), pointer  :: config           ! structure holding sections of configuration file
    integer, intent(in), optional :: fileunit         !fileunit for reading config file 

    type(ConfigSection), pointer :: ncconfig
    integer :: unit

    unit = 99
    if (present(fileunit)) then
       unit = fileunit
    endif

!TODO - glimmer_readconfig and printconfig?
    ! read configuration file
    call glide_readconfig(model,config)
    call glide_printconfig(model)

    ! Read alternate sigma levels from config file, if necessary
    call glide_read_sigma(model,config)

    ! read isostasy configuration file
    call isos_readconfig(model%isos,config)
    call isos_printconfig(model%isos)
    ! read mapping from config file
    ! **** Use of dew and dns here is an ugly fudge that
    ! **** allows the use of old [GLINT projection] config section
    ! **** for backwards compatibility. It will be deleted soon.
    ! **** (You have been warned!)
    ! **** N.B. Here, dew and dns are unscaled - i.e. real distances in m

    call glimmap_readconfig(model%projection,config, &
         model%numerics%dew, &
         model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile).eq.'') then
       ncconfig => config
    else
       call ConfigRead(process_path(model%funits%ncfile), ncconfig, unit)
    end if
    call glimmer_nc_readparams(model,ncconfig)

  end subroutine glissade_config

!=======================================================================

!TODO - Modify for glissade_global_type
!       Remove extraneous use statements

  subroutine glissade_initialise(model)

    ! initialise Glissade model instance

!TODO - Are all of these needed?
    use parallel
    use glide_stop, only: register_model
    use glide_setup
    use glimmer_ncio
    use glide_velo  !TODO - Remove this
!!    use glide_thck
!!    use glide_temp
    use glissade_temp
    use glimmer_log
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground
    use glam_strs2, only : glam_velo_init

    use glide_velo_higher  !TODO - Remove this when removing calc_run_ho_diagnostic option

!!    use remap_glamutils, only : horizontal_remap_init
!!    use fo_upwind_advect, only : fo_upwind_advect_init
!!    use glam_Basal_Proc, only : Basal_Proc_init

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance

!lipscomb - TODO - build glimmer_vers file or put this character elsewhere?
    character(len=100), external :: glimmer_version_char

    call write_log(trim(glimmer_version_char()))

!SCALING - to be removed?
    ! initialise scales
    call glimmer_init_scales

!TODO - Is this call needed?
    !Initialize the NAN representation, hack to get smart compilers like gfortran to divide by zero
    call initnan 

!SCALING - to be removed?
    ! scale parameters
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

!TODO - When glide is separate from the HO driver, only glide_init_temp to be called here;
!       glissade_init_temp will be called from HO driver.

    call glissade_init_temp(model)  ! temperature lives at layer centers

!TODO - this call needed for SIA only
!!    call init_thck(model)

!TODO - Not sure backstress is ever used
!       Comment out for now; can uncomment if ever needed.
!!    call glide_initialise_backstress(model%geometry%thck,&
!!                                     model%climate%backstressmap,&
!!                                     model%climate%backstress, &
!!                                     model%climate%stressin, &
!!                                     model%climate%stressout)

    if (model%options%gthf > 0) then
       call not_parallel(__FILE__,__LINE__)
       call glide_lithot_io_createall(model)
       call init_lithot(model)
    end if

    if (model%options%which_ho_diagnostic == HO_DIAG_PP ) then

        call glam_velo_init(model%general%ewn,    model%general%nsn,  &
                            model%general%upn,                        &
                            model%numerics%dew,   model%numerics%dns, &
                            model%numerics%sigma)
    end if

!TODO - Eliminate old remapping once the new remapping is certified to be working.
!       Note: old_remapping now suppported only if call_inc_remap_driver = .true.

    if ((model%options%whichevol == EVOL_INC_REMAP ) .or. &
       (model%options%whichevol == EVOL_NO_THICKNESS)) then
      if (old_remapping) then

        if (model%options%whichtemp == TEMP_REMAP_ADV) then ! Use IR to advect temperature

           !old horizontal remapping is serialized, so it needs to be initialized with full grid size
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       global_ewn,  global_nsn,   &
                                       model%options%periodic_ew, model%options%periodic_ns, &
                                       model%general%upn,  model%numerics%sigma )

        else  ! Use IR to transport thickness only
           !old horizontal remapping is serialized, so it needs to be initialized with full grid size
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       global_ewn,  global_nsn, &
                                       model%options%periodic_ew, model%options%periodic_ns)
        endif ! whichtemp

      endif 
    endif 

!TODO - Will this module ever be supported?
    ! *mb* added; initialization of basal proc. module
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then        
!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!                              model%numerics%ntem)
!!    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
!!    stop
!!    end if      

!TODO - Change to glissade_set_mask?

!TODO - Verify exact restart.
    ! calculate mask
    if (model%options%hotstart /= 1) then  ! setting the mask destroys exact restart
       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)
    endif

!TODO- Not sure why this needs to be called here.
    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

!TODO - inline calclsrf?
    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

!TODO - Pretty sure these are SIA only
    ! initialise thckwk variables; used in timeders subroutine
!!    model%thckwk%olds(:,:,1) = model%geometry%thck(:,:)
!!    model%thckwk%olds(:,:,2) = model%geometry%usrf(:,:)

!TODO - Unclear on how this is used - Is it needed for parallel code?
    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)
    
  end subroutine glissade_initialise
  
!=======================================================================

  subroutine glissade_tstep(model,time)

    ! Perform time-step of an ice model instance with glissade dycore

    use parallel
    use glimmer_global, only : rk
    use glide_thck  !TODO - Remove this (after moving geometry_derivs)
    use glide_velo  !TODO - Remove this (make sure wvel is not needed)
    use glide_setup

    use glissade_temp, only: glissade_temp_driver
    use glissade_velo, only: glissade_velo_driver
    use glissade_transport, only: glissade_transport_driver,  &
                                  nghost_transport, ntracer_transport
    use glide_mask
    use glide_deriv, only : df_field_2d_staggered
    use glimmer_paramets, only: tim0, vis0, vis0_glam
    use glimmer_physcon, only: scyr
    use glide_thckmask  !TODO - Remove this?
    use glide_grids
    use glide_ground, only: glide_marinlim
    use stress_hom, only : glide_calcstrsstr
    use isostasy

!TODO - To be removed
    use glam, only: inc_remap_driver
    use glide_velo_higher, only: run_ho_diagnostic

    implicit none

    type(glide_global_type), intent(inout) :: model   ! model instance
    real(rk),  intent(in)   :: time         !*FD Current time in years

    integer :: sc  ! subcycling index

!TODO - These two parameters to be removed
 
    ! Temporary parameter to call inc_remap_driver in glam.F90 (or not)
    logical, parameter :: call_inc_remap_driver = .false.

    ! Temporary parameter to call run_ho_diagnostic in glide_velo_higher.F90 (or not)
    logical, parameter :: call_run_ho_diagnostic = .false.

    ! temporary variables needed to reset geometry for the EVOL_NO_THICKNESS option
    real(dp), dimension(model%general%ewn,model%general%nsn) :: thck_old
    real(dp), dimension(model%general%ewn-1,model%general%nsn-1) :: stagthck_old

    ! Update internal clock
    model%numerics%time = time  
    model%temper%newtemps = .false.

    model%thckwk%oldtime = model%numerics%time - (model%numerics%dt * tim0/scyr)

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     

!HALO - Make sure these geometry derivs are computed everywhere they are needed
!       (all locally owned velocity points?)
!TODO - I suggest explicit calls to the appropriate subroutines in glide_derivs (dfdx_2d, etc.)
!       Then we would not need to use the geometry_derivs subroutine in glide_thck.
!       Currently, this subroutine computes stagthck, staglsrf, stagtopg, dusrfdew/ns, dthckdew/ns,
!        dlsrfdew/ns, d2usrfdew/ns2, d2thckdew/ns2 (2nd derivs are HO only)   

!TODO - Remove this call. 
!        Note that there is another call to geometry_derivs in glam.F90.
!       None of the derivatives below are needed by glissade_temp_driver.
!       (Some are passed to subroutines but are never used; need to modify those subroutines.)
    call geometry_derivs(model)
    
    !EIB! from gc2 - think this was all replaced by geometry_derivs??
!TODO - The subroutine geometry_derivs calls subroutine stagthickness to compute stagthck.
!       Similarly for dthckdew/ns and dusrfdew/ns
!       No need to call the next three subroutines as well as geometry_derivs

    call stagvarb(model%geometry%thck, model%geomderv%stagthck,&
                  model%general%ewn,   model%general%nsn)

    call df_field_2d_staggered(model%geometry%usrf, &
                               model%numerics%dew,      model%numerics%dns, &
                               model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                               model%geometry%thck,     model%numerics%thklim )

    call df_field_2d_staggered(model%geometry%thck, &
                               model%numerics%dew,      model%numerics%dns, &
                               model%geomderv%dthckdew, model%geomderv%dthckdns, &
                               model%geometry%thck,     model%numerics%thklim )

    !EIB!
    
!TODO - Pretty sure that glide_maskthck is SIA only
!       totpts and empty are used only in glide_thck.
!       model%geometry%mask (not to be confused with model%geometry%thkmask) is used only in glide_thck
!       I don't see where dom is used.

    !TREY This sets local values of dom, mask, totpts, and empty
    !EIB! call veries between lanl and gc2, this is lanl version
    !magi a hack, someone explain what whichthck=5 does
    !call glide_maskthck(0, &       
!!    call glide_maskthck (model%geometry% thck,      &
!!                         model%climate%  acab,      &
!!                         .true.,                    &
!!                         model%numerics% thklim,    &
!!                         model%geometry% dom,       &
!!                         model%geometry% mask,      &
!!                         model%geometry% totpts,    &
!!                         model%geometry% empty)

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    if (model%options%gthf > 0) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! Vertical diffusion and strain heating only; no advection
    ! ------------------------------------------------------------------------ 

!SCALING - Make sure inequality makes sense with scaling removed
! I think this is OK, since these timesteps have not been scaled by tim0
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

!TODO - Remove model derived type from argument list
      call t_startf('glissade_temp_driver')
       call glissade_temp_driver(model, model%options%whichtemp)
      call t_stopf('glissade_temp_driver')

       model%temper%newtemps = .true.

    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
!TODO - SIA only (glam dycore computes beta rather than btrc)
!!    call calc_btrc(model,model%options%whichbtrc,model%velocity%btrc)

!TODO - Will this ever be supported?
    ! ------------------------------------------------------------------------ 
    ! Calculate basal shear strength from Basal Proc module, if necessary
    ! ------------------------------------------------------------------------    
!!    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
!!        model%options%which_bmod == BAS_PROC_FASTCALC) then
!        call Basal_Proc_driver (model%general%ewn,model%general%nsn,model%general%upn,       &
!                                model%numerics%ntem,model%velocity%uvel(model%general%upn,:,:), &
!                                model%velocity%vvel(model%general%upn,:,:), &
!                                model%options%which_bmod,model%temper%bmlt,model%basalproc)
!!    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
!!    stop
!!    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 

    select case(model%options%whichevol)

    case(EVOL_INC_REMAP, EVOL_NO_THICKNESS) 

       ! Use incremental remapping scheme for advecting ice thickness ---
       ! (and temperature too, if whichtemp = TEMP_REMAP_ADV)
       ! MJH: I put the no thickness evolution option here so that it is still possible 
       ! (but not required) to use IR to advect temperature when thickness evolution is turned off.

       if (model%options%whichevol .eq. EVOL_NO_THICKNESS) then
          ! store old thickness
          thck_old = model%geometry%thck
          stagthck_old = model%geomderv%stagthck
       endif

!WHL - Moved inc_remap_driver code from glam.F90 to this level
!TODO - Remove glam.F90.

      call t_startf('inc_remap_driver')

       if (call_inc_remap_driver) then

          call inc_remap_driver(model)   ! in glam.F90 - to be removed

       else   ! use inlined code from inc_remap_driver

        ! Compute the new geometry derivatives for this time step

!HALO - Would be better to have one derivative call per field.
!       Also could consider computing derivatives in glam_strs2.F90, so as
!        not to have to pass them through the interface.
!       Do halo updates as needed (e.g., thck) before computing derivatives.
!       If we need a derivative on the staggered grid (for all locally owned cells),
!        then we need one layer of halo scalars before calling the derivative routine.

          call geometry_derivs(model)

!HALO - Pretty sure this is not needed
          call geometry_derivs_unstag(model)

          if (main_task) then
             print *, ' '
             print *, 'Compute higher-order ice velocities, time =', model%numerics%time
          endif

         call t_startf('ho_velo_diagnostic')

!TODO - Remove run_ho_diagnostic; it has been superseded by glissade_velo_driver.

          if (call_run_ho_diagnostic) then

             call run_ho_diagnostic(model)   ! in glide_velo_higher.F90

          else

!TODO - Replace model derived type with explicit arguments?
             call glissade_velo_driver(model)

          endif

         call t_stopf('ho_velo_diagnostic')

          if (main_task) then
             print *, ' '
             print *, 'Compute dH/dt'
          endif

!WHL - Only the new remapping scheme supported here

          if( model%numerics%tend > model%numerics%tstart) then

!HALO - Need halo updates here for thck, temp (and any other advected tracers), uvel and vvel.
!       If nhalo >= 2, then no halo updates should be needed inside glissade_transport_driver.

!PW FOLLOWING NECESSARY?
!HALO - These halo updates could be moved up a level to the new glissade driver.

           ! Halo updates for velocities, thickness and tracers
            call t_startf('new_remap_halo_upds')
             call staggered_parallel_halo(model%velocity%uvel)
             call staggered_parallel_halo(model%velocity%vvel)
             call parallel_halo(model%geometry%thck)
             if (model%options%whichtemp == TEMP_REMAP_ADV) then
                !If advecting other tracers, add parallel_halo update here
                call parallel_halo(model%temper%temp)
             endif
            call t_stopf('new_remap_halo_upds')

            call t_startf('glissade_transport_driver')

             model%numerics%dt = model%numerics%dt / model%numerics%subcyc

             do sc = 1 , model%numerics%subcyc

                if (model%options%whichtemp == TEMP_REMAP_ADV) then  ! Use IR to transport thickness, temperature
                                                                     ! (and other tracers, if present)

                   call glissade_transport_driver(model%numerics%dt * tim0,                             &
                                                  model%numerics%dew * len0, model%numerics%dns * len0, &
                                                  model%general%ewn,         model%general%nsn,         &
                                                  model%general%upn-1,       model%numerics%sigma,      &
                                                  nghost_transport,          ntracer_transport,         &
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
                                                  nghost_transport,          1,                         &
                                                  model%velocity%uvel(:,:,:) * vel0,                &
                                                  model%velocity%vvel(:,:,:) * vel0,                &
                                                  model%geometry%thck(:,:))

                endif  ! whichtemp

             enddo     ! subcycling

             model%numerics%dt = model%numerics%dt * model%numerics%subcyc

            call t_stopf('glissade_transport_driver')

          endif   ! tend > tstart

       endif      ! call inc_remap_driver

      call t_stopf('inc_remap_driver')

!TODO - Would this be the appropriate place to add/remove ice at the upper and lower surfaces?
!       Should probably do this before vertical remapping.
!       (Vertical remapping currently lives in glissade_transport_driver)
!       Note that vertical remapping is needed to return to standard sigma levels,
!        as assumed by both the temperature and velocity solvers.

        !Update halos of modified fields
      call t_startf('after_remap_haloupds')

!HALO - Move these updates to the new glissade driver.

       call parallel_halo(model%geometry%thck)
       call parallel_halo(model%temper%temp)
      call t_stopf('after_remap_haloupds')

!HALO - I think that these are not needed, provided that glide_stress loops over locally owned cells only.
       !Halo updates required for inputs to glide_stress?
       call staggered_parallel_halo(model%geomderv%dusrfdew)
       call staggered_parallel_halo(model%geomderv%dusrfdns)
       call staggered_parallel_halo(model%geomderv%dthckdew)
       call staggered_parallel_halo(model%geomderv%dthckdns)

!HALO - Pretty sure these can be removed
       ! call parallel_halo(model%geometry%thck) in inc_remap_driver
       ! call staggered_parallel_halo(model%velocity%uvel) in inc_remap_driver
       ! call staggered_parallel_halo(model%velocity%vvel) in inc_remap_driver

!HALO - I think this update is not needed, provided that glide_calcstrsstr loops over locally owned cells.
       call parallel_halo(model%stress%efvs)

       !Tau is calculated in glide_stress and initialized in glide_types.
       call glide_calcstrsstr( model )       !*sfp* added for populating stress tensor w/ HO fields
       !Includes halo updates of 
       ! model%stress%tau%xx, model%stress%tau%yy, model%stress%tau%xy,
       ! model%stress%tau%scalar, model%stress%tau%xz, model%stress%tau%yz
       !at end of call
!HALO - If the stress%tau halo updates are needed, they should go here (in glissade.F90)
!       But I think they are not needed.

       if (model%options%whichevol .eq. EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck = thck_old
          model%geomderv%stagthck = stagthck_old
       endif

end select

!HALO - I suggest putting the various geometry halo updates here, after thickness and temperature evolution.
!       This would include thck at a minimum.
!       Also include topg if basal topography is evolving.

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------

    !Halo updates required for inputs to glide_set_mask?
    ! call parallel_halo(model%geometry%thck) in inc_remap_driver
    call parallel_halo(model%geometry%topg)

!TODO - May want to write a new subroutine, glissade_set_mask, that loops over locally owned cells
!        and is followed by a halo update (for thkmask) in the driver.
!       For the serial SIA code, we probably shouldn't change the loops.
 
!Note: The call to glide_set_mask is needed before glide_marinlim.

       call glide_set_mask(model%numerics,                                &
                           model%geometry%thck,  model%geometry%topg,     &
                           model%general%ewn,    model%general%nsn,       &
                           model%climate%eus,    model%geometry%thkmask,  &
                           model%geometry%iarea, model%geometry%ivol)

    !Includes a halo update of model%geometry%thkmask at end of call
!HALO - Halo update of thkmask should go here

!HALO - It should be possible to compute stagthck without a halo update,
!       provided that thck has been updated.

    call staggered_parallel_halo(model%geomderv%stagthck)
    ! call parallel_halo(model%geometry%thkmask) in earlier glide_set_mask call

!TODO - Not sure we need to update ubas, vbas, or surfvel,
!       because these are already part of the 3D uvel and vvel arrays
    call staggered_parallel_halo(model%velocity%surfvel)
    call staggered_parallel_halo(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)

!TODO - Not sure if we will support calc_gline_flux.  
!       It computes a diagnostic flux, but doesn't compute it accurately.
!       Comment out for now.

    !calculate the grounding line flux after the mask is correct
    !Halo updates required for inputs to calc_gline_flux?
!!    call calc_gline_flux(model%geomderv%stagthck,model%velocity%surfvel, &
!!                         model%geometry%thkmask,model%ground%gline_flux, model%velocity%ubas, &
!!                         model%velocity%vbas, model%numerics%dew)
    !Includes a halo update of model%ground%gline_flux at end of call
!HALO - Halo update of gline_flux (if needed) should go here.

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

!HALO - Look at marinlim more carefully and see which fields need halo updates before it is called.

    call parallel_halo(model%isos%relx)

!HALO - not sure if needed for glide_marinlim.
    call parallel_halo(model%temper%flwa)

!HALO - Not sure backstress is ever used
    call parallel_halo(model%climate%backstress)
    ! call parallel_halo(model%geometry%usrf) not actually used

!TODO - glissade_marinlim?
    call glide_marinlim(model%options%whichmarn, &
         model%geometry%thck,      &
         model%isos%relx,      &
         model%geometry%topg,   &
         model%temper%flwa,   &
         model%numerics%sigma,   &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving,  &
!TODO - The remaining arguments may not be needed
         model%climate%backstress, &
         model%climate%tempanmly, &
         model%numerics%dew,    &
         model%numerics%dns, &
         model%climate%backstressmap, &
         model%climate%stressout, &
         model%climate%stressin, &
         model%ground, &
         model%general%nsn, &
         model%general%ewn)
         ! model%geometry%usrf) not used in routine

    !Includes halo updates of model%geometry%thck and model%climate%calving
    ! and of model%climate%backstress), when needed
!HALO - Those updates should go here instead.
!       Note: thck will have changed, and updated halo values are needed below in calc_lsrf.

    !Note that halo updates for model%geometry%thkmask not needed in current implementation of case 3
    ! model%climate%eus needed only for disabled case 6
    ! model%ground components needed only for disabled case 6

!TODO - Write a better comment here.  The mask needs to be recalculated after marinlim.
    !issues with ice shelf, calling it again fixes the mask

!HALO - Halo updates are not needed here if glide_set_mask loops over locally owned cells.

!TODO - glissade_set_mask?
    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

!HALO - glide_set_mask includes a halo update of model%geometry%thkmask at end of call
!       That update should be moved here if needed later (but may not be needed).

    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%iarea,  model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

!HALO - Need a global sum here (currently done inside calc_iareaf_iareag)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------

!TODO - Are we supporting an isostasy calculation in the parallel model?
!       While this may be a low priority in the near term, we should do so eventually.

    if (model%isos%do_isos) then
       !JEFF the isos_icewaterload() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       if (model%numerics%time.ge.model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model)
          model%isos%new_load = .true.
       end if
    end if
    
    ! basal shear stress calculations

!HALO - If these values are needed, it should be possible to compute them without halo updates,
!        provided that thck and usrf have been updated in halos.
!       But I think they are not needed, because calc_basal_shear is SIA-only.
    ! call staggered_parallel_halo(model%geomderv%stagthck) prior to calc_gline_flux
    ! call staggered_parallel_halo(model%geomderv%dusrfdew) prior to glide_calcstrsstr
    ! call staggered_parallel_halo(model%geomderv%dusrfdns) prior to glide_calcstrsstr

!HALO - Does the last part of the time step require any halo info?

      ! calculate isostatic adjustment and upper and lower ice surface

!TODO - Need to support an isostasy calculation in the parallel model?

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

    if (model%isos%do_isos) then
       !JEFF the isos_isostasy() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       call isos_isostasy(model)
    end if

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------

! Calculate time-derivatives of thickness and upper surface elevation ------------
    !Halo updates required for inputs to glide_calcsrf?
    ! call parallel_halo(model%geometry%thck) within glide_marinlim
    ! call parallel_halo(model%geometry%topg) before previous call to glide_set_mask

!HALO - Verify that both thck and topg are up to date before this call.
!       Note that glide_calclsrf loops over all cells (not just locally owned)

!TODO - Inline calclsrf
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    !If input halos are up to date, halo update for model%geometry%lsrf should not be necessary

    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)
    !If input halos are up to date, halo update for model%geometry%usrf should not be necessary

!TODO - Is this the right place to increment the timecounter?  
!CESM Glimmer code has this after the netCDF write.

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glissade_tstep

!=======================================================================

!TODO - Merge with glissade_tstep or move to simple_driver?
!       Not sure these calls are needed for HO solver.

  subroutine glissade_post_tstep(model, no_write)

    !* This routine does the parallel routines and output that was in _p3()
    !* _p3() is executed in serial on main node.  These are parallel operations.
    !* Jeff Nichols, created for Bill Lipscomb September 2011

    use parallel
    use glide_setup

!TODO - not needed?
    use glide_velo, only: gridwvel
    use glide_thck, only: timeders

    implicit none
    type(glide_global_type), intent(inout) :: model   ! model instance

    logical, optional, intent(in) :: no_write
    logical nw

!TODO - Determine correct location for these calls (related to exact restart of serial model).
!       Remove if not needed here.

    ! These three calls are in glide_temp in the full-temperature section replaced by Bill's new temperature code.
    ! Commenting out until I hear otherwise.
   call t_startf('postp3_timeders')
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
   call t_stopf('postp3_timeders')

    ! Calculate the vertical velocity of the grid ------------------------------------

   call t_startf('postp3_gridwvel')
    call gridwvel(model%numerics%sigma,  &
            model%numerics%thklim, &
            model%velocity%uvel,   &
            model%velocity%vvel,   &
            model%geomderv,        &
            model%geometry%thck,   &
            model%velocity%wgrd)
   call t_stopf('postp3_gridwvel')

    !---------------------------------------------------------------------
    ! write to netCDF file
    ! ------------------------------------------------------------------------

    if (present(no_write)) then
       nw = no_write
    else
       nw = .false.
    end if
   
!TODO - Change to glimmer_io_writeall?
!       Move to simple_driver?
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

  end subroutine glissade_post_tstep

!=======================================================================

end module glissade
