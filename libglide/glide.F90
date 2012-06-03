!CLEANUP - glide.F90
! Moved higher-order computations to a new module, glissade.F90.
! Simplified glide.F90 to include only SIA computations.
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide.f90 - part of the Glimmer-CISM ice model           + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

!=======================================================================

module glide

  ! Driver for Glide (serial, SIA) dynamical core
  
  use glide_types
  use glide_stop
  use glide_nc_custom
  use glide_io
  use glide_lithot_io
  use glide_lithot
  use glide_profile
  use glimmer_config
  use glimmer_global

  implicit none

#ifdef GLC_DEBUG
    use glimmer_paramets, only: itest, jtest, thk0
#endif

  integer, private, parameter :: dummyunit=99

contains

!=======================================================================

!TODO - Rename to glimmer_config and put in separate module?

  subroutine glide_config(model,config,fileunit)

    ! Read glide configuration from file and print it to the log

    use glide_setup
    use isostasy
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames

    implicit none

    type(glide_global_type), intent(inout) :: model  ! model instance
    type(ConfigSection), pointer  :: config          ! structure holding sections of configuration file
    integer, intent(in), optional :: fileunit        ! fileunit for reading config file 

    type(ConfigSection), pointer :: ncconfig
    integer :: unit

    unit = 99
    if (present(fileunit)) then
       unit = fileunit
    endif

    ! read configuration file
    call glide_readconfig (model,config)
    call glide_printconfig(model)

    ! read alternate sigma levels from config file, if necessary
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

    call glimmap_readconfig(model%projection,   config,   &
                            model%numerics%dew, model%numerics%dns)

    ! netCDF I/O
    if (trim(model%funits%ncfile).eq.'') then
       ncconfig => config
    else
       call ConfigRead(process_path(model%funits%ncfile), ncconfig, unit)
    end if

    call glimmer_nc_readparams(model, ncconfig)

  end subroutine glide_config

!=======================================================================

!TODO - Move some of this to a glimmer_initialise module? Keep only the SIA stuff here.

  subroutine glide_initialise(model)

    ! Initialise Glide model instance

    use glide_setup
    use glimmer_ncio
    use glide_velo
    use glide_thck
    use glide_temp
    use glimmer_log
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use fo_upwind_advect, only : fo_upwind_advect_init

!TODO - Can we eliminate this?  Currently it's needed.
    use parallel, only: distributed_grid

    type(glide_global_type), intent(inout) :: model     ! model instance

!TODO - build glimmer_vers file or put this character elsewhere?
    character(len=100), external :: glimmer_version_char

    call write_log(trim(glimmer_version_char()))

!SCALING - Remove scaling
    ! initialise scales
    call glimmer_init_scales

!TODO - Is this call needed?
    !Initialize the NAN representation, hack to get smart compilers like gfortran to divide by zero
    call initnan 

!SCALING - Remove this call?
    ! scale parameters
    call glide_scale_params(model)

    ! set up coordinate systems

!TODO - Can we eliminate this call?  Currently get a fatal error if it's commented out.
    call distributed_grid(model%general%ewn,model%general%nsn)

    model%general%ice_grid = coordsystem_new(0.d0,               0.d0, &
                                             model%numerics%dew, model%numerics%dns, &
                                             model%general%ewn,  model%general%nsn)

    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.d0, model%numerics%dns/2.d0, &
                                              model%numerics%dew,      model%numerics%dns,      &
                                              model%general%ewn-1,     model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

!TODO - May be able to eliminate the bed softness parameter 
!       and set btrc to model%velowo%btrac_const in glide_velo
    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

    ! set uniform basal heat flux 
    model%temper%bheatflx = model%paramets%geot

    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)

    ! read first time slice
    call glide_io_readall(model,model)

    ! write projection info to log
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
       call isos_relaxed(model)
    end select

    ! open all output files
    call openall_out(model)

    ! create glide variables
    call glide_io_createall(model)

    ! initialise velocity calc
    call init_velo(model)

!TODO - Make sure this does not overwrite temperature for hotstarts
    ! initialize temperature field - this needs to happen after input file is
    ! read so we can assign artm (which could possibly be read in) if temp has not been input.
       call glide_init_temp(model)

    ! initialise thickness evolution calc
    call init_thck(model)

    if (model%options%gthf > 0) then
       call glide_lithot_io_createall(model)
       call init_lithot(model)
    end if

    ! *sfp** added for summer modeling school
    if (model%options%whichevol == EVOL_FO_UPWIND ) then
        call fo_upwind_advect_init( model%general%ewn, model%general%nsn )
    endif

!TODO - Will this module ever be supported?
    ! *mb* added; initialization of basal proc. module
    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
        model%options%which_bmod == BAS_PROC_FASTCALC) then
        
!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!                              model%numerics%ntem)
!TODO - Remove 'stop' and exit cleanly (in glide_setup?)
    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
    stop

    end if      

!TODO - Modify glide_set_mask for SIA?
!     - Mask names are perverse: 
!       model%geometry%thkmask is set in glide_mask.F90, whereas model%geometry%mask is set in glide_thckmask.F90

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

!TODO - Inline calclsrf?
    ! calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)

    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise thckwk variables; used in timeders subroutine
    model%thckwk%olds(:,:,1) = model%geometry%thck(:,:)
    model%thckwk%olds(:,:,2) = model%geometry%usrf(:,:)

    ! initialise standard glide profiling
    call glide_prof_init(model)

!TODO - Unclear on how this is used - Is it needed for serial code?
    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)
    
  end subroutine glide_initialise

!=======================================================================

  subroutine glide_tstep_p1(model,time)

    ! Perform first part of time-step of an ice model instance:
    ! temperature advection, vertical conduction, and internal dissipation.

    use glimmer_global, only : rk
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    use glide_deriv, only : df_field_2d_staggered
    use glimmer_paramets, only: tim0
    use glimmer_physcon, only: scyr
    use glide_thckmask
    use glide_grids

    type(glide_global_type), intent(inout) :: model     ! model instance
    real(rk),  intent(in)   :: time                     ! current time in years

    ! Update internal clock
    model%numerics%time = time  
    model%temper%newtemps = .false.

    model%thckwk%oldtime = model%numerics%time - (model%numerics%dt * tim0/scyr)

    call glide_prof_start(model,model%glide_prof%geomderv)

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives
    ! ------------------------------------------------------------------------     

!TODO - I suggest explicit calls to the appropriate subroutines in glide_derivs (dfdx_2d, etc.)
!       Then we would not need to use the geometry_derivs subroutine in glide_thck.
!       Currently, this subroutine computes stagthck, staglsrf, stagtopg, dusrfdew/ns, dthckdew/ns,
!        dlsrfdew/ns, d2usrfdew/ns2, d2thckdew/ns2 (2nd derivs are HO only)   

    call geometry_derivs(model)
    
    !EIB! from gc2 - think this was all replaced by geometry_derivs??
!TODO - The subroutine geometry_derivs calls subroutine stagthickness to compute stagthck.
!       Similarly for dthckdew/ns and dusrfdew/ns
!       No need to call the next three subroutines as well as geometry_derivs?

    call stagvarb(model%geometry%thck, model%geomderv%stagthck,&
                  model%general%ewn,   model%general%nsn)

    call df_field_2d_staggered(model%geometry%usrf,                              &
                               model%numerics%dew,      model%numerics%dns,      &
                               model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                               model%geometry%thck,     model%numerics%thklim )

    call df_field_2d_staggered(model%geometry%thck,                              &
                               model%numerics%dew,      model%numerics%dns,      &
                               model%geomderv%dthckdew, model%geomderv%dthckdns, &
                               model%geometry%thck,     model%numerics%thklim )

    call glide_prof_stop(model,model%glide_prof%geomderv)

    call glide_prof_start(model,model%glide_prof%ice_mask1)

    !TREY This sets local values of dom, mask, totpts, and empty
    !EIB! call veries between lanl and gc2, this is lanl version
    !magi a hack, someone explain what whichthck=5 does
    !call glide_maskthck(0, &       
    call glide_maskthck( model%geometry% thck,      &
                         model%climate%  acab,      &
                         .true.,                    &
                         model%numerics% thklim,    &
                         model%geometry% dom,       &
                         model%geometry% mask,      &
                         model%geometry% totpts,    &
                         model%geometry% empty)

    call glide_prof_stop(model,model%glide_prof%ice_mask1)

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 

    if (model%options%gthf > 0) then
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glen's A, if necessary
    ! ------------------------------------------------------------------------ 

    call glide_prof_start(model,model%glide_prof%temperature)

!SCALING - Make sure inequality makes sense with scaling removed
! I think this is OK, since these timesteps have not been scaled by tim0

    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

!TODO - Remove model derived type from argument list
       ! temperature advection, vertical conduction, and internal dissipation
       call glide_temp_driver(model, model%options%whichtemp, model%options%which_ho_diagnostic)

       model%temper%newtemps = .true.

    end if
    call glide_prof_stop(model,model%glide_prof%temperature)

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
!TODO - Remove model derived type from argument list

    call calc_btrc(model,                    &
                   model%options%whichbtrc,  &
                   model%velocity%btrc)

!TODO - Delete?
    ! ------------------------------------------------------------------------ 
    ! Calculate basal shear strength from Basal Proc module, if necessary
    ! ------------------------------------------------------------------------    
    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
        model%options%which_bmod == BAS_PROC_FASTCALC) then
!        call Basal_Proc_driver (model%general%ewn,model%general%nsn,model%general%upn,       &
!                                model%numerics%ntem,model%velocity%uvel(model%general%upn,:,:), &
!                                model%velocity%vvel(model%general%upn,:,:), &
!                                model%options%which_bmod,model%temper%bmlt,model%basalproc)
      write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
      stop
    end if

  end subroutine glide_tstep_p1

!=======================================================================

  subroutine glide_tstep_p2(model)

    ! Perform second part of time-step of an ice model instance:
    ! thickness evolution by one of several methods.

    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    use isostasy
    use fo_upwind_advect, only: fo_upwind_advect_driver
    use glide_ground, only: glide_marinlim

    type(glide_global_type), intent(inout) :: model    ! model instance

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 

    call glide_prof_start(model,model%glide_prof%ice_evo)

!TODO - Remove model derived type from argument lists

    select case(model%options%whichevol)

    case(EVOL_PSEUDO_DIFF) ! Use precalculated uflx, vflx -----------------------------------

       call thck_lin_evolve(model,model%temper%newtemps)

    case(EVOL_ADI) ! Use explicit leap frog method with uflx,vflx -------------------

       call stagleapthck(model,model%temper%newtemps)

    case(EVOL_DIFFUSION) ! Use non-linear calculation that incorporates velocity calc -----

       call thck_nonlin_evolve(model,model%temper%newtemps)

    case(EVOL_FO_UPWIND) ! Use first order upwind scheme for mass transport

       call fo_upwind_advect_driver( model )

    end select

    call glide_prof_stop(model,model%glide_prof%ice_evo)

    ! ------------------------------------------------------------------------ 
    ! get new mask
    ! Note: A call to glide_set_mask is needed before glide_marinlim.
    ! ------------------------------------------------------------------------ 

    call glide_prof_start(model,model%glide_prof%ice_mask2)

!TODO - Calculate area and vol separately from glide_set_mask?
    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    call glide_prof_stop(model,model%glide_prof%ice_mask2)

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

!TODO - Are all these arguments needed?
!       Old Glimmer code includes only arguments through model%climate%calving.

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
! The remaining arguments may not be needed?
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

!TODO - Write a better comment here.  The mask needs to be recalculated after marinlim.
    !issues with ice shelf, calling it again fixes the mask

    call glide_set_mask(model%numerics,                                &
                        model%geometry%thck,  model%geometry%topg,     &
                        model%general%ewn,    model%general%nsn,       &
                        model%climate%eus,    model%geometry%thkmask,  &
                        model%geometry%iarea, model%geometry%ivol)

    call calc_iareaf_iareag(model%numerics%dew,    model%numerics%dns,     &
                            model%geometry%iarea,  model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------

    call glide_prof_start(model,model%glide_prof%isos_water)

!TODO - Remove model derived type from argument list
    if (model%isos%do_isos) then
       if (model%numerics%time .ge. model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model)
          model%isos%new_load = .true.
       end if
    end if

    call glide_prof_stop(model,model%glide_prof%isos_water)
    
    ! ------------------------------------------------------------------------
    ! basal shear stress calculation
    ! ------------------------------------------------------------------------

    call calc_basal_shear(model%geomderv%stagthck,                          &
                          model%geomderv%dusrfdew, model%geomderv%dusrfdns, &
                          model%stress%tau_x,      model%stress%tau_y)

  end subroutine glide_tstep_p2

!=======================================================================

!TODO - no_write has been moved to tstep_postp3.  Move back?

  subroutine glide_tstep_p3(model, no_write)

    ! Perform third part of time-step of an ice model instance:
    ! calculate isostatic adjustment and upper and lower ice surface

    use isostasy
    use glide_setup

!TODO - Move back to tstep_p3?
    use glide_velo, only: gridwvel
    use glide_thck, only: timeders    

    implicit none
    type(glide_global_type), intent(inout) :: model     ! model instance
    
!TODO - Move back to tstep_p3?
    logical, optional, intent(in) :: no_write
    logical nw

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 

!TODO - Remove model derived type from argument list
    call glide_prof_start(model,model%glide_prof%isos)
    if (model%isos%do_isos) then
       call isos_isostasy(model)
    end if
    call glide_prof_stop(model,model%glide_prof%isos)

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------

!TODO - Inline calclsrf
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, &
                        model%climate%eus,   model%geometry%lsrf)

    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

!TODO - Is this the right place to increment the timecounter?  
!CESM Glimmer code has this after the netCDF write.

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glide_tstep_p3

!=======================================================================

!TODO - Move this back to tstep_p3?

  subroutine glide_tstep_postp3(model, no_write)

    use glide_setup
    use glide_velo, only: gridwvel
    use glide_thck, only: timeders

    implicit none
    type(glide_global_type), intent(inout) :: model     ! model instance

    logical, optional, intent(in) :: no_write
    logical nw

    ! For exact restart, compute wgrd here and write it to the hotstart file.
    ! (This is easier than writing thckwk quantities to the restart file.)

   call t_startf('postp3_timeders')
    call timeders(model%thckwk,            &
                  model%geometry%thck,     &
                  model%geomderv%dthckdtm, &
                  model%geometry%mask,     &
                  model%numerics%time,     &
                  1)

    call timeders(model%thckwk,            &
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
       nw=no_write
    else
       nw=.false.
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

  end subroutine glide_tstep_postp3

!=======================================================================

end module glide
