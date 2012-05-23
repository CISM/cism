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

!TODO - Make a new model called glissade.F90 patterned after glide.F90.
!       Remove HO and parallel stuff from glide.F90, remove SIA-only stuff from glissade.F90.

module glide
  !*FD the top-level GLIDE module

  use glide_types
  use glide_stop
  use glide_nc_custom
  use glide_io
  use glide_lithot_io
  use glide_lithot
  use glide_profile
  use glimmer_config
  use glimmer_global
  use glam, only : old_remapping

#ifdef GLC_DEBUG
    use glimmer_paramets, only: itest, jtest, thk0
#endif

  integer, private, parameter :: dummyunit=99

contains

  subroutine glide_config(model,config,fileunit)

    !*FD read glide configuration from file and print it to the log
    use glide_setup
    use isostasy
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames
    implicit none

    type(glide_global_type) :: model          !*FD model instance
    type(ConfigSection), pointer  :: config   !*FD structure holding sections of configuration file
    integer, intent(in), optional :: fileunit !*FD fileunit for reading config file 

    type(ConfigSection), pointer :: ncconfig
    integer :: unit

    unit = 99
    if (present(fileunit)) then
       unit = fileunit
    endif

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

  end subroutine glide_config

  subroutine glide_initialise(model)
    !*FD initialise GLIDE model instance
    use parallel
    use glide_setup
    use glimmer_ncio
    use glide_velo
    use glide_velo_higher
    use glide_thck
    use glide_temp
    use glissade_temp
    use glimmer_log
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    use glide_ground
    use glam_strs2, only : glam_velo_fordsiapstr_init
    use remap_glamutils, only : horizontal_remap_init
    use fo_upwind_advect, only : fo_upwind_advect_init
    use glam_Basal_Proc, only : Basal_Proc_init

    implicit none
    type(glide_global_type) :: model        !*FD model instance

!lipscomb - TO DO - build glimmer_vers file or put this character elsewhere?
    character(len=100), external :: glimmer_version_char

#ifdef GLC_DEBUG
    integer :: i, j, k
#endif

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales
    call initnan !Initialize the NAN representation, hack to get smart compilers like gfortran to divide by zero

    ! scale parameters
    call glide_scale_params(model)
    ! set up coordinate systems
    ! time to change to the parallel values of ewn and nsn
    call distributed_grid(model%general%ewn,model%general%nsn)
    model%general%ice_grid = coordsystem_new(0.d0, 0.d0, &
         model%numerics%dew, model%numerics%dns, &
         model%general%ewn, model%general%nsn)
    model%general%velo_grid = coordsystem_new(model%numerics%dew/2.,model%numerics%dns/2., &
         model%numerics%dew,model%numerics%dns, &
         model%general%ewn-1,model%general%nsn-1)

    ! allocate arrays
    call glide_allocarr(model)

!TODO - May be able to eliminate the bed softness parameter 
!       and set btrc to model%velowo%btrac_const in glide_velo
    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

!TODO - Can remove; PBJ only
    !Initialize boundary condition fields to be NaN everywhere 
    model%geometry%marine_bc_normal = NaN
    
    ! set uniform basal heat flux 
    model%temper%bheatflx = model%paramets%geot

    ! load sigma file
    call glide_load_sigma(model,dummyunit)

    ! open all input files
    call openall_in(model)
    ! and read first time slice
    call glide_io_readall(model,model)
    ! Write projection info to log
    call glimmap_printproj(model%projection)

#ifdef GLC_DEBUG
       write(6,*) 'Opened input files'
       write(6,*) 'i, j, thck, thck(m):', itest, jtest, &
               model%geometry%thck(itest,jtest), model%geometry%thck(itest,jtest)*thk0
#endif

    ! read lithot if required
    if (model%options%gthf.gt.0) then
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
    ! create glide variables
    call glide_io_createall(model)

    ! initialise glide components

!TODO - Most of what's done in init_velo is needed for SIA only
!       Can remove this call provided velowk is not used elsewhere (e.g., to call wvelintg)
       
    call init_velo(model)

!TODO - When glide is separate from the HO driver, only glide_init_temp to be called here;
!       glissade_init_temp will be called from HO driver.

    ! MJH: Initialize temperature field - this needs to happen after input file is
    ! read so we can assign artm (which could possibly be read in) if temp has not been input.
    if (model%options%whichtemp == TEMP_GLIMMER) then
       call glide_init_temp(model)  ! temperature lives at layer intefaces
    else
       call glissade_init_temp(model)  ! temperature lives at layer centers
    endif 

!TODO - this call needed for SIA only
    call init_thck(model)

!TODO - Not sure backstress is ever used
!       Probably safe to comment out the call; can uncomment if ever needed.
    call glide_initialise_backstress(model%geometry%thck,&
                                     model%climate%backstressmap,&
                                     model%climate%backstress, &
                                     model%climate%stressin, &
                                     model%climate%stressout)

    if (model%options%gthf.gt.0) then
       call not_parallel(__FILE__,__LINE__)
       call glide_lithot_io_createall(model)
       call init_lithot(model)
    end if

!TODO - Move the HO code below to the HO driver, leaving only the old Glide code

    if (model%options%which_ho_diagnostic == HO_DIAG_PP ) then

        call glam_velo_fordsiapstr_init(model%general%ewn,    model%general%nsn,  &
                                        model%general%upn,                        &
                                        model%numerics%dew,   model%numerics%dns, &
                                        model%numerics%sigma)
    end if

!TODO - Eliminate old remapping once the new remapping is certified to be working.

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

!TODO - Eliminate this option once the parallel upwind transport is certified to be working.
    ! *sfp** added for summer modeling school
    if (model%options%whichevol== EVOL_FO_UPWIND ) then

        ! JEFF - Review for distributed. OK for parallel Trilinos. call not_parallel(__FILE__,__LINE__)
        call fo_upwind_advect_init( model%general%ewn, model%general%nsn )

    endif

!TODO - Will this module ever be supported?
    ! *mb* added; initialization of basal proc. module
    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
        model%options%which_bmod == BAS_PROC_FASTCALC) then
        
!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!                              model%numerics%ntem)
    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
    stop

    end if      

!TODO - Compute and advect the ice age.
!       Initialization should take place in glide_types
    ! initialise ice age
    !! This is a placeholder; currently the ice age is not computed.  
    !! Currently the ice age is only computed for remapping transport
    !! (whichevol = 3 or 4)
    model%geometry%age(:,:,:) = 0._dp

!TODO - Can remove these lines
!MJH moved flwa init to glide_init_temp/glissade_init_temp
!    if (model%options%hotstart.ne.1) then
!       ! initialise Glen's flow parameter A using an isothermal temperature distribution
!       call glide_temp_driver(model,0,0)
!    end if       

!TODO - Decide how to modify glide_set_mask.  Might be able to remove some masks.
!       Also, may want to rename to cism_set_mask

!TODO - Mask names are perverse: 
!       model%geometry%thkmask is set in glide_mask.F90, whereas model%geometry%mask is set in glide_thckmask.F90

    ! calculate mask
    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)

    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

!TODO - Remove this call?  Steve thinks it is PBJ only.
!       The variable marine_bc_normal is never used.
 
    !calculate the normal at the marine margin
    call glide_marine_margin_normal(model%geometry%thck, model%geometry%thkmask, model%geometry%marine_bc_normal)

!TODO - Not sure why this was commented out, but OK to remove since it is called just above.
!     - Still have to verify exact restart

    ! calculate mask
    !if (model%options%hotstart.ne.1) then  ! setting the mask destroys exact restart
    !   call glide_set_mask(model)
    !end if

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

!TODO - Pretty sure these are SIA only
    ! initialise thckwk variables; used in timeders subroutine
    model%thckwk%olds(:,:,1) = model%geometry%thck(:,:)
    model%thckwk%olds(:,:,2) = model%geometry%usrf(:,:)

!TODO - SIA only?
    ! initialise standard glide profiling
    call glide_prof_init(model)

!TODO - Unclear on how this is used - Is it needed for parallel code?
    ! register the newly created model so that it can be finalised in the case
    ! of an error without needing to pass the whole thing around to every
    ! function that might cause an error
    call register_model(model)
    
#ifdef GLC_DEBUG
       write(6,*) ' '
       write(6,*) 'End of glide_init'
       i = itest
       j = jtest
       write(6,*) 'i, j, thck =', i, j, model%geometry%thck(i,j)
       write(6,300) k, model%geometry%thck(i,j)
       write(6,*) ' '
       write(6,*) 'k, temperature'
       do k = 1, model%general%upn
            write(6,300) k, model%temper%temp(k,i,j)
       enddo
  300  format(i3, Z24.20)
#endif

  end subroutine glide_initialise
  
!TODO - glide_tstep_p1, p2, and p3 will be combined in glissade.F90.
!       Could be part of a single glissade_tstep subroutine

  subroutine glide_tstep_p1(model,time)
    !*FD Performs first part of time-step of an ice model instance.
    !*FD calculate velocity and temperature
    use parallel
    use glimmer_global, only : rk
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glissade_temp
    use glide_mask
    use glide_deriv, only : df_field_2d_staggered
    use glimmer_paramets, only: tim0
    use glimmer_physcon, only: scyr
    use glide_thckmask
    use glide_grids
    use glam_Basal_Proc, only : Basal_Proc_driver

    implicit none

    type(glide_global_type) :: model        !*FD model instance
    real(rk),  intent(in)   :: time         !*FD Current time in years

    ! Update internal clock
    model%numerics%time = time  
    model%temper%newtemps = .false.

    model%thckwk%oldtime = model%numerics%time - (model%numerics%dt * tim0/scyr)

#ifdef GLC_DEBUG
       write(6,*) ' '
       write(6,*) 'time =', model%numerics%time
       write(6,*) 'tinc =', model%numerics%tinc
       write(6,*) 'oldtime =', model%thckwk%oldtime
#endif

    call glide_prof_start(model,model%glide_prof%geomderv)

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     

!HALO - If these geometry derivs are computed only on locally owned nodes,
!       then some derivs may need halo updates.
!TODO - I suggest explicit calls to the appropriate subroutines in glide_derivs (dfdx_2d, etc.)
!       Then we would not need to use the geometry_derivs subroutine in glide_thck.

    call geometry_derivs(model)
    
    !EIB! from gc2 - think this was all replaced by geometry_derivs??
    call stagvarb(model%geometry% thck, &
         model%geomderv% stagthck,&
         model%general%  ewn, &
         model%general%  nsn)

    call df_field_2d_staggered(model%geometry%usrf, &
         model%numerics%dew, model%numerics%dns, &
         model%geomderv%dusrfdew, & 
         model%geomderv%dusrfdns, &
        model%geometry%thck,     &
        model%numerics%thklim )

    call df_field_2d_staggered(model%geometry%thck, &
         model%numerics%dew, model%numerics%dns, &
         model%geomderv%dthckdew, & 
         model%geomderv%dthckdns, &
        model%geometry%thck,     &
        model%numerics%thklim )

    !EIB!
    
    call glide_prof_stop(model,model%glide_prof%geomderv)

    call glide_prof_start(model,model%glide_prof%ice_mask1)

!TODO - Pretty sure that glide_maskthck is SIA only
!       totpts and empty are used only in glide_thck.
!       model%geometry%mask (not to be confused with model%geometry%thkmask) is used only in glide_thck
!       I don't see where dom is used.

    !TREY This sets local values of dom, mask, totpts, and empty
    !EIB! call veries between lanl and gc2, this is lanl version
    !magi a hack, someone explain what whichthck=5 does
    !call glide_maskthck(0, &       
    call glide_maskthck( &       
         model%geometry% thck,      &
         model%climate%  acab,      &
         .true.,                    &
         model%numerics%thklim,     &
         model%geometry% dom,       &
         model%geometry% mask,      &
         model%geometry% totpts,    &
         model%geometry% empty)

    call glide_prof_stop(model,model%glide_prof%ice_mask1)

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    if (model%options%gthf.gt.0) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

!TODO - Can have one call in SIA driver, the other call in the HO driver

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
    call glide_prof_start(model,model%glide_prof%temperature)
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

       if (model%options%whichtemp == TEMP_GLIMMER) then 

         ! Glimmer temp driver, including temperature advection

        call t_startf('glide_temp_driver')
         call glide_temp_driver(model, model%options%whichtemp, model%options%which_ho_diagnostic)
        call t_stopf('glide_temp_driver')

       else  ! either TEMP_REMAP_ADV, or a simple option like TEMP_SURFACE_AIR_TEMP

         ! Vert diffusion and strain heating only; no advection
         ! Remapping routine is used to advect temperature in glide_tstep_p2

        call t_startf('glissade_temp_driver')
         call glissade_temp_driver(model, model%options%whichtemp)
        call t_stopf('glissade_temp_driver')

       endif

       model%temper%newtemps = .true.

    end if
    call glide_prof_stop(model,model%glide_prof%temperature)

!TODO - SIA only (glam dycore computes beta rather than btrc)

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
    call calc_btrc(model,model%options%whichbtrc,model%velocity%btrc)

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

!----------------------------------------------------------------------------- 

  subroutine glide_tstep_p2(model)
    !*FD Performs second part of time-step of an ice model instance.
    !*FD write data and move ice
    use parallel
    use glide_thck
    use glide_velo
    use glide_ground
    use glide_setup
    use glide_temp
    use glide_mask
    use isostasy

    ! driver module/subroutines for Payne/Price HO dynamics and LANL inc. remapping for dH/dt 
    use glam, only: inc_remap_driver
    use fo_upwind_advect, only: fo_upwind_advect_driver
    use stress_hom, only : glide_calcstrsstr

    implicit none

    type(glide_global_type) :: model        !*FD model instance

    ! temporary variables needed to reset geometry for the EVOL_NO_THICKNESS option
    real (kind = dp), dimension(model%general%ewn,model%general%nsn) :: thck_old
    real (kind = dp), dimension(model%general%ewn-1,model%general%nsn-1) :: stagthck_old

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
    call glide_prof_start(model,model%glide_prof%ice_evo)

!TODO - The first 3 cases are SIA only.
!       For the glissade driver we will support only the INC_REMAP and EVOL_NO_THICKNESS options 

    select case(model%options%whichevol)
    case(EVOL_PSEUDO_DIFF) ! Use precalculated uflx, vflx -----------------------------------

       call not_parallel(__FILE__,__LINE__)
       call thck_lin_evolve(model,model%temper%newtemps)

    case(EVOL_ADI) ! Use explicit leap frog method with uflx,vflx -------------------

       call not_parallel(__FILE__,__LINE__)
       call stagleapthck(model,model%temper%newtemps)

    case(EVOL_DIFFUSION) ! Use non-linear calculation that incorporates velocity calc -----

       call not_parallel(__FILE__,__LINE__)
       call thck_nonlin_evolve(model,model%temper%newtemps)

    case(EVOL_INC_REMAP, EVOL_NO_THICKNESS) ! Use incremental remapping scheme for advecting ice thickness ---
                         ! (and temperature too, if whichtemp = TEMP_REMAP_ADV)
                         ! MJH: I put the no thickness evolution option here so that it is still possible 
                         ! (but not required) to use IR to advect temperature when thickness evolution is turned off.

       if (model%options%whichevol .eq. EVOL_NO_THICKNESS) then
          ! store old thickness
          thck_old = model%geometry%thck
          stagthck_old = model%geomderv%stagthck
       endif

!TODO - Suggest moving inc_remap_driver code to this level, to separate velocity solve from thck/temp evolution
      call t_startf('inc_remap_driver')
       call inc_remap_driver( model )
      call t_stopf('inc_remap_driver')

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

!TODO - Remove this call once the parallel upwind scheme has been tested.

    case(EVOL_FO_UPWIND) ! Use first order upwind scheme for mass transport
       call not_parallel(__FILE__,__LINE__)
       call fo_upwind_advect_driver( model )

end select

!HALO - I suggest putting the various geometry halo updates here, after thickness and temperature evolution.
!       This would include thck at a minimum.
!       Also include topg if basal topography is evolving.

    call glide_prof_stop(model,model%glide_prof%ice_evo)

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------
    call glide_prof_start(model,model%glide_prof%ice_mask2)

    !Halo updates required for inputs to glide_set_mask?
    ! call parallel_halo(model%geometry%thck) in inc_remap_driver
    call parallel_halo(model%geometry%topg)

!TODO - May want to write a new subroutine, glissade_set_mask, that loops over locally owned cells
!        and is followed by a halo update (for thkmask) in the driver.
!       For the serial SIA code, we probably shouldn't change the loops.
 
!Note: The call to glide_set_mask is needed before glide_marinlim.

    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)
    !Includes a halo update of model%geometry%thkmask at end of call
!HALO - Halo update of thkmask should go here

    call glide_prof_stop(model,model%glide_prof%ice_mask2)

!TODO - This call was used by the PBJ core and probably is not needed.
    call glide_marine_margin_normal(model%geometry%thck, model%geometry%thkmask, model%geometry%marine_bc_normal)
    !Includes a halo update of model%geometry%marine_bc_normal at end of call

!TODO - Not sure if we will support calc_gline_flux.  
!       It computes a diagnostic flux, but doesn't compute it accurately.
!       Comment out for now?

    !Halo updates required for inputs to calc_gline_flux?
!HALO - It should be possible to compute stagthck without a halo update,
!       provided that thck has been updated.

    call staggered_parallel_halo(model%geomderv%stagthck)
    ! call parallel_halo(model%geometry%thkmask) in earlier glide_set_mask call

!TODO - Not sure we need to update ubas, vbas, or surfvel,
!       because these are already part of the 3D uvel and vvel arrays
    call staggered_parallel_halo(model%velocity%surfvel)
    call staggered_parallel_halo(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)
    ! call parallel_halo(model%ground%gline_flux) - halo update not required for correctness

!TODO - Not sure this info is currently used.  Comment out the call to calc_gline_flux for now?
    !calculate the grounding line flux after the mask is correct
    call calc_gline_flux(model%geomderv%stagthck,model%velocity%surfvel, &
                         model%geometry%thkmask,model%ground%gline_flux, model%velocity%ubas, &
                         model%velocity%vbas, model%numerics%dew)
    !Includes a halo update of model%ground%gline_flux at end of call
!HALO - Halo update of gline_flux should go here.

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

!HALO - Look at marinlim more carefully and see which fields need halo updates before it is called.

!HALO - Not sure we need halo values for isostasy
    call parallel_halo(model%isos%relx)

!HALO - not sure if needed for glide_marinlim.
    call parallel_halo(model%temper%flwa)

!HALO - Not sure backstress is ever used
    call parallel_halo(model%climate%backstress)

    ! call parallel_halo(model%geometry%usrf) not actually used

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

    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)

!HALO - glide_set_mask includes a halo update of model%geometry%thkmask at end of call
!       That update should be moved here if needed later (but may not be needed).

    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call

    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

!HALO - Need a global sum here (currently done inside calc_iareaf_iareag)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------
    call glide_prof_start(model,model%glide_prof%isos_water)

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

    call glide_prof_stop(model,model%glide_prof%isos_water)
    
    ! basal shear stress calculations

    !Halo updates required for inputs to glide_set_mask?
!HALO - If these values are needed, it should be possible to compute them without halo updates,
!       provided that thck and usrf have been updated in halos
    ! call staggered_parallel_halo(model%geomderv%stagthck) prior to calc_gline_flux
    ! call staggered_parallel_halo(model%geomderv%dusrfdew) prior to glide_calcstrsstr
    ! call staggered_parallel_halo(model%geomderv%dusrfdns) prior to glide_calcstrsstr

!TODO - This call is useful only for SIA diagnostics; not needed for glissade driver.
    call calc_basal_shear(model%geomderv%stagthck, &
                          model%geomderv%dusrfdew, &
                          model%geomderv%dusrfdns, &
                          model%stress%tau_x, &
                          model%stress%tau_y)
    !Includes a halo update of model%stress%tau_x and model%stress%tau_y at end of call
    !HALO - That update will be removed

  end subroutine glide_tstep_p2

!HALO - Does tstep_p3 require any halo info?

  subroutine glide_tstep_p3(model, no_write)
    !*FD Performs third part of time-step of an ice model instance:
    !*FD calculate isostatic adjustment and upper and lower ice surface
    use parallel
    use isostasy
    use glide_setup
    use glide_velo, only: gridwvel
    use glide_thck, only: timeders    

    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    logical, optional, intent(in) :: no_write
    logical nw
#ifdef GLC_DEBUG
    integer :: i, j, k, upn 
       upn = model%general%upn

       i = itest
       j = jtest
       write(6,*) ' '
       write(6,*) 'Starting tstep_p3, i, j, thck =', i, j, model%geometry%thck(i,j)
#endif

!TODO - Are we supporting an isostasy calculation in the parallel model?
!       While this may be a low priority in the near term, we should do so eventually.

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
    call glide_prof_start(model,model%glide_prof%isos)
    if (model%isos%do_isos) then
       !JEFF the isos_isostasy() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       call isos_isostasy(model)
    end if
    call glide_prof_stop(model,model%glide_prof%isos)

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------

! Calculate time-derivatives of thickness and upper surface elevation ------------
    !Halo updates required for inputs to glide_calcsrf?
    ! call parallel_halo(model%geometry%thck) within glide_marinlim
    ! call parallel_halo(model%geometry%topg) before previous call to glide_set_mask

!HALO - Verify that both thck and topg are up to date.
!       Note that glide_calclsrf loops over all cells (not just locally owned)

    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    !If input halos are up to date, halo update for model%geometry%lsrf should not be necessary

    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)
    !If input halos are up to date, halo update for model%geometry%usrf should not be necessary

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glide_tstep_p3

  !-------------------------------------------------------------------

  subroutine glide_tstep_postp3(model, no_write)
    !* This routine does the parallel routines and output that was in _p3()
    !* _p3() is executed in serial on main node.  These are parallel operations.
    !* Jeff Nichols, created for Bill Lipscomb September 2011
    use parallel

    use glide_setup
    use glide_velo, only: gridwvel
    use glide_thck, only: timeders

    implicit none
    type(glide_global_type) :: model        !*FD model instance

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

#ifdef GLC_DEBUG
    i = itest
    j = jtest
    write(6,*) ' '
    write(6,*) 'Before restart write, i, j, thck =', i, j, model%geometry%thck(i,j)
    write(6,300) k, model%geometry%thck(i,j)
    write(6,*) ' '
    write(6,*) 'k, temperature'
    do k = 1, upn
         write(6,300) k, model%temper%temp(k,i,j)
    enddo
300  format(i3, Z24.20)
#endif

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
       if (model%options%gthf.gt.0) then
         call t_startf('glide_lithot_io_writeall')
          call glide_lithot_io_writeall(model,model)
         call t_stopf('glide_lithot_io_writeall')
       end if
    end if

  end subroutine glide_tstep_postp3

  !-------------------------------------------------------------------

end module glide
