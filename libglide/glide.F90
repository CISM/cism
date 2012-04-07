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

    ! initialise bed softness to uniform parameter
    model%velocity%bed_softness = model%velowk%btrac_const

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
    call init_velo(model)

    ! MJH: Initialize temperature field - this needs to happen after input file is
    ! read so we can assign artm (which could possibly be read in) if temp has not been input.
    if (model%options%whichtemp == TEMP_REMAP_ADV) then
       call glissade_init_temp(model)
    else
       call glide_init_temp(model)
    endif 

    call init_thck(model)

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

    if (model%options%which_ho_diagnostic == HO_DIAG_PP ) then

        call glam_velo_fordsiapstr_init(model%general%ewn,    model%general%nsn,  &
                                        model%general%upn,                        &
                                        model%numerics%dew,   model%numerics%dns, &
                                        model%numerics%sigma)
    end if

    if ((model%options%whichevol == EVOL_INC_REMAP ) .or. &
       (model%options%whichevol == EVOL_NO_THICKNESS)) then
      if (old_remapping) then

        if (model%options%whichtemp == TEMP_REMAP_ADV) then ! Use IR to advect temperature

#ifdef JEFFORIG
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       model%general%ewn,  model%general%nsn,   &
                                       model%options%periodic_ew, model%options%periodic_ns, &
                                       model%general%upn,  model%numerics%sigma )    
#endif
           !JEFF In distributed, horizontal remapping is serialized, so it needs to be initialized with full grid size
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       global_ewn,  global_nsn,   &
                                       model%options%periodic_ew, model%options%periodic_ns, &
                                       model%general%upn,  model%numerics%sigma )

        else  ! Use IR to transport thickness only
#ifdef JEFFORIG
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       model%general%ewn,  model%general%nsn, &
                                       model%options%periodic_ew, model%options%periodic_ns))
#endif
           !JEFF In distributed, horizontal remapping is serialized, so it needs to be initialized with full grid size
           call horizontal_remap_init( model%remap_wk,    &
                                       model%numerics%dew, model%numerics%dns,  &
                                       global_ewn,  global_nsn, &
                                       model%options%periodic_ew, model%options%periodic_ns)
        endif ! whichtemp

      endif 
    endif 

    ! *sfp** added for summer modeling school
    if (model%options%whichevol== EVOL_FO_UPWIND ) then

        ! JEFF - Review for distributed. OK for parallel Trilinos. call not_parallel(__FILE__,__LINE__)
        call fo_upwind_advect_init( model%general%ewn, model%general%nsn )

    endif

    ! *mb* added; initialization of basal proc. module
    if (model%options%which_bmod == BAS_PROC_FULLCALC .or. &
        model%options%which_bmod == BAS_PROC_FASTCALC) then
        
!        call Basal_Proc_init (model%general%ewn, model%general%nsn,model%basalproc,     &
!                              model%numerics%ntem)
    write(*,*)"ERROR: Basal processes module is not supported in this release of CISM."
    stop

    end if      

    ! initialise ice age
    !! This is a placeholder; currently the ice age is not computed.  
    !! lipscomb - to do - Compute and advect the ice age.
    !! Currently the ice age is only computed for remapping transport
    !! (whichevol = 3 or 4)
    model%geometry%age(:,:,:) = 0._dp

!MJH moved flwa init to glide_init_temp/glissade_init_temp
!    if (model%options%hotstart.ne.1) then
!       ! initialise Glen's flow parameter A using an isothermal temperature distribution
!       call glide_temp_driver(model,0,0)
!    end if       

    ! calculate mask
    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)

    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)

    !calculate the normal at the marine margin
    call glide_marine_margin_normal(model%geometry%thck, model%geometry%thkmask, model%geometry%marine_bc_normal)

    ! calculate mask
    !if (model%options%hotstart.ne.1) then  ! setting the mask destroys exact restart
    !   call glide_set_mask(model)
    !end if

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise thckwk variables; used in timeders subroutine
    model%thckwk%olds(:,:,1) = model%geometry%thck(:,:)
    model%thckwk%olds(:,:,2) = model%geometry%usrf(:,:)

    ! initialise profile
#ifdef PROFILING
    call glide_prof_init(model)
#endif

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

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%geomderv)
#endif

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
    
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%geomderv)
#endif

#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask1)
#endif
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
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask1)
#endif

    ! ------------------------------------------------------------------------ 
    ! calculate geothermal heat flux
    ! ------------------------------------------------------------------------ 
    if (model%options%gthf.gt.0) then
       call not_parallel(__FILE__,__LINE__)
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%temperature)
#endif
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then

       if (model%options%whichtemp == TEMP_REMAP_ADV) then 

         ! Vert diffusion and strain heating only; no advection
         ! Remapping routine is used to advect temperature in glide_tstep_p2

         call glissade_temp_driver(model)

       else

         ! standard Glide driver, including temperature advection

         call glide_temp_driver(model, model%options%whichtemp, model%options%which_ho_diagnostic)

       endif

       model%temper%newtemps = .true.

    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%temperature)
#endif

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
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_evo)
#endif

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

       call inc_remap_driver( model )

       !Halo updates required for inputs to glide_stress?
       call staggered_parallel_halo(model%geomderv%dusrfdew)
       call staggered_parallel_halo(model%geomderv%dusrfdns)
       call staggered_parallel_halo(model%geomderv%dthckdew)
       call staggered_parallel_halo(model%geomderv%dthckdns)
       ! call parallel_halo(model%geometry%thck) in inc_remap_driver
       ! call staggered_parallel_halo(model%velocity%uvel) in inc_remap_driver
       ! call staggered_parallel_halo(model%velocity%vvel) in inc_remap_driver
       call parallel_halo(model%stress%efvs)

       !Tau is calculated in glide_stress and initialized in glide_types.
       call glide_calcstrsstr( model )       !*sfp* added for populating stress tensor w/ HO fields
       !Includes halo updates of 
       ! model%stress%tau%xx, model%stress%tau%yy, model%stress%tau%xy,
       ! model%stress%tau%scalar, model%stress%tau%xz, model%stress%tau%yz
       !at end of call

       if (model%options%whichevol .eq. EVOL_NO_THICKNESS) then
          ! restore old thickness
          model%geometry%thck = thck_old
          model%geomderv%stagthck = stagthck_old
       endif

    case(EVOL_FO_UPWIND) ! Use first order upwind scheme for mass transport
       call not_parallel(__FILE__,__LINE__)
       call fo_upwind_advect_driver( model )

end select

#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_evo)
#endif

    ! ------------------------------------------------------------------------
    ! get new mask
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask2)
#endif

    !Halo updates required for inputs to glide_set_mask?
    ! call parallel_halo(model%geometry%thck) in inc_remap_driver
    call parallel_halo(model%geometry%topg)

    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)
    !Includes a halo update of model%geometry%thkmask at end of call

#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask2)
#endif

    !Halo updates required for inputs to glide_set_mask?
    ! call parallel_halo(model%geometry%thck) in inc_remap_driver
    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call

    call glide_marine_margin_normal(model%geometry%thck, model%geometry%thkmask, model%geometry%marine_bc_normal)
    !Includes a halo update of model%geometry%marine_bc_normal at end of call

    !Halo updates required for inputs to calc_gline_flux?
    call staggered_parallel_halo(model%geomderv%stagthck)
    call staggered_parallel_halo(model%velocity%surfvel)
    ! call parallel_halo(model%geometry%thkmask) in earlier glide_set_mask call
    call staggered_parallel_halo(model%velocity%ubas)
    call staggered_parallel_halo(model%velocity%vbas)
    ! call parallel_halo(model%ground%gline_flux) - halo update not required for correctness

    !calculate the grounding line flux after the mask is correct
    call calc_gline_flux(model%geomderv%stagthck,model%velocity%surfvel, &
                         model%geometry%thkmask,model%ground%gline_flux, model%velocity%ubas, &
                         model%velocity%vbas, model%numerics%dew)
    !Includes a halo update of model%ground%gline_flux at end of call

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 

    !Halo updates required for inputs to glide_marinlim?
    !(Note: case specific, so push into glide_marinlim?)
    ! call parallel_halo(model%geometry%thck) in inc_remap_driver
    call parallel_halo(model%isos%relx)
    ! call parallel_halo(model%geometry%topg) before previous set_glide_mask
    call parallel_halo(model%temper%flwa)
    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call
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
    !Note that halo updates for 
    ! model%geometry%thkmask not needed in current implementation of case 3
    ! model%climate%eus needed only for disabled case 6
    ! model%ground components needed only for disabled case 6

    !issues with ice shelf, calling it again fixes the mask

    !Halo updates required for inputs to glide_set_mask?
    ! call parallel_halo(model%geometry%thck) within glide_marinlim
    ! call parallel_halo(model%geometry%topg) before previous call to glide_set_mask

    call glide_set_mask(model%numerics, model%geometry%thck, model%geometry%topg, &
                        model%general%ewn, model%general%nsn, model%climate%eus, &
                        model%geometry%thkmask, model%geometry%iarea, model%geometry%ivol)
    !Includes a halo update of model%geometry%thkmask at end of call

    !Halo updates required for inputs to calc_iareaf_iareag?
    ! call parallel_halo(model%geometry%thkmask) in previous glide_set_mask call

    call calc_iareaf_iareag(model%numerics%dew,model%numerics%dns, &
                            model%geometry%iarea, model%geometry%thkmask, &
                            model%geometry%iareaf, model%geometry%iareag)
    !No output requiring halo update (though global_sum called in routine)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos_water)
#endif
    if (model%isos%do_isos) then
       !JEFF the isos_icewaterload() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       if (model%numerics%time.ge.model%isos%next_calc) then
          model%isos%next_calc = model%isos%next_calc + model%isos%period
          call isos_icewaterload(model)
          model%isos%new_load = .true.
       end if
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos_water)
#endif
    
    ! basal shear stress calculations

    !Halo updates required for inputs to glide_set_mask?
    ! call staggered_parallel_halo(model%geomderv%stagthck) prior to calc_gline_flux
    ! call staggered_parallel_halo(model%geomderv%dusrfdew) prior to glide_calcstrsstr
    ! call staggered_parallel_halo(model%geomderv%dusrfdns) prior to glide_calcstrsstr

    call calc_basal_shear(model%geomderv%stagthck, &
                          model%geomderv%dusrfdew, &
                          model%geomderv%dusrfdns, &
                          model%stress%tau_x, &
                          model%stress%tau_y)
    !Includes a halo update of model%stress%tau_x and model%stress%tau_y at end of call

  end subroutine glide_tstep_p2

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

    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos)
#endif
    if (model%isos%do_isos) then
       !JEFF the isos_isostasy() is passed the entire model, so I don't know what gathered variables it needs.
       call not_parallel(__FILE__, __LINE__)

       call isos_isostasy(model)
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos)
#endif

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------

! Calculate time-derivatives of thickness and upper surface elevation ------------
    !Halo updates required for inputs to glide_calcsrf?
    ! call parallel_halo(model%geometry%thck) within glide_marinlim
    ! call parallel_halo(model%geometry%topg) before previous call to glide_set_mask

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

    ! These three calls are in glide_temp in the full-temperature section replaced by Bill's new temperature code.
    ! Commenting out until I hear otherwise.
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

    ! Calculate the vertical velocity of the grid ------------------------------------

    call gridwvel(model%numerics%sigma,  &
            model%numerics%thklim, &
            model%velocity%uvel,   &
            model%velocity%vvel,   &
            model%geomderv,        &
            model%geometry%thck,   &
            model%velocity%wgrd)

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
       call glide_io_writeall(model,model)
       if (model%options%gthf.gt.0) then
          call glide_lithot_io_writeall(model,model)
       end if
    end if
  end subroutine glide_tstep_postp3

  !-------------------------------------------------------------------

end module glide
