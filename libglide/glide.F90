! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide.f90 - part of the GLIMMER ice model                + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide
  !*FD the top-level GLIDE module

  use glide_types
  use glide_stop
  use glide_nc_custom
  use glide_io
  use glide_lithot
  use glide_profile
  use glimmer_config
  integer, private, parameter :: dummyunit=99

contains

  subroutine glide_config(model,config)
    !*FD read glide configuration from file and print it to the log
    use glide_setup
    use isostasy
    use glimmer_ncparams
    use glimmer_config
    use glimmer_map_init
    use glimmer_filenames
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file

    type(ConfigSection), pointer :: ncconfig
   
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
       call ConfigRead(process_path(model%funits%ncfile),ncconfig)
    end if
    call glimmer_nc_readparams(model,ncconfig)
  end subroutine glide_config

  subroutine glide_initialise(model)
    !*FD initialise GLIDE model instance
    use glide_setup
    use glimmer_ncio
    use glide_velo
    use glide_thck
    use glide_temp
    use glimmer_log
    use glide_mask
    use glimmer_scales
    use glide_mask
    use isostasy
    use glimmer_map_init
    implicit none
    type(glide_global_type) :: model        !*FD model instance

    character(len=100), external :: glimmer_version_char

    call write_log(trim(glimmer_version_char()))

    ! initialise scales
    call glimmer_init_scales

    ! scale parameters
    call glide_scale_params(model)
    ! set up coordinate systems
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
    
    ! load sigma file
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
    case(1) ! Supplied topography is relaxed
       model%isos%relx = model%geometry%topg
    case(2) ! Supplied topography is in equilibrium
       call isos_relaxed(model)
    end select

    ! set uniform basal heat flux
    model%temper%bheatflx = model%paramets%geot

    ! open all output files
    call openall_out(model)
    ! create glide variables
    call glide_io_createall(model)

    ! initialise glide components
    call init_velo(model)
    call init_temp(model)
    call init_thck(model)
    if (model%options%gthf.gt.0) then
       call init_lithot(model)
    end if

    if (model%options%hotstart.ne.1) then
       ! initialise Glen's flow parameter A using an isothermal temperature distribution
       call timeevoltemp(model,0)
    end if

    ! calculate mask
    call glide_set_mask(model)

    ! and calculate lower and upper ice surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! initialise profile
#ifdef PROFILING
    call glide_prof_init(model)
#endif
  end subroutine glide_initialise
  
  subroutine glide_tstep_p1(model,time)
    !*FD Performs first part of time-step of an ice model instance.
    !*FD calculate velocity and temperature
    use glimmer_global, only : rk
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    implicit none

    type(glide_global_type) :: model        !*FD model instance
    real(rk),  intent(in)   :: time         !*FD Current time in years

    ! Update internal clock
    model%numerics%time=time  
    model%temper%newtemps = .false.

    ! ------------------------------------------------------------------------ 
    ! Calculate various derivatives...
    ! ------------------------------------------------------------------------     
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%geomderv)
#endif
    call stagvarb(model%geometry% thck, &
         model%geomderv% stagthck,&
         model%general%  ewn, &
         model%general%  nsn)

    call geomders(model%numerics, &
         model%geometry% usrf, &
         model%geomderv% stagthck,&
         model%geomderv% dusrfdew, &
         model%geomderv% dusrfdns)

    call geomders(model%numerics, &
         model%geometry% thck, &
         model%geomderv% stagthck,&
         model%geomderv% dthckdew, &
         model%geomderv% dthckdns)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%geomderv)
#endif

#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_mask1)
#endif
    call glide_maskthck(0, &                                    !magi a hack, someone explain what whichthck=5 does
         model%geometry% thck,      &
         model%climate%  acab,      &
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
       call calc_lithot(model)
    end if

    ! ------------------------------------------------------------------------ 
    ! Calculate temperature evolution and Glenn's A, if necessary
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%temperature)
#endif
    if ( model%numerics%tinc >  mod(model%numerics%time,model%numerics%ntem)) then
       call timeevoltemp(model, model%options%whichtemp)
       model%temper%newtemps = .true.
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%temperature)
#endif

    ! ------------------------------------------------------------------------ 
    ! Calculate basal traction factor
    ! ------------------------------------------------------------------------ 
    call calc_btrc(model,model%options%whichbtrc,model%velocity%btrc)

  end subroutine glide_tstep_p1


  subroutine glide_tstep_p2(model,no_write)
    !*FD Performs second part of time-step of an ice model instance.
    !*FD write data and move ice
    use glide_thck
    use glide_velo
    use glide_setup
    use glide_temp
    use glide_mask
    use isostasy
    implicit none

    type(glide_global_type) :: model        !*FD model instance
    logical,optional :: no_write

    logical nw

    ! ------------------------------------------------------------------------ 
    ! write to netCDF file
    ! ------------------------------------------------------------------------ 

    if (present(no_write)) then
       nw=no_write
    else
       nw=.false.
    end if

    if (.not. nw) call glide_io_writeall(model,model)

    ! ------------------------------------------------------------------------ 
    ! Calculate flow evolution by various different methods
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%ice_evo)
#endif
    select case(model%options%whichevol)
    case(0) ! Use precalculated uflx, vflx -----------------------------------

       call thck_lin_evolve(model,model%temper%newtemps)

    case(1) ! Use explicit leap frog method with uflx,vflx -------------------

       call stagleapthck(model,model%temper%newtemps)

    case(2) ! Use non-linear calculation that incorporates velocity calc -----

       call thck_nonlin_evolve(model,model%temper%newtemps)

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
    call glide_set_mask(model)
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%ice_mask2)
#endif

    ! ------------------------------------------------------------------------ 
    ! Remove ice which is either floating, or is present below prescribed
    ! depth, depending on value of whichmarn
    ! ------------------------------------------------------------------------ 
    call glide_marinlim(model%options%  whichmarn, &
         model%geometry% thck,      &
         model%isos% relx,      &
         model%geometry%topg,   &
         model%geometry%thkmask,    &
         model%numerics%mlimit,     &
         model%numerics%calving_fraction, &
         model%climate%eus,         &
         model%climate%calving)

    ! ------------------------------------------------------------------------
    ! update ice/water load if necessary
    ! ------------------------------------------------------------------------
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos_water)
#endif
    if (model%isos%do_isos) then
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
    call calc_basal_shear(model)
  end subroutine glide_tstep_p2

  subroutine glide_tstep_p3(model)
    !*FD Performs third part of time-step of an ice model instance:
    !*FD calculate isostatic adjustment and upper and lower ice surface
    use isostasy
    use glide_setup
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    
    ! ------------------------------------------------------------------------ 
    ! Calculate isostasy
    ! ------------------------------------------------------------------------ 
#ifdef PROFILING
    call glide_prof_start(model,model%glide_prof%isos)
#endif
    if (model%isos%do_isos) then
       call isos_isostasy(model)
    end if
#ifdef PROFILING
    call glide_prof_stop(model,model%glide_prof%isos)
#endif

    ! ------------------------------------------------------------------------
    ! calculate upper and lower ice surface
    ! ------------------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

    ! increment time counter
    model%numerics%timecounter = model%numerics%timecounter + 1

  end subroutine glide_tstep_p3

  !-------------------------------------------------------------------

!MH!  subroutine glide_write_mod_rst(rfile)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart_common
!MH!
!MH!#ifdef RESTARTS
!MH!    use glide_types
!MH!    use isostasy_types
!MH!#endif
!MH!
!MH!    type(restart_file) :: rfile      !*FD Open restart file 
!MH!
!MH!#ifdef RESTARTS
!MH!    call glide_types_modrsw(rfile)
!MH!    call isostasy_types_modrsw(rfile)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_write_mod_rst
!MH!
!MH!  !-------------------------------------------------------------------
!MH!
!MH!  subroutine glide_read_mod_rst(rfile)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart_common
!MH!
!MH!#ifdef RESTARTS
!MH!    use glide_types
!MH!    use isostasy_types
!MH!#endif
!MH!
!MH!    type(restart_file) :: rfile      !*FD Open restart file 
!MH!
!MH!#ifdef RESTARTS
!MH!    call glide_types_modrsr(rfile)
!MH!    call isostasy_types_modrsr(rfile)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_read_mod_rst
!MH!
!MH!  !-------------------------------------------------------------------
!MH!
!MH!  subroutine glide_write_restart(model,rfile)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart
!MH!    use glimmer_restart_common
!MH!    implicit none
!MH!
!MH!    type(glide_global_type) :: model !*FD model instance
!MH!    type(restart_file) :: rfile      !*FD Open restart file     
!MH!
!MH!#ifdef RESTARTS
!MH!    call glimmer_write_mod_rst(rfile)
!MH!    call glide_write_mod_rst(rfile)
!MH!    call rsw_glide_global_type(rfile,model)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_write_restart

!MH!  !-------------------------------------------------------------------
!MH!
!MH!  subroutine glide_read_restart(model,rfile,prefix)
!MH!
!MH!    use glimmer_log
!MH!    use glimmer_restart
!MH!    use glimmer_restart_common
!MH!    use glimmer_ncdf
!MH!    use glimmer_ncio
!MH!    implicit none
!MH!
!MH!    type(glide_global_type) :: model !*FD model instance
!MH!    type(restart_file) :: rfile      !*FD Open restart file 
!MH!    character(*),optional,intent(in) :: prefix !*FD prefix for new output files
!MH!
!MH!    character(40) :: pf
!MH!
!MH!    if (present(prefix)) then
!MH!       pf = prefix
!MH!    else
!MH!       pf = 'RESTART_'
!MH!    end if
!MH!
!MH!#ifdef RESTARTS
!MH!    call glimmer_read_mod_rst(rfile)
!MH!    call glide_read_mod_rst(rfile)
!MH!    call rsr_glide_global_type(rfile,model)
!MH!    call nc_repair_outpoint(model%funits%out_first)
!MH!    call nc_repair_inpoint(model%funits%in_first)
!MH!    call nc_prefix_outfiles(model%funits%out_first,trim(pf))
!MH!    call openall_out(model)
!MH!    call glide_io_createall(model)
!MH!    call glide_nc_fillall(model)
!MH!#else
!MH!    call write_log('No restart code available - rebuild GLIMMER with --enable-restarts',GM_FATAL)
!MH!#endif
!MH!
!MH!  end subroutine glide_read_restart

end module glide
