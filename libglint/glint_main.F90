! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_main.f90 - part of the Glimmer-CISM ice model      + 
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

module glint_main

  !*FD  This is the main glimmer module, which contains the top-level 
  !*FD  subroutines and derived types comprising the glimmer ice model.

  use glimmer_global
  use glint_type
  use glint_global_grid
  use glint_constants
  use glimmer_anomcouple
  use glimmer_paramets, only: idiag, jdiag
  use glide_diagnostics
  use glimmer_paramets, only: itest, jtest, jjtest, stdout  

  ! ------------------------------------------------------------
  ! GLIMMER_PARAMS derived type definition
  ! This is where default values are set.
  ! ------------------------------------------------------------

  type glint_params 

     !*FD Derived type containing parameters relevant to all instances of 
     !*FD the model - i.e. those parameters which pertain to the global model. 

     ! Global grids used ----------------------------------------

     type(global_grid) :: g_grid      !*FD The main global grid, used for 
                                      !*FD input and most outputs
     type(global_grid) :: g_grid_orog !*FD Global grid used for orography output.

     ! Ice model instances --------------------------------------

     integer                                   :: ninstances = 1       !*FD Number of ice model instances
     type(glint_instance),pointer,dimension(:) :: instances  => null() !*FD Array of glimmer\_instances

     ! Global model parameters ----------------------------------

     integer  :: tstep_mbal = 1        !*FD Mass-balance timestep (hours)
     integer  :: start_time            !*FD Time of first call to glint (hours)
     integer  :: time_step             !*FD Calling timestep of global model (hours)

     ! Parameters that can be set by the GCM calling Glint

     logical  :: gcm_smb = .false.     !*FD If true, receive surface mass balance from the GCM 
     logical  :: gcm_restart = .false. !*FD If true, hotstart the model from a GCM restart file
     character(fname_length) :: gcm_restart_file   !*FD Name of restart file
     integer  :: gcm_fileunit = 99     !*FD Fileunit specified by GCM for reading config files
   
     ! Averaging parameters -------------------------------------

     integer  :: av_start_time = 0   !*FD Holds the value of time from 
                                     !*FD the last occasion averaging was restarted (hours)
     integer  :: av_steps      = 0   !*FD Holds the number of times glimmer has 
                                     !*FD been called in current round of averaging.
     integer  :: next_av_start = 0   !*FD Time when we expect next averaging to start
     logical  :: new_av     = .true. !*FD Set to true if the next correct call starts a new averaging round

     ! Averaging arrays -----------------------------------------

     real(rk),pointer,dimension(:,:) :: g_av_precip  => null()  !*FD globally averaged precip
     real(rk),pointer,dimension(:,:) :: g_av_temp    => null()  !*FD globally averaged temperature 
     real(rk),pointer,dimension(:,:) :: g_max_temp   => null()  !*FD global maximum temperature
     real(rk),pointer,dimension(:,:) :: g_min_temp   => null()  !*FD global minimum temperature
     real(rk),pointer,dimension(:,:) :: g_temp_range => null()  !*FD global temperature range
     real(rk),pointer,dimension(:,:) :: g_av_zonwind => null()  !*FD globally averaged zonal wind 
     real(rk),pointer,dimension(:,:) :: g_av_merwind => null()  !*FD globally averaged meridional wind 
     real(rk),pointer,dimension(:,:) :: g_av_humid   => null()  !*FD globally averaged humidity (%)
     real(rk),pointer,dimension(:,:) :: g_av_lwdown  => null()  !*FD globally averaged downwelling longwave (W/m$^2$)
     real(rk),pointer,dimension(:,:) :: g_av_swdown  => null()  !*FD globally averaged downwelling shortwave (W/m$^2$)
     real(rk),pointer,dimension(:,:) :: g_av_airpress => null() !*FD globally averaged surface air pressure (Pa)
     real(rk),pointer,dimension(:,:,:) :: g_av_qsmb => null()   ! globally averaged surface mass balance (kg m-2 s-1)
     real(rk),pointer,dimension(:,:,:) :: g_av_tsfc => null()   ! globally averaged surface temperature (deg C)
     real(rk),pointer,dimension(:,:,:) :: g_av_topo => null()   ! globally averaged surface elevation   (m)

     ! Fractional coverage information --------------------------

     real(rk),pointer,dimension(:,:) :: total_coverage  => null()     !*FD Fractional coverage by 
                                                                      !*FD all ice model instances.
     real(rk),pointer,dimension(:,:) :: cov_normalise   => null()     !*FD Normalisation values 
                                                                      !*FD for coverage calculation.
     real(rk),pointer,dimension(:,:) :: total_cov_orog  => null()     !*FD Fractional coverage by 
                                                                      !*FD all ice model instances (orog).
     real(rk),pointer,dimension(:,:) :: cov_norm_orog   => null()     !*FD Normalisation values 
                                                                      !*FD for coverage calculation (orog).
     logical                         :: coverage_calculated = .false. !*FD Have we calculated the
                                                                      !*FD coverage map yet?
     ! File information -----------------------------------------

     character(fname_length) :: paramfile      !*FD Name of global parameter file

     ! Accumulation/averaging flags -----------------------------

     logical :: need_winds=.false. !*FD Set if we need the winds to be accumulated/downscaled
     logical :: enmabal=.false.    !*FD Set if we're using the energy balance mass balance model anywhere

     ! Anomaly coupling for global climate ------------------------------------------

     type(anomaly_coupling) :: anomaly_params !*FD Parameters for anomaly coupling

  end type glint_params

  ! Private names -----------------------------------------------

  private glint_allocate_arrays
  private glint_readconfig,calc_bounds,check_init_args

    !---------------------------------------------------------------------------------------
    ! Some notes on coupling to the Community Earth System Model (CESM).  These may be applicable
    ! for coupling to other GCMs:
    !
    ! If gcm_smb is true, then Glint receives three fields from the CESM coupler on a global grid
    ! in each of several elevation classes:
    !   qsmb = surface mass balance (kg/m^2/s)
    !   tsfc = surface ground temperature (deg C)
    !   topo = surface elevation (m)
    ! Both qsmb and tsfc are computed in the CESM land model.
    ! Five fields are returned to CESM on the global grid for each elevation class:
    !   gfrac = fractional ice coverage
    !   gtopo = surface elevation (m)
    !   grofi = ice runoff (i.e., calving) (kg/m^2/s)
    !   grofl = liquid runoff (i.e., basal melting; the land model handles sfc runoff) (kg/m^2/s)
    !   ghflx = heat flux from the ice interior to the surface (W/m^2)
    ! The land model has the option to update its ice coverage and surface elevation, given
    ! the fields returned from Glint.
    !
    ! If gcm_smb is false, then the CESM wrapper that calls Glint receives three fields of the same name
    ! from the coupler, but the meanings of qsmb and tsfc are different:
    !   qsmb = net precipitation (kg/m^2/s)
    !   tsfc = 2-m air temperature (deg C)
    !   topo = surface elevation (m)
    ! These fields are received for a single elevation class.  The precip and 2-m temperature are
    ! sent to Glint as inputs to a daily PDD scheme.  Glint does not return values for gfrac, 
    ! gtopo, etc., and the land model will not update its ice coverage and surface elevation.
    ! 
    ! Glimmer-CISM will normally be coupled to CESM with gcm_smb = T.
    ! The PDD option is included for comparison to other models with PDD schemes.
    !---------------------------------------------------------------------------------------

contains

  subroutine initialise_glint(params,                         &
                              lats,         longs,            &
                              time_step,    paramfile,        &
                              latb,         lonb,             &
                              orog,         albedo,           &
                              ice_frac,     veg_frac,         &
                              snowice_frac, snowveg_frac,     &
                              snow_depth,                     &
                              orog_lats,    orog_longs,       &
                              orog_latb,    orog_lonb,        &
                              output_flag,  daysinyear,       &
                              snow_model,   ice_dt,           &
                              extraconfigs, start_time,       &
                              gcm_nec,      gcm_smb,          &
                              gfrac,        gtopo,            &
                              grofi,        grofl,            &
                              ghflx,        gmask,            &
                              gcm_restart,  gcm_restart_file, &
                              gcm_debug,    gcm_fileunit)

    !*FD Initialises the model

    use glimmer_config
    use glint_initialise
    use glimmer_log
    use glimmer_filenames
    implicit none

    ! Subroutine argument declarations --------------------------------------------------------

    type(glint_params),              intent(inout) :: params      !*FD parameters to be set
    real(rk),dimension(:),           intent(in)    :: lats,longs  !*FD location of gridpoints 
                                                                  !*FD in global data.
    integer,                         intent(in)    :: time_step   !*FD Timestep of calling model (hours)
    character(*),dimension(:),       intent(in)    :: paramfile   !*FD array of configuration filenames.
    real(rk),dimension(:),  optional,intent(in)    :: latb        !*FD Locations of the latitudinal 
                                                                  !*FD boundaries of the grid-boxes.
    real(rk),dimension(:),  optional,intent(in)    :: lonb        !*FD Locations of the longitudinal
                                                                  !*FD boundaries of the grid-boxes.
    real(rk),dimension(:,:),optional,intent(out)   :: orog        !*FD Initial global orography
    real(rk),dimension(:,:),optional,intent(out)   :: albedo      !*FD Initial albedo
    real(rk),dimension(:,:),optional,intent(out)   :: ice_frac    !*FD Initial ice fraction 
    real(rk),dimension(:,:),optional,intent(out)   :: veg_frac    !*FD Initial veg fraction
    real(rk),dimension(:,:),optional,intent(out)   :: snowice_frac !*FD Initial snow-covered ice fraction
    real(rk),dimension(:,:),optional,intent(out)   :: snowveg_frac !*FD Initial snow-covered veg fraction
    real(rk),dimension(:,:),optional,intent(out)   :: snow_depth  !*FD Initial snow depth 
    real(rk),dimension(:),  optional,intent(in)    :: orog_lats   !*FD Latitudinal location of gridpoints 
                                                                  !*FD for global orography output.
    real(rk),dimension(:),  optional,intent(in)    :: orog_longs  !*FD Longitudinal location of gridpoints 
                                                                  !*FD for global orography output.
    real(rk),dimension(:),  optional,intent(in)    :: orog_latb   !*FD Locations of the latitudinal 
                                                                  !*FD boundaries of the grid-boxes (orography).
    real(rk),dimension(:),  optional,intent(in)    :: orog_lonb   !*FD Locations of the longitudinal
                                                                  !*FD boundaries of the grid-boxes (orography).
    logical,                optional,intent(out)   :: output_flag !*FD Flag to show output set (provided for
                                                                  !*FD consistency)
    integer,                optional,intent(in)    :: daysinyear  !*FD Number of days in the year
    logical,                optional,intent(out)   :: snow_model  !*FD Set if the mass-balance scheme has a snow-depth model
    integer,                optional,intent(out)   :: ice_dt      !*FD Ice dynamics time-step in hours
    type(ConfigData),dimension(:),optional ::  extraconfigs !*FD Additional configuration information - overwrites
                                                                  !*FD config data read from files
    integer,                optional,intent(in)    :: start_time  !*FD Time of first call to glint (hours)
    integer,                  optional,intent(in)  :: gcm_nec     !*FD number of elevation classes for GCM input
    logical,                  optional,intent(in)  :: gcm_smb     !*FD true if getting sfc mass balance from a GCM
    real(rk),dimension(:,:,:),optional,intent(out) :: gfrac       !*FD ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(out) :: gtopo       !*FD surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofi       !*FD ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofl       !*FD liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: ghflx       !*FD heat flux (W/m^2, positive down)
    integer, dimension(:,:),  optional,intent(in)  :: gmask       !*FD mask = 1 where global data are valid
    logical,                  optional,intent(in)  :: gcm_restart ! logical flag to hotstart from a GCM restart file
    character(*),             optional, intent(in) :: gcm_restart_file ! hotstart filename for a GCM restart
                                                                  ! (currently assumed to be CESM)
    logical,                  optional,intent(in)  :: gcm_debug   ! logical flag from GCM to output debug information
    integer,                  optional,intent(in)  :: gcm_fileunit! fileunit for reading config files

    ! Internal variables -----------------------------------------------------------------------

    type(ConfigSection), pointer :: global_config, instance_config, section  ! configuration stuff
    character(len=100) :: message                                            ! For log-writing
    character(fname_length),dimension(:),pointer :: config_fnames=>null()    ! array of config filenames
    type(ConfigSection), pointer :: econf
    integer :: i
    real(rk),dimension(:,:),allocatable :: orog_temp, if_temp, vf_temp, sif_temp,  &
                                           svf_temp,  sd_temp, alb_temp      ! Temporary output arrays
    integer,dimension(:),allocatable :: mbts,idts ! Array of mass-balance and ice dynamics timesteps
    logical :: anomaly_check ! Set if we've already initialised anomaly coupling

    real(rk),dimension(:,:,:),allocatable ::   &
               gfrac_temp, gtopo_temp, grofi_temp, grofl_temp, ghflx_temp    ! Temporary output arrays
    integer :: n
    integer :: nec       ! number of elevation classes
    real(rk):: timeyr    ! time in years
    integer :: j, ii, jj

    if (present(gcm_debug)) then
       GLC_DEBUG = gcm_debug
    endif

    if (GLC_DEBUG) then
       write(stdout,*) 'Starting initialise_glint'
    endif

    ! Initialise start time and calling model time-step ----------------------------------------
    ! We ignore t=0 by default 

    params%time_step = time_step

    if (present(start_time)) then
       params%start_time = start_time
    else
       params%start_time = time_step
    end if

    params%next_av_start = params%start_time

    ! Initialisation for runs where the surface mass balance is received from a GCM
 
    params%gcm_smb = .false.
    if (present(gcm_smb)) then
       params%gcm_smb = gcm_smb
    endif

    params%gcm_restart = .false.
    if (present(gcm_restart)) then
       params%gcm_restart = gcm_restart
    endif

    params%gcm_restart_file = ''
    if (present(gcm_restart_file)) then
       params%gcm_restart_file = gcm_restart_file
    endif

    params%gcm_fileunit = 99
    if (present(gcm_fileunit)) then
       params%gcm_fileunit = gcm_fileunit
    endif

    nec = 1
    if (present(gcm_nec)) then
       nec = gcm_nec
    endif

    if (GLC_DEBUG) then
       write(stdout,*) 'time_step =', params%time_step
       write(stdout,*) 'start_time =', params%start_time
       write(stdout,*) 'next_av_start =', params%next_av_start
    endif

    ! Initialise year-length -------------------------------------------------------------------

    if (present(daysinyear)) then
       call glint_set_year_length(daysinyear)
    end if

    if (GLC_DEBUG) then
       write(stdout,*) 'Initialize global grid'
       write(stdout,*) 'present =', present(gmask)
    endif

    ! Initialise main global grid --------------------------------------------------------------

    if (present(gmask)) then
       call new_global_grid(params%g_grid, longs, lats, lonb=lonb, latb=latb, nec=nec, mask=gmask)
    else
       call new_global_grid(params%g_grid, longs, lats, lonb=lonb, latb=latb, nec=nec)
    endif

    if (GLC_DEBUG) then
       write (stdout,*) ' ' 
       write (stdout,*) 'time_step (hr)  =', params%time_step
       write (stdout,*) 'start_time (hr) =', params%start_time
       write (stdout,*) 'Called new_global_grid '
       write (stdout,*) 'g_grid%nx =', params%g_grid%nx
       write (stdout,*) 'g_grid%ny =', params%g_grid%ny
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lons =', params%g_grid%lons
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lats =', params%g_grid%lats
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lon_bound =', params%g_grid%lon_bound
       write (stdout,*) ' '
       write (stdout,*) 'g_grid%lat_bound =', params%g_grid%lat_bound
       do j = 5, 10
          write (stdout,*)
          write (stdout,*) 'j, g_grid%mask =', j, params%g_grid%mask(:,j)
       enddo
    endif

    ! Initialise orography grid ------------------------------------

    call check_init_args(orog_lats, orog_longs, orog_latb, orog_lonb)

    if (present(orog_lats) .and. present(orog_longs)) then
       call new_global_grid(params%g_grid_orog, orog_longs, orog_lats,  &
                            lonb=orog_lonb, latb=orog_latb)
    else
       call copy_global_grid(params%g_grid, params%g_grid_orog)
    end if

    ! Allocate arrays -----------------------------------------------

!lipscomb - TO DO - The following arrays may not be needed for gcm_smb runs
    call glint_allocate_arrays(params)

    if (params%gcm_smb) call glint_allocate_arrays_gcm(params)

    ! Initialise arrays ---------------------------------------------

    params%g_av_precip  = 0.0
    params%g_av_temp    = 0.0
    params%g_max_temp   = -1000.0
    params%g_min_temp   = 1000.0
    params%g_temp_range = 0.0
    params%g_av_zonwind = 0.0
    params%g_av_merwind = 0.0
    params%g_av_humid   = 0.0
    params%g_av_lwdown  = 0.0
    params%g_av_swdown  = 0.0
    params%g_av_airpress = 0.0

    if (params%gcm_smb) then
       params%g_av_qsmb    = 0.0
       params%g_av_tsfc    = 0.0
       params%g_av_topo    = 0.0
    endif

    ! ---------------------------------------------------------------
    ! Zero coverage maps and normalisation fields for main grid and
    ! orography grid
    ! ---------------------------------------------------------------

    params%total_coverage = 0.0
    params%total_cov_orog = 0.0

    params%cov_normalise = 0.0
    params%cov_norm_orog = 0.0

    if (GLC_DEBUG) then
       write(stdout,*) 'Read paramfile'
       write(stdout,*) 'paramfile =', paramfile
    endif

    ! ---------------------------------------------------------------
    ! Determine how many instances there are, according to what
    ! configuration files we've been provided with
    ! ---------------------------------------------------------------

    if (size(paramfile) == 1) then
       ! Load the configuration file into the linked list
       call ConfigRead(process_path(paramfile(1)), global_config, params%gcm_fileunit)    
       call glint_readconfig(global_config, params%ninstances, config_fnames, paramfile) ! Parse the list
    else
       params%ninstances = size(paramfile)
       allocate(config_fnames(params%ninstances))
       config_fnames = paramfile
    end if

    allocate(params%instances(params%ninstances))
    allocate(mbts(params%ninstances), idts(params%ninstances))

    if (GLC_DEBUG) then
       write(stdout,*) 'Number of instances =', params%ninstances
       write(stdout,*) 'Read config files and initialize each instance'
    endif
    ! ---------------------------------------------------------------
    ! Read config files, and initialise instances accordingly
    ! ---------------------------------------------------------------

    call write_log('Reading instance configurations')
    call write_log('-------------------------------')

    anomaly_check=.false.

    do i=1,params%ninstances

       call ConfigRead(process_path(config_fnames(i)),instance_config, params%gcm_fileunit)
       if (present(extraconfigs)) then
          if (size(extraconfigs)>=i) then
             call ConfigCombine(instance_config,extraconfigs(i))
          end if
       end if
 
      call glint_i_initialise(instance_config,     params%instances(i),     &
                               params%g_grid,      params%g_grid_orog,      &
                               mbts(i),            idts(i),                 &
                               params%need_winds,  params%enmabal,          &
                               params%start_time,  params%time_step,        &
                               params%gcm_restart, params%gcm_restart_file, &
                               params%gcm_fileunit )

       params%total_coverage = params%total_coverage + params%instances(i)%frac_coverage
       params%total_cov_orog = params%total_cov_orog + params%instances(i)%frac_cov_orog

       where (params%total_coverage > 0.0) params%cov_normalise = params%cov_normalise + 1.0
       where (params%total_cov_orog > 0.0) params%cov_norm_orog = params%cov_norm_orog + 1.0

       ! Write initial diagnostics for this instance
       timeyr = real(params%start_time/8760.)
       if (GLC_DEBUG) then
          write(stdout,*) 'Write model diagnostics, time =', timeyr
       endif
       call glide_write_diag(params%instances(i)%model, timeyr, idiag, jdiag)

       ! Initialise anomaly coupling
       if (.not.anomaly_check) then 
          call anomaly_init(params%anomaly_params, instance_config)
          if (params%anomaly_params%enabled .and. &
               (params%anomaly_params%nx/=params%g_grid%nx .or. &
                params%anomaly_params%ny/=params%g_grid%ny) ) then
             call write_log("Anomaly coupling grids have different "// &
                  "sizes to GLINT coupling grids",GM_FATAL,__FILE__,__LINE__)
          end if
          if (params%anomaly_params%enabled) anomaly_check=.true.
       end if

    end do

    ! Check that all mass-balance time-steps are the same length and 
    ! assign that value to the top-level variable

    params%tstep_mbal = check_mbts(mbts)
    if (present(ice_dt)) then
       ice_dt = check_mbts(idts)
    end if

    if (GLC_DEBUG) then
       write(stdout,*) 'tstep_mbal =', params%tstep_mbal
       write(stdout,*) 'start_time =', params%start_time
       write(stdout,*) 'time_step =',  params%time_step
       if (present(ice_dt)) write(stdout,*) 'ice_dt =', ice_dt
    endif

    ! Check time-steps divide into one another appropriately.

    if (.not.(mod(params%tstep_mbal,params%time_step)==0)) then
       print*,params%tstep_mbal,params%time_step
       call write_log('The mass-balance timestep must be an integer multiple of the forcing time-step', &
            GM_FATAL,__FILE__,__LINE__)
    end if

    ! Check we don't have coverage greater than one at any point.

    where (params%total_coverage > 1.0) params%total_coverage = 1.0
    where (params%total_cov_orog > 1.0) params%total_cov_orog = 1.0
    params%coverage_calculated=.true.

    ! Zero optional outputs, if present

    if (present(orog)) orog = 0.0
    if (present(albedo)) albedo = 0.0
    if (present(ice_frac)) ice_frac = 0.0
    if (present(veg_frac)) veg_frac = 0.0
    if (present(snowice_frac)) snowice_frac = 0.0
    if (present(snowveg_frac)) snowveg_frac = 0.0
    if (present(snow_depth)) snow_depth = 0.0

    if (present(gfrac)) gfrac = 0.0
    if (present(gtopo)) gtopo = 0.0
    if (present(grofi)) grofi = 0.0
    if (present(grofl)) grofl = 0.0
    if (present(ghflx)) ghflx = 0.0

    ! Allocate arrays

    allocate(orog_temp(params%g_grid_orog%nx, params%g_grid_orog%ny))
    allocate(alb_temp (params%g_grid%nx, params%g_grid%ny))
    allocate(if_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(vf_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(sif_temp (params%g_grid%nx, params%g_grid%ny))
    allocate(svf_temp (params%g_grid%nx, params%g_grid%ny))
    allocate(sd_temp  (params%g_grid%nx, params%g_grid%ny))

    if (params%gcm_smb) then
       allocate(gfrac_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
       allocate(gtopo_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
       allocate(grofi_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
       allocate(grofl_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
       allocate(ghflx_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
    endif

    if (GLC_DEBUG) then
       write(stdout,*) 'Upscale and splice the initial fields'
    endif

    ! Get initial fields from instances, splice together and return

    do i=1,params%ninstances

       call get_i_upscaled_fields(params%instances(i),   &
                                  orog_temp,   alb_temp, &
                                  if_temp,     vf_temp,  &
                                  sif_temp,    svf_temp, &
                                  sd_temp)

       if (present(orog)) &
            orog = splice_field(orog, orog_temp, params%instances(i)%frac_cov_orog, &
            params%cov_norm_orog)

       if (present(albedo)) &
            albedo = splice_field(albedo, alb_temp, params%instances(i)%frac_coverage, &
            params%cov_normalise)

       if (present(ice_frac)) &
            ice_frac = splice_field(ice_frac, if_temp, params%instances(i)%frac_coverage, &
            params%cov_normalise)

       if (present(veg_frac)) &
            veg_frac = splice_field(veg_frac, vf_temp, params%instances(i)%frac_coverage, &
            params%cov_normalise)

       if (present(snowice_frac)) &
            snowice_frac = splice_field(snowice_frac,sif_temp,params%instances(i)%frac_coverage, &
            params%cov_normalise)

       if (present(snowveg_frac)) &
            snowveg_frac = splice_field(snowveg_frac,svf_temp,params%instances(i)%frac_coverage, &
            params%cov_normalise)

       if (present(snow_depth)) &
            snow_depth = splice_field(snow_depth,sd_temp,params%instances(i)%frac_coverage, &
            params%cov_normalise)

       if (params%gcm_smb) then

!lipscomb - TO DO - These temp arrays are not currently upscaled correctly
          call get_i_upscaled_fields_gcm(params%instances(i), params%g_grid%nec,  &
                                         params%instances(i)%lgrid%size%pt(1),    &
                                         params%instances(i)%lgrid%size%pt(2),    &
                                         params%g_grid%nx,    params%g_grid%ny,   &
                                         gfrac_temp,          gtopo_temp,         &
                                         grofi_temp,          grofl_temp,         &
                                         ghflx_temp)

          do n = 1, params%g_grid%nec

             if (present(gfrac))    &
                gfrac(:,:,n) = splice_field(gfrac(:,:,n),                      &
                                            gfrac_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

             if (present(gtopo))    &
                gtopo(:,:,n) = splice_field(gtopo(:,:,n),                      &
                                            gtopo_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

             if (present(grofi))    &
                grofi(:,:,n) = splice_field(grofi(:,:,n),                      &
                                            grofi_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

             if (present(grofl))    &
                grofl(:,:,n) = splice_field(grofl(:,:,n),                      &
                                            grofl_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

             if (present(ghflx))    &
                ghflx(:,:,n) = splice_field(ghflx(:,:,n),                      &
                                            ghflx_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

          enddo  ! nec

       endif     ! gcm_smb

    end do

    ! Deallocate

    deallocate(orog_temp, alb_temp, if_temp, vf_temp, sif_temp, svf_temp,sd_temp)
    if (params%gcm_smb) deallocate(gfrac_temp, gtopo_temp, grofi_temp, grofl_temp, ghflx_temp)

    ! Sort out snow_model flag

    if (present(snow_model)) then
       snow_model = .false.
       do i=1, params%ninstances
          snow_model = (snow_model .or. glint_has_snow_model(params%instances(i)))
       end do
    end if

    ! Set output flag

    if (present(output_flag)) output_flag = .true.

  end subroutine initialise_glint

  !================================================================================

  subroutine glint(params,         time,            &
                   rawtemp,        rawprecip,       &
                   orog,                            &
                   zonwind,        merwind,         &
                   humid,          lwdown,          &
                   swdown,         airpress,        &
                   output_flag,                     &
                   orog_out,       albedo,          &
                   ice_frac,       veg_frac,        &
                   snowice_frac,   snowveg_frac,    &
                   snow_depth,                      &
                   water_in,       water_out,       &
                   total_water_in, total_water_out, &
                   ice_volume,     ice_tstep,       &
                   qsmb,           tsfc,            &
                   topo,           gfrac,           &
                   gtopo,          grofi,           &
                   grofl,          ghflx)

    !*FD Main Glimmer subroutine.
    !*FD
    !*FD This should be called daily or hourly, depending on
    !*FD the mass-balance scheme being used. It does all necessary 
    !*FD spatial and temporal averaging, and calls the dynamics 
    !*FD part of the model when required. 
    !*FD
    !*FD Input fields should be taken as means over the period since the last call.
    !*FD See the user documentation for more information.
    !*FD
    !*FD Note that the total ice volume returned is the total at the end of the time-step;
    !*FD the water fluxes are valid over the duration of the timestep. Thus the difference
    !*FD between \texttt{total\_water\_in} and \texttt{total\_water\_out} should be equal
    !*FD to the change in \texttt{ice\_volume}, after conversion between m$^3$ and kg.

    use glimmer_utils
    use glint_interp
    use glint_timestep
    use glimmer_log
    implicit none

    ! Subroutine argument declarations -------------------------------------------------------------

    type(glint_params),              intent(inout) :: params          !*FD parameters for this run
    integer,                         intent(in)    :: time            !*FD Current model time        (hours)
    real(rk),dimension(:,:),target,  intent(in)    :: rawtemp         !*FD Surface temperature field (deg C)
    real(rk),dimension(:,:),target,  intent(in)    :: rawprecip       !*FD Precipitation rate        (mm/s)
    real(rk),dimension(:,:),         intent(in)    :: orog            !*FD The large-scale orography (m)
    real(rk),dimension(:,:),optional,intent(in)    :: zonwind,merwind !*FD Zonal and meridional components 
                                                                      !*FD of the wind field         (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: humid           !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown          !*FD Downwelling longwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown          !*FD Downwelling shortwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress        !*FD surface air pressure (Pa)
    logical,                optional,intent(out)   :: output_flag     !*FD Set true if outputs set
    real(rk),dimension(:,:),optional,intent(inout) :: orog_out        !*FD The fed-back, output orography (m)
    real(rk),dimension(:,:),optional,intent(inout) :: albedo          !*FD surface albedo
    real(rk),dimension(:,:),optional,intent(inout) :: ice_frac        !*FD grid-box ice-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: veg_frac        !*FD grid-box veg-fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snowice_frac    !*FD grid-box snow-covered ice fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snowveg_frac    !*FD grid-box snow-covered veg fraction
    real(rk),dimension(:,:),optional,intent(inout) :: snow_depth      !*FD grid-box mean snow depth (m water equivalent)
    real(rk),dimension(:,:),optional,intent(inout) :: water_in        !*FD Input water flux          (mm)
    real(rk),dimension(:,:),optional,intent(inout) :: water_out       !*FD Output water flux         (mm)
    real(rk),               optional,intent(inout) :: total_water_in  !*FD Area-integrated water flux in (kg)
    real(rk),               optional,intent(inout) :: total_water_out !*FD Area-integrated water flux out (kg)
    real(rk),               optional,intent(inout) :: ice_volume      !*FD Total ice volume (m$^3$)
    logical,                optional,intent(out)   :: ice_tstep       !*FD Set when an ice-timestep has been done, and
                                                                      !*FD water balance information is available
    real(rk),dimension(:,:,:),optional,intent(in)    :: qsmb          ! flux of glacier ice (kg/m^2/s)
    real(rk),dimension(:,:,:),optional,intent(in)    :: tsfc          ! surface ground temperature (deg C)
    real(rk),dimension(:,:,:),optional,intent(in)    :: topo          ! surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(inout) :: gfrac         ! ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(inout) :: gtopo         ! surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(inout) :: grofi         ! ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(inout) :: grofl         ! liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(inout) :: ghflx         ! heat flux (W/m^2, positive down)

    ! Internal variables ----------------------------------------------------------------------------

    integer :: i, n
    real(rk),dimension(:,:),allocatable :: albedo_temp, if_temp, vf_temp, sif_temp, svf_temp,  &
                                           sd_temp, wout_temp, orog_out_temp, win_temp
    real(rk) :: twin_temp,twout_temp,icevol_temp
    type(output_flags) :: out_f
    logical :: icets
    character(250) :: message
    real(rk),dimension(size(rawprecip,1),size(rawprecip,2)),target :: anomprecip
    real(rk),dimension(size(rawtemp,1),  size(rawtemp,2)),  target :: anomtemp
    real(rk),dimension(:,:),pointer :: precip
    real(rk),dimension(:,:),pointer :: temp
    real(rk) :: yearfrac
    real(rk) :: timeyr   ! time in years
    integer :: j, ig, jg

    real(rk), dimension(:,:,:), allocatable ::   &
       gfrac_temp    ,&! gfrac for a single instance
       gtopo_temp    ,&! gtopo for a single instance
       grofi_temp    ,&! grofi for a single instance
       grofl_temp    ,&! grofl for a single instance
       ghflx_temp      ! ghflx for a single instance

    if (GLC_DEBUG) then
!       write (stdout,*) 'In subroutine glint, current time (hr) =', time
!       write (stdout,*) 'av_start_time =', params%av_start_time
!       write (stdout,*) 'next_av_start =', params%next_av_start
!       write (stdout,*) 'new_av =', params%new_av
!       write (stdout,*) 'tstep_mbal =', params%tstep_mbal
    endif

    ! Check we're expecting a call now --------------------------------------------------------------

    if (params%new_av) then
       if (time == params%next_av_start) then
          params%av_start_time = time
          params%new_av = .false.
       else
          write(message,*) 'Unexpected calling of GLINT at time ', time
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
    else
       if (mod(time-params%av_start_time,params%time_step)/=0) then
          write(message,*) 'Unexpected calling of GLINT at time ', time
          call write_log(message,GM_FATAL,__FILE__,__LINE__)
       end if
    end if

    ! Check input fields are correct ----------------------------------------------------------------

    call check_input_fields(params, humid, lwdown, swdown, airpress, zonwind, merwind)

    ! Reset output flag

    if (present(output_flag)) output_flag = .false.
    if (present(ice_tstep))   ice_tstep = .false.

    ! Sort out anomaly coupling

    if (params%anomaly_params%enabled) then
       yearfrac = real(mod(time,days_in_year),rk)/real(days_in_year,rk)
       call anomaly_calc(params%anomaly_params, yearfrac, rawtemp, rawprecip, anomtemp, anomprecip)
       precip => anomprecip
       temp   => anomtemp
    else
       precip => rawprecip
       temp   => rawtemp
    end if

    ! Do averaging and so on...

    call accumulate_averages(params,           &
                             temp,    precip,  &
                             zonwind, merwind, &
                             humid,   lwdown,  &
                             swdown,  airpress)

    if (params%gcm_smb) call accumulate_averages_gcm(params, qsmb, tsfc, topo)

    ! Increment step counter

    params%av_steps = params%av_steps + 1

    ! ---------------------------------------------------------
    ! If this is a mass balance timestep, prepare global fields, and do a timestep
    ! for each model instance
    ! ---------------------------------------------------------

    if (time-params%av_start_time+params%time_step > params%tstep_mbal) then

       write(message,*) &
            'Incomplete forcing of GLINT mass-balance time-step detected at time ', time
       call write_log(message,GM_FATAL,__FILE__,__LINE__)

    else if (time-params%av_start_time+params%time_step == params%tstep_mbal) then

       ! Set output_flag

       ! At present, outputs are done for each mass-balance timestep, since
       ! that involved least change to the code. However, it might be good
       ! to change the output to occur with user-specified frequency.

       if (present(output_flag)) output_flag = .true.

       ! Allocate output fields

       if (present(orog_out)) then
          allocate(orog_out_temp(size(orog_out,1),size(orog_out,2)))
       else
          allocate(orog_out_temp(params%g_grid_orog%nx, params%g_grid_orog%ny))
       end if
       allocate(albedo_temp(size(orog,1),size(orog,2)))
       allocate(if_temp(size(orog,1),size(orog,2)))
       allocate(vf_temp(size(orog,1),size(orog,2)))
       allocate(sif_temp(size(orog,1),size(orog,2)))
       allocate(svf_temp(size(orog,1),size(orog,2)))
       allocate(sd_temp(size(orog,1),size(orog,2)))
       allocate(wout_temp(size(orog,1),size(orog,2)))
       allocate(win_temp(size(orog,1),size(orog,2)))
       if (params%gcm_smb) then
          allocate(gfrac_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
          allocate(gtopo_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
          allocate(grofi_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
          allocate(grofl_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
          allocate(ghflx_temp(params%g_grid%nx,params%g_grid%ny,params%g_grid%nec))
       endif

       ! Populate output flag derived type

       call populate_output_flags(out_f,                           &
                                  orog_out,       albedo,          &
                                  ice_frac,       veg_frac,        &
                                  snowice_frac,   snowveg_frac,    &
                                  snow_depth,                      &
                                  water_in,       water_out,       &
                                  total_water_in, total_water_out, &
                                  ice_volume)

       ! Zero outputs if present

       if (present(orog_out))        orog_out        = 0.0
       if (present(albedo))          albedo          = 0.0
       if (present(ice_frac))        ice_frac        = 0.0
       if (present(veg_frac))        veg_frac        = 0.0
       if (present(snowice_frac))    snowice_frac    = 0.0
       if (present(snowveg_frac))    snowveg_frac    = 0.0
       if (present(snow_depth))      snow_depth      = 0.0
       if (present(water_out))       water_out       = 0.0
       if (present(water_in))        water_in        = 0.0
       if (present(total_water_in))  total_water_in  = 0.0
       if (present(total_water_out)) total_water_out = 0.0
       if (present(ice_volume))      ice_volume      = 0.0

       if (present(gfrac)) gfrac = 0.0
       if (present(gtopo)) gtopo = 0.0
       if (present(grofi)) grofi = 0.0
       if (present(grofl)) grofl = 0.0
       if (present(ghflx)) ghflx = 0.0

       ! Calculate averages by dividing by number of steps elapsed
       ! since last model timestep.

       call calculate_averages(params)
       if (params%gcm_smb) call calculate_averages_gcm(params)

       ! Calculate total accumulated precipitation - multiply
       ! by time since last model timestep

       params%g_av_precip = params%g_av_precip*params%tstep_mbal*hours2seconds

       ! Calculate temperature half-range

       params%g_temp_range=(params%g_max_temp-params%g_min_temp)/2.0

       if (GLC_DEBUG) then
          i = itest
          j = jjtest
          write(stdout,*) 'Take a mass balance timestep, time (hr) =', time
          write(stdout,*) 'av_steps =', real(params%av_steps,rk)
          write(stdout,*) 'tstep_mbal (hr) =', params%tstep_mbal
          write(stdout,*) 'i, j =', i, j
          if (params%gcm_smb) then
             do n = 1, params%g_grid%nec
                write (stdout,*) ' '
                write (stdout,*) 'n =', n
                write (stdout,*) 'g_av_qsmb (kg m-2 s-1) =', params%g_av_qsmb(i,j,n)
                write (stdout,*) 'g_av_tsfc (Celsius) =',    params%g_av_tsfc(i,j,n)
                write (stdout,*) 'g_av_topo (m) =',          params%g_av_topo(i,j,n)
             enddo
          endif
          write(stdout,*) 'call glint_i_tstep'
       endif

       ! Calculate total surface mass balance - multiply by time since last model timestep
       ! Note on units: We want g_av_qsmb to have units of m per time step.
       ! Divide by 1000 to convert from mm to m.
       ! Multiply by 3600 to convert from 1/s to 1/hr.  (tstep_mbal has units of hours)

       if (params%gcm_smb) &
          params%g_av_qsmb(:,:,:) = params%g_av_qsmb(:,:,:) * params%tstep_mbal * hours2seconds / 1000._rk

       ! Do a timestep for each instance

       do i=1,params%ninstances

          if (params%gcm_smb) then

             !lipscomb - TO DO - Make some of these arguments optional?

             call glint_i_tstep(time,      params%instances(i),          &
                  params%g_av_temp,        params%g_temp_range,          &
                  params%g_av_precip,      params%g_av_zonwind,          &
                  params%g_av_merwind,     params%g_av_humid,            &
                  params%g_av_lwdown,      params%g_av_swdown,           &
                  params%g_av_airpress,                                  &
                  orog,                    orog_out_temp,                &
                  albedo_temp,             if_temp,                      &
                  vf_temp,                 sif_temp,                     &
                  svf_temp,                sd_temp,                      &
                  win_temp,                wout_temp,                    &
                  twin_temp,               twout_temp,                   &
                  icevol_temp,             out_f,                        &
                  .true.,                  icets,                        &
                  gcm_smb_in = params%gcm_smb,                           &
                  qsmb_g = params%g_av_qsmb, tsfc_g = params%g_av_tsfc,  &
                  topo_g = params%g_av_topo, gmask  = params%g_grid%mask,&
                  gfrac  = gfrac_temp,       gtopo  = gtopo_temp,        &
                  grofi  = grofi_temp,       grofl  = grofl_temp,        &
                  ghflx  = ghflx_temp  )

          else 

             call glint_i_tstep(time,                    params%instances(i),       &
                                params%g_av_temp,        params%g_temp_range,       &
                                params%g_av_precip,      params%g_av_zonwind,       &
                                params%g_av_merwind,     params%g_av_humid,         &
                                params%g_av_lwdown,      params%g_av_swdown,        &
                                params%g_av_airpress,                               &
                                orog,                    orog_out_temp,             &
                                albedo_temp,             if_temp,                   &
                                vf_temp,                 sif_temp,                  &
                                svf_temp,                sd_temp,                   &
                                win_temp,                wout_temp,                 &
                                twin_temp,               twout_temp,                &
                                icevol_temp,             out_f,                     &
                                .true.,                  icets)

           endif

          if (GLC_DEBUG) then
             write(stdout,*) 'Finished glc_glint_ice tstep, instance =', i
             write(stdout,*) 'Upscale fields to global grid'
          endif

          ! Add this contribution to the output orography

          if (present(orog_out)) orog_out=splice_field(orog_out,orog_out_temp, &
               params%instances(i)%frac_cov_orog,params%cov_norm_orog)

          if (present(albedo)) albedo=splice_field(albedo,albedo_temp, &
               params%instances(i)%frac_coverage,params%cov_normalise)

          if (present(ice_frac)) ice_frac=splice_field(ice_frac,if_temp, &
               params%instances(i)%frac_coverage,params%cov_normalise)

          if (present(veg_frac)) veg_frac=splice_field(veg_frac,vf_temp, &
               params%instances(i)%frac_coverage,params%cov_normalise)

          if (present(snowice_frac))snowice_frac=splice_field(snowice_frac,sif_temp, &
               params%instances(i)%frac_coverage,params%cov_normalise)

          if (present(snowveg_frac)) snowveg_frac=splice_field(snowveg_frac, &
               svf_temp,params%instances(i)%frac_coverage, params%cov_normalise)

          if (present(snow_depth)) snow_depth=splice_field(snow_depth, &
               sd_temp,params%instances(i)%frac_coverage,params%cov_normalise)

          if (present(water_in)) water_in=splice_field(water_in,win_temp, &
               params%instances(i)%frac_coverage,params%cov_normalise)

          if (present(water_out)) water_out=splice_field(water_out, &
               wout_temp, params%instances(i)%frac_coverage,params%cov_normalise)

          ! Add total water variables to running totals

          if (present(total_water_in))  total_water_in  = total_water_in  + twin_temp
          if (present(total_water_out)) total_water_out = total_water_out + twout_temp
          if (present(ice_volume))      ice_volume      = ice_volume      + icevol_temp

          ! Set flag
          if (present(ice_tstep)) then
             ice_tstep = (ice_tstep .or. icets)
          end if

          ! Upscale the output to elevation classes on the global grid

          if (params%gcm_smb) then

             call get_i_upscaled_fields_gcm(params%instances(i), params%g_grid%nec, &
                                            params%instances(i)%lgrid%size%pt(1),   &
                                            params%instances(i)%lgrid%size%pt(2),   &
                                            params%g_grid%nx,    params%g_grid%ny,  &
                                            gfrac_temp,          gtopo_temp,        &
                                            grofi_temp,          grofl_temp,        &
                                            ghflx_temp )

             if (GLC_DEBUG) then
                ig = itest
                jg = jjtest
                write(stdout,*) ' '
                write(stdout,*) 'After upscaling:'
                do n = 1, params%g_grid%nec
                  write(stdout,*) ' '
                  write(stdout,*) 'n =', n
                  write(stdout,*) 'gfrac(n) =', gfrac(ig,jg,n)
                  write(stdout,*) 'gtopo(n) =', gtopo(ig,jg,n)
!!                  write(stdout,*) 'grofi(n) =', grofi(ig,jg,n)
!!                  write(stdout,*) 'grofl(n) =', grofl(ig,jg,n)
!!                  write(stdout,*) 'ghflx(n) =', ghflx(ig,jg,n)
                enddo
             endif

          ! Add this contribution to the global output

             do n = 1, params%g_grid%nec

                gfrac(:,:,n) = splice_field(gfrac(:,:,n),                      &
                                            gfrac_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

                gtopo(:,:,n) = splice_field(gtopo(:,:,n),                      &
                                            gtopo_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

                grofi(:,:,n) = splice_field(grofi(:,:,n),                      &
                                            grofi_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

                grofl(:,:,n) = splice_field(grofl(:,:,n),                      &
                                            grofl_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

                ghflx(:,:,n) = splice_field(ghflx(:,:,n),                      &
                                            ghflx_temp(:,:,n),                 &
                                            params%instances(i)%frac_coverage, &
                                            params%cov_normalise)

             enddo   ! nec

          endif   ! gcm_smb

          ! write ice sheet diagnostics
          if (mod(params%instances(i)%model%numerics%timecounter,  &
                  params%instances(i)%model%numerics%ndiag)==0)  then
             timeyr = time / (days_in_year*24.d0) 
             if (GLC_DEBUG) then
                write(stdout,*) 'Write diagnostics, time (yr)=', timeyr     
             endif
             call glide_write_diag(params%instances(i)%model, timeyr, idiag, jdiag)
          endif

       enddo

       ! Scale output water fluxes to be in mm/s

       if (present(water_in))  water_in  = water_in/ &
                                           (params%tstep_mbal*hours2seconds)

       if (present(water_out)) water_out = water_out/ &
                                           (params%tstep_mbal*hours2seconds)

       ! ---------------------------------------------------------
       ! Reset averaging fields, flags and counters
       ! ---------------------------------------------------------

       params%g_av_temp    = 0.0
       params%g_av_precip  = 0.0
       params%g_av_zonwind = 0.0
       params%g_av_merwind = 0.0
       params%g_av_humid   = 0.0
       params%g_av_lwdown  = 0.0
       params%g_av_swdown  = 0.0
       params%g_av_airpress = 0.0
       params%g_temp_range = 0.0
       params%g_max_temp   = -1000.0
       params%g_min_temp   = 1000.0
       if (params%gcm_smb) then
          params%g_av_qsmb    = 0.0
          params%g_av_tsfc    = 0.0
          params%g_av_topo    = 0.0
       endif

       params%av_steps      = 0
       params%new_av        = .true.
       params%next_av_start = time+params%time_step

       deallocate(albedo_temp,if_temp,vf_temp,sif_temp,svf_temp,sd_temp,wout_temp,win_temp,orog_out_temp)
       if (params%gcm_smb) deallocate(gfrac_temp, gtopo_temp, grofi_temp, grofl_temp, ghflx_temp)

    endif

  end subroutine glint

  !===================================================================

  subroutine end_glint(params,close_logfile)

    !*FD perform tidying-up operations for glimmer
    use glint_initialise
    use glimmer_log
    implicit none

    type(glint_params),intent(inout) :: params          ! parameters for this run
    logical, intent(in), optional    :: close_logfile   ! if true, then close the log file
                                                        ! (GCM may do this elsewhere)                                  

    integer i

    ! end individual instances

    do i=1,params%ninstances
       call glint_i_end(params%instances(i))
    enddo

    if (present(close_logfile)) then
       if (close_logfile) call close_log
    else
          call close_log
    endif

  end subroutine end_glint

  !=====================================================

  integer function glint_coverage_map(params, coverage, cov_orog)

    !*FD Retrieve ice model fractional 
    !*FD coverage map. This function is provided so that glimmer may
    !*FD be restructured without altering the interface.
    !*RV Three return values are possible:
    !*RV \begin{description}
    !*RV \item[0 ---] Successful return
    !*RV \item[1 ---] Coverage map not calculated yet (fail)
    !*RV \item[2 ---] Coverage array is the wrong size (fail)
    !*RV \end{description}

    implicit none

    type(glint_params),intent(in) :: params         !*FD ice model parameters
    real(rk),dimension(:,:),intent(out) :: coverage !*FD array to hold coverage map
    real(rk),dimension(:,:),intent(out) :: cov_orog !*FD Orography coverage

    if (.not.params%coverage_calculated) then
       glint_coverage_map = 1
       return
    endif

    if ( size(coverage,1) /= params%g_grid%nx .or. &
         size(coverage,2) /= params%g_grid%ny) then
       glint_coverage_map = 2
       return
    endif

    glint_coverage_map = 0
    coverage = params%total_coverage
    cov_orog = params%total_cov_orog

  end function glint_coverage_map

  !=====================================================

  !----------------------------------------------------------------------
  ! PRIVATE INTERNAL GLIMMER SUBROUTINES FOLLOW.............
  !----------------------------------------------------------------------

  subroutine glint_allocate_arrays(params)

    !*FD allocates glimmer arrays

    implicit none

    type(glint_params),intent(inout) :: params !*FD ice model parameters

    allocate(params%g_av_precip (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_temp   (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_max_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_min_temp  (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_temp_range(params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_zonwind(params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_merwind(params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_humid  (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_lwdown (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_swdown (params%g_grid%nx, params%g_grid%ny))
    allocate(params%g_av_airpress(params%g_grid%nx,params%g_grid%ny))

    allocate(params%total_coverage(params%g_grid%nx, params%g_grid%ny))
    allocate(params%cov_normalise (params%g_grid%nx, params%g_grid%ny))
    allocate(params%total_cov_orog(params%g_grid_orog%nx, params%g_grid_orog%ny))
    allocate(params%cov_norm_orog (params%g_grid_orog%nx, params%g_grid_orog%ny))

  end subroutine glint_allocate_arrays

  !========================================================

  subroutine glint_allocate_arrays_gcm(params)

    !*FD allocates glimmer arrays for GCM input fields

    implicit none

    type(glint_params),intent(inout) :: params !*FD ice model parameters

    ! input fields from GCM

    allocate(params%g_av_qsmb (params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(params%g_av_tsfc (params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))
    allocate(params%g_av_topo (params%g_grid%nx, params%g_grid%ny, params%g_grid%nec))

  end subroutine glint_allocate_arrays_gcm

  !========================================================

  function splice_field(global, local, coverage, normalise)

    !*FD Splices an upscaled field into a global field

    real(rk),dimension(:,:),intent(in) :: global    !*FD Field to receive the splice
    real(rk),dimension(:,:),intent(in) :: local     !*FD The field to be spliced in
    real(rk),dimension(:,:),intent(in) :: coverage  !*FD The coverage fraction
    real(rk),dimension(:,:),intent(in) :: normalise !*FD The normalisation field

    real(rk),dimension(size(global,1),size(global,2)) :: splice_field

    where (coverage==0.0)
       splice_field=global
    elsewhere
       splice_field=(global*(1-coverage/normalise))+(local*coverage/normalise)
    end where

  end function splice_field

  !========================================================

  subroutine glint_readconfig(config, ninstances, fnames, infnames)

    !*FD Determine whether a given config file is a
    !*FD top-level glint config file, and return parameters
    !*FD accordingly.

    use glimmer_config
    use glimmer_log
    implicit none

    ! Arguments -------------------------------------------

    type(ConfigSection),      pointer :: config !*FD structure holding sections of configuration file
    integer,              intent(out) :: ninstances !*FD Number of instances to create
    character(fname_length),dimension(:),pointer :: fnames !*FD list of filenames (output)
    character(fname_length),dimension(:) :: infnames !*FD list of filenames (input)

    ! Internal variables ----------------------------------

    type(ConfigSection), pointer :: section
    character(len=100) :: message
    integer :: i

    if (associated(fnames)) nullify(fnames)

    call GetSection(config,section,'GLINT')
    if (associated(section)) then
       call GetValue(section,'n_instance',ninstances)
       allocate(fnames(ninstances))
       do i=1,ninstances
          call GetSection(section%next,section,'GLINT instance')
          if (.not.associated(section)) then
             write(message,*) 'Must specify ',ninstances,' instance config files'
             call write_log(message,GM_FATAL,__FILE__,__LINE__)
          end if
          call GetValue(section,'name',fnames(i))
       end do
    else
       ninstances=1
       allocate(fnames(1))
       fnames=infnames
    end if

    ! Print some configuration information

!!$    call write_log('GLINT global')
!!$    call write_log('------------')
!!$    write(message,*) 'number of instances :',params%ninstances
!!$    call write_log(message)
!!$    call write_log('')

  end subroutine glint_readconfig

  !========================================================

  subroutine calc_bounds(lon, lat, lonb, latb)

    !*FD Calculates the boundaries between
    !*FD global grid-boxes. Note that we assume that the boundaries lie 
    !*FD half-way between the 
    !*FD points, both latitudinally and longitudinally, although 
    !*FD this isn't strictly true for a Gaussian grid.

    use glimmer_map_trans, only: loncorrect

    implicit none

    real(rk),dimension(:),intent(in) :: lon,lat    !*FD locations of global grid-points (degrees)
    real(rk),dimension(:),intent(out) :: lonb,latb !*FD boundaries of grid-boxes (degrees)

    real(rk) :: dlon

    integer :: nxg,nyg,i,j

    nxg=size(lon) ; nyg=size(lat)

    ! Latitudes first - we assume the boundaries of the first and 
    ! last boxes coincide with the poles. Not sure how to
    ! handle it if they don't...

    latb(1)=90.0
    latb(nyg+1)=-90.0

    do j=2,nyg
       latb(j)=lat(j-1)-(lat(j-1)-lat(j))/2.0
    enddo

    ! Longitudes

    if (lon(1)<lon(nxg)) then
       dlon=lon(1)-lon(nxg)+360.0
    else
       dlon=lon(1)-lon(nxg)
    endif
    lonb(1)=lon(nxg)+dlon/2
    lonb(1)=loncorrect(lonb(1),0.0_rk)      

    lonb(nxg+1)=lonb(1)

    do i=2,nxg
       if (lon(i)<lon(i-1)) then
          dlon=lon(i)-lon(i-1)+360.0
       else
          dlon=lon(i)-lon(i-1)
       endif
       lonb(i)=lon(i-1)+dlon/2
       lonb(i)=loncorrect(lonb(i),0.0_rk)      
    enddo

  end subroutine calc_bounds


  !========================================================

  integer function check_mbts(timesteps)

    !*FD Checks to see that all mass-balance time-steps are
    !*FD the same. Flags a fatal error if not, else assigns that
    !*FD value to the output

    use glimmer_log

    implicit none

    integer,dimension(:) :: timesteps !*FD Array of mass-balance timsteps

    integer :: n,i

    n=size(timesteps)
    if (n==0) then
       check_mbts=0
       return
    endif

    check_mbts=timesteps(1)

    do i=2,n
       if (timesteps(i)/=check_mbts) then
          call write_log('All instances must have the same mass-balance and ice timesteps', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    enddo

  end function check_mbts

  !========================================================

  subroutine check_init_args(orog_lats, orog_longs, orog_latb, orog_lonb)

    !*FD Checks which combination arguments have been supplied to
    !*FD define the global grid, and rejects unsuitable combinations

    use glimmer_log

    real(rk),dimension(:),optional,intent(in) :: orog_lats 
    real(rk),dimension(:),optional,intent(in) :: orog_longs 
    real(rk),dimension(:),optional,intent(in) :: orog_latb
    real(rk),dimension(:),optional,intent(in) :: orog_lonb 

    integer :: args
    integer,dimension(5) :: allowed=(/0,3,7,11,15/)

    args=0

    if (present(orog_lats))  args=args+1
    if (present(orog_longs)) args=args+2
    if (present(orog_latb))  args=args+4
    if (present(orog_lonb))  args=args+8

    if (.not.any(args==allowed)) then
       call write_log('Unexpected combination of arguments to initialise_glint', &
            GM_FATAL,__FILE__,__LINE__)
    end if

  end subroutine check_init_args

  !========================================================

  subroutine check_input_fields(params, humid, lwdown, swdown, airpress, zonwind, merwind)

    use glimmer_log

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run
    real(rk),dimension(:,:),optional,intent(in)    :: humid    !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown   !*FD Downwelling longwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown   !*FD Downwelling shortwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress !*FD surface air pressure (Pa)
    real(rk),dimension(:,:),optional,intent(in)    :: zonwind  !*FD Zonal component of the wind field (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: merwind  !*FD Meridional component of the wind field (m/s)

    if (params%enmabal) then
       if (.not.(present(humid).and.present(lwdown).and. &
            present(swdown).and.present(airpress).and. &
            present(zonwind).and.present(merwind))) &
            call write_log('Necessary fields not supplied for Energy Balance Mass Balance model',GM_FATAL, &
            __FILE__,__LINE__)
    end if

    if (params%need_winds) then
       if (.not.(present(zonwind).and.present(merwind))) &
          call write_log('Need to supply zonal and meridional wind fields to GLINT',GM_FATAL, &
            __FILE__,__LINE__)
    end if

  end subroutine check_input_fields

  !========================================================

  subroutine accumulate_averages(params, temp, precip, zonwind, merwind, humid, lwdown, swdown, airpress)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run
    real(rk),dimension(:,:),         intent(in)    :: temp     !*FD Surface temperature field (celsius)
    real(rk),dimension(:,:),         intent(in)    :: precip   !*FD Precipitation rate        (mm/s)
    real(rk),dimension(:,:),optional,intent(in)    :: zonwind  !*FD Zonal component of the wind field (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: merwind  !*FD Meridional component of the wind field (m/s)
    real(rk),dimension(:,:),optional,intent(in)    :: humid    !*FD Surface humidity (%)
    real(rk),dimension(:,:),optional,intent(in)    :: lwdown   !*FD Downwelling longwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: swdown   !*FD Downwelling shortwave (W/m$^2$)
    real(rk),dimension(:,:),optional,intent(in)    :: airpress !*FD surface air pressure (Pa)

    params%g_av_temp    = params%g_av_temp    + temp
    params%g_av_precip  = params%g_av_precip  + precip

    if (params%need_winds) params%g_av_zonwind = params%g_av_zonwind + zonwind
    if (params%need_winds) params%g_av_merwind = params%g_av_merwind + merwind

    if (params%enmabal) then
       params%g_av_humid    = params%g_av_humid    + humid
       params%g_av_lwdown   = params%g_av_lwdown   + lwdown
       params%g_av_swdown   = params%g_av_swdown   + swdown
       params%g_av_airpress = params%g_av_airpress + airpress
    endif

    ! Ranges of temperature

    where (temp > params%g_max_temp) params%g_max_temp = temp
    where (temp < params%g_min_temp) params%g_min_temp = temp

  end subroutine accumulate_averages

  !========================================================

  subroutine accumulate_averages_gcm(params, qsmb, tsfc, topo)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run
    real(rk),dimension(:,:,:),optional,intent(in)  :: qsmb     ! flux of glacier ice (kg/m^2/s)
    real(rk),dimension(:,:,:),optional,intent(in)  :: tsfc     ! surface ground temperature (C)
    real(rk),dimension(:,:,:),optional,intent(in)  :: topo     ! surface elevation (m)

    if (present(qsmb)) params%g_av_qsmb(:,:,:) = params%g_av_qsmb(:,:,:) + qsmb(:,:,:)
    if (present(tsfc)) params%g_av_tsfc(:,:,:) = params%g_av_tsfc(:,:,:) + tsfc(:,:,:)
    if (present(topo)) params%g_av_topo(:,:,:) = params%g_av_topo(:,:,:) + topo(:,:,:)

!lipscomb - TO DO - No need to accumulate topo?

  end subroutine accumulate_averages_gcm

  !========================================================

  subroutine calculate_averages(params)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run

    params%g_av_temp    = params%g_av_temp   /real(params%av_steps)
    params%g_av_precip  = params%g_av_precip /real(params%av_steps)
    if (params%need_winds) params%g_av_zonwind = params%g_av_zonwind/real(params%av_steps)
    if (params%need_winds) params%g_av_merwind = params%g_av_merwind/real(params%av_steps)
    if (params%enmabal) then
       params%g_av_humid    = params%g_av_humid   /real(params%av_steps)
       params%g_av_lwdown   = params%g_av_lwdown  /real(params%av_steps)
       params%g_av_swdown   = params%g_av_swdown  /real(params%av_steps)
       params%g_av_airpress = params%g_av_airpress/real(params%av_steps)
    endif

  end subroutine calculate_averages

  !========================================================

  subroutine calculate_averages_gcm(params)

    type(glint_params),              intent(inout) :: params   !*FD parameters for this run

!lipscomb - TO DO - Do not average topo?  Remove 'rk' here?
    params%g_av_qsmb(:,:,:) = params%g_av_qsmb(:,:,:) / real(params%av_steps,rk)
    params%g_av_tsfc(:,:,:) = params%g_av_tsfc(:,:,:) / real(params%av_steps,rk)
    params%g_av_topo(:,:,:) = params%g_av_topo(:,:,:) / real(params%av_steps,rk)

  end subroutine calculate_averages_gcm

  !========================================================

end module glint_main

