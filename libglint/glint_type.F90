! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_type.f90 - part of the Glimmer-CISM ice model      + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010
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

#define NCO outfile%nc
#define NCI infile%nc

module glint_type

  !*FD contains type definitions for GLINT

  use glimmer_global
  use glint_interp
  use glide_types
  use glint_mbal_coupling

  implicit none

  type glint_instance

     !*FD Derived type holding information about ice model instance. 

     type(coordsystem_type)           :: lgrid              !*FD Local grid for interfacing with glide
     type(downscale)                  :: downs              !*FD Downscaling parameters.
     type(upscale)                    :: ups                !*FD Upscaling parameters
     type(upscale)                    :: ups_orog           !*FD Upscaling parameters for orography (to cope
                                                            !*FD with need to convert to spectral form).
     type(glide_global_type)          :: model              !*FD The instance and all its arrays.
     character(fname_length)          :: paramfile          !*FD The name of the configuration file.
     integer                          :: ice_tstep          !*FD Ice timestep in hours
     integer                          :: mbal_tstep         !*FD Mass-balance timestep in hours
     integer                          :: mbal_accum_time    !*FD Accumulation time for mass-balance (hours)
                                                            !*FD (defaults to ice time-step)
     integer                          :: ice_tstep_multiply=1 !*FD Ice time multiplier (non-dimensional)
     integer                          :: n_icetstep         !*FD Number of ice time-steps per mass-balance accumulation
     real(rk)                         :: glide_time         !*FD Time as seen by glide (years)
     integer                          :: next_time          !*FD The next time we expect to be called (hours)

     ! Climate inputs from global model --------------------------

     real(sp),dimension(:,:),pointer :: artm        => null() !*FD Annual mean air temperature
     real(sp),dimension(:,:),pointer :: arng        => null() !*FD Annual air temperature half-range
     real(sp),dimension(:,:),pointer :: prcp        => null() !*FD Precipitation (mm or m)
     real(sp),dimension(:,:),pointer :: snowd       => null() !*FD Snow depth (m)
     real(sp),dimension(:,:),pointer :: siced       => null() !*FD Superimposed ice depth (m)
     real(rk),dimension(:,:),pointer :: xwind       => null() !*FD $x$-component of surface winds (m/s)
     real(rk),dimension(:,:),pointer :: ywind       => null() !*FD $y$-component of surface winds (m/s)
     real(rk),dimension(:,:),pointer :: humid       => null() !*FD Surface humidity (%)
     real(rk),dimension(:,:),pointer :: lwdown      => null() !*FD Downwelling longwave (W/m^2)
     real(rk),dimension(:,:),pointer :: swdown      => null() !*FD Downwelling shortwave (W/m^2)
     real(rk),dimension(:,:),pointer :: airpress    => null() !*FD Surface air pressure (Pa)
     real(dp),dimension(:,:),pointer :: global_orog => null() !*FD Global orography (m)
     real(sp),dimension(:,:),pointer :: local_orog  => null() !*FD Local orography (m)

     ! Locally calculated climate/mass-balance fields ------------

     real(sp),dimension(:,:),pointer :: ablt => null() !*FD Annual ablation
     real(sp),dimension(:,:),pointer :: acab => null() !*FD Annual mass-balance

     ! Arrays to accumulate mass-balance quantities --------------

     type(glint_mbc) :: mbal_accum

     ! Fractional coverage information ---------------------------

     real(rk) ,dimension(:,:),pointer :: frac_coverage => null() 
     !*FD Fractional coverage of each global gridbox by the projected grid.
     real(rk) ,dimension(:,:),pointer :: frac_cov_orog => null() 
     !*FD Fractional coverage of each global gridbox by the projected grid (orography).

     ! Output masking --------------------------------------------

     integer, dimension(:,:),pointer :: out_mask => null() 

     !*FD Array indicating whether a point should be considered or ignored 
     !*FD when upscaling data for output. 1 means use, 0 means ignore.

     ! Climate options -------------------------------------------

     integer :: whichacab = 1

     !*FD Which mass-balance scheme: 
     !*FD \begin{description}
     !*FD \item[0] Receive surface mass balance from climate model
     !*FD \item[1] PDD mass-balance model
     !*FD \item[2] Accumulation only 
     !*FD \item[3] RAPID energy balance model
     !*FD \end{description}

     integer :: whichprecip = 1

     !*FD Source of precipitation:
     !*FD \begin{description}
     !*FD \item[1] Use large-scale precip as is.
     !*FD \item[2] Use parameterization of \emph{Roe and Lindzen} 
     !*FD \end{description}

     integer :: use_mpint = 0
   
     !*FD Flag to control if mean-preserving interpolation is used

     ! Climate parameters ----------------------------------------------------------

     real(sp) :: ice_albedo   =   0.4 !*FD Ice albedo. (fraction)
     real(sp) :: lapse_rate   =   8.0 !*FD Uniform lapse rate in deg C/km 
     !*FD (N.B. This should be \emph{positive} for 
     !*FD temperature falling with height!)
     real(sp) :: data_lapse_rate = 8.0 !*FD Implied lapse rate in large-scale data (used for
     !*FD tuning). Set equal to lapse\_rate if not supplied.

     ! Counter for averaging temperature input --------------------------------------

     integer  :: av_count = 0 !*FD Counter for averaging temperature input

     ! Pointers to file input and output

     type(glimmer_nc_output),pointer :: out_first => null() !*FD first element of linked list defining netCDF outputs
     type(glimmer_nc_input), pointer :: in_first => null()  !*FD first element of linked list defining netCDF inputs

  end type glint_instance

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type output_flags
     !*FD A derived type used internally to communicate the outputs which need
     !*FD to be upscaled, thus avoiding unnecessary calculation

     logical :: orog         !*FD Set if we need to upscale the orography
     logical :: albedo       !*FD Set if we need to upscale the albedo
     logical :: ice_frac     !*FD Set if we need to upscale the ice fraction
     logical :: veg_frac     !*FD Set if we need to upscale the veg fraction
     logical :: snowice_frac !*FD Set if we need to upscale the snow-covered ice fraction
     logical :: snowveg_frac !*FD Set if we need to upscale the snow-covered veg fraction
     logical :: snow_depth   !*FD Set if we need to upscale the snow depth
     logical :: water_in     !*FD Set if we need to upscale the input water flux
     logical :: water_out    !*FD Set if we need to upscale the output water flux
     logical :: total_win    !*FD Set if we need to sum the total water taken up by ice sheet
     logical :: total_wout   !*FD Set if we need to sum the total ablation by the ice sheet
     logical :: ice_vol      !*FD Set if we need to calculate the total ice volume
  end type output_flags

contains

  subroutine glint_i_allocate(instance,nxg,nyg,nxgo,nygo)

    !*FD Allocate top-level arrays in
    !*FD the model instance, and ice model arrays.

    implicit none

    type(glint_instance),intent(inout) :: instance  !*FD Instance whose elements are to be allocated.
    integer,             intent(in)    :: nxg       !*FD Longitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nyg       !*FD Latitudinal size of global grid (grid-points).
    integer,             intent(in)    :: nxgo      !*FD Longitudinal size of global orog grid (grid-points).
    integer,             intent(in)    :: nygo      !*FD Latitudinal size of global orog grid (grid-points).

    integer ewn,nsn

    ewn=get_ewn(instance%model)
    nsn=get_nsn(instance%model)

    ! First deallocate if necessary
    ! Downscaled global arrays

    if (associated(instance%artm))          deallocate(instance%artm)
    if (associated(instance%arng))          deallocate(instance%arng)
    if (associated(instance%prcp))          deallocate(instance%prcp)
    if (associated(instance%snowd))         deallocate(instance%snowd)
    if (associated(instance%siced))         deallocate(instance%siced)
    if (associated(instance%xwind))         deallocate(instance%xwind)
    if (associated(instance%ywind))         deallocate(instance%ywind)
    if (associated(instance%humid))         deallocate(instance%humid)
    if (associated(instance%lwdown))        deallocate(instance%lwdown)
    if (associated(instance%swdown))        deallocate(instance%swdown)
    if (associated(instance%airpress))      deallocate(instance%airpress)
    if (associated(instance%global_orog))   deallocate(instance%global_orog) 
    if (associated(instance%local_orog))    deallocate(instance%local_orog)

    ! Local climate arrays

    if (associated(instance%ablt))          deallocate(instance%ablt)
    if (associated(instance%acab))          deallocate(instance%acab)

    ! Fractional coverage

    if (associated(instance%frac_coverage)) deallocate(instance%frac_coverage)
    if (associated(instance%frac_cov_orog)) deallocate(instance%frac_cov_orog)

    ! Output mask

    if (associated(instance%out_mask))      deallocate(instance%out_mask)

    ! Then reallocate and zero...
    ! Global input fields

    allocate(instance%artm(ewn,nsn));        instance%artm        = 0.0
    allocate(instance%arng(ewn,nsn));        instance%arng        = 0.0
    allocate(instance%prcp(ewn,nsn));        instance%prcp        = 0.0
    allocate(instance%snowd(ewn,nsn));       instance%snowd       = 0.0
    allocate(instance%siced(ewn,nsn));       instance%siced       = 0.0
    allocate(instance%xwind(ewn,nsn));       instance%xwind       = 0.0
    allocate(instance%ywind(ewn,nsn));       instance%ywind       = 0.0
    allocate(instance%humid(ewn,nsn));       instance%humid       = 0.0
    allocate(instance%lwdown(ewn,nsn));      instance%lwdown      = 0.0
    allocate(instance%swdown(ewn,nsn));      instance%swdown      = 0.0
    allocate(instance%airpress(ewn,nsn));    instance%airpress    = 0.0
    allocate(instance%global_orog(ewn,nsn)); instance%global_orog = 0.0
    allocate(instance%local_orog(ewn,nsn));  instance%local_orog  = 0.0

    ! Local fields

    allocate(instance%ablt(ewn,nsn)); instance%ablt = 0.0
    allocate(instance%acab(ewn,nsn)); instance%acab = 0.0

    ! Fractional coverage map

    allocate(instance%frac_coverage(nxg,nyg)); instance%frac_coverage = 0.0
    allocate(instance%frac_cov_orog(nxgo,nygo)); instance%frac_cov_orog = 0.0

    ! Output mask

    allocate(instance%out_mask(ewn,nsn)); instance%out_mask = 1

  end subroutine glint_i_allocate

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_readconfig(instance,config)

    !*FD read glint configuration

    use glimmer_config
    use glimmer_log
    use glint_constants, only: years2hours

    implicit none

    ! Arguments

    type(ConfigSection), pointer       :: config      !*FD structure holding sections of configuration file   
    type(glint_instance),intent(inout) :: instance    !*FD The instance being initialised.

    ! Internals

    type(ConfigSection), pointer :: section
    real(rk) :: mbal_time_temp ! Accumulation time in years

    mbal_time_temp = -1.0

    call GetSection(config,section,'GLINT climate')
    if (associated(section)) then
       call GetValue(section,'precip_mode',instance%whichprecip)
       call GetValue(section,'acab_mode',instance%whichacab)
       call GetValue(section,'ice_albedo',instance%ice_albedo)
       call GetValue(section,'lapse_rate',instance%lapse_rate)
       instance%data_lapse_rate=instance%lapse_rate
       call GetValue(section,'data_lapse_rate',instance%data_lapse_rate)
       call GetValue(section,'mbal_accum_time',mbal_time_temp)
       call GetValue(section,'ice_tstep_multiply',instance%ice_tstep_multiply)
       call GetValue(section,'mean_preserving',instance%use_mpint)
    end if

    if (mbal_time_temp>0.0) then
       instance%mbal_accum_time = mbal_time_temp * years2hours
    else
       instance%mbal_accum_time = -1
    end if

    call glint_nc_readparams(instance,config)

  end subroutine glint_i_readconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  
  subroutine glint_nc_readparams(instance,config)
    !*FD read netCDF I/O related configuration file
    !*FD based on glimmer_ncparams
    use glimmer_config
    use glimmer_ncparams, only: handle_output, handle_input, configstring
    implicit none
    type(glint_instance)         :: instance  !*FD GLINT instance
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output
    type(glimmer_nc_input), pointer :: input

    ! Initialise local pointers 
    output => null()
    input => null()

    ! setup outputs
    call GetSection(config,section,'GLINT output')
    do while(associated(section))
       output => handle_output(section,output,0.0,configstring)
       if (.not.associated(instance%out_first)) then
          instance%out_first => output
       end if
       call GetSection(section%next,section,'GLINT output')
    end do

    ! setup inputs
    call GetSection(config,section,'GLINT input')
    do while(associated(section))
       input => handle_input(section,input)
       if (.not.associated(instance%in_first)) then
          instance%in_first => input
       end if
       call GetSection(section%next,section,'GLINT input')
    end do
    
    output => null()
    input => null()

  end subroutine glint_nc_readparams

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_printconfig(instance)

    use glimmer_log
    use glint_constants, only: hours2years

    implicit none

    ! Argument

    type(glint_instance),intent(inout) :: instance    !*FD The instance to be printed

    ! Internal

    character(len=100) :: message

    call write_log('GLINT climate')
    call write_log('-------------')
    write(message,*) 'precip_mode ',instance%whichprecip
    call write_log(message)
    write(message,*) 'acab_mode   ',instance%whichacab
    call write_log(message)
    write(message,*) 'ice_albedo  ',instance%ice_albedo
    call write_log(message)
    write(message,*) 'lapse_rate  ',instance%lapse_rate
    call write_log(message)
    write(message,*) 'data_lapse_rate',instance%data_lapse_rate
    call write_log(message)
    if (instance%mbal_accum_time==-1) then
       call write_log('Mass-balance accumulation time == ice time-step')
    else
       write(message,*) 'Mass-balance accumulation time:',instance%mbal_accum_time * hours2years,' years'
       call write_log(message)
    end if
    write(message,*) 'ice_tstep_multiply',instance%ice_tstep_multiply
    call write_log(message)
    select case(instance%use_mpint)
    case(1)
       write(message,*) 'Using mean-preserving interpolation'
       call write_log(message)
    case(0)
       write(message,*) 'Using normal interpolation'
       call write_log(message)
    case default
       write(message,*) 'Unrecognised value of instance%use_mpint'
       call write_log(message,GM_FATAL)
    end select

    call write_log('')

  end subroutine glint_i_printconfig

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_i_upscaled_fields(instance,                    &
                                   orog,         albedo,        &
                                   ice_frac,     veg_frac,      &
                                   snowice_frac, snowveg_frac,  &
                                   snow_depth)

    !*FD Upscales and returns certain fields
    !*FD 
    !*FD \begin{itemize}
    !*FD \item \texttt{orog} --- the orographic elevation (m)
    !*FD \item \texttt{albedo} --- the albedo of ice/snow (this is only a notional value --- need to do
    !*FD some work here)
    !*FD \item \texttt{ice\_frac} --- The fraction covered by ice
    !*FD \item \texttt{veg\_frac} --- The fraction of exposed vegetation
    !*FD \item \texttt{snowice\_frac} --- The fraction of snow-covered ice
    !*FD \item \texttt{snowveg\_frac} --- The fraction of snow-covered vegetation
    !*FD \item \texttt{snow_depth} --- The mean snow-depth over those parts covered in snow (m w.e.)
    !*FD \end{itemize}

    use glimmer_paramets

    ! Arguments ----------------------------------------------------------------------------------------

    type(glint_instance),   intent(in)  :: instance      !*FD the model instance
    real(rk),dimension(:,:),intent(out) :: orog          !*FD the orographic elevation (m)
    real(rk),dimension(:,:),intent(out) :: albedo        !*FD the albedo of ice/snow
    real(rk),dimension(:,:),intent(out) :: ice_frac      !*FD The fraction covered by ice
    real(rk),dimension(:,:),intent(out) :: veg_frac      !*FD The fraction of exposed vegetation
    real(rk),dimension(:,:),intent(out) :: snowice_frac  !*FD The fraction of snow-covered ice
    real(rk),dimension(:,:),intent(out) :: snowveg_frac  !*FD The fraction of snow-covered vegetation
    real(rk),dimension(:,:),intent(out) :: snow_depth    !*FD The mean snow-depth over those 
    !*FD parts covered in snow (m w.e.)

    ! Internal variables -------------------------------------------------------------------------------

    real(rk),dimension(:,:),pointer :: temp => null()

    ! --------------------------------------------------------------------------------------------------
    ! Orography

    call mean_to_global(instance%ups_orog, &
         instance%model%geometry%usrf, &
         orog,    &
         instance%out_mask)
    orog=thk0*orog

    call coordsystem_allocate(instance%lgrid,temp)

    ! Ice-no-snow fraction
    where (instance%mbal_accum%snowd==0.0.and.instance%model%geometry%thck>0.0)
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
         temp, &
         ice_frac,    &
         instance%out_mask)

    ! Ice-with-snow fraction
    where (instance%mbal_accum%snowd>0.0.and.instance%model%geometry%thck>0.0)
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
         temp, &
         snowice_frac,    &
         instance%out_mask)

    ! Veg-with-snow fraction (if ice <10m thick)
    where (instance%mbal_accum%snowd>0.0.and.instance%model%geometry%thck<=(10.0/thk0))
       temp=1.0
    elsewhere
       temp=0.0
    endwhere
    call mean_to_global(instance%ups, &
         temp, &
         snowveg_frac,    &
         instance%out_mask)

    ! Remainder is veg only
    veg_frac=1.0-ice_frac-snowice_frac-snowveg_frac

    ! Snow depth

    call mean_to_global(instance%ups, &
         instance%mbal_accum%snowd, &
         snow_depth,    &
         instance%out_mask)

    ! Albedo

    where ((ice_frac+snowice_frac)>0.0)
       albedo=instance%ice_albedo
    elsewhere
       albedo=0.0
    endwhere

    deallocate(temp)
    temp => null()

  end subroutine get_i_upscaled_fields

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine get_i_upscaled_fields_gcm(instance,    nec,      &
                                       nxl,         nyl,      &
                                       nxg,         nyg,      &
                                       gfrac,       gtopo,    &
                                       grofi,       grofl,    &
                                       ghflx)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! The upscaled fields are passed to the GCM land surface model, which has the option
    !  of updating the fractional area and surface elevation of glaciated gridcells.

    use glimmer_paramets, only: thk0
    use glimmer_log

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance),     intent(in)  :: instance      ! the model instance
    integer,                  intent(in)  :: nec           ! number of elevation classes
    integer,                  intent(in)  :: nxl,nyl       ! local grid dimensions 
    integer,                  intent(in)  :: nxg,nyg       ! global grid dimensions 

    real(dp),dimension(nxg,nyg,nec),intent(out) :: gfrac   ! ice-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: grofi   ! ice runoff (calving) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: grofl   ! liquid runoff (basal melt) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: ghflx   ! heat flux (m)
 
    ! Internal variables ----------------------------------------------------------------------
 
    real(dp),dimension(nxl,nyl) :: local_field
    real(dp),dimension(nxl,nyl) :: local_topo   ! local surface elevation (m)
    real(dp),dimension(nxl,nyl) :: local_thck   ! local ice thickness (m)
    real(dp), parameter :: min_thck = 1.0_dp    ! min thickness (m) for setting gfrac = 1

    integer :: i, j            ! indices
 
    integer :: il, jl, ig, jg
    character(len=100) :: message

!lipscomb - TO DO - Read topomax from data file at initialization
    real(dp), dimension(0:nec) :: topomax   ! upper elevation limit of each class

    ! Given the value of nec, specify the upper and lower elevation boundaries of each class.
    ! Note: These must be consistent with the values in the GCM.  Better to pass as an argument.
    if (nec == 1) then
       topomax = (/ 0._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 3) then
       topomax = (/ 0._dp,  1000._dp,  2000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 5) then
       topomax = (/ 0._dp,   500._dp,  1000._dp,  1500._dp,  2000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 10) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   700._dp,  1000._dp,  1300._dp,  &
                            1600._dp,  2000._dp,  2500._dp,  3000._dp, 10000._dp /)
    else
       if (GLC_DEBUG) then 
          write(message,'(a6,i3)') 'nec =', nec
          call write_log(trim(message), GM_DIAGNOSTIC)
       endif
       call write_log('ERROR: Current supported values of nec (no. of elevation classes) are 1, 3, 5, or 10', &
                       GM_FATAL,__FILE__,__LINE__)
    endif

    local_topo(:,:) = thk0 * instance%model%geometry%usrf(:,:)
    local_thck(:,:) = thk0 * instance%model%geometry%thck(:,:)
        
    if (GLC_DEBUG) then
       ig = itest
       jg = jjtest
       il = itest_local
       jl = jtest_local
       write(stdout,*) 'In get_i_upscaled_fields_gcm'
       write(stdout,*) 'il, jl =', il, jl
       write(stdout,*) 'ig, jg =', ig, jg
       write(stdout,*) 'nxl, nyl =', nxl,nyl
       write(stdout,*) 'nxg, nyg =', nxg,nyg
       write(stdout,*) 'topo =', local_topo(il,jl) 
       write(stdout,*) 'thck =', local_thck(il,jl) 
       write(stdout,*) 'local out_mask =', instance%out_mask(il,jl)
    endif

    ! temporary field: = 1 where ice thickness exceeds threshold, else = 0

    do j = 1, nyl
    do i = 1, nxl
       if (local_thck(i,j) > min_thck) then
          local_field(i,j) = 1._dp
       else
          local_field(i,j) = 0._dp
       endif
    enddo
    enddo

    ! ice fraction

    call mean_to_global_mec(instance%ups,                       &
                            nxl,                nyl,            &
                            nxg,                nyg,            &
                            nec,                topomax,        &
                            local_field,        gfrac,          &
                            local_topo,         instance%out_mask)

    ! surface elevation

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_topo,          gtopo,     &
                            local_topo,          instance%out_mask)

!lipscomb - TO DO - Copy the appropriate fields into local_field array

    ! ice runoff

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         grofi,     &
                            local_topo,          instance%out_mask)

    ! liquid runoff

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         grofl,     &
                            local_topo,          instance%out_mask)

    ! heat flux

    local_field(:,:) = 0._dp

    call mean_to_global_mec(instance%ups,                   &
                            nxl,                 nyl,       &
                            nxg,                 nyg,       &
                            nec,                 topomax,   &
                            local_field,         ghflx,     &
                            local_topo,          instance%out_mask)
 
    if (GLC_DEBUG) then
!       write(stdout,*) ' '
!       write(stdout,*) 'global ifrac:'
!       do n = 1, nec
!          write(stdout,*) n, gfrac(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global gtopo:'
!       do n = 1, nec
!          write(stdout,*) n, gtopo(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofi:'
!       do n = 1, nec
!          write(stdout,*) n, grofi(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofl:'
!       do n = 1, nec
!          write(stdout,*) n, grofl(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global ghflx:'
!       do n = 1, nec
!          write(stdout,*) n, ghflx(ig, jg, n)
!       enddo
    endif

  end subroutine get_i_upscaled_fields_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical function glint_has_snow_model(instance)

    type(glint_instance),            intent(in)  :: instance

    glint_has_snow_model=mbal_has_snow_model(instance%mbal_accum%mbal)

  end function glint_has_snow_model

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine populate_output_flags(out_f,orog_out,albedo,ice_frac,veg_frac,snowice_frac, &
       snowveg_frac,snow_depth,water_in,water_out,total_water_in,total_water_out,ice_volume)

    type(output_flags),intent(inout) :: out_f
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


    out_f%orog         = present(orog_out)
    out_f%albedo       = present(albedo)
    out_f%ice_frac     = present(ice_frac)
    out_f%veg_frac     = present(veg_frac)
    out_f%snowice_frac = present(snowice_frac)
    out_f%snowveg_frac = present(snowveg_frac)
    out_f%snow_depth   = present(snow_depth)
    out_f%water_out    = present(water_out)
    out_f%water_in     = present(water_in)
    out_f%total_win    = present(total_water_in)
    out_f%total_wout   = present(total_water_out)
    out_f%ice_vol      = present(ice_volume)

  end subroutine populate_output_flags

end module glint_type
