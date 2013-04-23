! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_setup.f90 - part of the Glimmer-CISM ice model     + 
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

#include "glide_mask.inc"

module glide_setup

  ! general routines for initialisation, etc, called from top-level glimmer subroutines

  use glimmer_global, only: dp

  implicit none

  private
  public :: glide_readconfig, glide_printconfig, glide_scale_params, &
            glide_load_sigma, glide_read_sigma, glide_calc_sigma

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine glide_readconfig(model,config)

    ! read GLIDE configuration file
    ! Note: sigma coordinates are handled by a subsequent call to glide_read_sigma
 
    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section

    ! read grid size  parameters
    call GetSection(config,section,'grid')
    if (associated(section)) then
       call handle_grid(section, model)
    end if

    ! read time parameters
    call GetSection(config,section,'time')
    if (associated(section)) then
       call handle_time(section, model)
    end if

    ! read options parameters
    call GetSection(config,section,'options')
    if (associated(section)) then
       call handle_options(section, model)
    end if

    !read options for higher-order computation
    call GetSection(config,section,'ho_options')
    if (associated(section)) then
        call handle_ho_options(section, model)
    end if

     !read options for computation using an external dycore -- Doug Ranken 04/20/12
    call GetSection(config,section,'external_dycore_options')
    if (associated(section)) then
        call handle_dycore_options(section, model)
    end if

    ! read parameters
    call GetSection(config,section,'parameters')
    if (associated(section)) then
       call handle_parameters(section, model)
    end if

    ! read GTHF 
    call GetSection(config,section,'GTHF')
    if (associated(section)) then
!!       model%options%gthf = 1       
       model%options%gthf = GTHF_COMPUTE
       call handle_gthf(section, model)
    end if

    !WHL - disabled for now
    ! read till parameters
!!    call GetSection(config,section,'till_options')
!!    if (associated(section)) then
!!       call handle_till_options(section, model)
!!    end if

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file (specific to the glide & glide_litho sections - other sections,
    ! e.g. glint, isos, are handled separately by their setup routines).
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13

    call define_glide_restart_variables(model%options)

  end subroutine glide_readconfig

!-------------------------------------------------------------------------

  subroutine glide_printconfig(model)

    !*FD print model configuration to log
    use glimmer_log
    use glide_types
    implicit none
    type(glide_global_type)  :: model !*FD model instance

    call write_log_div
    call print_grid(model)
    call print_time(model)
    call print_options(model)
    call print_parameters(model)
    call print_gthf(model)
!!    call print_till_options(model)  ! disabled for now

  end subroutine glide_printconfig

!-------------------------------------------------------------------------
    
!TODO - Remove most of the scaling?
!       (Probably will need to keep scyr)

  subroutine glide_scale_params(model)
    !*FD scale parameters
    use glide_types
    use glimmer_physcon,  only: scyr

!SCALING - Can delete the following when scaling is removed
    use glimmer_physcon,  only: gn
    use glimmer_paramets, only: thk0,tim0,len0, vel0, vis0, acc0

    implicit none

    type(glide_global_type)  :: model !*FD model instance

!SCALING - Some of this code is still needed, even after thk0, len0, etc. are set to 1.0.
!          In particular, keep the conversions with scyr.

!TODO - Change ntem and nvel to dttem and dtvel?  Is nvel used?
    model%numerics%ntem = model%numerics%ntem * model%numerics%tinc   ! keep
    model%numerics%nvel = model%numerics%nvel * model%numerics%tinc   ! keep

    model%numerics%dt     = model%numerics%tinc * scyr / tim0   ! keep scyr?  (or let dt be in yr)
    model%numerics%dttem  = model%numerics%ntem * scyr / tim0   ! keep scyr?  (or let dt be in yr)
    model%numerics%thklim = model%numerics%thklim  / thk0       ! can remove scaling here 

    model%numerics%dew = model%numerics%dew / len0         ! remove
    model%numerics%dns = model%numerics%dns / len0         ! remove

    model%numerics%mlimit = model%numerics%mlimit / thk0   ! remove

    model%numerics%periodic_offset_ew = model%numerics%periodic_offset_ew / thk0
    model%numerics%periodic_offset_ns = model%numerics%periodic_offset_ns / thk0

    model%velowk%trc0   = vel0 * len0 / (thk0**2)          ! keep scyr
    model%velowk%btrac_const = model%paramets%btrac_const/model%velowk%trc0/scyr  ! remove? 
    model%velowk%btrac_max = model%paramets%btrac_max/model%velowk%trc0/scyr      ! remove?
    model%velowk%btrac_slope = model%paramets%btrac_slope*acc0/model%velowk%trc0  ! remove?

  end subroutine glide_scale_params

!-------------------------------------------------------------------------

  subroutine glide_read_sigma(model,config)

    ! read sigma levels from configuration file, if present
    ! called immediately after glide_readconfig

    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file
        
    ! local variables
    type(ConfigSection), pointer :: section

    ! read sigma levels
    ! NOTE: If a sigma section is present in the config file, then we will always
    !       use the sigma levels set there, regardless of the value of model%options%which_sigma

    call GetSection(config,section,'sigma')
    if (associated(section)) then
       call handle_sigma(section, model)
    end if

  end subroutine glide_read_sigma

!-------------------------------------------------------------------------

  subroutine glide_load_sigma(model,unit)

    ! Compute sigma coordinates or read them from a file
    ! Note: This subroutine is called from glide_initialise or glissade_initialise.
    !       If sigma levels are provided in the config file, then they are read
    !        in by glide_read_sigma, and model%options%which_sigma is set to
    !        SIGMA_CONFIG, in which case this subroutine does nothing.

    use glide_types
    use glimmer_log
    use glimmer_filenames
    use glimmer_global, only: dp
    use parallel

    implicit none

    ! Arguments
    type(glide_global_type),intent(inout) :: model !*FD Ice model to use
    integer,               intent(in)    :: unit   !*FD Logical file unit to use. 
                                                   !*FD (Must not already be in use)

    ! Internal variables

    integer :: up,upn
    logical :: there
    real(dp) :: level

    ! Beginning of code

    upn=model%general%upn

    select case(model%options%which_sigma)

!!    case(0)    
    case(SIGMA_COMPUTE_GLIDE)   !  compute standard Glide sigma levels

       do up = 1,upn
          level = real(up-1,kind=dp) / real(upn-1,kind=dp)
!!          model%numerics%sigma(up) = glide_find_level(level, model%options%which_sigma_builtin, up, upn)
          model%numerics%sigma(up) = glide_calc_sigma(level, 2.d0)
       end do

       call write_log('Computing Glide sigma levels')

!!    case(1)
    case(SIGMA_EXTERNAL)        ! read from external file

       if (main_task) inquire (exist=there, file=process_path(model%funits%sigfile))
       call broadcast(there)
       if (.not.there) then
          call write_log('Sigma levels file: '//trim(process_path(model%funits%sigfile))// &
               ' does not exist',GM_FATAL)
       end if
       call write_log('Reading sigma file: '//process_path(model%funits%sigfile))
       if (main_task) then
          open(unit,file=process_path(model%funits%sigfile))
          read(unit,'(f9.7)',err=10,end=10) (model%numerics%sigma(up), up=1,upn)
          close(unit)
       end if
       call broadcast(model%numerics%sigma)

!!    case(2)
    case(SIGMA_CONFIG)          ! read from config file

       ! sigma levels have already been read from glide_read_sigma

       call write_log('Getting sigma levels from main configuration file')

    case(SIGMA_COMPUTE_EVEN)

       do up = 1,upn
          model%numerics%sigma(up) = real(up-1,kind=dp) / real(upn-1,kind=dp)
       enddo

       call write_log('Computing evenly spaced sigma levels')

    case(SIGMA_COMPUTE_PATTYN)

       do up = 1,upn
          if (up == 1) then
             model%numerics%sigma(up) = 0.d0
          else if (up == upn) then
             model%numerics%sigma(up) = 1.d0
          else
             level = real(up-1,kind=dp) / real(upn-1,kind=dp)
             model%numerics%sigma(up) = glide_calc_sigma_pattyn(level)
          end if
       enddo

       call write_log('Computing Pattyn sigma levels')

    end select


    !WHL - adding logic so that Glam will always use evenly spaced levels,
    !      overriding other values of which_sigma (including sigma levels in config file)

    if (model%options%whichdycore == DYCORE_GLAM) then   ! evenly spaced levels are required

       do up = 1,upn
          model%numerics%sigma(up) = real(up-1,kind=dp) / real(upn-1,kind=dp)
       enddo

       call write_log('Using evenly spaced sigma levels for Glam as required')

    endif

    ! Compute stagsigma (= sigma values at layers midpoints)

    model%numerics%stagsigma(1:upn-1) =   &
            (model%numerics%sigma(1:upn-1) + model%numerics%sigma(2:upn)) / 2.0_dp

    ! Compute stagwbndsigma, adding the boundaries to stagsigma

    model%numerics%stagwbndsigma(1:upn-1) = model%numerics%stagsigma(1:upn-1)
    model%numerics%stagwbndsigma(0) = 0.d0
    model%numerics%stagwbndsigma(upn) = 1.d0        

    call print_sigma(model)

    return

10  call write_log('something wrong with sigma coord file',GM_FATAL)
    
  end subroutine glide_load_sigma

!--------------------------------------------------------------------------------

!WHL - I removed the call to this function.
!      I added the 'even' and 'Pattyn' options to the whichsigma options.

!!  function glide_find_level(level, scheme, up, upn)

  !Returns the sigma coordinate of one level using a specific builtin scheme

!!    use glide_types
!!    use glimmer_global, only: dp
!!    real(dp) :: level
!!    integer  :: scheme, up, upn
!!    real(dp) :: glide_find_level

!!    select case(scheme)
!!      case (SIGMA_BUILTIN_DEFAULT)
!!        glide_find_level = glide_calc_sigma(level,2D0)
!!      case (SIGMA_BUILTIN_EVEN)
!!        glide_find_level = level
!!      case (SIGMA_BUILTIN_PATTYN)
!!        if (up == 1) then
!!          glide_find_level = 0
!!        else if (up == upn) then
!!          glide_find_level = 1
!!        else
!!           glide_find_level = glide_calc_sigma_pattyn(level)
!!        end if
!!    end select
     
!!  end function glide_find_level

!--------------------------------------------------------------------------------

  function glide_calc_sigma(x,n)

     use glimmer_global, only: dp
     implicit none
     real(dp) :: glide_calc_sigma, x, n
      
     glide_calc_sigma = (1-(x+1)**(-n)) / (1-2**(-n))

  end function glide_calc_sigma

!--------------------------------------------------------------------------------

  function glide_calc_sigma_pattyn(x)

     ! Implements an alternate set of sigma levels that encourages better
     ! convergence for higher-order velocities

     use glimmer_global, only:dp
     implicit none
     real(dp) :: glide_calc_sigma_pattyn, x

     glide_calc_sigma_pattyn =   &
         (-2.5641025641d-4)*(41d0*x)**2+3.5256410256d-2*(41d0*x)-8.0047080075d-13

  end function glide_calc_sigma_pattyn

!--------------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! grid sizes

  subroutine handle_grid(section, model)
    use glimmer_config
    use glide_types
    use glimmer_filenames
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'ewn',model%general%ewn)
    call GetValue(section,'nsn',model%general%nsn)
    call GetValue(section,'upn',model%general%upn)
    call GetValue(section,'dew',model%numerics%dew)
    call GetValue(section,'dns',model%numerics%dns)
    call GetValue(section,'sigma_file',model%funits%sigfile)
      !WHL - removed this option
!!    call GetValue(section,'sigma_builtin',model%options%which_sigma_builtin)

    ! We set this flag to one to indicate we've got a sigfile name.
    ! A warning/error is generated if sigma levels are specified in some other way
    ! and mangle the name
    if (trim(model%funits%sigfile) /= '') then
       model%funits%sigfile = filenames_inputname(model%funits%sigfile)
       model%options%which_sigma = SIGMA_EXTERNAL
    end if

  end subroutine handle_grid

!--------------------------------------------------------------------------------

  subroutine print_grid(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=512) :: message


    call write_log('Grid specification')
    call write_log('------------------')
    write(message,*) 'ewn             : ',model%general%ewn
    call write_log(trim(message))
    write(message,*) 'nsn             : ',model%general%nsn
    call write_log(trim(message))
    write(message,*) 'upn             : ',model%general%upn
    call write_log(trim(message))
    write(message,*) 'EW grid spacing : ',model%numerics%dew
    call write_log(trim(message))
    write(message,*) 'NS grid spacing : ',model%numerics%dns
    call write_log(trim(message))
    write(message,*) 'sigma file      : ',trim(model%funits%sigfile)
    call write_log(trim(message))
    call write_log('')
  end subroutine print_grid

!--------------------------------------------------------------------------------

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

!TODO - To handle timesteps both greater and less than one year, we may want to
!        define ice_dt_option and ice_dt_count in place of the current dt.
!       For instance, ice_dt_option could be either 'nyears' or 'steps_per_year'.
!       For timesteps < 1 year, we would use ice_dt_option = 'steps_per_year'.
!       This would ensure that the ice sheet dynamic timestep divides evenly
!        into the mass balance timestep (= 1 year) when running with Glint.
    call GetValue(section,'tstart',model%numerics%tstart)
    call GetValue(section,'tend',model%numerics%tend)
    call GetValue(section,'dt',model%numerics%tinc)
    call GetValue(section,'subcyc',model%numerics%subcyc)
    call GetValue(section,'ntem',model%numerics%ntem)
    call GetValue(section,'nvel',model%numerics%nvel)
    call GetValue(section,'profile',model%numerics%profile_period)

    call GetValue(section,'dt_diag',model%numerics%dt_diag)
    call GetValue(section,'idiag',model%numerics%idiag_global)
    call GetValue(section,'jdiag',model%numerics%jdiag_global)

    !WHL - ndiag replaced by dt_diag, but retained (for now) for backward compatibility
    call GetValue(section,'ndiag',model%numerics%ndiag)

  end subroutine handle_time
  
!--------------------------------------------------------------------------------

  subroutine print_time(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Time steps')
    call write_log('----------')
    write(message,*) 'start time (yr)     : ',model%numerics%tstart
    call write_log(message)
    write(message,*) 'end time (yr)       : ',model%numerics%tend
    call write_log(message)
    write(message,*) 'time step (yr)      : ',model%numerics%tinc
    call write_log(message)
    write(message,*) 'thermal dt factor   : ',model%numerics%ntem
    call write_log(message)
    if ( (model%numerics%ntem < 1.0d0) .or. & 
       (floor(model%numerics%ntem) /= model%numerics%ntem) ) then
       call write_log('ntem is a multiplier on the basic time step.  It should be a positive integer.  Aborting.',GM_FATAL)
    endif
    write(message,*) 'velo dt factor      : ',model%numerics%nvel
    call write_log(message)
    write(message,*) 'profile frequency   : ',model%numerics%profile_period
    call write_log(message)

    if (model%numerics%dt_diag > 0.d0) then
       write(message,*) 'diagnostic time (yr): ',model%numerics%dt_diag
       call write_log(message)
       if (mod(model%numerics%dt_diag, model%numerics%tinc) > 1.e-11) then
          write(message,*) 'Warning: diagnostic interval does not divide evenly into ice timestep dt'
          call write_log(message)
       endif
    endif

    !WHL - ndiag replaced by dt_diag, but retained (for now) for backward compatibility
    if (model%numerics%ndiag > 0) then
       write(message,*) 'diag time (steps)   : ',model%numerics%ndiag
       call write_log(message)
    endif

    write(message,*) 'idiag_global        : ',model%numerics%idiag_global
    call write_log(message)
    write(message,*) 'jdiag_global        : ',model%numerics%jdiag_global
    call write_log(message)
    call write_log('')

  end subroutine print_time

!--------------------------------------------------------------------------------

  ! options
  subroutine handle_options(section, model)

    use glimmer_config
    use glide_types

    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'dycore',model%options%whichdycore)
    call GetValue(section,'evolution',model%options%whichevol)
    call GetValue(section,'temperature',model%options%whichtemp)
    call GetValue(section,'temp_init',model%options%temp_init)
    call GetValue(section,'flow_law',model%options%whichflwa)
    call GetValue(section,'slip_coeff',model%options%whichbtrc)
    call GetValue(section,'basal_water',model%options%whichbwat)
    call GetValue(section,'basal_mass_balance',model%options%basal_mbal)
    call GetValue(section,'gthf',model%options%gthf)
    call GetValue(section,'marine_margin',model%options%whichmarn)
    call GetValue(section,'vertical_integration',model%options%whichwvel)
    call GetValue(section,'topo_is_relaxed',model%options%whichrelaxed)
    call GetValue(section,'periodic_ew',model%options%periodic_ew)

    !WHL - added this one
    call GetValue(section,'sigma',model%options%which_sigma)

    !TODO - Not sure if this is still needed
    call GetValue(section,'ioparams',model%funits%ncfile)

    ! Both terms 'hotstart' and 'restart' are supported in the config file, 
    ! but if they are both supplied for some reason, then restart will be used.
    ! 'restart' is the preferred term moving forward.  
    ! 'hotstart' is retained for backward compatability.
    call GetValue(section,'hotstart',model%options%is_restart)
    call GetValue(section,'restart',model%options%is_restart)

    ! These are not currently supported
    !call GetValue(section, 'use_plume',model%options%use_plume)
    !call GetValue(section,'basal_proc',model%options%which_bproc)

!WHL - Removed these (no longer used)
!!    call GetValue(section,'periodic_ns',model%options%periodic_ns)  ! not used
!!    call GetValue(section,'diagnostic_run',model%options%diagnostic_run)

  end subroutine handle_options

!--------------------------------------------------------------------------------

  !Higher order options
  subroutine handle_ho_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type) :: model
    
    call GetValue(section, 'which_ho_efvs',      model%options%which_ho_efvs)
    call GetValue(section, 'which_disp',         model%options%which_disp)
    call GetValue(section, 'which_ho_babc',      model%options%which_ho_babc)
    call GetValue(section, 'which_ho_resid',     model%options%which_ho_resid)
    call GetValue(section, 'which_ho_nonlinear', model%options%which_ho_nonlinear)
    call GetValue(section, 'which_ho_sparse',    model%options%which_ho_sparse)

!WHL  - Added which_ho_approx option for glissade dycore. Commented out for now
!    call GetValue(section, 'which_ho_approx',    model%options%which_ho_approx)

    !WHL - Removed HO options that are not used
!!    call GetValue(section, 'which_bmelt',        model%options%which_bmelt)
!!    call GetValue(section, 'which_ho_source',    model%options%which_ho_source)
!!    call GetValue(section, 'include_thin_ice',   model%options%ho_include_thinice)
!!    call GetValue(section, 'guess_specified',    model%velocity%is_velocity_valid)
!!    call GetValue(section, 'which_ho_sparse_fallback', model%options%which_ho_sparse_fallback)

  end subroutine handle_ho_options

!--------------------------------------------------------------------------------

  ! Handles external dycore options -- Doug Ranken 03/26/12
  subroutine handle_dycore_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type) :: model
    
    call GetValue(section, 'external_dycore_type', model%options%external_dycore_type)
    call GetValue(section, 'dycore_input_file',  model%options%dycore_input_file)
    print *,"In handle_dycore_options, type_code, input file = ", &
             model%options%external_dycore_type,model%options%dycore_input_file 
  end subroutine handle_dycore_options

!--------------------------------------------------------------------------------

  subroutine print_options(model)

    !TODO - Make sure these are consistent with options in glide_types
    use glide_types
    use glimmer_log

    use parallel

    implicit none

    type(glide_global_type)  :: model
    character(len=500) :: message

    ! basic Glide/Glimmer options

    character(len=*), dimension(0:2), parameter :: dycore = (/ &
         'glide              ', &  ! Glimmer SIA
         'glam               ', &  ! Payne-Price finite difference
         'glissade           ' /)  ! prototype finite element

    character(len=*), dimension(-1:4), parameter :: evolution = (/ &
         'no thickness evolution                ', &
         'pseudo-diffusion                      ', &
         'ADI scheme                            ', &
         'iterated diffusion                    ', &
         'remap thickness                       ', &   
         '1st order upwind                      ' /)

    !TODO : Remove option 3 (after cleaning up config files)
    character(len=*), dimension(0:3), parameter :: temperature = (/ &
         'isothermal         ', &
         'full prognostic    ', &
         'constant in time   ', &
         'remapping advection' /)

    character(len=*), dimension(0:2), parameter :: temp_init = (/ &
         'set to 0 C             ', &
         'set to surface air temp', &
         'linear vertical profile' /)

    !WHL - Changed order
    character(len=*), dimension(0:2), parameter :: flow_law = (/ &
         'const 1e-16 Pa^-n a^-1      ', &
         'Paterson and Budd (T = -5 C)', &
         'Paterson and Budd           ' /)

    !WHL - Changed order
    !TODO - Rename to something like which_btrc?
    character(len=*), dimension(0:5), parameter :: slip_coeff = (/ &
         'no basal sliding       ', &
         'constant basal traction', &
         'constant where bwat > 0', &
         'constant where T = Tpmp', &
         'linear function of bmlt', &
         'tanh function of bwat  ' /)

    !WHL - Changed order
    character(len=*), dimension(0:3), parameter :: basal_water = (/ &
         'none                     ', &
         'local water balance      ', &
         'local + steady-state flux', &
         'Constant value (= 10 m)  ' /)
!!         'From basal proc model    '/) ! not supported

      !WHL - basal proc model is disabled for now.
!!    character(len=*), dimension(0:2), parameter :: which_bproc = (/ &
!!         'Basal proc mod disabled '  , &
!!         'Basal proc, high res.   '   , &
!!         'Basal proc, fast calc.  ' /)
    character(len=*), dimension(0:0), parameter :: which_bproc = (/ &
         'Basal process model disabled ' /)

    character(len=*), dimension(0:1), parameter :: b_mbal = (/ &
         'not in continuity eqn', &
         'in continuity eqn    ' /)

    ! NOTE: Set gthf = 1 in the config file to read the geothermal heat flux from an input file.
    !       Otherwise it will be overwritten, even if the 'bheatflx' field is present.

    character(len=*), dimension(0:2), parameter :: gthf = (/ &
         'uniform geothermal flux         ', &
         'read flux from file, if present ', &
         'compute flux from diffusion eqn ' /)

    !WHL - changed order
    character(len=*), dimension(0:5), parameter :: marine_margin = (/ &
         'do nothing at marine margin     ', &
         'remove all floating ice         ', &
         'remove fraction of floating ice ', &
         'relaxed bedrock threshold       ', &
         'present bedrock threshold       ', &
         'Huybrechts grounding line scheme' /)

    character(len=*), dimension(0:1), parameter :: vertical_integration = (/ &
         'standard     ', &
         'obey upper BC' /)

    ! higher-order options

    !WHL - Changed order
    character(len=*), dimension(0:2), parameter :: ho_whichefvs = (/ &
         'constant value                 ', & 
         'multiple of flow factor        ', &
         'nonlinear, from eff strain rate' /)

    character(len=*), dimension(0:1), parameter :: dispwhich = (/ &
         '0-order SIA                       ', &
         '1-st order model (Blatter-Pattyn) ' /)
!!         '1-st order depth-integrated (SSA) ' /)  ! not supported

    character(len=*), dimension(0:7), parameter :: ho_whichbabc = (/ &
         'constant B^2                           ', &
         'simple pattern of B^2                  ', &
         'till yield stress (Picard)             ', &
         'circular ice shelf                     ', &
         'no slip (using large B^2)              ', &
         'B^2 passed from CISM                   ', &
         'no slip (Dirichlet implementation)     ', &
         'till yield stress (Newton)             ' /)

    !WHL - Removed this option
    !      Glide always uses 0-order SIA, glam/glissade always uses 1st order,
    !      and SSA has never been supported
!!    character(len=*), dimension(0:2), parameter :: bmeltwhich = (/ &
!!         '0-order SIA                       ', &
!!         '1-st order model (Blatter-Pattyn) ', &
!!         '1-st order depth-integrated (SSA) ' /)

    !WHL - Removed this option
!!    character(len=*), dimension(0:2), parameter :: ho_whichsource = (/ &
!!         'vertically averaged     ', &
!!         'vertically explicit     ', &
!!         'shelf front disabled    '/)

    character(len=*), dimension(0:1), parameter :: which_ho_nonlinear = (/ &
         'use standard Picard iteration  ', &
         'use JFNK                       '/)

    character(len=*), dimension(0:3), parameter :: ho_whichresid = (/ &
         'max value               ', &
         'max value ignoring ubas ', &
         'mean value              ', &
         'L2 norm of Ax-b = resid ' /)

    character(len=*), dimension(0:4), parameter :: ho_whichsparse = (/ &
         'BiCG with LU preconditioner                ', &
         'GMRES with LU preconditioner               ', &
         'PCG with incomplete Cholesky preconditioner', &
         'Standalone PCG solver (structured)         ', &
         'Standalone Trilinos interface              '/)

!WHL - added glissade options for solving different Stokes approximations
!      commented out for now
!!    character(len=*), dimension(0:2), parameter :: ho_whichapprox = (/ &
!!         'SIA only (glissade dycore)         ', &
!!         'SSA only (glissade dycore)         ', &
!!         'Blatter-Pattyn HO (glissade dycore)' /)


    call write_log('GLIDE options')
    call write_log('-------------')

    write(message,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    call write_log(message)

    if (model%options%whichdycore < 0 .or. model%options%whichdycore >= size(dycore)) then
       call write_log('Error, dycore option out of range',GM_FATAL)
    end if
    write(message,*) 'Dycore                  : ',model%options%whichdycore,dycore(model%options%whichdycore)
    call write_log(message)

    !WHL - Option 3 (TEMP_REMAP_ADV) is now deprecated.  
    ! If this has been set, then change to option 1 (TEMP_PROGNOSTIC), which applies to any dycore.
    !TODO - Remove this option after cleaning up config files.
    if (model%options%whichtemp == TEMP_REMAP_ADV) model%options%whichtemp = TEMP_PROGNOSTIC

    if (model%options%whichtemp < 0 .or. model%options%whichtemp >= size(temperature)) then
       call write_log('Error, temperature option out of range',GM_FATAL)
    end if
    write(message,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    call write_log(message)

    ! Forbidden options to use with the Glide dycore
    if (model%options%whichdycore == DYCORE_GLIDE) then

       if (model%options%whichevol==EVOL_INC_REMAP) then
          call write_log('Error, incremental remapping evolution is not supported for the Glide dycore', GM_FATAL)
       endif

       if (tasks > 1) then
          call write_log('Error, Glide dycore not supported for runs with more than one processor', GM_FATAL)
       end if

       if (model%options%whichevol==EVOL_ADI) then
          call write_log('Warning, exact restarts are not currently possible with ADI evolution', GM_WARNING)
       endif

    else   ! glam/glissade dycore

       if (model%options%whichevol==EVOL_PSEUDO_DIFF .or.  &
           model%options%whichevol==EVOL_ADI         .or.  &
           model%options%whichevol==EVOL_DIFFUSION) then
          call write_log('Error, Glide thickness evolution options can be used only with Glide dycore', GM_FATAL)
       endif

    endif

    if (model%options%whichdycore /= DYCORE_GLISSADE) then 

       if (model%options%which_ho_sparse == HO_SPARSE_PCG_STRUC) then
          call write_log('Error, structured PCG solver requires glissade dycore')
       endif

    endif

    ! Forbidden options associated with Glam dycore
   
!WHL - commented out for now
!!    if (model%options%whichdycore == DYCORE_GLAM) then
!!       if (model%options%which_ho_approx == HO_APPROX_SIA .or.   &
!!           model%options%which_ho_approx == HO_APPROX_SSA) then 
!!          call write_log('Error, Glide dycore must use higher-order Blatter-Pattyn approximation', GM_FATAL)
!!       endif
!!    endif

    if (model%options%temp_init < 0 .or. model%options%temp_init >= size(temp_init)) then
       call write_log('Error, temp_init option out of range',GM_FATAL)
    end if
    ! Note: If reading temperature from an input or restart file, the temp_init option is overridden,
    !        in which case it could be confusing here to write the option to the log file.
    !       The method actually used is written to the log file by glide_init_temp. 

    if (model%options%whichflwa < 0 .or. model%options%whichflwa >= size(flow_law)) then
       call write_log('Error, flow_law out of range',GM_FATAL)
    end if
    write(message,*) 'flow law                : ',model%options%whichflwa,flow_law(model%options%whichflwa)
    call write_log(message)

    if (model%options%whichbwat < 0 .or. model%options%whichbwat >= size(basal_water)) then
       call write_log('Error, basal_water out of range',GM_FATAL)
    end if
    write(message,*) 'basal_water             : ',model%options%whichbwat,basal_water(model%options%whichbwat)
    call write_log(message)

    if (model%options%whichmarn < 0 .or. model%options%whichmarn >= size(marine_margin)) then
       call write_log('Error, marine_margin out of range',GM_FATAL)
    end if
    write(message,*) 'marine_margin           : ', model%options%whichmarn, marine_margin(model%options%whichmarn)
    call write_log(message)

    if (model%options%whichbtrc < 0 .or. model%options%whichbtrc >= size(slip_coeff)) then
       call write_log('Error, slip_coeff out of range',GM_FATAL)
    end if

    write(message,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
    call write_log(message)

    if (model%options%whichevol < -1 .or. model%options%whichevol >= (size(evolution)-1)) then
       call write_log('Error, evolution out of range',GM_FATAL)
    end if

    write(message,*) 'evolution               : ', model%options%whichevol, evolution(model%options%whichevol)
    call write_log(message)

    if (model%options%whichwvel < 0 .or. model%options%whichwvel >= size(vertical_integration)) then
       call write_log('Error, vertical_integration out of range',GM_FATAL)
    end if

    write(message,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    call write_log(message)

    if (model%options%basal_mbal < 0 .or. model%options%basal_mbal >= size(b_mbal)) then
       call write_log('Error, basal_mass_balance out of range',GM_FATAL)
    end if

    write(message,*) 'basal_mass_balance      : ',model%options%basal_mbal,b_mbal(model%options%basal_mbal)
    call write_log(message)

    if (model%options%gthf < 0 .or. model%options%gthf >= size(gthf)) then
       print*, 'gthf =', model%options%gthf
       call write_log('Error, geothermal flux option out of range',GM_FATAL)
    end if

    write(message,*) 'geothermal heat flux    : ',model%options%gthf,gthf(model%options%gthf)
    call write_log(message)

    if (model%options%whichrelaxed==1) then
       call write_log('First topo time slice is relaxed')
    end if

    if (model%options%periodic_ew) then
       if (model%options%whichevol == EVOL_ADI) then
          call write_log('Periodic boundary conditions not implemented in ADI scheme',GM_FATAL)
       end if
       call write_log('Periodic EW lateral boundary condition')
       call write_log('  Slightly cheated with how temperature is implemented.',GM_WARNING)
    end if

    if (model%options%is_restart == RESTART_TRUE) then
       call write_log('Restarting model from a previous run')
    end if

!!     This option is not currently supported
!!    if (model%options%which_bproc < 0 .or. model%options%which_bproc >= size(which_bproc)) then
!!       call write_log('Error, basal_proc out of range',GM_FATAL)
!!    end if
!!    write(message,*) 'basal_proc              : ',model%options%which_bproc,which_bproc(model%options%which_bproc)
!!    call write_log(message)

    !HO options

    if (model%options%whichdycore /= DYCORE_GLIDE) then   ! glam/glissade higher-order

       call write_log(' ')
       call write_log('Higher-order options:')
       call write_log('----------')

       write(message,*) 'ho_whichefvs            : ',model%options%which_ho_efvs,  &
                         ho_whichefvs(model%options%which_ho_efvs)
       call write_log(message)
       if (model%options%which_ho_efvs < 0 .or. model%options%which_ho_efvs >= size(ho_whichefvs)) then
          call write_log('Error, HO effective viscosity input out of range', GM_FATAL)
       end if

       write(message,*) 'dispwhich               : ',model%options%which_disp,  &
                         dispwhich(model%options%which_disp)
       call write_log(message)
       if (model%options%which_disp < 0 .or. model%options%which_disp >= size(dispwhich)) then
          call write_log('Error, which dissipation input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichbabc            : ',model%options%which_ho_babc,  &
                         ho_whichbabc(model%options%which_ho_babc)
       call write_log(message)
       if (model%options%which_ho_babc < 0 .or. model%options%which_ho_babc >= size(ho_whichbabc)) then
          call write_log('Error, HO basal BC input out of range', GM_FATAL)
       end if

       write(message,*) 'which_ho_nonlinear      : ',model%options%which_ho_nonlinear,  &
                         which_ho_nonlinear(model%options%which_ho_nonlinear)
       call write_log(message)
       if (model%options%which_ho_nonlinear < 0 .or. model%options%which_ho_nonlinear >= size(which_ho_nonlinear)) then
          call write_log('Error, HO nonlinear solution input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichresid           : ',model%options%which_ho_resid,  &
                         ho_whichresid(model%options%which_ho_resid)
       call write_log(message)
       if (model%options%which_ho_resid < 0 .or. model%options%which_ho_resid >= size(ho_whichresid)) then
          call write_log('Error, HO residual input out of range', GM_FATAL)
       end if

       write(message,*) 'ho_whichsparse          : ',model%options%which_ho_sparse,  &
                         ho_whichsparse(model%options%which_ho_sparse)
       call write_log(message)
       if (model%options%which_ho_sparse < 0 .or. model%options%which_ho_sparse >= size(ho_whichsparse)) then
          call write_log('Error, HO sparse solver input out of range', GM_FATAL)
       end if

!WHL - commented out for now
!!       write(message,*) 'ho_whichapprox          : ',model%options%which_ho_approx,  &
!!                         ho_whichapprox(model%options%which_ho_approx)
!!       call write_log(message)
!!       if (model%options%which_ho_approx < 0 .or. model%options%which_ho_approx >= size(ho_whichapprox)) then
!!          call write_log('Error, Stokes approximation out of range', GM_FATAL)
!!       end if

         !WHL - Removed this option
!!       write(message,*) 'bmeltwhich              : ',model%options%which_bmelt,  &
!!                         bmeltwhich(model%options%which_bmelt)
!!       call write_log(message)
!!       if (model%options%which_bmelt < 0 .or. model%options%which_bmelt >= size(bmeltwhich)) then
!!          call write_log('Error, which bmelt input out of range', GM_FATAL)
!!       end if

         !WHL - Removed this option
!!       if (model%options%which_ho_source < 0 .or. model%options%which_ho_source >= size(ho_whichsource)) then
!!          call write_log('Error, which_ho_source out of range', GM_FATAL)
!!       end if
!!       write(message,*) 'ice_shelf_source_term   :',model%options%which_ho_source, &
!!                         ho_whichsource(model%options%which_ho_source)
!!       call write_log(message)

    endif   ! whichdycore

  end subroutine print_options

!--------------------------------------------------------------------------------

  ! parameters
  subroutine handle_parameters(section, model)

    use glimmer_config
    use glide_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model
    real, pointer, dimension(:) :: tempvar => NULL()
    integer :: loglevel

    loglevel = GM_levels-GM_ERROR

    call GetValue(section,'log_level',loglevel)
    call glimmer_set_msg_level(loglevel)
    call GetValue(section,'ice_limit',        model%numerics%thklim)
    call GetValue(section,'marine_limit',     model%numerics%mlimit)
    call GetValue(section,'calving_fraction', model%numerics%calving_fraction)
    call GetValue(section,'geothermal',       model%paramets%geot)
    call GetValue(section,'flow_factor',      model%paramets%flow_factor)
    call GetValue(section,'default_flwa',     model%paramets%default_flwa)
    call GetValue(section,'hydro_time',       model%paramets%hydtim)

    ! NOTE: bpar is used only for BTRC_TANH_BWAT
    !       btrac_max and btrac_slope are used (with btrac_const) for BTRC_LINEAR_BMLT
    !       btrac_const is used for several options

    call GetValue(section,'basal_tract_const', model%paramets%btrac_const)
    call GetValue(section,'basal_tract_max',   model%paramets%btrac_max)
    call GetValue(section,'basal_tract_slope', model%paramets%btrac_slope)

    !WHL - Changed this so that bpar can be read correctly from config file.
    !      Parameter is now called 'basal_tract_tanh' instead of 'basal_tract'.
    call GetValue(section,'basal_tract_tanh',  tempvar, 5)
    if (associated(tempvar)) then
!!       model%paramets%btrac_const = tempvar(1)  ! old code
       model%paramets%bpar(:) = tempvar(:)
       deallocate(tempvar)
    end if

!!    call GetValue(section,'sliding_constant',  model%climate%slidconst)  ! not currently used

    !WHL - added for ismip-hom
    call GetValue(section,'periodic_offset_ew',model%numerics%periodic_offset_ew)
    call GetValue(section,'periodic_offset_ns',model%numerics%periodic_offset_ns)

  end subroutine handle_parameters

!--------------------------------------------------------------------------------

  subroutine print_parameters(model)

    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log(' ')
    call write_log('Parameters')
    call write_log('----------')

    write(message,*) 'ice limit (m)                 : ',model%numerics%thklim
    call write_log(message)

    write(message,*) 'marine depth limit (m)        : ',model%numerics%mlimit
    call write_log(message)

    if (model%options%whichmarn == MARINE_FLOAT_FRACTION) then
       write(message,*) 'ice fraction lost due to calving : ', model%numerics%calving_fraction
       call write_log(message)
    end if

    write(message,*) 'geothermal flux  (W/m2)       : ', model%paramets%geot
    call write_log(message)

    write(message,*) 'flow enhancement factor       : ', model%paramets%flow_factor
    call write_log(message)

    write(message,*) 'basal hydro time constant (yr): ', model%paramets%hydtim
    call write_log(message)

    if (model%options%whichbtrc == BTRC_CONSTANT      .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_BWAT .or.  &
        model%options%whichbtrc == BTRC_LINEAR_BMLT   .or.  &
        model%options%whichbtrc == BTRC_CONSTANT_TPMP) then
       write(message,*) 'basal traction param (m/yr/Pa): ', model%paramets%btrac_const
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_TANH_BWAT) then
       write(message,*) 'basal traction tanh factors: ',model%paramets%bpar(1)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(2)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(3)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(4)
       call write_log(message)
       write(message,*) '                             ',model%paramets%bpar(5)
       call write_log(message)
    end if

    if (model%options%whichbtrc == BTRC_LINEAR_BMLT) then
       write(message,*) 'basal traction max            : ',model%paramets%btrac_max
       call write_log(message)
       write(message,*) 'basal traction slope          : ',model%paramets%btrac_slope
       call write_log(message)
    end if

    !TODO - Should this be in a different subroutine?
    if (model%numerics%idiag_global < 1 .or. model%numerics%idiag_global > model%general%ewn     &
                                        .or.                                                     &
        model%numerics%jdiag_global < 1 .or. model%numerics%jdiag_global > model%general%nsn) then
        call write_log('Error, global diagnostic point (idiag, jdiag) is out of bounds', GM_FATAL)
    endif

!WHL - added for ismip-hom
    if (model%numerics%periodic_offset_ew /= 0.d0) then
       write(message,*) 'periodic offset_ew (m)  :  ',model%numerics%periodic_offset_ew
       call write_log(message)
    endif

    if (model%numerics%periodic_offset_ns /= 0.d0) then
       write(message,*) 'periodic offset_ns (m)  :  ',model%numerics%periodic_offset_ns
       call write_log(message)
    endif
   
    call write_log('')
  end subroutine print_parameters

!--------------------------------------------------------------------------------

  ! Sigma levels
  subroutine handle_sigma(section, model)

    use glimmer_config
    use glide_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    if (model%options%which_sigma==SIGMA_EXTERNAL) then
       call write_log('Sigma levels specified twice - use only'// &
            ' config file or separate file, not both',GM_FATAL)
    else
       model%options%which_sigma = SIGMA_CONFIG   ! override any previously set value
       call GetValue(section,'sigma_levels',model%numerics%sigma,model%general%upn)
    end if

  end subroutine handle_sigma

!--------------------------------------------------------------------------------

  subroutine print_sigma(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message,temp
    integer :: i

    call write_log('Sigma levels:')
    call write_log('------------------')
    message=''
    do i=1,model%general%upn
       write(temp,'(f6.3)') model%numerics%sigma(i)
       message=trim(message)//trim(temp)
    enddo
    call write_log(trim(message))
    call write_log('')
    
  end subroutine print_sigma

!--------------------------------------------------------------------------------

  ! geothermal heat flux calculations
  subroutine handle_gthf(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'num_dim',model%lithot%num_dim)
    call GetValue(section,'nlayer',model%lithot%nlayer)
    call GetValue(section,'surft',model%lithot%surft)
    call GetValue(section,'rock_base',model%lithot%rock_base)
    call GetValue(section,'numt',model%lithot%numt)
    call GetValue(section,'rho',model%lithot%rho_r)
    call GetValue(section,'shc',model%lithot%shc_r)
    call GetValue(section,'con',model%lithot%con_r)
  end subroutine handle_gthf

!--------------------------------------------------------------------------------

  subroutine print_gthf(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message
    
!!    if (model%options%gthf > 0) then
    if (model%options%gthf == GTHF_COMPUTE) then
       call write_log('GTHF configuration')
       call write_log('------------------')
       if (model%lithot%num_dim==1) then
          call write_log('solve 1D diffusion equation')
       else if (model%lithot%num_dim==3) then          
          call write_log('solve 3D diffusion equation')
       else
          call write_log('Wrong number of dimensions.',GM_FATAL,__FILE__,__LINE__)
       end if
       write(message,*) 'number of layers                     : ',model%lithot%nlayer
       call write_log(message)
       write(message,*) 'initial surface temperature          : ',model%lithot%surft
       call write_log(message)
       write(message,*) 'rock base                            : ',model%lithot%rock_base
       call write_log(message)
       write(message,*) 'density of rock layer                : ',model%lithot%rho_r
       call write_log(message)
       write(message,*) 'specific heat capacity of rock layer : ',model%lithot%shc_r
       call write_log(message)
       write(message,*) 'thermal conductivity of rock layer   : ',model%lithot%con_r
       call write_log(message)
       write(message,*) 'number of time steps for spin-up     : ',model%lithot%numt
       call write_log(message)
       call write_log('')
    end if
  end subroutine print_gthf

!--------------------------------------------------------------------------------

!WHL - These options are disabled for now.
!Till options
!!  subroutine handle_till_options(section,model)
!!    use glimmer_config
!!    use glide_types
!!    implicit none
!!    type(ConfigSection), pointer :: section
!!    type(glide_global_type) :: model

!!    if (model%options%which_bproc==1) then
!!        call GetValue(section, 'fric',  model%basalproc%fric)
!!        call GetValue(section, 'etillo',  model%basalproc%etillo)
!!        call GetValue(section, 'No',  model%basalproc%No)
!!        call GetValue(section, 'Comp',  model%basalproc%Comp)
!!        call GetValue(section, 'Cv',  model%basalproc%Cv)
!!        call GetValue(section, 'Kh',  model%basalproc%Kh)
!!    else if (model%options%which_bproc==2) then
!!        call GetValue(section, 'aconst',  model%basalproc%aconst)
!!        call GetValue(section, 'bconst',  model%basalproc%bconst)
!!    end if
!!    if (model%options%which_bproc > 0) then
!!        call GetValue(section, 'Zs',  model%basalproc%Zs)
!!        call GetValue(section, 'tnodes',  model%basalproc%tnodes)
!!        call GetValue(section, 'till_hot', model%basalproc%till_hot)
!!    end if  
!!  end subroutine handle_till_options    

!WHL - These options are disabled for now.
!!  subroutine print_till_options(model)
!!    use glide_types
!!    use glimmer_log
!!    implicit none
!!    type(glide_global_type)  :: model
!!    character(len=100) :: message

!!    if (model%options%which_bproc > 0) then 
!!        call write_log('Till options')
!!        call write_log('----------')
!!        if (model%options%which_bproc==1) then
!!            write(message,*) 'Internal friction           : ',model%basalproc%fric
!!            call write_log(message)
!!            write(message,*) 'Reference void ratio        : ',model%basalproc%etillo
!!            call write_log(message)
!!            write(message,*) 'Reference effective Stress  : ',model%basalproc%No
!!            call write_log(message)
!!            write(message,*) 'Compressibility             : ',model%basalproc%Comp
!!            call write_log(message)
!!            write(message,*) 'Diffusivity                 : ',model%basalproc%Cv
!!            call write_log(message)
!!            write(message,*) 'Hyd. conductivity           : ',model%basalproc%Kh
!!            call write_log(message)
!!        end if
!!        if (model%options%which_bproc==2) then
!!            write(message,*) 'aconst  : ',model%basalproc%aconst
!!            call write_log(message)
!!            write(message,*) 'bconst  : ',model%basalproc%aconst
!!            call write_log(message)
!!        end if
!!        write(message,*) 'Solid till thickness : ',model%basalproc%Zs
!!        call write_log(message)
!!        write(message,*) 'Till nodes number : ',model%basalproc%tnodes
!!        call write_log(message)
!!        write(message,*) 'till_hot  :',model%basalproc%till_hot
!!        call write_log(message)
!!    end if
!!  end subroutine print_till_options

!--------------------------------------------------------------------------------

  subroutine define_glide_restart_variables(options)
    !*FD This subroutine analyzes the glide/glissade options input by the user in the config file
    !*FD and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.

    use glide_types
    use glide_io, only: glide_add_to_restart_variable_list

!    use glide_lithot_io, only: glide_lithot_add_to_restart_variable_list
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_options), intent (in) :: options  !*FD Derived type holding all model options

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    !This was the restart list as of 1/11/13 using the old hot=1 systme in glide_vars.def:
    !restart_variable_list=' lat  relx  tauf  thk  thkmask  topg  bheatflx  bmlt  bwat  uvel  vvel  wgrd  flwa  temp  litho_temp  age '

    ! Start with a few variables that we always want - prognostic variables and b.c.
    ! topg - needed to reconstruct all other geometry fields
    ! thk - prognostic variable
    ! temp - prognostic variable
    ! Note: the conversion from temp/flwa to tempstag/flwastag (if necessary) happens in glide_io.F90
    ! bheatflx, artm, acab - boundary conditions.  Of course if these fields are 0 they don't need 
    !        to be in the restart file, but without adding a check for that we cannot assume any of them are.
    !        There are some options where artm would not be needed.  Logic could be added to make that distinction.
    call glide_add_to_restart_variable_list('topg thk temp bheatflx artm acab')

    ! add dycore specific restart variables
    select case (options%whichdycore)

      case (DYCORE_GLIDE)
        ! thkmask - TODO is this needed?
        ! wgrd & wvel - temp driver calculates weff = f(wgrd, wvel) so both are needed by temp code.
        !               It looks possible to calculate wvel on a restart from wgrd because wvel does not 
        !               appear to require a time derivative (see subroutine wvelintg).  
        !               wgrd does require time derivatives and therefore should be
        !               calculated at the end of each time step and stored as a restart variable
        !               so that the time derivatives do not need to be restart variables.
        !               For now I am calculating wvel at the same time (end of glide time step) 
        !               and then saving both as restart variables.  This has the advantage of
        !               them being on consistent time levels in the output file.  
        !               (If we waited to calculate wvel in the temp driver, we would not need to
        !                add it as a restart variable, been then in the output wgrd and wvel would
        !                be based on different time levels.)
        ! flwa - in principal this could be reconstructed from temp.  However in the current 
        !        implementation of glide the flwa calculation occurs after temp evolution but 
        !        before thk evolution.  This means flwa is calculated from the current temp and 
        !        the old thk.  The old thk is not available on a restart (just the current thk).
        !        (thk is needed to calculate flwa for 1) a mask for where ice is, 2) correction for pmp.)
        call glide_add_to_restart_variable_list('thkmask wgrd wvel flwa')
    
        ! slip option for SIA
        select case (options%whichbtrc)
          case (0)
            ! no restart variable needed when no-slip is chosen
          case default
            call glide_add_to_restart_variable_list('btrc')
        end select

        if (options%whichevol == EVOL_PSEUDO_DIFF) then
          ! uvel,vvel are needed for linear diffusion evolution.  Their values
          ! cannot be reconstructed from the current thk or diffu because their
          ! values were originaly calculated from the old thk, which is not available
          ! on a restart.
          call glide_add_to_restart_variable_list('uvel vvel')
        endif

      case (DYCORE_GLAM, DYCORE_GLISSADE)
        ! uvel,vvel - these are needed for an exact restart because we can only 
        !             recalculate them to within the picard/jfnk convergence tolerance.
        ! beta - b.c. needed for runs with sliding - could add logic to only include in that case
        ! flwa is not needed for glissade.
        ! TODO not sure if thkmask is needed for HO
        ! TODO not sure if wgrd is needed for HO
        call glide_add_to_restart_variable_list('uvel vvel beta thkmask wgrd')

    end select

    ! ==== Other non-dycore specific options ====

    ! basal water option
    select case (options%whichbwat)
      case (BWATER_NONE, BWATER_CONST)
        ! no restart needed
      case default
        ! restart needs to know bwat value
        call glide_add_to_restart_variable_list('bwat')
    end select

    !WHL - This should be handled by glide_add_to_restart_variable_list
    ! if the GTHF calculation is enabled, litho_temp needs to be a restart variable
!!    if (options%gthf == GTHF_COMPUTE) then
!!        call glide_lithot_add_to_restart_variable_list('litho_temp')
!!    endif

    ! basal processes module - requires tauf for a restart
!!    if (options%which_bproc /= BAS_PROC_DISABLED ) then
!!        call glide_add_to_restart_variable_list('tauf')
!!    endif

    ! TODO bmlt was set as a restart variable, but I'm not sure when or if it is needed.

    ! TODO age should be a restart variable if it is an input variable.  
    ! Same goes for b.c. (bheatflxm, artm, acab) and any other tracers that get introduced.
    ! These could be included all the time (as I have down above for b.c.), or 
    ! we could add logic to only include them when they were in the input file.
    ! To do this, this subroutine would have to be moved to after where input files are read,
    ! glide_io_readall(), but before the output files are created, glide_io_createall()

    ! TODO lat is only needed for some climate drivers.  It is not needed for simple_glide.
    ! Need to add logic that will add it only when those drivers are used.

  end subroutine define_glide_restart_variables

!--------------------------------------------------------------------------------

end module glide_setup

!--------------------------------------------------------------------------------
