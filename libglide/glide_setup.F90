! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_setup.f90 - part of the Glimmer-CISM ice model     + 
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

module glide_setup

  !*FD Contains general routines for initialisation, etc, called
  !*FD from the top-level glimmer subroutines.

  private
  public :: glide_readconfig, glide_printconfig, glide_scale_params, &
       glide_calclsrf, glide_marinlim, glide_load_sigma, glide_maskthck, &
       glide_read_sigma, glide_calc_sigma

contains
  
  subroutine glide_readconfig(model,config)
    !*FD read GLIDE configuration file
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
    ! read parameters
    call GetSection(config,section,'parameters')
    if (associated(section)) then
       call handle_parameters(section, model)
    end if
    ! read GTHF 
    call GetSection(config,section,'GTHF')
    if (associated(section)) then
       model%options%gthf = 1
       call handle_gthf(section, model)
    end if
  end subroutine glide_readconfig

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
  end subroutine glide_printconfig
    
  subroutine glide_scale_params(model)
    !*FD scale parameters
    use glide_types
    use glimmer_physcon,  only: scyr, gn
    use glimmer_paramets, only: thk0,tim0,len0, tau0, vel0, vis0, acc0
    implicit none
    type(glide_global_type)  :: model !*FD model instance

    tau0 = (vel0/(vis0*len0))**(1.0/gn)

    model%numerics%ntem = model%numerics%ntem * model%numerics%tinc
    model%numerics%nvel = model%numerics%nvel * model%numerics%tinc

    model%numerics%dt     = model%numerics%tinc * scyr / tim0
    model%numerics%dttem  = model%numerics%ntem * scyr / tim0 
    model%numerics%thklim = model%numerics%thklim  / thk0

    model%numerics%dew = model%numerics%dew / len0
    model%numerics%dns = model%numerics%dns / len0

    model%numerics%mlimit = model%numerics%mlimit / thk0

    model%velowk%trc0   = vel0 * len0 / (thk0**2)
    model%velowk%btrac_const = model%paramets%btrac_const/model%velowk%trc0/scyr
    model%velowk%btrac_max = model%paramets%btrac_max/model%velowk%trc0/scyr
    model%velowk%btrac_slope = model%paramets%btrac_slope*acc0/model%velowk%trc0

  end subroutine glide_scale_params

  subroutine glide_read_sigma(model,config)
    !*FD read sigma levels from configuration file
    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file
        
    ! local variables
    type(ConfigSection), pointer :: section

    ! read grid size  parameters
    call GetSection(config,section,'sigma')
    if (associated(section)) then
       call handle_sigma(section, model)
    end if

  end subroutine glide_read_sigma

!-------------------------------------------------------------------------

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glide_maskthck(whichthck,crita,critb,dom,pointno,totpts,empty)
    
    !*FD Calculates the contents of the mask array.

    use glimmer_global, only : dp, sp 

    implicit none

    !-------------------------------------------------------------------------
    ! Subroutine arguments
    !-------------------------------------------------------------------------

    integer,                intent(in)  :: whichthck  !*FD Option determining
                                                      !*FD which method to use.
    real(dp),dimension(:,:),intent(in)  :: crita      !*FD Ice thickness
    real(sp),dimension(:,:),intent(in)  :: critb      !*FD Mass balance
    integer, dimension(:),  intent(out) :: dom        
    integer, dimension(:,:),intent(out) :: pointno    !*FD Output mask
    integer,                intent(out) :: totpts     !*FD Total number of points
    logical,                intent(out) :: empty      !*FD Set if no mask points set.

    !-------------------------------------------------------------------------
    ! Internal variables
    !-------------------------------------------------------------------------

    integer,dimension(size(crita,2),2) :: band
    logical,dimension(size(crita,2))   :: full
    integer :: covtot 
    integer :: ew,ns,ewn,nsn

    !-------------------------------------------------------------------------

    ewn=size(crita,1) ; nsn=size(crita,2)

    pointno = 0
    covtot  = 0 

    !-------------------------------------------------------------------------

    do ns = 1,nsn

      full(ns) = .false.

      do ew = 1,ewn
        if ( thckcrit(crita(max(1,ew-1):min(ewn,ew+1),max(1,ns-1):min(nsn,ns+1)),critb(ew,ns)) ) then

          covtot = covtot + 1
          pointno(ew,ns) = covtot 

          if ( .not. full(ns) ) then
            band(ns,1) = ew
            full(ns)   = .true.
          else
            band(ns,2) = ew
          end if
               
        end if
      end do
    end do
  
    totpts = covtot
                                             
    dom(1:2) = (/ewn,1/); empty = .true.

    do ns = 1,nsn
           
      if (full(ns)) then

        if (empty) then
          empty  = .false.
          dom(3) = ns
        end if
        dom(4) = ns
        dom(1) = min0(dom(1),band(ns,1))
        dom(2) = max0(dom(2),band(ns,2))
      end if
    end do

  contains

    logical function thckcrit(ca,cb)

      implicit none

      real(dp),dimension(:,:),intent(in) :: ca 
      real(sp),               intent(in) :: cb

      select case (whichthck)
      case(5)

        ! whichthck=5 is not a 'known case'

        if ( ca(2,2) > 0.0d0 .or. cb > 0.0) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      case default

        ! If the thickness in the region under consideration
        ! or the mass balance is positive, thckcrit is .true.

        if ( any((ca(:,:) > 0.0d0)) .or. cb > 0.0 ) then
          thckcrit = .true.
        else
          thckcrit = .false.
        end if

      end select

    end function thckcrit

  end subroutine glide_maskthck

  subroutine glide_calclsrf(thck,topg,eus,lsrf)

    !*FD Calculates the elevation of the lower surface of the ice, 
    !*FD by considering whether it is floating or not.

    use glimmer_global, only : dp
    use glimmer_physcon, only : rhoi, rhoo

    implicit none

    real(dp), intent(in),  dimension(:,:) :: thck !*FD Ice thickness
    real(dp), intent(in),  dimension(:,:) :: topg !*FD Bedrock topography elevation
    real, intent(in)                      :: eus  !*FD global sea level
    real(dp), intent(out), dimension(:,:) :: lsrf !*FD Lower ice surface elevation

    real(dp), parameter :: con = - rhoi / rhoo

    where (topg-eus < con * thck)
      lsrf = con * thck
    elsewhere
      lsrf = topg
    end where
  end subroutine glide_calclsrf

!-------------------------------------------------------------------------

  subroutine glide_marinlim(which,thck,relx,topg,mask,mlimit,calving_fraction,eus,ablation_field)

    !*FD Removes non-grounded ice, according to one of two altenative
    !*FD criteria, and sets upper surface of non-ice-covered points 
    !*FD equal to the topographic height, or sea-level, whichever is higher.

    use glimmer_global, only : dp, sp
    use glimmer_physcon, only : rhoi, rhoo
    use glide_mask
    implicit none

    !---------------------------------------------------------------------
    ! Subroutine arguments
    !---------------------------------------------------------------------

    integer,                intent(in)    :: which   !*FD Option to choose ice-removal method
                                                     !*FD \begin{description}
                                                     !*FD \item[0] Set thickness to zero if 
                                                     !*FD relaxed bedrock is below a given level.
                                                     !*FD \item[1] Set thickness to zero if
                                                     !*FD ice is floating.
                                                     !*FD \end{description}
    real(dp),dimension(:,:),intent(inout) :: thck    !*FD Ice thickness (scaled)
    real(dp),dimension(:,:),intent(in)    :: relx    !*FD Relaxed topography (scaled)
    real(dp),dimension(:,:),intent(in)    :: topg    !*FD Present bedrock topography (scaled)
    integer, dimension(:,:),pointer       :: mask    !*FD grid type mask
    real(dp)                              :: mlimit  !*FD Lower limit on topography elevation for
                                                     !*FD ice to be present (scaled). Used with 
                                                     !*FD $\mathtt{which}=0$.
    real(dp), intent(in) :: calving_fraction         !*FD fraction of ice lost when calving Used with 
                                                     !*FD $\mathtt{which}=3$.
    real, intent(inout) :: eus                       !*FD eustatic sea level
    real(sp),dimension(:,:),intent(inout) :: ablation_field 

    integer ew,ns
    real(dp), parameter :: con = - rhoi / rhoo
    !---------------------------------------------------------------------

    ablation_field=0.0

    select case (which)
        
    case(1) ! Set thickness to zero if ice is floating
      where (is_float(mask))
        ablation_field=thck
        thck = 0.0d0
      end where

    case(2) ! Set thickness to zero if relaxed bedrock is below a 
       ! given level
       where (relx <= mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where

    case(3) ! remove fraction of ice when floating
       do ns = 2,size(thck,2)-1
          do ew = 2,size(thck,1)-1
             if (is_calving(mask(ew,ns))) then
                ablation_field(ew,ns)=(1.0-calving_fraction)*thck(ew,ns)
                thck(ew,ns) =  calving_fraction*thck(ew,ns)
                !mask(ew,ns) = glide_mask_ocean
             end if
          end do
       end do

    case(4) ! Set thickness to zero at marine edge if present
            ! bedrock is below a given level
       where (is_marine_ice_edge(mask).and.topg<mlimit+eus)
          ablation_field=thck
          thck = 0.0d0
       end where

    end select

  end subroutine glide_marinlim

  subroutine glide_load_sigma(model,unit)

    !*FD Loads a file containing
    !*FD sigma vertical coordinates.
    use glide_types
    use glimmer_log
    use glimmer_filenames
    implicit none

    ! Arguments
    type(glide_global_type),intent(inout) :: model !*FD Ice model to use
    integer,               intent(in)    :: unit  !*FD Logical file unit to use. 
                                                  !*FD The logical file unit specified 
                                                  !*FD must not already be in use

    ! Internal variables

    integer :: up,upn
    logical :: there

    ! Beginning of code

    upn=model%general%upn

    select case(model%options%which_sigma)
    case(0)
       call write_log('Calculating sigma levels')
       do up=1,upn
          model%numerics%sigma(up) = glide_calc_sigma(real(up-1)/real(upn-1),2.)
       end do
    case(1)
       there = .false.
       inquire (exist=there,file=process_path(model%funits%sigfile))
       if (.not.there) then
          call write_log('Sigma levels file: '//trim(process_path(model%funits%sigfile))// &
               ' does not exist',GM_FATAL)
       end if
       call write_log('Reading sigma file: '//process_path(model%funits%sigfile))
       open(unit,file=process_path(model%funits%sigfile))
       read(unit,'(f5.2)',err=10,end=10) (model%numerics%sigma(up), up=1,upn)
       close(unit)
    case(2)
       call write_log('Using sigma levels from main configuration file')
    end select
    call print_sigma(model)
    return

10  call write_log('something wrong with sigma coord file',GM_FATAL)
    
  end subroutine glide_load_sigma

  function glide_calc_sigma(x,n)
      implicit none
      real :: glide_calc_sigma,x,n
      
      glide_calc_sigma = (1-(x+1)**(-n))/(1-2**(-n))
    end function glide_calc_sigma


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

    ! We set this flag to one to indicate we've got a sigfile name.
    ! A warning/error is generated if sigma levels are specified in some other way
    ! and mangle the name
    if (trim(model%funits%sigfile)/='') then
       model%funits%sigfile = filenames_inputname(model%funits%sigfile)
       model%options%which_sigma=1
    end if

  end subroutine handle_grid

  subroutine print_grid(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message


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

  ! time
  subroutine handle_time(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'tstart',model%numerics%tstart)
    call GetValue(section,'tend',model%numerics%tend)
    call GetValue(section,'dt',model%numerics%tinc)
    call GetValue(section,'ntem',model%numerics%ntem)
    call GetValue(section,'nvel',model%numerics%nvel)
    call GetValue(section,'profile',model%numerics%profile_period)
    call GetValue(section,'ndiag',model%numerics%ndiag)

  end subroutine handle_time
  
  subroutine print_time(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Time steps')
    call write_log('----------')
    write(message,*) 'start time        : ',model%numerics%tstart
    call write_log(message)
    write(message,*) 'end time          : ',model%numerics%tend
    call write_log(message)
    write(message,*) 'main time step    : ',model%numerics%tinc
    call write_log(message)
    write(message,*) 'thermal dt factor : ',model%numerics%ntem
    call write_log(message)
    write(message,*) 'velo dt factor    : ',model%numerics%nvel
    call write_log(message)
    write(message,*) 'profile frequency : ',model%numerics%profile_period
    call write_log(message)
    write(message,*) 'diag frequency    : ',model%numerics%ndiag
    call write_log(message)
    call write_log('')
  end subroutine print_time

  ! options
  subroutine handle_options(section, model)
    use glimmer_config
    use glide_types
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    call GetValue(section,'ioparams',model%funits%ncfile)
    call GetValue(section,'temperature',model%options%whichtemp)
    call GetValue(section,'flow_law',model%options%whichflwa)
    call GetValue(section,'basal_water',model%options%whichbwat)
    call GetValue(section,'marine_margin',model%options%whichmarn)
    call GetValue(section,'slip_coeff',model%options%whichbtrc)
    call GetValue(section,'evolution',model%options%whichevol)
    call GetValue(section,'vertical_integration',model%options%whichwvel)
    call GetValue(section,'topo_is_relaxed',model%options%whichrelaxed)
    call GetValue(section,'hotstart',model%options%hotstart)
    call GetValue(section,'periodic_ew',model%options%periodic_ew)
    call GetValue(section,'basal_mass_balance',model%options%basal_mbal)
  end subroutine handle_options
  
  subroutine print_options(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    ! local variables
    character(len=*), dimension(0:1), parameter :: temperature = (/ &
         'isothermal', &
         'full      '/)
    character(len=*), dimension(0:2), parameter :: flow_law = (/ &
         'Patterson and Budd               ', &
         'Patterson and Budd (temp=-10degC)', &
         'const 1e-16a^-1Pa^-n             ' /)
    character(len=*), dimension(0:2), parameter :: basal_water = (/ &
         'local water balance', &
         'local + const flux ', &
         'none               ' /)
    character(len=*), dimension(0:4), parameter :: marine_margin = (/ &
         'ignore            ', &
         'no ice shelf      ', &
         'threshold         ', &
         'const calving rate', &
         'edge threshold    '/)
    character(len=*), dimension(0:5), parameter :: slip_coeff = (/ &
         'zero           ', &
         'const          ', &
         'const if bwat>0', &
         '~basal water   ', &
         '~basal melt    ', &
         'const if T>Tpmp'/)
    character(len=*), dimension(0:2), parameter :: evolution = (/ &
         'pseudo-diffusion', &
         'ADI scheme      ', &
         'diffusion       ' /)
    character(len=*), dimension(0:1), parameter :: vertical_integration = (/ &
         'standard     ', &
         'obey upper BC' /)
    character(len=*), dimension(0:1), parameter :: b_mbal = (/ &
         'not in continutity eqn', &
         'in continutity eqn    ' /)

    call write_log('GLIDE options')
    call write_log('-------------')
    write(message,*) 'I/O parameter file      : ',trim(model%funits%ncfile)
    call write_log(message)
    if (model%options%whichtemp.lt.0 .or. model%options%whichtemp.ge.size(temperature)) then
       call write_log('Error, temperature out of range',GM_FATAL)
    end if
    write(message,*) 'temperature calculation : ',model%options%whichtemp,temperature(model%options%whichtemp)
    call write_log(message)
    if (model%options%whichflwa.lt.0 .or. model%options%whichflwa.ge.size(flow_law)) then
       call write_log('Error, flow_law out of range',GM_FATAL)
    end if
    write(message,*) 'flow law                : ',model%options%whichflwa,flow_law(model%options%whichflwa)
    call write_log(message)
    if (model%options%whichbwat.lt.0 .or. model%options%whichbwat.ge.size(basal_water)) then
       call write_log('Error, basal_water out of range',GM_FATAL)
    end if
    write(message,*) 'basal_water             : ',model%options%whichbwat,basal_water(model%options%whichbwat)
    call write_log(message)
    if (model%options%whichmarn.lt.0 .or. model%options%whichmarn.ge.size(marine_margin)) then
       call write_log('Error, marine_margin out of range',GM_FATAL)
    end if
    write(message,*) 'marine_margin           : ', model%options%whichmarn, marine_margin(model%options%whichmarn)
    call write_log(message)
    if (model%options%whichbtrc.lt.0 .or. model%options%whichbtrc.ge.size(slip_coeff)) then
       call write_log('Error, slip_coeff out of range',GM_FATAL)
    end if
    write(message,*) 'slip_coeff              : ', model%options%whichbtrc, slip_coeff(model%options%whichbtrc)
    call write_log(message)
    if (model%options%whichevol.lt.0 .or. model%options%whichevol.ge.size(evolution)) then
       call write_log('Error, evolution out of range',GM_FATAL)
    end if
    write(message,*) 'evolution               : ', model%options%whichevol, evolution(model%options%whichevol)
    call write_log(message)
    if (model%options%whichwvel.lt.0 .or. model%options%whichwvel.ge.size(vertical_integration)) then
       call write_log('Error, vertical_integration out of range',GM_FATAL)
    end if
    write(message,*) 'vertical_integration    : ',model%options%whichwvel,vertical_integration(model%options%whichwvel)
    call write_log(message)
    if (model%options%basal_mbal.lt.0 .or. model%options%basal_mbal.ge.size(b_mbal)) then
       call write_log('Error, basal_mass_balance out of range',GM_FATAL)
    end if
    write(message,*) 'basal_mass_balance      : ',model%options%basal_mbal,b_mbal(model%options%basal_mbal)
    call write_log(message)
    if (model%options%whichrelaxed.eq.1) then
       call write_log('First topo time slice is relaxed')
    end if
    if (model%options%periodic_ew.eq.1) then
       if (model%options%whichevol.eq.1) then
          call write_log('Periodic boundary conditions not implemented in ADI scheme',GM_FATAL)
       end if
       call write_log('Periodic EW lateral boundary condition')
       call write_log('  Slightly cheated with how temperature is implemented.',GM_WARNING)
    end if
    if (model%options%hotstart.eq.1) then
       call write_log('Hotstarting model')
    end if
    call write_log('')
  end subroutine print_options

  ! parameters
  subroutine handle_parameters(section, model)
    use glimmer_config
    use glide_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model
    real, pointer, dimension(:) :: temp => NULL()
    integer :: loglevel

    loglevel = GM_levels-GM_ERROR

    call GetValue(section,'log_level',loglevel)
    call glimmer_set_msg_level(loglevel)
    call GetValue(section,'ice_limit',model%numerics%thklim)
    call GetValue(section,'marine_limit',model%numerics%mlimit)
    call GetValue(section,'calving_fraction',model%numerics%calving_fraction)
    call GetValue(section,'geothermal',model%paramets%geot)
    call GetValue(section,'flow_factor',model%paramets%fiddle)
    call GetValue(section,'hydro_time',model%paramets%hydtim)
    call GetValue(section,'basal_tract',temp,5)
    if (associated(temp)) then
       model%paramets%btrac_const=temp(1)
       deallocate(temp)
    end if
    call GetValue(section,'basal_tract_const',model%paramets%btrac_const)
    call GetValue(section,'basal_tract_max',model%paramets%btrac_max)
    call GetValue(section,'basal_tract_slope',model%paramets%btrac_slope)
    call GetValue(section,'basal_water_smoothing',model%paramets%bwat_smooth)
  end subroutine handle_parameters

  subroutine print_parameters(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message

    call write_log('Parameters')
    call write_log('----------')
    write(message,*) 'ice limit             : ',model%numerics%thklim
    call write_log(message)
    write(message,*) 'marine depth limit    : ',model%numerics%mlimit
    call write_log(message)
    if (model%options%whichmarn.eq.3) then
       write(message,*) 'ice fraction lost due to calving :', model%numerics%calving_fraction
       call write_log(message)
    end if
    write(message,*) 'geothermal heat flux  : ',model%paramets%geot
    call write_log(message)
    write(message,*) 'flow enhancement      : ',model%paramets%fiddle
    call write_log(message)
    write(message,*) 'basal hydro time const: ',model%paramets%hydtim
    call write_log(message)
    if (model%options%whichbtrc.eq.1 .or. model%options%whichbtrc.eq.2 .or. model%options%whichbtrc.eq.4) then
       write(message,*) 'basal traction param  : ',model%paramets%btrac_const
       call write_log(message)
    end if
    if (model%options%whichbtrc.eq.4) then
       write(message,*) 'basal traction max  : ',model%paramets%btrac_max
       call write_log(message)
       write(message,*) 'basal traction slope  : ',model%paramets%btrac_slope
       call write_log(message)
    end if
    if (model%options%whichbtrc.eq.3) then
       write(message,*) 'basal traction factors: ',model%paramets%bpar(1)
       call write_log(message)
       write(message,*) '                        ',model%paramets%bpar(2)
       call write_log(message)
       write(message,*) '                        ',model%paramets%bpar(3)
       call write_log(message)
       write(message,*) '                        ',model%paramets%bpar(4)
       call write_log(message)
       write(message,*) '                        ',model%paramets%bpar(5)
       call write_log(message)
    end if
    write(message,*) 'basal water field smoothing strength: ',model%paramets%bwat_smooth
    call write_log(message)
    call write_log('')
  end subroutine print_parameters

  ! Sigma levels
  subroutine handle_sigma(section, model)
    use glimmer_config
    use glide_types
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    type(glide_global_type)  :: model

    if (model%options%which_sigma==1) then
       call write_log('Sigma levels specified twice - use only'// &
            ' config file or separate file, not both',GM_FATAL)
    else
       model%options%which_sigma = 2
       call GetValue(section,'sigma_levels',model%numerics%sigma,model%general%upn)
    end if

  end subroutine handle_sigma

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
       write(temp,'(F5.2)') model%numerics%sigma(i)
       message=trim(message)//trim(temp)
    enddo
    call write_log(trim(message))
    call write_log('')
    
  end subroutine print_sigma

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

  subroutine print_gthf(model)
    use glide_types
    use glimmer_log
    implicit none
    type(glide_global_type)  :: model
    character(len=100) :: message
    
    if (model%options%gthf.gt.0) then
       call write_log('GTHF configuration')
       call write_log('------------------')
       if (model%lithot%num_dim.eq.1) then
          call write_log('solve 1D diffusion equation')
       else if (model%lithot%num_dim.eq.3) then          
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

end module glide_setup
