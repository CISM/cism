! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eismint3_forcing.f90 - part of the GLIMMER ice model     + 
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

module eismint3_forcing

  use glimmer_global, only: sp
  use eismint3_types

  implicit none

  !*FD Provides climate forcing for EISMINT3 Greeland scenario 1.
  !*FD N.B. Input precip should be water-equivalent not ice equivalent.

contains

  subroutine eismint3_initialise(climate,config,model)

    !*FD initialise EISMINT3 climate

    use glide_types
    use glimmer_config
    use eismint3_io
    use glide_io
    use glide_setup
    use glimmer_log
    use physcon, only: rhoi,rhow
    use glimmer_paramets, only: thk0

    implicit none

    type(eismint3_climate) :: climate      !*FD structure holding EISMINT3 climate
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file   
    type(glide_global_type) :: model       !*FD model instance

    call eismint3_readconfig(climate,config)
    call eismint3_printconfig(climate)

    call coordsystem_allocate(model%general%ice_grid,climate%prcp)
    call coordsystem_allocate(model%general%ice_grid,climate%artm)
    call coordsystem_allocate(model%general%ice_grid,climate%arng)
    call coordsystem_allocate(model%general%ice_grid,climate%acab)
    call coordsystem_allocate(model%general%ice_grid,climate%usrf)
    call coordsystem_allocate(model%general%ice_grid,climate%ablt)
    call coordsystem_allocate(model%general%ice_grid,climate%presusurf)
    call coordsystem_allocate(model%general%ice_grid,climate%presprcp)
    call coordsystem_allocate(model%general%ice_grid,climate%presartm)
    call coordsystem_allocate(model%general%ice_grid,climate%landsea)

    call eismint3_io_readall(climate,model)
    call glimmer_pdd_init(climate%pdd_scheme,config)

    select case(climate%loadthk)
    case(0)
       ! Calculate initial thickness.

       where (climate%presusurf>0.0)
          climate%landsea=.true.
       elsewhere
          climate%landsea=.false.
       end where

       call eismint3_temp(climate%artm,climate%arng,climate%presusurf,model%climate%lati,0.0_sp)
       call glimmer_pdd_mbal(climate%pdd_scheme,climate%artm,climate%arng, &
            climate%presprcp,climate%ablt,climate%acab,climate%landsea)

       ! Convert to ice-equivalent depth

       climate%acab=climate%acab*(rhow/rhoi)

       ! Put it into glide

       call glide_set_thk(model,max(0.0,climate%acab*model%numerics%tinc))
    case(1)
       ! do nothing
    case(2)
       call write_log('Unknown value for climate%loadthk',GM_FATAL)
    end select

    ! Calculate the initial ice surface
    ! The reference present-day surface (for precip enhancement calc) is already loaded from file
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus,model%geometry%lsrf)
    model%geometry%usrf = model%geometry%thck + model%geometry%lsrf

    ! Work out present-day temperature for precip enhancement calculation

    call eismint3_temp(climate%presartm,climate%arng,climate%presusurf,model%climate%lati,0.0_sp)

  end subroutine eismint3_initialise

  subroutine eismint3_readconfig(climate,config)
    !*FD read configuration
    use glimmer_log
    use glimmer_config
    implicit none
    type(eismint3_climate) :: climate         !*FD structure holding climate info
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    call GetSection(config,section,'EISMINT-3')
    if (associated(section)) then
       call GetValue(section,'load_thk',climate%loadthk)
       call GetValue(section,'precip_enhance',climate%pfac)
       call GetValue(section,'temp_perturb',climate%temp_perturb)
    end if

  end subroutine eismint3_readconfig

  subroutine eismint3_printconfig(climate)
    !*FD print eismint3 climate configuration
    use glimmer_log
    implicit none
    type(eismint3_climate) :: climate   !*FD structure holding climate info
    character(len=100) :: message

    call write_log_div
    call write_log('EISMINT-3 Greenland configuration')
    call write_log('------------------------------------')
    write(message,*) 'load thickness            : ',climate%loadthk
    call write_log(message)
    write(message,*) 'precip enhancement factor : ',climate%pfac
    call write_log(message)
    write(message,*) 'temperature perturbation  : ',climate%temp_perturb,' degC'
    call write_log(message)
    call write_log('')

  end subroutine eismint3_printconfig

  subroutine eismint3_clim(climate,model)

    use glide_types
    use glide_io
    use physcon, only: rhoi,rhow

    type(eismint3_climate) :: climate      !*FD structure holding EISMINT3 climate
    type(glide_global_type) :: model       !*FD model instance

    call glide_get_usurf(model,climate%usrf)
    where (climate%usrf>0.0)
       climate%landsea=.true.
    elsewhere
       climate%landsea=.false.
    end where

    call eismint3_temp(climate%artm,climate%arng,climate%usrf,model%climate%lati,climate%temp_perturb)
    call eismint3_prcp(climate%prcp,climate%artm,climate%presartm,climate%presprcp,climate%pfac)
    call glimmer_pdd_mbal(climate%pdd_scheme,climate%artm,climate%arng,climate%prcp,climate%ablt,climate%acab,climate%landsea)

    where (.not.climate%landsea) climate%acab=0.0

    ! Convert to ice-equivalent depth
    climate%acab=climate%acab*(rhow/rhoi)

    call glide_set_acab(model,climate%acab)
    call glide_set_artm(model,climate%artm)

  end subroutine eismint3_clim

  subroutine eismint3_prcp(prcp,artm,presartm,presprcp,pfac)

    real(sp),dimension(:,:),intent(out) :: prcp
    real(sp),dimension(:,:),intent(in)  :: artm,presartm,presprcp
    real(sp) :: pfac

    prcp = presprcp * pfac**(artm - presartm)

  end subroutine eismint3_prcp

  subroutine eismint3_temp(artm,arng,usrf,lati,perturb)

    real(sp),dimension(:,:),intent(out) :: artm,arng
    real(sp),dimension(:,:),intent(in)  :: usrf,lati
    real(sp),               intent(in)  :: perturb

    artm=49.13-0.007992*max(usrf,20*(lati-65.0))-0.7576*lati+perturb
    arng=30.78-0.006277*usrf-0.3262*lati+perturb-artm

  end subroutine eismint3_temp

end module eismint3_forcing
