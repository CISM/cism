! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_climate.f90 - part of the Glimmer-CISM ice model   + 
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

module glint_climate

  use glint_type

  implicit none

  !*FD Subroutines that do various things to prepare the glint climate

  interface glint_lapserate
     module procedure glint_lapserate_dp, glint_lapserate_sp
  end interface

contains

  subroutine glint_calc_precip(instance)

    use glint_precip_param
    use glimmer_log

    !*FD Process precip if necessary

    type(glint_instance) :: instance

    select case (instance%whichprecip)
    case(1)
       ! Do nothing to the precip field
    case(2)
       ! Use the Roe/Lindzen parameterisation
       call glint_precip(instance%prcp, &
            instance%xwind, &
            instance%ywind, &
            instance%artm, &
            instance%local_orog, &
            real(instance%lgrid%delta%pt(1),rk), &
            real(instance%lgrid%delta%pt(2),rk), &
            fixed_a=.true.)
    case default
       call write_log('Invalid value of whichprecip',GM_FATAL,__FILE__,__LINE__)
    end select

    ! Convert from mm to m - very important!

    instance%prcp=instance%prcp*0.001

  end subroutine glint_calc_precip

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling(instance,                  &
                               g_temp,     g_temp_range,  &
                               g_precip,   g_orog,        &
                               g_zonwind,  g_merwind,     &
                               g_humid,    g_lwdown,      &
                               g_swdown,   g_airpress,    &
                               orogflag)

    use glint_interp

    !*FD Downscale relevant fields

    type(glint_instance) :: instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)
    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)
    logical,                intent(in)   :: orogflag

    call interp_to_local(instance%lgrid,g_temp,      instance%downs,localsp=instance%artm)
    call interp_to_local(instance%lgrid,g_temp_range,instance%downs,localsp=instance%arng,z_constrain=.true.)
    call interp_to_local(instance%lgrid,g_precip,    instance%downs,localsp=instance%prcp,z_constrain=.true.)

    if (instance%whichacab==3) then
       call interp_to_local(instance%lgrid,g_humid,   instance%downs,localrk=instance%humid,z_constrain=.true.)
       call interp_to_local(instance%lgrid,g_lwdown,  instance%downs,localrk=instance%lwdown)
       call interp_to_local(instance%lgrid,g_swdown,  instance%downs,localrk=instance%swdown)
       call interp_to_local(instance%lgrid,g_airpress,instance%downs,localrk=instance%airpress,z_constrain=.true.)
    end if

    if (orogflag) call interp_to_local(instance%lgrid,g_orog,instance%downs,localdp=instance%global_orog,z_constrain=.true.)

    if (instance%whichprecip==2.or.instance%whichacab==3) &
         call interp_wind_to_local(instance%lgrid,g_zonwind,g_merwind,instance%downs,instance%xwind,instance%ywind)

  end subroutine glint_downscaling

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_downscaling_gcm (instance,            &
                                    qsmb_g,     tsfc_g,  &
                                    topo_g,     gmask)
 
    use glimmer_paramets, only: thk0

    use glint_type, only: glint_instance
    use glint_interp, only: interp_to_local

    ! Downscale fields from the global grid (with multiple elevation classes)
    ! to the local grid.
    ! 
    ! This routine is used for downscaling when the surface mass balance is
    ! computed in the GCM land surface model.

    type(glint_instance), intent(inout) :: instance
    real(dp),dimension(:,:,:),intent(in) :: qsmb_g       ! Surface mass balance (m)
    real(dp),dimension(:,:,:),intent(in) :: tsfc_g       ! Surface temperature (C)
    real(dp),dimension(:,:,:),intent(in) :: topo_g       ! Surface elevation (m)
    integer ,dimension(:,:),  intent(in),optional :: gmask ! = 1 where global data are valid
                                                           ! = 0 elsewhere

    real(dp), parameter :: maskval = 0.0_dp    ! value written to masked out gridcells

    integer ::       &
       nec,          &      ! number of elevation classes
       i, j, n,      &      ! indices 
       nxl, nyl             ! local grid dimensions
 
    real(dp), dimension(:,:,:), allocatable ::   &
       qsmb_l,    &! interpolation of global mass balance to local grid
       tsfc_l,    &! interpolation of global sfc temperature to local grid
       topo_l      ! interpolation of global topography in each elev class to local grid

    real(dp) :: fact, usrf

    real(dp), parameter :: lapse = 0.0065_dp   ! atm lapse rate, deg/m
                                               ! used only for extrapolating temperature outside
                                               !  the range provided by the climate model
    nec = size(qsmb_g,3)
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)

    allocate(qsmb_l(nxl,nyl,nec))
    allocate(tsfc_l(nxl,nyl,nec))
    allocate(topo_l(nxl,nyl,nec))

!   Downscale global fields for each elevation class to local grid.

    if (present(gmask)) then   ! set local field = maskval where the global field is masked out
                               ! (i.e., where instance%downs%lmask = 0).
       do n = 1, nec
          call interp_to_local(instance%lgrid, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          call interp_to_local(instance%lgrid, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
          call interp_to_local(instance%lgrid, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n), &
                               gmask = gmask, maskval=maskval)
       enddo

    else    ! global field values are assumed to be valid everywhere

       do n = 1, nec
          call interp_to_local(instance%lgrid, qsmb_g(:,:,n), instance%downs, localdp=qsmb_l(:,:,n))
          call interp_to_local(instance%lgrid, tsfc_g(:,:,n), instance%downs, localdp=tsfc_l(:,:,n))
          call interp_to_local(instance%lgrid, topo_g(:,:,n), instance%downs, localdp=topo_l(:,:,n))
       enddo

    endif

    if (GLC_DEBUG) then
       write (stdout,*) ' ' 
       write (stdout,*) 'Interpolate fields to local grid'
       write (stdout,*) 'Global cell =', itest, jjtest
       do n = 1, nec
          write(stdout,*) n, topo_g(itest,jjtest, n)
       enddo

       do j = 1, nyl
       do i = 1, nxl
           if ( (instance%downs%xloc(i,j,1) == itest .and. instance%downs%yloc(i,j,1) == jjtest) .or.  &
                (instance%downs%xloc(i,j,2) == itest .and. instance%downs%yloc(i,j,2) == jjtest) .or.  &
                (instance%downs%xloc(i,j,3) == itest .and. instance%downs%yloc(i,j,3) == jjtest) .or.  &
                (instance%downs%xloc(i,j,4) == itest .and. instance%downs%yloc(i,j,4) == jjtest) ) then
               write(stdout,*) i, j, thk0 * instance%model%geometry%usrf(i,j)
           endif
       enddo
       enddo
    
       i = itest_local
       j = jtest_local
       write (stdout,*) ' ' 
       write (stdout,*) 'Interpolated to local cells: i, j =', i, j
       do n = 1, nec
          write (stdout,*) ' '
          write (stdout,*) 'n =', n
          write (stdout,*) 'qsmb_l =', qsmb_l(i,j,n)
          write (stdout,*) 'tsfc_l =', tsfc_l(i,j,n)
          write (stdout,*) 'topo_l =', topo_l(i,j,n)
       enddo
    endif

!   Interpolate tsfc and qsmb to local topography using values in the neighboring 
!    elevation classes.
!   If the local topography is outside the bounds of the global elevations classes,
!    extrapolate the temperature using the prescribed lapse rate.

    do j = 1, nyl
    do i = 1, nxl

       usrf = instance%model%geometry%usrf(i,j) * thk0   ! actual sfc elevation (m)

       if (usrf <= topo_l(i,j,1)) then
          instance%acab(i,j) = qsmb_l(i,j,1)
          instance%artm(i,j) = tsfc_l(i,j,1) + lapse*(topo_l(i,j,1)-usrf)
       elseif (usrf > topo_l(i,j,nec)) then
          instance%acab(i,j) = qsmb_l(i,j,nec)
          instance%artm(i,j) = tsfc_l(i,j,nec) - lapse*(usrf-topo_l(i,j,nec))
       else
          do n = 2, nec
             if (usrf > topo_l(i,j,n-1) .and. usrf <= topo_l(i,j,n)) then
                fact = (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
                instance%acab(i,j) = fact*qsmb_l(i,j,n-1) + (1._dp-fact)*qsmb_l(i,j,n)
                instance%artm(i,j) = fact*tsfc_l(i,j,n-1) + (1._dp-fact)*tsfc_l(i,j,n)
                exit
             endif
          enddo
       endif   ! usrf

       if (GLC_DEBUG) then
          if (i==itest_local .and. j==jtest_local) then
             n = 4  
             write (stdout,*) ' '
             write (stdout,*) 'Interpolated values, i, j, n =', i, j, n
             write (stdout,*) 'usrf =', usrf
             write (stdout,*) 'acab =', instance%acab(i,j)
             write (stdout,*) 'artm =', instance%artm(i,j)
             write (stdout,*) 'topo(n-1) =', topo_l(i,j,n-1)
             write (stdout,*) 'topo(n) =', topo_l(i,j,n)
             write (stdout,*) 'qsmb(n-1) =', qsmb_l(i,j,n-1)
             write (stdout,*) 'qsmb(n) =', qsmb_l(i,j,n)
             write (stdout,*) 'tsfc(n-1) =', tsfc_l(i,j,n-1)
             write (stdout,*) 'tsfc(n) =', tsfc_l(i,j,n)
             write (stdout,*) 'fact = ', (topo_l(i,j,n) - usrf) / (topo_l(i,j,n) - topo_l(i,j,n-1)) 
          endif
       endif

    enddo  ! i
    enddo  ! j

  end subroutine glint_downscaling_gcm

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_lapserate_dp(temp,topo,lr)

    !*FD Corrects the temperature field
    !*FD for height, using a constant lapse rate.
    !*FD
    !*FD This the double-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(dp),dimension(:,:), intent(inout) :: temp !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(rk),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(rk),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp=temp-(lr*topo/1000.0)                     ! The lapse rate calculation.

  end subroutine glint_lapserate_dp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_lapserate_sp(temp,topo,lr)

    !*FD Corrects the temperature field
    !*FD for height, using a constant lapse rate.
    !*FD
    !*FD This the single-precision version, aliased as \texttt{glimmer\_lapserate}.

    implicit none

    real(sp),dimension(:,:),intent(inout) :: temp  !*FD temperature at sea-level in $^{\circ}$C
                                                   !*FD used for input and output
    real(rk),dimension(:,:), intent(in)    :: topo !*FD topography field (m above msl)
    real(rk),                intent(in)    :: lr   !*FD Lapse rate ($^{\circ}\mathrm{C\,km}^{-1}$).
                                                   !*FD
                                                   !*FD NB: the lapse rate is positive for 
                                                   !*FD falling temp with height\ldots

    temp=temp-(lr*topo/1000.0)                     ! The lapse rate calculation.

  end subroutine glint_lapserate_sp

end module glint_climate
