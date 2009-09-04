! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_climate.f90 - part of the GLIMMER ice model        + 
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

  subroutine glint_downscaling(instance,g_temp,g_temp_range,g_precip,g_orog,g_zonwind, &
       g_merwind,g_humid,g_lwdown,g_swdown,g_airpress,orogflag)

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
