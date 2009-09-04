! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_constants.f90 - part of the GLIMMER ice model      + 
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

module glint_constants

  use glimmer_global

  implicit none

  ! ------------------------------------------------------------
  ! global parameters
  ! ------------------------------------------------------------

  real(rk),parameter :: pi=3.141592654          !*FD The value of pi
  real(rk),parameter :: days2hours=24.0
  real(rk),parameter :: hours2seconds=3600.0    !*FD Hours to seconds conversion factor

  integer, parameter :: default_diy=360                    !*FD Default number of days in year
  integer, parameter :: default_y2h=days2hours*default_diy !*FD Default years to hours conversion

  ! Constants set at run-time

  integer  :: days_in_year=default_diy        !*FD The number of days in a year  
  real(rk) :: years2hours =default_y2h        !*FD Years to hours conversion factor
  real(rk) :: hours2years =1.0_rk/default_y2h !*FD Hours to years conversion factor

  private :: default_diy,default_y2h

contains

  subroutine glint_set_year_length(daysinyear)

    integer, intent(in) :: daysinyear

    days_in_year=daysinyear
    years2hours=days2hours*days_in_year 
    hours2years=1.0_rk/years2hours      

  end subroutine glint_set_year_length

end module glint_constants
