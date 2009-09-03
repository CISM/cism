! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  test_ts.f90 - part of the GLIMMER ice model              + 
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

program test_ts
  !*FD testing the time series module
  use glimmer_ts
  implicit none

  character(len=50) :: fname
  integer numv
  real, dimension(:), allocatable :: val1,val2
  real time
  type(glimmer_tseries) :: ts

  write(*,*) 'Enter name containing time series'
  read(*,*) fname
  write(*,*) 'Enter number of values per time'
  read(*,*) numv
  
  allocate(val1(numv))
  allocate(val2(numv))
  
  call glimmer_read_ts(ts,fname,numv)
  
  do
     write(*,*) 'Enter new time'
     read(*,*) time
     call glimmer_ts_step(ts,time,val1)
     call glimmer_ts_linear(ts,time,val2)
     write(*,*) val1,val2
  end do
end program test_ts
