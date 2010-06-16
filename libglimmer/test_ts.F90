! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  test_ts.f90 - part of the Glimmer-CISM ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009
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
