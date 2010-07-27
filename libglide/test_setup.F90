! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  test_setup.f90 - part of the Glimmer-CISM ice model      +
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

program testsetup
  !*FD testing setup module
  !*FD written by Magnus Hagdorn, June 2004
  use glide_types
  use glide_setup
  use glimmer_config
  use glimmer_log
  implicit none

  type(glide_global_type) :: model        ! model instance
  character(len=20) :: fname   ! name of paramter file
  type(ConfigSection), pointer :: config  ! configuration stuff

  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname

  ! start logging
  call open_log(unit=50)
  
  ! read configuration
  call ConfigRead(fname,config)  

  call glide_readconfig(model,config)
  call glide_printconfig(model)
end program testsetup
