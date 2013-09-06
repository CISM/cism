!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_driver.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
program cism_driver

  use parallel
  use glimmer_commandline
  use cism_front_end_module 

  call parallel_initialise
  call glimmer_GetCommandline()

  call cism_front_end()
!  call parallel_finalise
end program cism_driver

! need in LD_LIBRARY_PATH: /opt/cray/xt-tools/lgdb/1.4/lib/alps
