! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  test_commandline.f90 - part of the Glimmer-CISM ice model+ 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2009
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

program test_commandline
  use glimmer_commandline
  implicit none

  call glimmer_GetCommandline()
  call glimmer_PrintCommandline()
end program test_commandline
