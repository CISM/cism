! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eis_types.f90 - part of the Glimmer-CISM ice model       + 
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

module eis_types
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  
  use eis_ela
  use eis_temp
  use eis_slc
  use eis_cony

  type eis_climate_type
     !*FD holds EIS climate
     type(eis_ela_type)  :: ela  !*FD ELA forcing
     type(eis_temp_type) :: temp !*FD temperature forcing
     type(eis_slc_type)  :: slc  !*FD sea level forcing
     type(eis_cony_type) :: cony !*FD continentality forcing
  end type eis_climate_type
end module eis_types
