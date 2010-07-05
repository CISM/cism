! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                               +
! +  erosion_trans_type.f90 - part of the Glimmer-CISM ice model  + 
! +                                                               +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2005, 2006, 2010
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
! module defining er_transport_type for conservation of 2nd moment advection scheme

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module erosion_transport_type
  use advect_2ndmo
  use glimmer_global, only : rk
  type er_transport_type
     type(advect_type) :: mo_seds1
     type(advect_type) :: mo_seds2
  end type er_transport_type
end module erosion_transport_type
