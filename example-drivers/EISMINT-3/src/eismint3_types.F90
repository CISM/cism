! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eismint3_types.f90 - part of the GLIMMER ice model       + 
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

module eismint3_types

  use glimmer_global, only: sp
  use glimmer_pdd

  implicit none

  type eismint3_climate
     real(sp),dimension(:,:),pointer :: prcp !*FD Annual accumulation
     real(sp),dimension(:,:),pointer :: acab !*FD Mass-balance
     real(sp),dimension(:,:),pointer :: artm !*FD Surface temp
     real(sp),dimension(:,:),pointer :: arng !*FD Surface temp half-range
     real(sp),dimension(:,:),pointer :: usrf !*FD Surface elevation
     real(sp),dimension(:,:),pointer :: ablt !*FD Ablation
     real(sp),dimension(:,:),pointer :: presusurf !*FD Present-day upper-surface elevation
     real(sp),dimension(:,:),pointer :: presartm  !*FD Present-day surface temperature
     real(sp),dimension(:,:),pointer :: presprcp  !*FD Present-day precip (water-equivalent)
     logical,dimension(:,:),pointer :: landsea !*FD Land-sea mask
     type(glimmer_pdd_params) :: pdd_scheme
     integer :: loadthk=0 !*FD Load thickness from file or start from scratch
     real(sp) :: pfac=1.0533 !*FD Precip enhancement factor (default is supposed EISMINT value)
     real(sp) :: temp_perturb = 0.0 !*FD Climate temperature perturbation
  end type eismint3_climate

end module eismint3_types
