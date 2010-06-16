! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy_types.f90 - part of the Glimmer-CISM ice model  + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
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

module isostasy_types
  !*FD types for isostasy model

  use glimmer_global, only : dp
  
  type isostasy_elastic
     !*FD Holds data used by isostatic adjustment calculations
     
     real(kind=dp) :: d = 0.24e25                !*FD flexural rigidity
     real(kind=dp) :: lr                         !*FD radius of relative stiffness
     real(kind=dp) :: a                          !*FD radius of disk
     real(kind=dp) :: c1,c2,cd3,cd4              !*FD coefficients
     real(kind=dp), dimension(:,:), pointer :: w !*FD matrix operator for lithosphere deformation
     integer :: wsize                            !*FD size of operator (0:rbel_wsize, 0:rbel_wsize), operator is axis symmetric
  end type isostasy_elastic

  type isos_type
     !*FD contains isostasy configuration
     logical :: do_isos = .False.    
     !*FD set to .True. if isostatic adjustment should be handled
     integer :: lithosphere = 0      
     !*FD method for calculating equilibrium bedrock depression:
     !*FD \begin{description} 
     !*FD \item[0] local lithosphere, equilibrium bedrock depression is found using Archimedes' principle
     !*FD \item[1] elastic lithosphere, flexural rigidity is taken into account
     !*FD \end{description}
     integer :: asthenosphere = 0
     !*FD method for approximating the mantle
     !*FD \begin{description} 
     !*FD \item[0] fluid mantle, isostatic adjustment happens instantaneously
     !*FD \item[1] relaxing mantle, mantle is approximated by a half-space
     !*FD \end{description}
     real :: relaxed_tau = 4000.
     !*FD characteristic time constant of relaxing mantle
     real :: period = 500.
     !*FD lithosphere update period
     real :: next_calc
     !*FD when to upate lithosphere
     logical :: new_load=.false.
     !*FD set to true if there is a new surface load
     type(isostasy_elastic) :: rbel
     !*FD structure holding elastic lithosphere setup

     real(dp),dimension(:,:),pointer :: relx => null()
     !*FD The elevation of the relaxed topography, by \texttt{thck0}.
     real(dp),dimension(:,:),pointer :: load => null()
     !*FD the load imposed on lithosphere
     real(dp),dimension(:,:),pointer :: load_factors => null()
     !*FD temporary used for load calculation
  end type isos_type  

contains

   subroutine isos_allocate(isos, ewn, nsn)
    !*FD allocate data for isostasy calculations
    implicit none
    type(isos_type) :: isos                !*FD structure holding isostasy configuration
    integer, intent(in) :: ewn, nsn        !*FD size of grid

    allocate(isos%relx(ewn,nsn))
    allocate(isos%load(ewn,nsn))
    allocate(isos%load_factors(ewn,nsn))
  end subroutine isos_allocate

  subroutine isos_deallocate(isos)
    !*FD deallocate data for isostasy calculations
    implicit none
    type(isos_type) :: isos                !*FD structure holding isostasy configuration
    
    deallocate(isos%relx)
    deallocate(isos%load)
    deallocate(isos%load_factors)
  end subroutine isos_deallocate
end module isostasy_types
