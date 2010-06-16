! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy.f90 - part of the Glimmer-CISM ice model        + 
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

module isostasy
  
  !*FD calculate isostatic adjustment due to changing surface loads

  use glimmer_global, only : dp
  use isostasy_setup
  use isostasy_types
  use isostasy_el

  private :: relaxing_mantle
  
contains
  subroutine init_isostasy(model)
    !*FD initialise isostasy calculations
    use glide_types
    use glimmer_physcon,  only: scyr
    use glimmer_paramets, only: tim0
    implicit none
    type(glide_global_type) :: model

    if (model%isos%lithosphere .eq. 1) then
       call init_elastic(model%isos%rbel,model%numerics%dew)
    end if
    model%isos%next_calc = model%numerics%tstart

    ! scale tau
    model%isos%relaxed_tau = model%isos%relaxed_tau * scyr / tim0

  end subroutine init_isostasy
  
  subroutine isos_icewaterload(model)
    !*FD calculate surface load factors due to water and ice distribution
    use glimmer_physcon
    use glide_types
    implicit none
    type(glide_global_type) :: model

    real(kind=dp) :: ice_mass, water_depth, water_mass
    integer :: ew,ns
  
     do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          ice_mass = rhoi * model%geometry%thck(ew,ns)
          if (model%geometry%topg(ew,ns)-model%climate%eus.lt.0) then             ! check if we are below sea level
             water_depth = model%climate%eus - model%geometry%topg(ew,ns)
             water_mass = rhoo * water_depth
             ! Just the water load due to changes in sea-level
             model%isos%load_factors(ew,ns) = rhoo* model%climate%eus/rhom
             ! Check if ice is not floating
             if ( ice_mass .gt. water_mass ) then
                model%isos%load_factors(ew,ns) = model%isos%load_factors(ew,ns) + (ice_mass - water_mass)/rhom
             end if
          else                                       ! bedrock is above sea level
             model%isos%load_factors(ew,ns) = ice_mass/rhom
          end if
       end do
    end do
  end subroutine isos_icewaterload

  subroutine isos_isostasy(model)
    !*FD calculate isostatic adjustment due to changing surface loads
    use glide_types
    implicit none
    type(glide_global_type) :: model

    ! update load if necessary
    if (model%isos%new_load) then
       call isos_lithosphere(model,model%isos%load,model%isos%load_factors)
       ! update bed rock with (non-viscous) fluid mantle
       if (model%isos%asthenosphere .eq. 0) then
          model%geometry%topg = model%isos%relx - model%isos%load
       end if
       model%isos%new_load = .false.
    end if
    ! update bed rock with relaxing mantle
    if (model%isos%asthenosphere .eq. 1) then
       call relaxing_mantle(model)
    end if
  end subroutine isos_isostasy

  subroutine isos_lithosphere(model,load,load_factors)
    use glide_types
    implicit none
    type(glide_global_type) :: model
    real(kind=dp), dimension(:,:), intent(out) :: load !*FD loading effect due to load_factors
    real(kind=dp), dimension(:,:), intent(in)  :: load_factors !*FD load mass divided by mantle density

    if (model%isos%lithosphere .eq. 0) then
       ! local lithosphere
       load = load_factors
    else if (model%isos%lithosphere .eq. 1) then
       call calc_elastic(model%isos%rbel,load,load_factors)
    end if
  end subroutine isos_lithosphere

  subroutine isos_relaxed(model)
    !*FD Calculate the relaxed topography, assuming the isostatic depression
    !*FD is the equilibrium state for the current topography.
    use glide_types
    implicit none
    type(glide_global_type) :: model

    ! Calculate the load
    call isos_icewaterload(model)
    ! Apply lithosphere model
    call isos_lithosphere(model,model%isos%load,model%isos%load_factors)
    ! Add to present topography to get relaxed topography
    model%isos%relx = model%geometry%topg + model%isos%load

  end subroutine isos_relaxed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine relaxing_mantle(model)
    !*FD approximate mantle with a relaxing half-space: dh/dt=-1/tau*(w-h)
    use glide_types
    implicit none
    type(glide_global_type) :: model
    
    integer :: ew,ns
    real(kind=dp) :: ft1, ft2

    ft1 = exp(-model%numerics%dt/model%isos%relaxed_tau)
    ft2 = 1. - ft1
    
    do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          model%geometry%topg(ew,ns) = ft2*(model%isos%relx(ew,ns)-model%isos%load(ew,ns)) + ft1*model%geometry%topg(ew,ns)
       end do
    end do
  end subroutine relaxing_mantle

end module isostasy
