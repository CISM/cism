!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   isostasy.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module isostasy

  !-------------------------------------------------------------------------
  ! Some notes on the isostasy calculation (WHL, May 2017):
  !
  ! For the most part, the isostasy has not changed since the original Glimmer release.
  ! The major change in CISM2.1 is to enable the elastic lithosphere calculation
  !  in simulations with more than one task. This is done in a simple way, by gathering
  !  load factors to the main task, doing a serial calculation, and then scattering the
  !  resulting load back to the local tasks.
  !
  ! The following config settings are relevant to the isostasy:
  ! (1) To run with isostasy, set isostasy = 1 in the [options] section.
  !     The default is 0 (no isostasy).
  ! (2) To run with an elastic lithosphere, set lithosphere = 1 in the [isostasy] section.
  !     This is now the default value, so it no longer needs to be set explicitly in the config file.
  !     Note that on multiple tasks, this calculation requires a gather/scatter that
  !     does not scale well.  It seems sufficiently fast, though, on a 4-km mesh.
  !     The alternative is a local lithosphere (lithosphere = 0) that is less realistic.
  ! (3) To run with a relaxing asthenosphere, set asthenosphere = 1 in the [isostasy] section.
  !     This is now the default and does not need to be set explicitly.
  !     The alternative is a fluid asthenosphere (asthenosphere = 0) with instantaneous
  !     isostatic adjustment, which is less realistic.
  ! (4) The flexural rigidity of the elastic lithosphere is controlled by the parameter 'flexural_rigidity',
  !     which can be set in the [isostasy] section. The default is 0.24e25 N m.
  ! (5) The period for recomputing the load in the elastic lithosphere calculation is controlled
  !     by the parameter 'lithosphere_period', which can be set in the [isostasy] section.
  !     The default is 500 yr.  As long as the load is not recomputed too often, the isostasy
  !     calculation should have minimal cost compared to the whole simulation
  !     (at least on grids of moderate resolution, ~4 km).
  ! (6) The adjustment time scale in the relaxing asthenosphere calculation is controlled
  !     by the parameter relaxed_tau, which can be set in the [isostasy] section.
  !     The default is 4000 yr.
  !
  ! Finally, a few words on the 'whichrelaxed' parameter.  This used to be called 'topo_is_relaxed'
  ! and was in the [options] section; now it is called 'whichrelaxed' and is in the [isostasy] section.
  ! There are three possible values:
  !
  ! - whichrelaxed = 0, the default setting. In this case, both topg and relx, if present, are read
  !   from the input file. The model topography is initialized as topg.  The relx field is interpreted
  !   as the topography we would have eventually (after the asthenosphere fully relaxes) with zero load.
  !   The asthenosphere calculation continually adjusts the topography toward a state with topg = relx - load.
  !   NOTE: If relx is not present in the input file, the model will be initialized with relx = 0
  !         everywhere, which may be OK for idealized problems but will be wrong for real ice sheets.
  !
  ! - whichrelaxed = 1. In this case, the input 'topg' field is interpreted as the relaxed field.
  !   That is, the model sets relx = topg at initialization.  Then topg will be correct if there is no load
  !   (e.g., prior to ice sheet inception), but in general will be wrong. If relx is different from
  !   the initial topography, it is better to input each field separately with whichrelaxed = 0.
  !
  ! - whichrelaxed = 2. In this case, the input 'topg' field is interpreted as the equilibrium topography.
  !   The field 'relx' (i.e., the steady-state topography with zero load) is computed at initialization
  !   as relx = topg + load. This setting could be useful if we happen to know the equilibrium value
  !   of topg and want to compute relx. But if the model is stopping and restarting, the interpretation
  !   of topg as the equilibrium topography will usually be wrong on restart.
  !
  ! In general, the preferred setting is whichrelaxed = 0, with topg and relx read in separately
  ! from the input file. The other settings have specialized uses but may be inappropriate for production.
  !-------------------------------------------------------------------------


  !> calculate isostatic adjustment due to changing surface loads

  use glimmer_global, only : dp
  use isostasy_elastic

  implicit none

  private :: relaxing_mantle

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------

  subroutine init_isostasy(model)

    !> initialise isostasy calculations
    use parallel
    use glide_types
    use glimmer_physcon,  only: scyr
    use glimmer_paramets, only: tim0
    implicit none

    type(glide_global_type) :: model

    if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then

       call init_elastic(model%isostasy%rbel, model%numerics%dew)

    end if

    !-----------------------------------------------------------------
    ! Based on the update period, determine how frequently the lithosphere load should be updated.
    ! The load is updated every nlith timesteps.
    ! An integer is used instead of a real number to decide when to update, in order to avoid roundoff issues.
    ! NOTE: The ratio isostasy%period/tinc is rounded to the nearest integer.
    !       Use numerics%tinc because it has units of years (like isostasy%period), whereas numerics%dt has model timeunits.
    !-----------------------------------------------------------------

    if (model%isostasy%period > 0.0d0) then
       model%isostasy%nlith = nint(model%isostasy%period / model%numerics%tinc)
    else
       model%isostasy%nlith = 0  ! never update
    endif

    model%isostasy%relaxed_tau = model%isostasy%relaxed_tau * scyr / tim0

  end subroutine init_isostasy

!-------------------------------------------------------------------------
  
  subroutine isos_icewaterload(model)

    !> calculate surface load factors due to water and ice distribution

    use glimmer_physcon
    use glide_types
    implicit none

    type(glide_global_type) :: model

    real(dp) :: ice_mass, water_depth, water_mass
    integer :: ew,ns
  
     do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          ice_mass = rhoi * model%geometry%thck(ew,ns)

          if (model%geometry%topg(ew,ns) - model%climate%eus < 0.d0) then   ! check if we are below sea level

             water_depth = model%climate%eus - model%geometry%topg(ew,ns)
             water_mass = rhoo * water_depth

             ! Just the water load due to changes in sea-level
             model%isostasy%load_factors(ew,ns) = rhoo* model%climate%eus/rhom

             ! Check if ice is not floating
             if ( ice_mass > water_mass ) then
                model%isostasy%load_factors(ew,ns) = model%isostasy%load_factors(ew,ns) + (ice_mass - water_mass)/rhom
             end if

          else                                       ! bedrock is above sea level

             model%isostasy%load_factors(ew,ns) = ice_mass/rhom

          end if

       end do
    end do

  end subroutine isos_icewaterload

!-------------------------------------------------------------------------

  subroutine isos_compute(model)

    !> calculate isostatic adjustment due to changing surface loads

    use glide_types
    implicit none

    type(glide_global_type) :: model

    ! update load if it is time to do so
    if (model%isostasy%new_load) then

       call isos_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

       ! update bedrock if the mantle is fluid (non-viscous)
       if (model%isostasy%asthenosphere == ASTHENOSPHERE_FLUID) then
          model%geometry%topg = model%isostasy%relx - model%isostasy%load
       end if

       model%isostasy%new_load = .false.

    end if

    ! update bedrock if the mantle is relaxing
    if (model%isostasy%asthenosphere == ASTHENOSPHERE_RELAXING) then
       call relaxing_mantle(model)
    end if

  end subroutine isos_compute

!-------------------------------------------------------------------------

  subroutine isos_lithosphere(model,load,load_factors)

    use glide_types
    implicit none

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(out) :: load !> loading effect due to load_factors
    real(dp), dimension(:,:), intent(in)  :: load_factors !> load mass divided by mantle density

    if (model%isostasy%lithosphere == LITHOSPHERE_LOCAL) then

       load = load_factors

    else if (model%isostasy%lithosphere == LITHOSPHERE_ELASTIC) then

       call calc_elastic(&
            model%isostasy%rbel,  &
            load_factors,         &
            load,                 &
            model%numerics%idiag, &
            model%numerics%jdiag, &
            model%numerics%idiag_local, &
            model%numerics%jdiag_local, &
            model%numerics%rdiag_local)

    end if

  end subroutine isos_lithosphere

!-------------------------------------------------------------------------

  subroutine isos_relaxed(model)

    !> Calculate the relaxed topography, assuming the isostatic depression
    !> is the equilibrium state for the current topography.

    use glide_types
    implicit none
    type(glide_global_type) :: model

    ! Calculate the load
    call isos_icewaterload(model)

    ! Apply lithosphere model
    call isos_lithosphere(model, model%isostasy%load, model%isostasy%load_factors)

    ! Add to present topography to get relaxed topography
    model%isostasy%relx = model%geometry%topg + model%isostasy%load

  end subroutine isos_relaxed

!-------------------------------------------------------------------------
! private subroutines
!-------------------------------------------------------------------------

  subroutine relaxing_mantle(model)

    !> approximate mantle with a relaxing half-space: dh/dt=-1/tau*(w-h)
    use glide_types
    implicit none
    type(glide_global_type) :: model
    
    integer :: ew,ns
    real(dp) :: ft1, ft2

    ft1 = exp(-model%numerics%dt/model%isostasy%relaxed_tau)
    ft2 = 1.d0 - ft1

    do ns=1,model%general%nsn
       do ew=1,model%general%ewn
          model%geometry%topg(ew,ns) = ft2 * (model%isostasy%relx(ew,ns) - model%isostasy%load(ew,ns)) &
                                     + ft1 *  model%geometry%topg(ew,ns)
       end do
    end do

  end subroutine relaxing_mantle

!-------------------------------------------------------------------------

end module isostasy

!-------------------------------------------------------------------------
