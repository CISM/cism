!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_diagnostics.F90 - part of the Community Ice Sheet Model (CISM)  
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

!TODO - Calculations of iarea, iareaf and areag in calc_iareaf_iareag() and glide_set_mask() could be replaced by values computed here.  
!       These could be saved to the model derived type (model%geometry%iarea, etc.) for output.

module glide_diagnostics

  ! subroutines for computing various useful diagnostics
  ! Author: William Lipscomb, LANL 
 
  use glimmer_global, only: dp
  use glimmer_log
  use glide_types
  use parallel

  implicit none

contains

  subroutine glide_write_diagnostics (model,  time,    &
                                      tstep_count)

    ! Short driver subroutine to decide whether it's time to write diagnostics.
    ! If so, it calls glide_write_diag.   

    ! input/output arguments

    ! Note: model is intent(inout) so that some global diagnostics can be computed below
    type(glide_global_type), intent(inout) :: model    ! model instance

    real(dp), intent(in) :: time          ! current time in years

    integer,  intent(in) :: tstep_count   ! current timestep

    ! local arguments

    logical, parameter :: verbose_diagnostics = .false.

    ! debug
    if (main_task .and. verbose_diagnostics) then
       print*, '	'
       print*, 'In glide_write_diagnostics'
       print*, 'time =', time
       print*, 'dt_diag =', model%numerics%dt_diag
       print*, 'ndiag =', model%numerics%ndiag
       print*, 'tstep_count =', tstep_count
    endif

    if (model%numerics%ndiag > 0) then

       if (mod(tstep_count, model%numerics%ndiag) == 0)  then    ! time to write
          call glide_write_diag(model, time)
       endif

    endif    ! ndiag > 0

  end subroutine glide_write_diagnostics
 
!--------------------------------------------------------------------------

  subroutine glide_init_diag (model)

    ! Initialize model diagnostics for glide or glissade.
    ! (1) Set ndiag based on dt_diag.  (Diagnostics are written every ndiag steps.)
    ! (2) Find the local rank and indices of the global diagnostic point


    implicit none

    ! input/output arguments

    type(glide_global_type), intent(inout) :: model    ! model instance

    ! local variables

    character(len=100) :: message

    !-----------------------------------------------------------------
    ! Given dt_diag, compute the interval ndiag of diagnostic output.
    ! (Output is written every ndiag timesteps.)
    ! NOTE: The ratio dt_diag/tinc is rounded to the nearest integer.
    !-----------------------------------------------------------------

    if (model%numerics%dt_diag > 0.0d0) then   ! dt_diag was specified in the config file; use it to compute ndiag

       ! Note: tinc and dt_diag have units of years, whereas dt has model timeunits
       model%numerics%ndiag = nint(model%numerics%dt_diag / model%numerics%tinc)
       model%numerics%ndiag = max(model%numerics%ndiag, 1)  ! cannot write more often than once per timestep

    endif

    !-----------------------------------------------------------------
    ! Find the local rank and indices of the global diagnostic point
    !-----------------------------------------------------------------

    call parallel_localindex(model%numerics%idiag,       model%numerics%jdiag,       &
                             model%numerics%idiag_local, model%numerics%jdiag_local, &
                             model%numerics%rdiag_local)

    !WHL - debug
    if (main_task) then
       write(6,'(a25,2i6)') 'Global idiag, jdiag:     ',   &
                             model%numerics%idiag, model%numerics%jdiag
       write(6,'(a25,3i6)') 'Local idiag, jdiag, task:',   &
                             model%numerics%idiag_local,  &
                             model%numerics%jdiag_local,  &
                             model%numerics%rdiag_local
    endif

    if (main_task) then

       write(message,'(a25,2i6)') 'Global idiag, jdiag:     ',   &
                                   model%numerics%idiag, model%numerics%jdiag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,3i6)') 'Local idiag, jdiag, task:',   &
                                   model%numerics%idiag_local,  &
                                   model%numerics%jdiag_local,  &
                                   model%numerics%rdiag_local
       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif  ! main_task

  end subroutine glide_init_diag

!--------------------------------------------------------------------------

  subroutine glide_write_diag (model,       time)

    ! Write global diagnostics
    ! Also write local diagnostics for a selected grid cell
 
    use parallel

    use glimmer_paramets, only: thk0, len0, vel0, tim0, unphys_val
    use glimmer_physcon, only: scyr, rhoi, shci
 
    implicit none
 
    ! input/output arguments

    ! Note: model is intent(inout) so that some global diagnostics can be computed here 
    type(glide_global_type), intent(inout) :: model ! model instance

    real(dp),  intent(in) :: time                   ! current time in years

    ! local variables

    real(dp) ::                         &
         minthck,                       &    ! ice thickness threshold (m) for global diagnostics
         tot_area,                      &    ! total ice area (m^2)
         tot_area_ground,               &    ! total area of grounded ice (m^2)
         tot_area_float,                &    ! total area of floating ice (m^2)
         area_cell,                     &    ! cell area
         tot_volume,                    &    ! total ice volume (m^3)
         tot_mass,                      &    ! total ice mass (kg)
         tot_mass_above_flotation,      &    ! total ice mass above flotation (kg)
         thck_floating,                 &    ! thickness of floating ice
         thck_above_flotation,          &    ! thickness above flotation
         tot_energy,                    &    ! total ice energy (J)
         tot_smb_flux,                  &    ! total surface mass balance flux (kg/s)
         tot_bmb_flux,                  &    ! total basal mass balance flux (kg/s)
         tot_calving_flux,              &    ! total calving flux (kg/s)
         tot_gl_flux,                   &    ! total grounding line flux (kg/s)
         tot_acab,                      &    ! total surface accumulation/ablation rate (m^3/yr)
         tot_bmlt,                      &    ! total basal melt rate (m^3/yr)
         tot_calving,                   &    ! total calving rate (m^3/yr)
         tot_dmass_dt,                  &    ! rate of change of total mass (kg/s)
         err_dmass_dt,                  &    ! mass conservation error (kg/s)
                                             ! given by dmass_dt - (tot_acab - tot_bmlt - tot_calving)
         mean_thck,                     &    ! mean ice thickness (m)
         mean_temp,                     &    ! mean ice temperature (deg C)
         mean_acab,                     &    ! mean surface accumulation/ablation rate (m/yr)
         mean_bmlt,                     &    ! mean basal melt (m/yr)
         mean_calving,                  &    ! mean calving (m/yr)
         max_thck, max_thck_global,     &    ! max ice thickness (m)
         max_temp, max_temp_global,     &    ! max ice temperature (deg C)
         min_temp, min_temp_global,     &    ! min ice temperature (deg C)
         max_spd_sfc, max_spd_sfc_global,   &    ! max surface ice speed (m/yr)
         max_spd_bas, max_spd_bas_global,   &    ! max basal ice speed (m/yr)
         spd,                           &    ! speed
         thck_diag, usrf_diag,          &    ! local column diagnostics
         topg_diag, relx_diag,          &    
         load_diag,                     &
         artm_diag, acab_diag,          &
         bmlt_diag, bwat_diag,          &
         bheatflx_diag, level,          &
         factor                              ! unit conversion factor

    integer, dimension(model%general%ewn,model%general%nsn) ::  &
         ice_mask,     &! = 1 where ice is present with thck > minthck, else = 0
         floating_mask  ! = 1 where ice is present and floating, else = 0

    real(dp), dimension(model%general%upn) ::  &
         temp_diag,                     &    ! Note: sfc temp not included if temps are staggered
                                             !       (use artm instead)
         spd_diag

    real(dp), dimension(model%lithot%nlayer) ::  &
         lithtemp_diag                       ! lithosphere column diagnostics

    integer :: i, j, k, ktop, kbed,               &
               imax, imin,                        &
               jmax, jmin,                        &
               kmax, kmin,                        &
               imax_global, imin_global,          &
               jmax_global, jmin_global,          &
               kmax_global, kmin_global,          &
               procnum,                           &
               ewn, nsn, upn,                     &    ! model%numerics%ewn, etc.
               nlith,                             &    ! model%lithot%nlayer
               velo_ew_ubound, velo_ns_ubound          ! upper bounds for velocity variables
 
    character(len=100) :: message
    
    real(dp), dimension(:,:), allocatable :: &
         cell_area     ! grid cell areas (scaled model units)
                       ! optionally, divide by scale factor^2 to account for grid distortion

    real(dp), parameter ::   &
       eps = 1.0d-11,         & ! small number
       eps_thck = 1.0d-11       ! threshold thickness (m) for writing diagnostics

    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn

    allocate(cell_area(ewn,nsn))
    cell_area(:,:) = model%numerics%dew * model%numerics%dns

    ! Note: If projection%stere%compute_area_factor = .true., then area factors will differ from 1.
    !       Then the total ice area and volume computed below will be corrected for area distortions,
    !        giving a better estimate of the true ice area and volume.
    !       However, applying scale factors will give a mass conservation error (total dmass_dt > 0)
    !        in the diagnostics, because horizontal transport does not account for area factors.
    !        Transport conserves mass only under the assumption of rectangular grid cells.

    if (associated(model%projection%stere)) then   ! divide cell area by area_factor^2
       do j = 1, nsn
          do i = 1, ewn
             if (model%projection%stere%area_factor(i,j) > 0.0d0) then
                cell_area(i,j) = cell_area(i,j) / model%projection%stere%area_factor(i,j)**2
             endif
          enddo
       enddo
    endif

    nlith = model%lithot%nlayer

    if (uhalo > 0) then
       velo_ns_ubound = nsn-uhalo
       velo_ew_ubound = ewn-uhalo
    else
       ! for uhalo==0 (as is the case for the glide dycore), the velocity grid has one less
       ! point than the main grid, so we need to subtract one to avoid out-of-bounds problems
       velo_ns_ubound = nsn-uhalo-1
       velo_ew_ubound = ewn-uhalo-1
    end if

    ! Set the minimum ice thickness for including cells in diagnostics
    if (model%options%diag_minthck == DIAG_MINTHCK_ZERO) then
       minthck = eps_thck  ! slightly > 0
    elseif (model%options%diag_minthck == DIAG_MINTHCK_THKLIM) then
       minthck = model%numerics%thklim*thk0
    endif

    !-----------------------------------------------------------------
    ! Compute some masks that are useful for diagnostics
    !-----------------------------------------------------------------

    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j)*thk0 > minthck) then
             ice_mask(i,j) = 1
             if (model%geometry%topg(i,j) - model%climate%eus < (-rhoi/rhoo)*model%geometry%thck(i,j)) then
                floating_mask(i,j) = 1
             else
                floating_mask(i,j) = 0
             endif
          else
             ice_mask(i,j) = 0
             floating_mask(i,j) = 0
          endif
       enddo
    enddo

    !-----------------------------------------------------------------
    ! Compute and write global diagnostics
    !-----------------------------------------------------------------
 
    call write_log('----------------------------------------------------------')
    call write_log(' ')
    write(message,'(a25,f24.16)') 'Diagnostic output, time =', time
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    call write_log(' ')

    ! total ice area (m^2)

    tot_area = 0.d0
    tot_area_ground = 0.d0
    tot_area_float = 0.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
             tot_area = tot_area + cell_area(i,j)
             if (floating_mask(i,j) == 1) then
                tot_area_float = tot_area_float + cell_area(i,j)
             else
                tot_area_ground = tot_area_ground + cell_area(i,j)
             endif
          endif
       enddo
    enddo

    tot_area = tot_area * len0**2
    tot_area = parallel_reduce_sum(tot_area)

    tot_area_ground = tot_area_ground * len0**2
    tot_area_ground = parallel_reduce_sum(tot_area_ground)

    tot_area_float = tot_area_float * len0**2
    tot_area_float = parallel_reduce_sum(tot_area_float)

    ! total ice volume (m^3)
 
    tot_volume = 0.d0
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
             tot_volume = tot_volume + model%geometry%thck(i,j) * cell_area(i,j)
          endif
       enddo
    enddo
    tot_volume = tot_volume * thk0 * len0**2
    tot_volume = parallel_reduce_sum(tot_volume)

    ! total ice mass (kg)
    tot_mass = tot_volume * rhoi

    ! total ice mass above flotation (kg)
    tot_mass_above_flotation = 0.d0

    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
             if (floating_mask(i,j) == 0) then  ! grounded ice
                if (model%geometry%topg(i,j) - model%climate%eus < 0.0d0) then  ! grounded below sea level
                   thck_floating = (-rhoo/rhoi) * (model%geometry%topg(i,j) - model%climate%eus)  ! thickness of ice that is exactly floating
                   thck_above_flotation = model%geometry%thck(i,j) - thck_floating
                   tot_mass_above_flotation = tot_mass_above_flotation    &
                                            + thck_above_flotation * cell_area(i,j)
                else   ! grounded above sea level
                   tot_mass_above_flotation = tot_mass_above_flotation    &
                                            + model%geometry%thck(i,j) * cell_area(i,j)
                endif
             endif
          endif
       enddo
    enddo

    tot_mass_above_flotation = tot_mass_above_flotation * thk0 * len0**2  ! convert to m^3
    tot_mass_above_flotation = tot_mass_above_flotation * rhoi            ! convert from m^3 to kg
    tot_mass_above_flotation = parallel_reduce_sum(tot_mass_above_flotation)

    ! total ice energy relative to T = 0 deg C (J)
 
    tot_energy = 0.d0
    if (size(model%temper%temp,1) == upn+1) then  ! temps are staggered in vertical
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             if (ice_mask(i,j) == 1) then
                do k = 1, upn-1
                   tot_energy = tot_energy +   &
                                model%geometry%thck(i,j) * model%temper%temp(k,i,j) * cell_area(i,j)   &
                                *(model%numerics%sigma(k+1) - model%numerics%sigma(k))
                enddo
             endif
          enddo
       enddo
    
    else   ! temps are unstaggered in vertical
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             if (ice_mask(i,j) == 1) then
                ! upper half-layer, T = upper sfc temp
                tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(1,i,j) * cell_area(i,j)    &
                             * 0.5d0 * model%numerics%sigma(2)
                do k = 2, upn-1
                   tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(k,i,j) * cell_area(i,j)  &
                             * 0.5d0*(model%numerics%sigma(k+1) - model%numerics%sigma(k-1))
                enddo
                ! lower half-layer, T = lower sfc temp
                tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(upn,i,j) * cell_area(i,j)  &
                             * 0.5d0 * (1.0d0 - model%numerics%sigma(upn-1))
             endif
          enddo
       enddo
    endif

    tot_energy = tot_energy * thk0 * len0**2 * rhoi * shci
    tot_energy = parallel_reduce_sum(tot_energy)

    ! mean thickness

    if (tot_area > eps) then
       mean_thck = tot_volume/tot_area
    else
       mean_thck = 0.d0
    endif

    ! mean temperature
 
    if (tot_volume > eps) then
       mean_temp = tot_energy/ (rhoi*shci*tot_volume)
    else
       mean_temp = 0.d0
    endif
 
    ! copy some global scalars to the geometry derived type
    ! Note: These have SI units (e.g, m^2 for area, m^3 for volume)

    model%geometry%iarea  = tot_area
    model%geometry%iareag = tot_area_ground
    model%geometry%iareaf = tot_area_float
    model%geometry%ivol   = tot_volume
    model%geometry%imass  = tot_mass
    model%geometry%imass_above_flotation  = tot_mass_above_flotation

    ! For Glissade only, compute a global mass budget and check mass conservation

    if (model%options%whichdycore == DYCORE_GLISSADE) then

       ! total surface accumulation/ablation rate (m^3/yr ice)
 
       tot_acab = 0.d0
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             tot_acab = tot_acab + model%climate%acab_applied(i,j) * cell_area(i,j)
          enddo
       enddo

       tot_acab = tot_acab * scyr * thk0 / tim0 * len0**2   ! convert to m^3/yr
       tot_acab = parallel_reduce_sum(tot_acab)

       ! total surface mass balance flux (kg/s)
       tot_smb_flux = tot_acab * rhoi / scyr   ! convert m^3/yr to kg/s

       ! mean accumulation/ablation rate (m/yr)
       ! Note: This will be only approximate if some ice has melted completely during the time step
       if (tot_area > eps) then
          mean_acab = tot_acab/tot_area    ! divide by total area to get m/yr
       else
          mean_acab = 0.d0
       endif

       ! total basal melting rate (positive for ice loss)
       tot_bmlt = 0.d0
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             tot_bmlt = tot_bmlt + model%basal_melt%bmlt_applied(i,j) * cell_area(i,j)
          enddo
       enddo

       tot_bmlt = tot_bmlt * scyr * thk0 / tim0 * len0**2   ! convert to m^3/yr
       tot_bmlt = parallel_reduce_sum(tot_bmlt)

       ! total basal mass balance (kg/s, positive for freeze-on, negative for melt)
       tot_bmb_flux = -tot_bmlt * rhoi / scyr   ! convert m^3/yr to kg/s

       ! mean basal melt rate (m/yr)
       ! Note: This will be only approximate if some ice has melted completely during the time step
       if (tot_area > eps) then
          mean_bmlt = tot_bmlt/tot_area    ! divide by total area to get m/yr
       else
          mean_bmlt = 0.d0
       endif

       ! total calving rate (m^3/yr ice)
       ! Note: calving%calving_rate has units of m/yr ice

       tot_calving = 0.d0
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             tot_calving = tot_calving + model%calving%calving_rate(i,j) * (cell_area(i,j)*len0**2)  ! m^3/yr ice
          enddo
       enddo
       tot_calving = parallel_reduce_sum(tot_calving)

       ! total calving mass balance flux (kg/s, negative for ice loss by calving)
       tot_calving_flux = -tot_calving * rhoi / scyr   ! convert m^3/yr to kg/s

       ! mean calving rate (m/yr)
       ! Note: This will be only approximate if some ice has melted completely during the time step
       if (tot_area > eps) then
          mean_calving = tot_calving/tot_area    ! divide by total area to get m/yr
       else
          mean_calving = 0.d0
       endif

       ! total grounding line mass balance flux (< 0 by definition)
       ! Note: At this point, gl_flux_east and gl_flux_north are already dimensionalized in kg/m/s,
       !       so tot_gl_flux will have units of kg/s

       tot_gl_flux = 0.d0
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             tot_gl_flux = tot_gl_flux - abs(model%geometry%gl_flux_east(i,j)) * model%numerics%dns*len0 &
                                       - abs(model%geometry%gl_flux_north(i,j)) * model%numerics%dew*len0
          enddo
       enddo
       tot_gl_flux = parallel_reduce_sum(tot_gl_flux)

       ! total rate of change of ice mass (kg/s)
       ! Note: dthck_dt has units of m/s
       tot_dmass_dt = 0.d0
       do j = lhalo+1, nsn-uhalo
          do i = lhalo+1, ewn-uhalo
             tot_dmass_dt = tot_dmass_dt + model%geometry%dthck_dt(i,j) * cell_area(i,j)
          enddo
       enddo
       tot_dmass_dt = tot_dmass_dt * rhoi * len0**2  ! convert to kg/s
       tot_dmass_dt = parallel_reduce_sum(tot_dmass_dt)

       ! mass conservation error
       ! Note: For most runs, this should be close to zero.

       err_dmass_dt = tot_dmass_dt - (tot_smb_flux + tot_bmb_flux + tot_calving_flux)

       ! uncomment to convert total fluxes from kg/s to Gt/yr
!!!    tot_smb_flux = tot_smb_flux * scyr/1.0d12
!!!    tot_bmb_flux = tot_bmb_flux * scyr/1.0d12
!!!    tot_calving_flux = tot_calving_flux * scyr/1.0d12
!!!    tot_gl_flux = tot_gl_flux * scyr/1.0d12
!!!    tot_dmass_dt = tot_dmass_dt * scyr/1.0d12
!!!    err_dmass_dt = err_dmass_dt * scyr/1.0d12

       ! copy some global scalars to the geometry derived type
       ! Note: These have SI units (e.g, m^2 for area, m^3 for volume)
       model%geometry%total_smb_flux = tot_smb_flux
       model%geometry%total_bmb_flux = tot_bmb_flux
       model%geometry%total_calving_flux = tot_calving_flux
       model%geometry%total_gl_flux = tot_gl_flux

    endif  ! Glissade dycore

    ! write global sums and means to log file

    write(message,'(a25,e24.16)') 'Total ice area (km^2)    ',   &
                                   tot_area*1.0d-6           ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Grounded ice area (km^2) ',   &
                                   tot_area_ground*1.0d-6    ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Floating ice area (km^2) ',   &
                                   tot_area_float*1.0d-6     ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Total ice volume (km^3)  ',   &
                                   tot_volume*1.0d-9         ! convert to km^3
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Total ice mass (kg)      ',   &
                                   tot_mass                  ! kg
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Mass above flotation (kg)',   &
                                   tot_mass_above_flotation  ! kg
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,e24.16)') 'Total ice energy (J)     ', tot_energy
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    if (model%options%whichdycore == DYCORE_GLISSADE) then

       if (model%options%dm_dt_diag == DM_DT_DIAG_KG_S) then

          write(message,'(a25,e24.16)') 'Total SMB flux (kg/s)    ', tot_smb_flux
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total BMB flux (kg/s)    ', tot_bmb_flux
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total calving flux (kg/s)', tot_calving_flux
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total dmass/dt (kg/s)    ', tot_dmass_dt
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'dmass/dt error (kg/s)    ', err_dmass_dt
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total gr line flux (kg/s)', tot_gl_flux
          call write_log(trim(message), type = GM_DIAGNOSTIC)

       elseif (model%options%dm_dt_diag == DM_DT_DIAG_GT_Y) then

          factor = scyr / 1.0d12

          write(message,'(a25,e24.16)') 'Total SMB flux (Gt/y)    ', tot_smb_flux * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total BMB flux (Gt/y)    ', tot_bmb_flux * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total calving flux (Gt/y)', tot_calving_flux * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total dmass/dt (Gt/y)    ', tot_dmass_dt * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'dmass/dt error (Gt/y)    ', err_dmass_dt * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          write(message,'(a25,e24.16)') 'Total gr line flux (Gt/y)', tot_gl_flux * factor
          call write_log(trim(message), type = GM_DIAGNOSTIC)

       endif

!       write(message,'(a25,e24.16)') 'Mean accum/ablat (m/yr)  ', mean_acab
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

!       write(message,'(a25,e24.16)') 'Mean basal melt (m/yr)   ', mean_bmlt
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

!       write(message,'(a25,e24.16)') 'Mean calving (m/yr)      ', mean_calving
!       call write_log(trim(message), type = GM_DIAGNOSTIC)

    endif  ! Glissade dycore

    write(message,'(a25,f24.16)') 'Mean thickness (m)       ', mean_thck
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    write(message,'(a25,f24.16)') 'Mean temperature (C)     ', mean_temp
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! Find various global maxes and mins

    ! max thickness

    imax = 0
    jmax = 0
    max_thck = unphys_val   ! = an arbitrary large negative number
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (model%geometry%thck(i,j) > max_thck) then
             max_thck = model%geometry%thck(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo

    imax_global = 0
    jmax_global = 0
    max_thck_global = parallel_reduce_max(max_thck)
    if (max_thck == max_thck_global) then  ! max_thck lives on this processor
       imax_global = (imax - lhalo) + global_col_offset
       jmax_global = (jmax - lhalo) + global_row_offset
    endif
    imax_global = parallel_reduce_max(imax_global)
    jmax_global = parallel_reduce_max(jmax_global)

    write(message,'(a25,f24.16,2i6)') 'Max thickness (m), i, j  ',   &
                                       max_thck_global*thk0, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max temperature

    ktop = lbound(model%temper%temp,1)
    kbed = ubound(model%temper%temp,1)

    imax = 0
    jmax = 0
    kmax = 0
    max_temp =  unphys_val
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
            do k = ktop, kbed
                if (model%temper%temp(k,i,j) > max_temp) then
                   max_temp = model%temper%temp(k,i,j)
                   imax = i
                   jmax = j
                   kmax = k
                endif
             enddo
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_temp, xout=max_temp_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global)
    kmax_global = kmax
    call broadcast(imax_global, procnum)
    call broadcast(jmax_global, procnum)
    call broadcast(kmax_global, procnum)

    write(message,'(a25,f24.16,3i6)') 'Max temperature, i, j, k ',   &
                    max_temp_global, imax_global, jmax_global, kmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! min temperature

    imin = 0
    jmin = 0
    kmin = 0
    min_temp = 999.d0  ! arbitrary large positive number
    do j = lhalo+1, nsn-uhalo
       do i = lhalo+1, ewn-uhalo
          if (ice_mask(i,j) == 1) then
             do k = ktop, kbed
                if (model%temper%temp(k,i,j) < min_temp) then
                   min_temp = model%temper%temp(k,i,j)
                   imin = i
                   jmin = j
                   kmin = k
                endif
             enddo
          endif
       enddo
    enddo

    call parallel_reduce_minloc(xin=min_temp, xout=min_temp_global, xprocout=procnum)
    call parallel_globalindex(imin, jmin, imin_global, jmin_global)
    kmin_global = kmin
    call broadcast(imin_global, procnum)
    call broadcast(jmin_global, procnum)
    call broadcast(kmin_global, procnum)

    write(message,'(a25,f24.16,3i6)') 'Min temperature, i, j, k ',   &
                    min_temp_global, imin_global, jmin_global, kmin_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max surface speed
    imax = 0
    jmax = 0
    max_spd_sfc = unphys_val

    do j = lhalo+1, velo_ns_ubound
       do i = lhalo+1, velo_ew_ubound
          spd = sqrt(model%velocity%uvel(1,i,j)**2   &
                   + model%velocity%vvel(1,i,j)**2)
          if (model%geomderv%stagthck(i,j)*thk0 > minthck .and. spd > max_spd_sfc) then
             max_spd_sfc = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_spd_sfc, xout=max_spd_sfc_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global)
    call broadcast(imax_global, procnum)
    call broadcast(jmax_global, procnum)

    write(message,'(a25,f24.16,2i6)') 'Max sfc spd (m/yr), i, j ',   &
                    max_spd_sfc_global*vel0*scyr, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! max basal speed

    imax = 0
    jmax = 0
    max_spd_bas = unphys_val
    do j = lhalo+1, velo_ns_ubound
       do i = lhalo+1, velo_ew_ubound
          spd = sqrt(model%velocity%uvel(upn,i,j)**2   &
                   + model%velocity%vvel(upn,i,j)**2)
          if (model%geomderv%stagthck(i,j)*thk0 > minthck  .and. spd > max_spd_bas) then
             max_spd_bas = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo

    call parallel_reduce_maxloc(xin=max_spd_bas, xout=max_spd_bas_global, xprocout=procnum)
    call parallel_globalindex(imax, jmax, imax_global, jmax_global)
    call broadcast(imax_global, procnum)
    call broadcast(jmax_global, procnum)

    write(message,'(a25,f24.16,2i6)') 'Max base spd (m/yr), i, j',   &
                    max_spd_bas_global*vel0*scyr, imax_global, jmax_global
    call write_log(trim(message), type = GM_DIAGNOSTIC)

    ! local diagnostics

    ! initialize to unphysical negative values
    usrf_diag     = unphys_val
    thck_diag     = unphys_val
    topg_diag     = unphys_val
    relx_diag     = unphys_val
    load_diag     = unphys_val
    artm_diag     = unphys_val
    acab_diag     = unphys_val
    bmlt_diag     = unphys_val
    bwat_diag     = unphys_val
    bheatflx_diag = unphys_val
    temp_diag(:)  = unphys_val
    spd_diag (:)  = unphys_val
    lithtemp_diag(:) = unphys_val    

    ! Set local diagnostic values, and communicate them to main_task
       
    if (model%numerics%idiag_local >= 1 .and. model%numerics%idiag_local <= ewn  &
                                        .and.                                    &
        model%numerics%jdiag_local >= 1 .and. model%numerics%jdiag_local <= nsn) then

       if (this_rank == model%numerics%rdiag_local) then

          i = model%numerics%idiag_local
          j = model%numerics%jdiag_local
          usrf_diag = model%geometry%usrf(i,j)*thk0
          thck_diag = model%geometry%thck(i,j)*thk0
          topg_diag = model%geometry%topg(i,j)*thk0
          if (model%options%isostasy == ISOSTASY_COMPUTE) then
             relx_diag = model%isostasy%relx(i,j)*thk0
             load_diag = model%isostasy%load(i,j)*thk0
          endif
          artm_diag = model%climate%artm(i,j)
          acab_diag = model%climate%acab_applied(i,j) * thk0*scyr/tim0
          bmlt_diag = model%basal_melt%bmlt_applied(i,j) * thk0*scyr/tim0
          bwat_diag = model%temper%bwat(i,j) * thk0
          bheatflx_diag = model%temper%bheatflx(i,j)
       
          temp_diag(:) = model%temper%temp(1:upn,i,j)          
          spd_diag(:) = sqrt(model%velocity%uvel(1:upn,i,j)**2   &
                           + model%velocity%vvel(1:upn,i,j)**2) * vel0*scyr
          if (model%options%gthf == GTHF_COMPUTE) &
               lithtemp_diag(:) = model%lithot%temp(i,j,:)
       endif

       usrf_diag = parallel_reduce_max(usrf_diag)
       thck_diag = parallel_reduce_max(thck_diag)
       topg_diag = parallel_reduce_max(topg_diag)
       if (model%options%isostasy == ISOSTASY_COMPUTE) then
          relx_diag = parallel_reduce_max(relx_diag)
          load_diag = parallel_reduce_max(load_diag)
       endif
       artm_diag = parallel_reduce_max(artm_diag)
       acab_diag = parallel_reduce_max(acab_diag)
       bmlt_diag = parallel_reduce_max(bmlt_diag)
       bwat_diag = parallel_reduce_max(bwat_diag)
       bheatflx_diag = parallel_reduce_max(bheatflx_diag)

       do k = 1, upn
          temp_diag(k) = parallel_reduce_max(temp_diag(k))
          spd_diag(k)  = parallel_reduce_max(spd_diag(k))
       enddo

       do k = 1, nlith
          lithtemp_diag(k) = parallel_reduce_max(lithtemp_diag(k))
       enddo

       call write_log(' ')
       write(message,'(a39,2i6)')  &
            'Grid point diagnostics: (i,j) =', model%numerics%idiag, &
                                               model%numerics%jdiag
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       write(message,'(a39,3i6)')  &
            'Local (i,j,rank) =             ', model%numerics%idiag_local, &
                                               model%numerics%jdiag_local, &
                                               model%numerics%rdiag_local
       call write_log(trim(message), type = GM_DIAGNOSTIC)
       call write_log(' ')
 
       write(message,'(a25,f24.16)') 'Upper surface (m)        ', usrf_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Thickness (m)            ', thck_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Bedrock topo (m)         ', topg_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       if (model%options%isostasy == ISOSTASY_COMPUTE) then
          write(message,'(a25,f24.16)') 'Relaxed bedrock (m)   ', relx_diag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
          write(message,'(a25,f24.16)') 'Load deflection (m)   ', load_diag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
       endif

       write(message,'(a25,f24.16)') 'Sfc mass balance (m/yr)  ', acab_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal mass balance (m/yr)', -bmlt_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal water depth (m)    ', bwat_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)

       write(message,'(a25,f24.16)') 'Basal heat flux (W/m^2)  ', bheatflx_diag
       call write_log(trim(message), type = GM_DIAGNOSTIC)
 
       ! Vertical profile of ice speed and temperature

       call write_log(' ')
       write(message,'(a55)') ' Sigma       Ice speed (m/yr)       Ice temperature (C)'
       call write_log(trim(message), type = GM_DIAGNOSTIC)
 
       if (size(model%temper%temp,1) == upn+1) then   ! temperatures staggered in vertical
                                                      ! (at layer midpoints)

           ! upper surface 
           write (message,'(f6.3,2f24.16)') model%numerics%sigma(1), spd_diag(1), artm_diag
           call write_log(trim(message), type = GM_DIAGNOSTIC)

           ! internal
           do k = 1, upn-1

              ! speed at top of layer
              if (k > 1) then
                 write (message,'(f6.3,f24.16)') model%numerics%sigma(k), spd_diag(k)
                 call write_log(trim(message), type = GM_DIAGNOSTIC)
              endif

              ! temp at layer midpoint
              write (message,'(f6.3,24x,f24.16)') model%numerics%stagsigma(k), temp_diag(k)
              call write_log(trim(message), type = GM_DIAGNOSTIC)

           enddo

           ! lower surface
           write (message,'(f6.3,2f24.16)') model%numerics%sigma(upn), spd_diag(upn), temp_diag(upn)
           call write_log(trim(message), type = GM_DIAGNOSTIC)
           
       else    ! temperatures unstaggered in vertical (at layer interfaces)
 
           do k = 1, upn
             write (message,'(f6.3,2f24.16)') model%numerics%sigma(k), spd_diag(k), temp_diag(k)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

       endif  ! temps staggered

       ! Vertical profile of upper lithosphere temperature

       if (model%options%gthf == GTHF_COMPUTE) then

          call write_log(' ')
          write(message,'(a41)') '  Level (m)          Lithosphere temp (C)'
          call write_log(trim(message), type = GM_DIAGNOSTIC)

          level = 0.d0
          do k = 1, nlith
             level = level + model%lithot%deltaz(nlith)
             write (message,'(f10.0,6x,f24.16)') level, lithtemp_diag(k)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

       endif  ! gthf_compute

    endif     ! idiag_local and jdiag_local in bounds

    call write_log(' ')

  end subroutine glide_write_diag
     
!==============================================================

end module glide_diagnostics
