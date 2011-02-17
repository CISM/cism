! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                             +
! +  glide_diagnostics.f90 - part of the Glimmer-CISM ice model + 
! +                                                             +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Copyright (C) 2010
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
 
module glide_diagnostics
  !*FD subroutines for computing various useful diagnostics
  ! Author: William Lipscomb, LANL 
 
contains
  subroutine glide_write_diag (model, time, idiag, jdiag)
    !*FD Writes diagnostic output
 
    use glimmer_log
    use glimmer_paramets, only: thk0, len0, vel0, tim0
    use glimmer_physcon, only: scyr, rhoi, shci
    use glide_types
 
    implicit none
 
    type(glide_global_type), intent(in) :: model    ! model instance
    real(dp),  intent(in)   :: time                 ! current time in years
    integer, intent(in), optional :: idiag, jdiag
 
    real(dp) ::          &
         tot_area,       &    ! total ice area (km^2)
         tot_volume,     &    ! total ice volume (km^3)
         tot_energy,     &    ! total ice energy (J)
         max_thck,       &    ! max ice thickness (m)
         min_temp,       &    ! min ice temperature (deg C)
         mean_thck,      &    ! mean ice thickness (m)
         mean_temp,      &    ! mean ice temperature (deg C)
         max_spd_sfc,    &    ! max surface ice speed (m/yr)
         max_spd_bas,    &    ! max basal ice speed (m/yr)
         thck,           &    ! thickness
         spd                  ! speed

    integer :: i, j, k,            &
               imin, jmin, kmin,   &
               imax, jmax, kmax,   &
               ewn, nsn, upn       ! model%numerics%ewn, etc.
 
    character(len=100) :: message
 
    real(dp), parameter ::   &
         eps = 1.0e-11_dp     ! small number
 
    ewn = model%general%ewn
    nsn = model%general%nsn
    upn = model%general%upn
 
    ! If running HO 
    !-----------------------------------------------------------------
    ! Compute and write global diagnostics
    !-----------------------------------------------------------------
 
    print*, ' '
    print*, 'Writing diagnostics to log file, time(yr) =', time
 
    call write_log(' ')
    write(message,'(a32,f10.2)') 'Global diagnostic output, time =', time
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    call write_log(' ')
 
    ! total ice area
 
    tot_area = 0.d0
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             tot_area = tot_area + 1.0d0
          endif
       enddo
    enddo
    tot_area = tot_area * model%numerics%dew * model%numerics%dns * len0**2
    write(message,'(a25,e24.16)') 'Total ice area (km^2)     ',   &
                                   tot_area*1.0d-6   ! convert to km^2
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! total ice volume
 
    tot_volume = 0.d0
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             tot_volume = tot_volume + model%geometry%thck(i,j)
          endif
       enddo
    enddo
    tot_volume = tot_volume * model%numerics%dew * model%numerics%dns  &
               * thk0 * len0**2
    write(message,'(a25,e24.16)') 'Total ice volume (km^3)  ',   &
                                   tot_volume*1.0d-9   ! convert to km^3
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! total ice energy relative to T = 0 deg C
 
    tot_energy = 0.d0

    if (size(model%temper%temp,1) == upn+1) then  ! staggered temps
       do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             do k = 1, upn-1
                tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(k,i,j) *   &
                             (model%numerics%sigma(k+1) - model%numerics%sigma(k))
             enddo
          endif
       enddo
    enddo
    tot_energy = tot_energy * model%numerics%dew * model%numerics%dns  &
               * thk0 * len0**2
    write(message,'(a25,e24.16)') 'Total ice energy (J)     ', tot_energy
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    
    ! total sum of volume * age
 
    ! tot_age = 0._dp ! JCC - Removed in parallel
    else   ! unstaggered temps
       do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > model%numerics%thklim) then
             ! upper half-layer, T = upper sfc temp
             tot_energy = tot_energy +   &
                          model%geometry%thck(i,j) * model%temper%temp(1,i,j) *   &
                          0.5d0 * model%numerics%sigma(2)
             do k = 2, upn-1
                tot_energy = tot_energy +   &
                             model%geometry%thck(i,j) * model%temper%temp(k,i,j) *   &
                             0.5d0*(model%numerics%sigma(k+1) - model%numerics%sigma(k-1))
             enddo
             ! lower half-layer, T = lower sfc temp
             tot_energy = tot_energy +   &
                          model%geometry%thck(i,j) * model%temper%temp(upn,i,j) *   &
                          0.5d0 * (1.0d0 - model%numerics%sigma(upn-1))
          endif
       enddo
       enddo
    endif

    tot_energy = tot_energy * model%numerics%dew * model%numerics%dns  &
               * thk0 * len0**2 * rhoi * shci
    write(message,'(a25,e24.16)') 'Total ice energy (J)     ', tot_energy
    call write_log(trim(message), type = GM_DIAGNOSTIC)
    
    ! mean thickness
 
    if (tot_area > eps) then
       mean_thck = tot_volume/tot_area
    else
       mean_thck = 0.d0
    endif
    write(message,'(a25,f24.16)') 'Mean thickness (m)       ', mean_thck
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! max thickness
 
    imax = 0
    jmax = 0
    max_thck = 0.d0
    do j = 1, nsn
       do i = 1, ewn
          if (model%geometry%thck(i,j) > max_thck) then
             max_thck = model%geometry%thck(i,j)
             imax = i
             jmax = j
          endif
       enddo
    enddo
    write(message,'(a25,f24.16,2i4)') 'Max thickness (m), i, j  ',   &
                                       max_thck*thk0, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! mean temperature
 
    if (tot_volume > eps) then
       mean_temp = tot_energy/ (rhoi*shci*tot_volume)
    else
       mean_temp = 0.d0
    endif
    write(message,'(a25,f24.16)') 'Mean temperature (C)     ', mean_temp
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! min temperature
 
    min_temp =  9999.d0
    do j = 1, nsn
       do i = 1, ewn
          do k = 1, upn-1
             if (model%temper%temp(k,i,j) < min_temp) then
                min_temp = model%temper%temp(k,i,j)
                imin = i
                jmin = j
                kmin = k
             endif
          enddo
       enddo
    enddo
    write(message,'(a25,f24.16,3i4)') 'Min temperature, i, j, k ',   &
                                       min_temp, imin, jmin, kmin
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! mean age
 
!     if (tot_age > eps) then ! JCC - Removed in parallel
!        mean_age = tot_age/tot_volume
!     else
!        mean_age = 0._dp
!     endif
!     write(message,'(a25,f24.16)') 'Mean ice age (yr)        ', mean_age
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! max surface speed
    ! Note: uvel and vvel not defined at i = ewn, j = nsn
    imax = 0
    jmax = 0
    max_spd_sfc = 0.d0
    do j = 1, nsn-1
       do i = 1, ewn-1

          !Note: HO velocities are stored in a separate type, velocity_hom
          if (model%options%which_ho_diagnostic  /= 0) then
             spd = sqrt(model%velocity_hom%uvel(1,i,j)**2   &
                      + model%velocity_hom%vvel(1,i,j)**2)
          else
             spd = sqrt(model%velocity%uvel(1,i,j)**2   &
                      + model%velocity%vvel(1,i,j)**2)
          endif
          if (model%geometry%thck(i,j) > eps .and. spd > max_spd_sfc) then
             max_spd_sfc = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo
    write(message,'(a25,f24.16,2i4)') 'Max sfc spd (m/yr), i, j ',   &
                                       max_spd_sfc*vel0*scyr, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
    ! max basal speed
 
    imax = 0
    jmax = 0
    max_spd_bas = 0.d0
    do j = 1, nsn-1
       do i = 1, ewn-1
          if (model%options%which_ho_diagnostic  /= 0) then
             spd = sqrt(model%velocity_hom%uvel(upn,i,j)**2   &
                      + model%velocity_hom%vvel(upn,i,j)**2)
          else
             spd = sqrt(model%velocity%uvel(upn,i,j)**2   &
                      + model%velocity%vvel(upn,i,j)**2)
          endif
          if (model%geometry%thck(i,j) > eps .and. spd > max_spd_bas) then
             max_spd_bas = spd
             imax = i
             jmax = j
          endif
       enddo
    enddo
    write(message,'(a25,f24.16,2i4)') 'Max base spd (m/yr), i, j',   &
                                       max_spd_bas*vel0*scyr, imax, jmax
    call write_log(trim(message), type = GM_DIAGNOSTIC)
 
 
    ! local diagnostics
 
    if (present(idiag) .and. present(jdiag)) then
       if (idiag >= 1 .and. idiag <=model%general%ewn  &
                      .and.                            &
           jdiag >= 1 .and. jdiag <=model%general%nsn) then
          call write_log(' ')
          write(message,'(a30,2i4)')  &
               'Grid point diagnostics: i, j =', idiag, jdiag
          call write_log(trim(message), type = GM_DIAGNOSTIC)
          call write_log(' ')
 
          i = idiag
          j = jdiag
          write(message,'(a25,f24.16)')   &
               'Thickness (m)            ',model%geometry%thck(i,j)*thk0
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          write(message,'(a25,f24.16)')   &
               'Sfc air temperature (C)  ',model%climate%artm(i,j)
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          write(message,'(a25,f24.16)')   &
               'Basal temperature (C)    ',model%temper%temp(upn,i,j)
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          spd = sqrt(model%velocity%ubas(i,j)**2   &
                   + model%velocity%vbas(i,j)**2) * vel0*scyr
          write(message,'(a25,f24.16)')   &
               'Basal speed (m/yr)       ', spd
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          call write_log(' ')
          write(message,'(a50)')   &
               'Level     Speed (m/yr)             Temperature (C)'
          call write_log(trim(message), type = GM_DIAGNOSTIC)
 
          do k = 1, upn

             if (model%options%which_ho_diagnostic  /= 0) then
                spd = sqrt(model%velocity_hom%uvel(k,i,j)**2   &
                         + model%velocity_hom%vvel(k,i,j)**2)
             else
                spd = sqrt(model%velocity%uvel(k,i,j)**2   &
                         + model%velocity%vvel(k,i,j)**2)
             endif
             spd = spd * vel0*scyr
             write (message,'(i4,2f24.16)')  &
                k, spd, model%temper%temp(k,i,j)
             call write_log(trim(message), type = GM_DIAGNOSTIC)
          enddo

#ifdef GLC_DEBUG
             do k = 1, upn-1
                spd = sqrt(model%velocity%uvel(k,i,j)**2   &
                         + model%velocity%vvel(k,i,j)**2) * vel0*scyr

                write (message,'(i4,Z24.20,Z24.20)')  &
                   k, spd, model%temper%temp(k,i,j)
                call write_log(trim(message), type = GM_DIAGNOSTIC)
             enddo
#endif
 
       endif  ! idiag and jdiag in bounds
    endif     ! idiag and jdiag present
 
  end subroutine glide_write_diag
     
!==============================================================

end module glide_diagnostics
