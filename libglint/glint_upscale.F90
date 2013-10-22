!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_upscale.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

  module glint_upscale

  ! This module contains subroutines for upscaling fields from the local to the global grid.
  ! Much of the actual work is done at a lower level, in glint_interp.F90.

  use glint_type
  use glint_constants
  use glimmer_global, only: dp
  implicit none

  private
  public glint_upscaling, glint_upscaling_gcm,  &
         glint_accumulate_output_gcm

contains

  subroutine glint_upscaling(instance,                    &
                             orog,         albedo,        &
                             ice_frac,     veg_frac,      &
                             snowice_frac, snowveg_frac,  &
                             snow_depth)

    !*FD Upscales and returns certain fields
    !*FD Output fields are only valid on the main task
    !*FD 
    !*FD \begin{itemize}
    !*FD \item \texttt{orog} --- the orographic elevation (m)
    !*FD \item \texttt{albedo} --- the albedo of ice/snow (this is only a notional value --- need to do
    !*FD some work here)
    !*FD \item \texttt{ice\_frac} --- The fraction covered by ice
    !*FD \item \texttt{veg\_frac} --- The fraction of exposed vegetation
    !*FD \item \texttt{snowice\_frac} --- The fraction of snow-covered ice
    !*FD \item \texttt{snowveg\_frac} --- The fraction of snow-covered vegetation
    !*FD \item \texttt{snow_depth} --- The mean snow-depth over those parts covered in snow (m w.e.)
    !*FD \end{itemize}

    use glimmer_paramets
    use glimmer_coordinates, only: coordsystem_allocate

    ! Arguments ----------------------------------------------------------------------------------------

    type(glint_instance),   intent(in)  :: instance      !*FD the model instance

    real(dp),dimension(:,:),intent(out) :: orog          !*FD the orographic elevation (m)
    real(dp),dimension(:,:),intent(out) :: albedo        !*FD the albedo of ice/snow
    real(dp),dimension(:,:),intent(out) :: ice_frac      !*FD The fraction covered by ice
    real(dp),dimension(:,:),intent(out) :: veg_frac      !*FD The fraction of exposed vegetation
    real(dp),dimension(:,:),intent(out) :: snowice_frac  !*FD The fraction of snow-covered ice
    real(dp),dimension(:,:),intent(out) :: snowveg_frac  !*FD The fraction of snow-covered vegetation
    real(dp),dimension(:,:),intent(out) :: snow_depth    !*FD The mean snow-depth over those 
    !*FD parts covered in snow (m w.e.)

    ! Internal variables -------------------------------------------------------------------------------

    real(dp),dimension(:,:),pointer :: temp => null()

    ! --------------------------------------------------------------------------------------------------
    ! Orography

    call local_to_global_avg(instance%ups_orog, &
                             instance%model%geometry%usrf, &
                             orog,    &
                             instance%out_mask)
    orog=thk0*orog

    call coordsystem_allocate(instance%lgrid,temp)

    ! Ice-no-snow fraction
    where (instance%mbal_accum%snowd == 0.d0 .and. instance%model%geometry%thck > 0.d0)
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere

    call local_to_global_avg(instance%ups, &
                             temp, &
                             ice_frac,    &
                             instance%out_mask)

    ! Ice-with-snow fraction
    where (instance%mbal_accum%snowd > 0.d0 .and. instance%model%geometry%thck > 0.d0)
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere
    call local_to_global_avg(instance%ups, &
                             temp, &
                             snowice_frac,    &
                             instance%out_mask)

    ! Veg-with-snow fraction (if ice <10m thick)
    where (instance%mbal_accum%snowd > 0.d0 .and. instance%model%geometry%thck <= (10.d0/thk0))
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere
    call local_to_global_avg(instance%ups, &
                             temp, &
                             snowveg_frac,    &
                             instance%out_mask)

    ! Remainder is veg only
    veg_frac = 1.d0 - ice_frac - snowice_frac - snowveg_frac

    ! Snow depth

    call local_to_global_avg(instance%ups, &
                             instance%mbal_accum%snowd, &
                             snow_depth,    &
                             instance%out_mask)

    ! Albedo

    where ((ice_frac+snowice_frac) > 0.d0)
       albedo = instance%ice_albedo
    elsewhere
       albedo = 0.d0
    endwhere

    deallocate(temp)
    temp => null()

  end subroutine glint_upscaling

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_upscaling_gcm(instance,    nec,      &
                                 nxl,         nyl,      &
                                 nxg,         nyg,      &
                                 box_areas,             &
                                 gfrac,       gtopo,    &
                                 grofi,       grofl,    &
                                 ghflx,                 &
                                 init_call)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! Output fields are only valid on the main task.
    ! The upscaled fields are passed to the GCM land surface model, which has the option
    !  of updating the fractional area and surface elevation of glaciated gridcells.

    use glimmer_paramets, only: thk0, GLC_DEBUG
    use glimmer_log
    use parallel, only: tasks, main_task

!WHL - debug
    use glimmer_paramets, only: tim0

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance), intent(inout) :: instance      ! the model instance
    integer,              intent(in)    :: nec           ! number of elevation classes
    integer,              intent(in)    :: nxl,nyl       ! local grid dimensions 
    integer,              intent(in)    :: nxg,nyg       ! global grid dimensions 

    real(dp),dimension(nxg,nyg),    intent(in)  :: box_areas ! global grid cell areas (m^2)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gfrac   ! ice-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: ghflx   ! heat flux (m)
    real(dp),dimension(nxg,nyg),    intent(out) :: grofi   ! ice runoff (calving) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg),    intent(out) :: grofl   ! liquid runoff (basal melt) flux (kg/m^2/s)
 
    logical, intent(in), optional :: init_call   ! true if called during initialization

    ! Internal variables ----------------------------------------------------------------------
 
    !TODO - Put this parameter elsewhere? 
    real(dp), parameter :: min_thck = 0.d0    ! min thickness (m) for setting gfrac = 1

    integer :: i, j, n      ! indices
    integer :: il, jl, ig, jg

    character(len=100) :: message

    real(dp) :: dew, dns    ! gridcell dimensions

    real(dp) ::   &
       usrf,               &! surface elevation (m)
       thck                 ! ice thickness (m)

    real(dp), dimension(nxl,nyl) ::  &
       area_l               ! local gridcell area

    real(dp), dimension(nxg,nyg) ::  &
       area_g               ! global gridcell area

    real(dp), dimension(nxl,nyl,nec) ::  &
       area_frac_l,        &! area*frac per elevation class on local grid
       area_topo_l,        &! area*topo per elevation class on local grid
       area_hflx_l          ! area*hflx per elevation class on local grid

    real(dp), dimension(nxl,nyl) ::  &
       area_rofi_l,        &! area*rofi on local grid
       area_rofl_l          ! area*rofl on local grid

    real(dp), dimension(nxg,nyg,nec) ::  &
       area_frac_g,        &! area*frac per elevation class on global grid
       area_topo_g,        &! area*topo per elevation class on global grid
       area_hflx_g          ! area*hflx per elevation class on global grid

    real(dp), dimension(nxg,nyg) ::  &
       area_rofi_g,        &! area*rofi on global grid
       area_rofl_g          ! area*rofl on global grid

    !TODO - Pass in topomax as an argument instead of hardwiring it here
    real(dp), dimension(0:nec) :: topomax   ! upper elevation limit of each class

    logical :: first_call   ! if calling the first time, then do not average the accumulated fluxes
                            ! use values from restart file if available

    first_call = .false.
    if (present(init_call)) then
       if (init_call) first_call = .true.
    endif 

    dew = get_dew(instance%model)
    dns = get_dns(instance%model)

    ! Given the value of nec, specify the upper and lower elevation boundaries of each class.
    ! TODO: These must be consistent with the values in the GCM.  Better to pass as an argument.

    if (nec == 1) then
       topomax = (/ 0._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 3) then
       topomax = (/ 0._dp,  1000._dp,  2000._dp, 10000._dp, 10000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 5) then
       topomax = (/ 0._dp,   500._dp,  1000._dp,  1500._dp,  2000._dp, 10000._dp,  &
                           10000._dp, 10000._dp, 10000._dp, 10000._dp, 10000._dp /)
    elseif (nec == 10) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   700._dp,  1000._dp,  1300._dp,  &
                            1600._dp,  2000._dp,  2500._dp,  3000._dp, 10000._dp /)
    elseif (nec == 36) then
       topomax = (/ 0._dp,   200._dp,   400._dp,   600._dp,   800._dp,  &
                 1000._dp,  1200._dp,  1400._dp,  1600._dp,  1800._dp,  &
                 2000._dp,  2200._dp,  2400._dp,  2600._dp,  2800._dp,  &
                 3000._dp,  3200._dp,  3400._dp,  3600._dp,  3800._dp,  &
                 4000._dp,  4200._dp,  4400._dp,  4600._dp,  4800._dp,  &
                 5000._dp,  5200._dp,  5400._dp,  5600._dp,  5800._dp,  &
                 6000._dp,  6200._dp,  6400._dp,  6600._dp,  6800._dp,  &
                 7000._dp, 10000._dp /)
    else
       if (GLC_DEBUG .and. main_task) then
          write(message,'(a6,i3)') 'nec =', nec
          call write_log(trim(message), GM_DIAGNOSTIC)
       end if
       call write_log('ERROR: Current supported values of nec (no. of elevation classes) are 1, 3, 5, 10, or 36', &
                       GM_FATAL,__FILE__,__LINE__)
    endif

    ! The following output only works correctly if running with a single task
    ig = iglint_global    ! defined in glint_type
    jg = jglint_global
    il = instance%model%numerics%idiag
    jl = instance%model%numerics%jdiag
    if (GLC_DEBUG .and. tasks==1) then
       write(stdout,*) 'In glint_upscaling_gcm'
       write(stdout,*) 'il, jl =', il, jl
       write(stdout,*) 'ig, jg =', ig, jg
       write(stdout,*) 'nxl, nyl =', nxl,nyl
       write(stdout,*) 'nxg, nyg =', nxg,nyg
       write(stdout,*) 'av_count_output =', instance%av_count_output
       write(stdout,*) 'local out_mask =', instance%out_mask(il,jl)
       write(stdout,*) 'max, min out_mask =', maxval(instance%out_mask), minval(instance%out_mask)       
    end if

    ! Initialize some fields

    gfrac(:,:,:) = 0.d0
    gtopo(:,:,:) = 0.d0
    ghflx(:,:,:) = 0.d0
    grofi(:,:)   = 0.d0
    grofl(:,:)   = 0.d0

    area_l(:,:) = dew*dns   ! if all grid cells are identical rectangles
    area_g(:,:) = 0.d0

    area_frac_l(:,:,:) = 0.d0
    area_frac_g(:,:,:) = 0.d0

    area_topo_l(:,:,:) = 0.d0
    area_topo_g(:,:,:) = 0.d0

    area_hflx_l(:,:,:) = 0.d0
    area_hflx_g(:,:,:) = 0.d0

    area_rofi_l(:,:) = 0.d0
    area_rofi_g(:,:) = 0.d0

    area_rofl_l(:,:) = 0.d0
    area_rofl_g(:,:) = 0.d0

    ! Compute time-average fluxes (unless called during initialization)

    if (first_call) then
       ! do nothing; use values from restart file if restarting
    else
       if (instance%av_count_output > 0) then
          instance%rofi_tavg(:,:) = instance%rofi_tavg(:,:) / real(instance%av_count_output,dp)
          instance%rofl_tavg(:,:) = instance%rofl_tavg(:,:) / real(instance%av_count_output,dp)
          instance%hflx_tavg(:,:) = instance%hflx_tavg(:,:) / real(instance%av_count_output,dp)
       else
          instance%rofi_tavg(:,:) = 0.d0
          instance%rofl_tavg(:,:) = 0.d0
          instance%hflx_tavg(:,:) = 0.d0
       endif
    endif

    ! Reset the logical variable for averaging output

    instance%new_tavg_output = .true.
    
    ! Loop over local grid cells
    ! Note: The area calculation is not strictly needed on a uniform grid, but is included
    !       to support the general case of spatially varying grid dimensions.

    do j = 1, nyl
    do i = 1, nxl

       usrf = thk0 * instance%model%geometry%usrf(i,j)
       thck = thk0 * instance%model%geometry%thck(i,j)

!WHL - debug
!       if (i==il .and. j==jl) then
!          write(stdout,*) 'usrf =', usrf
!          write(stdout,*) 'thck =', thck
!          write(stdout,*) 'hflx =', instance%hflx_tavg(i,j)
!          write(stdout,*) 'rofi =', instance%rofi_tavg(i,j)
!          write(stdout,*) 'rofl =', instance%rofl_tavg(i,j)
!       endif

       if (thck > min_thck) then   ! this cell is ice-covered

          do n = 1, nec
             if (usrf >= topomax(n-1) .and. usrf < topomax(n)) then  ! local cell is in elev class n
                area_frac_l(i,j,n) = dew*dns           ! Assume fractional coverage = 1.0
                area_topo_l(i,j,n) = dew*dns * usrf
                area_hflx_l(i,j,n) = dew*dns * instance%hflx_tavg(i,j)

                !WHL - debug
!                if (i==il .and. j==jl) then
!                   print*, ' '
!                   print*, 'n =', n
!                   nloc = n
!                   print*, 'dew*dns =', dew*dns
!                   print*, 'dew*dns*usrf =', area_topo_l(i,j,n)
!                   print*, 'dew*dns*hflx =', area_hflx_l(i,j,n)
!                endif

                exit
             endif
          enddo   ! nec
       endif      ! thck > min_thck

       ! Runoff fluxes, without nec loop
       ! Note: These fluxes can be nonzero for cells that no longer have ice.

       area_rofi_l(i,j) = dew*dns * instance%rofi_tavg(i,j)
       area_rofl_l(i,j) = dew*dns * instance%rofl_tavg(i,j)

       !WHL - debug
!       if (i==il .and. j==jl) then
!          print*, ' '
!          print*, 'dew*dns*rofi =', area_rofi_l(i,j)
!          print*, 'dew*dns*rofl =', area_rofl_l(i,j)
!       endif

    enddo         ! i
    enddo         ! j

    ! Map the area-weighted local values to the global grid.

    ! Gridcell area

    call local_to_global_sum(instance%ups, &
                             area_l,       &
                             area_g,       &
                             instance%out_mask)

!WHL - debug
!    print*, ' '
!    print*, 'il, jl, area_l:', il, jl, area_l(il,jl)
!    print*, 'ig, jg, area_g:', ig, jg, area_g(ig,jg)
!    print*, 'ig, jg, box_area:', ig, jg, box_areas(ig,jg)
!    print*, 'area_g/area_l:', area_g(ig,jg)/area_l(il,jl)

    ! Loop over elevation classes for gfrac, gtopo, ghflx

    do n = 1, nec

!WHL - debug
!       print*, ' '
!       print*, 'n =', n

       ! area-weighted gfrac

       call local_to_global_sum(instance%ups,         &
                                area_frac_l(:,:,n),   &
                                area_frac_g(:,:,n),   &
                                instance%out_mask)

!WHL - debug
!       if (n==nloc) then
!          print*, ' '
!          print*, 'il, jl, area_frac_l(il,jl,n):', il, jl, area_frac_l(il,jl,n)
!          print*, 'ig, jg, area_frac_g(ig,jg,n):', ig, jg, area_frac_g(ig,jg,n)
!          print*, 'area_frac_g/area_frac_l:', area_frac_g(ig,jg,n)/area_frac_l(il,jl,n)
!       endif

       ! area-weighted gtopo

       call local_to_global_sum(instance%ups,         &
                                area_topo_l(:,:,n),   &
                                area_topo_g(:,:,n),   &
                                instance%out_mask)

!WHL - debug
!       if (n==nloc) then
!          print*, ' '
!          print*, 'il, jl, area_topo_l(il,jl,n):', il, jl, area_topo_l(il,jl,n)
!          print*, 'ig, jg, area_topo_g(ig,jg,n):', ig, jg, area_topo_g(ig,jg,n)
!          print*, 'area_topo_g/area_topo_l:', area_topo_g(ig,jg,n)/area_topo_l(il,jl,n)
!       endif

       ! area-weighted ghflx

       call local_to_global_sum(instance%ups,         &
                                area_hflx_l(:,:,n),   &
                                area_hflx_g(:,:,n),   &
                                instance%out_mask)

       ! Compute the mean ice fraction, surface elevation, and surface heat flux in each global grid cell

       do j = 1, nyg
       do i = 1, nxg

          ! Find fraction of gridcell area covered by ice in this elevation class
          if (area_g(i,j) > 0.d0) then      
             gfrac(i,j,n) = area_frac_g(i,j,n) / area_g(i,j)
          endif

          ! Find mean surface elevation
          if (area_frac_g(i,j,n) > 0.d0) then
             gtopo(i,j,n) = area_topo_g(i,j,n) / area_frac_g(i,j,n)
          endif

          ! Find mean heat flux
          if (area_frac_g(i,j,n) > 0.d0) then
             ghflx(i,j,n) = area_hflx_g(i,j,n) / area_frac_g(i,j,n)
          endif

!WHL - debug
!          if (i==ig .and. j==jg) then
!             print*, ' '
!             print*, 'ig, jg, n =', i, j, n
!             print*, 'area_g =', area_g(i,j)
!             print*, 'area_frac_g =', area_frac_g(i,j,n)
!             print*, 'area_topo_g =', area_topo_g(i,j,n)
!             print*, 'area_hflx_g =', area_hflx_g(i,j,n)
!             print*, 'gfrac =', gfrac(i,j,n)
!             print*, 'gtopo =', gtopo(i,j,n)
!             print*, 'ghflx =', ghflx(i,j,n)
!          endif

       enddo  ! i
       enddo  ! j

    enddo     ! nec

    ! area-weighted grofi

    call local_to_global_sum(instance%ups,       &
                             area_rofi_l(:,:),   &
                             area_rofi_g(:,:),   &
                             instance%out_mask)

    ! area-weighted grofl

    call local_to_global_sum(instance%ups,       &
                             area_rofl_l(:,:),   &
                             area_rofl_g(:,:),   &
                             instance%out_mask)

    ! Compute the mean ice and liquid runoff in each global grid cell

    do j = 1, nyg
    do i = 1, nxg

       ! Find mean ice runoff from calving 
       ! Note: Here we divide by box_areas (the area of the global grid cell), which in general is not equal
       !       to area_g (the sum over area of the local grid cells associated with the global cell).
       !       We do this to ensure conservation of ice mass (=flux*area) when multiplying later by the
       !       area of the global grid cell.

       if (area_g(i,j) > 0.d0) then
          grofi(i,j) = area_rofi_g(i,j) / box_areas(i,j)
       endif

       ! Find mean liquid runoff from internal/basal melting 
       if (area_g(i,j) > 0.d0) then      
          grofl(i,j) = area_rofl_g(i,j) / box_areas(i,j)
       endif

!WHL - debug
!       if (i==ig .and. j==jg) then
!          print*, ' '
!          print*, 'i, j =', i, j
!          print*, 'area_g =', area_g(i,j)
!          print*, 'box_areas =', box_areas(i,j)
!          print*, 'area_rofi_g =', area_rofi_g(i,j)
!          print*, 'area_rofl_g =', area_rofl_g(i,j)
!          print*, 'grofi(kg/m2/s) =', grofi(i,j)
!          print*, 'grofi(m ice/yr) =', grofi(i,j)/rhoi*scyr
!          print*, 'grofl(kg/m2/s) =', grofl(i,j)
!          print*, 'grofl(m ice/yr) =', grofl(i,j)/rhoi*scyr
!       endif

    enddo   ! i
    enddo   ! j

  end subroutine glint_upscaling_gcm

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_accumulate_output_gcm(model,            &
                                         av_count_output,  &
                                         new_tavg_output,  &
                                         rofi_tavg,        &
                                         rofl_tavg,        &
                                         hflx_tavg)

    ! Given the calving, basal melting, and conductive heat flux fields from the dycore,
    ! accumulate contributions to the rofi, rofl, and hflx fields to be sent to the coupler.

    use glimmer_paramets, only: thk0, tim0

    use glimmer_scales, only: scale_acab  ! for testing

!WHL - debug - Set to inout if specifying the model fields for testing
    type(glide_global_type), intent(in)  :: model
!    type(glide_global_type), intent(inout)  :: model

    integer,  intent(inout) :: av_count_output     ! step counter 
    logical,  intent(inout) :: new_tavg_output     ! if true, start new averaging
    real(dp), dimension(:,:), intent(inout) :: rofi_tavg    ! solid ice runoff (kg m-2 s-1)
    real(dp), dimension(:,:), intent(inout) :: rofl_tavg    ! liquid runoff from basal/interior melting (kg m-2 s-1)
    real(dp), dimension(:,:), intent(inout) :: hflx_tavg    ! conductive heat flux at top surface (W m-2)

!WHL - debug - Uncomment if specifying the model fields for testing

    ! calving field for testing (1 m, scaled to model units)
!    model%climate%calving(:,:) = 1.0d0 / thk0
          
    ! bmlt field for testing (1 m/yr, scaled to model units)
!    model%temper%bmlt(:,:) = 1.0d0 / scale_acab

    ! things to do the first time

    if (new_tavg_output) then

       new_tavg_output = .false.
       av_count_output  = 0

       ! Initialise
       rofi_tavg(:,:) = 0.d0
       rofl_tavg(:,:) = 0.d0
       hflx_tavg(:,:) = 0.d0

    end if

    av_count_output = av_count_output + 1

    !--------------------------------------------------------------------
    ! Accumulate solid runoff (calving)
    !--------------------------------------------------------------------
                       
    ! Note on units: model%climate%calving has dimensionless ice thickness units
    !                Multiply by thk0 to convert to meters of ice
    !                Multiply by rhoi to convert to kg/m^2 water equiv.
    !                Divide by (dt*tim0) to convert to kg/m^2/s

    ! Convert to kg/m^2/s
    rofi_tavg(:,:) = rofi_tavg(:,:)  &
                   + model%climate%calving(:,:) * thk0 * rhoi / (model%numerics%dt * tim0)

    !--------------------------------------------------------------------
    ! Accumulate liquid runoff (basal melting)
    !--------------------------------------------------------------------
    !TODO - Add internal melting for enthalpy case
                       
    ! Note on units: model%temper%bmlt has dimensionless units of ice thickness per unit time
    !                Multiply by thk0/tim0 to convert to meters ice per second
    !                Multiply by rhoi to convert to kg/m^2/s water equiv.

    ! Convert to kg/m^2/s
    rofl_tavg(:,:) = rofl_tavg(:,:)  &
                   + model%temper%bmlt(:,:) * thk0/tim0 * rhoi

    !--------------------------------------------------------------------
    ! Accumulate basal heat flux
    !--------------------------------------------------------------------

    ! Note on units: model%temper%ucondflx has units of W/m^2, positive down
    !                Flip the sign so that hflx_tavg is positive up.

    hflx_tavg(:,:) = hflx_tavg(:,:) &
                   - model%temper%ucondflx(:,:)

  end subroutine glint_accumulate_output_gcm

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glint_upscale

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
