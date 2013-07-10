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

    call mean_to_global(instance%ups_orog, &
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

    call mean_to_global(instance%ups, &
                        temp, &
                        ice_frac,    &
                        instance%out_mask)

    ! Ice-with-snow fraction
    where (instance%mbal_accum%snowd > 0.d0 .and. instance%model%geometry%thck > 0.d0)
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        snowice_frac,    &
                        instance%out_mask)

    ! Veg-with-snow fraction (if ice <10m thick)
    where (instance%mbal_accum%snowd > 0.d0 .and. instance%model%geometry%thck <= (10.d0/thk0))
       temp = 1.d0
    elsewhere
       temp = 0.d0
    endwhere
    call mean_to_global(instance%ups, &
                        temp, &
                        snowveg_frac,    &
                        instance%out_mask)

    ! Remainder is veg only
    veg_frac = 1.d0 - ice_frac - snowice_frac - snowveg_frac

    ! Snow depth

    call mean_to_global(instance%ups, &
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
                                 gfrac,       gtopo,    &
                                 grofi,       grofl,    &
                                 ghflx)

    ! Upscale fields from the local grid to the global grid (with multiple elevation classes).
    ! Output fields are only valid on the main task.
    ! The upscaled fields are passed to the GCM land surface model, which has the option
    !  of updating the fractional area and surface elevation of glaciated gridcells.

    use glimmer_paramets, only: thk0, GLC_DEBUG
    use glimmer_log
    use parallel, only: tasks, main_task

    ! Arguments ----------------------------------------------------------------------------
 
    type(glint_instance), intent(inout) :: instance      ! the model instance
    integer,              intent(in)    :: nec           ! number of elevation classes
    integer,              intent(in)    :: nxl,nyl       ! local grid dimensions 
    integer,              intent(in)    :: nxg,nyg       ! global grid dimensions 

    !TODO - Should these be inout?
    !       Remove hardwired dimensions?
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gfrac   ! ice-covered fraction [0,1]
    real(dp),dimension(nxg,nyg,nec),intent(out) :: gtopo   ! surface elevation (m)
    real(dp),dimension(nxg,nyg,nec),intent(out) :: ghflx   ! heat flux (m)
    real(dp),dimension(nxg,nyg),    intent(out) :: grofi   ! ice runoff (calving) flux (kg/m^2/s)
    real(dp),dimension(nxg,nyg),    intent(out) :: grofl   ! liquid runoff (basal melt) flux (kg/m^2/s)
 
    ! Internal variables ----------------------------------------------------------------------
 
!    real(dp),dimension(nxl,nyl) :: local_field, local_topo, local_thck

    !TODO - Put this parameter elsewhere? 
    real(dp), parameter :: min_thck = 0.d0    ! min thickness (m) for setting gfrac = 1

    integer :: i, j, n      ! indices
    integer :: il, jl, ig, jg

    character(len=100) :: message

    real(dp) :: dew, dns    ! gridcell dimensions

    real(dp) ::   &
       usrf,               &! surface elevation (m)
       thck,               &! ice thickness (m)
       ucondflx,           &! conductive heat flux at upper surface (W m-2)
                            ! defined as positive down
       calving,            &! time-average calving rate (kg m-2 s-1)
       bmlting              ! time-average rate of basal/internal melting (kg m-2 s-1)
                            ! (melting at top surface is handled by GCM)

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

!WHL - debug
    print*, 'Start upscaling'

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
!!    if (GLC_DEBUG .and. tasks==1) then
       ig = iglint_global    ! defined in glint_type
       jg = jglint_global
       il = instance%model%numerics%idiag
       jl = instance%model%numerics%jdiag
       write(stdout,*) 'In glint_upscaling_gcm'
       write(stdout,*) 'il, jl =', il, jl
       write(stdout,*) 'ig, jg =', ig, jg
       write(stdout,*) 'nxl, nyl =', nxl,nyl
       write(stdout,*) 'nxg, nyg =', nxg,nyg
       write(stdout,*) 'av_count_output =', instance%av_count_output
       write(stdout,*) 'local out_mask =', instance%out_mask(il,jl)
       write(stdout,*) 'max, min out_mask =', maxval(instance%out_mask), minval(instance%out_mask)
       
!!    end if

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

    !TODO - Call a short subroutine to compute these?
    ! Compute time-average fluxes

    if (instance%av_count_output > 0) then
       instance%rofi_tavg(:,:) = instance%rofi_tavg(:,:) / real(instance%av_count_output,dp)
       instance%rofl_tavg(:,:) = instance%rofl_tavg(:,:) / real(instance%av_count_output,dp)
       instance%hflx_tavg(:,:) = instance%hflx_tavg(:,:) / real(instance%av_count_output,dp)
    else
       instance%rofi_tavg(:,:) = 0.d0
       instance%rofl_tavg(:,:) = 0.d0
       instance%hflx_tavg(:,:) = 0.d0
    endif

    ! Reset the logical variable for averaging output

    instance%new_tavg_output = .true.
    
    ! Loop over local grid cells
    ! Note: The area calculation is not strictly needed on a uniform grid, but is included
    !       to suppoort the general case of spatially varying grid dimensions.

    do j = 1, nyl
    do i = 1, nxl

       usrf = thk0 * instance%model%geometry%usrf(i,j)
       thck = thk0 * instance%model%geometry%thck(i,j)

!WHL - debug
       if (i==il .and. j==jl) then
          write(stdout,*) 'usrf =', usrf
          write(stdout,*) 'thck =', thck
          write(stdout,*) 'hflx =', instance%hflx_tavg(i,j)
          write(stdout,*) 'rofi =', instance%rofi_tavg(i,j)
          write(stdout,*) 'rofl =', instance%rofl_tavg(i,j)
       endif

       if (thck > min_thck) then   ! this cell is ice-covered

          do n = 1, nec
             if (usrf >= topomax(n-1) .and. usrf < topomax(n)) then  ! local cell is in elev class n
                area_frac_l(i,j,n) = dew*dns           ! Assume fractional coverage = 1.0
                area_topo_l(i,j,n) = dew*dns * usrf
                area_hflx_l(i,j,n) = dew*dns * instance%hflx_tavg(i,j)

                !WHL - debug
                if (i==il .and. j==jl) then
                   print*, ' '
                   print*, 'n =', n
                   print*, 'dew*dns =', dew*dns
                   print*, 'dew*dns*usrf =', area_topo_l(i,j,n)
                   print*, 'dew*dns*hflx =', area_hflx_l(i,j,n)
                endif

                exit
             endif
          enddo   ! nec
       endif      ! thck > min_thck

       ! Runoff fluxes, without nec loop
       ! Note: These fluxes can be nonzero for cells that no longer have ice.

       area_rofi_l(i,j) = dew*dns * instance%rofi_tavg(i,j)
       area_rofl_l(i,j) = dew*dns * instance%rofl_tavg(i,j)

       !WHL - debug
       if (i==il .and. j==jl) then
          print*, ' '
          print*, 'dew*dns*rofi =', area_rofi_l(i,j)
          print*, 'dew*dns*rofl =', area_rofl_l(i,j)
       endif

    enddo         ! i
    enddo         ! j

    ! Map the area-weighted local values to the global grid.

    ! Gridcell area

    call mean_to_global(instance%ups, &
                        area_l,       &
                        area_g,       &
                        instance%out_mask)

    ! Loop over elevation classes for gfrac, gtopo, ghflx

    do n = 1, nec

       ! area-weighted gfrac

       call mean_to_global(instance%ups,         &
                           area_frac_l(:,:,n),   &
                           area_frac_g(:,:,n),   &
                           instance%out_mask)

       ! area-weighted gtopo

       call mean_to_global(instance%ups,         &
                           area_topo_l(:,:,n),   &
                           area_topo_g(:,:,n),   &
                           instance%out_mask)

       ! area-weighted ghflx

       call mean_to_global(instance%ups,         &
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
          if (i==ig .and. j==jg) then
             print*, ' '
             print*, 'i, j, n =', i, j, n
             print*, 'area_g =', area_g(i,j)
             print*, 'area_frac_g =', area_frac_g(i,j,n)
             print*, 'area_topo_g =', area_topo_g(i,j,n)
             print*, 'area_hflx_g =', area_hflx_g(i,j,n)
             print*, 'gfrac =', gfrac(i,j,n)
             print*, 'gtopo =', gtopo(i,j,n)
             print*, 'ghflx =', ghflx(i,j,n)
          endif

       enddo  ! i
       enddo  ! j

    enddo     ! nec

!TODO - Modify mean_to_global subroutine so that it upscales correctly
!       regardless of whether ice is present in a given local cell

    ! area-weighted grofi

    call mean_to_global(instance%ups,       &
                        area_rofi_l(:,:),   &
                        area_rofi_g(:,:),   &
                        instance%out_mask)

    ! area-weighted grofl

    call mean_to_global(instance%ups,       &
                        area_rofl_l(:,:),   &
                        area_rofl_g(:,:),   &
                        instance%out_mask)


    ! Compute the mean ice and liquid runoff in each global grid cell

    do j = 1, nyg
    do i = 1, nxg

       ! Find mean ice runoff from calving 
       if (area_g(i,j) > 0.d0) then
          grofi(i,j) = area_rofi_g(i,j) / area_g(i,j)
       endif

       ! Find mean liquid runoff from internal/basal melting 
       if (area_g(i,j) > 0.d0) then      
          grofl(i,j) = area_rofl_g(i,j) / area_g(i,j)
       endif

!WHL - debug
       if (i==ig .and. j==jg) then
          print*, ' '
          print*, 'i, j =', i, j
          print*, 'area_g =', area_g(i,j)
          print*, 'area_rofi_g =', area_rofi_g(i,j)
          print*, 'area_rofl_g =', area_rofl_g(i,j)
          print*, 'grofi =', grofi(i,j)
          print*, 'grofl =', grofl(i,j)
       endif

    enddo   ! i
    enddo   ! j

!WHL - Old code follows

    if (GLC_DEBUG .and. main_task) then

!       write(stdout,*) ' '
!       write(stdout,*) 'global gfrac:'
!       do n = 1, nec
!          write(stdout,*) n, gfrac(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global gtopo:'
!       do n = 1, nec
!          write(stdout,*) n, gtopo(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global ghflx:'
!       do n = 1, nec
!          write(stdout,*) n, ghflx(ig, jg, n)
!       enddo

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofi:'
!       write(stdout,*) grofi(ig, jg)

!       write(stdout,*) ' '
!       write(stdout,*) 'global grofl:'
!       write(stdout,*) grofl(ig, jg)

    end if

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

    type(glide_global_type), intent(in)  :: model

    integer,  intent(inout) :: av_count_output     ! step counter 
    logical,  intent(inout) :: new_tavg_output     ! if true, start new averaging
    real(dp), dimension(:,:), intent(inout) :: rofi_tavg    ! solid ice runoff (kg m-2 s-1)
    real(dp), dimension(:,:), intent(inout) :: rofl_tavg    ! liquid runoff from basal/interior melting (kg m-2 s-1)
    real(dp), dimension(:,:), intent(inout) :: hflx_tavg    ! conductive heat flux at top surface (W m-2)

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

    ! Accumulate
    !TODO - Make sure this works with the real Glide fields
                                                                                                                               
!!    if (associated(model%climate%calving)) &
!!         rofi_tavg(:,:) = rofi_tavg(:,:) + model%climate%calving(:,:)  !TODO - Check units

!!    if (associated(model%temper%bmlt)) &
!!         rofl_tavg(:,:) = rofl_tavg(:,:) + model%temper%bmlt(:,:)

    !TODO - Associate and compute uncondflx for SIA dycore
    !       Currently computed only for HO dycore
 
!!    if (associated(model%temper%ucondflx)) &
!!         hflx_tavg(:,:) = hflx_tavg(:,:) + model%temper%ucondflx(:,:)

!WHL - debug - fake values for now

    rofi_tavg(:,:) = rofi_tavg(:,:) + 1.d0
    rofl_tavg(:,:) = rofl_tavg(:,:) + 1.d-1
    hflx_tavg(:,:) = hflx_tavg(:,:) - 1.d0

  end subroutine glint_accumulate_output_gcm

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  end module glint_upscale

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
