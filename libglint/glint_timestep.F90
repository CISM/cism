! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_timestep.f90 - part of the GLIMMER ice model       + 
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

module glint_timestep
  !*FD timestep of a GLINT instance

  use glint_type
  use glint_constants
  private
  public glint_i_tstep

contains

  subroutine glint_i_tstep(time,            instance,       &
                           g_temp,          g_temp_range,   &
                           g_precip,        g_zonwind,      &
                           g_merwind,       g_humid,        &
                           g_lwdown,        g_swdown,       &
                           g_airpress,      g_orog,         &
                           g_orog_out,      g_albedo,       &
                           g_ice_frac,      g_veg_frac,     &
                           g_snowice_frac,  g_snowveg_frac, &
                           g_snow_depth,                    &
                           g_water_in,      g_water_out,    &
                           t_win,           t_wout,         &
                           ice_vol,         out_f,          &
                           orogflag,        ice_tstep,      &
                           ccsm_smb_in,                     &
                           qice_g,          tsfc_g,         &
                           topo_g,          gmask,          &
                           gfrac,           gtopo,          &
                           grofi,           grofl,          &
                           ghflx )

    !*FD Performs time-step of an ice model instance. 
    !*FD Note that input quantities here are accumulated/average totals since the
    !*FD last call.
    !
    ! to do - 
    ! This code will need to be altered to take account of the 
    ! energy-balance mass-balance model when it is completed.
    ! Some of these arguments (e.g., zonwind, merwind, humid, lwdown, swdown, airpress)
    ! should probably be made optional.
    
    use glide
    use glide_setup
    use glide_io
    use glimmer_paramets
    use glint_io
    use glint_mbal_io
    use glint_climate
    !EIB! use glint_routing
    use glimmer_routing
    use glimmer_log
    use glimmer_physcon, only: rhow,rhoi
    use glide_mask, only: glide_mask_ocean

    implicit none

    ! ------------------------------------------------------------------------  
    ! Arguments
    ! ------------------------------------------------------------------------  

    integer,                intent(in)   :: time         !*FD Current time in hours
    type(glint_instance), intent(inout)  :: instance     !*FD Model instance
    real(rk),dimension(:,:),intent(in)   :: g_temp       !*FD Global mean surface temperature field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_temp_range !*FD Global surface temperature half-range field ($^{\circ}$C)
    real(rk),dimension(:,:),intent(in)   :: g_precip     !*FD Global precip field total (mm)

    real(rk),dimension(:,:),intent(in)   :: g_zonwind    !*FD Global mean surface zonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_merwind    !*FD Global mean surface meridonal wind (m/s)
    real(rk),dimension(:,:),intent(in)   :: g_humid      !*FD Global surface humidity (%)
    real(rk),dimension(:,:),intent(in)   :: g_lwdown     !*FD Global downwelling longwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_swdown     !*FD Global downwelling shortwave (W/m^2)
    real(rk),dimension(:,:),intent(in)   :: g_airpress   !*FD Global surface air pressure (Pa)

    real(rk),dimension(:,:),intent(in)   :: g_orog       !*FD Input global orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_orog_out   !*FD Output orography (m)
    real(rk),dimension(:,:),intent(out)  :: g_albedo     !*FD Output surface albedo 
    real(rk),dimension(:,:),intent(out)  :: g_ice_frac   !*FD Output ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_veg_frac   !*FD Output veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowice_frac !*FD Output snow-ice fraction
    real(rk),dimension(:,:),intent(out)  :: g_snowveg_frac !*FD Output snow-veg fraction
    real(rk),dimension(:,:),intent(out)  :: g_snow_depth !*FD Output snow depth (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_in   !*FD Input water flux (m)
    real(rk),dimension(:,:),intent(out)  :: g_water_out  !*FD Output water flux (m)
    real(rk),               intent(out)  :: t_win        !*FD Total water input (kg)
    real(rk),               intent(out)  :: t_wout       !*FD Total water output (kg)
    real(rk),               intent(out)  :: ice_vol      !*FD Output ice volume (m$^3$)
    type(output_flags),     intent(in)   :: out_f        !*FD Flags to tell us whether to do output   
    logical,                intent(in)   :: orogflag     !*FD Set if we have new global orog
    logical,                intent(out)  :: ice_tstep    !*FD Set if we have done an ice time step
    logical,                  optional,intent(in)  :: ccsm_smb_in ! true if getting sfc mass balance from CCSM
    real(rk),dimension(:,:,:),optional,intent(in)  :: qice_g    ! Depth of new ice (m)
    real(rk),dimension(:,:,:),optional,intent(in)  :: tsfc_g    ! Surface temperature (C)
    real(rk),dimension(:,:,:),optional,intent(in)  :: topo_g    ! Surface elevation (m)
    integer, dimension(:,:),  optional,intent(in)  :: gmask     ! = 1 where global data are valid, else = 0
    real(rk),dimension(:,:,:),optional,intent(out) :: gfrac     ! ice fractional area [0,1]
    real(rk),dimension(:,:,:),optional,intent(out) :: gtopo     ! surface elevation (m)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofi     ! ice runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: grofl     ! liquid runoff (kg/m^2/s = mm H2O/s)
    real(rk),dimension(:,:,:),optional,intent(out) :: ghflx     ! heat flux (W/m^2, positive down)

    ! ------------------------------------------------------------------------  
    ! Internal variables
    ! ------------------------------------------------------------------------  

    real(rk),dimension(:,:),pointer :: upscale_temp => null() ! temporary array for upscaling
    real(rk),dimension(:,:),pointer :: routing_temp => null() ! temporary array for flow routing
    real(rk),dimension(:,:),pointer :: accum_temp   => null() ! temporary array for accumulation
    real(rk),dimension(:,:),pointer :: ablat_temp   => null() ! temporary array for ablation
    integer, dimension(:,:),pointer :: fudge_mask   => null() ! temporary array for fudging
    real(sp),dimension(:,:),pointer :: thck_temp    => null() ! temporary array for volume calcs
    real(sp),dimension(:,:),pointer :: calve_temp   => null() ! temporary array for calving flux
    real(rk) :: start_volume,end_volume,flux_fudge
    integer :: i

#ifdef GLC_DEBUG
    integer :: j, ii, jj, nx, ny, il, jl
#endif

    logical :: ccsm_smb   ! true if getting sfc mass balance from CCSM

    if (present(ccsm_smb_in)) then
       ccsm_smb = ccsm_smb_in
    else
       ccsm_smb = .false.
    endif

    ! Zero outputs

    if (present(gfrac)) gfrac(:,:,:) = 0._rk
    if (present(gtopo)) gtopo(:,:,:) = 0._rk
    if (present(grofi)) grofi(:,:,:) = 0._rk
    if (present(grofl)) grofl(:,:,:) = 0._rk
    if (present(ghflx)) ghflx(:,:,:) = 0._rk

#ifdef GLC_DEBUG
    write(stdout,*) ' '
    write(stdout,*) 'In glint_i_tstep, time =', time
    write(stdout,*) 'next_time =', instance%next_time
#endif

    ! Check whether we're doing anything this time.

    if (time /= instance%next_time) then
       return
    else
       instance%next_time = instance%next_time + instance%mbal_tstep
    end if

    ! Assume we always need this, as it's too complicated to work out when we do and don't

    call coordsystem_allocate(instance%lgrid, thck_temp)
    call coordsystem_allocate(instance%lgrid, calve_temp)
    ice_tstep = .false.

    ! Downscale input fields from global to local grid
    ! This subroutine computes instance%acab and instance%artm, the key inputs to GLIDE.

    if (ccsm_smb) then
       call glint_downscaling_ccsm (instance,              &
                                    qice_g,      tsfc_g,   &
                                    topo_g,      gmask)
    else
       call glint_downscaling(instance,                  &
                              g_temp,     g_temp_range,  &
                              g_precip,   g_orog,        &
                              g_zonwind,  g_merwind,     &
                              g_humid,    g_lwdown,      &
                              g_swdown,   g_airpress,    &
                              orogflag)
    endif

    ! ------------------------------------------------------------------------  
    ! Sort out some local orography and remove bathymetry. This relies on the 
    ! point 1,1 being underwater. However, it's a better method than just 
    ! setting all points < 0.0 to zero
    ! ------------------------------------------------------------------------  

    call glide_get_usurf(instance%model, instance%local_orog)
    call glint_remove_bath(instance%local_orog,1,1)


    if (.not. ccsm_smb) then

       ! ------------------------------------------------------------------------  
       ! Adjust the surface temperatures using the lapse-rate, by reducing to
       ! sea-level and then back up to high-res orography
       ! ------------------------------------------------------------------------  

       call glint_lapserate(instance%artm, real(instance%global_orog,rk), real(-instance%data_lapse_rate,rk))
       call glint_lapserate(instance%artm, real(instance%local_orog,rk),  real(instance%lapse_rate,rk))

       ! Process the precipitation field if necessary ---------------------------
       ! and convert from mm/s to m/s

       call glint_calc_precip(instance)

    endif   ! not ccsm_smb

    ! Get ice thickness ----------------------------------------

!lipscomb -  Note that thck_temp is single-precision. [to do - does this matter for exact restart?] 

    call glide_get_thk(instance%model,thck_temp)

    ! Do accumulation --------------------------------------------------------

    if (ccsm_smb) then
       call glint_accumulate_ccsm(instance%mbal_accum, time, instance%acab, instance%artm)
    else
       call glint_accumulate(instance%mbal_accum, time, instance%artm, instance%arng, instance%prcp, &
                             instance%snowd, instance%siced, instance%xwind, instance%ywind, &
                             instance%local_orog, real(thck_temp,rk), instance%humid,    &
                             instance%swdown, instance%lwdown, instance%airpress)
    endif

    ! Initialise water budget quantities to zero. These will be over-ridden if
    ! there's an ice-model time-step

    t_win=0.0       ; t_wout=0.0
    g_water_out=0.0 ; g_water_in=0.0

#ifdef GLC_DEBUG
    write(stdout,*) ' '
    write(stdout,*) 'Check for ice dynamics timestep'
    write(stdout,*) 'time =', time
    write(stdout,*) 'start_time =', instance%mbal_accum%start_time
    write(stdout,*) 'mbal_step =', instance%mbal_tstep
    write(stdout,*) 'mbal_accum_time =', instance%mbal_accum_time
#endif

    ! ------------------------------------------------------------------------  
    ! ICE TIMESTEP begins HERE ***********************************************
    ! ------------------------------------------------------------------------  

    if (time - instance%mbal_accum%start_time + instance%mbal_tstep == instance%mbal_accum_time) then

       if (instance%mbal_accum_time < instance%ice_tstep) then 
          instance%next_time = instance%next_time + instance%ice_tstep - instance%mbal_tstep
       end if

       ice_tstep = .true.

       ! Prepare arrays for water budgeting

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
          call coordsystem_allocate(instance%lgrid, accum_temp)
          call coordsystem_allocate(instance%lgrid, ablat_temp)
          accum_temp = 0.0
          ablat_temp = 0.0
       end if

       ! Calculate the initial ice volume (scaled and converted to water equivalent)
       call glide_get_thk(instance%model, thck_temp)
       thck_temp = thck_temp*real(rhoi/rhow)
       start_volume = sum(thck_temp)

       ! ---------------------------------------------------------------------
       ! do the different parts of the glint timestep
       ! ---------------------------------------------------------------------

       do i = 1, instance%n_icetstep

#ifdef GLC_DEBUG
          write (stdout,*) 'GLIDE timestep, iteration =', i
#endif

          ! Calculate the initial ice volume (scaled and converted to water equivalent)
          call glide_get_thk(instance%model,thck_temp)
          thck_temp = thck_temp*real(rhoi/rhow)

          ! Get latest upper-surface elevation (needed for masking)
          call glide_get_usurf(instance%model, instance%local_orog)
          call glint_remove_bath(instance%local_orog,1,1)

          ! Get the mass-balance, as m water/year 
          call glint_get_mbal(instance%mbal_accum, instance%artm, instance%prcp, instance%ablt, &
                              instance%acab, instance%snowd, instance%siced, instance%mbal_accum_time)

          ! Mask out non-accumulation in ice-free areas

          where(thck_temp<=0.0 .and. instance%acab<0.0)
             instance%acab = 0.0
             instance%ablt = instance%prcp
          end where

          ! Set acab to zero for ocean cells (bed below sea level, no ice present)

          where (instance%model%geometry%thkmask == glide_mask_ocean)
             instance%acab = 0.0
          endwhere

          ! Put climate inputs in the appropriate places, with conversion ----------

          ! Note on units: 
          ! For this subroutine, input acab is in m/yr; this value is multiplied 
          !  by tim0/(scyr*thk0) and copied to data%climate%acab.
          ! Input artm is in deg C; this value is copied to data%climate%artm (no unit conversion).

          call glide_set_acab(instance%model, instance%acab*real(rhow/rhoi))
          call glide_set_artm(instance%model, instance%artm)

#ifdef GLC_DEBUG
             il = itest_local
             jl = jtest_local
             write (stdout,*) ' '
             write (stdout,*) 'After glide_set_acab, glide_set_artm: i, j =', il, jl
             write (stdout,*) 'acab (m/y), artm (C) =', instance%acab(il,jl), instance%artm(il,jl)
#endif

          ! Adjust glint acab and ablt for output
 
          where (instance%acab < -thck_temp .and. thck_temp > 0.0)
             instance%acab = -thck_temp
             instance%ablt =  thck_temp
          end where

          instance%glide_time=instance%glide_time+instance%model%numerics%tinc
          call glide_tstep_p1(instance%model,instance%glide_time)
          call glide_tstep_p2(instance%model)
          !EIB! difference in call between gc2 and lanl
          call glide_tstep_p3(instance%model)
          !EIB! call glide_tstep_p3(instance%model,no_write=.true.)

          ! Add the calved ice to the ablation field

!lipscomb - to do - Calving still needs to be added to the ice runoff.
!                   Also add basal melting (bmlt) to the liquid runoff.

          call glide_get_calving(instance%model, calve_temp)
          calve_temp = calve_temp * real(rhoi/rhow)

          instance%ablt = instance%ablt + calve_temp/instance%model%numerics%tinc
          instance%acab = instance%acab - calve_temp/instance%model%numerics%tinc

          ! Accumulate for water-budgeting
          if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then
             accum_temp = accum_temp + instance%prcp*instance%model%numerics%tinc
             ablat_temp = ablat_temp + instance%ablt*instance%model%numerics%tinc
          endif

          ! Save GLIDE output until now
          call glide_io_writeall(instance%model,instance%model)
          call glint_io_writeall(instance,instance%model)

       end do

       ! Calculate flux fudge factor --------------------------------------------

       if (out_f%water_out .or. out_f%total_wout .or. out_f%water_in .or. out_f%total_win) then

          call coordsystem_allocate(instance%lgrid,fudge_mask)

          call glide_get_thk(instance%model,thck_temp)
          end_volume=sum(thck_temp)

          where (thck_temp>0.0)
             fudge_mask = 1
          elsewhere
             fudge_mask = 0
          endwhere

          flux_fudge = (start_volume + sum(accum_temp) - sum(ablat_temp) - end_volume) / sum(fudge_mask)

          ! Apply fudge_factor

          where(thck_temp > 0.0)
             ablat_temp = ablat_temp + flux_fudge
          endwhere
          
          deallocate(fudge_mask)
          fudge_mask => null()

       endif

       ! Upscale water flux fields ----------------------------------------------
       ! First water input (i.e. mass balance + ablation)

       if (out_f%water_in) then
          call coordsystem_allocate(instance%lgrid, upscale_temp)

          where (thck_temp > 0.0)
             upscale_temp = accum_temp
          elsewhere
             upscale_temp = 0.0
          endwhere

          call mean_to_global(instance%ups,   &
                              upscale_temp,   &
                              g_water_in,     &
                              instance%out_mask)
          deallocate(upscale_temp)
          upscale_temp => null()
       endif

       ! Now water output (i.e. ablation) - and do routing

!lipscomb - to do - Modify (or skip?) the following for SMB option

       if (out_f%water_out) then
          call coordsystem_allocate(instance%lgrid, upscale_temp)
          call coordsystem_allocate(instance%lgrid, routing_temp)

          where (thck_temp > 0.0)
             upscale_temp = ablat_temp
          elsewhere
             upscale_temp = 0.0
          endwhere

          call glide_get_usurf(instance%model, instance%local_orog)
          call flow_router(instance%local_orog, &
                           upscale_temp, &
                           routing_temp, &
                           instance%out_mask, &
                           real(instance%lgrid%delta%pt(1),rk), &
                           real(instance%lgrid%delta%pt(2),rk))

          call mean_to_global(instance%ups,   &
                              routing_temp,   &
                              g_water_out,    &
                              instance%out_mask)

          deallocate(upscale_temp,routing_temp)
          upscale_temp => null()
          routing_temp => null()

       endif

       ! Sum water fluxes and convert if necessary ------------------------------

       if (out_f%total_win) then
          t_win  = sum(accum_temp) * instance%lgrid%delta%pt(1)* &
                                     instance%lgrid%delta%pt(2)
       endif

       if (out_f%total_wout) then
          t_wout = sum(ablat_temp) * instance%lgrid%delta%pt(1)* &
                                     instance%lgrid%delta%pt(2)
       endif

    end if

    ! Output instantaneous values

    call glint_mbal_io_writeall(instance, instance%model,       &
                                outfiles = instance%out_first,  &
                                time = real(time*hours2years,sp))

    ! ------------------------------------------------------------------------ 
    ! Upscaling of output
    ! ------------------------------------------------------------------------ 

    ! We now upscale all fields at once...

    call get_i_upscaled_fields(instance, g_orog_out, g_albedo, g_ice_frac, g_veg_frac, &
                               g_snowice_frac, g_snowveg_frac, g_snow_depth)

    ! Calculate ice volume ---------------------------------------------------

    if (out_f%ice_vol) then
       call glide_get_thk(instance%model, thck_temp)
       ice_vol = sum(thck_temp) * instance%lgrid%delta%pt(1)* &
                                  instance%lgrid%delta%pt(2)
    endif

    ! Tidy up ----------------------------------------------------------------

    if (associated(accum_temp)) then 
       deallocate(accum_temp)
       accum_temp => null()
    end if

    if (associated(ablat_temp)) then
       deallocate(ablat_temp)
       ablat_temp => null()
    end if

    if (associated(calve_temp)) then
       deallocate(calve_temp)
       calve_temp => null()
    end if

    deallocate(thck_temp)
    thck_temp => null()

  end subroutine glint_i_tstep

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_remove_bath(orog,x,y)

    !*FD Sets ocean areas to zero height, working recursively from
    !*FD a known ocean point.

    real(sp),dimension(:,:),intent(inout) :: orog !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y  !*FD Location of starting point (index)

    integer :: nx,ny

    nx=size(orog,1) ; ny=size(orog,2)

    if (orog(x,y).lt.0.0) orog(x,y)=0.0
    call glint_find_bath(orog,x,y,nx,ny)

  end subroutine glint_remove_bath

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  recursive subroutine glint_find_bath(orog,x,y,nx,ny)

    !*FD Recursive subroutine called by {\tt glimmer\_remove\_bath}.

    real(sp),dimension(:,:),intent(inout) :: orog  !*FD Orography --- used for input and output
    integer,                intent(in)    :: x,y   !*FD Starting point
    integer,                intent(in)    :: nx,ny !*FD Size of array {\tt orography}

    integer,dimension(4) :: xi=(/ -1,1,0,0 /)
    integer,dimension(4) :: yi=(/ 0,0,-1,1 /)
    integer :: ns=4,i

    do i=1,ns
       if (x+xi(i).le.nx.and.x+xi(i).gt.0.and. &
            y+yi(i).le.ny.and.y+yi(i).gt.0) then
          if (orog(x+xi(i),y+yi(i)).lt.0.0) then
             orog(x+xi(i),y+yi(i))=0.0
             call glint_find_bath(orog,x+xi(i),y+yi(i),nx,ny)
          endif
       endif
    enddo

  end subroutine glint_find_bath

end module glint_timestep


