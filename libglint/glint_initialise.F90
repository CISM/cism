! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                            +
! +  glint_initialise.f90 - part of the Glimmer-CISM ice model + 
! +                                                            +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

module glint_initialise

  !*FD Initialise GLINT model instance

  use glint_type

  private
  public glint_i_initialise, glint_i_end, calc_coverage

contains

  subroutine glint_i_initialise(config,           instance,         &
                                grid,             grid_orog,        &
                                mbts,             idts,             &
                                need_winds,       enmabal,          &
                                force_start,      force_dt,         &
                                gcm_restart,      gcm_restart_file, &
                                gcm_config_unit)

    !*FD Initialise a GLINT ice model instance

    use glimmer_config
    use glint_global_grid
    use glint_io
    use glint_mbal_io
    use glimmer_ncio
    use glide
    use glimmer_log
    use glint_constants
    use glimmer_restart_gcm
    implicit none

    ! Arguments
    type(ConfigSection), pointer         :: config           !*FD structure holding sections of configuration file   
    type(glint_instance),  intent(inout) :: instance         !*FD The instance being initialised.
    type(global_grid),     intent(in)    :: grid             !*FD Global grid to use
    type(global_grid),     intent(in)    :: grid_orog        !*FD Global grid to use for orography
    integer,               intent(out)   :: mbts             !*FD mass-balance time-step (hours)
    integer,               intent(out)   :: idts             !*FD ice dynamics time-step (hours)
    logical,               intent(inout) :: need_winds       !*FD Set if this instance needs wind input
    logical,               intent(inout) :: enmabal          !*FD Set if this instance uses the energy balance
                                                             !    mass-bal model
    integer,               intent(in)    :: force_start      !*FD glint forcing start time (hours)
    integer,               intent(in)    :: force_dt         !*FD glint forcing time step (hours)
    logical,     optional, intent(in)    :: gcm_restart      !*FD logical flag to read from a hotstart file
    character(*),optional, intent(in)    :: gcm_restart_file !*FD hotstart filename for restart
    integer,     optional, intent(in)    :: gcm_config_unit  !*FD fileunit for reading config files

    ! Internal
    real(sp),dimension(:,:),allocatable :: thk
    integer :: config_fileunit, restart_fileunit

    config_fileunit = 99
    if (present(gcm_config_unit)) then
       config_fileunit = gcm_config_unit
    endif

    ! initialise model

    call glide_config(instance%model, config, config_fileunit)

    ! if this is a continuation run, then set up to read restart
    ! (currently assumed to be a CESM restart file)

    if (present(gcm_restart)) then

      if (gcm_restart) then

         if (present(gcm_restart_file)) then

            ! read the hotstart file
            call glimmer_read_restart_gcm(instance%model, gcm_restart_file)
            instance%model%options%hotstart = 1
 
         else

            call write_log('Missing gcm_restart_file when gcm_restart is true',&
                           GM_FATAL,__FILE__,__LINE__)

         endif

      endif
    endif

    call glide_initialise(instance%model)
    instance%ice_tstep = get_tinc(instance%model)*years2hours
    instance%glide_time = instance%model%numerics%tstart
    idts = instance%ice_tstep

    ! read glint configuration

    call glint_i_readconfig(instance, config)    
    call glint_i_printconfig(instance)    
 
    ! create glint variables for the glide output files
    call glint_io_createall(instance%model, data=instance)

    ! create instantaneous glint variables
    call openall_out(instance%model, outfiles=instance%out_first)
    call glint_mbal_io_createall(instance%model, data=instance, outfiles=instance%out_first)

    ! fill dimension variables
    call glide_nc_fillall(instance%model)
    call glide_nc_fillall(instance%model, outfiles=instance%out_first)

    ! Check we've used all the config sections

    call CheckSections(config)

    ! New grid and downscaling

    instance%lgrid = coordsystem_new(0.d0, 0.d0, &
                                     get_dew(instance%model), &
                                     get_dns(instance%model), &
                                     get_ewn(instance%model), &
                                     get_nsn(instance%model))

    call new_downscale(instance%downs, instance%model%projection, grid, &
                       instance%lgrid, mpint=(instance%use_mpint==1))

    ! Allocate arrays appropriately

    call glint_i_allocate(instance, grid%nx, grid%ny, grid_orog%nx, grid_orog%ny)

    ! Read data and initialise climate

    call glint_i_readdata(instance)

    ! New upscaling

    call new_upscale(instance%ups, grid, instance%model%projection, &
                     instance%out_mask,  instance%lgrid) ! Initialise upscaling parameters
    call new_upscale(instance%ups_orog, grid_orog, instance%model%projection, &
                     instance%out_mask, instance%lgrid) ! Initialise upscaling parameters

    ! Calculate coverage map

    call calc_coverage(instance%lgrid, &
                       instance%ups,   &             
                       grid,           &
                       instance%out_mask, &
                       instance%frac_coverage)

    ! Calculate coverage map for orog

    call calc_coverage(instance%lgrid, &               
                       instance%ups_orog,  &             
                       grid_orog,     &
                       instance%out_mask, &
                       instance%frac_cov_orog)

    ! initialise the mass-balance accumulation

    call glint_mbc_init(instance%mbal_accum, &
                        instance%lgrid, &
                        config,         &
                        instance%whichacab, &
                        instance%snowd, &
                        instance%siced, &
                        instance%lgrid%size%pt(1), &
                        instance%lgrid%size%pt(2), &
                        real(instance%lgrid%delta%pt(1),rk))

    instance%mbal_tstep=instance%mbal_accum%mbal%tstep
    mbts=instance%mbal_tstep

    instance%next_time = force_start-force_dt+instance%mbal_tstep

    if (GLC_DEBUG) then
       write (6,*) 'Called glint_mbc_init'
       write (6,*) 'mbal tstep =', mbts
       write (6,*) 'next_time =', instance%next_time
       write (6,*) 'start_time =', instance%mbal_accum%start_time
    endif

    ! Mass-balance accumulation length

    if (instance%mbal_accum_time == -1) then
       instance%mbal_accum_time = max(instance%ice_tstep,instance%mbal_tstep)
       if (GLC_DEBUG) then
!Set mbal_accum_time = mbal_tstep
! lipscomb - TO DO - Make it easy to run Glimmer/Glint for ~5 days, e.g. for CESM smoke tests,
!         with all major components exercised. 
!!          instance%mbal_accum_time = instance%mbal_tstep
!!          write(6,*) 'WARNING: Seting mbal_accum_time =', instance%mbal_accum_time
       endif
    end if

    if (instance%mbal_accum_time < instance%mbal_tstep) then
       call write_log('Mass-balance accumulation timescale must be as '//&
            'long as mass-balance time-step',GM_FATAL,__FILE__,__LINE__)
    end if

    if (mod(instance%mbal_accum_time,instance%mbal_tstep)/=0) then
       call write_log('Mass-balance accumulation timescale must be an '// &
            'integer multiple of the mass-balance time-step',GM_FATAL,__FILE__,__LINE__)
    end if

    if (.not.(mod(instance%mbal_accum_time,instance%ice_tstep)==0.or.&
         mod(instance%ice_tstep,instance%mbal_accum_time)==0)) then
       call write_log('Mass-balance accumulation timescale and ice dynamics '//&
            'timestep must divide into one another',GM_FATAL,__FILE__,__LINE__)
    end if

    if (instance%ice_tstep_multiply/=1 .and. mod(instance%mbal_accum_time,int(years2hours))/=0.0) then
       call write_log('For ice time-step multiplication, mass-balance accumulation timescale '//&
            'must be an integer number of years',GM_FATAL,__FILE__,__LINE__)
    end if

    ! Initialise some other stuff

    if (instance%mbal_accum_time>instance%ice_tstep) then
       instance%n_icetstep = instance%ice_tstep_multiply*instance%mbal_accum_time/instance%ice_tstep
    else
       instance%n_icetstep = instance%ice_tstep_multiply
    end if

!This was commented out because it destroys exact restart
!lipscomb - TO DO - Find another way to set thk to snowd?
    ! Copy snow-depth to thickness if no thickness is present

!!    allocate(thk(get_ewn(instance%model),get_nsn(instance%model)))
!!    call glide_get_thk(instance%model,thk)
!!    where (instance%snowd>0.0 .and. thk==0.0)
!!       thk=instance%snowd
!!    elsewhere
!!       thk=thk
!!    endwhere
!!    call glide_set_thk(instance%model,thk)
!!    deallocate(thk)

    call glide_io_writeall(instance%model, instance%model)
    call glint_io_writeall(instance, instance%model)
    call glint_mbal_io_writeall(instance, instance%model, outfiles=instance%out_first)

    if (instance%whichprecip == 2) need_winds=.true.
    if (instance%whichacab == 3) then
       need_winds = .true.
       enmabal = .true.
    end if

  end subroutine glint_i_initialise

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_end(instance)

    !*FD Performs tidying-up for an ice model. 

    use glide
    use glimmer_ncio
    implicit none
    type(glint_instance),  intent(inout) :: instance    !*FD The instance being initialised.

    call glide_finalise(instance%model)
    call closeall_out(instance%model,outfiles=instance%out_first)
    instance%out_first => null()

  end subroutine glint_i_end

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine glint_i_readdata(instance)
    !*FD read data from netCDF file and initialise climate
    use glimmer_log
    use glimmer_ncdf
    use glint_climate
    use glint_io
    use glide_setup
    use glide_temp
    implicit none

    type(glint_instance),intent(inout)   :: instance    !*FD Instance whose elements are to be allocated.

    ! read data
    call glint_io_readall(instance,instance%model)

    call glide_calclsrf(instance%model%geometry%thck,instance%model%geometry%topg, &
         instance%model%climate%eus,instance%model%geometry%lsrf)
    instance%model%geometry%usrf = instance%model%geometry%thck + instance%model%geometry%lsrf

  end subroutine glint_i_readdata

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_coverage(lgrid,ups,grid,mask,frac_coverage)

    !*FD Calculates the fractional
    !*FD coverage of the global grid-boxes by the ice model
    !*FD domain.

    use glimmer_map_types
    use glimmer_coordinates
    use glint_global_grid

    ! Arguments

    type(coordsystem_type), intent(in)  :: lgrid         !*FD Local grid
    type(upscale),          intent(in)  :: ups           !*FD Upscaling used
    type(global_grid),      intent(in)  :: grid          !*FD Global grid used
    integer, dimension(:,:),intent(in)  :: mask          !*FD Mask of points for upscaling
    real(rk),dimension(:,:),intent(out) :: frac_coverage !*FD Map of fractional 
                                                         !*FD coverage of global by local grid-boxes.
    ! Internal variables

    integer,dimension(grid%nx,grid%ny) :: tempcount
    integer :: i,j

    ! Beginning of code

    tempcount=0

    do i=1,lgrid%size%pt(1)
       do j=1,lgrid%size%pt(2)
          tempcount(ups%gboxx(i,j),ups%gboxy(i,j))=tempcount(ups%gboxx(i,j),ups%gboxy(i,j))+mask(i,j)
       enddo
    enddo

    do i=1,grid%nx
       do j=1,grid%ny
          if (tempcount(i,j)==0) then
             frac_coverage(i,j)=0.0
          else
             frac_coverage(i,j)=(tempcount(i,j)*lgrid%delta%pt(1)*lgrid%delta%pt(2))/ &
                  (lon_diff(grid%lon_bound(i+1),grid%lon_bound(i))*D2R*EQ_RAD**2*    &
                  (sin(grid%lat_bound(j)*D2R)-sin(grid%lat_bound(j+1)*D2R)))
          endif
       enddo
    enddo

    ! Fix points that should be 1.0 by checking their surroundings

    ! Interior points first

    do i=2,grid%nx-1
       do j=2,grid%ny-1
          if ((frac_coverage(i,j).ne.0).and. &
               (frac_coverage(i+1,j).ne.0).and. &
               (frac_coverage(i,j+1).ne.0).and. &
               (frac_coverage(i-1,j).ne.0).and. &
               (frac_coverage(i,j-1).ne.0)) &
               frac_coverage(i,j)=1.0
       enddo
    enddo

    ! top and bottom edges

    do i=2,grid%nx/2
       if ((frac_coverage(i,1).ne.0).and. &
           (frac_coverage(i+1,1).ne.0).and. &
           (frac_coverage(i,2).ne.0).and. &
           (frac_coverage(i-1,1).ne.0).and. &
           (frac_coverage(i+grid%nx/2,1).ne.0)) &
            frac_coverage(i,1)=1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
       if ((frac_coverage(i,1).ne.0).and. &
           (frac_coverage(i+1,1).ne.0).and. &
           (frac_coverage(i,2).ne.0).and. &
           (frac_coverage(i-1,1).ne.0).and. &
           (frac_coverage(i-grid%nx/2,1).ne.0)) &
            frac_coverage(i,1)=1.0
    enddo

    do i=2,grid%nx/2
       if ((frac_coverage(i,grid%ny).ne.0).and. &
           (frac_coverage(i+1,grid%ny).ne.0).and. &
           (frac_coverage(i+grid%nx/2,grid%ny).ne.0).and. &
           (frac_coverage(i-1,grid%ny).ne.0).and. &
           (frac_coverage(i,grid%ny-1).ne.0)) &
            frac_coverage(i,grid%ny)=1.0
    enddo

    do i=grid%nx/2+1,grid%nx-1
       if ((frac_coverage(i,grid%ny).ne.0).and. &
           (frac_coverage(i+1,grid%ny).ne.0).and. &
           (frac_coverage(i-grid%nx/2,grid%ny).ne.0).and. &
           (frac_coverage(i-1,grid%ny).ne.0).and. &
           (frac_coverage(i,grid%ny-1).ne.0)) &
            frac_coverage(i,grid%ny)=1.0
    enddo

    ! left and right edges

    do j=2,grid%ny-1
       if ((frac_coverage(1,j).ne.0).and. &
           (frac_coverage(2,j).ne.0).and. &
           (frac_coverage(1,j+1).ne.0).and. &
           (frac_coverage(grid%nx,j).ne.0).and. &
           (frac_coverage(1,j-1).ne.0)) &
            frac_coverage(1,j)=1.0
       if ((frac_coverage(grid%nx,j).ne.0).and. &
           (frac_coverage(1,j).ne.0).and. &
           (frac_coverage(grid%nx,j+1).ne.0).and. &
           (frac_coverage(grid%nx-1,j).ne.0).and. &
           (frac_coverage(grid%nx,j-1).ne.0)) &
            frac_coverage(grid%nx,j)=1.0
    enddo

    ! corners

    if ((frac_coverage(1,1).ne.0).and. &
        (frac_coverage(2,1).ne.0).and. &
        (frac_coverage(1,2).ne.0).and. &
        (frac_coverage(grid%nx,1).ne.0).and. &
        (frac_coverage(grid%nx/2+1,1).ne.0)) &
         frac_coverage(1,1)=1.0

    if ((frac_coverage(1,grid%ny).ne.0).and. &
        (frac_coverage(2,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx/2+1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx,grid%ny).ne.0).and. &
        (frac_coverage(1,grid%ny-1).ne.0)) &
         frac_coverage(1,grid%ny)=1.0

    if ((frac_coverage(grid%nx,1).ne.0).and. &
        (frac_coverage(1,1).ne.0).and. &
        (frac_coverage(grid%nx,2).ne.0).and. &
        (frac_coverage(grid%nx-1,1).ne.0).and. &
        (frac_coverage(grid%nx/2,1).ne.0)) &
         frac_coverage(grid%nx,1)=1.0

    if ((frac_coverage(grid%nx,grid%ny).ne.0).and. &
        (frac_coverage(1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx/2,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx-1,grid%ny).ne.0).and. &
        (frac_coverage(grid%nx,grid%ny-1).ne.0)) &
         frac_coverage(grid%nx,grid%ny)=1.0

    ! Finally fix any rogue points > 1.0

    where (frac_coverage>1.0) frac_coverage=1.0

  end subroutine calc_coverage

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function lon_diff(a,b)

    implicit none

    real(rk),intent(in) :: a,b
    real(rk) :: aa,bb

    aa=a ; bb=b

    do
       if (aa>bb) exit
       aa=aa+360.0
    enddo

    lon_diff=aa-bb

  end function lon_diff

end module glint_initialise
