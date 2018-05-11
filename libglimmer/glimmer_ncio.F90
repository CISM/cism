!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_ncio.F90 - part of the Community Ice Sheet Model (CISM)  
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

#define NCO outfile%nc
#define NCI infile%nc

module glimmer_ncio
  !> module for common netCDF I/O
  !> written by Magnus Hagdorn, 2004

  use glimmer_ncdf

  implicit none

  integer,parameter,private :: msglen=512
  
contains

  !*****************************************************************************
  ! netCDF output
  !*****************************************************************************  
  subroutine openall_out(model,outfiles)

    !> open all netCDF files for output
    use glide_types
    use glimmer_ncdf
    use glimmer_filenames, only: process_path
    use parallel, only: parallel_open

    implicit none

    type(glide_global_type) :: model
    type(glimmer_nc_output), pointer, optional :: outfiles
    
    ! local variables
    type(glimmer_nc_output), pointer :: oc
    integer :: status

    if (present(outfiles)) then
       oc => outfiles
    else
       oc => model%funits%out_first
    end if

    do while(associated(oc))

       if (oc%append) then   ! assume the file exists, and reopen it

          call glimmer_nc_openappend(oc,model)

       elseif (model%options%is_restart == RESTART_TRUE) then   ! reopen the file if it exists

          status = parallel_open(process_path(oc%nc%filename),NF90_WRITE,oc%nc%id)

          if (status == NF90_NOERR) then  ! file exists and is now open; append it

             oc%append = .true.
             call glimmer_nc_openappend(oc, model, already_open_in=.true.)

          else  ! file does not exist; create it

             call glimmer_nc_createfile(oc, model)

          endif

       else  ! assume the file does not exist; create it

          call glimmer_nc_createfile(oc, model)

       end if

       oc => oc%next

    end do

  end subroutine openall_out

  !------------------------------------------------------------------------------

  subroutine closeall_out(model,outfiles)

    !> close all netCDF files for output
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    type(glimmer_nc_output),pointer,optional :: outfiles

    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc => model%funits%out_first
    end if

    do while(associated(oc))
       oc => delete(oc)
    end do
    if (.not.present(outfiles)) model%funits%out_first=>NULL()

  end subroutine closeall_out

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_openappend(outfile, model, &
                                   already_open_in)

    !> open netCDF file for appending
    use parallel
    use glimmer_log
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    implicit none

    type(glimmer_nc_output), pointer :: outfile       !> structure containing output netCDF descriptor
    type(glide_global_type) :: model                  !> the model instance
    logical, intent(in), optional :: already_open_in  !> if true, then the file is already open 

    ! local variables
    integer :: status, timedimid, ntime
    character(len=msglen) :: message
    logical :: already_open   ! if true, the file is already open

    if (present(already_open_in)) then
       already_open = already_open_in
    else
       already_open = .false.
    endif

    ! open the existing netCDF file, if not already open
    if (.not. already_open) then
       status = parallel_open(process_path(NCO%filename),NF90_WRITE,NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
    endif

    call write_log_div
    write(message,*) 'Reopening file ',trim(process_path(NCO%filename)),' for output; '
    call write_log(trim(message))

    ! Find out when last time slice was written
    status = parallel_inq_dimid(NCO%id,'time',timedimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCO%id,timedimid,len=ntime)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Set timecounter
    outfile%timecounter = ntime+1

    write(message,*) '  Write every ', outfile%freq, ' years'
    call write_log(trim(message))

    ! Get time-related varids
    status = parallel_inq_varid(NCO%id,glimmer_nc_internal_time_varname,NCO%internal_timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_varid(NCO%id,glimmer_nc_time_varname,NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inq_varid(NCO%id,glimmer_nc_tstep_count_varname,NCO%tstep_count_var)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! Put dataset into define mode
    status = parallel_redef(NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! setting the size of the level and staglevel dimension
    NCO%nlevel = model%general%upn
    NCO%nstaglevel = model%general%upn-1
    NCO%nstagwbndlevel = model%general%upn ! MJH this is the max index, not the size

  end subroutine glimmer_nc_openappend

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_createfile(outfile, model, baseline_year)

    !> create a new netCDF file
    use parallel
    use glimmer_log
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_filenames
    implicit none

    type(glimmer_nc_output), pointer :: outfile     !> structure containing output netCDF descriptor
    type(glide_global_type) :: model                !> the model instance
    integer, intent(in), optional :: baseline_year  !> baseline year to use for time units - i.e., the year to use in the string,
                                                    !> 'common_year since YYYY-01-01'
                                                    !> if not provided, use year 1 (0001)

    ! local variables
    integer, parameter :: time_units_len = 128
    integer status
    integer mapid
    integer :: sub_baseline_year  ! local version of baseline_year
    character(len=time_units_len) :: time_units
    character(len=msglen) message

    if (present(baseline_year)) then
       sub_baseline_year = baseline_year
    else
       sub_baseline_year = 1
    end if

    ! create new netCDF file
    !WHL - Changed the following line to support large netCDF output files
!!    status = parallel_create(process_path(NCO%filename),NF90_CLOBBER,NCO%id)
    status = parallel_create(process_path(NCO%filename), ior(NF90_CLOBBER,NF90_64BIT_OFFSET), NCO%id)
    call nc_errorhandle(__FILE__,__LINE__,status)
    call write_log_div
    write(message,*) 'Opening file ', trim(process_path(NCO%filename)), ' for output; '
    call write_log(trim(message))

    if (outfile%write_init) then
       write(message,*) '  Write output at start of run and every ', outfile%freq, ' years'
    else
       write(message,*) '  Write output every ', outfile%freq, ' years'
    endif
    call write_log(trim(message))

    if (outfile%end_write < glimmer_nc_max_time) then
       write(message,*) '  Stop writing at ', outfile%end_write
       call write_log(trim(message))
    end if
    NCO%define_mode=.TRUE.

    ! writing meta data
    status = parallel_put_att(NCO%id, NF90_GLOBAL, 'Conventions', "CF-1.3")
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'title',trim(outfile%metadata%title))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'institution',trim(outfile%metadata%institution))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'source',trim(outfile%metadata%source))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'history',trim(outfile%metadata%history))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'references',trim(outfile%metadata%references))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'comment',trim(outfile%metadata%comment))
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NF90_GLOBAL,'configuration',trim(outfile%metadata%config))
    call nc_errorhandle(__FILE__,__LINE__,status)
  
    ! defining time dimension and variable
    status = parallel_def_dim(NCO%id,'time',NF90_UNLIMITED,NCO%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)

    !     time -- Model time
    ! (see note in glimmer_ncdf regarding the reason for having multiple time-related variables)
    call write_log('Creating variables internal_time, time, and tstep_count')

    status = parallel_def_var(NCO%id,glimmer_nc_internal_time_varname,&
         outfile%default_xtype,(/NCO%timedim/),NCO%internal_timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'long_name', &
         'Model time - internal representation')
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'standard_name', 'time')
    ! CISM currently assumes a noleap calendar - exactly 365 days. For now, we hard-code
    ! this assumption in the units (by hard-coding that we're using units of common_year:
    ! CF/Udunits defines common_year to be 365 days, whereas year means 365.242198781
    ! days) and the calendar attribute.
    !
    ! For internal time, the baseline year is meaningless - so arbitrarily use year 1.
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'units', 'common_year since 1-1-1 0:0:0')
    status = parallel_put_att(NCO%id, NCO%internal_timevar, 'calendar', 'noleap')

    status = parallel_def_var(NCO%id,glimmer_nc_time_varname,&
         outfile%default_xtype,(/NCO%timedim/),NCO%timevar)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%timevar, 'long_name', 'Model time')
    status = parallel_put_att(NCO%id, NCO%timevar, 'standard_name', 'time')
    ! CISM currently assumes a noleap calendar - exactly 365 days. For now, we hard-code
    ! this assumption in the units (by hard-coding that we're using units of common_year:
    ! CF/Udunits defines common_year to be 365 days, whereas year means 365.242198781
    ! days) and the calendar attribute.
    !
    ! For time units, we write the year in YYYY format (but allowing for more digits if
    ! sub_baseline_year is greater than 9999).
    write(time_units,'(a,i0.4,a)') 'common_year since ', sub_baseline_year, '-01-01 0:0:0'
    status = parallel_put_att(NCO%id, NCO%timevar, 'units', time_units)
    status = parallel_put_att(NCO%id, NCO%timevar, 'calendar', 'noleap')

    status = parallel_def_var(NCO%id,glimmer_nc_tstep_count_varname,&
         NF90_INT,(/NCO%timedim/),NCO%tstep_count_var)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_put_att(NCO%id, NCO%tstep_count_var, 'long_name', &
         'Time step count')
    status = parallel_put_att(NCO%id, NCO%tstep_count_var, 'units', '-')

    ! adding projection info
    if (glimmap_allocated(model%projection)) then
       status = parallel_def_var(NCO%id,glimmer_nc_mapvarname,NF90_CHAR,mapid)
       call nc_errorhandle(__FILE__,__LINE__,status)
       call glimmap_CFPutProj(NCO%id,mapid,model%projection)
    end if

    ! setting the size of the level and staglevel dimension
    NCO%nlevel = model%general%upn
    NCO%nstaglevel = model%general%upn-1
    NCO%nstagwbndlevel = model%general%upn ! MJH this is the max index, not the size

  end subroutine glimmer_nc_createfile

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_checkwrite(outfile,model,forcewrite,time,external_time)

    !> check if we should write to file
    use parallel
    use glimmer_log
    use glide_types
    use glimmer_filenames
    implicit none
    type(glimmer_nc_output), pointer :: outfile    
    type(glide_global_type) :: model
    logical forcewrite
    real(dp),optional :: time  ! time in years (written to 'internal_time')
    real(dp),optional :: external_time  ! time in years (written to 'time') (if not present, uses the same time as internal_time)
    ! external_time only has an effect if it's present in the first call to this routine for a given time

    character(len=msglen) :: message
    integer status
    real(dp) :: sub_time  ! local version of time (years)
    real(dp) :: sub_external_time  ! local version of external_time (years)
    integer :: nfreq      ! freq/tinc; write output every nfreq timesteps

    ! Check for optional time argument
    if (present(time)) then
       sub_time=time
    else
       sub_time=model%numerics%time
    end if

    if (present(external_time)) then
       sub_external_time = external_time
    else
       sub_external_time = sub_time
    end if

    ! check if we are still in define mode and if so leave it
    if (NCO%define_mode) then
       status = parallel_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    if (sub_time > NCO%processsed_time) then
       if (NCO%just_processed) then
          ! finished writing during last time step, need to increase counter...
          
          outfile%timecounter = outfile%timecounter + 1
          status = parallel_sync(NCO%id)
          call nc_errorhandle(__FILE__,__LINE__,status)
          NCO%just_processed = .FALSE.
       end if
    end if

    ! Compute the desired integer frequency for writing output (every nfreq timesteps), rounding if needed.
    ! Note: Both outfile%freq and model%general%tinc have units of years.
    !       If tinc does not divide evenly into freq, then output will be written at regular intervals,
    !        but not exactly at the user-desired frequency. For example, suppose tinc = 0.3 yr and freq = 1.0 yr.
    !       Then output will be written every 3 timesteps, since nint(1.0/0.3) = 3.
    !
    nfreq = nint(outfile%freq / model%numerics%tinc)

    if (nfreq == 0) then  ! freq < tinc/2
       nfreq = 1
       write(message,*) 'WARNING: output file frequency is smaller than timestep; writing output every timestep'
       call write_log(trim(message))
    endif

    ! Write output if any of the following is true:
    ! (1) forcewrite = T
    ! (2) tstep_count = 0 & write_init = T
    ! (3) tstep_count > 0 & mod(tstep_count,nfreq) = 0
    ! Note: write_init = T by default, but can be turned off in the config file (e.g., for restart files)

    if ( forcewrite .or.  &
        (model%numerics%tstep_count == 0 .and. outfile%write_init) .or.  &
        (model%numerics%tstep_count > 0 .and. mod(model%numerics%tstep_count, nfreq) == 0) ) then

       if (sub_time <= outfile%end_write .and. .not.NCO%just_processed) then
          call write_log_div
          write(message,*) 'Writing to file ', trim(process_path(NCO%filename)), ' at time ', sub_time
          call write_log(trim(message))

          NCO%processsed_time = sub_time

          ! write time and tstep_count
          status = parallel_put_var(NCO%id,NCO%internal_timevar,sub_time,(/outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_var(NCO%id,NCO%timevar,sub_external_time,(/outfile%timecounter/))
          call nc_errorhandle(__FILE__,__LINE__,status)
          status = parallel_put_var(NCO%id,NCO%tstep_count_var,model%numerics%tstep_count,(/outfile%timecounter/))
          NCO%just_processed = .TRUE.         
       end if

    end if

  end subroutine glimmer_nc_checkwrite

  !*****************************************************************************
  ! netCDF input
  !*****************************************************************************  

  subroutine openall_in(model)

    !> open all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none

    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ! open input files
    ic=>model%funits%in_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do

    ! open forcing files
    ic=>model%funits%frc_first
    do while(associated(ic))
       call glimmer_nc_openfile(ic,model)
       ic=>ic%next
    end do

  end subroutine openall_in

  !------------------------------------------------------------------------------

  subroutine closeall_in(model)

    !> close all netCDF files for input
    use glide_types
    use glimmer_ncdf
    implicit none
    type(glide_global_type) :: model
    
    ! local variables
    type(glimmer_nc_input), pointer :: ic

    ! Input files
    ic=>model%funits%in_first
    do while(associated(ic))
       ic=>delete(ic)
    end do
    model%funits%in_first=>NULL()

    ! Forcing files
    ic=>model%funits%frc_first
    do while(associated(ic))
       ic=>delete(ic)
    end do
    model%funits%frc_first=>NULL()

  end subroutine closeall_in

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_openfile(infile,model)

    !> open an existing netCDF file
    use glide_types
    use glimmer_map_CFproj
    use glimmer_map_types
    use glimmer_log
    use glimmer_paramets, only: len0
    use glimmer_filenames
    use parallel
    implicit none

    type(glimmer_nc_input), pointer :: infile    !> structure containg input netCDF descriptor
    type(glide_global_type) :: model             !> the model instance

    ! local variables
    integer dimsize, dimid, varid
    real, dimension(2) :: delta
    integer status    
    character(len=msglen) message
    
    real,parameter :: small = 1.e-6

    ! open netCDF file
    status = parallel_open(process_path(NCI%filename),NF90_NOWRITE,NCI%id)
    if (status /= NF90_NOERR) then
       call write_log('Error opening file '//trim(process_path(NCI%filename))//': '//nf90_strerror(status),&
            type=GM_FATAL,file=__FILE__,line=__LINE__)
    end if
    call write_log_div
    call write_log('opening file '//trim(process_path(NCI%filename))//' for input')

    ! getting projection, if none defined already
    if (.not.glimmap_allocated(model%projection)) model%projection = glimmap_CFGetProj(NCI%id)

    ! getting time dimension
    status = parallel_inq_dimid(NCI%id, 'time', NCI%timedim)
    call nc_errorhandle(__FILE__,__LINE__,status)
    ! get id of time variable
    status = parallel_inq_varid(NCI%id,glimmer_nc_internal_time_varname,NCI%internal_timevar)
    ! BACKWARDS_COMPATIBILITY(wjs, 2017-04-28) Older files may not have 'internal_time',
    ! so if we can't find that variable, fall back on 'time'.
    if (status /= NF90_NOERR) then
       status = parallel_inq_varid(NCI%id,glimmer_nc_time_varname,NCI%internal_timevar)
    end if
    call nc_errorhandle(__FILE__,__LINE__,status)

    ! getting length of time dimension and allocating memory for array containing times
    status = parallel_inquire_dimension(NCI%id,NCI%timedim,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    allocate(infile%times(dimsize))
    infile%nt=dimsize
    status = parallel_get_var(NCI%id,NCI%internal_timevar,infile%times)

    ! getting tstep_count
    status = parallel_inq_varid(NCI%id,glimmer_nc_tstep_count_varname,NCI%tstep_count_var)
    ! BACKWARDS_COMPATIBILITY(wjs, 2017-05-17) Older files may not have 'tstep_count', so
    ! only read it if present.
    if (status == NF90_NOERR) then
       allocate(infile%tstep_counts(infile%nt))
       status = parallel_get_var(NCI%id,NCI%tstep_count_var,infile%tstep_counts)
       call nc_errorhandle(__FILE__,__LINE__,status)
       infile%tstep_counts_read = .true.
    else
       infile%tstep_counts_read = .false.
    end if

    ! setting the size of the level and staglevel dimension
    NCI%nlevel = model%general%upn
    NCI%nstaglevel = model%general%upn-1
    NCI%nstagwbndlevel = model%general%upn !MJH This is the max index, not size

    ! checking if dimensions and grid spacing are the same as in the configuration file
    ! x1
    status = parallel_inq_dimid(NCI%id,'x1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize /= global_ewn) then
       write(message,*) 'Dimension x1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config dimension: ', dimsize, global_ewn
       call write_log(message,type=GM_FATAL)
    end if
    status = parallel_inq_varid(NCI%id,'x1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)

!WHL - mod to prevent code from crashing due to small roundoff error
!    if (abs(delta(2)-delta(1) - model%numerics%dew*len0) > small) then
    if (abs( (delta(2)-delta(1) - model%numerics%dew*len0) / (model%numerics%dew*len0) ) > small) then
       write(message,*) 'deltax1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config deltax: ', delta(2)-delta(1),model%numerics%dew*len0
       call write_log(message,type=GM_FATAL)
    end if

    ! x0
    !status = nf90_inq_dimid(NCI%id,'x0',dimid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (dimsize /= model%general%ewn-1) then
    !   write(message,*) 'Dimension x0 of file ',trim(process_path(NCI%filename)),' does not match with config dimension: ', &
    !        dimsize, model%general%ewn-1
    !   call write_log(message,type=GM_FATAL)
    !end if
    !status = nf90_inq_varid(NCI%id,'x0',varid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_get_var(NCI%id,varid,delta)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (abs(delta(2)-delta(1) - model%numerics%dew*len0) > small) then
    !   write(message,*) 'deltax0 of file '//trim(process_path(NCI%filename))//' does not match with config deltax: ', &
    !        delta(2)-delta(1),model%numerics%dew*len0
    !   call write_log(message,type=GM_FATAL)
    !end if

    ! y1
    status = parallel_inq_dimid(NCI%id,'y1',dimid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_inquire_dimension(NCI%id,dimid,len=dimsize)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (dimsize /= global_nsn) then
       write(message,*) 'Dimension y1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config dimension: ', dimsize, global_nsn
       call write_log(message,type=GM_FATAL)
    end if
    status = parallel_inq_varid(NCI%id,'y1',varid)
    call nc_errorhandle(__FILE__,__LINE__,status)
    status = parallel_get_var(NCI%id,varid,delta)
    call nc_errorhandle(__FILE__,__LINE__,status)


!WHL - mod to prevent code from crashing due to small roundoff error
!    if (abs(delta(2)-delta(1) - model%numerics%dns*len0) > small) then
    if (abs( (delta(2)-delta(1) - model%numerics%dns*len0) / (model%numerics%dns*len0) ) > small) then
       write(message,*) 'deltay1 of file '//trim(process_path(NCI%filename))// &
            ' does not match with config deltay: ', delta(2)-delta(1),model%numerics%dns*len0
       call write_log(message,type=GM_FATAL)
    end if
    
    ! y0
    !status = nf90_inq_dimid(NCI%id,'y0',dimid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_inquire_dimension(NCI%id,dimid,len=dimsize)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (dimsize /= model%general%nsn-1) then
    !   write(message,*) 'Dimension y0 of file '//trim(process_path(NCI%filename))//' does not match with config dimension: ',&
    !        dimsize, model%general%nsn-1
    !   call write_log(message,type=GM_FATAL)
    !end if
    !status = nf90_inq_varid(NCI%id,'y0',varid)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !status = nf90_get_var(NCI%id,varid,delta)
    !call nc_errorhandle(__FILE__,__LINE__,status)
    !if (abs(delta(2)-delta(1) - model%numerics%dns*len0) > small) then
    !   write(message,*) 'deltay0 of file '//trim(process_path(NCI%filename))//' does not match with config deltay: ',&
    !        delta(2)-delta(1),model%numerics%dns*len0
    !   call write_log(message,type=GM_FATAL)
    !end if
  
  ! Check that the number of vertical layers is the same, though it's asking for trouble
  ! to check whether the spacing is the same (don't want to put that burden on setup,
  ! plus f.p. compare has been known to cause problems here)
  status = parallel_inq_dimid(NCI%id,'level',dimid)
  ! If we couldn't find the 'level' dimension fail with a warning.
  ! We don't want to throw an error, as input files are only required to have it if they
  ! include 3D data fields.
  if (status == NF90_NOERR) then
        status = parallel_inquire_dimension(NCI%id, dimid, len=dimsize)
        call nc_errorhandle(__FILE__, __LINE__, status)
        if (dimsize /= model%general%upn .and. dimsize  /=  1) then
            write(message,*) 'Dimension level of file '//trim(process_path(NCI%filename))//&
                ' does not match with config dimension: ', &
                dimsize, model%general%upn
            call write_log(message,type=GM_FATAL)
        end if
  else
        call write_log("Input file contained no level dimension.  This is not necessarily a problem.", type=GM_WARNING)
  end if
  
  end subroutine glimmer_nc_openfile

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_checkread(infile,model,time)

    !> check if we should read from file
    !>
    !> If we're reading a restart file, then also sets model%numerics%tstart,
    !  model%numerics%time and model%numerics%tstep_count
    use glimmer_log
    use glide_types
    use glimmer_filenames

    implicit none

    type(glimmer_nc_input), pointer :: infile  !> structure containg output netCDF descriptor
    type(glide_global_type) :: model    !> the model instance
    real(dp),optional :: time           !> Optional alternative time

    character(len=msglen) :: message

    integer :: pos       ! to identify CISM standalone restart files
    integer :: pos_cesm  ! to identify CESM restart files

    real(dp) :: restart_time   ! time of restart (yr)

    ! Note: infile%nt = number of time slices in the file
    !       infile%current_time = current time index

    ! Parse the filename to see if it is a restart file (standalone or CESM)
    pos = index(infile%nc%filename,'.restart.') ! CISM naming convention for restart files
    pos_cesm = index(infile%nc%filename,'.r.')  ! CESM naming convention for restart files

    ! If a standalone file, then set current_time to the latest time slice
    ! (Not necessary for CESM restart files, which contain just one time slice)
    if (pos /= 0) then
       infile%current_time = infile%nt
    endif

    if (infile%current_time <= infile%nt) then
       if (.not.NCI%just_processed) then

          call write_log_div

          ! Reset model%numerics%tstart if reading a restart file
          !write(message,*) 'Check for restart:', trim(infile%nc%filename)
          !call write_log(message)

          if (pos /= 0 .or. pos_cesm /= 0) then   ! get the start time based on the current time slice

             restart_time = infile%times(infile%current_time)      ! years
             model%numerics%tstart = restart_time
             model%numerics%time = restart_time

             if (infile%tstep_counts_read) then
                model%numerics%tstep_count = infile%tstep_counts(infile%current_time)
             else
                ! BACKWARDS_COMPATIBILITY(wjs, 2017-05-17) Older files may not have
                ! 'tstep_count', so compute it ourselves here. We don't want to use this
                ! formulation in general because it is prone to roundoff errors.
                model%numerics%tstep_count = nint(model%numerics%time/model%numerics%tinc)
             end if

             write(message,*) 'Restart: New tstart, tstep_count =', model%numerics%tstart, model%numerics%tstep_count
             call write_log(message)

          endif  ! pos/=0 or pos_cesm/=0

          write(message,*) 'Reading time slice ',infile%current_time,'(',infile%times(infile%current_time),') from file ', &
               trim(process_path(NCI%filename)), ' at time ', sub_time(model, time)
          call write_log(message)
          NCI%just_processed = .TRUE.
          NCI%processsed_time = sub_time(model, time)

       end if  ! not just processed
    end if  ! current_time <= nt

    if (sub_time(model, time) > NCI%processsed_time) then
       if (NCI%just_processed) then
          ! finished reading during last time step, need to increase counter...
          infile%current_time = infile%current_time + 1
          NCI%just_processed = .FALSE.
       end if
    end if

  contains

    real(dp) function sub_time(model, time)
      ! Get the current time applicable to this subroutine. 
      ! If time is present, use that; otherwise use model%numerics%time
      !
      ! We need this function to avoid code duplication. We canNOT simply set a local
      ! sub_time variable variable at the start of glimmer_nc_checkread, because model
      ! %numerics%time can be updated in the midst of this routine... so we need to
      ! determine sub_time when it's actually needed, with this function.
      use glide_types
      implicit none
      type(glide_global_type) :: model    !> the model instance
      real(dp),optional :: time           !> Optional alternative time

      if (present(time)) then
         sub_time = time
      else
         sub_time = model%numerics%time
      end if

    end function sub_time

  end subroutine glimmer_nc_checkread

  !------------------------------------------------------------------------------

  subroutine check_for_tempstag(whichdycore, nc)

      ! Check for the need to output tempstag and update the output variables if needed.
      !
      ! For the glam/glissade dycore, the vertical temperature grid has an extra level.
      ! In that case, the netCDF output file should include a variable
      ! called tempstag(0:nz) instead of temp(1:nz). This subroutine is added for
      ! convenience to allow the variable "temp" to be specified in the config
      ! file in all cases and have it converted to "tempstag" when appropriate.
      ! MJH

      use glimmer_log
      use glide_types

      implicit none
      integer, intent(in) :: whichdycore
      type(glimmer_nc_stat) :: nc

      ! Locals
      integer :: i

      ! Check if tempstag should be output

      ! If both temp and tempstag are specified, temp will get converted to tempstag
      ! and then there will be two tempstags in the list, but that is ok because
      ! the parser ignores duplicate entries in the varlist.  
      ! (The check for the existence of variables looks like:    pos = index(NCO%vars,' acab ')  )

      !print *, "Original varstring:", varstring

      if (whichdycore/=DYCORE_GLIDE) then 
          ! We want temp to become tempstag
          i = index(nc%vars, " temp ")
          if (i > 0) then
            ! temp was specified - change it to tempstag
            ! If temp is listed more than once, this just changes the first instance
            nc%vars = nc%vars(1:i-1) // " tempstag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses temperature on a staggered vertical grid.' // &
              '  The netCDF output variable "temp" has been changed to "tempstag".' )
          endif 
          ! Now check if flwa needs to be changed to flwastag
          i = index(nc%vars, " flwa ") ! Look for flwa
          if (i > 0) then
            ! flwa was specified - change to flwastag
            nc%vars = nc%vars(1:i-1) // " flwastag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses flwa on a staggered vertical grid.' // &
            '  The netCDF output variable "flwa" has been changed to "flwastag".' )
          endif
          ! Now check if dissip needs to be changed to dissipstag
          i = index(nc%vars, " dissip ") ! Look for dissip
          if (i > 0) then
            ! dissip was specified - change to dissipstag
            nc%vars = nc%vars(1:i-1) // " dissipstag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses dissip on a staggered vertical grid.' // &
            '  The netCDF output variable "dissip" has been changed to "dissipstag".' )
          endif
      else  ! glide dycore
          ! We want tempstag to become temp
          i = index(nc%vars, " tempstag ")
          if (i > 0) then
            !Change tempstag to temp
            nc%vars = nc%vars(1:i-1) // " temp " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "tempstag" should not be used with the Glide dycore.' // &
              '  The netCDF output variable "tempstag" has been changed to "temp".' )
          endif
          ! We want flwastag to become flwa
          i = index(nc%vars, " flwastag ")
          if (i > 0) then
            !Change flwastag to flwa
            nc%vars = nc%vars(1:i-1) // " flwa " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "flwastag" should not be used with the Glide dycore.' // &
              '  The netCDF output variable "flwastag" has been changed to "flwa".' )
          endif
          ! We want dissipstag to become dissip
          i = index(nc%vars, " dissipstag ")
          if (i > 0) then
            !Change dissipstag to dissip
            nc%vars = nc%vars(1:i-1) // " dissip " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "dissipstag" should not be used with the Glide dycore.' // &
              '  The netCDF output variable "dissipstag" has been changed to "dissip".' )
          endif
      endif  ! whichdycore

      ! Copy any changes to vars_copy
      nc%vars_copy = nc%vars

  end subroutine check_for_tempstag

!------------------------------------------------------------------------------


end module glimmer_ncio

!------------------------------------------------------------------------------
