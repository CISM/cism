! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                            +
! +  glimmer_ncparams.f90 - part of the Glimmer-CISM ice model + 
! +                                                            +
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2007, 2008, 2009, 2010
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

#define NCO outfile%nc
#define NCI infile%nc

module glimmer_ncparams
  !*FD read netCDF I/O related configuration files
  !*FD written by Magnus Hagdorn, May 2004
  use glimmer_ncdf, only: glimmer_nc_meta
    
  private
  public :: glimmer_nc_readparams, default_metadata, handle_output, handle_input, configstring

  type(glimmer_nc_meta),save :: default_metadata
  character(10000) :: configstring

  character(310) :: restart_variable_list=''    !*FD list of variables needed for a restart
!TODO change 310 to a variable - see glimmer_ncdf.F90 in the definition for type glimmer_nc_stat for other instances of this value.


contains
    subroutine glimmer_nc_readparams(model,config)
    !*FD read netCDF I/O related configuration file
    use glide_types
    use glimmer_config
    implicit none
    type(glide_global_type)      :: model  !*FD model instance
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output => null()
    type(glimmer_nc_input), pointer :: input => null()

    ! get config string
    call ConfigAsString(config,configstring)

    ! get default meta data
    call GetSection(config,section,'CF default')
    if (associated(section)) then
       call handle_metadata(section, default_metadata, .true.)
    end if

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file and processed by *_readconfig (e.g. glide, isos).
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/11/13
    call define_restart_variables(model%options, model%isos)

    ! setup outputs
    call GetSection(config,section,'CF output')
    do while(associated(section))
       output => handle_output(section,output,model%numerics%tstart,configstring)
       if (.not.associated(model%funits%out_first)) then
          model%funits%out_first => output
       end if
       call GetSection(section%next,section,'CF output')
    end do

    ! setup inputs
    call GetSection(config,section,'CF input')
    do while(associated(section))
       input => handle_input(section,input)
       if (.not.associated(model%funits%in_first)) then
          model%funits%in_first => input
       end if
       call GetSection(section%next,section,'CF input')
    end do
    
    output => null()
    input => null()



  end subroutine glimmer_nc_readparams

  !==================================================================================
  ! private procedures
  !==================================================================================

  subroutine handle_metadata(section,metadata, default)
    use glimmer_ncdf
    use glimmer_config
    !use glimmer_global, only: glimmer_version !EIB! glimmer_verision not module in gc2
    implicit none
    type(ConfigSection), pointer :: section
    type(glimmer_nc_meta) ::metadata
    logical :: default

    !EIB! from gc2, may have been replaced by glimmer_version about, or vice versa??

    character(len=100), external :: glimmer_version_char

    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    if (.not.default) then
       metadata%title = trim(default_metadata%title)
       metadata%institution = trim(default_metadata%institution)
       metadata%references = trim(default_metadata%references)
       metadata%comment = trim(default_metadata%comment)
    end if

    call GetValue(section,'title',metadata%title)
    call GetValue(section,'institution',metadata%institution)
    call GetValue(section,'references',metadata%references)
    call GetValue(section,'comment',metadata%comment)

    if (default) then
       call date_and_time(date,time)
       !EIB!metadata%source = 'Generated by '//trim(glimmer_version)
       metadata%source = 'Generated by '//trim(glimmer_version_char())
       write(metadata%history,fmt="(a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6,' : ',a)") date(1:4),date(5:6),date(7:8),&
            !EIB!time(1:2),time(3:4),time(5:10),trim(glimmer_version)
            time(1:2),time(3:4),time(5:10),trim(glimmer_version_char())
    else
       metadata%source = trim(default_metadata%source)
       metadata%history = trim(default_metadata%history)
    end if
  end subroutine handle_metadata
  

  function handle_output(section, output, start_yr, configstring)
    use glimmer_ncdf
    use glimmer_config
    use glimmer_log
!WHLTSTEP
    use glimmer_global, only: dp
    implicit none

    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output
    type(glimmer_nc_output), pointer :: handle_output
!WHLTSTEP - changed start_yr to dp
!    real, intent(in) :: start_yr
    real(dp), intent(in) :: start_yr
    character(*),intent(in) :: configstring
    character(10) :: mode_str,xtype_str

    handle_output=>add(output)
    
    handle_output%next_write = start_yr
    mode_str=''
    xtype_str = 'real'

    ! get filename
    call GetValue(section,'name',handle_output%nc%filename)
    call GetValue(section,'start',handle_output%next_write)
    call GetValue(section,'stop',handle_output%end_write)
    call GetValue(section,'frequency',handle_output%freq)
    call GetValue(section,'variables',handle_output%nc%vars)
    call GetValue(section,'mode',mode_str)
    call GetValue(section,'xtype',xtype_str) !EIB! from gc2, needed?

    ! handle mode field
    if (trim(mode_str)=='append'.or.trim(mode_str)=='APPEND') then
       handle_output%append = .true.
    else
       handle_output%append = .false.
    end if

    !EIB! from gc2
    ! handle xtype field
    if (trim(xtype_str)=='real'.or.trim(xtype_str)=='REAL') then
       handle_output%default_xtype = NF90_REAL
    else if (trim(xtype_str)=='double'.or.trim(xtype_str)=='DOUBLE') then
       handle_output%default_xtype = NF90_DOUBLE
    else
       call write_log('Error, unknown xtype, must be real or double [netCDF output]',GM_FATAL)
    end if  
    !EIB!

    ! Determine if this is a restart, and if so set necessary flags and variable names
    call setup_restart_file(handle_output)

    ! add config data
    handle_output%metadata%config=trim(configstring)

    ! Make copy of variables for future reference
    handle_output%nc%vars_copy=handle_output%nc%vars

    ! get metadata
    call handle_metadata(section, handle_output%metadata,.false.)
    if (handle_output%nc%filename(1:1)==' ') then
       call write_log('Error, no file name specified [netCDF output]',GM_FATAL)
    end if  
  end function handle_output
  

  function handle_input(section, input)
    use glimmer_ncdf
    use glimmer_config
    use glimmer_log
    use glimmer_filenames, only : filenames_inputname !EIB! not in lanl, which is newer?
    implicit none
    type(ConfigSection), pointer :: section
    type(glimmer_nc_input), pointer :: input
    type(glimmer_nc_input), pointer :: handle_input

    handle_input=>add(input)
    
    ! get filename
    call GetValue(section,'name',handle_input%nc%filename)
    call GetValue(section,'time',handle_input%get_time_slice)
    
    handle_input%current_time = handle_input%get_time_slice

    if (handle_input%nc%filename(1:1)==' ') then
       call write_log('Error, no file name specified [netCDF input]',GM_FATAL)
    end if
    
    !EIB! from gc2
    handle_input%nc%filename = trim(filenames_inputname(handle_input%nc%filename))

  end function handle_input


  subroutine define_restart_variables(options, isos)
    !*FD This subroutine analyzes the options input by the user in the config file
    !*FD and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.

    use glide_types
    use isostasy_types
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glide_options), intent (in) :: options  !*FD Derived type holding all model options
    type(isos_type), intent (in) :: isos  !*FD Derived type holding isotasy options and variables (we only need the options here)
    !character(*), intent (inout) :: restart_variable_list  !*FD list of variables needed to perform an exact restart - module variable

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    !This was the restart list as of 1/11/13 using the old hot=1 systme in glide_vars.def:
    !restart_variable_list=' lat  relx  tauf  thk  thkmask  topg  bheatflx  bmlt  bwat  uvel  vvel  wgrd  flwa  temp  litho_temp  age '

    ! Start with a space and a few variables that we always want - prognostic variables and b.c.
    ! topg - needed to reconstruct all other geometry fields
    ! thk - prognostic variable
    ! temp - prognostic variable
    ! flwa - in principal this could be reconstructed from temp.  However in the current 
    !        implementation of glide the flwa calculation occurs after temp evolution but 
    !        before thk evolution.  This means flwa is calculated from the current temp and 
    !        the old thk.  The old thk is not available on a restart (just the current thk).
    !        (thk is needed to calculate flwa for 1) a mask for where ice is, 2) correction for pmp.)
    ! Note: the conversion from temp/flwa to tempstag/flwastag (if necessary) happens in glide_io.F90
    ! bheatflx, artm, acab - boundary conditions.  Of course if these fields are 0 they don't need 
    !        to be in the restart file, but without adding a check for that we cannot assume any of them are.
    !        There are some options where artm would not be needed.  Logic could be added to make that distinction.
    restart_variable_list = 'topg thk temp flwa bheatflx artm acab'

    ! add dycore specific restart variables
    select case (options%whichdycore)

      case (DYCORE_GLIDE)
        ! thkmask - TODO is this needed?
        ! wgrd - this is needed by the temperature code on a restart and cannot
        !        be re-calculated because it requires time derivatives that are stored in thckwrk
        !        TODO sorting out wgrd and wvel calculations
        restart_variable_list = trim(restart_variable_list) // ' thkmask wgrd'

        ! slip option for SIA
        select case (options%whichbtrc)
          case (0)
            ! no restart variable needed when no-slip is chosen
          case default
            restart_variable_list = trim(restart_variable_list) // ' btrc'
        end select

        if (options%whichevol == EVOL_PSEUDO_DIFF) then
          ! uvel,vvel are needed for linear diffusion evolution.  Their values
          ! cannot be reconstructed from the current thk or diffu because their
          ! values were originaly calculated from the old thk, which is not available
          ! on a restart.
          restart_variable_list = trim(restart_variable_list) // ' uvel vvel'
        endif

      case (DYCORE_GLAM, DYCORE_GLISSADE)
        ! uvel,vvel - these are needed for an exact restart because we can only 
        !             recalculate them to within the picard/jfnk convergence tolerance.
        ! beta - b.c. needed for runs with sliding - could add logic to only include in that case
        ! TODO not sure if thkmask is needed for HO
        ! TODO not sure if wgrd is needed for HO
        restart_variable_list = trim(restart_variable_list) // ' uvel vvel beta thkmask wgrd'

    end select

    ! ==== Other non-dycore specific options ====

    ! basal water option
    select case (options%whichbwat)
      case (BWATER_NONE, BWATER_CONST)
        ! no restart needed
      case default
        ! restart needs to know bwat value
        restart_variable_list = trim(restart_variable_list) // ' bwat'
    end select

    ! if the GTHF calculation is enabled, litho_temp needs to be a restart variable
    if (options%gthf /= 0) then
        restart_variable_list = trim(restart_variable_list) // ' litho_temp'
    endif

    ! basal processes module - requires tauf for a restart
    if (options%which_bmod /= BAS_PROC_DISABLED ) then
        restart_variable_list = trim(restart_variable_list) // ' tauf'
    endif

    ! TODO not sure if wvel is needed when vertical_integration=1.  Restarts are not
    ! exact with that option, so more work needs to be done to figure that out.

    ! TODO bmlt was set as a restart variable, but I'm not sure when or if it is needed.

    ! restart variables needed for isostasy calculation
    if (isos%do_isos) then
        restart_variable_list = trim(restart_variable_list) // ' relx'
        ! TODO I suspect that relx is only needed when asthenosphere=1 (relaxing mantle), but I'm not sure - this should be tested when isostasy implementation is finalized/tested.
    endif

    ! TODO age should be a restart variable if it is an input variable.  
    ! Same goes for b.c. (bheatflxm, artm, acab) and any other tracers that get introduced.
    ! These could be included all the time (as I have down above for b.c.), or 
    ! we could add logic to only include them when they were in the input file.
    ! To do this, this subroutine would have to be moved to after where input files are read,
    ! glide_io_readall(), but before the output files are created, glide_io_createall()

    ! TODO lat is only needed for some climate drivers.  It is not needed for simple_glide.
    ! Need to add logic that will add it only when those drivers are used.

    ! start & end end with a space
    restart_variable_list = ' ' // trim(restart_variable_list) // ' '

  end subroutine define_restart_variables



  subroutine setup_restart_file(handle_output)
    !*FD This subroutine determines if a given output file is meant to be a restart file
    !*FD and if so sets variables and stats accordingly.  
    !*FD An output file is set as a restart file if the words 'restart' or 'hot' are included
    !*FD in the variable list.   The inclusion of one of these words ensures the file is setup
    !*FD properly to allow an exact restart.     MJH 1/11/2013

    use glide_types

    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(glimmer_nc_output), pointer :: handle_output !*FD Derived type holding output information about this netcdf output file

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------
    integer :: pos
    !------------------------------------------------------------------------------------

    ! Expanding restart variables: if 'restart' or 'hot' is present, we remove that
    ! word from the variable list, and flip the restartfile flag.  
    ! In CISM 2.0, 'restart' is the preferred name to represent restart variables, 
    ! but 'hot' is supported for backward compatibility.  Thus, we check for both.
    handle_output%nc%vars = ' '//trim(handle_output%nc%vars)//' '
    pos = index(handle_output%nc%vars,' restart ') 
    if (pos.ne.0) then
       handle_output%nc%vars = handle_output%nc%vars(:pos)//handle_output%nc%vars(pos+8:)
       handle_output%nc%restartfile = .true.
    end if
    pos = index(handle_output%nc%vars,' hot ') 
    if (pos.ne.0) then
       handle_output%nc%vars = handle_output%nc%vars(:pos)//handle_output%nc%vars(pos+4:)
       handle_output%nc%restartfile = .true.
    end if

    ! Now apply necessary changes if the file is a restart file.
    if (handle_output%nc%restartfile) then
       ! Expand the restart variable list 
       handle_output%nc%vars = trim(handle_output%nc%vars) // restart_variable_list ! (a module variable)
       ! Set the xtype to be double (required for an exact restart)
       handle_output%default_xtype = NF90_DOUBLE
    end if

    ! Note: the variables are actually created in the file in glide_io.F90.  
    ! That is also where conversion from temp/flwa get converted to tempstag/flwastag
  end subroutine setup_restart_file


end module glimmer_ncparams
