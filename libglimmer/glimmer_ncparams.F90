!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_ncparams.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glimmer_ncparams

  ! read netCDF I/O related configuration files
  ! written by Magnus Hagdorn, May 2004

  use glimmer_ncdf, only: glimmer_nc_meta

  implicit none
    
  private

  public :: glimmer_nc_readparams, default_metadata, handle_output, handle_input, configstring

  type(glimmer_nc_meta),save :: default_metadata
  character(10000) :: configstring


contains

  !------------------------------------------------------------------------------

  subroutine glimmer_nc_readparams(model,config)

    ! read netCDF I/O related configuration file
    use glide_types
    use glimmer_config
    use glimmer_global, only : fname_length

    implicit none

    type(glide_global_type)      :: model  ! model instance
    type(ConfigSection), pointer :: config ! structure holding sections of configuration file
    
    ! local variables
    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output => null()
    type(glimmer_nc_input), pointer :: input => null()
    type(glimmer_nc_input), pointer :: forcing => null()

    integer :: pos
    integer :: ierr
    character(len=fname_length) :: restart_filename
    character(len=256) :: message

    ! Note on restart files:
    ! If a file is listed in the 'CF restart' section, then it is added to the glimmer_nc_output data structure
    !  and written at the specified frequency.
    ! If model%options%is_restart = RESTART_TRUE, then the file listed in 'CF restart' (provided it exists) 
    !  is added to the glimmer_nc_input data structure, overriding any file listed in the 'CF input' section.
    !  The latest time slice will be read in.
    ! Thus when restarting the model, it is only necessary to set restart = RESTART_TRUE (i.e, restart = 1)
    !  in the config file; it is not necesssary to change filenames in 'CF input' or 'CF restart'.
    ! At most one file should be listed in the 'CF restart' section, and it should contain the string 'restart'
    ! If model%options%is_restart = RESTART_TRUE and there is no 'CF restart' section, then the model will restart
    !  from the file and time slice specified in the 'CF input' section. (This is the old Glimmer behavior.)

    ! get config string
    call ConfigAsString(config,configstring)

    ! get default meta data
    call GetSection(config,section,'CF default')
    if (associated(section)) then
       call handle_metadata(section, default_metadata, .true.)
    end if

    ! set up outputs
    call GetSection(config,section,'CF output')
    do while(associated(section))
       output => handle_output(section,output,configstring)
       if (.not.associated(model%funits%out_first)) then
          model%funits%out_first => output
       end if
       call GetSection(section%next,section,'CF output')
    end do

    ! set up restart output
    ! If there is a 'CF restart' section, the file listed there is added to the output list.
    ! Note: There should be only one 'CF restart' section. Duplicate sections will be ignored.
    !TODO: Check that there is only one 'CF restart' section, and abort if there are more than one.
    call GetSection(config,section,'CF restart')
    if (associated(section)) then
       output => handle_output(section,output,configstring)
       if (.not.associated(model%funits%out_first)) then
          model%funits%out_first => output
       end if

       ! Make sure the filename contains 'restart'
       pos = index(output%nc%filename,'restart')
       if (pos == 0) then
          call write_log ('Error, filename in CF restart section should include "restart"', GM_FATAL)
       endif
    endif

    ! set up inputs
    call GetSection(config,section,'CF input')
    do while(associated(section))
       input => handle_input(section,input)
       if (.not.associated(model%funits%in_first)) then
          model%funits%in_first => input
       end if
       call GetSection(section%next,section,'CF input')
    end do

    ! set up restart input
    if (model%options%is_restart == RESTART_TRUE) then

       ! If there is a 'CF restart' section, the model will restart from the file listed there (if it exists).
       ! Else the model will start from the input file in the 'CF input' section.

       call GetSection(config,section,'CF restart')

       if (associated(section)) then

          ! get filename
          call GetValue(section, 'name', restart_filename)

          ! check whether a file with this name exists
          open(unit=1, file=trim(restart_filename), status='old', iostat=ierr)

          if (ierr == 0) then   ! file exists; set input pointers to the file in 'CF restart'

             ! nullify the input data structure set above
             input => null()
             model%funits%in_first => null()

             ! set new pointers
             input => handle_input(section,input)
             model%funits%in_first => input

             ! Make sure the filename contains 'restart'
             pos = index(input%nc%filename,'restart')
             if (pos == 0) then
                call write_log ('Error, filename in CF restart section should include "restart"', GM_FATAL)
             endif

             ! Make sure there is only one 'CF restart' section
             if (associated(section%next)) then
                call write_log ('Error, there should not be more than one CF restart section', GM_FATAL)
             endif

             write(message,*) 'Starting from restart file:', trim(input%nc%filename)
             call write_log(message)

          else   ! file does not exist; do not reset input pointers

             write(message,*) 'Cannot find restart file:', trim(restart_filename)
             call write_log(message)
             write(message,*) 'Starting from input file:', trim(input%nc%filename)
             call write_log(message)

          endif  ! ierr = 0 (restart file exists)

          close(unit=1)

       endif   ! associated(section)

    endif  ! model%options%is_restart = RESTART_TRUE

    ! setup forcings
    call GetSection(config,section,'CF forcing')
    do while(associated(section))
       forcing => handle_forcing(section,forcing)
       if (.not.associated(model%funits%frc_first)) then
          model%funits%frc_first => forcing
       end if
       call GetSection(section%next,section,'CF forcing')
    end do

    output => null()
    input => null()
    forcing => null()

  end subroutine glimmer_nc_readparams

  !==================================================================================
  ! private procedures
  !==================================================================================

  subroutine handle_metadata(section,metadata, default)
    use glimmer_ncdf
    use glimmer_config
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
  
  !------------------------------------------------------------------------------

  function handle_output(section, output, configstring)

    use glimmer_ncdf
    use glimmer_config
    use glimmer_log
    use glimmer_global, only: dp

    implicit none

    type(ConfigSection), pointer :: section
    type(glimmer_nc_output), pointer :: output
    type(glimmer_nc_output), pointer :: handle_output
    character(*),intent(in) :: configstring
    character(10) :: mode_str, xtype_str

    handle_output => add(output)
    
    mode_str = ''
    xtype_str = 'real'

    ! get filename and other info from config file
    call GetValue(section, 'name', handle_output%nc%filename)
    call GetValue(section, 'stop', handle_output%end_write)
    call GetValue(section, 'frequency', handle_output%freq)
    call GetValue(section, 'variables', handle_output%nc%vars)
    call GetValue(section, 'write_init', handle_output%write_init)
    call GetValue(section, 'mode', mode_str)
    call GetValue(section, 'xtype', xtype_str)

    ! handle mode field
    if (trim(mode_str)=='append' .or. trim(mode_str)=='APPEND') then
       handle_output%append = .true.
    else
       handle_output%append = .false.
    end if

    ! handle xtype field
    if (trim(xtype_str)=='real' .or. trim(xtype_str)=='REAL') then
       handle_output%default_xtype = NF90_REAL
    else if (trim(xtype_str)=='double' .or. trim(xtype_str)=='DOUBLE') then
       handle_output%default_xtype = NF90_DOUBLE
    else
       call write_log('Error, unknown xtype, must be real or double [netCDF output]',GM_FATAL)
    end if  

    ! add config data
    handle_output%metadata%config = trim(configstring)

    ! Make copy of variables for future reference
    handle_output%nc%vars_copy = handle_output%nc%vars

    ! get metadata
    call handle_metadata(section, handle_output%metadata, .false.)
    if (handle_output%nc%filename(1:1) == ' ') then
       call write_log('Error, no file name specified [netCDF output]',GM_FATAL)
    end if

  end function handle_output
  
  !------------------------------------------------------------------------------

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
    call GetValue(section, 'name', handle_input%nc%filename)
    call GetValue(section, 'time', handle_input%get_time_slice)
    
    handle_input%current_time = handle_input%get_time_slice

    if (handle_input%nc%filename(1:1) == ' ') then
       call write_log('Error, no file name specified [netCDF input]',GM_FATAL)
    end if
    
    handle_input%nc%filename = trim(filenames_inputname(handle_input%nc%filename))

  end function handle_input

  !------------------------------------------------------------------------------

  function handle_forcing(section, forcing)

    use glimmer_ncdf
    use glimmer_config
    use glimmer_log
    use glimmer_filenames, only : filenames_inputname
    implicit none
    type(ConfigSection), pointer :: section
    type(glimmer_nc_input), pointer :: forcing
    type(glimmer_nc_input), pointer :: handle_forcing

    handle_forcing=>add(forcing)
    
    ! get filename
    call GetValue(section,'name',handle_forcing%nc%filename)
    call GetValue(section,'time',handle_forcing%get_time_slice)  ! MJH don't think we'll use 'time' keyword in the forcing config section
    
    handle_forcing%current_time = handle_forcing%get_time_slice

    if (handle_forcing%nc%filename(1:1)==' ') then
       call write_log('Error, no file name specified [netCDF forcing]',GM_FATAL)
    end if

    handle_forcing%nc%filename = trim(filenames_inputname(handle_forcing%nc%filename))

  end function handle_forcing

  !------------------------------------------------------------------------------

end module glimmer_ncparams
