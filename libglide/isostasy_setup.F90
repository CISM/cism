! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  isostasy_setup.f90 - part of the Glimmer-CISM ice model  + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010
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

module isostasy_setup

  !*FD routines for setting up isostasy module

  use isostasy_types
  use glimmer_global, only : dp

  implicit none

contains
  
  subroutine isos_readconfig(isos,config)
    !*FD read isostasy configuration
    use glimmer_config
    implicit none
    type(isos_type) :: isos                !*FD structure holding isostasy configuration
    type(ConfigSection), pointer :: config !*FD structure holding sections of configuration file

    ! local variables
    type(ConfigSection), pointer :: section

    ! read isostasy section
    call GetSection(config,section,'isostasy')
    if (associated(section)) then
       isos%do_isos = .true.
       call GetValue(section,'lithosphere',isos%lithosphere)
       call GetValue(section,'asthenosphere',isos%asthenosphere)
       call GetValue(section,'relaxed_tau',isos%relaxed_tau)       
       call GetValue(section,'update',isos%period)      
    end if
    ! read elastic section
    call GetSection(config,section,'elastic lithosphere')
    if (associated(section)) then
       call GetValue(section,'flexural_rigidity',isos%rbel%d)
    end if

    ! Construct the list of necessary restart variables based on the config options 
    ! selected by the user in the config file (specific to the isos section - other sections,
    ! e.g. glint, glide, are handled separately by their setup routines).
    ! This is done regardless of whether or not a restart ouput file is going 
    ! to be created for this run, but this information is needed before setting up outputs.   MJH 1/17/13
    call define_isos_restart_variables(isos)

  end subroutine isos_readconfig

  subroutine isos_printconfig(isos)
    !*FD print isostasy configuration to log

    use parallel, only: tasks
    use glimmer_log
    implicit none
    type(isos_type) :: isos                !*FD structure holding isostasy configuration

    character(len=100) :: message

    if (isos%do_isos) then
       call write_log('Isostasy setup')
       call write_log('--------------')
       if (isos%lithosphere==0) then
          call write_log('using local lithosphere approximation')
       else if (isos%lithosphere==1) then
          if (tasks > 1) then
             call write_log('Error, elastic lithosphere not supported for multiple processors',GM_FATAL)
          endif
          call write_log('using elastic lithosphere approximation')
          write(message,*) ' flexural rigidity: ', isos%rbel%d
          call write_log(message)
          write(message,*) ' update period: ',isos%period
          call write_log(message)
       else
          call write_log('unknown lithosphere')
          call close_log
          stop
       end if

       if (isos%asthenosphere==0) then
          call write_log('using fluid mantle')
       else if (isos%asthenosphere==1) then
          call write_log('using relaxing mantle')
          write(message,*) ' characteristic time constant: ',isos%relaxed_tau
          call write_log(message)
       else
          call write_log('unknown mantle')
          call close_log
          stop
       end if
       call write_log('')
    end if
  end subroutine isos_printconfig

  subroutine define_isos_restart_variables(isos)
    !*FD This subroutine analyzes the isostasy options input by the user in the config file
    !*FD and determines which variables are necessary for an exact restart.  MJH 1/11/2013

    ! Please comment thoroughly the reasons why a particular variable needs to be a restart variable for a given config.

    use isostasy_types
    use glide_io, only: glide_add_to_restart_variable_list
    implicit none

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------
    type(isos_type), intent (in) :: isos  !*FD Derived type holding isotasy options and variables (we only need the options here)

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    !------------------------------------------------------------------------------------

    ! restart variables needed for isostasy calculation
    if (isos%do_isos) then
        call glide_add_to_restart_variable_list('relx')
        ! TODO I suspect that relx is only needed when asthenosphere=1 (relaxing mantle), but I'm not sure - this should be tested when isostasy implementation is finalized/tested.
    endif
  end subroutine define_isos_restart_variables

end module isostasy_setup
