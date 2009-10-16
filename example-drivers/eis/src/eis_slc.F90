! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  eis_slc.f90 - part of the GLIMMER ice model              + 
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

module eis_slc
  !*FD climate forcing similar to the old Edinburgh Ice Sheet model
  !*FD Magnus Hagdorn, June 2004
  !*FD this modules handles the sea level component

  use glimmer_ts
  use glimmer_global, only : fname_length

  type eis_slc_type
     !*FD parameters for the EIS sea level forcing
     character(len=fname_length) :: fname=''     !*FD name of file containing ELA ts
     type(glimmer_tseries) :: slc_ts             !*FD ELA time series 
  end type eis_slc_type
  
contains  
  subroutine eis_slc_config(config,slc)
    !*FD get SLC configuration from config file
    use glimmer_config
    use glimmer_filenames, only : filenames_inputname
    implicit none
    type(eis_slc_type)           :: slc     !*FD slc data
    type(ConfigSection), pointer :: config  !*FD structure holding sections of configuration file   
    ! local variables
    type(ConfigSection), pointer :: section

    slc%fname=''

    call GetSection(config,section,'EIS SLC')
    if (associated(section)) then
       call GetValue(section,'slc_file',slc%fname)
       if (trim(slc%fname).ne.'') then
          slc%fname = trim(filenames_inputname(slc%fname))
       end if
    end if
  end subroutine eis_slc_config

  subroutine eis_slc_printconfig(slc)
    !*FD print configuration to log
    use glimmer_log
    use glimmer_global, only : msg_length
    implicit none
    type(eis_slc_type)           :: slc     !*FD slc data
    ! local variables
    character(len=msg_length) :: message

    call write_log('EIS SLC')
    call write_log('-------')
    write(message,*) 'SLC file: ',trim(slc%fname)
    call write_log(message)
    call write_log('')
  end subroutine eis_slc_printconfig

  subroutine eis_init_slc(slc)
    !*FD initialise SLC forcing
    use glimmer_paramets, only: thk0
    implicit none
    type(eis_slc_type)           :: slc     !*FD slc data
    
    call glimmer_read_ts(slc%slc_ts,slc%fname)
    ! scale parameters
    slc%slc_ts%values = slc%slc_ts%values/thk0
  end subroutine eis_init_slc

  subroutine eis_eus(slc,model,time)
    !*FD calculate mass balance
    use glide_types
    use glimmer_global, only : rk
    implicit none
    type(eis_slc_type)        :: slc   !*FD slc data
    type(glide_global_type)   :: model !*FD model instance
    real(kind=rk), intent(in) :: time  !*FD current time

    call glimmer_ts_linear(slc%slc_ts,real(time),model%climate%eus)
  end subroutine eis_eus
end module eis_slc
