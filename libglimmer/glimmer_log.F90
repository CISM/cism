! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_log.f90 - part of the GLIMMER ice model          + 
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

module glimmer_log
  !*FD module providing file logging and error/message handling
  !*FD Six levels of message/error are defined:
  !*FD \begin{itemize}
  !*FD \item Diagnostic messages
  !*FD \item Timestep enumeration and related information
  !*FD \item Information messages
  !*FD \item Warning messages
  !*FD \item Error messages
  !*FD \item Fatal error messages
  !*FD \end{itemize}
  !*FD These are numbered 1--6, with increasing severity, and the level of
  !*FD message output may be set to output all messages, only those above a particular 
  !*FD severity, or none at all. It should be noted that even if all messages are
  !*FD turned off, the model will still halt if it encounters a fatal
  !*FD error!
  !*FD 
  !*FD The other point to note is that when calling the messaging routines,
  !*FD the numerical identifier of a message level should be replaced by the
  !*FD appropriate parameter:
  !*FD \begin{itemize}
  !*FD \item \texttt{GM\_DIAGNOSTIC}
  !*FD \item \texttt{GM\_TIMESTEP}
  !*FD \item \texttt{GM\_INFO}
  !*FD \item \texttt{GM\_WARNING}
  !*FD \item \texttt{GM\_ERROR}
  !*FD \item \texttt{GM\_FATAL}
  !*FD \end{itemize}

  use glimmer_global, only : fname_length,dirsep

  integer,parameter :: GM_DIAGNOSTIC = 1 !*FD Numerical identifier for diagnostic messages.
  integer,parameter :: GM_TIMESTEP   = 2 !*FD Numerical identifier for timestep messages.
  integer,parameter :: GM_INFO       = 3 !*FD Numerical identifier for information messages.
  integer,parameter :: GM_WARNING    = 4 !*FD Numerical identifier for warning messages.
  integer,parameter :: GM_ERROR      = 5 !*FD Numerical identifier for (non-fatal) error messages.
  integer,parameter :: GM_FATAL      = 6 !*FD Numerical identifier for fatal error messages.

  integer, parameter            :: GM_levels = 6
  logical, private, dimension(GM_levels) :: gm_show = .false.

  character(len=*), parameter, dimension(0:GM_levels), private :: msg_prefix = (/ &
       '* UNKNOWN      ', &
       '*              ', &
       '*              ', &
       '               ', &
       '* WARNING:     ', &
       '* ERROR:       ', &
       '* FATAL ERROR :' /)


  character(len=fname_length),private :: glimmer_logname !*FD name of log file
  integer,private :: glimmer_unit=6                      !*FD log unit

contains
  function logname(fname)
    !*FD derives name of log file from file name by stripping directories and appending .log
    implicit none
    character(len=*), intent(in) :: fname
    character(len=fname_length) :: logname
    
    character(len=*), parameter :: suffix='.log'
    integer i
    i = scan(fname,dirsep,.True.)
    if (i.ne.0) then
       logname = trim(fname(i+1:))//suffix
    else
       logname = trim(fname)//suffix
    end if
  end function logname
  
  subroutine open_log(unit,fname)
    !*FD opens log file
    implicit none
    integer, optional          :: unit   !*FD file unit to use
    character(len=*), optional :: fname  !*FD name of log file

    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    if (present(unit)) then
       glimmer_unit = unit
    end if
    if (present(fname)) then
       glimmer_logname = adjustl(trim(fname))
    else
       glimmer_logname = 'glide.log'
    end if

    if (glimmer_unit.ne.6) then
       open(unit=glimmer_unit,file=glimmer_logname,status='unknown')
    end if

    call date_and_time(date,time)
    call write_log_div
    write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Started logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
  end subroutine open_log

  subroutine write_log(message,type,file,line)
    !*FD write to log
    implicit none
    integer,intent(in),optional          :: type    !*FD Type of error to be generated (see list above).
    character(len=*),intent(in)          :: message !*FD message to be written
    character(len=*),intent(in),optional :: file    !*FD the name of the file which triggered the message
    integer,intent(in),optional          :: line    !*FD the line number at the which the message was triggered

    ! local variables
    character(len=250) :: msg
    integer :: local_type
    character(len=6) :: line_num

    local_type = 0
    if (present(type)) then
       if (type.ge.1 .or. type.le.GM_levels) then
          local_type = type
       end if
    else
       local_type = GM_INFO
    end if

    ! constructing message
    if (present(file) .and. present(line)) then
       write(line_num,'(I6)')line
       write(msg,*) trim(msg_prefix(local_type))//' (',trim(file),':',trim(adjustl(line_num)),') '//message
    else
       write(msg,*) trim(msg_prefix(local_type))//' '//message
    end if
    ! messages are always written to file log
    write(glimmer_unit,*) trim(msg)
    ! and maybe to std out
    if (local_type.ne.0) then
       if (gm_show(local_type)) write(*,*) trim(msg)
    end if
    ! stop logging if we encountered a fatal error
    if (local_type.eq.GM_FATAL) then
       call close_log
       stop
    end if
  end subroutine write_log

  subroutine write_log_div
    !*FD start a new section
    implicit none
    write(glimmer_unit,*) '*******************************************************************************'
  end subroutine write_log_div

  subroutine close_log
    !*FD close log file
    implicit none
    ! local variables
    character(len=8) :: date
    character(len=10) :: time

    call date_and_time(date,time)
    call write_log_div
    write(unit=glimmer_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") ' Finished logging at ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    call write_log_div
    
    close(glimmer_unit)
  end subroutine close_log

  subroutine sync_log
    !*FD synchronise log to disk
    implicit none
    close(glimmer_unit)
    open(unit=glimmer_unit,file=glimmer_logname, position="append", status='old')
  end subroutine sync_log

  subroutine glimmer_set_msg_level(level)
    !*FD Sets the output message level.
    integer, intent(in) :: level !*FD The message level (6 is all messages; 0 is no messages). 
    integer :: i

    do i=1,GM_levels
       if (i>(GM_levels-level)) then
          gm_show(i)=.true.
       else
          gm_show(i)=.false.
       endif
    enddo

  end subroutine glimmer_set_msg_level

  function glimmer_get_logunit()
    !*FD return glimmer log unit
    implicit none
    integer glimmer_get_logunit

    glimmer_get_logunit = glimmer_unit
  end function glimmer_get_logunit
  
end module glimmer_log
