! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  profile.f90 - part of the GLIMMER ice model              + 
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

module profile
  !*FD Magnus Hagdorn
  !*FD January 2005
  !*FD module for profiling programs

  integer, private :: current_unit = 200
  integer, private,parameter :: max_prof = 100

  type profile_type
     integer :: profile_unit=0  !*FD file unit to be written to
     real :: start_time         !*FD CPU time at start of log
     integer :: nump=0          !*FD number of profiles

     real, dimension(max_prof) :: pstart,ptotal      !*FD for each log store start and totals
     character(len=50), dimension(max_prof) :: pname !*FD name for each profile
  end type profile_type

  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_PROFILE
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_PROFILE
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_PROFILE
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_PROFILE
!MH!#endif
  
  subroutine profile_init(prof,name)
    !*FD initialise a profile
    implicit none
    type(profile_type), intent(out) :: prof !*FD structure storing profile definitions
    character(len=*), intent(in) :: name    !*FD name of file
    ! local variables
    character(len=8)  :: date
    character(len=10) :: time    

    prof%profile_unit = current_unit
    current_unit = current_unit + 1
    call cpu_time(prof%start_time)
    call date_and_time (date, time)
    open(unit=prof%profile_unit,file=name,status='unknown')
    write(unit=prof%profile_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") '# Started profile on ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
  end subroutine profile_init
  
  function profile_register(prof,msg)
    !*FD register a new series of meassurements
    use glimmer_log
    implicit none
    type(profile_type) :: prof !*FD structure storing profile definitions
    character(len=*), intent(in) :: msg
    integer profile_register

    prof%nump = prof%nump+1
    if (prof%nump .gt. max_prof) then
       call write_log('Maximum number of profiles reached',type=GM_FATAL, &
            file=__FILE__,line=__LINE__)
    end if
    profile_register = prof%nump
    prof%pname(prof%nump) = trim(msg)
  end function profile_register

  subroutine profile_start(prof,profn)
    !*FD start profiling
    implicit none
    type(profile_type) :: prof !*FD structure storing profile definitions
    integer, intent(in) :: profn
    
    call cpu_time(prof%pstart(profn))
  end subroutine profile_start

  subroutine profile_stop(prof,profn)
    !*FD stop profiling
    implicit none
    type(profile_type)  :: prof !*FD structure storing profile definitions
    integer, intent(in) :: profn
    
    real t
    call cpu_time(t)
    prof%ptotal(profn) = prof%ptotal(profn) + t-prof%pstart(profn)
  end subroutine profile_stop

  subroutine profile_log(prof,profn,msg)
    !*FD log a message to profile
    implicit none
    type(profile_type)           :: prof !*FD structure storing profile definitions
    integer, intent(in) :: profn
    character(len=*), intent(in), optional :: msg     !*FD message to be written to profile

    real t

    call cpu_time(t)
    if (present(msg)) then
       write(prof%profile_unit,*) t-prof%start_time,prof%ptotal(profn),profn,trim(msg)//' '//trim(prof%pname(profn))
    else
       write(prof%profile_unit,*) t-prof%start_time,prof%ptotal(profn),profn,trim(prof%pname(profn))
    end if
    prof%ptotal(profn) = 0.
    prof%pstart(profn) = 0.
  end subroutine profile_log

  subroutine profile_close(prof)
    !*FD close profile
    implicit none
    type(profile_type), intent(in) :: prof !*FD structure storing profile definitions
    ! local variables
    character(len=8)  :: date
    character(len=10) :: time    
    real t

    call cpu_time(t)
    call date_and_time (date, time)
    write(prof%profile_unit,*) '# total elapse cpu time: ',t-prof%start_time
    write(unit=prof%profile_unit,fmt="(a,a4,'-',a2,'-',a2,' ',a2,':',a2,':',a6)") '# Finished profile on ',&
         date(1:4),date(5:6),date(7:8),time(1:2),time(3:4),time(5:10)
    close(prof%profile_unit)
  end subroutine profile_close
end module profile
