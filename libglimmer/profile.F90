! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  profile.f90 - part of the Glimmer-CISM ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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

!> module for profiling programs
!! \author Magnus Hagdorn
!! \date January 2005
module profile

  integer, private :: current_unit = 200
  integer, private,parameter :: max_prof = 100

  !> the profiling type
  type profile_type
     integer :: profile_unit=0  !< file unit to be written to
     real :: start_time         !< CPU time at start of log
     integer :: nump=0          !< number of profiles

     real, dimension(max_prof) :: pstart,ptotal      !< for each log store start and totals
     character(len=50), dimension(max_prof) :: pname !< name for each profile
  end type profile_type

contains

  !> initialise a profile
  subroutine profile_init(prof,name)
    implicit none
    type(profile_type), intent(out) :: prof !< structure storing profile definitions
    character(len=*), intent(in) :: name    !< name of file
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
  
  !> register a new series of meassurements
  function profile_register(prof,msg)
    use glimmer_log
    implicit none
    type(profile_type) :: prof !< structure storing profile definitions
    character(len=*), intent(in) :: msg !< the message to be associated
    integer profile_register

    prof%nump = prof%nump+1
    if (prof%nump .gt. max_prof) then
       call write_log('Maximum number of profiles reached',type=GM_FATAL, &
            file=__FILE__,line=__LINE__)
    end if
    profile_register = prof%nump
    prof%pname(prof%nump) = trim(msg)
  end function profile_register

  !> start profiling
  subroutine profile_start(prof,profn)
    implicit none
    type(profile_type) :: prof !< structure storing profile definitions
    integer, intent(in) :: profn !< the profile ID
    
    call cpu_time(prof%pstart(profn))
  end subroutine profile_start

  !> stop profiling
  subroutine profile_stop(prof,profn)
    implicit none
    type(profile_type)  :: prof !< structure storing profile definitions
    integer, intent(in) :: profn !< the profile ID
    
    real t
    call cpu_time(t)
    prof%ptotal(profn) = prof%ptotal(profn) + t-prof%pstart(profn)
  end subroutine profile_stop

  !> log a message to profile
  subroutine profile_log(prof,profn,msg)
    implicit none
    type(profile_type)           :: prof !< structure storing profile definitions
    integer, intent(in) :: profn !< the profile ID
    character(len=*), intent(in), optional :: msg     !< message to be written to profile

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

  !> close profile
  subroutine profile_close(prof)
    implicit none
    type(profile_type), intent(in) :: prof !< structure storing profile definitions
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
