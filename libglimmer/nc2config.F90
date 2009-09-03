! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  nc2config.f90 - part of the GLIMMER ice model            + 
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

#ifdef HAVE_2003ARGS
#define NARGS   command_argument_count
#define GETARG  get_command_argument
#else
#define NARGS   iargc
#define GETARG  getarg
#endif

program nc2config

  !*FD Program to extract the config file data from a
  !*FD glimmer netcdf output file, and write it out, either
  !*FD to file or the screen. Command-line argument and 
  !*FD option support is provided, where available

  use netcdf
  use glimmer_global, only: endline

  implicit none

  character(100) :: infile, outfile
  logical :: stdout = .true.
  character(10000) :: config
  integer :: status,ncid,attlen,unit,i,ellen
  integer numargs,nfiles
  integer, external :: iargc
  character(len=100) :: argument
  integer, dimension(100) :: argumentIdx
  
  ! get number of arguments
  numargs = NARGS()
  if (numargs.gt.0) then
     i=0
     nfiles = 0
     ! loop over command line arguments
     do while (i.lt.numargs)
        i = i + 1
        call GETARG(i,argument)
        ! check if it is an option
        if (argument(1:1).eq.'-') then
           select case (trim(argument))
           case ('-h')
              call usage()
              stop
           case ('-o')
              i = i+1
              if (i.gt.numargs) then
                 write(*,*) 'Error, expect name of output file to follow -o option'
                 call usage()
                 stop
              end if
              call GETARG(i,outfile)
              stdout = .false.
           case default
              write(*,*) 'Unkown option ',trim(argument)
              call usage()
              stop
           end select
        else
           ! it's not an option
           nfiles = nfiles+1
           argumentIdx(nfiles) = i
        end if
     end do
     if (nfiles.gt.1) then
        call GETARG(argumentIdx(1),infile)
     else
        write(*,*) 'No input file specified'
        call usage()
        stop
     end if
  else
     ! These are the default inputs
     write(*,*) 'Enter name of input file:'
     write(*,'(a)') infile
     write(*,*) 'Enter name of output file:'
     write(*,'(a)') outfile
     stdout=.false.
  end if

  ! Open file and look for appropriate attribute

  status=nf90_open(infile,0,ncid)
  if (status/=NF90_NOERR) call netcdf_error(status)
  status=nf90_inquire_attribute(ncid,NF90_GLOBAL,'configuration',len=attlen)
  if (status==NF90_ENOTATT) then
     write(0,*)'ERROR: No configuration data found in file ',trim(infile)
     stop
  else if (status/=NF90_NOERR) then
     call netcdf_error(status)
  end if
  status=nf90_get_att(ncid,NF90_GLOBAL,'configuration',config)
  if (status/=NF90_NOERR) call netcdf_error(status)
  config=config(:attlen)

  ! Format and print output

  ellen=len(endline)

  if (stdout) then
     unit=6
  else
     unit=20
     open(unit,file=outfile)
  end if

  i=1
  do 
     if (trim(config(i:))=='') exit
     if (config(i:i+ellen-1)==endline) then
        write(unit,'(A)')''
        i=i+ellen
     else if (config(i:i)=='[') then
        write(unit,'(A)')''
        write(unit,'(A)',advance='no')config(i:i)
        i=i+1
     else
        write(unit,'(A)',advance='no')config(i:i)
        i=i+1
     end if
  end do

contains

  subroutine usage

    print*,'usage: nc2config [-h  -o outfile] [infile]'
    print*,'    -h: display this message'
    print*,'    -o: specify output file. If absent, stdout is used'

  end subroutine usage

  subroutine netcdf_error(status)

    integer :: status

    print*,nf90_strerror(status)
    stop

  end subroutine netcdf_error

end program nc2config
