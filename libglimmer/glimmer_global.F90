! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_global.f90 - part of the Glimmer-CISM ice model  + 
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

!> Module holding global variables for Glimmer. Holds real-type
!! kind values, and other global code parameters.
module glimmer_global


  integer,parameter :: sp = kind(1.0)  !< Single precision --- Fortran single-precision real-type kind value.
  
  ! Note that if the code is being compiled with forced typing (e.g. with 
  ! the -r8 flag), then this parameter may need to be set in agreement with 
  ! that.

  integer,parameter :: dp = kind(1.0d0) !< Double precision --- Fortran double-precision real-type kind value
  
  ! Note that if the code is being compiled with forced typing (e.g. with
  ! the -r8 flag), then this parameter may need to be set in agreement
  ! with that

#ifdef GLIMMER_SP

  integer,parameter :: rk=sp !< Precision of glimmer module --- the general Fortran real-type kind value for the Glimmer module and its interfaces.

#else

  integer,parameter :: rk=dp !< Precision of glimmer module --- the general Fortran real-type kind value for the Glimmer module and its interfaces.

#endif


  integer,parameter :: fname_length=200 !< Specifies the length of character string variables used to hold filenames.
  integer,parameter :: msg_length=500  !< lenght of message buffers

  character, parameter :: dirsep = '/' !< directory separator

  character, parameter :: linefeed = achar(10)          !< ASCII linefeed
  character, parameter :: char_ret = achar(13)          !< ASCII carriage-return
  character(2), parameter :: cr_lf = char_ret//linefeed !< default newline appropriate for UNIX-type systems
  character, parameter :: endline = linefeed

end module glimmer_global
