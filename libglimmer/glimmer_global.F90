! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_global.f90 - part of the GLIMMER ice model       + 
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

module glimmer_global

  !*FD Module holding global variables for Glimmer. Holds real-type
  !*FD kind values, and other global code parameters.

  use glimmer_vers

  integer,parameter :: sp = kind(1.0) 
  
  !*FD Single precision --- Fortran single-precision real-type kind 
  !*FD value. Used internally.
  !*FD
  !*FD Note that if the code is being compiled with forced typing (e.g. with 
  !*FD the -r8 flag), then this parameter may need to be set in agreement with 
  !*FD that.

  integer,parameter :: dp = kind(1.0d0) 
  
  !*FD Double precision --- Fortran double-precision real-type kind 
  !*FD value. Used internally.
  !*FD
  !*FD Note that if the code is being compiled with forced typing (e.g. with
  !*FD the -r8 flag), then this parameter may need to be set in agreement
  !*FD with that

#ifdef GLIMMER_SP

  integer,parameter :: rk=sp

#else

  integer,parameter :: rk=dp

#endif

  !*FD Precision of glimmer module --- the general Fortran real-type kind value 
  !*FD for the Glimmer module and its interfaces.
  !*FD
  !*FD Note that if the code is being compiled with forced typing (e.g. with 
  !*FD the -r8 flag), then this parameter must be set in agreement with that. 

  integer,parameter :: fname_length=100

  !*FD Specifies the length of character string variables used to
  !*FD hold filenames.

  character, parameter :: dirsep = '/'
  !*FD directory separator

  character, parameter :: linefeed = achar(10)
  character, parameter :: char_ret = achar(13)
  character(2), parameter :: cr_lf = char_ret//linefeed
  character, parameter :: endline = linefeed
  !*FD ASCII linefeed and carriage-return characters,
  !*FD and set up default newline appropriate for UNIX-type systems

end module glimmer_global
