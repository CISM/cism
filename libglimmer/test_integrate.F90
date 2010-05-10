! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  test_integrate.f90 - part of the GLIMMER ice model       + 
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

module test_int
  use glimmer_global, only : dp, sp
contains
  
  ! integrate $\int_0^\pi\a*sin(x+b)dx$

  real(dp) function dfp(x,p)
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(in), dimension(:) :: p
    
    dfp = p(1)*sin(p(2)*x)
  end function dfp
  
  real(dp) function df(x)
    implicit none
    real(dp), intent(in) :: x

    real(dp), dimension(2) :: p = (/3.d0,0.5d0/)
    
    df = p(1)*sin(p(2)*x)
  end function df

  real(sp) function sfp(x,p)
    implicit none
    real(sp), intent(in) :: x
    real(sp), intent(in), dimension(:) :: p
    
    sfp =  p(1)*sin(p(2)*x)
  end function sfp
  
  real(sp) function sf(x)
    implicit none
    real(sp), intent(in) :: x

    real(sp), dimension(2) :: p = (/3.,0.5/)
    
    sf = p(1)*sin(p(2)*x)
  end function sf
end module test_int

program test_integrate
  !*FD test numerical integration schemes
  use test_int
  use glimmer_physcon, only : pi
  use glimmer_global, only : sp,dp
  use glimmer_integrate
  implicit none

  real(dp), dimension(2) :: dprms = (/3.d0,0.5d0/)
  real(sp), dimension(2) :: sprms = (/3.0,0.50/)

  real(dp) :: dupper = 4.d0*pi
  real(sp) :: supper = 4.0*pi

  write(*,*) 'dp p ',romberg_int(dfp,0.d0,dupper,dprms)
  write(*,*) 'sp p ',romberg_int(sfp,0.,supper,sprms)
  write(*,*) 'dp   ',romberg_int(df,0.d0,dupper)
  write(*,*) 'sp   ',romberg_int(sf,0.,supper)

end program test_integrate



