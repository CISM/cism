! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_integrate.f90 - part of the GLIMMER ice model    + 
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

module glimmer_integrate
  !*FD integration of functions

  use glimmer_global, only : dp,sp

  private :: dp, sp

  interface romberg_int
     module procedure sromberg_int, sromberg_int_prms, dromberg_int, dromberg_int_prms
  end interface

contains

  recursive real(sp) function sromberg_int(fct,lgr,rgr)

    !*FD Function to perform Romberg Integration on function \texttt{fct}, between
    !*FD limits \texttt{lgr} and \texttt{rgr}. The precision of the routine is 
    !*FD determined by the value of \texttt{ord}, an internal variable. 
    !*FD
    !*FD This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
    !*FD (Comm. ACM, vol. 4, issue 6, June 1961).

    implicit none

    interface 
       function fct(x)
         use glimmer_global, only : sp
         implicit none
         real(sp), intent(in) :: x
         real(sp) :: fct
       end function fct
    end interface
    
    real(sp),intent(in) :: lgr    !*FD Lower bound
    real(sp),intent(in) :: rgr    !*FD Upper bound
    integer,parameter :: ord = 6

    real(sp),dimension(ord+1) :: t
    real(sp) :: l,u,m
    integer :: f,h,j,n

    l=rgr-lgr
    t(1)=(fct(lgr)+fct(rgr))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    sromberg_int=t(1)*l

  end function sromberg_int


  recursive real(sp) function sromberg_int_prms(fct,lgr,rgr,params)

    !*FD Function to perform Romberg Integration on function \texttt{fct}, between
    !*FD limits \texttt{lgr} and \texttt{rgr}. \texttt{params} is an array 
    !*FD passed to the function containing parameters. The precision of the routine is 
    !*FD determined by the value of \texttt{ord}, an internal variable. 
    !*FD
    !*FD This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
    !*FD (Comm. ACM, vol. 4, issue 6, June 1961).

    implicit none

    interface 
       function fct(x,params)
         use glimmer_global, only : sp
         implicit none
         real(sp), intent(in) :: x
         real(sp), intent(in), dimension(:) :: params
         real(sp) :: fct
       end function fct
    end interface
    
    real(sp),intent(in) :: lgr    !*FD Lower bound
    real(sp),intent(in) :: rgr    !*FD Upper bound
    real(sp),intent(in),dimension(:) :: params !*FD parameters for function
    integer,parameter :: ord = 6

    real(sp),dimension(ord+1) :: t
    real(sp) :: l,u,m
    integer :: f,h,j,n

    l=rgr-lgr
    t(1)=(fct(lgr,params)+fct(rgr,params))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m,params)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    sromberg_int_prms=t(1)*l

  end function sromberg_int_prms

  recursive real(dp) function dromberg_int(fct,lgr,rgr)

    !*FD Function to perform Romberg Integration on function \texttt{fct}, between
    !*FD limits \texttt{lgr} and \texttt{rgr}. The precision of the routine is 
    !*FD determined by the value of \texttt{ord}, an internal variable. 
    !*FD
    !*FD This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
    !*FD (Comm. ACM, vol. 4, issue 6, June 1961).

    implicit none

    interface 
       function fct(x)
         use glimmer_global, only : dp
         implicit none
         real(dp), intent(in) :: x
         real(dp) :: fct
       end function fct
    end interface
    
    real(dp),intent(in) :: lgr    !*FD Lower bound
    real(dp),intent(in) :: rgr    !*FD Upper bound
    integer,parameter :: ord = 6

    real(dp),dimension(ord+1) :: t
    real(dp) :: l,u,m
    integer :: f,h,j,n

    l=rgr-lgr
    t(1)=(fct(lgr)+fct(rgr))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    dromberg_int=t(1)*l

  end function dromberg_int


  recursive real(dp) function dromberg_int_prms(fct,lgr,rgr,params)

    !*FD Function to perform Romberg Integration on function \texttt{fct}, between
    !*FD limits \texttt{lgr} and \texttt{rgr}. \texttt{params} is an array 
    !*FD passed to the function containing parameters. The precision of the routine is 
    !*FD determined by the value of \texttt{ord}, an internal variable. 
    !*FD
    !*FD This routine is an implementation of ACM algorithm 60, by F. L. Bauer.
    !*FD (Comm. ACM, vol. 4, issue 6, June 1961).

    implicit none

    interface 
       function fct(x,params)
         use glimmer_global, only : dp
         implicit none
         real(dp), intent(in) :: x
         real(dp), intent(in), dimension(:) :: params
         real(dp) :: fct
       end function fct
    end interface
    
    real(dp),intent(in) :: lgr    !*FD Lower bound
    real(dp),intent(in) :: rgr    !*FD Upper bound
    real(dp),intent(in),dimension(:) :: params !*FD parameters for function
    integer,parameter :: ord = 6

    real(dp),dimension(ord+1) :: t
    real(dp) :: l,u,m
    integer :: f,h,j,n

    l=rgr-lgr
    t(1)=(fct(lgr,params)+fct(rgr,params))/2.0
    n=1

    do h=1,ord
       u=0
       m=l/(2*n)

       do j=1,2*n-1,2
          u=u+fct(lgr+j*m,params)
       end do

       t(h+1)=((u/n)+t(h))/2.0
       f=1
       
       do j=h,1,-1
          f=f*4
          t(j)=t(j+1)+(t(j+1)-t(j))/(f-1)
       end do

       n=2*n

    end do

    dromberg_int_prms=t(1)*l

  end function dromberg_int_prms

end module glimmer_integrate
