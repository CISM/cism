! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_smb.f90 - part of the GLIMMER ice model            + 
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

module smb_dummy

  ! This module provides a dummy, hopefully warning-free interface
  ! in place of the Energy-balance mass-balance scheme. If either
  ! subroutine is called, a fatal error is flagged.

  use glimmer_global

  implicit none

  type smb_params
     integer       :: dummyint
     real(rk)      :: dummyreal
     character(40) :: dummypath
  end type smb_params

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_SMB_DUMMY
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_SMB_DUMMY
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_SMB_DUMMY
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_SMB_DUMMY
!MH!#endif

  subroutine SMBInitWrapper(params,nx,ny,dxr,tstep,path)

    use glimmer_log

    type(smb_params) :: params
    integer :: nx,ny,dxr,tstep
    character(*) :: path

    ! Fatal error

    call write_log('Glimmer not compiled with EBMB mass-balance scheme',GM_FATAL, &
         __FILE__,__LINE__)

    ! Need these lines to avoid warnings, though they are never executed

    params%dummyint=nx
    params%dummyint=ny
    params%dummyint=dxr
    params%dummyint=tstep
    params%dummypath=path

  end subroutine SMBInitWrapper

  !---------------------------------------------------------------------------------------------

  subroutine SMBStepWrapper(params,temp,thck,artm,prcp,U10m,V10m,humidity,SWdown,LWdown,Psurf)

    use glimmer_log

    type(smb_params)        :: params
    real(rk),dimension(:,:) :: temp,thck,artm,prcp,U10m,V10m,humidity,SWdown,LWdown,Psurf

    ! Fatal error

    call write_log('Glimmer not compiled with EBMB mass-balance scheme',GM_FATAL, &
         __FILE__,__LINE__)

    ! Need this line to avoid warnings, though it is never executed

    params%dummyreal=sum(temp+thck+artm+prcp+U10m+V10m+humidity+SWdown+LWdown+Psurf)

  end subroutine SMBStepWrapper

end module smb_dummy
