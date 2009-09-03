! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_coordinates.f90 - part of the GLIMMER ice model  + 
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

module glimmer_coordinates

  !*FD module for handling regular coordinate systems
  use glimmer_global, only: dp, sp

  type coord_point
     !*FD  point type
     real(kind=dp), dimension(2) :: pt
  end type coord_point

  type coord_ipoint
     !*FD integer point type
     integer, dimension(2) :: pt
  end type coord_ipoint

  type coordsystem_type
     !*FD type describing coordinate systems
     type(coord_point) :: origin   !*FD origin of coordinate space
     type(coord_point) :: delta    !*FD stepsize in x and y direction
     type(coord_point) :: delta_r  !*FD reciprocal stepsize in x and y direction
     type(coord_ipoint) :: size    !*FD extent in x and y direction
  end type coordsystem_type
  
  ! interface of creating new coord system
  interface coordsystem_new
     module procedure coordsystem_new_real, coordsystem_new_pt
  end interface

  interface coordsystem_allocate
     module procedure coordsystem_allocate_d, coordsystem_allocate_s, coordsystem_allocate_i, coordsystem_allocate_l, &
          coordsystem_allocate_d2, coordsystem_allocate_s2, coordsystem_allocate_i2
  end interface

#ifdef DEBUG_COORDS
  character(len=100), private :: message
#endif
  !NO_RESTART message
  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_COORDINATES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_COORDINATES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_COORDINATES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_COORDINATES
!MH!#endif

  subroutine coordsystem_print(coord, unit)
    !*FD print coordsystem info to unit
    implicit none
    type(coordsystem_type), intent(in) :: coord  !*FD coordinate system
    integer,intent(in) :: unit                   !*FD unit to be printed to
    write(unit,*) 'Origin  ',coord%origin%pt
    write(unit,*) 'Delta   ',coord%delta%pt
    write(unit,*) '1/Delta ',coord%delta_r%pt
    write(unit,*) 'Size    ',coord%size%pt
  end subroutine coordsystem_print

  function coordsystem_new_real(ox, oy, dx, dy, sx, sy)
    !*FD create new coordinate system from individual variables
    implicit none
    real(kind=dp), intent(in) :: ox, oy !*FD coordinates of origin
    real(kind=dp), intent(in) :: dx, dy !*FD offsets
    integer, intent(in) :: sx, sy       !*FD x and y dimension
    type(coordsystem_type) :: coordsystem_new_real
    
    ! origin
    coordsystem_new_real%origin%pt(1) = ox
    coordsystem_new_real%origin%pt(2) = oy
    ! deltas
    coordsystem_new_real%delta%pt(1) = dx
    coordsystem_new_real%delta%pt(2) = dy
    coordsystem_new_real%delta_r%pt(1) = 1.d0/dx
    coordsystem_new_real%delta_r%pt(2) = 1.d0/dy
    ! size
    coordsystem_new_real%size%pt(1) = sx
    coordsystem_new_real%size%pt(2) = sy
  end function coordsystem_new_real

  function coordsystem_new_pt(o, d, s)
    !*FD create new coordinate system from points
    implicit none
    type(coord_point), intent(in) :: o  !*FD coordinates of origin
    type(coord_point), intent(in) :: d  !*FD offsets
    type(coord_ipoint), intent(in) :: s !*FD x and y dimension
    type(coordsystem_type) :: coordsystem_new_pt

    ! origin
    coordsystem_new_pt%origin = o
    ! deltas
    coordsystem_new_pt%delta = d
    coordsystem_new_pt%delta_r%pt(:) = 1.d0/d%pt(:)
    ! size
    coordsystem_new_pt%size = s
  end function coordsystem_new_pt

  function coordsystem_get_coord(coord,node)
    !*FD get coordinates of node
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord  !*FD coordinate system
    type(coord_ipoint), intent(in) :: node       !*FD node

    type(coord_point) :: coordsystem_get_coord
  
#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,node)) then
       write(message,*) 'node (',node%pt,') not inside coord system'
       call coordsystem_print(coord,glimmer_get_logunit())
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    coordsystem_get_coord%pt(:) = coord%origin%pt(:) + (node%pt(:) - 1)*coord%delta%pt(:)
  end function coordsystem_get_coord

  function coordsystem_get_node(coord,point)
    !*FD get index of nearest node given coords of a point
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    type(coord_point), intent(in) :: point      !*FD point
    
    type(coord_ipoint) :: coordsystem_get_node
    
    coordsystem_get_node%pt(:) = 1+floor(0.5+(point%pt(:)-coord%origin%pt(:))*coord%delta_r%pt(:))
    if (coordsystem_get_node%pt(1).eq.coord%size%pt(1)+1) coordsystem_get_node%pt(1) = coord%size%pt(1)
    if (coordsystem_get_node%pt(2).eq.coord%size%pt(2)+1) coordsystem_get_node%pt(2) = coord%size%pt(2)

#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,coordsystem_get_node)) then
       write(message,*) 'point (',point%pt,') not inside coord system'
       call coordsystem_print(coord,glimmer_get_logunit())
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif
  end function coordsystem_get_node

  function coordsystem_get_llnode(coord,point)
    !*FD get index of lower-left node of cell into which point falls
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    type(coord_point), intent(in) :: point      !*FD point
    
    type(coord_ipoint) :: coordsystem_get_llnode

    coordsystem_get_llnode%pt(:) = 1+floor((point%pt(:)-coord%origin%pt(:))*coord%delta_r%pt(:))
  end function coordsystem_get_llnode

  function coordsystem_node_inside(coord,node)
    !*FD return true iff node is inside coord system
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    type(coord_ipoint), intent(in) :: node      !*FD node

    logical coordsystem_node_inside
    
    coordsystem_node_inside = (all(node%pt.ge.1) .and. all(node%pt.le.coord%size%pt))
  end function coordsystem_node_inside

  function coordsystem_point_inside(coord,point)
    !*FD return true iff point is inside coord system
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    type(coord_point), intent(in) :: point      !*FD point
    logical coordsystem_point_inside
    integer i

    coordsystem_point_inside = .true.
    do i=1,2
       coordsystem_point_inside = (point%pt(i).ge.coord%origin%pt(i)) .and. &
            (point%pt(i).le.coord%origin%pt(i)+coord%size%pt(i)*coord%delta%pt(i))
       if (.not.coordsystem_point_inside) then
          exit
       end if
    end do
  end function coordsystem_point_inside
    
  function coordsystem_linearise2d(coord,node)
    !*FD linearise node, given coord
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    type(coord_ipoint), intent(in) :: node      !*FD node
    integer coordsystem_linearise2d

    coordsystem_linearise2d = -1

#ifdef DEBUG_COORDS
    if (.not.coordsystem_node_inside(coord,node)) then
       write(message,*) 'node (',node%pt,') not inside coord system'
       call write_log(message,GM_ERROR,__FILE__,__LINE__)
       return
    end if
#endif
    
    coordsystem_linearise2d = node%pt(1) + (node%pt(2)-1)*coord%size%pt(1)
  end function coordsystem_linearise2d

  function coordsystem_delinearise2d(coord, ind)
    !*FD expand linearisation
    use glimmer_log
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    integer, intent(in) :: ind                  !*FD index
    type(coord_ipoint) :: coordsystem_delinearise2d

#ifdef DEBUG_COORDS
    if (ind.lt.1 .or. ind.gt.coord%size%pt(1)*coord%size%pt(2)) then
       write(message,*) 'index ',ind,' outside coord system'
       call write_log(message,GM_FATAL,__FILE__,__LINE__)
    end if
#endif

    coordsystem_delinearise2d%pt(1) = mod(ind-1,coord%size%pt(1)) + 1
    coordsystem_delinearise2d%pt(2) = (ind-1)/coord%size%pt(1) + 1
  end function coordsystem_delinearise2d

  subroutine coordsystem_allocate_d(coord, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    real(kind=dp), dimension(:,:), pointer :: field !*FD unallocated field

    allocate(field(coord%size%pt(1),coord%size%pt(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_d
  
  subroutine coordsystem_allocate_s(coord, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    real(kind=sp), dimension(:,:), pointer :: field !*FD unallocated field

    allocate(field(coord%size%pt(1),coord%size%pt(2)))
    field = 0.e0
  end subroutine coordsystem_allocate_s

  subroutine coordsystem_allocate_i(coord, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    integer, dimension(:,:), pointer :: field !*FD unallocated field

    allocate(field(coord%size%pt(1),coord%size%pt(2)))
    field = 0
  end subroutine coordsystem_allocate_i

  subroutine coordsystem_allocate_l(coord, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    logical, dimension(:,:), pointer :: field !*FD unallocated field

    allocate(field(coord%size%pt(1),coord%size%pt(2)))
    field = .FALSE.
  end subroutine coordsystem_allocate_l

  subroutine coordsystem_allocate_d2(coord, nup, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    integer, intent(in) :: nup
    real(kind=dp), dimension(:,:,:), pointer :: field !*FD unallocated field

    allocate(field(nup,coord%size%pt(1),coord%size%pt(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_d2

  subroutine coordsystem_allocate_s2(coord, nup, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    integer, intent(in) :: nup
    real(kind=sp), dimension(:,:,:), pointer :: field !*FD unallocated field

    allocate(field(nup,coord%size%pt(1),coord%size%pt(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_s2

  subroutine coordsystem_allocate_i2(coord, nup, field)
    !*FD allocate memory to pointer field
    implicit none
    type(coordsystem_type), intent(in) :: coord !*FD coordinate system
    integer, intent(in) :: nup
    integer, dimension(:,:,:), pointer :: field !*FD unallocated field

    allocate(field(nup,coord%size%pt(1),coord%size%pt(2)))
    field = 0.d0
  end subroutine coordsystem_allocate_i2

end module glimmer_coordinates
