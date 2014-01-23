!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_higher.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains routines for computing the ice velocity
!  using a variational finite-element approach.
!
! See these papers for details:
!
! J.K. Dukowicz, S.F. Price and W.H. Lipscomb, 2010: Consistent
!    approximations and boundary conditions for ice-sheet dynamics
!    using a principle of least action.  J. Glaciology, 56 (197),
!    480-495.
!
! M. Perego, M. Gunzburger, and J. Burkardt, 2012: Parallel
!    finite-element implementation for higher-order ice-sheet models.
!    J. Glaciology, 58 (207), 76-88.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

  module glissade_velo_higher

    use glimmer_global, only: dp
    use glimmer_physcon, only: gn, rhoi, rhoo, grav, scyr, pi
    use glimmer_paramets, only: thk0, len0, tim0, tau0, vel0, vis0, evs0
    use glimmer_paramets, only: vel_scale, len_scale   ! used for whichefvs = HO_EFVS_FLOWFACT
    use glimmer_log, only: write_log
    use glimmer_sparse_type
    use glimmer_sparse
     
    use glide_types  ! for HO_EFVS and other options

    use cism_sparse_pcg, only: pcg_solver_structured

    use parallel

    implicit none

    private
    !IK, 9/18/13: added staggered_scalar here so it can be called from
    !felix_dycore_interface.F90 
    public :: glissade_velo_higher_init, glissade_velo_higher_solve, staggered_scalar

    !----------------------------------------------------------------
    ! Here are some definitions:
    !
    ! The horizontal mesh is composed of cells and vertices.
    ! All cells are assumed to be quadrilaterals, but the code can be
    !  generalized later to triangles (e.g., for MPAS mesh).
    ! Each cell can be extruded to form a column with a specified number of layers.
    ! 
    ! An element is a layer of a cell, and a node is a corner of an element.
    ! Elements and nodes live in 3D space, whereas cells and vertices live in
    !  the horizontal plane.
    !
    ! Locally owned cells and vertices have indices (nhalo+1:nx-nhalo, nhalo+1,ny-nhalo).
    ! Active cells are cells that (1) contain ice and (2) border locally owned vertices.
    ! Active vertices are all vertices of active cells.
    ! Active nodes are all nodes in the columns associated with active vertices.
    !
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Finite element properties
    !
    ! Assume 3D hexahedral elements.
    !----------------------------------------------------------------

    integer, parameter ::        &
       nNodesPerElement = 8,     & ! 8 nodes for hexahedral elements
       nQuadPoints = 8,          & ! number of quadrature points per element
                                   ! These live at +- 1/sqrt(3) for reference hexahedron
       nNodesPerElement_2d = 4,  & ! 4 nodes for faces of hexahedral elements
       nQuadPoints_2d = 4          ! number of quadrature points per element face
                                   ! These live at +- 1/sqrt(3) for reference square

    real(dp), parameter ::     &
       rsqrt3 = 1.d0/sqrt(3.d0)    ! for quadrature points
         
    !----------------------------------------------------------------
    ! Arrays used for finite-element calculations
    !
    ! Most integals are done over 3D hexahedral elements.
    ! Surface integrals are done over 2D faces of these elements. 
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement, nQuadPoints) ::   & 
       phi,            &    ! trilinear basis function, evaluated at quad pts
       dphi_dxr,       &    ! dphi/dx for reference element, evaluated at quad pts
       dphi_dyr,       &    ! dphi/dy for reference element, evaluated at quad pts
       dphi_dzr             ! dphi/dy for reference element, evaluated at quad pts

    real(dp), dimension(nQuadPoints) ::  &
       xqp, yqp, zqp,  &    ! quad pt coordinates in reference element
       wqp                  ! quad pt weights

    real(dp), dimension(nNodesPerElement_2d, nQuadPoints_2d) ::   & 
       phi_2d,         &    ! bilinear basis function, evaluated at quad pts
       dphi_dxr_2d,    &    ! dphi/dx for reference square, evaluated at quad pts
       dphi_dyr_2d          ! dphi/dy for reference square, evaluated at quad pts

    real(dp), dimension(nQuadPoints_2d) ::  &
       xqp_2d, yqp_2d, &    ! quad pt coordinates in reference square
       wqp_2d               ! quad pt weights

    integer, dimension(nNodesPerElement, nNodesPerElement) ::  &
       ishift, jshift, kshift   ! matrices describing relative indices of nodes in an element

    integer, dimension(-1:1,-1:1,-1:1) :: &
       indxA                 ! maps relative (x,y,z) coordinates to an index between 1 and 27
                             ! index order is (i,j,k)

    real(dp), dimension(3,3) ::  &
       identity3                ! 3 x 3 identity matrix

    real(dp), parameter ::   &
       eps10 = 1.d-10           ! small number

!    real(dp), parameter ::   &
!       eps09 = 1.d-09           ! small number

    real(dp) :: vol0    ! volume scale = dx * dy * (1000 m)

    logical, parameter ::  &
       check_symmetry = .true.   ! if true, then check symmetry of assembled matrix

    logical, parameter ::  &
       remove_small_values = .false.   ! if true, then remove small values from assembled matrix
                                       ! (resulting from taking the difference of two large terms)
    logical, parameter ::  &
       zero_uvel = .false.,  &         ! if true, zero out uvel everywhere
       zero_vvel = .false.             ! if true, zero out vvel everywhere

!WHL - debug
!    logical :: verbose = .false.      ! for debug print statements
    logical :: verbose = .true.  
    logical :: verbose_init = .false.   
!    logical :: verbose_init = .true.   
    logical :: verbose_Jac = .false.
!    logical :: verbose_Jac = .true.
    logical :: verbose_residual = .false.
!    logical :: verbose_residual = .true.
!    logical :: verbose_state = .false.
    logical :: verbose_state = .true.
    logical :: verbose_load = .false.
!    logical :: verbose_load = .true.
!    logical :: verbose_shelf = .false.
    logical :: verbose_shelf = .true.
    logical :: verbose_matrix = .false.
!    logical :: verbose_matrix = .true.
    logical :: verbose_basal = .false.
!    logical :: verbose_basal = .true.
    logical :: verbose_umc = .false.
!    logical :: verbose_umc = .true.
    logical :: verbose_slapsolve = .false.
!    logical :: verbose_slapsolve = .true.

    logical :: verbose_efvs = .false.
!    logical :: verbose_efvs = .true.

!WHL - debug
    logical :: trial_efvs = .false.   ! if true, compute what nonlinear efvs would be (if not constant)

    integer :: itest, jtest    ! coordinates of diagnostic point
    integer :: rtest           ! task number for processor containing diagnostic point

    integer, parameter :: ktest = 1     ! vertical level of diagnostic point
    integer, parameter :: ptest = 1     ! diagnostic quadrature point

    !TODO - Compute ntest on the fly.
    integer, parameter :: ntest = 1

!    integer, parameter :: rtest = 0    ! rank for any single-process run
!    integer, parameter :: rtest = 1

!    integer, parameter :: &
!       itest = 9, jtest = 19, ktest = 1, ptest = 1  ! for dome, global (i,j) = (7,17), 1 proc
!       itest = 26, jtest = 19, ktest = 1, ptest = 1  ! for dome, global (i,j) = (24,17), 1 proc
!       integer, parameter :: ntest = 2371  ! nodeID for dome global (24,17,1)    

!        itest = 8, jtest = 8, ktest = 1, ptest = 1    ! for block test, global (i,j) = (6,6)

!        itest = 3, jtest = 3, ktest = 1, ptest = 1    ! for ishom, global (i,j) = (1,1)
!        itest = 7, jtest = 7, ktest = 1, ptest = 1    ! for ishom, global (i,j) = (5,5)
!        itest = 6, jtest = 6, ktest = 1, ptest = 1    ! for ishom, global (i,j) = (4,4)


!       itest = 24, jtest = 6, ktest = 1, ptest = 1  ! for confined/linear (south-flowing) shelf, global (i,j) = (22,4), 1 proc
!!       itest = 24, jtest = 42, ktest = 1, ptest = 1  ! for confined/linear (north-flowing) shelf, global (i,j) = (22,40), 1 proc
!       itest = 11, jtest = 11, ktest = 1, ptest = 1  ! for circular shelf, global (i,j) = (9,9), 1 proc
!       itest = 5, jtest = 42, ktest = 1, ptest = 1  ! for confined/linear (north-flowing) shelf, global (i,j) = (3,40), 1 proc

!       itest = 9, jtest = 9, ktest = 1, ptest = 1  ! for confined/linear shelf, global (i,j) = (7,7), 1 proc


!       integer, parameter :: ntest = 73    ! for 5x5 test, northeast corner node at upper surface
!       integer, parameter :: ntest = 48    ! for 4x4 test, northeast corner node at bed

!    integer, parameter :: ntest = 2372  ! nodeID for dome global (24,17,2)
!    integer, parameter :: ntest = 2380  ! nodeID for dome global (24,17,10)
    
!    integer, parameter :: ntest = 1    ! for ishom.a, global (i,j) = (1,1)
!    integer, parameter :: ntest = 882    ! for ishom.a.80km symmetry check

!    integer, parameter :: ntest = 101  ! nodeID for confined (south-flowing) shelf, global (22,4,1)    
!    integer, parameter :: ntest = 7701  ! nodeID for confined shelf global (22,40,1)    


!WHL = debug
    integer, parameter :: rowtest = -999
    integer, parameter :: coltest = -999

    integer, parameter :: &
      iAtest=1, jAtest=0, kAtest=0

    integer, parameter :: &
      Krowtest=1, Kcoltest=2

    logical, parameter :: write_matrix = .false.
!    logical, parameter :: write_matrix = .true.
!    character(*), parameter :: matrix_label = 'ishomC_block'   ! Change depending on the case we're running
!    character(*), parameter :: matrix_label = 'block'   ! Change depending on the case we're running
!    character(*), parameter :: matrix_label = 'ishomC_periodic'   ! Change depending on the case we're running
    character(*), parameter :: matrix_label = 'shelf'   ! Change depending on the case we're running

    contains

!****************************************************************************

  subroutine glissade_velo_higher_init

    !----------------------------------------------------------------
    ! Initial calculations for glissade higher-order solver.
    !----------------------------------------------------------------

    integer :: i, j, k, m, n, p

!WHL - debug
    real(dp) :: sumx, sumy, sumz

    !----------------------------------------------------------------
    ! Initialize some time-independent finite element arrays
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Trilinear basis set for reference hexahedron, x=(-1,1), y=(-1,1), z=(-1,1)             
    ! Indexing is counter-clockwise from SW corner, with 1-4 on lower surface
    !  and 5-8 on upper surface
    ! In the code we use "phi" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)*(1-z)/8             N4----N3
    ! N2 = (1+x)*(1-y)*(1-z)/8             |     |    Lower layer        
    ! N3 = (1+x)*(1+y)*(1-z)/8             |     |
    ! N4 = (1-x)*(1+y)*(1-z)/8             N1----N2

    ! N5 = (1-x)*(1-y)*(1+z)/8             N8----N7
    ! N6 = (1+x)*(1-y)*(1+z)/8             |     |    Upper layer
    ! N7 = (1+x)*(1+y)*(1+z)/8             |     |
    ! N8 = (1-x)*(1+y)*(1+z)/8             N5----N6
    !----------------------------------------------------------------
   
    ! Set coordinates and weights of quadrature points for reference hexahedral element.
    ! Numbering is counter-clockwise from southwest, lower face (1-4) followed by
    !  upper face (5-8).

    xqp(1) = -rsqrt3; yqp(1) = -rsqrt3; zqp(1) = -rsqrt3
    wqp(1) =  1.d0

    xqp(2) =  rsqrt3; yqp(2) = -rsqrt3; zqp(2) = -rsqrt3
    wqp(2) =  1.d0

    xqp(3) =  rsqrt3; yqp(3) =  rsqrt3; zqp(3) = -rsqrt3
    wqp(3) =  1.d0

    xqp(4) = -rsqrt3; yqp(4) =  rsqrt3; zqp(4) = -rsqrt3
    wqp(4) =  1.d0

    xqp(5) = -rsqrt3; yqp(5) = -rsqrt3; zqp(5) =  rsqrt3
    wqp(5) =  1.d0

    xqp(6) =  rsqrt3; yqp(6) = -rsqrt3; zqp(6) =  rsqrt3
    wqp(6) =  1.d0

    xqp(7) =  rsqrt3; yqp(7) =  rsqrt3; zqp(7) =  rsqrt3
    wqp(7) =  1.d0

    xqp(8) = -rsqrt3; yqp(8) =  rsqrt3; zqp(8) =  rsqrt3
    wqp(8) =  1.d0

    if (verbose_init) then
       print*, ' '
       print*, 'Hexahedral elements, quad points, x, y, z:'
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints
          print*, p, xqp(p), yqp(p), zqp(p)
          sumx = sumx + xqp(p); sumy = sumy + yqp(p); sumz = sumz + zqp(p)
       enddo
       print*, ' '
       print*, 'sums:', sumx, sumy, sumz
    endif

    ! Evaluate trilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints

       phi(1,p) = (1.d0 - xqp(p)) * (1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(2,p) = (1.d0 + xqp(p)) * (1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(3,p) = (1.d0 + xqp(p)) * (1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(4,p) = (1.d0 - xqp(p)) * (1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       phi(5,p) = (1.d0 - xqp(p)) * (1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0
       phi(6,p) = (1.d0 + xqp(p)) * (1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0
       phi(7,p) = (1.d0 + xqp(p)) * (1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0
       phi(8,p) = (1.d0 - xqp(p)) * (1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0

       dphi_dxr(1,p) = -(1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dxr(2,p) =  (1.d0 - yqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dxr(3,p) =  (1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dxr(4,p) = -(1.d0 + yqp(p)) * (1.d0 - zqp(p)) / 8.d0
       dphi_dxr(5,p) = -(1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dxr(6,p) =  (1.d0 - yqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dxr(7,p) =  (1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dxr(8,p) = -(1.d0 + yqp(p)) * (1.d0 + zqp(p)) / 8.d0

       dphi_dyr(1,p) = -(1.d0 - xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(2,p) = -(1.d0 + xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(3,p) =  (1.d0 + xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(4,p) =  (1.d0 - xqp(p)) * (1.d0 - zqp(p)) / 8.d0 
       dphi_dyr(5,p) = -(1.d0 - xqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dyr(6,p) = -(1.d0 + xqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dyr(7,p) =  (1.d0 + xqp(p)) * (1.d0 + zqp(p)) / 8.d0 
       dphi_dyr(8,p) =  (1.d0 - xqp(p)) * (1.d0 + zqp(p)) / 8.d0 

       dphi_dzr(1,p) = -(1.d0 - xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(2,p) = -(1.d0 + xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(3,p) = -(1.d0 + xqp(p)) * (1.d0 + yqp(p)) / 8.d0 
       dphi_dzr(4,p) = -(1.d0 - xqp(p)) * (1.d0 + yqp(p)) / 8.d0 
       dphi_dzr(5,p) =  (1.d0 - xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(6,p) =  (1.d0 + xqp(p)) * (1.d0 - yqp(p)) / 8.d0 
       dphi_dzr(7,p) =  (1.d0 + xqp(p)) * (1.d0 + yqp(p)) / 8.d0 
       dphi_dzr(8,p) =  (1.d0 - xqp(p)) * (1.d0 + yqp(p)) / 8.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          print*, 'n, phi, dphi_dxr, dphi_dyr, dphi_dzr:'
          do n = 1, 8
             print*, n, phi(n,p), dphi_dxr(n,p), dphi_dyr(n,p), dphi_dzr(n,p)
          enddo
          print*, ' '
          print*, 'sum(phi)', sum(phi(:,p))  ! verified that sum = 1
          print*, 'sum(dphi/dx)', sum(dphi_dxr(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy)', sum(dphi_dyr(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dz)', sum(dphi_dzr(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints

    ! Identity matrix
    identity3(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
    identity3(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
    identity3(3,:) = (/ 0.d0, 0.d0, 1.d0 /)

    ! Initialize some matrices that describe how the i, j and k indices of each node
    ! in each element are related to one another.

    ! The ishift matrix describes how the i indices of the 8 nodes are related to one another.
    ! E.g, if ishift (1,2) = 1, this means that node 2 has an i index
    ! one greater than the i index of node 1.

    ishift(1,:) = (/ 0,  1,  1,  0,  0,  1,  1,  0/)   
    ishift(2,:) = (/-1,  0,  0, -1, -1,  0,  0, -1/)   
    ishift(3,:) = ishift(2,:)
    ishift(4,:) = ishift(1,:)
    ishift(5,:) = ishift(1,:)
    ishift(6,:) = ishift(2,:)
    ishift(7,:) = ishift(2,:)
    ishift(8,:) = ishift(1,:)

    ! The jshift matrix describes how the j indices of the 8 nodes are related to one another.
    ! E.g, if jshift (1,4) = 1, this means that node 4 has a j index
    ! one greater than the j index of node 1.

    jshift(1,:) = (/ 0,  0,  1,  1,  0,  0,  1,  1/)   
    jshift(2,:) = jshift(1,:)
    jshift(3,:) = (/-1, -1,  0,  0, -1, -1,  0,  0/)   
    jshift(4,:) = jshift(3,:)
    jshift(5,:) = jshift(1,:)
    jshift(6,:) = jshift(1,:)
    jshift(7,:) = jshift(3,:)
    jshift(8,:) = jshift(3,:)

    ! The kshift matrix describes how the k indices of the 8 nodes are related to one another.
    ! E.g, if kshift (1,5) = -1, this means that node 5 has a k index
    ! one less than the k index of node 1.  (Assume that k increases downward.)

    kshift(1,:) = (/ 0,  0,  0,  0, -1, -1, -1, -1/)   
    kshift(2,:) = kshift(1,:)
    kshift(3,:) = kshift(1,:)
    kshift(4,:) = kshift(1,:)
    kshift(5,:) = (/ 1,  1,  1,  1,  0,  0,  0,  0/)
    kshift(6,:) = kshift(5,:)
    kshift(7,:) = kshift(5,:)
    kshift(8,:) = kshift(5,:)

    if (verbose_init) then
       print*, ' '
       print*, 'ishift:'
       do n = 1, 8
          write (6,'(8i4)') ishift(n,:)
       enddo
       print*, ' '
       print*, 'jshift:'
       do n = 1, 8
          write (6,'(8i4)') jshift(n,:)
       enddo
       print*, ' '
       print*, 'kshift:'
       do n = 1, 8
          write (6,'(8i4)') kshift(n,:)
       enddo
    endif

    !----------------------------------------------------------------
    ! Bilinear basis set for reference square, x=(-1,1), y=(-1,1)             
    ! Indexing is counter-clockwise from SW corner
    ! In the code we use "phi_2d" to denote these basis functions. 
    !
    ! N1 = (1-x)*(1-y)/4             N4----N3
    ! N2 = (1+x)*(1-y)/4             |     |
    ! N3 = (1+x)*(1+y)/4             |     |
    ! N4 = (1-x)*(1+y)/4             N1----N2
    !----------------------------------------------------------------

    ! Set coordinates and weights of quadrature points for reference square.
    ! Numbering is counter-clockwise from southwest

    xqp_2d(1) = -rsqrt3; yqp_2d(1) = -rsqrt3
    wqp_2d(1) =  1.d0

    xqp_2d(2) =  rsqrt3; yqp_2d(2) = -rsqrt3
    wqp_2d(2) =  1.d0

    xqp_2d(3) =  rsqrt3; yqp_2d(3) =  rsqrt3
    wqp_2d(3) =  1.d0

    xqp_2d(4) = -rsqrt3; yqp_2d(4) =  rsqrt3
    wqp_2d(4) =  1.d0

    if (verbose_init) then
       print*, ' '
       print*, ' '
       print*, 'Quadrilateral elements, quad points, x, y:'
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints_2d
          print*, p, xqp_2d(p), yqp_2d(p)
          sumx = sumx + xqp_2d(p); sumy = sumy + yqp_2d(p)
       enddo
       print*, ' '
       print*, 'sumx, sumy:', sumx, sumy
    endif

    ! Evaluate bilinear basis functions and their derivatives at each quad pt

    do p = 1, nQuadPoints_2d

       phi_2d(1,p) = (1.d0 - xqp_2d(p)) * (1.d0 - yqp_2d(p)) / 4.d0 
       phi_2d(2,p) = (1.d0 + xqp_2d(p)) * (1.d0 - yqp_2d(p)) / 4.d0
       phi_2d(3,p) = (1.d0 + xqp_2d(p)) * (1.d0 + yqp_2d(p)) / 4.d0 
       phi_2d(4,p) = (1.d0 - xqp_2d(p)) * (1.d0 + yqp_2d(p)) / 4.d0

       dphi_dxr_2d(1,p) = -(1.d0 - yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(2,p) =  (1.d0 - yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(3,p) =  (1.d0 + yqp_2d(p)) / 4.d0 
       dphi_dxr_2d(4,p) = -(1.d0 + yqp_2d(p)) / 4.d0

       dphi_dyr_2d(1,p) = -(1.d0 - xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(2,p) = -(1.d0 + xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(3,p) =  (1.d0 + xqp_2d(p)) / 4.d0 
       dphi_dyr_2d(4,p) =  (1.d0 - xqp_2d(p)) / 4.d0 

       if (verbose_init) then
          print*, ' '
          print*, 'Quad point, p =', p
          print*, 'n, phi_2d, dphi_dxr_2d, dphi_dyr_2d:'
          do n = 1, 4
             print*, n, phi_2d(n,p), dphi_dxr_2d(n,p), dphi_dyr_2d(n,p)
          enddo
          print*, 'sum(phi_2d)', sum(phi_2d(:,p))        ! verified that sum = 1
          print*, 'sum(dphi/dx_2d)', sum(dphi_dxr_2d(:,p))  ! verified that sum = 0 (within roundoff)
          print*, 'sum(dphi/dy_2d)', sum(dphi_dyr_2d(:,p))  ! verified that sum = 0 (within roundoff)
       endif

    enddo   ! nQuadPoints_2d

    ! Compute indxA; maps displacements i,j,k = (-1,0,1) onto an index from 1 to 27
    ! Numbering starts in SW corner of layers k-1, finishes in NE corner of layer k+1
    ! Diagonal term has index 14

    ! Layer k-1:           Layer k:            Layer k+1:
    !
    !   7    8    9          16   17   18        25   26   27 
    !   4    5    6          13   14   15        22   23   24
    !   1    2    3          10   11   12        19   20   21                                                                                               

    m = 0
    do k = -1,1
       do j = -1,1
          do i = -1,1
             m = m + 1
             indxA(i,j,k) = m
          enddo
       enddo
    enddo

  end subroutine glissade_velo_higher_init

!****************************************************************************

  subroutine glissade_velo_higher_solve(model,                &
                                        nx,     ny,     nz)

!TODO - Remove nx, ny, nz from argument list?
!       Would then have to allocate many local arrays.

!TODO - Something like this may be needed if building with Trilinos.
!!    use glam_strs2, only: linearSolveTime, totalLinearSolveTime

    use parallel, only: nhalo

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    !----------------------------------------------------------------
    ! Note that the glissade solver uses SI units.
    ! Thus we have grid cell dimensions and ice thickness in meters,
    !  velocity in m/s, and the rate factor in Pa^(-n) s(-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Note: nx and ny are the horizontal dimensions of scalar arrays (e.g., thck and temp)
    !       The velocity arrays have horizontal dimensions (nx-1, ny-1)
    !       nz is the number of levels at which uvel and vvel are computed
    !       The scalar variables generally live at layer midpoints and have
    !         vertical dimension nz-1.
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each direction
       nz                       ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)
 
    !----------------------------------------------------------------
    ! Local variables and pointers set to components of model derived type 
    !----------------------------------------------------------------

!    integer, pointer ::   &
!       nx, ny,               &  ! number of grid cells in each direction
!       nz                       ! number of vertical levels where velocity is computed
                                 ! (same as model%general%upn)
    real(dp) ::  &
       dx,  dy                  ! grid cell length and width (m)
                                ! assumed to have the same value for each grid cell

    real(dp), dimension(:), pointer :: &
       sigma                    ! vertical sigma coordinate, [0,1]

    real(dp)  ::   & 
       thklim, &                ! minimum ice thickness for active cells (m)
       eus                      ! eustatic sea level (m)
                                ! = 0. by default

    real(dp), dimension(:,:), pointer ::  &
       thck,                 &  ! ice thickness (m)
       usrf,                 &  ! upper surface elevation (m)
       topg,                 &  ! elevation of topography (m)
       beta                     ! basal traction parameter (Pa/(m/yr))

    real(dp), dimension(:,:,:), pointer ::  &
       flwa,   &                ! flow factor in units of Pa^(-n) yr^(-1)
       efvs,   &                ! effective viscosity (Pa yr)
       uvel, vvel               ! velocity components (m/yr)

    integer,  dimension(:,:), pointer ::   &
       kinbcmask                ! = 1 at vertices where u = v = 0 (Dirichlet BC)
                                ! = 0 elsewhere

    integer ::   &
       whichefvs, &             ! option for effective viscosity calculation 
                                ! (calculate it or make it uniform)
       whichresid, &            ! option for method of calculating residual
       whichsparse, &           ! option for method of doing elliptic solve
                                ! (BiCG, GMRES, standalone Trilinos, etc.)
       whichbabc,  &            ! option for basal boundary condition
       whichapprox, &           ! option for which Stokes approximation to use
                                ! 0 = SIA, 1 = SSA, 2 = Blatter-Pattyn HO
                                ! default = 2
       whichprecond             ! option for which preconditioner to use with 
                                !  structured PCG solver
                                ! 0 = none, 1 = diag, 2 = SIA

    !--------------------------------------------------------
    ! Local parameters
    !--------------------------------------------------------

    integer, parameter :: cmax = 300          ! max number of outer iterations

    !WHL - What I call resid_target is called minresid in glam_strs2

    real(dp), parameter :: resid_target = 1.0d-04   ! assume velocity fields have converged below this resid 
    real(dp), parameter :: NL_tol   = 1.0d-06   ! to have same criterion as JFNK

!WHL - UMC
    logical, parameter ::  &                 !TODO - Make this an input argument?
       unstable_manifold = .false. ! if true, do an unstable manifold correction
                                   ! (based on Hindmarsh & Payne, Ann Glac, 1996)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    real(dp), dimension(nx-1,ny-1) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex

    real(dp), dimension(nx-1,ny-1) ::   &
       stagusrf,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    real(dp), dimension(nx,ny) ::    &
       rmask                  ! = 1. where ice is present, else = 0.

    logical, dimension(nx,ny) ::     &
       active_cell,          &! true for active cells (thck > thklim and border locally owned vertices)
       floating_cell,        &! true for cells where ice is present and is floating
       ocean_cell             ! true for cells where topography is below sea level and ice is absent

    logical, dimension(nx-1,ny-1) :: &
       active_vertex          ! true for vertices of active cells

    real(dp), dimension(nz-1,nx,ny) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute effective viscosity
                          ! units: Pa yr^(1/n)

    real(dp), dimension(nz,nx-1,ny-1) ::   &
       usav, vsav,                 &! previous guess for velocity solution
       bu, bv,                     &! assembled load vector, divided into 2 parts
       resid_vec_u, resid_vec_v     ! residual vector, divided into 2 parts

    logical, dimension(nz,nx-1,ny-1) ::    &
       umask_dirichlet    ! Dirichlet mask for velocity (if true, u = v = 0)

    real(dp), dimension(27,nz,nx-1,ny-1) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 27 (node and its nearest neighbors in x, y and z direction) 
                          ! other dimensions = (k,i,j)

    real(dp) :: &
       resid_velo,          & ! quantity related to velocity convergence
       L2_norm,             & ! L2 norm of residual
       L2_target,           & ! nonlinear convergence target   !TODO - Make sure this is set correctly
       err,                 & ! solution error from sparse_easy_solve
       outer_it_criterion,  & ! current value of outer (nonlinear) loop converence criterion
       outer_it_target        ! target value for outer-loop convergence

    logical, save ::    &
       converged_soln = .false.    ! true if we get a converged solution for velocity

    integer ::  & 
       counter,         & ! outer (nonlinear) iteration counter
       niters             ! linear iteration count from sparse_easy_solve

    real(dp) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

!!    !Passed in as an argument instead
!!    real(dp), dimension(nx-1,ny-1) ::  &
!!       beta                   ! basal traction parameter

    ! The following are used only for the single-processor SLAP solver

    integer ::            &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1) ::  &
       NodeID                 ! ID for each node where we solve for velocity
                              ! For periodic BCs, halo node IDs will be copied
                              !  from the other side of the grid

    integer, dimension((nx-1)*(ny-1)*nz) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    type(sparse_matrix_type) ::  &
       matrix             ! sparse matrix for SLAP solver, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), allocatable ::   &   ! for SLAP solver
       rhs,             & ! right-hand-side (b) in Ax = b
       answer,          & ! answer (x) in Ax = b
       resid_vec          ! residual vector Ax - b

    integer ::    &
       matrix_order,    & ! order of matrix = number of rows
       nNonzero           ! upper bound for number of nonzero entries in sparse matrix

!WHL - debug
    integer :: i, j, k, m, n, r
    integer :: iA, jA, kA, colA
    real(dp) :: ds_dx, ds_dy
    real(dp) :: maxthck, maxusrf
    real(dp) :: sumuvel, sumvvel
    logical, parameter :: test_matrix = .false.
    integer, parameter :: test_order = 20
    integer :: rowi
    logical, parameter :: sia_test = .false.
!    logical, parameter :: sia_test = .true.

    integer :: nNonzeros    ! number of nonzero entries in structured matrices

!WHL - UMC
    real(dp), dimension(nz,nx-1,ny-1) ::   &
       uvel_old, vvel_old,         &! uvel and vvel from previous iteration
       ucorr_old, vcorr_old         ! correction vectors from previous iteration

    rtest = -999
    itest = 1
    jtest = 1
    if (this_rank == model%numerics%rdiag_local) then
       rtest = model%numerics%rdiag_local
       itest = model%numerics%idiag_local
       jtest = model%numerics%jdiag_local
    endif

    if (verbose .and. main_task) then
       print*, 'In glissade_velo_higher_solve'
       print*, 'itest, jtest, ktest =', itest, jtest, ktest
    endif

    !--------------------------------------------------------
    ! Assign local pointers to derived type components
    !--------------------------------------------------------

!    nx = model%general%ewn
!    ny = model%general%nsn
!    nz = model%general%upn

     dx = model%numerics%dew
     dy = model%numerics%dns

     thklim = model%numerics%thklim
     eus = model%climate%eus
 
     sigma => model%numerics%sigma(:)     
     thck => model%geometry%thck(:,:)
     usrf => model%geometry%usrf(:,:)
     topg => model%geometry%topg(:,:)

     flwa => model%temper%flwa(:,:,:)
     efvs => model%stress%efvs(:,:,:)
     beta => model%velocity%beta(:,:)

     uvel => model%velocity%uvel(:,:,:)
     vvel => model%velocity%vvel(:,:,:)

     kinbcmask => model%velocity%kinbcmask(:,:)

     whichefvs    = model%options%which_ho_efvs
     whichresid   = model%options%which_ho_resid
     whichsparse  = model%options%which_ho_sparse
     whichbabc    = model%options%which_ho_babc
     whichapprox  = model%options%which_ho_approx
     whichprecond = model%options%which_ho_precond
     
    !--------------------------------------------------------
    ! Convert input variables to appropriate units for this solver
    !TODO - Add some comments about units for this module
    !--------------------------------------------------------

    call glissade_velo_higher_scale_input(dx,     dy,      &
                                          thck,   usrf,    &
                                          topg,            &
                                          eus,    thklim,  &
                                          flwa,   efvs,    &
                                          beta,            &
                                          uvel,   vvel)

!WHL - hack for ishomC test
    if (trim(matrix_label) == 'ishomC_block') then
       print*, ' '
       print*, 'Re-initializing for ishomC block'
       call ishomC_block_init(model, dx, dy, thck, usrf, topg, beta)
    endif

    ! Set volume scale
    ! This is not strictly necessary, but dividing by this scale gives matrix coefficients 
    !  that are not quite so large.

    vol0 = 1.0d9    ! volume scale (m^3)

    if (whichapprox == HO_APPROX_SIA) then   ! SIA
!!          if (verbose .and. main_task) print*, 'Solving shallow-ice approximation'
       if (main_task) print*, 'Solving shallow-ice approximation'
       sia_factor = 1.d0
       ssa_factor = 0.d0
    elseif (whichapprox == HO_APPROX_SSA) then  ! SSA
!!          if (verbose .and. main_task) print*, 'Solving shallow-shelf approximation'
       if (main_task) print*, 'Solving shallow-shelf approximation'
       sia_factor = 0.d0
       ssa_factor = 1.d0
    else   ! Blatter-Pattyn higher-order 
!!          if (verbose .and. main_task) print*, 'Solving Blatter-Pattyn higher-order approximation'
       if (main_task) print*, 'Solving Blatter-Pattyn higher-order approximation'
       sia_factor = 1.d0
       ssa_factor = 1.d0
    endif


    if (test_matrix) then
       call solve_test_matrix(test_order, whichsparse)
    endif

    ! Make sure that the geometry and flow factor are correct in halo cells.
    ! These calls are commented out, since the halo updates are done in 
    !  module glissade.F90, before calling glissade_velo_higher_solve.

!    call parallel_halo(thck)
!    call parallel_halo(topg)
!    call parallel_halo(usrf)
!    call parallel_halo(flwa)

    !------------------------------------------------------------------------------
    ! Setup for higher-order solver: Compute nodal geometry, allocate storage, etc.
    ! These are quantities that do not change during the outer nonlinear loop. 
    !------------------------------------------------------------------------------

    maxthck = maxval(thck(:,:))
    maxthck = parallel_reduce_max(maxthck)
    maxusrf = maxval(usrf(:,:))
    maxusrf = parallel_reduce_max(maxusrf)

    if (verbose .and. main_task) then
       print*, ' '
       print*, 'nx, ny, nz:', nx, ny, nz
       print*, 'vol0:', vol0
       print*, ' '
       print*, 'thklim:', thklim
       print*, 'max thck:', maxthck
       print*, 'max usrf:', maxusrf
    endif

!WHL - debug
    if (verbose_state .and. this_rank==rtest) then

       print*, 'sigma coordinate:'
       do k = 1, nz
          print*, k, sigma(k)
       enddo

       print*, ' '
       print*, 'Thickness field, rank =', rtest
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(f6.0)',advance='no') thck(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'Topography field, rank =', rtest
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(f6.0)',advance='no') topg(i,j)
          enddo
          write(6,*) ' '
       enddo
       print*, ' '

       print*, 'Upper surface field, rank =', rtest
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(f6.0)',advance='no') usrf(i,j)
          enddo
          write(6,*) ' '
       enddo

       print*, ' '
       print*, 'flwa (Pa-3 yr-1), k = 1, rank =', rtest
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(e12.5)',advance='no') flwa(1,i,j)
          enddo
          write(6,*) ' '
       enddo

    endif
 
    !------------------------------------------------------------------------------
    ! Specify Dirichlet (u = 0) boundary conditions
    !------------------------------------------------------------------------------

    ! initialize
    umask_dirichlet(:,:,:) = .false.   

    if (whichbabc == HO_BABC_NO_SLIP) then
       ! Impose zero sliding everywhere at the bed.
       umask_dirichlet(nz,:,:) = .true.    ! u = v = 0 at bed
    endif

    ! use kinbcmask, typically read from file at initialization
    ! zero out entire column where kinbcmask = 1
    do j = 1,ny-1
       do i = 1, nx-1
          if (kinbcmask(i,j) == 1) then
             umask_dirichlet(:,i,j) = .true.
          endif
       enddo
    enddo

!WHL - debug
    if (verbose_state .and. this_rank==rtest) then

       print*, ' '
       print*, 'beta:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(e12.5)',advance='no') beta(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'umask_dirichlet, k = 1:'
       k = 1
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(L3)',advance='no') umask_dirichlet(k,i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'umask_dirichlet, k =', nz
       k = nz
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(L3)',advance='no') umask_dirichlet(k,i,j)
          enddo
          print*, ' '
       enddo

    endif   ! verbose_state

!WHLt2 - If iterating on the velocity during the remapping transport calculation,
!        we do not want these masks to change.
!TODO - Pull these out into a separate subroutine?

    !------------------------------------------------------------------------------
    ! Compute masks: 
    ! (1) real mask = 1.0 where ice is present, 0.0 elsewhere
    ! (2) floating mask = .true. where ice is present and is floating
    ! (3) ocean mask = .true. where topography is below sea level and ice is absent
    !------------------------------------------------------------------------------

    do j = 1, ny
       do i = 1, nx

          if (thck(i,j) > thklim) then
	     rmask(i,j) = 1.d0
          else
             rmask(i,j) = 0.d0
          endif

          if (topg(i,j) - eus < (-rhoi/rhoo)*thck(i,j)) then
             if (thck(i,j) > thklim) then
                floating_cell(i,j) = .true.
                ocean_cell(i,j) = .false.
             else
                floating_cell(i,j) = .false.
                ocean_cell(i,j) = .true.
             endif
          else
             floating_cell(i,j) = .false.
             ocean_cell(i,j) = .false.
          endif

       enddo
    enddo

!WHL - debug
    if (verbose_state .and. this_rank==rtest) then
       print*, ' '
       print*, 'ocean_cell, rank =', rtest
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(L3)',advance='no') ocean_cell(i,j)
          enddo
          print*, ' '
       enddo
    endif

!WHL - debug
    if (verbose_state .and. this_rank==rtest) then
       print*, ' '
       print*, 'floating_cell, rank =', rtest
       do j = ny, 1, -1
          do i = 1, nx
             write(6,'(L3)',advance='no') floating_cell(i,j)
          enddo
          print*, ' '
       enddo
    endif

    !------------------------------------------------------------------------------
    ! Compute ice thickness and upper surface on staggered grid
    ! (requires that thck and usrf are up to date in halo cells)
    !------------------------------------------------------------------------------

    call staggered_scalar(nx,           ny,         &
                          nhalo,        rmask,      &
                          thck,         stagthck)

    call staggered_scalar(nx,           ny,         &
                          nhalo,        rmask,      &
                          usrf,         stagusrf)

!WHL - debug
    if (verbose_state .and. this_rank==rtest) then
       print*, ' '
       print*, 'stagthck, rank =', rtest
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f6.0)',advance='no') stagthck(i,j)
          enddo
          print*, ' '
       enddo
    endif

    !------------------------------------------------------------------------------
    ! Compute vertices of each element.
    ! Identify the active cells (i.e., cells with thck > thklim,
    !  bordering a locally owned vertex) and active vertices (all vertices
    !  of active cells).
    !
    ! For the SLAP solver, count and assign a unique ID to each active node.
    ! TODO - Move SLAP-only calculation to another subroutine?
    !------------------------------------------------------------------------------

    call get_vertex_geometry(nx,          ny,         nz,  &   
                             nhalo,       sigma,           &
                             dx,          dy,              &
                             thck,        thklim,          &
                             stagusrf,    stagthck,        &
                             xVertex,     yVertex,         &
                             active_cell, active_vertex,   &
                             nNodesSolve, NodeID,          &
                             iNodeIndex,  jNodeIndex,  kNodeIndex)

!WHL - debug
    !TODO - Change to global node ID
    if (verbose_state .and. this_rank==rtest) then
       print*, ' '
       print*, 'NodeID before halo update, k = 1:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(i5)',advance='no') NodeID(1,i,j)
          enddo
          print*, ' '
       enddo
    endif

    ! Assign the appropriate ID to nodes in the halo.
    ! NOTE: This will work for single-processor runs with periodic BCs
    !       (e.g., ISMIP-HOM), but not for multiple processors.

    call staggered_parallel_halo(NodeID)

!WHL - debug
    !TODO - Change to global node ID
    if (verbose_state .and. this_rank==rtest) then
       print*, ' '
       print*, 'NodeID after halo update, k = 1:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(i5)',advance='no') NodeID(1,i,j)
          enddo
          print*, ' '
       enddo
    endif

    !------------------------------------------------------------------------------
    ! Compute the factor A^(-1/n) appearing in the expression for effective viscosity.
    ! This factor is often denoted as B in the literature.
    ! Note: The rate factor (flwa = A) is assumed to have units of Pa^(-n) yr^(-1).
    !       Thus flwafact = 0.5 * A^(-1/n) has units Pa yr^(1/n).
    !------------------------------------------------------------------------------

    flwafact(:,:,:) = 0.d0

    ! Loop over all cells that border locally owned vertices
    !TODO - Simply compute for all cells?  We should have flwa for all cells.

    do j = 1+nhalo, ny-nhalo+1
       do i = 1+nhalo, nx-nhalo+1
          if (active_cell(i,j)) then
             ! gn = exponent in Glen's flow law (= 3 by default)
             flwafact(:,i,j) = 0.5d0 * flwa(:,i,j)**(-1.d0/real(gn,dp))  
          endif
       enddo
    enddo

    if (verbose .and. this_rank==rtest) then
       n = ntest
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)
       print*, ' '
       print*, 'ntest, i, j, k:', n, i, j, k
    endif
 
    if (verbose .and. this_rank==rtest) then
       print*, ' '
       i = itest; j = jtest; k = ktest
       print*, 'itest, jtest, ktest:', i, j, k
       print*, 'NodeID =', NodeID(k,i,j)
!       print*, ' '
!       print*, 'thck(j+1):', thck(i:i+1,j+1)
!       print*, 'thck(j)  :', thck(i:i+1,j)
!       print*, 'stagthck :', stagthck(i,j)
       print*, ' '
       print*, 'usrf(j+1):', usrf(i:i+1,j+1)
       print*, 'usrf(j)  :', usrf(i:i+1,j)
       print*, 'stagusrf :', stagusrf(i,j)
       print*, ' '
       ds_dx = 0.5d0 * (stagusrf(i,j) + stagusrf(i,j-1) - stagusrf(i-1,j) - stagusrf(i-1,j-1)) / &
               (xVertex(i+1,j) - xVertex(i,j))
       ds_dy = 0.5d0 * (stagusrf(i,j) - stagusrf(i,j-1) + stagusrf(i-1,j) - stagusrf(i-1,j-1)) / &
               (yVertex(i,j+1) - yVertex(i,j))
       print*, 'Finite-difference ds/dx, ds_dy:', ds_dx, ds_dy
       print*, ' '                
    endif

    !------------------------------------------------------------------------------
    ! If using SLAP solver, then allocate space for the sparse matrix (A), rhs (b), 
    !  answer (x), and residual vector (Ax-b).
    !------------------------------------------------------------------------------

    if (whichsparse /= STANDALONE_PCG_STRUC .and.    &
        whichsparse /= STANDALONE_TRILINOS_SOLVER) then   

       matrix_order = 2*nNodesSolve    ! Is this exactly enough?
       nNonzero = matrix_order*54      ! 27 = node plus 26 nearest neighbors in hexahedral lattice   
                                       ! 54 = 2 * 27 (since solving for both u and v)

       allocate(matrix%row(nNonzero), matrix%col(nNonzero), matrix%val(nNonzero))
       allocate(rhs(matrix_order), answer(matrix_order), resid_vec(matrix_order))

       if (verbose_matrix) then
          print*, 'matrix_order =', matrix_order
          print*, 'nNonzero = ', nNonzero
       endif

    endif   ! SLAP solver
 
    if (verbose .and. this_rank==rtest) then
       print*, 'size(uvel, vvel) =', size(uvel), size(vvel)
       print*, 'sum (uvel, vvel) =', sum(uvel), sum(vvel)
       print*, 'size(usav, vsav) =', size(usav), size(vsav)
    endif

    !---------------------------------------------------------------
    ! Form and solve Ax = b
    !---------------------------------------------------------------

    if (main_task) then
       ! print some info to the screen to update on iteration progress
       print *, ' '
       print *, 'Running Glissade higher-order dynamics solver'
       print *, ' '
       if (whichresid == HO_RESID_L2NORM) then  ! use L2 norm of residual
          print *, 'iter #     resid (L2 norm)       target resid'
       else                                     ! residual based on velocity
          print *, 'iter #     velo resid            target resid'
       end if
       print *, ' '
    endif

    !TODO - Any ghost preprocessing needed, as in glam_strs2?

    !------------------------------------------------------------------------------
    ! set initial values 
    !------------------------------------------------------------------------------

!WHL - UMC
!TODO - Replace usav, vsav with uvel_old, vvel_old?

    counter = 0
    resid_velo = 1.d0
    L2_norm   = 1.0d20
    L2_target = 1.0d-4      !TODO - Is this the best value?  
                            ! Should it stay fixed during the outer loop, or should it evolve?

    outer_it_criterion = 1.0d10   ! guarantees at least one loop
    outer_it_target    = 1.0d-12 

    !------------------------------------------------------------------------------
    ! Assemble the load vector b
    ! This goes before the outer loop because the load vector
    !  does not change from one nonlinear iteration to the next.
    !------------------------------------------------------------------------------

    bu(:,:,:) = 0.d0
    bv(:,:,:) = 0.d0

    !------------------------------------------------------------------------------
    ! gravitational forcing
    !------------------------------------------------------------------------------

    call load_vector_gravity(nx,               ny,              &
                             nz,               nhalo,           &
                             nNodesPerElement,                  &
                             sigma,            active_cell,     &
                             xVertex,          yVertex,         &
                             stagusrf,         stagthck,        &
                             bu,               bv)

    !WHL - debug
!    print*, ' '
!    print*, 'After gravity, bv:'
!    do n = 1, nNodesSolve
!       i = iNodeIndex(n)
!       j = jNodeIndex(n)
!       k = kNodeIndex(n) 
!       print*, n, bv(k,i,j)
!    enddo

    !------------------------------------------------------------------------------
    ! lateral pressure at vertical ice edge
    !------------------------------------------------------------------------------

    call load_vector_lateral_bc(nx,               ny,              &
                                nz,               nhalo,           &
                                nNodesPerElement_2d,               &
                                sigma,                             &
                                floating_cell,    ocean_cell,      &
                                active_cell,                       &
                                xVertex,          yVertex,         &
                                stagusrf,         stagthck,        &
                                bu,               bv)

    !WHL - debug
!    print*, ' '
!    print*, 'After lateral BC, bv:'
!    do n = 1, nNodesSolve
!       i = iNodeIndex(n)
!       j = jNodeIndex(n)
!       k = kNodeIndex(n) 
!       print*, n, bv(k,i,j)
!    enddo

!WHL - debug - adjust bu and bv.
!!    bu(:,:,:) = bu(:,:,:) / 1.09d0
!!    bv(:,:,:) = bv(:,:,:) / 1.09d0

    !------------------------------------------------------------------------------
    ! main outer loop: iteration to solve the nonlinear problem
    !------------------------------------------------------------------------------

    do while (outer_it_criterion >= outer_it_target .and. counter < cmax)

       ! Advance the iteration counter

       counter = counter + 1

       if (verbose .and. main_task) then
          print*, ' '
          print*, 'Outer counter =', counter
          print*, 'whichresid =', whichresid
          print*, 'L2_norm, L2_target =', L2_norm, L2_target
          print*, 'resid_velo, resid_target =', resid_velo, resid_target
       endif

       ! save current velocity
       usav(:,:,:) = uvel(:,:,:)
       vsav(:,:,:) = vvel(:,:,:)

       !---------------------------------------------------------------------------
       ! Assemble the stiffness matrix A
       !---------------------------------------------------------------------------

       if (verbose_matrix .and. this_rank==rtest) print*, 'call assemble_stiffness_matrix'

       call assemble_stiffness_matrix(nx,               ny,              &
                                      nz,               nhalo,           &
                                      sigma,            active_cell,     &
                                      xVertex,          yVertex,         &
                                      uvel,             vvel,            &
                                      stagusrf,         stagthck,        &
                                      flwafact,                          &
                                      efvs,             whichefvs,       &
                                      Auu,              Auv,             &
                                      Avu,              Avv,             &
                                      sia_factor,       ssa_factor)
       
       if (verbose_matrix .and. this_rank==rtest) print*, 'Assembled the stiffness matrix'

!WHL - debug
       if (verbose_matrix .and. this_rank==rtest) then

          i = itest
          j = jtest
          k = ktest
          n = ntest
          r = rtest
          print*, ' '
          print*, 'Auu after assembly: rank, i, j, k, n =', r, i, j, k, n
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             m = indxA(iA,jA,kA)
             print*, 'iA, jA, kA, Auu:', iA, jA, kA, Auu(m, k, i, j)
          enddo
          enddo
          enddo

       endif

       !---------------------------------------------------------------------------
       ! Incorporate basal sliding boundary conditions
       ! Assume a linear sliding law for now

       ! Note: We could call this subroutine before the main outer loop if beta
       !       is assumed to be independent of velocity.  Putting the call here,
       !       however, allows for more general sliding laws.
       !---------------------------------------------------------------------------

!TODO - Test this subroutine

       if (whichbabc /= HO_BABC_NO_SLIP) then

       !WHL - debug
          if (verbose .and. main_task) then
             print*, ' '
             print*, 'max, min beta (Pa/(m/yr)) =', maxval(beta), minval(beta)
             print*, 'Call basal_sliding_bc'
          endif

          call basal_sliding_bc(nx,               ny,              &
                                nz,               nhalo,           &
                                nNodesPerElement_2d,               &
                                active_cell,      beta,            &
                                xVertex,          yVertex,         &
                                Auu,              Avv)

       endif

       if (verbose_matrix .and. this_rank==rtest) then
          i = itest
          j = jtest
          k = ktest
          n = ntest
!          print*, ' '
!          print*, 'Auu after basal sliding BC: rank, i, j, k, n =', r, i, j, k, n
          do kA = -1, 1
             do jA = -1, 1
                do iA = -1, 1
                   m = indxA(iA,jA,kA)
                   !             print*, 'iA, jA, kA, Auu:', iA, jA, kA, Auu(m, k, i, j)
                enddo
             enddo
          enddo
       endif

       !---------------------------------------------------------------------------
       ! Incorporate Dirichlet (u = 0) boundary conditions
       !---------------------------------------------------------------------------

       if (verbose .and. main_task) then
          print*, ' '
          print*, 'Call dirichlet_bc'
       endif

       call dirichlet_boundary_conditions(nx,              ny,                &
                                          nz,              nhalo,             &
                                          active_vertex,   umask_dirichlet,   &
                                          Auu,             Auv,               &
                                          Avu,             Avv,               &
                                          bu,              bv)

       !---------------------------------------------------------------------------
       ! Halo updates for matrices
       !
       ! These updates are not strictly necessary unless we're concerned about
       !  roundoff errors.
       ! But suppose we are comparing two entries that are supposed to be equal
       !  (e.g., to preserve symmetry), where entry 1 is owned by processor A and 
       !  entry 2 is owned by processor B.  
       ! Processor A might compute a local version of entry 2 in its halo, with 
       !  entry 2 = entry 1 locally.  But processor B's entry 2 might be different
       !  because of roundoff.  We need to make sure that processor B's value 
       !  is communicated to processor A.  If these values are slightly different, 
       !  they will be reconciled by the subroutine check_symmetry_assembled_matrix.
       !---------------------------------------------------------------------------
     
        call staggered_parallel_halo(Auu(:,:,:,:))
        call staggered_parallel_halo(Auv(:,:,:,:))
        call staggered_parallel_halo(Avu(:,:,:,:))
        call staggered_parallel_halo(Avv(:,:,:,:))

       !---------------------------------------------------------------------------
       ! Halo updates for load vectors
       !WHL - Not sure if these are necessary
       !---------------------------------------------------------------------------

    !WHL - debug
!    print*, ' '
!    print*, 'After Dirichlet conditions, bv:'
!    do n = 1, nNodesSolve
!       i = iNodeIndex(n)
!       j = jNodeIndex(n)
!       k = kNodeIndex(n) 
!       print*, n, bv(k,i,j)
!    enddo

!WHL - debug
    if (verbose_matrix .and. this_rank==rtest) then
!       print*, ' '
!       print*, 'Before halo update, bu(1,:,:), rank =', rtest
!       do j = ny-1, 1, -1
!          do i = 1, nx-1
!             write(6,'(e10.3)',advance='no') bu(1,:,j)
!          enddo
!          print*, ' '
!       enddo
    endif

        call staggered_parallel_halo(bu(:,:,:))
        call staggered_parallel_halo(bv(:,:,:))

!WHL - debug
    if (verbose_matrix .and. this_rank==rtest) then
!       print*, ' '
!       print*, 'After halo update, bu(1,:,:), rank =', rtest
!       do j = ny-1, 1, -1
!          do i = 1, nx-1
!             write(6,'(e10.3)',advance='no') bu(1,:,j)
!          enddo
!          print*, ' '
!       enddo
    endif

       !---------------------------------------------------------------------------
       ! Check symmetry of assembled matrix
       ! 
       ! There may be small differences from perfect symmetry due to roundoff errors.  
       ! If sufficiently small, these differences are fixed by averaging the two values 
       !  that should be symmetric.  Otherwise the code aborts.
       !
       ! Note: It may be OK to skip this check for production code.  However,
       !       small violations of symmetry are not tolerated well by some solvers.
       !       For example, the SLAP PCG solver with incomplete Cholesky preconditioning
       !       can crash if symmetry is not perfect. 
       !---------------------------------------------------------------------------

       if (verbose .and. main_task) print*, 'Check matrix symmetry'

       if (check_symmetry) then

          call check_symmetry_assembled_matrix(nx,          ny,      &
                                               nz,          nhalo,   &
                                               active_vertex,        &
                                               Auu,         Auv,     &
                                               Avu,         Avv)

       endif

       !WHL - debug
       !TODO - Determine how these very small values get in the matrix.

       !---------------------------------------------------------------------------
       ! 

       !---------------------------------------------------------------------------

       if (remove_small_values) then

          call remove_small_values_assembled_matrix(nx,          ny,      &
                                                    nz,          nhalo,   &
                                                    active_vertex,        &
                                                    Auu,         Auv,     &
                                                    Avu,         Avv)

       endif


       if (zero_uvel) then  ! enforce uvel = 0 everywhere

          print*, 'Zeroing out u components of velocity'

          ! zero out all matrices except Avv
          Auu(:,:,:,:) = 0.d0
          Auv(:,:,:,:) = 0.d0
          Avu(:,:,:,:) = 0.d0

          ! Put 1's on main diagnonal of Auu
          m = indxA(0,0,0)
          Auu(m,:,:,:) = 1.d0

          ! zero out u terms on RHS
          bu(:,:,:) = 0.d0

       elseif (zero_vvel) then   ! enforce vvel = 0 everywhere

          print*, 'Zeroing out v components of velocity'

          ! zero out all matrices except Auu
          Avv(:,:,:,:) = 0.d0
          Auv(:,:,:,:) = 0.d0
          Avu(:,:,:,:) = 0.d0

          ! Put 1's on main diagnonal of Avv
          m = indxA(0,0,0)
          Avv(m,:,:,:) = 1.d0

          ! zero out v terms on RHS
          bv(:,:,:) = 0.d0

       endif
           
!WHL - debug - Try stripping out columns from the matrix and see if it can still solve an SIA problem.

          if (sia_test) then

             if (verbose .and. main_task) then
                print*, ' '
                print*, 'SIA test: Removing Auv, Avu, and iA/jA columns'
             endif

             ! Set Auv = Avu = 0

             Auv(:,:,:,:) = 0.d0 
             Avu(:,:,:,:) = 0.d0 

             ! Remove terms from iA and jA columns

             do j = 1, ny-1
             do i = 1, nx-1
             do k = 1, nz
                do kA = -1,1
                do jA = -1,1
                do iA = -1,1
                   if (iA /= 0 .or. jA /= 0) then
                      m = indxA(iA,jA,kA)
                      Auu(m,k,i,j) = 0.d0
                      Avv(m,k,i,j) = 0.d0
                   endif
                enddo
                enddo
                enddo
             enddo
             enddo
             enddo

          endif   ! sia_test             

       !---------------------------------------------------------------------------
       ! Count nonzero elements in structured matrices (strictly diagnostic)
       !---------------------------------------------------------------------------

       nNonzeros = 0
       do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz
                do kA = -1, 1
                do jA = -1, 1
                do iA = -1, 1
                   m = indxA(iA,jA,kA)
                   if (Auu(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Auv(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Avu(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                   if (Avv(m,k,i,j) /= 0.d0) nNonzeros = nNonzeros + 1
                enddo 
                enddo
                enddo
             enddo  ! k
          endif     ! active_vertex
       enddo        ! i
       enddo        ! j

       if (verbose_matrix .and. this_rank==rtest) then
          print*, ' '
          print*, 'nNonzeros in structured matrices =', nNonzeros
       endif

!WHL - debug - print out some matrix values for test point

       if (verbose_matrix .and. this_rank==rtest) then
          i = itest
          j = jtest
          k = ktest
          n = ntest

          print*, ' '
      	  print*, 'i,j,k,n =', i, j, k, n
          print*, 'Auu sum =', sum(Auu(:,k,i,j))
          print*, 'Auv sum =', sum(Auv(:,k,i,j))
          print*, 'Avu sum =', sum(Avu(:,k,i,j))
          print*, 'Avv sum =', sum(Avv(:,k,i,j))

          print*, ' '
          print*, 'iA, jA, kA, Auu, Avv, Auv, Avu:'
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             m = indxA(iA,jA,kA)
             print*, iA, jA, kA, Auu(m,k,i,j), Avv(m,k,i,j), Auv(m,k,i,j), Avu(m,k,i,j) 
          enddo
          enddo
          enddo
     
          i = itest+1
          print*, ' '
      	  print*, 'i,j,k,n =', i, j, k, n
          print*, 'iA, jA, kA, Auu, Avv, Auv, Avu:'
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             m = indxA(iA,jA,kA)
             print*, iA, jA, kA, Auu(m,k,i,j), Avv(m,k,i,j), Auv(m,k,i,j), Avu(m,k,i,j) 
          enddo
          enddo
          enddo

          i = itest
          do kA = -1, 1
          do jA = -1, 1
          do iA = -1, 1
             if ( (k+kA >= 1 .and. k+kA <= nz)       &
                             .and.                   &
                  (i+iA >= 1 .and. i+iA <= nx-1)     &
                             .and.                   &
                  (j+jA >= 1 .and. j+jA <= ny-1) ) then
                colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
                if (colA > 0) then
                   m = indxA(iA,jA,kA)
!                   print*, ' '
!                   print*, 'iA, jA, kA, colA =', iA, jA, kA, colA
!                   print*, 'i, j, k for colA:', iNodeIndex(colA), jNodeIndex(colA), kNodeIndex(colA)
!                   print*, 'Auu =', Auu(m,k,i,j)
!                   print*, 'Auv =', Auv(m,k,i,j)
!                   print*, 'Avu =', Avu(m,k,i,j)
!                   print*, 'Avv =', Avv(m,k,i,j)
                endif
             endif
          enddo
          enddo
          enddo

          print*, 'i, j, k: ', i, j, k
          print*, 'bu =', bu(k,i,j)
          print*, 'bv =', bv(k,i,j)

       endif  ! verbose_matrix

       if (verbose_matrix .and. this_rank==rtest) then

          k = 1
          
          m = indxA(0,0,0)
          print*, ' '
          print*, 'Auu_diag, k =', k
          do j = ny-1, 1, -1
             do i = 1, nx-1
                write(6,'(e10.2)',advance='no'), Auu(m,k,i,j)
             enddo
             print*, ' '
          enddo

          print*, ' '
          print*, 'bu, k =', k
          do j = ny-1, 1, -1
             do i = 1, nx-1
                write(6,'(e10.2)',advance='no'), bu(k,i,j)
             enddo
             print*, ' '
          enddo

          print*, ' '
          print*, 'Avv_diag, k =', k
          do j = ny-1, 1, -1
             do i = 1, nx-1
                write(6,'(e10.2)',advance='no'), Avv(m,k,i,j)
             enddo
             print*, ' '
          enddo

          print*, ' '
          print*, 'bv, k =', k
          do j = ny-1, 1, -1
             do i = 1, nx-1
                write(6,'(e10.2)',advance='no'), bv(k,i,j)
             enddo
             print*, ' '
          enddo
          
       endif   ! verbose_matrix, this_rank==rtest

       if (write_matrix) then

          if (counter == 1) then    ! first outer iteration only
 
             call write_matrix_elements(nx,    ny,   nz,  &
                                        nNodesSolve, NodeID,       &
                                        iNodeIndex,  jNodeIndex,  &
                                        kNodeIndex,            &
                                        Auu,         Auv,   &
                                        Avu,         Avv,   &
                                        bu,          bv)

          endif

       endif


       if (whichsparse == STANDALONE_PCG_STRUC) then   ! standalone PCG for structured grid
                                                       ! works for both serial and parallel runs

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2 norm
          !------------------------------------------------------------------------

          if (verbose .and. this_rank==rtest) then
             print*, ' '
             print*, 'Compute residual vector'
          endif

          call compute_residual_vector(nx,          ny,            &
                                       nz,          nhalo,         &
                                       active_vertex,              &
                                       Auu,         Auv,           &
                                       Avu,         Avv,           &
                                       bu,          bv,            &
                                       uvel,        vvel,          &
                                       resid_vec_u, resid_vec_v,   &
                                       L2_norm)

          if (verbose .and. this_rank==rtest) then
             print*, ' '
             print*, 'counter, L2_norm =', counter, L2_norm
          endif

          ! preprocessing is not needed because the matrix, rhs, and solution already
          ! have the required data structure

!          if (verbose .and. this_rank==rtest) then
!             print*, ' '
!             print*, 'Calling structured PCG solver'
!          endif

          !------------------------------------------------------------------------
          ! Call linear PCG solver, compute uvel and vvel on local processor
          !------------------------------------------------------------------------

          call pcg_solver_structured(nx,           ny,            &
                                     nz,           nhalo,         &
                                     indxA,        active_vertex, &
                                     Auu,          Auv,           &
                                     Avu,          Avv,           &
                                     bu,           bv,            &
                                     uvel,         vvel,          &
                                     whichprecond, err,           &
                                     niters)

          if (verbose .and. this_rank==rtest) then
             print*, 'Solved the linear system, niters, err =', niters, err
          endif

          ! Halo updates for uvel and vvel
          !TODO - Put these after the 'endif'?

          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)

       elseif (whichsparse /= STANDALONE_TRILINOS_SOLVER) then   ! one-processor SLAP solve   
          
          !------------------------------------------------------------------------
          ! Given the stiffness matrices (Auu, etc.) and load vector (bu, bv) in
          !  structured format, form the global matrix and rhs in SLAP format.
          !------------------------------------------------------------------------

          if (verbose) print*, 'Form global matrix in sparse format'
 
          matrix%order = matrix_order
          matrix%nonzeros = nNonzero
          matrix%symmetric = .false.

          call slap_preprocess(nx,           ny,          &   
                               nz,           nNodesSolve, &
                               NodeID,                    &
                               iNodeIndex,   jNodeIndex,  &
                               kNodeIndex,                &
                               Auu,          Auv,         &
                               Avu,          Avv,         &
                               bu,           bv,          &
                               uvel,         vvel,        &
                               matrix_order, nNonzero,    &
                               matrix,       rhs,         &
                               answer)

          !------------------------------------------------------------------------
          ! Compute the residual vector and its L2_norm
          !------------------------------------------------------------------------

          call slap_compute_residual_vector(matrix,    answer,   rhs,  &
                                            resid_vec, L2_norm)

!WHL - bug check
          if (verbose) then
             print*, ' '
             print*, 'Before linear solve: n, row, col, val:'
             print*, ' '
             do n = 1, 10              
                print*, n, matrix%row(n), matrix%col(n), matrix%val(n)
             enddo

             print*, 'L2_norm of residual =', L2_norm
             print*, 'Call sparse_easy_solve, counter =', counter
          endif

          !------------------------------------------------------------------------
          ! Solve the linear matrix problem
          !------------------------------------------------------------------------

          call sparse_easy_solve(matrix, rhs,    answer,  &
                                 err,    niters, whichsparse)

          if (verbose .and. this_rank==rtest) then
             print*, 'Solved the linear system, niters, err =', niters, err
             print*, ' '
!!             print*, 'n, u, v (m/yr):', ntest, answer(2*ntest-1), answer(2*ntest)
!!           do n = 1, matrix%order 
!!              print*, n, answer(n)
!!           enddo
          endif

          !------------------------------------------------------------------------
          ! Put the velocity solution back into 3D arrays
          !------------------------------------------------------------------------

          call slap_postprocess(nNodesSolve,  answer,                   &
                                iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                uvel,         vvel)

          !------------------------------------------------------------------------
          ! Halo updates for uvel and vvel
          !------------------------------------------------------------------------

          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)

#ifdef TRILINOS
       else    ! solve with Trilinos

!WHL - Commented out for now until testing later with Trilinos
!!           call solvewithtrilinos(rhs, answer, linearSolveTime)
!!           totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
          !  write(*,*) 'Total linear solve time so far', totalLinearSolveTime                                           
#endif

       endif   ! whichsparse (STANDALONE_PCG_STRUC or not) 

       ! Any ghost postprocessing needed, as in glam_strs2?

!WHL - bug check
       if (verbose .and. this_rank==rtest) then

          i = itest
          j = jtest
          print*, ' '
          print*, 'After postprocess: rank, i, j =', this_rank, i, j
          print*, 'k, uvel, vvel (m/yr):'
          do k = 1, nz
             print*, k, uvel(k,i,j), vvel(k,i,j)               
          enddo

       endif

       !---------------------------------------------------------------------------
       ! Compute residual quantities based on the velocity solution
       !---------------------------------------------------------------------------

       call compute_residual_velocity(nhalo,  whichresid,   &
                                      uvel,   vvel,        &
                                      usav,   vsav,        &
                                      resid_velo)
  
       if (verbose_state .and. this_rank==rtest) then
!!       if verbose_residual .and. this_rank==rtest) then

          print*, ' '
          print*, 'whichresid, resid_velo =', whichresid, resid_velo

!          print*, ' '
!          print*, 'uvel - usav, k = 1:'
          do j = ny-1, 1, -1
             do i = 1, nx-1
!                write(6,'(f6.2)',advance='no') (uvel(1,i,j) - usav(1,i,j))
             enddo 
!         print*, ' '
          enddo

!          print*, ' '
!          print*, 'vvel - vsav, k = 1:'
          do j = ny-1, 1, -1
             do i = 1, nx-1
!                write(6,'(f6.2)',advance='no') (vvel(1,i,j) - vsav(1,i,j))
             enddo 
!             print*, ' '
          enddo

          print*, ' '
          print*, 'uvel (m/yr), k = 1:'
          do j = ny-1, 1, -1
             do i = 1, nx-1
                write(6,'(f12.0)',advance='no') uvel(1,i,j)
             enddo 
             print*, ' '
          enddo
          print*, ' '
          print*, 'vvel (m/yr), k = 1:'
          do j = ny-1, 1, -1
             do i = 1, nx-1
                write(6,'(f12.0)',advance='no') vvel(1,i,j)
             enddo 
             print*, ' '
          enddo

          print*, 'max(uvel, vvel) =', maxval(uvel), maxval(vvel)
       endif

       !WHL - UMC
       !---------------------------------------------------------------------------
       ! Optionally, do an unstable manifold correction to improve convergence
       ! of the Picard iteration.
       !---------------------------------------------------------------------------

       if (unstable_manifold) then

          !TODO - Replace usav, vsav with uvel_old, vvel_old?
          if (counter==1) then
             uvel_old(:,:,:) = uvel(:,:,:)
             vvel_old(:,:,:) = vvel(:,:,:)
             ucorr_old(:,:,:) = 0.d0
             vcorr_old(:,:,:) = 0.d0
          endif                      

          ! correct uvel and vvel as needed

          call unstable_manifold_correction(nx,        ny,        &
                                            nz,        nhalo,     &
                                            uvel,      vvel,      &
                                            uvel_old,  vvel_old,  &
                                            ucorr_old, vcorr_old)

          ! Halo updates for uvel and vvel
          ! TODO - Are these needed?

          call staggered_parallel_halo(uvel)
          call staggered_parallel_halo(vvel)

       endif

       !---------------------------------------------------------------------------
       ! Write diagnostics (iteration number, max residual, and location of max residual
       ! (send output to the screen or to the log file, per whichever line is commented out) 
       !---------------------------------------------------------------------------

       !TODO - Find out if this note from glam_strs2 still applies:
       ! "Can't use main_task flag because main_task is true for all processors in case of parallel_single"

       if (main_task) then
          if (whichresid == HO_RESID_L2NORM) then
             if (verbose) then
                print*, ' '
                print*, 'Using L2 norm convergence criterion'
                print*, 'iter#, L2 norm, L2 target:', counter, L2_norm, L2_target
             else   ! standard short message
                print '(i4,2g20.6)', counter, L2_norm, L2_target
             endif
             !write(message,'(i4,3g20.6)') counter, L2_norm, L2_target
             !call write_log (message)
          else
             if (verbose) then
                print*, ' '
                print*, 'Using velocity residual convergence criterion'
                print*, 'iter#, resid velo, resid target:', counter, resid_velo, resid_target
             else
                print '(i4,2g20.6)', counter, resid_velo, resid_target
             endif
             !write(message,'(i4,2g20.6)') counter, resid_velo, resid_target
             !call write_log (message)
          end if
       endif

       !---------------------------------------------------------------------------
       ! update the outer loop stopping criterion
       !---------------------------------------------------------------------------

       if (whichresid == HO_RESID_L2NORM) then
          outer_it_criterion = L2_norm
          outer_it_target = L2_target      ! L2_target is currently set to 1.d-4 and held constant
       else
          outer_it_criterion = resid_velo
          outer_it_target = resid_target   ! resid_target is currently a parameter = 1.d-4  
       end if

    enddo  ! while (outer_it_criterion >= outer_it_target .and. counter < cmax)

    if (counter < cmax) converged_soln = .true.

    if (verbose .and. converged_soln) then

       if (main_task) then
          print*, ' '
          print*, 'GLISSADE SOLUTION HAS CONVERGED, outer counter =', counter
       endif

    endif

    !------------------------------------------------------------------------------
    ! Clean up
    !------------------------------------------------------------------------------

    if (whichsparse /= STANDALONE_PCG_STRUC .and.   &
        whichsparse /= STANDALONE_TRILINOS_SOLVER) then
       deallocate(matrix%row, matrix%col, matrix%val)
       deallocate(rhs, answer, resid_vec)
    endif

    !------------------------------------------------------------------------------
    ! Convert output variables to appropriate units for Glimmer-CISM
    ! (generally dimensionless)
    !------------------------------------------------------------------------------

    call glissade_velo_higher_scale_output(thck,   usrf,   &
                                           topg,           &
                                           flwa,   efvs,   &
                                           beta,           &
                                           uvel,   vvel)

  end subroutine glissade_velo_higher_solve

!****************************************************************************

    subroutine glissade_velo_higher_scale_input(dx,     dy,     &
                                                thck,   usrf,   &
                                                topg,           &
                                                eus,    thklim, &
                                                flwa,   efvs,   &
                                                beta,           &
                                                uvel,   vvel)

    !--------------------------------------------------------
    ! Convert input variables (generally dimensionless)
    ! to appropriate units for the glissade_velo_higher solver.
    !--------------------------------------------------------

    real(dp), intent(inout) ::   &
       dx, dy                  ! grid cell length and width 

    real(dp), dimension(:,:), intent(inout) ::   &
       thck,                &  ! ice thickness
       usrf,                &  ! upper surface elevation
       topg                    ! elevation of topography

    real(dp), intent(inout) ::   &
       eus,  &                 ! eustatic sea level (= 0 by default)
       thklim                  ! minimum cell thickness for active cells

    real(dp), dimension(:,:,:), intent(inout) ::  &
       flwa,   &               ! flow factor in units of Pa^(-n) yr^(-1)
       efvs                    ! effective viscosity (Pa yr)

    real(dp), dimension(:,:), intent(inout)  ::  &
       beta                    ! basal traction parameter (Pa/(m/yr))

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel              ! velocity components (m/yr)

    ! grid cell dimensions: rescale from dimensionless to m
    dx = dx * len0
    dy = dy * len0

    ! ice geometry: rescale from dimensionless to m
    thck = thck * thk0
    usrf = usrf * thk0
    topg = topg * thk0
    eus  = eus  * thk0
    thklim = thklim * thk0

    ! rate factor: rescale from dimensionless to Pa^(-n) yr^(-1)
    flwa = flwa * (vis0*scyr)

    ! effective viscosity: rescale from dimensionless to Pa yr
    efvs = efvs * (evs0/scyr)

    ! beta: rescale from dimensionless to Pa/(m/yr)
    beta = beta * tau0/(vel0*scyr)

    ! ice velocity: rescale from dimensionless to m/yr
    uvel = uvel * (vel0*scyr)
    vvel = vvel * (vel0*scyr)

    end subroutine glissade_velo_higher_scale_input

!****************************************************************************

    subroutine glissade_velo_higher_scale_output(thck,   usrf,   &
                                                 topg,           &
                                                 flwa,   efvs,   &
                                                 beta,           &
                                                 uvel,   vvel)

    !--------------------------------------------------------
    ! Convert output variables to appropriate Glimmer-CISM units
    ! (generally dimensionless)
    !--------------------------------------------------------

    real(dp), dimension(:,:), intent(inout) ::  &
       thck,                 &  ! ice thickness
       usrf,                 &  ! upper surface elevation
       topg                     ! elevation of topography

    real(dp), dimension(:,:,:), intent(inout) ::  &
       flwa,   &                ! flow factor in units of Pa^(-n) yr^(-1)
       efvs                     ! effective viscosity (Pa yr)

    real(dp), dimension(:,:), intent(inout)  ::  &
       beta                     ! basal traction parameter (Pa/(m/yr))

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel               ! velocity components (m/yr)

    ! Convert geometry variables from m to dimensionless units
    thck = thck / thk0
    usrf = usrf / thk0
    topg = topg / thk0

    ! Convert flow factor from Pa^(-n) yr^(-1) to dimensionless units
    flwa = flwa / (vis0*scyr)

    ! Convert effective viscosity from Pa yr to dimensionless units
    efvs = efvs / (evs0/scyr)

    ! Convert beta from Pa/(m/yr) to dimensionless units
    beta = beta / (tau0/(vel0*scyr))

    ! Convert velocity from m/yr to dimensionless units
    uvel = uvel / (vel0*scyr)
    vvel = vvel / (vel0*scyr)

    end subroutine glissade_velo_higher_scale_output

!****************************************************************************

  subroutine staggered_scalar(nx,           ny,        &
                              nhalo,        rmask,     &
                              var,          stagvar)

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nhalo                    ! number of rows/columns of halo cells

    real(dp), dimension(nx,ny), intent(in) ::        &
       rmask                    ! = 1. where ice is present, else = 0.

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! scalar field defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field defined at cell corners

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j

    real(dp) :: sumvar, summask

    stagvar(:,:) = 0.d0

    ! Loop over all vertices of locally owned cells
    ! Average the input field over the neighboring cells with ice present (rmask = 1)

    do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = rmask(i,j+1)*var(i,j+1) + rmask(i+1,j+1)*var(i+1,j+1)  &
                 + rmask(i,j)  *var(i,j)   + rmask(i+1,j)  *var(i+1,j)	  
          summask = rmask(i,j+1) + rmask(i+1,j+1) + rmask(i,j) + rmask(i+1,j)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
    enddo  

  end subroutine staggered_scalar

!****************************************************************************

  subroutine get_vertex_geometry(nx,          ny,          nz,      &   
                                 nhalo,       sigma,                &
                                 dx,          dy,                   &
                                 thck,        thklim,               &
                                 stagusrf,    stagthck,             &
                                 xVertex,     yVertex,              &
                                 active_cell, active_vertex,        &
                                 nNodesSolve, NodeID,               & 
                                 iNodeIndex,  jNodeIndex,  kNodeIndex)
                            
    !----------------------------------------------------------------
    ! Compute coordinates for each vertex.
    ! Identify and count the active cells and vertices for the finite-element calculations.
    ! Active cells include all cells that contain ice (thck > thklin) and border locally owned vertices.
    ! Active vertices include all vertices of active cells.
    !
    ! Also compute some node indices needed for the SLAP single-processor solver.
    !TODO - Move SLAP part to different subroutine?
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &    ! number of grid cells in each direction
       nz,                   &    ! number of vertical levels where velocity is computed
       nhalo                      ! number of halo layers

    real(dp), dimension(nz), intent(in) :: &
       sigma                  ! sigma vertical coordinate

    real(dp), intent(in) ::  &
       dx,  dy                ! grid cell length and width (m)
                              ! assumed to have the same value for each grid cell

    real(dp), dimension(nx,ny), intent(in) ::  &
       thck                   ! ice thickness

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    real(dp), dimension(nx-1,ny-1), intent(out) :: &
       xVertex, yVertex       ! x and y coordinates of each vertex

    logical, dimension(nx,ny), intent(out) :: &
       active_cell            ! true for active cells 
                              ! (thck > thklim, bordering a locally owned vertex)

    logical, dimension(nx-1,ny-1), intent(out) :: &
       active_vertex          ! true for vertices of active cells

    ! The remaining input/output arguments are for the SLAP solver

    integer, intent(out) :: &
       nNodesSolve            ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(out) ::  &
       NodeID                 ! ID for each node where we solve for velocity

    integer, dimension((nx-1)*(ny-1)*nz), intent(out) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, k

    !----------------------------------------------------------------
    ! Compute the x and y coordinates of each vertex.
    ! By convention, vertex (i,j) lies at the NE corner of cell(i,j).
    !----------------------------------------------------------------

    xVertex(:,:) = 0.d0
    yVertex(:,:) = 0.d0
    do j = 1, ny-1
    do i = 1, nx-1
       xVertex(i,j) = dx * i
       yVertex(i,j) = dy * j
    enddo
    enddo

    if (verbose_init) then
       write(6,*) ' '
!       write(6,*) 'Vertex coordinates:'
       do j = 1, ny-1
!          write(6,*) ' '    
          do i = 1, nx-1
!             write(6,*) 'i, j, xV, yV:', i, j, xVertex(i,j), yVertex(i,j)
          enddo
       enddo
    endif

    ! Identify the active cells.
    ! Include all cells that border locally owned vertices and contain ice.

    active_cell(:,:) = .false.

    do j = 1+nhalo, ny-nhalo+1
    do i = 1+nhalo, nx-nhalo+1
       if (thck(i,j) > thklim) then
          active_cell(i,j) = .true.
       endif
    enddo
    enddo

    ! Identify the active vertices.
    ! Include all vertices of active cells.

    active_vertex(:,:) = .false.

    do j = 1+nhalo, ny-nhalo+1  ! include east and north layer of halo cells
    do i = 1+nhalo, nx-nhalo+1
       if (active_cell(i,j)) then
          active_vertex(i-1:i, j-1:j) = .true.  ! all vertices of this cell
       endif
    enddo
    enddo

    ! Identify and count the nodes where we must solve for the velocity.
    ! This indexing is used for pre- and post-processing of the assembled matrix
    !  when we call the SLAP solver (one processor only).
    ! It is not required by the structured solver.
    !TODO - Move to separate subroutine?

    nNodesSolve = 0
    NodeID(:,:,:) = 0
    iNodeIndex(:) = 0
    jNodeIndex(:) = 0
    kNodeIndex(:) = 0

    do j = nhalo+1, ny-nhalo    ! locally owned vertices only
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then   ! all nodes in column are active
          do k = 1, nz               
             nNodesSolve = nNodesSolve + 1   
             NodeID(k,i,j) = nNodesSolve   ! unique index for each node
             iNodeIndex(nNodesSolve) = i
             jNodeIndex(nNodesSolve) = j
             kNodeIndex(nNodesSolve) = k

             if (verbose .and. this_rank==rtest .and. nNodesSolve==ntest .and. k < nz) then
                print*, ' '
                print*, 'i, j, k, n:', i, j, k, nNodesSolve
                print*, 'sigma, stagusrf, stagthck:', sigma(k), stagusrf(i,j), stagthck(i,j)
                print*, 'dx, dy, dz:', xVertex(i,j) - xVertex(i-1,j), &
                                       yVertex(i,j) - yVertex(i,j-1), &
                                       (sigma(k+1) - sigma(k)) * stagthck(i,j)
             endif

           enddo   ! k
        endif      ! active vertex
    enddo          ! i
    enddo          ! j

    if (verbose .and. this_rank==rtest) then
       print*, ' '
       print*, 'nNodesSolve =', nNodesSolve
    endif

  end subroutine get_vertex_geometry

!****************************************************************************

  subroutine load_vector_gravity(nx,               ny,              &
                                 nz,               nhalo,           &
                                 nNodesPerElement,                  &
                                 sigma,            active_cell,     &
                                 xVertex,          yVertex,         &
                                 stagusrf,         stagthck,        &
                                 bu,               bv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    integer, intent(in) :: nNodesPerElement

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       s                  ! upper surface elevation at nodes

    real(dp)  ::   &
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    real(dp), dimension(nNodesPerElement) ::   &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis functions, evaluated at quad pts

    real(dp) ::    &
       ds_dx, ds_dy       ! horizontal gradient of upper surface elevation

    integer :: i, j, k, n, p

    integer :: iNode, jNode, kNode

    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       if (active_cell(i,j)) then

          do k = 1, nz-1    ! loop over elements in this column 
                            ! assume k increases from upper surface to bed

             if (verbose .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, 'i, j, k:', i, j, k
                print*, ' '
             endif

             ! compute spatial coordinates, velocity, and upper surface elevation for each node

             do n = 1, nNodesPerElement

                ! Determine (k,i,j) for this node
                ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
                ! Indices for other nodes are computed relative to this node.
                iNode = i + ishift(7,n)
                jNode = j + jshift(7,n)
                kNode = k + kshift(7,n)

                x(n) = xVertex(iNode,jNode)
                y(n) = yVertex(iNode,jNode)
                z(n) = stagusrf(iNode,jNode) - sigma(kNode)*stagthck(iNode,jNode)
                s(n) = stagusrf(iNode,jNode)
 
                if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!!                   print*, ' '
                   print*, 'i, j, k, n, x, y, z, s:', i, j, k, n, x(n), y(n), z(n), s(n)
!!                   print*, 'iNode, jNode, kNode, sigma:', iNode, jNode, kNode, sigma(kNode)
                endif

             enddo   ! nodes per element

             ! Loop over quadrature points for this element
   
             do p = 1, nQuadPoints

                ! Evaluate the derivatives of the element basis functions
                ! at this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

                call get_basis_function_derivatives(nNodesPerElement,                             &
                                                    x(:),          y(:),          z(:),           &
                                                    dphi_dxr(:,p), dphi_dyr(:,p), dphi_dzr(:,p),  &
                                                    dphi_dx(:),    dphi_dy(:),    dphi_dz(:),     &
                                                    detJ , i, j, k, p                      )

!WHL - debug
                if (detJ < 0.d0) then
                   print*, 'detJ < 0:, i, j, k, p, detJ =', i, j, k, p, detJ
                endif

                ! Evaluate ds_dx and ds_dy at this quadrature point
 
                ds_dx = 0.d0
                ds_dy = 0.d0

                do n = 1, nNodesPerElement
                   ds_dx = ds_dx + dphi_dx(n) * s(n)
                   ds_dy = ds_dy + dphi_dy(n) * s(n)
                enddo

                if (verbose_load .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                   print*, ' '
                   print*, 'Increment load vector, i, j, k, p =', i, j, k, p
                   print*, 'ds/dx, ds/dy =', ds_dx, ds_dy
                   print*, 'detJ/vol0 =', detJ/vol0
                   print*, 'detJ/vol0* (ds/dx, ds/dy) =', detJ/vol0*ds_dx, detJ/vol0*ds_dy
                endif

                ! Increment the load vector with the gravitational contribution from
                ! this quadrature point.

                do n = 1, nNodesPerElement

                   ! Determine (k,i,j) for this node
                   iNode = i + ishift(7,n)
                   jNode = j + jshift(7,n)
                   kNode = k + kshift(7,n)
         
                   bu(kNode,iNode,jNode) = bu(kNode,iNode,jNode) - rhoi*grav * wqp(p) * detJ/vol0 * ds_dx * phi(n,p)
                   bv(kNode,iNode,jNode) = bv(kNode,iNode,jNode) - rhoi*grav * wqp(p) * detJ/vol0 * ds_dy * phi(n,p)

                   if (verbose_load .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest .and. p==ptest) then
!                      print*, ' '
                      print*, 'n, phi(n), delta(bu), delta(bv):', n, phi(n,p), &
                               rhoi*grav*wqp(p)*detJ/vol0 * ds_dx * phi(n,p), &
                               rhoi*grav*wqp(p)*detJ/vol0 * ds_dy * phi(n,p)
                   endif

                enddo   ! nNodesPerElement

             enddo      ! nQuadPoints

          enddo         ! k

       endif            ! active cell

    enddo               ! i
    enddo               ! j

  end subroutine load_vector_gravity

!****************************************************************************

  subroutine load_vector_lateral_bc(nx,               ny,              &
                                    nz,               nhalo,           &
                                    nNodesPerElement_2d,               &
                                    sigma,                             &
                                    floating_cell,    ocean_cell,      &
                                    active_cell,                       &
                                    xVertex,          yVertex,         &
                                    stagusrf,         stagthck,        &
                                    bu,               bv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    integer, intent(in) :: nNodesPerElement_2d

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    !TODO - Use active_cell and not floating/ocean?
    logical, dimension(nx,ny), intent(in) ::  &
       active_cell,                 &! true if cell contains ice and borders a locally owned vertex
       floating_cell,               &! true if ice is present and is floating
       ocean_cell                    ! true if topography is below sea level and ice is absent

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    ! Sum over elements in active cells 
    ! Loop over cells that contain locally owned vertices
    !TODO - Make sure we don't step out of bounds
    !       With more care, we could skip some computations for vertices that are not locally owned.
    !TODO - Generalize for all grid cells at ice edge, including those not floating.
    !       if (active_cell(i,j)) then
    !          if (.not. active_cell(i-1,j)) call lateral_bc('west')
    !            etc.
!WHL - debug
!    print*, ' '
!    print*, 'In load_vector_lateral_bc'

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       !WHL - debug
       if (verbose_shelf .and. i==itest .and. j==jtest .and. this_rank==rtest) then
          print*, 'i, j =', i, j
          print*, 'active =', active_cell(i,j)
          print*, 'floating =', floating_cell(i,j)
          print*, 'ocean (i-1:i,j)  =', ocean_cell(i-1:i, j) 
          print*, 'ocean (i-1:i,j-1)=', ocean_cell(i-1:i, j-1) 
       endif

       if (active_cell(i,j)) then    ! ice is present
!!       if (floating_cell(i,j)) then   ! ice is present and is floating

!WHL - debug
!          print*, 'Floating:', i, j

          ! Note: The lateral shelf BC assumed below applies only to floating ice.
          !       It will not work for marine-based ice that borders the ocean but is not floating.
          !TODO - Make the BC more general.

!!          if (.not. active_cell(i-1,j)) then  ! compute lateral BC for west face
          if (ocean_cell(i-1,j)) then ! compute lateral BC for west face

!WHL - debug
!          print*, '   Ocean west:', i-1, j

             call lateral_shelf_bc(nx,                   ny,          &
                                   nz,                   sigma,       &
                                   nNodesPerElement_2d,  'west',     &
                                   i,                    j,           &
                                   stagusrf,             stagthck,    &
                                   xVertex,              yVertex,     &
                                   bu,                   bv)

          endif

!!          if (.not. active_cell(i+1,j)) then  ! compute lateral BC for east face
          if (ocean_cell(i+1,j)) then ! compute lateral BC for east face

!WHL - debug
!          print*, '   Ocean east:', i+1, j

             call lateral_shelf_bc(nx,                   ny,          &
                                   nz,                   sigma,       &
                                   nNodesPerElement_2d,  'east',     &
                                   i,                    j,           &
                                   stagusrf,             stagthck,    &
                                   xVertex,              yVertex,     &
                                   bu,                   bv)

          endif

          if (.not. active_cell(i,j-1)) then  ! compute lateral BC for south face
!!          if (ocean_cell(i,j-1)) then ! compute lateral BC for south face

!WHL - debug
!          print*, '   Ocean south:', i, j-1

             call lateral_shelf_bc(nx,                   ny,          &
                                   nz,                   sigma,       &
                                   nNodesPerElement_2d,  'south',     &
                                   i,                    j,           &
                                   stagusrf,             stagthck,    &
                                   xVertex,              yVertex,     &
                                   bu,                   bv)

          endif

          if (.not. active_cell(i,j+1)) then  ! compute lateral BC for north face
!!          if (ocean_cell(i,j+1)) then ! compute lateral BC for north face

!WHL - debug
!          print*, '   Ocean north:', i, j+1

             !WHL - debug
             if (verbose_shelf .and. i==itest .and. j==jtest .and. this_rank==rtest) then
                print*, 'Call lateral_shelf_bc, north face, i, j =', i, j
             endif

             call lateral_shelf_bc(nx,                   ny,          &
                                   nz,                   sigma,       &
                                   nNodesPerElement_2d,  'north',     &
                                   i,                    j,           &
                                   stagusrf,             stagthck,    &
                                   xVertex,              yVertex,     &
                                   bu,                   bv)

          endif

       endif      ! floating_cell

    enddo         ! i
    enddo         ! j

  end subroutine load_vector_lateral_bc

!****************************************************************************

  subroutine lateral_shelf_bc(nx,                  ny,           &
                              nz,                  sigma,        &
                              nNodesPerElement_2d, face,         &
                              iCell,               jCell,        &
                              stagusrf,            stagthck,     &
                              xVertex,             yVertex,      &
                              bu,                  bv)

!TODO -  Verify that the sum of contributions to bu or bv over an entire column
!        is equal to p_av times the surface area of the column.

    !----------------------------------------------------------------------------------
    ! Determine the contribution to the load vector from ice and water pressure at the
    !  vertical boundary between ice and ocean (or alternatively, from ice pressure alone
    !  at a vertical boundary between ice and air).
    !
    ! This subroutine computes the vertically averaged hydrostatic pressure at a vertical face
    !  associated with the grid cell column (iCell, jCell).
    !
    ! At a given point, this pressure is proportional to the difference between
    ! (1) the vertically averaged pressure exerted outward (toward the ocean) by the ice front
    ! (2) the vertically averaged pressure exerted by the ocean back toward the ice
    ! 
    ! (1) is given by p_out = 0.5*rhoi*grav*H
    ! (2) is given by p_in  = 0.5*rhoi*grav*H*(rhoi/rhoo) for a floating shelf
    !                       = 0.5*rhoo*grav*H*(1 - s/H)^2 for s <= H but ice not necessarily afloat
    !
    ! The second term goes to zero for a land-terminating cliff. 
    ! The two pressure terms are opposite in sign, so the net vertically averaged pressure,
    !  directed toward the ocean (or air), is given by
    ! 
    !                    p_av = 0.5*rhoi*grav*H 
    !                         - 0.5*rhoo*grav*H * (1 - min((s/H),1)^2
    ! 
    ! Here we sum over quadrature points for each ocean-bordering face of each element.
    ! The contribution from each quadrature point to node N is proportional to the product
    !
    !                    p_av(s,H) * detJ * phi(n,p)
    !
    ! where s and H are the surface elevation and ice thickness evaluated at that point,
    !  detJ is the determinant of the transformation linking the reference 2D element coordinates
    !  to the true coordinates at that point, and phi(n,p) is the basis function evaluated at that point.
    !
    !-----------------------------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       iCell, jCell,            &    ! i and j indices for this cell
       nNodesPerElement_2d           ! nodes per 2D (quadrilateral) element face

    character(len=*), intent(in) ::  &
       face                          ! 'north', 'south', 'east', or 'west'
 
    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex   ! x and y coordinates of vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

!    real(dp), parameter :: &
!       rhoeff = 0.5d0*rhoi*(1.d0 - rhoi/rhoo)  ! term in expression for vertical avg pressure

    real(dp), dimension(nNodesPerElement_2d) ::     &
       x, y,               & ! local coordinates of nodes
       s,                  & ! upper surface elevation at nodes
       h                     ! ice thickness at nodes

    integer, dimension(nNodesPerElement_2d) ::     &
       iNode, jNode, kNode   ! global indices of each node

    !TODO - These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at quad pts

    real(dp)  ::        &
       hqp,             & ! ice thickness at a given quadrature point (m)
       sqp,             & ! ice surface elevation at a given quadrature point (m)
       p_av,            & ! net outward pressure from ice, p_out - p_in
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    integer :: k, n, p

    ! Compute nodal geometry in a local xy reference system
    ! Note: The local y direction is really the vertical direction
    !       The local x direction depends on the face (N/S/E/W)
    ! The diagrams below show the node indexing convention, along with the true
    !  directions for each face.  The true directions are mapped to local (x,y).

    if (face=='west') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> -y

       iNode(1) = iCell-1
       jNode(1) = jCell

       iNode(2) = iCell-1
       jNode(2) = jCell-1

       x(1) = yvertex(iNode(1), jNode(1))
       x(2) = yvertex(iNode(2), jNode(2))

    elseif (face=='east') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> y

       iNode(1) = iCell
       jNode(1) = jCell-1

       iNode(2) = iCell
       jNode(2) = jCell

       x(1) = yvertex(iNode(1), jNode(1))
       x(2) = yvertex(iNode(2), jNode(2))

    elseif (face=='south') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> x

       iNode(1) = iCell-1
       jNode(1) = jCell-1

       iNode(2) = iCell
       jNode(2) = jCell-1

       x(1) = xvertex(iNode(1), jNode(1))
       x(2) = xvertex(iNode(2), jNode(2))

    elseif (face=='north') then

       !     4-----3       z
       !     |     |       ^
       !     |     |       |
       !     1-----2       ---> -x

       iNode(1) = iCell
       jNode(1) = jCell

       iNode(2) = iCell-1
       jNode(2) = jCell

       x(1) = xvertex(iNode(1), jNode(1))
       x(2) = xvertex(iNode(2), jNode(2))

    endif

    iNode(3) = iNode(2)
    jNode(3) = jNode(2)

    iNode(4) = iNode(1)
    jNode(4) = jNode(1)

    x(3) = x(2)
    x(4) = x(1)

    s(1) = stagusrf(iNode(1), jNode(1))
    s(2) = stagusrf(iNode(2), jNode(2))
    s(3) = s(2)
    s(4) = s(1)

    h(1) = stagthck(iNode(1), jNode(1))
    h(2) = stagthck(iNode(2), jNode(2))
    h(3) = h(2)
    h(4) = h(1)

    ! loop over element faces in column
    ! assume k increases from upper surface to bottom 

    do k = 1, nz-1

       ! Compute the local y coordinate (i.e., the actual z coordinate)
       y(1) = s(1) - sigma(k+1)*h(1)   ! lower left
       y(2) = s(2) - sigma(k+1)*h(2)   ! lower right
       y(3) = s(3) - sigma(k)  *h(3)   ! upper right
       y(4) = s(4) - sigma(k)  *h(4)   ! upper left

       ! Set the k index for each node
       kNode(1) = k+1
       kNode(2) = k+1
       kNode(3) = k
       kNode(4) = k

       ! loop over quadrature points

       do p = 1, nQuadPoints_2d  

          ! Compute basis function derivatives and det(J) for this quadrature point
          ! For now, pass in i, j, k, p for debugging
          !TODO - We don't actually need the derivatives, just detJ.  
          !       Should we modify the subroutine so that derivatives are optional?

          !WHL - debug
          if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Get detJ, i, j, k, p =', iCell, jCell, k, p
             print*, 'x =', x(:)
             print*, 'y =', y(:)
             print*, 'dphi_dxr_2d =', dphi_dxr_2d(:,p)
             print*, 'dphi_dyr_2d =', dphi_dyr_2d(:,p)
          endif

          call get_basis_function_derivatives_2d(nNodesPerElement_2d,                   &
                                                 x(:),              y(:),               & 
                                                 dphi_dxr_2d(:,p),  dphi_dyr_2d(:,p),   &   
                                                 dphi_dx_2d(:),     dphi_dy_2d(:),      &
                                                 detJ, iCell, jCell, k, p)

          ! For some faces, detJ is computed to be a negative number because the face is
          ! oriented opposite the x or y axis.  Fix this by taking the absolute value.

          !WHL - debug
          ! Should always have detJ < 0 for north face as it's been computed
          if (trim(face)=='north' .and. detJ > 0.d0) then
!             print*, 'north detJ > 0: i, j, k, p, detJ =', iCell, jCell, k, p, detJ
          endif

          detJ = abs(detJ)

          ! Evaluate the ice thickness and surface elevation at this quadrature point

          hqp = 0.d0
          sqp = 0.d0
          do n = 1, nNodesPerElement_2d
             hqp = hqp + phi_2d(n,p) * h(n)
             sqp = sqp + phi_2d(n,p) * s(n)
          enddo

          if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Increment shelf load vector, i, j, face, k, p =', iCell, jCell, trim(face), k, p
             print*, 'hqp, sqp =', hqp, sqp
             print*, 'detJ/vol0 =', detJ/vol0
             print*, 'grav =', grav
          endif

          ! Increment the load vector with the shelf water pressure contribution from 
          !  this quadrature point.
          ! Increment bu for east/west faces and bv for north/south faces.

          ! This formula works for ice that either is floating or is partially submerged without floating
!!          p_av = 0.5*rhoi*grav*hqp &     ! p_out
!!               - 0.5*rhoo*grav*hqp * (1.d0 - min(sqp/hqp,1.d0))**2   ! p_in

          ! This formula works for floating ice.
          p_av = 0.5*rhoi*grav*hqp * (1.d0 - rhoi/rhoo)

          if (trim(face) == 'west') then  ! net force in -x direction

             do n = 1, nNodesPerElement_2d

                bu(kNode(n),iNode(n),jNode(n)) = bu(kNode(n),iNode(n),jNode(n))    &
                                               - p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)

                if (verbose .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
                   print*, 'n, p, phi(n), delta(bu):', n, p, phi_2d(n,p), &
                           -p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
                endif

             enddo

          elseif (trim(face) == 'east') then  ! net force in x direction

             do n = 1, nNodesPerElement_2d
                bu(kNode(n),iNode(n),jNode(n)) = bu(kNode(n),iNode(n),jNode(n))    &
                                               + p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)

                if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
                   print*, 'n, p, phi(n), delta(bu):', n, p, phi_2d(n,p), &
                            p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
                endif

             enddo

          elseif (trim(face) == 'south') then  ! net force in -y direction

             do n = 1, nNodesPerElement_2d
                bv(kNode(n),iNode(n),jNode(n)) = bv(kNode(n),iNode(n),jNode(n))    &
                                               - p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)

                if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
                   print*, 'n, p, phi(n), delta(bv):', n, p, phi_2d(n,p), &
                           -p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
                endif

             enddo

          elseif (trim(face) == 'north') then  ! net force in y direction
 
             do n = 1, nNodesPerElement_2d
                bv(kNode(n),iNode(n),jNode(n)) = bv(kNode(n),iNode(n),jNode(n))    &
                                               + p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)

                if (verbose_shelf .and. this_rank==rtest .and. iCell==itest .and. jCell==jtest .and. k==ktest) then
                   print*, 'n, p, phi_2d(n,p), p_av, detJ, delta(bv):', n, p, phi_2d(n,p), p_av, detJ, &
                            p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
!WHL - debug
!!                if (verbose_shelf .and. this_rank==rtest .and. k==ktest) then
!!                   print*, 'i, j, n, p, p_av, delta(bv):', iCell, jCell, n, p, p_av, &
!!                            p_av * wqp_2d(p) * detJ/vol0 * phi_2d(n,p)
                endif

             enddo

          endif   ! face = N/S/E/W

       enddo   ! nQuadPoints_2d

    enddo   ! k (elements)

  end subroutine lateral_shelf_bc

!****************************************************************************

!  Note: This subroutine has to be called once per nonlinear iteration if
!         we're iterating on the effective viscosity.
!        Currently we recompute some geometric info that does not
!         change between nonlinear iterations.  Should we compute this info
!         once and store it?

  subroutine assemble_stiffness_matrix(nx,               ny,              &
                                       nz,               nhalo,           &
                                       sigma,            active_cell,     &
                                       xVertex,          yVertex,         &
                                       uvel,             vvel,            &
                                       stagusrf,         stagthck,        &
                                       flwafact,                          &
                                       efvs,             whichefvs,       &
                                       Auu,              Auv,             &
                                       Avu,              Avv,             &
                                       sia_factor,       ssa_factor)

    ! Assemble the stiffness matrix A and load vector b in the linear system Ax = b
 
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    real(dp), dimension(nz), intent(in) ::    &
       sigma                         ! sigma vertical coordinate

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       uvel, vvel         ! velocity components (m/yr)

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,       &  ! upper surface elevation on staggered grid (m)
       stagthck           ! ice thickness on staggered grid (m)

    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute the effective viscosity
                          ! units: Pa yr^(1/n)

    integer, intent(in) :: whichefvs      ! option for effective viscosity calculation 
                                          ! (calculate it or make it uniform)

    real(dp), dimension(nz-1,nx,ny), intent(out) ::  &
       efvs               ! effective viscosity (Pa yr)

    real(dp), dimension(27,nz,nx-1,ny-1), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    real(dp), dimension(nQuadPoints) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement,nQuadPoints) ::   &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis function, evaluated at quad pts

    !----------------------------------------------------------------
    ! Note: Kuu, Kuv, Kvu, and Kvv are 8x8 components of the stiffness matrix
    !       for the local element.  (The combined stiffness matrix is 16x16.)
    !
    ! Once these matrices are formed, their coefficients are summed into the assembled
    !  matrices Auu, Auv, Avu, Avv.  The A matrices each have as many rows as there are
    !  active nodes, but only 27 columns, corresponding to the 27 vertices that belong to
    !  the 8 elements sharing a given node.
    !
    ! The homegrown structured PCG solver works with the dense A matrices in the form
    ! computed here.  For a SLAP solver, the terms of the A matrices are put
    ! in a sparse matrix during preprocessing.
    !
    ! For a Trilinos solve, the terms of the element matrices should be sent to
    ! Trilinos a row at a time.  This capability has not been implemented as of Jan. 2013.
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement, nNodesPerElement) ::   &   !
       Kuu,          &    ! element stiffness matrix, divided into 4 parts as shown below
       Kuv,          &    !  
       Kvu,          &    !
       Kvv                !    Kuu  | Kuv
                          !    _____|____          
                          !         |
                          !    Kvu  | Kvv
                          !         
                          ! Kvu may not be needed if matrix is symmetric, but is included for now

    real(dp), dimension(nNodesPerElement) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       u, v,            & ! u and v at nodes
       s                  ! upper surface elevation at nodes

    real(dp), dimension(nQuadPoints) ::    &
       efvs_qp            ! effective viscosity at a quad pt

    real(dp) ::         &
       efvs_avg           ! efvs averaged over quad pts

    logical, parameter ::   &
       efvs_element_avg = .false.   ! if true, then average efvs over quad pts

    integer :: i, j, k, n, p
    integer :: iA, jA, kA

    integer :: iNode, jNode, kNode

    integer :: jj

    if (verbose_matrix .and. main_task) then
       print*, ' '
       print*, 'In assemble_stiffness_matrix'
    endif

    ! Initialize effective viscosity
    efvs(:,:,:) = 0.d0

    ! Initialize global stiffness matrix
    !TODO - Initialize at higher level and pass as inout?

    Auu(:,:,:,:) = 0.d0
    Auv(:,:,:,:) = 0.d0
    Avu(:,:,:,:) = 0.d0
    Avv(:,:,:,:) = 0.d0

    ! Sum over elements in active cells 
    ! Loop over all cells that border locally owned vertices.

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
     if (active_cell(i,j)) then

       do k = 1, nz-1    ! loop over elements in this column 
                         ! assume k increases from upper surface to bed

          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, 'i, j, k:', i, j, k
             print*, ' '
          endif

          ! Initialize element stiffness matrix
          Kuu(:,:) = 0.d0
          Kuv(:,:) = 0.d0
          Kvu(:,:) = 0.d0
          Kvv(:,:) = 0.d0
    
          ! compute spatial coordinates, velocity, and upper surface elevation for each node

          do n = 1, nNodesPerElement

             ! Determine (k,i,j) for this node
             ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
             ! Indices for other nodes are computed relative to this node.
             iNode = i + ishift(7,n)
             jNode = j + jshift(7,n)
             kNode = k + kshift(7,n)

             x(n) = xVertex(iNode,jNode)
             y(n) = yVertex(iNode,jNode)
             z(n) = stagusrf(iNode,jNode) - sigma(kNode)*stagthck(iNode,jNode)
             u(n) = uvel(kNode,iNode,jNode)
             v(n) = vvel(kNode,iNode,jNode)
             s(n) = stagusrf(iNode,jNode)

             if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, 'i, j, k, n, x, y, z:', i, j, k, n, x(n), y(n), z(n)
                print*, 's, u, v:', s(n), u(n), v(n)
             endif

          enddo   ! nodes per element

          ! Loop over quadrature points for this element
   
          !TODO - Think about how often to compute things.  Geometric quantities need to be
          !       computed only once per nonlinear solve, whereas efv smust be recomputed
          !       after each linear solve.

             if (verbose_state .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, ' '
                print*, 'Current uvel (m/yr): i, k =', itest, ktest
                do jj = jtest+1, jtest-1, -1
                   write(6, '(i4, 5e15.8)') jj, uvel(ktest,itest-1:itest+1,jj)
                enddo
                print*, 'Current uvel (m/yr): i, k =', itest, ktest+1
                do jj = jtest+1, jtest-1, -1
                   write(6, '(i4, 5e15.8)') jj, uvel(ktest+1,itest-1:itest+1,jj)
                enddo
                print*, ' '
                print*, 'Current vvel (m/yr): i, k =', itest, ktest
                do jj = jtest+1, jtest-1, -1
                   write(6, '(i4, 5e15.8)') jj, vvel(ktest,itest-1:itest+1,jj)
                enddo
                print*, 'Current vvel (m/yr): i, k =', itest, ktest+1
                do jj = jtest+1, jtest-1, -1
                   write(6, '(i4, 5e15.8)') jj, vvel(ktest+1,itest-1:itest+1,jj)
                enddo

                print*, ' '
                print*, 'Current uvel(lower, upper, m/yr):'
                write(6, '(2e15.8,a10,2e15.8)') u(4), u(3), '          ', u(8), u(7)
                write(6, '(2e15.8,a10,2e15.8)') u(1), u(2), '          ', u(5), u(6)
                print*, ' '
                print*, 'Current vvel(lower, upper, m/yr):'
                write(6, '(2e15.8,a10,2e15.8)') v(4), v(3), '          ', v(8), v(7)
                write(6, '(2e15.8,a10,2e15.8)') v(1), v(2), '          ', v(5), v(6)
             endif

          do p = 1, nQuadPoints

             ! Evaluate the derivatives of the element basis functions
             ! at this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

             call get_basis_function_derivatives(nNodesPerElement, &
                                                 x(:),          y(:),          z(:),           &
                                                 dphi_dxr(:,p), dphi_dyr(:,p), dphi_dzr(:,p),  &
                                                 dphi_dx(:,p),  dphi_dy(:,p),  dphi_dz(:,p),   &
                                                 detJ(p) , i, j, k, p                      )

             if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!                print*, ' '
!                print*, 'Derivatives of basis functions, p =', p
!                do n = 1, nNodesPerElement
!                   print*, 'n, dphi_dx, dphi_dy, dphi_dz:', n, dphi_dx(n,p), dphi_dy(n,p), dphi_dz(n,p)
!                enddo
             endif

          enddo   ! nQuadPoints

          do p = 1, nQuadPoints

!WHL - debug - Pass in i, j, k, and p for now

             call compute_effective_viscosity(whichefvs,        nNodesPerElement,             &
                                              dphi_dx(:,p),     dphi_dy(:,p),   dphi_dz(:,p), &
                                              u(:),             v(:),                         & 
                                              flwafact(k,i,j),  efvs_qp(p),                   &
                                              sia_factor,       ssa_factor,  i, j, k, p )

!             if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!                print*, ' '
!                print*, 'i, j, k, p, efvs (Pa yr):', i, j, k, p, efvs_qp(p)
!             endif

          enddo   ! nQuadPoints


          ! Compute average of effective viscosity over quad pts
          efvs_avg = 0.d0
          do p = 1, nQuadPoints
             efvs_avg = efvs_avg + efvs_qp(p)
          enddo
          efvs_avg = efvs_avg/nQuadPoints

          ! Fill efvs array for output
          efvs(k,i,j) = efvs_avg

          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'i, j, k, efvs_avg (Pa yr):', i, j, k, efvs_avg
          endif

          if (efvs_element_avg) then  ! apply element-average efvs instead of a different efvs for each quad pt
                                      ! less accurate but smoother?
             ! Assign average to each quad pt
             efvs_qp(:) = efvs_avg

          endif   ! efvs_element_avg

          ! Increment the element stiffness matrix with the contribution from
          ! each quadrature point.

          do p = 1, nQuadPoints

!WHL - debug - Pass in i, j, k, and p for now

             call element_matrix_blatter_pattyn(nNodesPerElement,                              & 
                                                wqp(p),           detJ(p),       efvs_qp(p),   &
                                                dphi_dx(:,p),     dphi_dy(:,p),  dphi_dz(:,p), &
                                                Kuu(:,:),         Kuv(:,:),                    &
                                                Kvu(:,:),         Kvv(:,:),                    &
                                                sia_factor,       ssa_factor,    i, j, k, p )

          enddo   ! nQuadPoints

!WHL - debug
          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Kvv: i, j, k =', i, j, k 
             do jj = 1, nNodesPerElement
                write(6,'(i4,8e18.11)') jj, Kvv(1:8,jj)
             enddo
          endif

          if (check_symmetry) then

             call check_symmetry_element_matrix(nNodesPerElement,  &
                                                Kuu, Kuv, Kvu, Kvv)

          endif

         ! If solving in parallel with Trilinos, we have the option at this point to 
         ! call a sum_into_global_matrix routine, passing one row at a time of the
         ! element matrix.  Trilinos should handle the rest.
         !

          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuu into Auu'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kuu,          Auu)

          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuv into Auv'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kuv,          Auv)

          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvu into Avu'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kvu,          Avu)

          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvv into Avv'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement,                       &
                                        i,            j,           k,           &
                                        Kvv,          Avv)

       enddo   ! nz  (loop over elements in this column)

     endif   ! active cell

    enddo      ! i
    enddo      ! j

  end subroutine assemble_stiffness_matrix

!****************************************************************************

!WHL - debug - pass in i, j, k, p for now
  subroutine get_basis_function_derivatives(nNodesPerElement, &
                                            xNode,       yNode,     zNode,    &
                                            dphi_dxr,    dphi_dyr,  dphi_dzr, &
                                            dphi_dx,     dphi_dy,   dphi_dz,  &
                                            detJ,        i, j, k, p)

    !------------------------------------------------------------------
    ! Evaluate the x, y and z derivatives of the element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element
 
    real(dp), dimension(nNodesPerElement), intent(in) :: &
                xNode, yNode, zNode,          &! nodal coordinates
                dphi_dxr, dphi_dyr, dphi_dzr   ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in reference element

    real(dp), dimension(nNodesPerElement), intent(out) :: &
                dphi_dx, dphi_dy, dphi_dz      ! derivatives of basis functions at quad pt
                                               !  wrt x, y and z in true Cartesian coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(3,3) ::  &
                Jac,      &! Jacobian matrix
                Jinv,     &! inverse Jacobian matrix
                cofactor   ! matrix of cofactors

    integer, intent(in) :: i, j, k, p

    integer :: n, row, col

!WHL - debug
    real(dp), dimension(3,3) :: prod     ! Jac * Jinv (should be identity matrix)

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !                 |                                                                          |
    !                 | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn}   sum_n{dphi_n/dxr * zn} |
    !   J(xr,yr,zr) = |                                                                          |
    !                 | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn}   sum_n{dphi_n/dyr * zn} |
    !                 |                                                                          |
    !                 | sum_n{dphi_n/dzr * xn}   sum_n{dphi_n/dzr * yn}   sum_n{dphi_n/dzr * zn} |
    !                 !                                                                          |
    !
    ! where (xn,yn,zn) are the true Cartesian nodal coordinates,
    !       (xr,yr,zr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'In get_basis_function_derivatives: i, j, k, p =', i, j, k, p
    endif

    do n = 1, nNodesPerElement
       if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'n, dphi_d(xyz)r:', n, dphi_dxr(n), dphi_dyr(n), dphi_dzr(n) 
       endif
       Jac(1,1) = Jac(1,1) + dphi_dxr(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr(n) * yNode(n)
       Jac(1,3) = Jac(1,3) + dphi_dxr(n) * zNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr(n) * yNode(n)
       Jac(2,3) = Jac(2,3) + dphi_dyr(n) * zNode(n)
       Jac(3,1) = Jac(3,1) + dphi_dzr(n) * xNode(n)
       Jac(3,2) = Jac(3,2) + dphi_dzr(n) * yNode(n)
       Jac(3,3) = Jac(3,3) + dphi_dzr(n) * zNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    cofactor(1,1) =   Jac(2,2)*Jac(3,3) - Jac(2,3)*Jac(3,2)
    cofactor(1,2) = -(Jac(2,1)*Jac(3,3) - Jac(2,3)*Jac(3,1))
    cofactor(1,3) =   Jac(2,1)*Jac(3,2) - Jac(2,2)*Jac(3,1)
    cofactor(2,1) = -(Jac(1,2)*Jac(3,3) - Jac(1,3)*Jac(3,2))
    cofactor(2,2) =   Jac(1,1)*Jac(3,3) - Jac(1,3)*Jac(3,1)
    cofactor(2,3) = -(Jac(1,1)*Jac(3,2) - Jac(1,2)*Jac(3,1))
    cofactor(3,1) =   Jac(1,2)*Jac(2,3) - Jac(1,3)*Jac(2,2)
    cofactor(3,2) = -(Jac(1,1)*Jac(2,3) - Jac(1,3)*Jac(2,1))
    cofactor(3,3) =   Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    detJ = Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'detJ1:', Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)
       print*, 'detJ2:', Jac(2,1)*cofactor(2,1) + Jac(2,2)*cofactor(2,2) + Jac(2,3)*cofactor(2,3)
       print*, 'detJ3:', Jac(3,1)*cofactor(3,1) + Jac(3,2)*cofactor(3,2) + Jac(3,3)*cofactor(3,3)
    endif

    if (abs(detJ) > 0.d0) then
       do col = 1, 3
          do row = 1, 3
             Jinv(row,col) = cofactor(col,row)
          enddo
       enddo
       Jinv(:,:) = Jinv(:,:) / detJ
    else
       !WHL - do a proper abort here
       print*, 'stopping, det J = 0'
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:) 
       stop
    endif

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Jacobian calc, p =', p
       print*, 'det J =', detJ
       print*, ' '
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, Jac(3,:)
       print*, ' '
       print*, 'cofactor matrix:'
       print*, cofactor(1,:)
       print*, cofactor(2,:)
       print*, cofactor(3,:)
       print*, ' '
       print*, 'Inverse matrix:'
       print*, Jinv(1,:)
       print*, Jinv(2,:)
       print*, Jinv(3,:)
       print*, ' '
       prod = matmul(Jac, Jinv)
       print*, 'Jac*Jinv:'
       print*, prod(1,:)
       print*, prod(2,:)
       print*, prod(3,:)
    endif

    ! bug check - Verify that J * Jinv = I

    prod = matmul(Jac,Jinv)
    do col = 1, 3
       do row = 1, 3
          if (abs(prod(row,col) - identity3(row,col)) > 1.d-12) then
             !TODO - do a proper abort here 
             print*, 'stopping, Jac * Jinv /= identity'
             print*, 'i, j, k, p:', i, j, k, p
             print*, 'Jac*Jinv:'
             print*, prod(1,:)
             print*, prod(2,:)
             print*, prod(3,:)
             stop
          endif
       enddo
    enddo

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           |          |            | 
    !   | dphi_n/dy | = Jinv * | dphi_n/dyr |
    !   |           |          |            |
    !   | dphi_n/dz |          | dphi_n/dzr |
    !
    !------------------------------------------------------------------

    dphi_dx(:) = 0.d0
    dphi_dy(:) = 0.d0
    dphi_dz(:) = 0.d0

    !TODO - Don't need first terms on RHS since it's always zero?
    do n = 1, nNodesPerElement
       dphi_dx(n) = dphi_dx(n) + Jinv(1,1)*dphi_dxr(n)  &
                               + Jinv(1,2)*dphi_dyr(n)  &
                               + Jinv(1,3)*dphi_dzr(n)
       dphi_dy(n) = dphi_dy(n) + Jinv(2,1)*dphi_dxr(n)  &
                               + Jinv(2,2)*dphi_dyr(n)  &
                               + Jinv(2,3)*dphi_dzr(n)
       dphi_dz(n) = dphi_dz(n) + Jinv(3,1)*dphi_dxr(n)  &
                               + Jinv(3,2)*dphi_dyr(n)  &
                               + Jinv(3,3)*dphi_dzr(n)

       if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'n, dphi_d(xyz):', n, dphi_dx(n), dphi_dy(n), dphi_dz(n) 
       endif

    enddo

    ! debug - Check that the sum of dphi_dx, etc. is close to zero  
    !TODO - Think about why this should be the case.
    !       Use a relative rather than an absolute threshold?

    if (abs( sum(dphi_dx)/maxval(dphi_dx) ) > 1.d-11) then
        print*, 'stopping, sum over basis functions of dphi_dx > 0'
        print*, 'dphi_dx =', dphi_dx(:)
        print*, 'sum =', sum(dphi_dx)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs( sum(dphi_dy)/maxval(dphi_dy) ) > 1.d-11) then
        print*, 'stopping, sum over basis functions of dphi_dy > 0'
        print*, 'dphi_dy =', dphi_dy(:)
        print*, 'sum =', sum(dphi_dy)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs( sum(dphi_dz)/maxval(dphi_dz) ) > 1.d-11) then
        print*, 'stopping, sum over basis functions of dphi_dz > 0'
        print*, 'dphi_dz =', dphi_dz(:)
        print*, 'sum =', sum(dphi_dz)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

  end subroutine get_basis_function_derivatives

!****************************************************************************

!WHL - Pass i, j, k, and p for now

  subroutine compute_effective_viscosity (whichefvs,   nNodesPerElement,       &
                                          dphi_dx,     dphi_dy,  dphi_dz,      &
                                          uvel,        vvel,                   &
                                          flwafact,    efvs,                   &
                                          sia_factor,  ssa_factor,  i, j, k, p )

    ! Compute effective viscosity at a quadrature point, based on the latest
    ! guess for the velocity field

    integer, intent(in) :: i, j, k, p

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: &
       whichefvs       ! method for computing effective viscosity
                       ! 0 = based on effective strain rate
                       ! 1 = constant value
                       !TODO - Let 0 = constant value, 1 = eff str rate?
                       ! Add an option representing a true constant value, 
                       ! rather than temperature-based?

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       dphi_dx, dphi_dy, dphi_dz         ! derivatives of basis functions,
                                         ! evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       uvel, vvel      ! current guess for velocity at each node of element (m/yr)

    real(dp), intent(in) ::  &
       flwafact        ! temperature-based flow factor for this element, 0.5 * A^{-1/n}
                       ! units: Pa yr^{1/n}

    real(dp), intent(out) ::   &
       efvs            ! effective viscosity at this quadrature point (Pa yr)
                       ! computed as 0.5 * A^{-1/n) * effstrain^{(1-n)/n)}
                       
    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !----------------------------------------------------------------
    ! Local parameters
    !----------------------------------------------------------------

!TODO - Test sensitivity of model convergence to this parameter
    real(dp), parameter ::   &
!!       effstrain_min = 1.d-20,          &! minimum value of effective strain rate, s^{-1}
!!       effstrain_min = 1.d-20*scyr,     &! minimum value of effective strain rate, yr^{-1}
                                         ! GLAM uses 1.d-20 s^{-1} for minimum effective strain rate
       effstrain_min = 1.d-8,     &! minimum value of effective strain rate, yr^{-1}
                                   ! Mauro suggests 1.d-8 yr^{-1}
       p_effstr = (1.d0 - real(gn,dp)) / real(gn,dp)    ! exponent (1-n)/n in effective viscosity relation
                                                               
    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::            &
       du_dx, dv_dx,       & ! strain rate components
       du_dy, dv_dy,       &
       du_dz, dv_dz,       &
       effstrain,          & ! effective strain rate, s^{-1}
       effstrainsq           ! square of effective strain rate
        
    integer :: n


!WHL - debug
    if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Compute efvs, i, j, k, p =', i, j, k, p
    endif

    select case(whichefvs)

    case(HO_EFVS_CONSTANT)

!WHL - debug
       if (trial_efvs) then

          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, 'Computing trial efvs'
             print*, ' '
          endif

          du_dx = 0.d0
          dv_dx = 0.d0
          du_dy = 0.d0
          dv_dy = 0.d0
          du_dz = 0.d0
          dv_dz = 0.d0

          ! Compute strain rate components at this quadrature point
          ! by summing over basis functions

          do n = 1, nNodesPerElement

             du_dx = du_dx + dphi_dx(n)*uvel(n)
             dv_dx = dv_dx + dphi_dx(n)*vvel(n)

             du_dy = du_dy + dphi_dy(n)*uvel(n)
             dv_dy = dv_dy + dphi_dy(n)*vvel(n)

             du_dz = du_dz + dphi_dz(n)*uvel(n)
             dv_dz = dv_dz + dphi_dz(n)*vvel(n)

          enddo   ! nNodesPerElement

          ! Compute effective strain rate at this quadrature point (PGB 2012, eq. 3 and 9)
          ! Units are yr^(-1)

          effstrainsq = effstrain_min**2                                      &
                      + ssa_factor * (du_dx**2 + dv_dy**2 + du_dx*dv_dy       &
                                      + 0.25d0 * (dv_dx + du_dy)**2)          &
                      + sia_factor * 0.25d0 * (du_dz**2 + dv_dz**2)

          effstrain = sqrt(effstrainsq)

          ! Compute effective viscosity (PGB 2012, eq. 4)
          ! Units: flwafact has units Pa yr^{1/n}
          !        effstrain has units yr^{-1}
          !        p_effstr = (1-n)/n 
          !                 = -2/3 for n=3
          ! Thus efvs has units Pa yr
 
          efvs = flwafact * effstrain**p_effstr

          if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!             print*, 'effstrain_min (yr-1)=', effstrain_min
!             print*, 'Trial du/dx, du/dy, du/dz (yr-1) =', du_dx, du_dy, du_dz
             print*, 'Trial dv/dx, dv/dy, dv/dz (yr-1) =', dv_dx, dv_dy, dv_dz
             print*, 'Trial flwafact, effstrain (yr-1), efvs(Pa yr) =', flwafact, effstrain, efvs
          endif

       endif  ! trial_efvs
!WHL - end debug (next the real thing)

       efvs = 1.d7      ! Steve Price recommends 10^7 Pa yr
                        ! (~3e14 Pa s)
!WHL - This is the glam-type scaling
!!       efvs = efvs * scyr/tim0 / tau0   ! tau0 = rhoi*grav*thk0

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'Set efvs = constant (Pa yr):', efvs
       endif


    case(HO_EFVS_FLOWFACT)      ! set the effective viscosity to a multiple of the flow factor, 0.5*A^(-1/n)
                 
                 !SCALING: Set the effective strain rate (s^{-1}) based on typical 
                 !         velocity and length scales, to agree with GLAM.
  
       effstrain = vel_scale/len_scale   ! typical strain rate, s^{-1}  
       efvs = flwafact * effstrain**p_effstr  

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, 'flwafact, effstrain (yr-1), efvs (Pa yr)=', flwafact, effstrain, efvs
       endif

    case(HO_EFVS_NONLINEAR)    ! calculate effective viscosity based on effective strain rate, n = 3

       du_dx = 0.d0
       dv_dx = 0.d0
       du_dy = 0.d0
       dv_dy = 0.d0
       du_dz = 0.d0
       dv_dz = 0.d0

       ! Compute strain rate components at this quadrature point
       ! by summing over basis functions

       do n = 1, nNodesPerElement

          du_dx = du_dx + dphi_dx(n)*uvel(n)
          dv_dx = dv_dx + dphi_dx(n)*vvel(n)

          du_dy = du_dy + dphi_dy(n)*uvel(n)
          dv_dy = dv_dy + dphi_dy(n)*vvel(n)

          du_dz = du_dz + dphi_dz(n)*uvel(n)
          dv_dz = dv_dz + dphi_dz(n)*vvel(n)

       enddo   ! nNodesPerElement

       ! Compute effective strain rate at this quadrature point (PGB 2012, eq. 3 and 9)
       ! Units are yr^(-1)

       effstrainsq = effstrain_min**2                                      &
                   + ssa_factor * (du_dx**2 + dv_dy**2 + du_dx*dv_dy       &
                                   + 0.25d0 * (dv_dx + du_dy)**2)          &
                   + sia_factor * 0.25d0 * (du_dz**2 + dv_dz**2)

       effstrain = sqrt(effstrainsq)

       ! Compute effective viscosity (PGB 2012, eq. 4)
       ! Units: flwafact has units Pa yr^{1/n}
       !        effstrain has units yr^{-1}
       !        p_effstr = (1-n)/n 
       !                 = -2/3 for n=3
       ! Thus efvs has units Pa yr
 
       efvs = flwafact * effstrain**p_effstr

       if (verbose_efvs .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
!          print*, 'effstrain_min (yr-1)=', effstrain_min
          print*, 'du/dx, du/dy, du/dz (yr-1) =', du_dx, du_dy, du_dz
          print*, 'dv/dx, dv/dy, dv/dz (yr-1) =', dv_dx, dv_dy, dv_dz
          print*, 'flwafact, effstrain (yr-1), efvs(Pa yr) =', flwafact, effstrain, efvs
       endif

    end select

  end subroutine compute_effective_viscosity

!****************************************************************************

!WHL - debug - Pass in i, j, k, and p for now

  subroutine element_matrix_blatter_pattyn(nNodesPerElement,                &
                                           wqp,         detJ,     efvs,     &
                                           dphi_dx,     dphi_dy,  dphi_dz,  &
                                           Kuu,         Kuv,                &
                                           Kvu,         Kvv,                &
                                           sia_factor,  ssa_factor,  ii, jj, k, p)

    !------------------------------------------------------------------
    ! Increment the stiffness matrices Kuu, Kuv, Kvu, Kvv with the
    ! contribution from a particular quadrature point, 
    ! based on the Blatter-Pattyn first-order equations.
    !
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

!WHL - debug - Pass in i, j, k, and p for now
    integer, intent(in) :: ii, jj, k, p

    integer, intent(in) :: nNodesPerElement  ! number of nodes per element

    real(dp), intent(in) ::    &
             wqp,        &! weight for this quadrature point
             detJ,       &! determinant of Jacobian for the transformation
                          !  between the reference element and true element
             efvs         ! effective viscosity at this quadrature point

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
             dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis functions,
                                         ! evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(inout) :: &
             Kuu, Kuv, Kvu, Kvv     ! components of element stiffness matrix

!WHL - debug (Make intent(in) again after testing)
!    real(dp), intent(in) ::  &
    real(dp) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

!WHL - debug
!    sia_factor = 0.d0
!    ssa_factor = 0.d0
 
!WHL - debug
    if (verbose_matrix .and. this_rank==rtest .and. ii==itest .and. jj==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Increment element matrix, i, j, k, p =', ii, jj, k, p
    endif

    ! Increment the element stiffness matrices for the first-order Blatter-Pattyn equations.
    ! The terms in parentheses can be derived from PGB 2012, eq. 13 and 15.
    ! The factor of 2 in front of efvs has been absorbed into the quantities in parentheses.
    
    do j = 1, nNodesPerElement      ! columns of K
       do i = 1, nNodesPerElement   ! rows of K

!       if (verbose_matrix .and. this_rank==rtest .and. ii==itest .and. jj==jtest .and. k==ktest .and. p==ptest) then
!          print*, 'efvs, wqp, detJ/vol0 =', efvs, wqp, detJ/vol0
!          print*, 'dphi_dz(1) =', dphi_dz(1)
!          print*, 'dphi_dx(1) =', dphi_dx(1)
!          print*, 'Kuu dphi/dz increment(1,1) =', efvs*wqp*detJ/vol0*dphi_dz(1)*dphi_dz(1)
!          print*, 'Kuu dphi/dx increment(1,1) =', efvs*wqp*detJ/vol0*4.d0*dphi_dx(1)*dphi_dx(1)
!       endif

       if (verbose_matrix .and. this_rank==rtest .and. ii==itest .and. jj==jtest .and. k==ktest & !.and. p==ptest &
                          .and. i==Krowtest .and. j==Kcoltest) then
          print*, 'irow, jcol =', i, j
!          print*, 'efvs, wqp, detJ/vol0 =', efvs, wqp, detJ/vol0
          print*, 'dphi_dx, dphi_dy, dphi_dz(irow) =', dphi_dx(i), dphi_dy(i), dphi_dz(i)
          print*, 'dphi_dx, dphi_dy, dphi_dz(jcol) =', dphi_dx(j), dphi_dy(j), dphi_dz(j)
          print*, 'Kuu SSA increment(irow,jcol) =', efvs*wqp*detJ/vol0*(4.d0*dphi_dx(i)*dphi_dx(j) + dphi_dy(j)*dphi_dy(i))
          print*, 'Kvv SSA increment(irow,jcol) =', efvs*wqp*detJ/vol0*(4.d0*dphi_dy(i)*dphi_dy(j) + dphi_dx(j)*dphi_dx(i))
          print*, 'SIA increment(irow,jcol) =', efvs*wqp*detJ/vol0*dphi_dz(i)*dphi_dz(j)
       endif


          !WHL - Note volume scaling such that detJ/vol0 is closer to unity

          Kuu(i,j) = Kuu(i,j) + efvs * wqp * detJ/vol0 *                                         &
                              ( ssa_factor * (4.d0*dphi_dx(i)*dphi_dx(j) + dphi_dy(i)*dphi_dy(j))  &
                              + sia_factor * (dphi_dz(i)*dphi_dz(j)) )

          Kuv(i,j) = Kuv(i,j) + efvs * wqp * detJ/vol0 *                                         &
                                ssa_factor * (2.d0*dphi_dx(i)*dphi_dy(j) + dphi_dy(i)*dphi_dx(j))

          Kvu(i,j) = Kvu(i,j) + efvs * wqp * detJ/vol0 *                                         &
                                ssa_factor * (2.d0*dphi_dy(i)*dphi_dx(j) + dphi_dx(i)*dphi_dy(j))

          Kvv(i,j) = Kvv(i,j) + efvs * wqp * detJ/vol0 *                                         &
                              ( ssa_factor * (4.d0*dphi_dy(i)*dphi_dy(j) + dphi_dx(i)*dphi_dx(j))  &
                              + sia_factor * (dphi_dz(i)*dphi_dz(j)) )

       enddo  ! i (rows)
    enddo     ! j (columns)

  end subroutine element_matrix_blatter_pattyn

!****************************************************************************

  subroutine element_to_global_matrix(nx,           ny,          nz,          &
                                      nNodesPerElement,                       &
                                      iElement,     jElement,    kElement,    &
                                      Kmat,         Amat)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! So this subroutine is called four times per element.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    integer, intent(in) ::   &
       nNodesPerElement         ! number of nodes per element

    integer, intent(in) ::   &
       iElement, jElement, kElement     ! i, j and k indices for this element

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: i, j, k, m
    integer :: iA, jA, kA
    integer :: n, nr, nc

!WHL - debug
    if (verbose_matrix .and. this_rank==rtest .and. iElement==itest .and. jElement==jtest .and. kElement==ktest) then
       print*, 'Element i, j, k:', iElement, jElement, kElement 
       print*, 'Rows of K:'
       do n = 1, nNodesPerElement
          write(6, '(8e12.4)') Kmat(n,:)
       enddo
    endif


!TODO - Switch loops or switch order of indices in K(nr,nc)?
!       Inner index nr is currently the outer loop

    do nr = 1, nNodesPerElement       ! rows of K

       ! Determine (k,i,j) for this node
       ! The reason for the '7' is that node 7, in the NE corner of the upper layer, has index (k,i,j).
       ! Indices for other nodes are computed relative to this node.
       i = iElement + ishift(7,nr)
       j = jElement + jshift(7,nr)
       k = kElement + kshift(7,nr)
      
       do nc = 1, nNodesPerElement    ! columns of K

          kA = kshift(nr,nc)           ! k index of A into which K(m,n) is summed
          iA = ishift(nr,nc)           ! similarly for i and j indices 
          jA = jshift(nr,nc)           ! these indices can take values -1, 0 and 1

          m = indxA(iA,jA,kA)
          Amat(m,k,i,j) = Amat(m,k,i,j) + Kmat(nr,nc)

!WHL - debug
          if (verbose_matrix .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest  &
                      .and. iA == iAtest .and. jA==jAtest .and. kA==kAtest) then
             print*, ' '
             print*, 'i, j, k, iA, jA, kA:', i, j, k, iA, jA, kA
             print*, 'i, j, k of element:', iElement, jElement, kElement
             print*, 'nr, nc, Kmat, new Amat:', nr, nc, Kmat(nr,nc), Amat(m,k,i,j)
          endif

       enddo     ! nc

    enddo        ! nr

  end subroutine element_to_global_matrix

!****************************************************************************
!WHL -  Added basal BC here.  Start with linear sliding law.
!       Similar to lateral shelf BC in terms of doing surface integrals.
!       Integrate over all faces that contain at least one node with a stress BC. (Not Dirichlet or free-slip)
!       Dirichlet BCs are enforced after matrix assembly. 


  subroutine basal_sliding_bc(nx,               ny,              &
                              nz,               nhalo,           &
                              nNodesPerElement_2d,               &
                              active_cell,      beta,            &
                              xVertex,          yVertex,         &
                              Auu,              Avv)

    !------------------------------------------------------------------------
    ! Increment the Auu and Avv matrices with basal traction terms.
    !
    ! For now, assume a linear sliding law:
    !   tau_x = -beta*u
    !   tau_y = -beta*v
    ! 
    ! where beta is defined at vertices.
    !------------------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
                                     ! Note: the number of elements per column is nz-1
       nhalo                         ! number of halo layers

    integer, intent(in) :: nNodesPerElement_2d

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if cell contains ice and borders a locally owned vertex

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       beta                          ! basal traction field (Pa/(m/yr)) at cell vertices

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       xVertex, yVertex     ! x and y coordinates of vertices

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::  &
       Auu, Avv             ! parts of stiffness matrix

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k, n, p

    integer :: iNode, jNode, kNode   ! global indices for nodes

    real(dp), dimension(nNodesPerElement_2d) ::   &
       x, y,        & ! Cartesian coordinates of nodes
       b              ! beta at nodes

    !TODO - These are not currently used except as dummy arguments
    real(dp), dimension(nNodesPerElement_2d) ::   &
       dphi_dx_2d, dphi_dy_2d    ! derivatives of basis functions, evaluated at quad pts

    real(dp) ::   &
       beta_qp,     & ! beta evaluated at quadrature point
       detJ           ! determinant of Jacobian for the transformation
                      !  between the reference element and true element

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d) ::   &
       Kuu, Kvv       ! components of element matrix associated with basal sliding
        
    ! Sum over elements in active cells 
    ! Loop over all cells that contain locally owned vertices
    !TODO - Make sure we don't step out of bounds
    !       With more care, we could skip some computations for vertices that are not locally owned.

    !WHL- debug
    if (verbose_basal .and. this_rank==rtest) then 
       print*, 'In basal_sliding_bc: itest, jtest, rank =', itest, jtest, rtest
    endif

    do j = nhalo+1, ny-nhalo+1
    do i = nhalo+1, nx-nhalo+1
       
       !TODO - Should we exclude cells that have Dirichlet basal BCs for all vertices?

       if (active_cell(i,j)) then   ! ice is present

          ! Set x and y for each node

          !     4-----3       y
          !     |     |       ^
          !     |     |       |
          !     1-----2       ---> x

          x(1) = xVertex(i-1,j-1)
          x(2) = xVertex(i,j-1)
          x(3) = xVertex(i,j)
          x(4) = xVertex(i-1,j)

          y(1) = yVertex(i-1,j-1)
          y(2) = yVertex(i,j-1)
          y(3) = yVertex(i,j)
          y(4) = yVertex(i-1,j)

          b(1) = beta(i-1,j-1)          
          b(2) = beta(i,j-1)
          b(3) = beta(i,j)
          b(4) = beta(i-1,j)

          k = nz     ! basal layer only

          ! loop over quadrature points

          do p = 1, nQuadPoints_2d  

             ! Compute basis function derivatives and det(J) for this quadrature point
             ! For now, pass in i, j, k, p for debugging
             !TODO - We don't actually need the derivatives, just detJ.  
             !       Should we modify the subroutine so that derivatives are optional?

             call get_basis_function_derivatives_2d(nNodesPerElement_2d,                  &
                                                    x(:),       y(:),                     & 
                                                    dphi_dxr_2d(:,p), dphi_dyr_2d(:,p),   &   
                                                    dphi_dx_2d(:),    dphi_dy_2d(:),      &
                                                    detJ, i, j, k, p)
          
             ! Evaluate beta at this quadrature point
 
             beta_qp = 0.d0
             do n = 1, nNodesPerElement_2d
                beta_qp = beta_qp + phi_2d(n,p) * b(n)
             enddo

             if (verbose_basal .and. this_rank==rtest .and. i==itest .and. j==jtest) then
                print*, ' '
                print*, 'Increment basal traction, i, j, p =', i, j, p
                print*, 'beta_qp =', beta_qp
                print*, 'detJ/vol0 =', detJ/vol0
             endif

             call element_matrix_basal_sliding(nNodesPerElement_2d,       &
                                               wqp_2d(p),   detJ,         & 
                                               beta_qp,     phi_2d(:,p),  &
                                               Kuu(:,:),                  &
                                               i, j, k, p)

             Kvv(:,:) = Kuu(:,:)   ! TODO: Is this true for more general sliding laws?

             ! Insert terms of element matrix into global matrices Auu and Avv

             call element_to_global_matrix_2d(nx,    ny,    nz,             &
                                              nNodesPerElement_2d,          &
                                              i,     j,    k,               &
                                              Kuu,   Auu)

             call element_to_global_matrix_2d(nx,    ny,    nz,             &
                                              nNodesPerElement_2d,          &
                                              i,     j,    k,               &
                                              Kvv,   Avv)

          enddo   ! nQuadPoints_2d

       endif      ! active_cell

    enddo         ! i
    enddo         ! j

  end subroutine basal_sliding_bc

!****************************************************************************
!WHL - Pass in i, j, k, p for now

  subroutine get_basis_function_derivatives_2d(nNodesPerElement_2d,     &
                                               xNode,       yNode,      &
                                               dphi_dxr,    dphi_dyr,   &
                                               dphi_dx,     dphi_dy,    &
                                               detJ, i, j, k, p)

!TODO - Change subroutine name?

    !------------------------------------------------------------------
    ! Evaluate the x and y derivatives of 2D element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 2D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement_2d   ! number of nodes per element
 
    real(dp), dimension(nNodesPerElement_2d), intent(in) :: &
                xNode, yNode,                   &! nodal coordinates
                dphi_dxr, dphi_dyr               ! derivatives of basis functions at quad pt
                                                 !  wrt x and y in reference element

    real(dp), dimension(nNodesPerElement_2d), intent(out) :: &
                dphi_dx , dphi_dy                ! derivatives of basis functions at quad pt
                                                 !  wrt x and y in true Cartesian coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(2,2) ::  &
                Jac,      &! Jacobian matrix
                Jinv       ! inverse Jacobian matrix

    integer, intent(in) :: i, j, k, p

    integer :: n, row, col

!WHL - debug
    real(dp), dimension(2,2) :: prod     ! Jac * Jinv (should be identity matrix)

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !              |                                                  |
    !              | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn}  |
    !   J(xr,yr) = |                                                  |
    !              | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn}  |
    !              |                                                  |
    !
    ! where (xn,yn) are the true Cartesian nodal coordinates,
    !       (xr,yr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'In get_basis_function_derivatives_2d: i, j, k, p =', i, j, k, p
    endif

    do n = 1, nNodesPerElement_2d
       if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'n, x, y:', n, xNode(n), yNode(n)
          print*, 'dphi_dxr, dphi_dyr:', dphi_dxr(n), dphi_dyr(n) 
       endif
       Jac(1,1) = Jac(1,1) + dphi_dxr(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr(n) * yNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr(n) * yNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    if (abs(detJ) > 0.d0) then
       Jinv(1,1) =  Jac(2,2)/detJ
       Jinv(1,2) = -Jac(1,2)/detJ
       Jinv(2,1) = -Jac(2,1)/detJ
       Jinv(2,2) =  Jac(1,1)/detJ
    else
       !WHL - do a proper abort here
       print*, 'stopping, det J = 0'
       print*, 'i, j, k, p:', i, j, k, p
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       stop
    endif

    if (verbose_Jac .and. this_rank==rtest .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Jacobian calc, p =', p
       print*, 'det J =', detJ
       print*, ' '
       print*, 'Jacobian matrix:'
       print*, Jac(1,:)
       print*, Jac(2,:)
       print*, ' '
       print*, 'Inverse matrix:'
       print*, Jinv(1,:)
       print*, Jinv(2,:)
       print*, ' '
       prod = matmul(Jac, Jinv)
       print*, 'Jac*Jinv:'
       print*, prod(1,:)
       print*, prod(2,:)
    endif

    ! bug check - Verify that J * Jinv = I

    prod = matmul(Jac,Jinv)
    do col = 1, 2
       do row = 1, 2
          if (abs(prod(row,col) - identity3(row,col)) > 1.d-12) then
             !TODO - do a proper abort here 
             print*, 'stopping, Jac * Jinv /= identity'
             print*, 'i, j, k, p:', i, j, k, p
             print*, 'Jac*Jinv:'
             print*, prod(1,:)
             print*, prod(2,:)
             stop
          endif
       enddo
    enddo

!TODO - Make this part optional?

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           | = Jinv * |            |
    !   | dphi_n/dy |          | dphi_n/dyr |
    !
    !------------------------------------------------------------------

    dphi_dx(:) = 0.d0
    dphi_dy(:) = 0.d0

    do n = 1, nNodesPerElement_2d
       dphi_dx(n) = dphi_dx(n) + Jinv(1,1)*dphi_dxr(n)  &
                               + Jinv(1,2)*dphi_dyr(n)
       dphi_dy(n) = dphi_dy(n) + Jinv(2,1)*dphi_dxr(n)  &
                               + Jinv(2,2)*dphi_dyr(n)
    enddo

    ! debug - Check that the sum of dphi_dx, etc. is close to zero  
    !TODO - Think about why this should be the case

    if (abs( sum(dphi_dx)/maxval(dphi_dx) ) > 1.d-11) then
        print*, 'stopping, sum over basis functions of dphi_dx > 0'
        print*, 'dphi_dx =', dphi_dx(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs( sum(dphi_dy)/maxval(dphi_dy) ) > 1.d-11) then
        print*, 'stopping, sum over basis functions of dphi_dy > 0'
        print*, 'dphi_dy =', dphi_dy(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

  end subroutine get_basis_function_derivatives_2d

!****************************************************************************

!WHL - debug - Pass in i, j, k, and p for now

  subroutine element_matrix_basal_sliding(nNodesPerElement_2d,   &
                                          wqp_2d,      detJ,     &
                                          beta,        phi_2d,   &
                                          Kuu,                   &
                                          i, j, k, p)

    !------------------------------------------------------------------
    ! Increment the matrix K corresponding to a 2D element face with the contribution
    ! from a particular quadrature point, given a linear sliding law.
    !
    ! For this sliding law, Kuu = Kvv, and Kvu = Kuv = 0.
    ! So it suffices here to compute Kuu.
    !
    ! TODO: Modify for more general sliding laws.
    !------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

!WHL - debug - Pass in i, j, k, and p for now
    integer, intent(in) :: i, j, k, p

    integer, intent(in) :: nNodesPerElement_2d  ! number of nodes per 2D element face

    real(dp), intent(in) ::    &
             wqp_2d,     &! weight for this quadrature point
             detJ,       &! determinant of Jacobian for the transformation
                          !  between the reference element and true element
             beta         ! basal traction factor beta at this quadrature point

    real(dp), dimension(nNodesPerElement_2d), intent(in) ::  &
             phi_2d       ! 2D basis functions phi evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d), intent(out) :: &
             Kuu          ! components of stiffness matrix for this element face

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: nc, nr

    if (verbose_basal .and. this_rank==rtest .and. i==itest .and. j==jtest) then
       print*, ' '
       print*, 'Increment element matrix for basal BC, p =', p
       print*, 'phi_2d:', phi_2d(1:4)
       print*, ' '
       print*, 'phi_2d(nr)*phi_2d(nc):'
       do nr = 1, nNodesPerElement_2d   ! rows of K
          print*, phi_2d(nr)*phi_2d(1:4)
       enddo       
    endif

    ! Increment the stiffness matrix for this quadrature point.
    
    Kuu(:,:) = 0.d0

    do nc = 1, nNodesPerElement_2d      ! columns of K
       do nr = 1, nNodesPerElement_2d   ! rows of K

          !WHL - Note volume scaling

          Kuu(nr,nc) = Kuu(nr,nc) + beta * wqp_2d * detJ/vol0 * phi_2d(nr)*phi_2d(nc)

       enddo  ! m (rows)
    enddo     ! n (columns)

  end subroutine element_matrix_basal_sliding

!****************************************************************************

  subroutine element_to_global_matrix_2d(nx,           ny,          nz,          &
                                         nNodesPerElement_2d,                    &
                                         iElement,     jElement,    kLayer,      &
                                         Kmat,         Amat)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! For basal BC, the contribution to Kuv and Kvu is typically zero.
    ! So this subroutine is called twice per 2D element face.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    integer, intent(in) ::   &
       nNodesPerElement_2d      ! number of nodes per 2D element face

    integer, intent(in) ::   &
       iElement, jElement,   &  ! i and j indices for this element
       kLayer                   ! k index for the layer corresponding to this 2D face
                                ! kLayer = nz for basal BC

    real(dp), dimension(nNodesPerElement_2d,nNodesPerElement_2d), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: i, j, k, m
    integer :: iA, jA, kA
    integer :: nr, nc

    if (verbose_basal .and. this_rank==rtest .and. iElement==itest .and. jElement==jtest) then
       print*, ' '
       print*, 'Basal BC: K matrix:'
       do nr = 1, nNodesPerElement_2d
          write(6, '(4e14.6)') Kmat(nr,:)
       enddo
    endif

!TODO - Switch loops or switch order of indices in K(m,n)?
!       Inner index m is currently the outer loop

    k = kLayer   ! all nodes lie in one horizontal layer, typically k = nz for basal BC
    kA = 0       ! k index of A into which K(nr,nc) is summed
                 ! = 0 since all nodes lie in the same layer

    do nr = 1, nNodesPerElement_2d     ! rows of K

       ! Determine (i,j) for this node
       ! The reason for the '3' is that node 3, in the NE corner of the cell, has horizontal indices (i,j).
       ! Indices for other nodes are computed relative to this node.

       i = iElement + ishift(3,nr)
       j = jElement + jshift(3,nr)
      
       do nc = 1, nNodesPerElement_2d ! columns of K

          iA = ishift(nr,nc)          ! iA index of A into which K(nr,nc) is summed
          jA = jshift(nr,nc)          ! similarly for jA
                                      ! these indices can take values -1, 0 and 1

          m = indxA(iA,jA,kA)
          Amat(m,k,i,j) = Amat(m,k,i,j) + Kmat(nr,nc)

       enddo     ! n

    enddo        ! m

  end subroutine element_to_global_matrix_2d

!****************************************************************************

  subroutine dirichlet_boundary_conditions(nx,  ny,  nz,    nhalo,            &
                                           active_vertex,   umask_dirichlet,  &
                                           Auu,             Auv,              &
                                           Avu,             Avv,              &
                                           bu,              bv)

    !----------------------------------------------------------------
    ! Modify the global matrix and RHS for Dirichlet boundary conditions.
    ! This subroutine assumes that u = v = 0 at Dirichlet points.
    ! For each such point, we find the corresponding row and column of the global matrix.
    ! We zero out each row and column, except for setting the diagonal term to 1.
    ! We set the corresponding rhs to 0, so that the solution must be 0.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::  &
       active_vertex       ! true for active vertices (vertices of active cells)

    logical, dimension(nz,nx-1,ny-1), intent(in) ::  &
       umask_dirichlet     ! Dirichlet mask for velocity (if true, u = 0)

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv,    &      ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       bu, bv              ! assembled load vector, divided into 2 parts

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
    
    integer :: i, j, k     ! Cartesian indices of nodes
    integer :: iA, jA, kA  ! i, j, and k offsets of neighboring nodes 
    integer :: m

!PARALLEL - Make sure not to step out of bounds here.
!           This requires nhalo >= 2.

! Loop over all vertices that border locally owned vertices.
! Locally owned vertices are (nhalo+1:nx-nhalo, nhalo+1:ny-nhalo)
!           
     do j = nhalo, ny-nhalo+1
        do i = nhalo, nx-nhalo+1
          if (active_vertex(i,j)) then
             do k = 1, nz
                if (umask_dirichlet(k,i,j)) then

                   !WHL - debug
!!                   print*, 'umask_dirichlet: i, j, k =', i, j, k

                   ! loop through potential nonzero matrix values in the row associated with this vertex
                   ! Note: row = NodeID(k,i,j) if we were forming a matrix with one row per node

                   do kA = -1,1
                   do jA = -1,1
                   do iA = -1,1

                      if ( (k+kA >= 1 .and. k+kA <= nz)         &
                                      .and.                     &
                           (i+iA >= 1 .and. i+iA <= nx-1)       &
                                      .and.                     &
                           (j+jA >= 1 .and. j+jA <= ny-1) ) then

                         ! zero out A(row,col) and A(col,row)
                         ! Note: col = NodeID(k+kA, i+iA, j+jA) if we were forming a matrix with one column per node

                         ! Note: If we were setting (u,v) to something other than (0,0),
                         !  we would have to move terms to the RHS instead of just zeroing them.

!                            Auu( kA,  iA,  jA, row) = 0.d0
!                            Auu(-kA, -iA, -jA, col) = 0.d0
!                            Avv( kA,  iA,  jA, row) = 0.d0
!                            Avv(-kA, -iA, -jA, col) = 0.d0
!                            Auv( kA,  iA,  jA, row) = 0.d0
!                            Auv(-kA, -iA, -jA, col) = 0.d0
!                            Avu( kA,  iA,  jA, row) = 0.d0
!                            Avu(-kA, -iA, -jA, col) = 0.d0

                         if (iA==0 .and. jA==0 .and. kA==0) then
                            !! uncomment to put 1 on the main diagonal
                            m = indxA(0,0,0)
                            Auu(m,k,i,j) = 1.d0
                            Auv(m,k,i,j) = 0.d0
                            Avu(m,k,i,j) = 0.d0
                            Avv(m,k,i,j) = 1.d0

                            !else leave the diagonal term unchanged
                            !WHL: For the dome problem, it seems to make no difference whether we put a '1'
                            !     on the diagonal or leave the diagonal term unchanged. Answers are BFB
                            !     for diagonal preconditioner, SIA preconditioner, and GMRES/ILU

                         else

                            ! zero out the associated term in this row
                            m = indxA(iA,jA,kA)
                            Auu(m, k, i, j) = 0.d0
                            Auv(m, k, i, j) = 0.d0
                            Avu(m, k, i, j) = 0.d0
                            Avv(m, k, i, j) = 0.d0

                            ! zero out the associated term in this column
                            m = indxA(-iA,-jA,-kA)
                            Auu(m, k+kA, i+iA, j+jA) = 0.d0
                            Auv(m, k+kA, i+iA, j+jA) = 0.d0
                            Avu(m, k+kA, i+iA, j+jA) = 0.d0
                            Avv(m, k+kA, i+iA, j+jA) = 0.d0

                         endif   ! not on the diagonal

                      endif     ! i+iA, j+jA, and k+kA in bounds

                   enddo        ! kA
                   enddo        ! iA
                   enddo        ! jA

                   ! force u = v = 0 for this node by zeroing the rhs  
                   bu(k,i,j) = 0.d0            
                   bv(k,i,j) = 0.d0

                endif    ! umask_dirichlet
             enddo       ! k
          endif          ! active_vertex
       enddo             ! i
    enddo                ! j

  end subroutine dirichlet_boundary_conditions

!****************************************************************************

  subroutine compute_residual_vector(nx,          ny,            &
                                     nz,          nhalo,         &
                                     active_vertex,              &
                                     Auu,         Auv,           &
                                     Avu,         Avv,           &
                                     bu,          bv,            &
                                     uvel,        vvel,          &
                                     resid_vec_u, resid_vec_v,   &
                                     L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in structured (x/y/z) format.

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 3 (node and its nearest neighbors in x, y and z direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv              ! assembled load (rhs) vector, divided into 2 parts

   real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel          ! u and v components of velocity (m/yr)

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::   &
       resid_vec_u,      & ! residual vector, divided into 2 parts
       resid_vec_v

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector

    integer :: i, j, k, iA, jA, kA, m 

    if (verbose_residual .and. this_rank==rtest) then
       print*, ' '
       i = itest
       j = jtest
       k = ktest
       print*, 'i, j, k, u answer, u rhs:', i, j, k, uvel(k,i,j), bu(k,i,j)
       print*, 'i, j, k, v answer, v rhs:', i, j, k, vvel(k,i,j), bv(k,i,j)
    endif

    ! Compute u and v components of A*x

    resid_vec_u(:,:,:) = 0.d0
    resid_vec_v(:,:,:) = 0.d0

!    do j = 1, ny-1
!    do i = 1, nx-1

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo
       !TODO - Use indirect addressing to avoid an 'if' here?
       if (active_vertex(i,j)) then

          do k = 1, nz

             do kA = -1,1
             do jA = -1,1
             do iA = -1,1

                if ( (k+kA >= 1 .and. k+kA <= nz)     &
                                .and.                     &
                     (i+iA >= 1 .and. i+iA <= nx-1)         &
                             .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   m = indxA(iA,jA,kA)

                   resid_vec_u(k,i,j) = resid_vec_u(k,i,j)                        & 
                                      + Auu(m,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                      + Auv(m,k,i,j)*vvel(k+kA,i+iA,j+jA)

                   resid_vec_v(k,i,j) = resid_vec_v(k,i,j)                        &
                                      + Avu(m,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                      + Avv(m,k,i,j)*vvel(k+kA,i+iA,j+jA)

                endif   ! in bounds

             enddo   ! kA
             enddo   ! iA
             enddo   ! jA

          enddo   ! k

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j


    if (verbose_residual .and. this_rank==rtest) then
!       print*, ' '
!       print*, 'Compute residual:'
!       print*, 'k, i, j, Axu, bu:'
    endif

    ! Subtract b to get A*x - b
    ! Sum up squared L2 norm as we go

    L2_norm = 0.d0

    ! Loop over locally owned vertices

    do j = nhalo+1, ny-nhalo
    do i = nhalo+1, nx-nhalo
       if (active_vertex(i,j)) then
          do k = 1, nz
             resid_vec_u(k,i,j) = resid_vec_u(k,i,j) - bu(k,i,j)
             resid_vec_v(k,i,j) = resid_vec_v(k,i,j) - bv(k,i,j)
             L2_norm = L2_norm + resid_vec_u(k,i,j)*resid_vec_u(k,i,j)  &
                               + resid_vec_v(k,i,j)*resid_vec_v(k,i,j)
          enddo  ! k
       endif     ! active vertex
    enddo        ! i
    enddo        ! j

    ! Take global sum, then take square root

    L2_norm = parallel_reduce_sum(L2_norm)
    L2_norm = sqrt(L2_norm)

  end subroutine compute_residual_vector

!****************************************************************************

  subroutine compute_residual_velocity(nhalo,  whichresid, &
                                       uvel,   vvel,        &
                                       usav,   vsav,        &
                                       resid_velo)

    integer, intent(in) ::   &
       nhalo,           & ! number of layers of halo cells
       whichresid         ! option for method to use when calculating residual

!TODO - Specify dimensions?
    real(dp), dimension(:,:,:), intent(in) ::  &
       uvel, vvel,      & ! current guess for velocity
       usav, vsav         ! previous guess for velocity

    real(dp), intent(out) ::    &
       resid_velo         ! quantity related to velocity convergence


    integer ::   &
       imaxdiff, jmaxdiff, kmaxdiff   ! location of maximum speed difference
                                      ! currently computed but not used
 
    integer :: i, j, k, count

    real(dp) ::   &
       speed,      &   ! current guess for ice speed
       oldspeed,   &   ! previous guess for ice speed
       diffspeed       ! abs(speed-oldspeed)


  ! Compute a residual quantity based on convergence of the velocity field.
  !TODO - Are all of these methods needed?

  ! options for residual calculation method, as specified in configuration file
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel )
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels
  ! case(2): use mean of abs( vel_old - vel ) / vel )
  ! case(3): use max of abs( vel_old - vel ) / vel ) (in addition to L2 norm)

    resid_velo = 0.d0
    imaxdiff = 0
    jmaxdiff = 0
    kmaxdiff = 0

    select case (whichresid)

    case(HO_RESID_MAXU_NO_UBAS)   ! max speed difference, excluding the bed

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,3)-nhalo
          do i = 1+nhalo, size(uvel,2)-nhalo
             do k = 1, size(uvel,1) - 1         ! ignore bed velocity
                speed = sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2)
                if (speed /= 0.d0) then
                   oldspeed = sqrt(usav(k,i,j)**2 + vsav(k,i,j)**2)
                   diffspeed = abs((oldspeed - speed)/speed)
                   if (diffspeed > resid_velo) then
                      resid_velo = diffspeed
                      imaxdiff = i
                      jmaxdiff = j
                      kmaxdiff = k
                   endif
                endif
             enddo
          enddo
       enddo
       
       ! take global max
       resid_velo = parallel_reduce_max(resid_velo)

    case(HO_RESID_MEANU)   ! mean relative speed difference

       count = 0

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,3)-nhalo
          do i = 1+nhalo, size(uvel,2)-nhalo
             do k = 1, size(uvel,1) - 1         ! ignore bed velocity
                speed = sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2)
                if (speed /= 0.d0) then
                   count = count+1
                   oldspeed = sqrt(usav(k,i,j)**2 + vsav(k,i,j)**2)
                   diffspeed = abs((oldspeed - speed)/speed)                
                   resid_velo = resid_velo + diffspeed
                endif
             enddo
          enddo
       enddo

       if (count > 0) resid_velo = resid_velo / count

!TODO - How to convert this to a global value?
!       Maybe remove this case
       call not_parallel(__FILE__, __LINE__)


   case default    ! max speed difference, including basal speeds
                   ! (case HO_RESID_MAXU or HO_RESID_L2NORM)

       ! Loop over locally owned vertices

       do j = 1+nhalo, size(uvel,3)-nhalo
          do i = 1+nhalo, size(uvel,2)-nhalo
             do k = 1, size(uvel,1)
                speed = sqrt(uvel(k,i,j)**2 + vvel(k,i,j)**2)
                if (speed /= 0.d0) then
                   oldspeed = sqrt(usav(k,i,j)**2 + vsav(k,i,j)**2)
                   diffspeed = abs((oldspeed - speed)/speed)
                   if (diffspeed > resid_velo) then
                      resid_velo = diffspeed
                      imaxdiff = i
                      jmaxdiff = j
                      kmaxdiff = k
                   endif
                endif
             enddo
          enddo
       enddo

       resid_velo = parallel_reduce_max(resid_velo)
       
  end select

  end subroutine compute_residual_velocity

!---------------------------------------------------------------------------

  subroutine unstable_manifold_correction(nx,        ny,        &
                                          nz,        nhalo,     &
                                          uvel,      vvel,      &
                                          uvel_old,  vvel_old,  &
                                          ucorr_old, vcorr_old)

    ! Correct the velocity to improve convergence of the Picard solution

    ! input/output arguments

    integer, intent(in) ::     &
         nx, ny, nz,           &! grid dimensions
         nhalo                  ! number of halo layers 

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) :: &
         uvel, vvel,           &! velocity: on input, from linear solution
                                !           on output, corrected based on UMC
         uvel_old, vvel_old,   &! old velocity solution
         ucorr_old, vcorr_old   ! old correction vector

    ! local variables

    real(dp), dimension(nz,nx-1,ny-1) ::   &
         ucorr, vcorr                 ! correction vectors, uvel-uvel_old and vvel-vvel_old

    real(dp) ::     &
         cnorm, cnorm_old,           &! norm of correction vectors
         cdot,                       &! dot product of correction vectors
         cdiff, cdiff_u, cdiff_v,    &! difference of correction vectors
         theta,                      &! angle theta between succesive correction vectors (radians)
         cos_theta,                  &! cosine of angle theta
         alpha                        ! factor multiplying correction vector for underrelaxation

    integer :: i, j, k

    ! local parameters

    logical, parameter :: umc_split = .true.

    real(dp), parameter ::    &  
       theta_umc_underrelax = 5.d0*pi/6.d0   ! threshold angle for underrelaxation
                                             ! value of 5*pi/6 suggested by Hindmarsh and Payne
    real(dp), parameter ::    &
       small_velo = 1.d-16 * vel0            ! as in GLAM code, but adjusted for scaling


    ! compute preliminary correction vector based on linear solution

    ucorr(:,:,:) = uvel(:,:,:) - uvel_old(:,:,:)
    vcorr(:,:,:) = vvel(:,:,:) - vvel_old(:,:,:)

    if (umc_split) then     ! compute u and v corrections independently (as in glam)

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             do k = 1, nz

                ! uvel
                !!    tmp_vel = uvel(k,i,j)
                cnorm = abs(ucorr(k,i,j))
                cnorm_old = abs(ucorr_old(k,i,j))
                cdot = ucorr(k,i,j)*ucorr_old(k,i,j)
                cdiff = ucorr(k,i,j) - ucorr_old(k,i,j)

                if (cnorm*cnorm_old > 0.d0) then
                   !!  theta = acos(cdot/(cnorm*cnorm_old + small_velo))
                   !WHL - debug = Make sure cos_theta is in range before taking acos.
                   cos_theta = cdot / (cnorm*cnorm_old + small_velo)
                   if (cos_theta < -1.d0) then
                      print*, ' '
                      print*, 'COSINE OUT OF RANGE: i, j, k =', i, j, k
                      print*, 'cdot, cnorm, cnorm_old, cos_theta:', cdot, cnorm, cnorm_old, cos_theta
                      print*, 'Setting cos_theta = -1'
                      cos_theta = -1.d0
                   elseif (cos_theta > 1.d0) then
                      print*, ' '
                      print*, 'COSINE OUT OF RANGE: i, j, k =', i, j, k
                      print*, 'cdot, cnorm, cnorm_old, cos_theta:', cdot, cnorm, cnorm_old, cos_theta
                      print*, 'Setting cos_theta = -1'
                      cos_theta = 1.d0
                   endif
                   theta = acos(cos_theta)
                else
                   theta = 0.d0
                endif

                if (theta > theta_umc_underrelax .and. cdiff /= 0.d0) then
                   alpha = cnorm / cdiff
                   uvel(k,i,j) = uvel_old(k,i,j) + alpha*ucorr(k,i,j)

                   if (verbose_umc .and. this_rank==rtest) then
                      print*, ' '
                      print*, 'Underrelax: i, j, k, theta (deg), alpha:', &
                           i, j, k, theta*180.d0/pi, alpha
!                      print*, 'cnorm, cnorm_old, cdot:', cnorm, cnorm_old, cdot
                   endif

                endif
                !!  uvel_old(k,i,j) = tmp_vel

                ! vvel
                !!  tmp_vel = vvel(k,i,j)
                cnorm = abs(vcorr(k,i,j))
                cnorm_old = abs(vcorr_old(k,i,j))
                cdot = vcorr(k,i,j)*vcorr_old(k,i,j)
                cdiff = vcorr(k,i,j) - vcorr_old(k,i,j)

                if (cnorm*cnorm_old > 0.d0) then
                   !!  theta = acos(cdot/(cnorm*cnorm_old + small_velo))
                   !WHL - debug = Make sure cos_theta is in range before taking acos.
                   cos_theta = cdot / (cnorm*cnorm_old + small_velo)
                   if (cos_theta < -1.d0) then
                      print*, ' '
                      print*, 'COSINE OUT OF RANGE: i, j, k =', i, j, k
                      print*, 'cdot, cnorm, cnorm_old, cos_theta:', cdot, cnorm, cnorm_old, cos_theta
                      print*, 'Setting cos_theta = -1'
                      cos_theta = -1.d0
                   elseif (cos_theta > 1.d0) then
                      print*, ' '
                      print*, 'COSINE OUT OF RANGE: i, j, k =', i, j, k
                      print*, 'cdot, cnorm, cnorm_old, cos_theta:', cdot, cnorm, cnorm_old, cos_theta
                      print*, 'Setting cos_theta = -1'
                      cos_theta = 1.d0
                   endif
                   theta = acos(cos_theta)
                else
                   theta = 0.d0
                endif

                if (theta > theta_umc_underrelax .and. cdiff /= 0.d0) then
                   alpha = cnorm / cdiff
                   vvel(k,i,j) = vvel_old(k,i,j) + alpha*vcorr(k,i,j)

                   if (verbose_umc .and. this_rank==rtest) then
                      print*, ' '
                      print*, 'Underrelax: i, j, k, theta (deg), alpha:', &
                           i, j, k, theta*180.d0/pi, alpha
!                      print*, 'cnorm, cnorm_old, cdot:', cnorm, cnorm_old, cdot
                   endif

                endif

                !!  vvel_old(k,i,j) = tmp_vel

             enddo  ! k
          enddo     ! i
       enddo        ! j

    else    ! not split between u and v

       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             do k = 1, nz

                cnorm = sqrt(ucorr(k,i,j)*ucorr(k,i,j)   &
                           + vcorr(k,i,j)*vcorr(k,i,j))

                cnorm_old = sqrt(ucorr_old(k,i,j)*ucorr_old(k,i,j)  &
                               + vcorr_old(k,i,j)*vcorr_old(k,i,j))

                cdot = ucorr(k,i,j)*ucorr_old(k,i,j) + vcorr(k,i,j)*vcorr_old(k,i,j)

                cdiff_u = ucorr(k,i,j) - ucorr_old(k,i,j)
                cdiff_v = vcorr(k,i,j) - vcorr_old(k,i,j)
                cdiff = sqrt(cdiff_u*cdiff_u + cdiff_v*cdiff_v)

                ! compute angle theta between these two correction vectors
                !
                !              |   (c_n, c_(n-1))    |
                ! theta= arccos| ------------------- |
                !              |  |c_n| * |c_(n-1)|  |

                if (cnorm*cnorm_old > 0.d0) then
                   cos_theta = cdot / (cnorm*cnorm_old + small_velo)
                   if (cos_theta < -1.d0) then
                      print*, ' '
                      print*, 'COSINE OUT OF RANGE: i, j, k =', i, j, k
                      print*, 'cdot, cnorm, cnorm_old, cos_theta:', cdot, cnorm, cnorm_old, cos_theta
                      print*, 'Setting cos_theta = -1'
                      cos_theta = -1.d0
                   elseif (cos_theta > 1.d0) then
                      print*, ' '
                      print*, 'COSINE OUT OF RANGE: i, j, k =', i, j, k
                      print*, 'cdot, cnorm, cnorm_old, cos_theta:', cdot, cnorm, cnorm_old, cos_theta
                      print*, 'Setting cos_theta = -1'
                      cos_theta = 1.d0
                   endif
                   theta = acos(cos_theta)
                else
                   theta = 0.d0
                endif

                if (theta > theta_umc_underrelax .and. cdiff /= 0.d0) then

                   !!  alpha = cnorm/cdiff  !TODO - Compare to line below
                   alpha = min(cnorm/cdiff,1.d0)

                   uvel(k,i,j) = uvel_old(k,i,j) + alpha*ucorr(k,i,j)
                   vvel(k,i,j) = vvel_old(k,i,j) + alpha*ucorr(k,i,j)

                   if (verbose_umc .and. this_rank==rtest) then
                      print*, ' '
                      print*, 'Underrelax: i, j, k, theta (deg), alpha:', &
                           i, j, k, theta*180.d0/pi, alpha
                      print*, 'cnorm, cnorm_old, cdot:', cnorm, cnorm_old, cdot
                   endif

                endif   ! theta > theta_umc_underrelax

             enddo    ! k
          enddo       ! i
       enddo          ! j

    endif   ! umc_split

    ! Copy new to old vectors

    uvel_old(:,:,:) = uvel(:,:,:)
    vvel_old(:,:,:) = vvel(:,:,:)

    ucorr_old(:,:,:) = ucorr(:,:,:)
    vcorr_old(:,:,:) = vcorr(:,:,:)

  end subroutine unstable_manifold_correction


!****************************************************************************
! The next three subroutines are used for the SLAP solver only.
!****************************************************************************

  subroutine slap_preprocess(nx,           ny,          &
                             nz,           nNodesSolve, &
                             NodeID,                    &
                             iNodeIndex,   jNodeIndex,  &
                             kNodeIndex,                &
                             Auu,          Auv,         &
                             Avu,          Avv,         &  
                             bu,           bv,          &
                             uvel,         vvel,        &
                             matrix_order, nNonzero,    &
                             matrix,       rhs,         &
                             answer)

    !----------------------------------------------------------------
    ! Using the intermediate matrices (Auu, Auv, Avu, Avv), load vectors (bu, bv),
    ! and velocity components (uvel, vvel), form the matrix and the rhs and answer
    ! vectors in the desired sparse matrix format.
    !
    ! The matrix is formed in ascending row order, so it can easily be transformed
    ! to compressed sparse row (CSR) format without further sorting.
    !
    ! Note: This works only for single-processor runs with the SLAP solver.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels at which velocity is computed
       nNodesSolve              ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       NodeID             ! ID for each node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (k,i,j) indices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts
                          
    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel         ! u and v components of velocity

    integer, intent(in) ::    &
       matrix_order,  &   ! order of matrix = number of rows
       nNonzero           ! upper bound for number of nonzero entries in sparse matrix

    type(sparse_matrix_type), intent(inout) ::  &    ! TODO: inout or out?
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), intent(out) ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, k, iA, jA, kA, m, mm, n, ct

    integer :: rowA, colA   ! row and column of A submatrices (order = nNodesSolve)
    integer :: row, col     ! row and column of sparse matrix (order = 2*nNodesSolve) 

    real(dp) :: val         ! value of matrix coefficient
    
    ! Set the nonzero coefficients of the sparse matrix 

    ct = 0

    do rowA = 1, nNodesSolve

       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)

       ! Load the nonzero values associated with Auu and Auv
       ! These are assigned a value of matrix%row = 2*rowA - 1

!TODO: Make sure that this procedure captures all the nonzero entries

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
             m = indxA(iA,jA,kA)

             ! Auu
             val = Auu(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA - 1                    
                matrix%val(ct) = val
             endif

             ! Auv 
             val = Auv(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

             if (verbose_slapsolve .and. 2*colA==rowtest .and. 2*rowA-1==coltest) then
               print*, ' '
               print*, 'rowA, colA, i, j, k =', rowA, colA, i, j, k
               print*, 'iA, jA, kA:', iA, jA, kA
               print*, 'Auv(iA,jA,kA,i,j,k) =', Auv(m,k,i,j)
               mm = indxA(-iA,-jA,-kA)
               print*, 'Avu(-iA,-jA,-kA,i+iA,j+jA,k+kA) =', Avu(mm, k+kA, i+iA, j+jA)
               print*, ' '
               print*, 'ct, row, col, val:', ct, matrix%row(ct), matrix%col(ct), matrix%val(ct)
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

       ! Load the nonzero values associated with Avu and Avv
       ! These are assigned a value of matrix%row = 2*rowA

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
             m  = indxA( iA, jA, kA)
             mm = indxA(-iA,-jA,-kA)

             ! Avu 
             val = Avu(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA - 1
                matrix%val(ct) = val
             endif

!WHL - debug
             if (verbose_slapsolve .and. 2*rowA==rowtest .and. 2*colA-1==coltest) then
               print*, ' '
               print*, 'rowA, colA, i, j, k =', rowA, colA, i, j, k
               print*, 'iA, jA, kA:', iA, jA, kA
               print*, 'Avu(iA,jA,kA,i,j,k) =', Avu(m,k,i,j)
               print*, 'Auv(-iA,-jA,-kA,i+iA,j+jA,k+kA) =', Auv(mm, k+kA, i+iA, j+jA)
               print*, ' '
               print*, 'ct, row, col, val:', ct, matrix%row(ct), matrix%col(ct), matrix%val(ct)
             endif

             ! Avv 
             val = Avv(m,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

    enddo           ! rowA 

    ! Set basic matrix parameters.

    matrix%order = matrix_order
    matrix%nonzeros = ct
    matrix%symmetric = .false.

    if (verbose_slapsolve) then
       print*, ' '
       print*, 'solver preprocess'
       print*, 'order, nonzeros =', matrix%order, matrix%nonzeros
    endif

    ! Initialize the answer vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       answer(2*n-1) = uvel(k,i,j)
       answer(2*n)   = vvel(k,i,j)

       if (verbose_slapsolve .and. n==ntest) then
          print*, ' '
          print*, 'n, initial uvel =', n, answer(2*n-1)
          print*, 'n, initial vvel =', n, answer(2*n)
       endif

    enddo

    ! Set the rhs vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodesSolve
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       rhs(2*n-1) = bu(k,i,j)
       rhs(2*n)   = bv(k,i,j)

       if (verbose_slapsolve .and. n==ntest) then
          print*, ' '
          print*, 'n, initial bu =', n, rhs(2*n-1)
          print*, 'n, initial bv =', n, rhs(2*n)
       endif

    enddo

  end subroutine slap_preprocess

!****************************************************************************

  subroutine slap_postprocess(nNodesSolve,  answer,                   &
                              iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                              uvel,         vvel)

  ! Extract the velocities from the solution vector.
  ! Note: This works only for single-processor runs with the SLAP solver.
                                            
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nNodesSolve     ! number of nodes where we solve for velocity

    real(dp), dimension(:), intent(in) ::  &
       answer             ! velocity solution vector

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of active nodes

    real(dp), dimension(:,:,:), intent(inout) ::   &
       uvel, vvel         ! u and v components of velocity

    integer :: i, j, k, n

       if (verbose_slapsolve) then
          print*, ' '
          print*, 'solver postprocess: n, i, j, k, u, v'
       endif

    do n = 1, nNodesSolve

       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       uvel(k,i,j) = answer(2*n-1)
       vvel(k,i,j) = answer(2*n)

    enddo

  end subroutine slap_postprocess

!****************************************************************************

  subroutine slap_compute_residual_vector(matrix,    answer,   rhs, &
                                          resid_vec, L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in triad (row/col/val) format.

    type(sparse_matrix_type), intent(in) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    real(dp), dimension(:), intent(in) ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    real(dp), dimension(:), intent(out) ::   &
       resid_vec          ! residual vector

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector

    integer :: i, j, n

    if (verbose_residual) then
       print*, ' '
       print*, 'Residual vector: order, nonzeros =', matrix%order, matrix%nonzeros 
       print*, 'ntest, u answer, u rhs:', ntest, answer(2*ntest-1), rhs(2*ntest-1)
       print*, 'ntest, v answer, v rhs:', ntest, answer(2*ntest),   rhs(2*ntest)
    endif

    resid_vec(:) = 0.d0

    do n = 1, matrix%nonzeros
       i = matrix%row(n)
       j = matrix%col(n)
       resid_vec(i) = resid_vec(i) + matrix%val(n)*answer(j)
    enddo

    L2_norm = 0.d0
    do i = 1, matrix%order
       resid_vec(i) = resid_vec(i) - rhs(i)
       L2_norm = L2_norm + resid_vec(i)*resid_vec(i)
    enddo

!PARALLEL - Parallel sum needed here?  Or will we never call this subroutine in parallel?
!    L2_norm = parallel_reduce_sum(L2_norm)

    L2_norm = sqrt(L2_norm)

  end subroutine slap_compute_residual_vector

!****************************************************************************
! The remaining subroutines are used for testing and bug-checking but may 
! not need to be called for production runs.
!****************************************************************************

  subroutine check_symmetry_element_matrix(nNodesPerElement,  &
                                           Kuu, Kuv, Kvu, Kvv)

    !------------------------------------------------------------------
    ! Check that the element stiffness matrix is symmetric.
    ! This is true provided that (1) Kuu = (Kuu)^T
    !                            (2) Kvv = (Kvv)^T
    !                            (3) Kuv = (Kvu)^T
    ! This check should not be needed for production runs with a well-tested code,
    !  but is included for now to help with debugging.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement  ! number of nodes per element

    real(dp), dimension(nNodesPerElement, nNodesPerElement), intent(in) ::   &
             Kuu, Kuv, Kvu, Kvv     ! component of element stiffness matrix
                                    !
                                    !    Kuu  | Kuv
                                    !    _____|____          
                                    !    Kvu  | Kvv
                                    !         |

    integer :: i, j

    real(dp) :: avg_val


    ! make sure Kuu = (Kuu)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kuu(i,j) - Kuu(j,i)) > eps10) then
             print*, 'Kuu is not symmetric'
             print*, 'i, j, Kuu(i,j), Kuu(j,i):', i, j, Kuu(i,j), Kuu(j,i)
             stop
          endif    
       enddo
    enddo

    ! check that Kvv = (Kvv)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kvv(i,j) - Kvv(j,i)) > eps10) then
             print*, 'Kvv is not symmetric'
             print*, 'i, j, Kvv(i,j), Kvv(j,i):', i, j, Kvv(i,j), Kvv(j,i)
             stop
          endif    
       enddo
    enddo

    ! Check that Kuv = (Kvu)^T

    do j = 1, nNodesPerElement
       do i = 1, nNodesPerElement
          if (abs(Kuv(i,j) - Kvu(j,i)) > eps10) then
             print*, 'Kuv .ne. (Kvu)^T'
             print*, 'i, j, Kuv(i,j), Kvu(j,i):', i, j, Kuv(i,j), Kvu(j,i)
             stop
          endif    
       enddo
    enddo

  end subroutine check_symmetry_element_matrix

!****************************************************************************

  subroutine check_symmetry_assembled_matrix(nx,  ny,  nz, nhalo,   &
                                             active_vertex,         &
                                             Auu, Auv, Avu, Avv)

    !------------------------------------------------------------------
    ! Check that the assembled stiffness matrix is symmetric.
    ! This is true provided that (1) Auu = (Auu)^T
    !                            (2) Avv = (Avv)^T
    !                            (3) Auv = (Avu)^T
    ! The A matrices are assembled in a dense fashion to save storage
    !  and preserve the i/j/k structure of the grid..
    !
    ! There can be small differences from perfect symmetry due to roundoff error.
    ! These differences are fixed provided they are small enough.
    !------------------------------------------------------------------    

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex            ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv, Avu, Avv       ! components of assembled stiffness matrix
                                !
                                !    Auu  | Auv
                                !    _____|____          
                                !         |
                                !    Avu  | Avv                                    

    integer :: i, j, k, iA, jA, kA, m, mm

    real(dp) :: val1, val2          ! values of matrix coefficients

    real(dp) :: maxdiff, diag_entry, avg_val

    ! Check matrix for symmetry

    ! Here we correct for small differences from symmetry due to roundoff error.
    ! The maximum departure from symmetry is set to be a small fraction (1.e-12) 
    !  of the diagonal entry for the row.
    ! If the departure from symmetry is larger than this, then the model prints a warning 
    !  and/or aborts.

    maxdiff = 0.d0

! Loop over locally owned vertices.
! Each active vertex is associate with 2*nz matrix rows belonging to this processor.
! Locally owned vertices are (nhalo+1:ny-nhalo, nhalo+1:nx-nhalo)
!           
    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz

                ! Check Auu and Auv for symmetry

                m = indxA(0,0,0)
                diag_entry = Auu(m,k,i,j)

                do jA = -1, 1
                do iA = -1, 1
                do kA = -1, 1

                   if (k+kA >= 1 .and. k+kA <=nz) then  ! to keep k index in bounds

                      m =  indxA( iA, jA, kA)
                      mm = indxA(-iA,-jA,-kA)

                      ! Check that Auu = Auu^T

                      val1 = Auu( m, k,    i,    j   )   ! value of Auu(row,col)
                      val2 = Auu(mm, k+kA, i+iA, j+jA)   ! value of Auu(col,row)

                      if (val2 /= val1) then
                          
                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps10*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Auu( m, k,   i,   j   ) = avg_val
                            Auu(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Auu, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'WARNING: Auu is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Auu(row,col), Auu(col,row), diff:', val1, val2, val2 - val1
!!                            stop  !TODO - Put in a proper abort
                         endif

                      endif   ! val2 /= val1
                
                      ! Check that Auv = (Avu)^T

                      val1 = Auv( m, k,    i,    j)      ! value of Auv(row,col)
                      val2 = Avu(mm, k+kA, i+iA, j+jA)   ! value of Avu(col,row)

                      if (val2 /= val1) then

                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps10*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Auv( m, k,   i,   j   ) = avg_val
                            Avu(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Auv, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'WARNING: Auv is not equal to (Avu)^T, i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Auv(row,col), Avu(col,row), diff:', val1, val2, val2 - val1
!!                            stop  !TODO - Put in a proper abort
                         endif

                      endif  ! val2 /= val1

                   endif     ! k+kA in bounds
            
                enddo        ! kA
                enddo        ! iA
                enddo        ! jA

                ! Now check Avu and Avv

                m = indxA(0,0,0)
                diag_entry = Avv(m,k,i,j)

                ! check that Avv = (Avv)^T

                do jA = -1, 1
                do iA = -1, 1
                do kA = -1, 1

                   if (k+kA >= 1 .and. k+kA <=nz) then  ! to keep k index in bounds

                      m  = indxA( iA, jA, kA)
                      mm = indxA(-iA,-jA,-kA)

                      val1 = Avv( m, k,    i,    j)      ! value of Avv(row,col)
                      val2 = Avv(mm, k+kA, i+iA, j+jA)   ! value of Avv(col,row)

                      if (val2 /= val1) then
                          
                         if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps10*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Avv( m, k,   i,   j   ) = avg_val
                            Avv(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Avv, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'WARNING: Avv is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Avv(row,col), Avv(col,row), diff:', val1, val2, val2 - val1
!!                            stop  !TODO - Put in a proper abort
                         endif

                      endif   ! val2 /= val1

                      ! Check that Avu = (Auv)^T

                      val1 = Avu( m, k,    i,    j)      ! value of Avu(row,col)
                      val2 = Auv(mm, k+kA, i+iA, j+jA)   ! value of Auv(col,row)

                      if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

                      if (val2 /= val1) then

                         ! if difference is small, then fix the asymmetry by averaging values
                         ! else print a warning and abort

                         if ( abs(val2-val1) < eps10*abs(diag_entry) ) then
                            avg_val = 0.5d0 * (val1 + val2)
                            Avu( m, k,   i,   j   ) = avg_val
                            Auv(mm, k+kA,i+iA,j+jA) = avg_val
!!                            print*, 'Avu, i, j, k, iA, jA, kA, val1, val2:', i, j, k, iA, jA, kA, val1, val2
                         else
                            print*, ' '
                            print*, 'Avu is not equal to (Auv)^T, i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                            print*, 'Avu(row,col), Auv(col,row), diff:', val1, val2, val2 - val1
                            stop  !TODO - Put in a proper abort
                         endif

                      endif  ! val2 /= val1

                   endif     ! k+kA in bounds

                enddo        ! kA
                enddo        ! iA
                enddo        ! jA

             enddo     ! k
          endif        ! active_vertex
       enddo           ! i
    enddo              ! j

    maxdiff = parallel_reduce_max(maxdiff)
    if (verbose_matrix .and. main_task) then
       print*, ' '
       print*, 'Max difference from symmetry =', maxdiff
    endif

  end subroutine check_symmetry_assembled_matrix

!****************************************************************************

  subroutine remove_small_values_assembled_matrix(nx,  ny,  nz, nhalo,   &
                                                  active_vertex,         &
                                                  Auu, Auv, Avu, Avv)

    !------------------------------------------------------------------
    ! Remove very small values from the assembled matrix.
    ! "Very small" means much smaller than the diagonal terms.
    !------------------------------------------------------------------    

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels where velocity is computed
       nhalo                    ! number of halo layers

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex            ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv, Avu, Avv       ! components of assembled stiffness matrix
                                !
                                !    Auu  | Auv
                                !    _____|____          
                                !         |
                                !    Avu  | Avv                                    

    integer :: i, j, k, iA, jA, kA, m, mm

    real(dp) :: val        ! value of matrix coefficients

    real(dp) :: diag_entry

    ! Look for small values and remove them, along with symmetric counterparts.

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo
          if (active_vertex(i,j)) then
             do k = 1, nz

                m = indxA(0,0,0)
                diag_entry = Auu(m,k,i,j)

                do jA = -1, 1
                do iA = -1, 1
                do kA = -1, 1

                   if (k+kA >= 1 .and. k+kA <= nz) then

                      m =  indxA( iA, jA, kA)
                      mm = indxA(-iA,-jA,-kA)

                      val = Auu( m, k,    i,    j   )   ! value of Auu(row,col)
                      if ( abs(val) < eps10*abs(diag_entry) ) then
                         Auu(m,k,i,j) = 0.d0             ! Auu(row,col)
                         Auu(mm,k+kA,i+iA,j+jA) = 0.d0   ! Auu(col,row)
                      endif

                      val = Auv( m, k,    i,    j   )   ! value of Auv(row,col)
                      if ( abs(val) < eps10*abs(diag_entry) ) then
                         Auv(m,k,i,j) = 0.d0             
                         Avu(mm,k+kA,i+iA,j+jA) = 0.d0
                      endif
                   endif

                enddo  ! kA
                enddo  ! iA
                enddo  ! jA

                m = indxA(0,0,0)
                diag_entry = Avv(m,k,i,j)

                do jA = -1, 1
                do iA = -1, 1
                do kA = -1, 1

                   if (k+kA >= 1 .and. k+kA <= nz) then
                      m =  indxA( iA, jA, kA)
                      mm = indxA(-iA,-jA,-kA)

                      val = Avv( m, k,    i,    j   )   ! value of Avv(row,col)
                      if ( abs(val) < eps10*abs(diag_entry) ) then
                         Avv(m,k,i,j) = 0.d0             ! Avv(row,col)
                         Avv(mm,k+kA,i+iA,j+jA) = 0.d0   ! Avv(col,row)
                      endif

                      val = Avu( m, k,    i,    j   )   ! value of Avu(row,col)
                      if ( abs(val) < eps10*abs(diag_entry) ) then
                         Avu(m,k,i,j) = 0.d0             ! Avu(row,col)
                         Auv(mm,k+kA,i+iA,j+jA) = 0.d0   ! Auv(col,row)
                      endif

                   endif

                enddo  ! kA
                enddo  ! iA
                enddo  ! jA

             enddo     ! k
          endif        ! active_vertex
       enddo           ! i
    enddo              ! j

    end subroutine remove_small_values_assembled_matrix

!****************************************************************************

  subroutine solve_test_matrix (matrix_order, whichsparse)

    ! solve a small test matrix

    integer, intent(in) :: &
       matrix_order,       & ! matrix order
       whichsparse           ! solution method (0=BiCG, 1=GMRES, 2=PCG_DIAG, 3=PCG_INCH, 5 = STANDALONE_PCG)

    logical :: verbose_test = .true.

    type(sparse_matrix_type) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types

    real(dp), dimension(:), allocatable ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer             ! answer (x) in Ax = b

    real(dp), dimension(:,:), allocatable :: Atest

    real(dp) :: err

    integer :: niters, nNonzero_max

    integer :: i, j, n

    print*, 'Solving test matrix, order =', matrix_order

    nNonzero_max = matrix_order*matrix_order    ! not sure how big this must be

    allocate(Atest(matrix_order,matrix_order))
    Atest(:,:) = 0.d0

    allocate(matrix%row(nNonzero_max), matrix%col(nNonzero_max), matrix%val(nNonzero_max))
    allocate(rhs(matrix_order), answer(matrix_order))

    rhs(:) = 0.d0
    answer(:) = 0.d0
    matrix%row(:) = 0
    matrix%col(:) = 0
    matrix%val(:) = 0.d0

    matrix%order = matrix_order
    matrix%symmetric = .false.

    if (matrix%order == 2) then    ! symmetric 2x2
       Atest(1,1:2) = (/3.d0, 2.d0 /)
       Atest(2,1:2) = (/2.d0, 6.d0 /)
       rhs(1:2) = (/2.d0, -8.d0 /)   ! answer = (2 -2) 

    elseif (matrix%order == 3) then

       ! symmetric
       Atest(1,1:3) = (/ 7.d0, -2.d0,  0.d0 /)
       Atest(2,1:3) = (/-2.d0,  6.d0, -2.d0 /)
       Atest(3,1:3) = (/ 0.d0, -2.d0,  5.d0 /)
       rhs(1:3)   =   (/ 3.d0,  8.d0,  1.d0 /)   ! answer = (1 2 1)

        ! non-symmetric
!       Atest(1,1:3) = (/3.d0,   1.d0,  1.d0 /)
!       Atest(2,1:3) = (/2.d0,   2.d0,  5.d0 /)
!       Atest(3,1:3) = (/1.d0,  -3.d0, -4.d0 /)
!       rhs(1:3)   =  (/ 6.d0,  11.d0, -9.d0 /)   ! answer = (1 2 1)

    else if (matrix%order == 4) then

        ! symmetric

       Atest(1,1:4) = (/ 2.d0, -1.d0,  0.d0,  0.d0 /)
       Atest(2,1:4) = (/-1.d0,  2.d0, -1.d0,  0.d0 /)
       Atest(3,1:4) = (/ 0.d0, -1.d0,  2.d0, -1.d0 /)
       Atest(4,1:4) = (/ 0.d0,  0.d0, -1.d0,  2.d0 /)
       rhs(1:4)    = (/  0.d0,  1.d0, -1.d0,  4.d0 /)   ! answer = (1 2 2 3)

        ! non-symmetric
!       Atest(1,1:4) = (/3.d0,  0.d0,  2.d0, -1.d0 /)
!       Atest(2,1:4) = (/1.d0,  2.d0,  0.d0,  2.d0 /)
!       Atest(3,1:4) = (/4.d0,  0.d0,  6.d0, -3.d0 /)
!       Atest(4,1:4) = (/5.d0,  0.d0,  2.d0,  0.d0 /)
!       rhs(1:4)    = (/ 6.d0,  7.d0, 13.d0,  9.d0 /)   ! answer = (1 2 2 1)

    elseif (matrix%order > 4) then
  
        Atest(:,:) = 0.d0
        do n = 1, matrix%order 
           Atest(n,n) = 2.d0
           if (n > 1) Atest(n,n-1) = -1.d0
           if (n < matrix%order) Atest(n,n+1) = -1.d0
        enddo

        rhs(1) = 1.d0
        rhs(matrix%order) = 1.d0
        rhs(2:matrix%order-1) = 0.d0              ! answer = (1 1 1 ... 1 1 1)
        
    endif

    if (verbose_test) then
       print*, ' '
       print*, 'Atest =', Atest
       print*, 'rhs =', rhs
    endif

    ! Put in SLAP triad format (column ascending order)

    n = 0
    do j = 1, matrix%order
       do i = 1, matrix%order
          if (Atest(i,j) /= 0.d0) then 
             n = n + 1
             matrix%row(n) = i
             matrix%col(n) = j
             matrix%val(n) = Atest(i,j)
          endif
       enddo
    enddo

    ! Set number of nonzero values
    matrix%nonzeros = n

    if (verbose_test) then
       print*, ' '
       print*, 'row,       col,       val:'
       do n = 1, matrix%nonzeros
          print*, matrix%row(n), matrix%col(n), matrix%val(n)
       enddo
       print*, 'Call sparse_easy_solve, whichsparse =', whichsparse
    endif

    ! Solve the linear matrix problem


    call sparse_easy_solve(matrix, rhs,    answer,  &
                           err,    niters, whichsparse)

    if (verbose_test) then
       print*, ' '
       print*, 'answer =', answer
       print*, 'err =', err       
       print*, 'niters =', niters
    endif

    stop

  end subroutine solve_test_matrix

!****************************************************************************

!WHL - This subroutine is currently not called.
!      It might be useful if the matrix contains many entries that are
!       close to but not quite equal to zero (e.g., due to roundoff error).
!      It may need to be modified to preserve symmetry (so we don't zero out 
!       an entry without zeroing out its symmetric partner).

  subroutine remove_small_matrix_terms(nx, ny, nz,   &
                                       Auu, Auv,     &
                                       Avu, Avv) 


    !------------------------------------------------------------------------
    ! Get rid of matrix entries that are very small compared to diagonal entries.
    !   
    ! This handles a potential numerical issue: Some entries should be zero because of 
    !  cancellation of terms, but are not exactly zero because of rounding errors.
    !------------------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz                            ! number of vertical layers at which velocity is computed

    real(dp), dimension(27,nz,nx-1,ny-1), intent(inout) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer :: i, j, k, m, iA, jA, kA

    real(dp) ::       &
       diag_entry,   &! mean size of diagonal entries
       min_entry      ! minimum size of entries kept in matrix

    do j = 1, ny-1
    do i = 1, nx-1
    do k = 1, nz

       ! Compute threshold value of matrix entries for Auu/Auv row

       m = indxA(0,0,0)
       diag_entry = Auu(m,k,i,j)
       min_entry = eps10 * diag_entry

       ! Remove very small values

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          m = indxA(iA,jA,kA)
          if (abs(Auu(m,k,i,j)) < min_entry) then
!             if (abs(Auu(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Auu term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Auu(m,k,i,j)
!             endif
             Auu(m,k,i,j) = 0.d0
          endif

          if (abs(Auv(m,k,i,j)) < min_entry) then
!             if (abs(Auv(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Auv term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Auv(m,k,i,j)
!             endif
             Auv(m,k,i,j) = 0.d0
          endif

       enddo
       enddo
       enddo

       ! Compute threshold value of matrix entries for Avu/Avv row

       m = indxA(0,0,0)
       diag_entry = Avv(m,k,i,j)
       min_entry = eps10 * diag_entry

       ! Remove very small values

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          m = indxA(iA,jA,kA)

          if (abs(Avu(m,k,i,j)) < min_entry) then
!             if (abs(Avu(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Avu term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Avu(m,k,i,j)
!             endif
             Avu(m,k,i,j) = 0.d0
          endif

          if (abs(Avv(m,k,i,j)) < min_entry) then
!             if (abs(Avv(m,k,i,j)) > 0.d0) then
!                print*, 'Dropping Avv term:, i, j, k, iA, jA, kA, val:', i, j, k, iA, jA, kA, Avv(m,k,i,j)
!             endif
             Avv(m,k,i,j) = 0.d0
          endif

       enddo
       enddo
       enddo

  enddo
  enddo
  enddo
     
  end subroutine remove_small_matrix_terms

!****************************************************************************

  subroutine write_matrix_elements(nx,    ny,   nz,     &
                                   nNodesSolve, NodeID, &
                                   iNodeIndex,  jNodeIndex,  &
                                   kNodeIndex,          &
                                   Auu,         Auv,    &
                                   Avu,         Avv,    &
                                   bu,          bv)

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels at which velocity is computed
       nNodesSolve              ! number of nodes where we solve for velocity

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       NodeID             ! ID for each node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(27,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its nearest neighbors in x, y and z direction 
                          ! other dimensions = (k,i,j) indices

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts


    ! Local variables

    integer :: rowA, colA
    integer :: i, j, k, m, iA, jA, kA

    real(dp), dimension(nNodesSolve, nNodesSolve) ::   &
       Auu_val, Auv_val, Avu_val, Avv_val   ! dense matrices

    real(dp), dimension(nNodesSolve) :: nonzeros

    Auu_val(:,:) = 0.d0
    Auv_val(:,:) = 0.d0
    Avu_val(:,:) = 0.d0
    Avv_val(:,:) = 0.d0

    do rowA = 1, nNodesSolve

       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)

       do kA = -1, 1
       do jA = -1, 1
       do iA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)         &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx-1)       &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
             m = indxA(iA,jA,kA)

             if (colA > 0) then 
                Auu_val(rowA, colA) = Auu(m,k,i,j)
                Auv_val(rowA, colA) = Auv(m,k,i,j)
                Avu_val(rowA, colA) = Avu(m,k,i,j)
                Avv_val(rowA, colA) = Avv(m,k,i,j)
             else
!                print*, 'row, i, j, k, iA, jA, kA, col:', rowA, i, j, k, iA, jA, kA, colA
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

    enddo           ! rowA 

    !WHL - bug check
    print*, ' '
    print*, 'nonzeros per row:'
    do rowA = 1, nNodesSolve
       nonzeros(rowA) = 0
       do colA = 1, nNodesSolve
          if (abs(Auu_val(rowA,colA)) > 1.d-11) then
             nonzeros(rowA) = nonzeros(rowA) + 1
          endif
       enddo
!       print*, rowA, nonzeros(rowA)
    enddo

    print*, 'Write matrix elements to file, label =', matrix_label

    ! Write matrices to file (one line of file corresponding to each row of matrix)

    open(unit=10, file='Auu.'//matrix_label, status='unknown')
    open(unit=11, file='Auv.'//matrix_label, status='unknown')
    open(unit=12, file='Avu.'//matrix_label, status='unknown')
    open(unit=13, file='Avv.'//matrix_label, status='unknown')

    do rowA = 1, nNodesSolve
!       print*, 'row =', rowA
       write(10,'(i6)',advance='no') rowA
       write(11,'(i6)',advance='no') rowA
       write(12,'(i6)',advance='no') rowA
       write(13,'(i6)',advance='no') rowA
       do colA = 1, nNodesSolve
          write(10,'(e16.8)',advance='no') Auu_val(rowA,colA)
          write(11,'(e16.8)',advance='no') Auv_val(rowA,colA)
          write(12,'(e16.8)',advance='no') Avu_val(rowA,colA)
          write(13,'(e16.8)',advance='no') Avv_val(rowA,colA)
       enddo
       write(10,*) ' '
       write(11,*) ' '
       write(12,*) ' '
       write(13,*) ' '
    enddo

    close(10)
    close(11)
    close(12)
    close(13)

    print*, 'Done writing matrix elements'

    ! write load vectors to file
    open(unit=14, file='bu.'//matrix_label, status='unknown')
    open(unit=15, file='bv.'//matrix_label, status='unknown')
    do rowA = 1, nNodesSolve
       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)
       write(14,'(i6, e16.8)') rowA, bu(k,i,j)
       write(15,'(i6, e16.8)') rowA, bv(k,i,j)
    enddo
    close(14)
    close(15)

  end subroutine write_matrix_elements
  
!****************************************************************************

!WHL TODO - Remove this subroutine when done with test case
  subroutine ishomC_block_init(model,  &
                               dx,   dy,  &
                               thck, usrf, topg, beta)

!WHL - Hack for removing periodic BC and isolating the ISHOM C geometry to a 4x4 block
!      at the center of the domain

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    real(dp), intent(in) :: dx, dy   ! gridcell dimensions

    real(dp), dimension(:,:), intent(inout) ::  &
       thck, usrf, topg, beta

    real(dp), parameter :: alpha = 0.1d0 * pi/180.d0   ! slope angle
    real(dp) :: L   ! domain size
    real(dp) :: omega  ! frequency of beta variation
    real(dp) :: xdiff, ydiff

    integer :: ilo, imid, ihi, jlo, jmid, jhi
    integer :: i, j

    imid = model%general%ewn / 2
    jmid = model%general%nsn / 2

    ilo = imid-1
    ihi = imid+2

    jlo = jmid-1
    jhi = jmid+2

    print*, 'ilo, ihi =', ilo, ihi
    print*, 'jlo, jhi =', jlo, jhi
    print*, 'tan(alpha) =', tan(alpha)

    ! Reset thck, usrf, topg
    thck(:,:) = 0.d0
    usrf(:,:) = 0.d0
    do j = jlo, jhi
       do i = ilo, ihi
          thck(i,j) = 1000.d0
          xdiff = dx * (i-ilo+0.5d0)
          usrf(i,j) = 1000.d0 - xdiff * tan(alpha)
       enddo
    enddo
    topg(:,:) = usrf(:,:) - thck(:,:)


    ! Reset beta

    beta(:,:) = 0.d0
!    L = model%general%ewn * dx ! for full domain
    L = 4 * dx    ! for 4x4 block

    omega = 2.d0*pi / L

    print*, 'L, omega =', L, omega
    do j = jlo-1, jhi
       do i = ilo-1, ihi
          xdiff = dx * (i-ilo+1)
          ydiff = dy * (j-jlo+1)
          beta(i,j) = 1000.d0 + 1000.d0 * sin(omega*xdiff) * sin(omega*ydiff)
       enddo
    enddo

  end subroutine ishomC_block_init

!****************************************************************************

  end module glissade_velo_higher

!****************************************************************************
