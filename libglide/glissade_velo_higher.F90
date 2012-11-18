! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! +                                                           + 
! +  glissade_velo_higher.F90                                 + 
! +                                                           + 
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

    use glimmer_global, only: sp, dp
    use glimmer_physcon, only: gn, rhoi, rhoo, grav, scyr
    use glimmer_paramets, only: thk0, len0, tim0, tau0   !TODO - remove thk0 and len0
    use glimmer_log, only: write_log
    use glimmer_sparse_type
    use glimmer_sparse

    use cism_sparse_pcg, only: pcg_solver_structured, &
                               cism_sparse_pcg_solve
    
    use parallel, only: main_task

!    use glide_deriv
!    use glide_grids, only: stagvarb

    implicit none

    private
    public :: glissade_velo_higher_init, glissade_velo_higher_solve

!TODO - Declare nhalo in parallel modules (in addition to or instead of lhalo and uhalo).    
    public :: nhalo, itest, jtest, ktest

    integer, parameter :: nhalo = 2

    !----------------------------------------------------------------
    ! Here are some definitions:
    !
    ! The horizontal mesh is composed of cells and vertices.
    ! All cells are assumed to be quadrilaterals, but the code can be
    !  generalized later to triangles (e.g., for MPAS mesh).
    ! Each cell can be extruded to form a column with a specified number of layers.
    ! 
    ! An element is a layer of a cell, and a node is a corner of an element.
    ! So elements and nodes live in 3D space, whereas cells and vertices live in
    !  the horizontal plane.
    !
    ! Active elements are those elements lying in cells with ice present
    !  (that is, with an ice thickness above a minimum threshold).
    ! Active nodes are the nodes of these active elements.
    !
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! basic cell properties (quadrilateral cells)
    !----------------------------------------------------------------

    integer :: nCells     ! total number of cells on this processor
    integer :: nVertices  ! total number of vertices associated with at least one cell

    integer, parameter ::      &
       nVerticesOnCell = 4,    & ! number of vertices bordering each cell   
       nVerticesOnVertex = 9     !TODO - not needed?

    integer, dimension(:,:), allocatable ::  &
       CellID                ! unique ID for each cell (i,j)

    integer, dimension(:,:), allocatable ::  &
       VertexID              ! unique ID for each vertex (i,j)

    real(dp), dimension(:,:), allocatable :: &
       xVertex, yVertex       ! x and y coordinates of each vertex

    !----------------------------------------------------------------
    ! Finite element properties
    !
    ! For now, assume 3D hexahedral elements.
    !----------------------------------------------------------------

    integer, parameter ::      &
       nNodesPerElement = 8,   &   ! 8 nodes for hexahedral elements
       nQuadPoints = 8             ! number of quadrature points per element
                                   ! These live at +- 1/sqrt(3) for reference hexahedron

    real(dp), parameter ::     &
       rsqrt3 = 1.d0/sqrt(3.d0)    ! for quadrature points
         
    !----------------------------------------------------------------
    ! Arrays used for finite-element calculations
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement, nQuadPoints) ::   & 
       phi,            &    ! basis function, evaluated at quad pts
       dphi_dxr,       &    ! dphi/dx for reference element, evaluated at quad pts
       dphi_dyr,       &    ! dphi/dy for reference element, evaluated at quad pts
       dphi_dzr             ! dphi/dy for reference element, evaluated at quad pts

    real(dp), dimension(nQuadPoints) ::  &
       xqp, yqp, zqp,  &    ! quad pt coordinates in reference element
       wqp                  ! quad pt weights

    integer, dimension(nNodesPerElement, nNodesPerElement) ::  &
       ishift, jshift, kshift   ! matrices describing relative indices of nodes in an element

    real(dp), dimension(3,3) ::  &
       identity3                ! 3 x 3 identity matrix

    real(dp), parameter :: eps12 = 1.d-12   ! small number

!WHL - TODO - Keep volume scale?
    real(dp) :: vol0    ! volume scale = dx * dy * (1000 m)

!WHL - debug
    logical :: verbose_init = .false.   ! for debug print statements
    logical :: verbose_Jac = .false.
    logical :: verbose = .true.  
!WHL - debug

    integer, parameter :: &
       itest = 24, jtest = 17, ktest = 1, ptest = 1
!       itest = 24, jtest = 17, ktest = 2, ptest = 1
!       itest = 24, jtest = 17, ktest = 10, ptest = 1
!       itest = 12, jtest = 5, ktest = 1, ptest = 1
!       itest = 12, jtest = 26, ktest = 1, ptest = 1

    integer, parameter :: ntest = 2371  ! nodeID for (24,17,1)    
!    integer, parameter :: ntest = 2372  ! nodeID for (24,17,2)
!    integer, parameter :: ntest = 2380  ! nodeID for (24,17,10)

!    integer, parameter :: ntest = 1     ! nodeID for (12,5,1)
!    integer, parameter :: ntest = 3841     ! nodeID for (12,26,1)

    contains

!****************************************************************************

  subroutine glissade_velo_higher_init (nx, ny,  &
                                        nhalo,   &
                                        dx, dy)

    !----------------------------------------------------------------
    ! Assign global indices to each cell and vertex.
    ! (Not sure where these will be needed.)
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &     ! number of grid cells in each direction
       nhalo                       ! number of halo layers

    real(dp), intent(in) ::  &
       dx,  dy                     ! grid cell dimensions

    integer :: i, j, n, p, nc, nv

!WHL - debug
    real(dp) :: sumx, sumy, sumz

    !----------------------------------------------------------------
    ! Assign indices for the cells and vertices
    ! TODO - Decide whether these indices should be local or global for parallel code,
    !        and rewrite as needed.
    !----------------------------------------------------------------

    ! Define global cell index for cells on this processor.
    ! Note: Cell 1 is at lower left, cell nCells is at upper right

!WHL - MPAS has indexToCellID and indexToVertexID.  

    ! Set volume scale
    ! This is not absolutely necessary, but dividing by this scale gives matrix coefficients 
    !  that are closer to unity.

    vol0 = dx * dy * 1000.d0    ! typical volume of ice column (m^3)

    if (verbose_init) then
       print*, ' '
       print*, 'In glissade_velo_higher_init'
       print*, 'vol0 =', vol0
    endif


!TODO - Deallocate CellID and VertexID at end of run? 
    ! Assign a unique index to each cell on this processor.

    allocate(CellID(nx,ny))
    nCells = 0
    do j = 1, ny   ! OK to count halo cells?  Or should we loop over local cells only?
    do i = 1, nx
       nCells = nCells + 1
       CellID(i,j) = nCells
    enddo
    enddo

    ! Assign a unique index to each vertex and compute its x and y coordinates.
    ! By convention, vertex (i,j) lies at the NE corner of cell(i,j).
    !TODO - These coordinates are relative to local processor; modify for parallel case?

    allocate(VertexID(nx-1, ny-1))
    allocate(xVertex(nx-1,ny-1), yVertex(nx-1,ny-1))

    nVertices = 0
    VertexID(:,:) = 0
    xVertex(:,:) = 0.d0
    yVertex(:,:) = 0.d0
    do j = 1, ny-1  ! Restrict loop to vertices of local cells?
    do i = 1, nx-1
       nVertices = nVertices + 1
       VertexID(i,j) = nVertices
       xVertex(i,j) = dx * i
       yVertex(i,j) = dy * j
    enddo
    enddo

    if (verbose_init) then

       ! Cells
       write(6,*) ' '
       write(6,*) 'Cell IDs:'
       do j = 1, ny
          write(6,*) ' '
          do i = 1, nx
             write(6,*) 'i, j, cellID:', i, j, CellID(i,j)
          enddo
       enddo
       
       ! Vertices
       write(6,*) ' '
       write(6,*) 'Vertex IDs:'
       do j = 1, ny-1
          write(6,*) ' '    
          do i = 1, nx-1
             write(6,*) 'i, j, vertexID, xV, yV:', i, j, VertexID(i,j), xVertex(i,j), yVertex(i,j)
          enddo
       enddo

    endif

    !----------------------------------------------------------------
    ! Initialize some time-independent finite element arrays
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
   
    ! Set coordinates and weights of quadrature points for reference hexahedral element.
    ! Numbering is counter-clockwise from southwest, lower face (1-4) followed by
    !  upper face (5-8).

    xqp(1) = -rsqrt3; yqp(1) = -rsqrt3; zqp(1) = -rsqrt3
    wqp(1) =  1.d0   !TODO - check that weight = 1 and not 1/8

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
       sumx = 0.d0; sumy = 0.d0; sumz = 0.d0
       do p = 1, nQuadPoints
          print*, xqp(p), yqp(p), zqp(p)
          sumx = sumx + xqp(p); sumy = sumy + yqp(p); sumz = sumz + zqp(p)
       enddo
       print*, ' '
       print*, 'sums:', sumx, sumy, sumz
    endif

    ! Evaluate basis functions and their derivatives at each quad pt

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
          do n = 1, 8
             print*, n, phi(n,p)
!             print*, n, dphi_dxr(n,p)
!             print*, n, dphi_dyr(n,p)
!             print*, n, dphi_dzr(n,p)
          enddo
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

  end subroutine glissade_velo_higher_init

!****************************************************************************

  subroutine glissade_velo_higher_solve(nx,       ny,           &
                                        nz,       sigma,        &
                                        nhalo,                  &
                                        thck,     usrf,         &
!!                                        stagthck, stagusrf,   &
                                        thklim,                 &
                                        flwa,                   &
                                        uvel,     vvel,         &
!                                        beta,                  &
!                                        whichbabc,             &
                                        whichefvs,              &
                                        whichresid,             &
                                        whichnonlinear,         &
                                        whichsparse,            &
                                        whichapprox)

!TODO - Something like these may be needed if building with Trilinos.
!!    use glam_strs2, only: linearSolveTime, totalLinearSolveTime

!TODO - Is this 'use' better than hardwiring the numbers?
!!    use glide_types, only: HO_APPROX_SIA, HO_APPROX_SSA 

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

!TODO - Compute nx, ny, and nlyr based on size(thck) and size(flwa)?

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz,                   &  ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)
                                ! Note: The number of layers for scalars is nz-1
       nhalo                    ! number of rows/columns of halo cells
 
    real(dp), dimension(:), intent(in) :: &
       sigma

    real (dp), dimension(:,:), intent(in) ::  &
       thck,                 &  ! ice thickness
       usrf                     ! upper surface elevation

!TODO - May want to compute locally instead of passing in
!    real (dp), dimension(:,:), intent(in) ::  &
!       stagusrf,            & ! upper surface averaged to vertices
!       stagthck               ! ice thickness averaged to vertices

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells

    real (dp), dimension(:,:,:), intent(in) ::  &
       flwa                   ! flow factor

    real (dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel             ! velocity components

    integer, intent(in) :: whichefvs      ! option for effective viscosity calculation 
                                          ! (calculate it or make it uniform)
    integer, intent(in) :: whichresid     ! option for method to use when calculating residual
    integer, intent(in) :: whichnonlinear ! options for which nonlinear method (Picard or JFNK)
    integer, intent(in) :: whichsparse    ! options for which method for doing elliptic solve
                                          ! (BiCG, GMRES, standalone Trilinos, etc.)
    integer, intent(in), optional :: whichapprox    ! options for which Stokes approximation to use
                                                    ! 0 = SIA, 1 = SSA, 2 = Blatter-Pattyn HO
                                                    ! default = 2

    !--------------------------------------------------------
    ! Local parameters
    !--------------------------------------------------------

    integer, parameter :: cmax = 300          ! max number of outer iterations

    !WHL - What I call resid_target is called minresid in glam_strs2
    !TODO - I think these values are tuned for old scaling parameters; adjust them?
    real(dp), parameter :: resid_target = 1.0d-04   ! assume velocity fields have converged below this resid 
    real(dp), parameter :: NL_tol   = 1.0d-06   ! to have same criterion as JFNK

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer ::            &
       nActiveCells,      &   ! no. of active cells (thck > thklim)
       nElements,         &   ! no. of elements (nz-1 per active cell)
       nNodes                 ! no. of nodes belonging to these elements

    real(dp), dimension(nx-1,ny-1) ::   &
       stagusrf,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    real(dp), dimension(nx,ny) ::    &
       rmask                  ! = 1. where ice is present, else = 0.

    logical, dimension(nx,ny) ::     &
       active_cell            ! true for active cells (thck > thklim)

    logical, dimension(nx-1,ny-1) :: &
       active_vertex          ! true for vertices of active cells

    integer, dimension(nz,nx-1,ny-1) ::  &
       NodeID             ! ID for each active node

    integer, dimension(nVertices*nz) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nVertices*nz) :: &
       xNode, ynode, znode   ! x, y and z coordinates for each node

    integer, dimension(nNodesPerElement,nz-1,nx,ny) ::  &
       NodeOnElementID     ! node ID for each node of each element

!    real(dp), dimension(nz-1,nx,ny) :: &
!       visc               ! effective viscosity of each element

    real(dp), dimension(nz-1,nx,ny) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute effective viscosity

    real(dp), dimension(:,:,:), allocatable ::   &
       usav, vsav         ! previous guess for velocity solution

    logical, dimension(:,:,:), allocatable ::    &
       umask_dirichlet    ! Dirichlet mask for velocity (if true, u = v = 0)

    real(dp), dimension(:,:,:,:,:,:), allocatable ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 3 (node and its 2 neighbors in z direction) 
                          ! 2nd dimension = 3 (node and its 2 neighbors in x direction) 
                          ! 3rd dimension = 3 (node and its 2 neighbors in y direction) 
                          ! 4th dimension = nNodes

    real(dp), dimension(:,:,:), allocatable ::  &
       bu, bv,          & ! assembled load vector, divided into 2 parts
       resid_vec_u,    & ! residual vector, divided into 2 parts
       resid_vec_v

    type(sparse_matrix_type) ::  &
       matrix             ! sparse matrix, defined in glimmer_sparse_types
                          ! includes nonzeroes, order, col, row, val 

    integer ::    &
       matrix_order,    & ! order of matrix = number of rows
       nNonzero           ! upper bound for number of nonzero entries in sparse matrix

    real(dp), dimension(:), allocatable ::   &
       rhs,             & ! right-hand-side (b) in Ax = b
       answer,          & ! answer (x) in Ax = b
       resid_vec          ! residual vector Ax - b

    real(dp) :: &
       resid_velo,          & ! quantity related to velocity convergence
       L2_norm,             & ! L2 norm of residual
       !WHL- What I call L2_target is called NL_target in glam_strs2
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

!WHL - debug
    integer :: i, j, k, n
    integer :: iA, jA, kA, colA
    real(dp) :: ds_dx, ds_dy
    logical, parameter :: test_matrix = .false.
    integer, parameter :: test_order = 20
    integer :: rowi
    logical, parameter :: sia_test = .false.

!TODO - Pass in whichapprox?

    if (present(whichapprox)) then
!!       if (whichapprox == HO_APPROX_SIA) then   ! SIA
       if (whichapprox == 0) then   ! SIA
          if (verbose) print*, 'Solving shallow-ice approximation'
          sia_factor = 1.d0
          ssa_factor = 0.d0
!!       elseif (whichapprox == HO_APPROX_SSA) then  ! SSA
       elseif (whichapprox == 1) then  ! SSA
          if (verbose) print*, 'Solving shallow-shelf approximation'
          sia_factor = 0.d0
          ssa_factor = 1.d0
       else   ! Blatter-Pattyn higher-order 
          if (verbose) print*, 'Solving Blatter-Pattyn higher-order approximation'
          sia_factor = 1.d0
          ssa_factor = 1.d0
       endif
    else   ! Blatter-Pattyn higher-order by default
       if (verbose) print*, 'Solving Blatter-Pattyn higher-order approximation'
       sia_factor = 1.d0
       ssa_factor = 1.d0
    endif

    if (test_matrix) then
       call solve_test_matrix(test_order, whichsparse)
    endif

    !------------------------------------------------------------------------------
    ! Setup for higher-order solver: Compute nodal geometry, allocate storage, etc.
    ! These are quantities that do not change during the outer nonlinear loop. 
    !------------------------------------------------------------------------------

    if (verbose) then
       print*, ' '
       print*, 'In glissade_velo_higher_solve'
       print*, 'nx, ny, nz:', nx, ny, nz
       print*, 'vol0:', vol0
       print*, 'sigma:', sigma
       print*, 'thklim:', thklim
       print*, 'max thck:', maxval(thck(:,:))
       print*, 'max usrf:', maxval(usrf(:,:))
       print*, ' '
    endif

!WHL - debug
    if (verbose) then
       print*, ' '
       print*, 'Thickness field:'
       do j = ny, 1, -1
          write(6,'(34f6.0)') thck(:,j)
       enddo
    endif
 
    if (verbose) then
        print*, ' '
        print*, 'i, j, usrf(i,j)'
        i = itest+1
        do j = ny, 1, -1
           print*, i, j, usrf(i,j)
        enddo
    endif

    ! Compute real mask: = 1 where ice is present, 0 elsewhere

    do j = 1, ny
       do i = 1, nx
          if (thck(i,j) > thklim) then
	     rmask(i,j) = 1.d0
          else
             rmask(i,j) = 0.d0
          endif
       enddo
    enddo

    ! Compute ice thickness on staggered grid

    call staggered_scalar(nx,           ny,         &
                          nhalo,        rmask,      &
                          thck,         stagthck)

    ! Compute upper surface elevation on staggered grid

    call staggered_scalar(nx,           ny,         &
                          nhalo,        rmask,      &
                          usrf,         stagusrf)

    ! Identify and count the active cells (i.e., cells with thck > thklim)
    ! Assign a unique ID to each node, and identify the nodes of each element.

    call get_nodal_geometry(nx,          ny,         nz, &   
                            nhalo,       sigma,    &
                            thck,        thklim,  &
                            stagusrf,    stagthck, &
                            nElements,   nNodes,          & 
                            active_cell, active_vertex,   &
                            NodeID,      NodeOnElementID,             &
                            iNodeIndex,  jNodeIndex,  kNodeIndex,   &
                            xNode,       yNode,       zNode)

    ! Compute flow factor appearing in the expression for effective viscosity.

    do j = 1, ny
       do i = 1, nx
          if (thck(i,j) > thklim) then
             ! gn = exponent in Glen's flow law (= 3 by default)
             flwafact(:,i,j) = 0.5d0 * flwa(:,i,j)**(-1.d0/real(gn,dp))  
          else
             flwafact(:,i,j) = 0.d0
          endif
       enddo
    enddo

    if (verbose) then
       n = ntest
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)
       print*, ' '
       print*, 'ntest, i, j, k:', n, i, j, k
    endif
 
    if (verbose) then
       print*, ' '
       i = itest; j = jtest; k = ktest
       print*, 'itest, jtest, ktest:', i, j, k
       print*, 'NodeID =', NodeID(k,i,j)
       print*, 'NodeOnElementID =', NodeOnElementID(:,k,i,j)
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
       print*, 'flwa, flwafact:', flwa(k,i,j), flwafact(k,i,j)
    endif

    ! Allocate various arrays
    ! Note: nx and ny are the horizontal dimensions of scalar arrays (e.g., thck and temp)
    !       The velocity arrays have dimensions (nx-1, ny-1)
    ! (TODO: Some are not needed if using standalone Trilinos solver?)

    allocate(umask_dirichlet(nz,nx-1,ny-1))
    allocate(usav(nz,nx-1,ny-1), vsav(nz,nx-1,ny-1))  
    allocate(Auu(-1:1,-1:1,-1:1,nz,nx-1,ny-1), Auv(-1:1,-1:1,-1:1,nz,nx-1,ny-1))
    allocate(Avu(-1:1,-1:1,-1:1,nz,nx-1,ny-1), Avv(-1:1,-1:1,-1:1,nz,nx-1,ny-1))
    allocate(bu(nz,nx-1,ny-1), bv(nz,nx-1,ny-1))

    if (whichsparse == STANDALONE_PCG_STRUC) then  ! allocate space for the residual vector

       allocate(resid_vec_u(nz,nx-1,ny-1), resid_vec_v(nz,nx-1,ny-1))

    else  ! allocate space for the sparse matrix (A), rhs (b), answer (x), and residual vector (Ax-b).
       
       matrix_order = 2*nNodes         ! Is this exactly enough?
       nNonzero = matrix_order*54      ! 27 = node plus 26 nearest neighbors in hexahedral lattice   
                                       ! 54 = 2 * 27 (since solving for both u and v)

       allocate(matrix%row(nNonzero), matrix%col(nNonzero), matrix%val(nNonzero))
       allocate(rhs(matrix_order), answer(matrix_order), resid_vec(matrix_order))

       if (verbose) then
          print*, 'matrix_order =', matrix_order
          print*, 'nNonzero = ', nNonzero
       endif

    endif   ! STANDALONE_PCG_STRUC or not
 

    print*, 'size(uvel, vvel) =', size(uvel), size(vvel)
    print*, 'sum(uvel, vvel) =', sum(uvel), sum(vvel)
    print*, 'size(usav, vsav) =', size(usav), size(vsav)

    !---------------------------------------------------------------
    ! Solve Ax = b
    !---------------------------------------------------------------

    if (main_task) then
       ! print some info to the screen to update on iteration progress
       print *, ' '
       print *, 'Running Glissade higher-order dynamics solver'
       print *, ' '
       if (whichresid == 3) then                 ! use L2 norm of residual
          print *, 'iter #     resid (L2 norm)       target resid'
       else                                     ! use max value of residual
          print *, 'iter #     velo resid            target resid'
       end if
       print *, ' '
    endif

    ! Any ghost preprocessing needed, as in glam_strs2?

    ! set initial values 

    counter = 0
    resid_velo = 1.d0
    L2_norm = 1.0d20
    L2_target = 1.0d-4      !TODO - Is this the best value?  
                            ! Should it stay fixed during the outer loop, or should it evolve?

    outer_it_criterion = 1.0d10   ! guarantees at least one loop
    outer_it_target = 1.0d-12 

    do while (outer_it_criterion >= outer_it_target .and. counter < cmax)

       ! Advance the iteration counter

       counter = counter + 1

       if (verbose) then
          print*, ' '
          print*, 'Outer counter =', counter
          print*, 'whichresid =', whichresid
          print*, 'L2_norm, L2_target =', L2_norm, L2_target
          print*, 'resid_velo, resid_target =', resid_velo, resid_target
       endif

       ! save current velocity
       usav(:,:,:) = uvel(:,:,:)
       vsav(:,:,:) = vvel(:,:,:)

       ! Assemble the stiffness matrix A and load vector b in the linear system Ax = b

       if (verbose) print*, 'call assemble_linear_system'

       !TODO - Remove unneeded arguments?

       call assemble_linear_system(nx, ny, nz,  nhalo,              &
                                   active_cell,                       &
                                   nCells,           nActiveCells,    &
                                   NodeID,           NodeOnElementID, &
                                   iNodeIndex,       jNodeIndex,      &
                                   kNodeIndex,                        &
                                   xNode,            yNode,           &
                                   zNode,                             &
                                   nElements,        nNodes,          &
                                   uvel,             vvel,            &
                                   stagusrf,         flwafact,        &
                                   whichefvs,                         &
                                   Auu,              Auv,             &
                                   Avu,              Avv,             &
                                   bu,               bv,              &
                                   sia_factor,       ssa_factor)
       
!WHL - debug
       if (verbose) then

          i = itest
          j = jtest
          k = ktest
          n = ntest
          print*, ' '
          print*, 'Auu after assembly: i, j, k, n =', i, j, k, n
          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1
             print*, 'iA, jA, kA, Auu:', iA, jA, kA, Auu(kA, iA, jA, k, i, j)
          enddo
          enddo
          enddo

          i = itest
          j = jtest
          k = ktest+1
          n = ntest+1
          print*, ' '
          print*, 'Auu after assembly: i, j, k, n =', i, j, k, n
          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1
             print*, 'iA, jA, kA, Auu:', iA, jA, kA, Auu(kA, iA, jA, k, i, j)
          enddo
          enddo
          enddo

       endif

       ! Incorporate Dirichlet (u = 0) basal boundary conditions

       ! For now, simply impose zero sliding everywhere at the bed.
       ! TODO: Move above the start of the nonlinear loop

       umask_dirichlet(nz,:,:) = .true.   ! u = v = 0 at bed

       if (verbose) print*, 'Call dirichlet_bc'

       call dirichlet_boundary_conditions(nx,  ny,         nz,                &
                                          active_vertex,   umask_dirichlet,   &
                                          nNodes,          NodeID,            &
                                          Auu,             Auv,               &
                                          Avu,             Avv,               &
                                          bu,              bv)

       !TODO - Incorporate a basal sliding BC

!WHL - Try stripping out columns from the matrix and see if it can still solve an SIA problem.

          if (sia_test) then

             print*, ' '
             print*, 'SIA test: Removing Auv, Avu, and iA/jA columns'

             ! Set Auv = Avu = 0

             Auv(:,:,:,:,:,:) = 0.d0 
             Avu(:,:,:,:,:,:) = 0.d0 

             ! Remove terms from iA and jA columns

             do j = 1, ny-1
             do i = 1, nx-1
             do k = 1, nz
                do jA = -1,1
                do iA = -1,1
                do kA = -1,1
                   if (iA /= 0 .or. jA /= 0) then
                      Auu(kA,iA,jA,k,i,j) = 0.d0
                      Avv(kA,iA,jA,k,i,j) = 0.d0
                   endif
                enddo
                enddo
                enddo
             enddo
             enddo
             enddo

          endif   ! sia_test             



    ! Check symmetry of assembled matrix
    ! This call could be skipped in a well-tested production code.

       if (verbose) print*, 'Check matrix symmetry'

       call check_symmetry_assembled_matrix(nx,          ny,         nz,      &
                                            nNodes,      NodeID,              &     
                                            iNodeIndex,  jNodeIndex,          &
                                            kNodeIndex,                       &
                                            Auu,         Auv,                 &
                                            Avu,         Avv)

       if (verbose) then
          i = itest
          j = jtest
          k = ktest
          n = ntest

          print*, ' '
      	  print*, 'i,j,k,n =', i, j, k, n
          print*, 'Auu sum =', sum(Auu(:,:,:,k,i,j))
          print*, 'Auv sum =', sum(Auv(:,:,:,k,i,j))
          print*, 'Avu sum =', sum(Avu(:,:,:,k,i,j))
          print*, 'Avv sum =', sum(Avv(:,:,:,k,i,j))

          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1
             if ( (k+kA >= 1 .and. k+kA <= nz)       &
                             .and.                   &
                  (i+iA >= 1 .and. i+iA <= nx-1)     &
                             .and.                   &
                  (j+jA >= 1 .and. j+jA <= ny-1) ) then
                colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node
                print*, ' '
                print*, 'iA, jA, kA, colA =', iA, jA, kA, colA
                print*, 'i, j, k for colA:', iNodeIndex(colA), jNodeIndex(colA), kNodeIndex(colA)
                print*, 'Auu =', Auu(kA,iA,jA,k,i,j)
                print*, 'Auv =', Auv(kA,iA,jA,k,i,j)
                print*, 'Avu =', Avu(kA,iA,jA,k,i,j)
                print*, 'Avv =', Avv(kA,iA,jA,k,i,j)
             endif
          enddo
          enddo
          enddo

       endif

       if (whichsparse == STANDALONE_PCG_STRUC) then   ! standalone PCG for structured grid

          ! preprocessing is not needed because the matrix, rhs, and solution already
          ! have the required data structure

          if (verbose) then
             print*, ' '
             print*, 'Calling structured PCG solver'
             print*, 'sum(uvel, vvel) =', sum(uvel), sum(vvel)
          endif

          call pcg_solver_structured(nx,        ny,            &
                                     nz,        active_vertex, &
                                     Auu,       Auv,           &
                                     Avu,       Avv,           &
                                     bu,        bv,            &
                                     uvel,      vvel,          &
                                     err,       niters,        &
                                     NodeID)   !WHL - Node ID is temporary for debugging

          ! Compute the residual vector and its L2 norm

          call compute_residual_vector_structured(nx,          ny,            &
                                                  nz,          active_vertex, &
                                                  Auu,         Auv,           &
                                                  Avu,         Avv,           &
                                                  bu,          bv,            &
                                                  uvel,        vvel,          &
                                                  resid_vec_u, resid_vec_v,   &
                                                  L2_norm)

       elseif (whichsparse /= STANDALONE_TRILINOS_SOLVER) then   
          
          ! Put the nonzero matrix elements into the SLAP triad format
          ! (For SLAP solver; may not be needed for Trilinos)

          if (verbose) print*, 'Form global matrix in sparse format'
 
          matrix%order = matrix_order
          matrix%nonzeros = nNonzero
          matrix%symmetric = .false.

          ! Given the intermediate stiffness matrices (Auu, etc.) and load vector (bu, bv),
          ! form the global matrix (in sparse matrix format) and rhs.

          call solver_preprocess(nx,           ny,           nz,           &
                                 nNodes,       NodeID,                     &
                                 iNodeIndex,   jNodeIndex,   kNodeIndex,   &
                                 Auu,          Auv,   &
                                 Avu,          Avv,   &
                                 bu,           bv,    &
                                 uvel,         vvel,  &
                                 matrix_order, nNonzero,  &
                                 matrix,       rhs,   &
                                 answer)

          !TODO - glam_strs2 calls res_vect here--why?
          !!!call res_vect( matrix, vk_1, rhs, size(rhs), g_flag, L2square, whichsparse )
          !! L2norm  = L2square
          !! F(1:pcgsize(1)) = vk_1(:)

!WHL - bug check
          if (verbose) then
             print*, ' '
             print*, 'Before linear solve: n, row, col, val:'
             print*, ' '
             do n = 1, 10              
                print*, n, matrix%row(n), matrix%col(n), matrix%val(n)
             enddo
          endif

          ! Solve the linear matrix problem

          if (whichsparse == STANDALONE_PCG_UNSTRUC) then

!TODO - Eliminate this solve at some point?  It doesn't add new functionality relative to SLAP.

             if (verbose) print*, 'Call standalone cism pcg solver, counter =', counter

             call cism_sparse_pcg_solve(matrix,  rhs,    answer,   &
                                        err,     niters, verbose)

          elseif (whichsparse /= STANDALONE_TRILINOS_SOLVER) then

             if (verbose) then
                print*, 'Call sparse_easy_solve, counter =', counter
             endif

             call sparse_easy_solve(matrix, rhs,    answer,  &
                                    err,    niters, whichsparse)

             ! Put the velocity solution back into 3D arrays

          endif   ! STANDALONE_PCG_UNSTRUC or not

          if (verbose) then
             print*, 'Solved the linear system, niters, err =', niters, err
             print*, ' '
             print*, 'n, u, v (m/yr):', ntest, answer(2*ntest-1), answer(2*ntest)
!!           do n = 1, matrix%order 
!!              print*, n, answer(n)
!!           enddo
          endif

          call solver_postprocess(nNodes,       answer,                   &
                                  iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                  uvel,         vvel)

          ! Compute the residual vector and its L2 norm

          call compute_residual_vector(matrix,    answer,   rhs,  &
                                       resid_vec, L2_norm)

#ifdef TRILINOS
       else

!WHL - Commented out for now until testing later with Trilinos
!!           call solvewithtrilinos(rhs, answer, linearSolveTime)
!!           totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
          !  write(*,*) 'Total linear solve time so far', totalLinearSolveTime                                           
#endif

       endif   ! whichsparse (STANDALONE_PCG_STRUC or not) 

       ! Any ghost postprocessing needed, as in glam_strs2?

       ! Halo updates for velocity (to be tested)
       !TODO - Are all four of these calls needed?
!       call staggered_parallel_halo(uvel)
!       call staggered_parallel_halo(vvel)
!       call horiz_bcs_stag_vector_ew(uvel)
!       call horiz_bcs_stag_vector_ns(vvel)

!WHL - bug check
       if (verbose) then
          print*, ' '
          print*, 'After postprocess: n, i, j, k, uvel, vvel'
          do n = ntest-10, ntest+10              
             i = iNodeIndex(n)
             j = jNodeIndex(n)
             k = kNodeIndex(n)
             print*, n, i, j, k, uvel(k,i,j), vvel(k,i,j)
          enddo
       endif

       ! Compute residual quantities based on the velocity solution

       call compute_residual_velocity(nhalo,  whichresid,   &
                                      uvel,   vvel,        &
                                      usav,   vsav,        &
                                      resid_velo)
  
       if (verbose) then
          print*, ' '
          print*, 'L2_resid =', L2_norm
          print*, 'resid_velo =', resid_velo
       endif

       ! Write diagnostics (iteration number, max residual, and location of max residual
       ! (send output to the screen or to the log file, per whichever line is commented out) 

       !TODO - Find out if this note from glam_strs2 still applies:
       ! "Can't use main_task flag because main_task is true for all processors in case of parallel_single"

       if (main_task) then
          if (whichresid == 3 )then
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

       ! update the outer loop stopping criterion
       if (whichresid == 3) then
          outer_it_criterion = L2_norm
          outer_it_target = L2_target      ! L2_target is currently set to 1.d-4 and held constant
       else
          outer_it_criterion = resid_velo
          outer_it_target = resid_target   ! resid_target is currently a parameter = 1.d-4  
       end if

    enddo  ! while (outer_it_criterion >= outer_it_target .and. counter < cmax)

    if (counter < cmax) converged_soln = .true.

    if (verbose .and. converged_soln) then
       print*, ' '
       print*, 'GLISSADE SOLUTION HAS CONVERGED, outer counter =', counter

       i = itest
       j = jtest
       k = ktest
       print*, 'i, j, k, uvel, vvel:', i, j, k, uvel(k,i,j), vvel(k,i,j)
    endif

    ! Clean up

    deallocate(umask_dirichlet)
    deallocate(usav, vsav)
    deallocate(Auu, Auv, Avu, Avv)
    deallocate(bu, bv)
    if (whichsparse == STANDALONE_PCG_STRUC) then
       deallocate(resid_vec_u, resid_vec_v)
    else
       deallocate(matrix%row, matrix%col, matrix%val)
       deallocate(rhs, answer, resid_vec)
    endif

  end subroutine glissade_velo_higher_solve

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

    do j = nhalo, ny-nhalo
       do i = nhalo, nx-nhalo
          sumvar = rmask(i,j+1)*var(i,j+1) + rmask(i+1,j+1)*var(i+1,j+1)  &
                 + rmask(i,j)  *var(i,j)   + rmask(i+1,j)  *var(i+1,j)	  
          summask = rmask(i,j+1) + rmask(i+1,j+1) + rmask(i,j) + rmask(i+1,j)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
    enddo  

  end subroutine staggered_scalar

!****************************************************************************

  subroutine get_nodal_geometry(nx,          ny,          nz,       &   
                                nhalo,       sigma,    &
                                thck,        thklim,   &
                                stagusrf,    stagthck, &
                                nElements,   nNodes,   &
                                active_cell, active_vertex,             &
                                NodeID,      NodeOnElementID,           &
                                iNodeIndex,  jNodeIndex,  kNodeIndex,   &
                                xNode,       yNode,       zNode)
                            
    ! Identify and count the active nodes for the calculations.
    ! All nodes of all elements lying within active cells (i.e., cells with ice present)
    !  are active.
    ! Active nodes may be either free (to be solved for) or constrained (e.g., Dirichlet).

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical levels where velocity is computed
       nhalo                    

    real(dp), dimension(nz), intent(in) :: &
       sigma                  ! sigma vertical coordinate

    real (dp), dimension(nx,ny), intent(in) ::  &
       thck                   ! ice thickness

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells

    real (dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    integer, intent(out)  ::            &
       nElements,         &   ! no. of elements (nz-1 per active cell)
       nNodes                 ! no. of nodes belonging to these elements

    logical, dimension(nx,ny), intent(out) :: &
       active_cell            ! true for active cells (thck > thklim)

    logical, dimension(nx-1,ny-1), intent(out) :: &
       active_vertex          ! true for vertices of active cells

    integer, dimension(nz,nx-1,ny-1), intent(out) ::  &
       NodeID             ! ID for each active node

    integer, dimension(nNodesPerElement,nz-1,nx,ny), intent(out) ::  &
       NodeOnElementID     ! node ID for each node of each element

    integer, dimension(nVertices*nz), intent(out) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nVertices*nz), intent(out) :: &
       xNode, yNode, zNode   ! x, y and z coordinates for each node


    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer ::   &
       nActiveCells,   &      ! number of active cells (thck > thklim)
       nActiveVertices        ! number of active vertices

    integer :: i, j, k

    if (verbose) then
       print*, 'In get_nodal_geometry, thklim =', thklim
       print*, 'Active cells:'
    endif

    nActiveCells = 0
    active_cell(:,:) = .false.
    do j = 1+nhalo, ny-nhalo
    do i = 1+nhalo, nx-nhalo
       if (thck(i,j) >= thklim) then
          active_cell(i,j) = .true.
          nActiveCells = nActiveCells + 1
          if (verbose) then
!!             print*, i, j, thck(i,j)
          endif
       endif
    enddo
    enddo

    ! Count the active elements (nz-1 elements per active cell)

    nElements = nActiveCells * (nz-1)

    if (verbose) then
       print*, 'nActiveCells =', nActiveCells
       print*, 'nElements =', nElements
!!       print*, 'Active vertices:'
    endif


    ! Identify the active vertices (i.e., all vertices of active cells)

    active_vertex(:,:) = .false.

    do j = 1+nhalo, ny-nhalo   ! loop over local cells
    do i = 1+nhalo, nx-nhalo
       if (active_cell(i,j)) then
          active_vertex(i-1:i, j-1:j) = .true.  ! all vertices of this cell
       endif
    enddo
    enddo

    ! Identify and count the active nodes, and compute their x, y and z coordinates.

    nNodes = 0
    iNodeIndex(:) = 0
    jNodeIndex(:) = 0
    kNodeIndex(:) = 0    
    xNode(:) = 0.d0
    yNode(:) = 0.d0
    zNode(:) = 0.d0

    nActiveVertices = 0

    do j = nhalo, ny-nhalo    ! include S edge
    do i = nhalo, nx-nhalo    ! include W edge
       if (active_vertex(i,j)) then
          nActiveVertices = nActiveVertices + 1
!!          if (verbose) print*, i, j
          do k = 1, nz    ! all nodes in column are active
             nNodes = nNodes + 1   
             NodeID(k,i,j) = nNodes   ! unique index for each active node
             iNodeIndex(nNodes) = i
             jNodeIndex(nNodes) = j
             kNodeIndex(nNodes) = k
             xNode(nNodes) = xVertex(i,j)
             yNode(nNodes) = yVertex(i,j)
             zNode(nNodes) = stagusrf(i,j) - sigma(k)*stagthck(i,j)

             if (verbose .and. nNodes==ntest) then
                print*, ' '
                print*, 'i, j, k, n:', i, j, k, nNodes
                print*, 'sigma, stagusrf, stagthck:', sigma(k), stagusrf(i,j), stagthck(i,j)
                print*, 'x, y, z:', xNode(nNodes), yNode(nNodes), zNode(nNodes)
                print*, 'dx, dy, dz:', xVertex(i,j) - xVertex(i-1,j), &
                                       yVertex(i,j) - yVertex(i,j-1), &
                                       (sigma(k+1) - sigma(k)) * stagthck(i,j)
             endif

           enddo   ! k
        endif      ! active cell
    enddo          ! i
    enddo          ! j

    ! Identify the nodes of each element
    ! Assume that k increases from top to bottom 
    ! (Note: The numbering of nodes within each element is from bottom to top--
    !  this can be confusing)

    if (verbose) then
       print*, ' '
       print*, 'nActiveVertices =', nActiveVertices
       print*, 'nNodes =', nNodes
    endif

    NodeOnElementID(:,:,:,:) = 0  ! indices are (n,k,i,j)

    do j = 1+nhalo, ny-nhalo    ! loop over locally owned cells 
    do i = 1+nhalo, nx-nhalo
       if (active_cell(i,j)) then
          do k = 1, nz-1
             NodeOnElementID(1,k,i,j) = NodeID(k+1,i-1,j-1)
             NodeOnElementID(2,k,i,j) = NodeID(k+1,i,  j-1)
             NodeOnElementID(3,k,i,j) = NodeID(k+1,i,  j)
             NodeOnElementID(4,k,i,j) = NodeID(k+1,i-1,j)
             NodeOnElementID(5,k,i,j) = NodeID(k,  i-1,j-1)
             NodeOnElementID(6,k,i,j) = NodeID(k,  i,  j-1)
             NodeOnElementID(7,k,i,j) = NodeID(k,  i,  j)
             NodeOnElementID(8,k,i,j) = NodeID(k,  i-1,j)
          enddo  ! k
       endif     ! active cell
    enddo        ! i
    enddo        ! j

  end subroutine get_nodal_geometry

!****************************************************************************

  subroutine assemble_linear_system(nx, ny, nz,  nhalo,                &
                                    active_cell,                       &
                                    nCells,           nActiveCells,    &
                                    NodeID,           NodeOnElementID, &
                                    iNodeIndex,       jNodeIndex,     kNodeIndex,   &     
                                    xNode,            yNode,          zNode,        &
                                    nElements,        nNodes,          &
                                    uvel,             vvel,            &
                                    stagusrf,         flwafact,        &
                                    whichefvs,                         &
                                    Auu,              Auv,             &
                                    Avu,              Avv,             &
                                    bu,               bv,              &
                                    sia_factor,       ssa_factor)

    ! Assemble the stiffness matrix A and load vector b in the linear system Ax = b
 
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nz,                      &    ! number of vertical layers at which velocity is computed
       nhalo,                   &    ! number of halo layers
       nCells,                  &    ! number of horizontal cells
       nActiveCells                  ! number of cells with ice present

    logical, dimension(nx,ny), intent(in) ::  &
       active_cell                   ! true if element contains ice

    integer, dimension(nVertices*nz), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

    real(dp), dimension(nVertices*nz), intent(in) ::   &
       xNode, yNode, zNode  ! x, y and z coordinates of nodes

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
        NodeID             ! ID for each active node

    integer, dimension(nNodesPerElement, nz-1, nx, ny), intent(in) :: &
        NodeOnElementID    ! ID for each node of each element

!TODO - Are these needed?
    integer, intent(in) ::   &
       nElements,    &    ! number of elements
       nNodes             ! number of nodes associated with these elements

    real (dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       uvel, vvel             ! velocity components

    real(dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusrf           ! upper surface elevation on staggered grid

    real(dp), dimension(nz-1,nx,ny), intent(in) ::  &
       flwafact           ! temperature-based flow factor, 0.5 * A^(-1/n), 
                          ! used to compute the effective viscosity

    integer, intent(in) :: whichefvs      ! option for effective viscosity calculation 
                                          ! (calculate it or make it uniform)

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::  &
       bu, bv             ! load vector, divided into u and v components

    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    real(dp) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement) ::   &
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
    ! At the next step, the terms of the dense A matrices are put in a sparse matrix
    !  format suitable for iterative solutions.
    !----------------------------------------------------------------

    real(dp), dimension(nNodesPerElement, nNodesPerElement) ::   &   !
       Kuu,          &    ! element stiffness matrix, divided into 4 parts as shown below
       Kuv,          &    !  
       Kvu,          &    !
       Kvv                !    Kuu  | Kuv
                          !    _____|____          
                          !    Kvu  | Kvv
                          !         |
                          ! Kvu may not be needed if matrix is symmetric, but is included for now

    integer, dimension(nNodesPerElement) :: nid        ! ID for each node of an element

    real(dp), dimension(nNodesPerElement) ::     &
       x, y, z,         & ! Cartesian coordinates of nodes
       u, v,            & ! u and v at nodes
       s                  ! upper surface elevation at nodes

    real(dp) ::    &
       efvs               ! effective viscosity

    integer :: i, j, k, n, p

    if (verbose) then
       print*, ' '
       print*, 'In assemble_linear_system'
    endif

    ! Initialize global stiffness matrix and load vector

    Auu(:,:,:,:,:,:) = 0.d0
    Auv(:,:,:,:,:,:) = 0.d0
    Avu(:,:,:,:,:,:) = 0.d0
    Avv(:,:,:,:,:,:) = 0.d0

    bu(:,:,:) = 0.d0
    bv(:,:,:) = 0.d0

    ! Sum over elements

    do j = nhalo+1, ny-nhalo    ! loop over the cells owned by this processor
    do i = nhalo+1, nx-nhalo
       
     if (active_cell(i,j)) then

       do k = 1, nz-1    ! loop over elements in this column
                         ! assume k increases from upper surface to bed

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, 'i, j, k:', i, j, k
             print*, ' '
          endif

          ! Initialize element stiffness matrix
          Kuu(:,:) = 0.d0
          Kuv(:,:) = 0.d0
          Kvu(:,:) = 0.d0
          Kvv(:,:) = 0.d0
    
          ! get nodeID, spatial coordinates, velocity, and upper surface elevation for each node
          !TODO - Compute s(n) only once per column?

          do n = 1, nNodesPerElement

             nid(n) = NodeOnElementID(n,k,i,j)
             x(n) = xNode(nid(n))
             y(n) = yNode(nid(n))
             z(n) = zNode(nid(n))
             u(n) = uvel(kNodeIndex(nid(n)), iNodeIndex(nid(n)), jNodeIndex(nid(n)))
             v(n) = vvel(kNodeIndex(nid(n)), iNodeIndex(nid(n)), jNodeIndex(nid(n)))
             s(n) = stagusrf(iNodeIndex(nid(n)), jNodeIndex(nid(n)))

             if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, 'n, nid, x, y, z:', n, nid(n), x(n), y(n), z(n)
                print*, 's, u, v:', s(n), u(n), v(n)
             endif

          enddo   ! nodes per element

          ! Loop over quadrature points for this element
   
          !TODO - Think about how often to compute things.  Geometric quantities need to be
          !       computed only once per nonlinear solve, whereas efvs must be recomputed
          !       after each linear solve.

          do p = 1, nQuadPoints

          ! Evaluate the derivatives of the element basis functions
          ! at this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

             call get_basis_function_derivatives(nNodesPerElement, &
                                                 x(:),          y(:),          z(:),           &
                                                 dphi_dxr(:,p), dphi_dyr(:,p), dphi_dzr(:,p),  &
                                                 dphi_dx(:),    dphi_dy(:),    dphi_dz(:),     &
                                                 detJ , i, j, k, p                      )

             if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
!!                print*, ' '
!!                print*, 'Derivatives of basis functions, p =', p
                do n = 1, nNodesPerElement
!!                   print*, 'n, dphi_dx, dphi_dy, dphi_dz:', n, dphi_dx(n), dphi_dy(n), dphi_dz(n)
                enddo
             endif


!WHL - debug - Pass in i, j, k, and p for now

             call compute_effective_viscosity(whichefvs,        nNodesPerElement,             &
                                              dphi_dx(:),       dphi_dy(:),     dphi_dz(:),   &
                                              u(:),             v(:),                         & 
                                              flwafact(k,i,j),  efvs,                         &
                                              sia_factor,       ssa_factor,  i, j, k, p )

             if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
!                print*, ' '
!                print*, 'efvs:', efvs
             endif

             ! Increment the element stiffness matrix with the contribution from
             ! this quadrature point.

!WHL - debug - Pass in i, j, k, and p for now

             !TODO - Don't have to pass in nNodesPerElement

             call increment_element_matrix(nNodesPerElement,                          &
                                           wqp(p),           detJ,        efvs,       &
                                           dphi_dx(:),       dphi_dy(:),  dphi_dz(:), &
                                           Kuu(:,:),         Kuv(:,:),                &
                                           Kvu(:,:),         Kvv(:,:),                &
                                           sia_factor,       ssa_factor,  i, j, k, p )

             ! Increment the global load vector (i.e., the right-hand side of Ax = b)
             ! with the contribution from this quadrature point.

             !TO DO - Currently this is just the gravitational driving stress.
             !        Need to add the pressure term for floating ice.

             !TODO - Don't have to pass in nNodesPerElement

!WHL - debug - Pass in i, j, k, and p for now

             call increment_load_vector(nNodesPerElement, nid(:),         &
                                        iNodeIndex,       jNodeIndex,   kNodeIndex,  &
                                        wqp(p),           detJ,           &
                                        s(:),             phi(:,p),       &
                                        dphi_dx(:),       dphi_dy(:),     &
                                        bu(:,:,:),        bv(:,:,:), i, j, k, p)

             if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, ' '
                do n = 1, nNodesPerElement
                   print*, 'nid, bu, bv:', nid(n), bu(k,i,j), bv(k,i,j)
                enddo
             endif

          enddo   ! nQuadPoints

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Kuu:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kuu(n,:)
             enddo
             print*, 'Kuv:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kuv(n,:)
             enddo
             print*, 'Kvu:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kvu(n,:)
             enddo
             print*, 'Kvv:'
             do n = 1, nNodesPerElement
                write(6, '(8e12.4)') Kvv(n,:)
             enddo
          endif

          call check_symmetry_element_matrix(nNodesPerElement,  &
                                             Kuu, Kuv, Kvu, Kvv)

         ! If solving in parallel with Trilinos, we have the option at this point to 
         ! call a sum_into_global_matrix routine, passing one row at a time of the
         ! element matrix.  Trilinos will handle the rest.
         !
         ! For the serial SLAP solve, we first form a dense intermediate matrix, 
         ! and later put that matrix in sparse format.

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuu into Auu'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement, nNodes,  nid,         &
                                        iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                        Kuu,              Auu, i, j, k)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuv into Auv'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement, nNodes,  nid,         &
                                        iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                        Kuv,              Auv, i, j, k)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvu into Avu'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement, nNodes,  nid,         &
                                        iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                        Kvu,              Avu, i, j, k)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvv into Avv'
          endif

          call element_to_global_matrix(nx,           ny,          nz,          &
                                        nNodesPerElement, nNodes,  nid,         &
                                        iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                        Kvv,              Avv, i, j, k)

       enddo   ! nz  (loop over elements in this column)

     endif   ! active cell

    enddo      ! i
    enddo      ! j

  end subroutine assemble_linear_system

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
                dphi_dx , dphi_dy, dphi_dz     ! derivatives of basis functions at quad pt
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
    real (dp), dimension(3,3) :: prod     ! Jac * Jinv (should be identity matrix)

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

!TODO - Check this code carefully.

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, 'In get_basis_function_derivatives: i, j, k, p =', i, j, k, p
    endif

    do n = 1, nNodesPerElement
       if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'n, x, y, z:', n, xNode(n), yNode(n), zNode(n)
          print*, 'dphi_dxyz:', dphi_dxr(n), dphi_dyr(n), dphi_dzr(n) 
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

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
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

    if (verbose_Jac .and. i==itest .and. j==jtest .and. k==ktest) then
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
    enddo

    ! debug - Check that the sum of dphi_dx, etc. is close to zero  
    !TODO - Think about why this should be the case

    if (abs(sum(dphi_dx)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dx > 0'
        print*, 'dphi_dx =', dphi_dx(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs(sum(dphi_dy)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dy > 0'
        print*, 'dphi_dy =', dphi_dy(:)
        print*, 'i, j, k, p =', i, j, k, p
        stop
    endif

    if (abs(sum(dphi_dz)) > 1.d-12) then
        print*, 'stopping, sum over basis functions of dphi_dz > 0'
        print*, 'dphi_dz =', dphi_dz(:)
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

!TODO - Change to different variable name?
    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       uvel, vvel      ! current guess for velocity field at each node of element

    real(dp), intent(in) ::  &
       flwafact        ! temperature-based flow factor for this element, 0.5 * A^(-1/n)

    real(dp), intent(out) ::   &
       efvs            ! effective viscosity at this quadrature point

    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !----------------------------------------------------------------
    ! Local parameters
    !----------------------------------------------------------------

    real(dp), parameter ::   &
       effstrainsq_min = (1.0d-20 *tim0)**2,  &! minimum value of squared effective strain rate
       p_effstr = (1.d0 - real(gn,dp)) / (2.d0*real(gn,dp))    ! exponent in efvs relation

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::            &
       du_dx, dv_dx,       & ! strain rate components
       du_dy, dv_dy,       &
       du_dz, dv_dz,       &
       effstrainsq           ! square of effective strain rate
        
    integer :: n

    select case(whichefvs)

    case(0)      ! calculate effective viscosity based on effective strain rate

       du_dx = 0.d0
       dv_dx = 0.d0
       du_dy = 0.d0
       dv_dy = 0.d0
       du_dz = 0.d0
       dv_dz = 0.d0

       ! Compute strain rate components at this quadrature point

       do n = 1, nNodesPerElement

          du_dx = du_dx + dphi_dx(n)*uvel(n)
          dv_dx = dv_dx + dphi_dx(n)*vvel(n)

          du_dy = du_dy + dphi_dy(n)*uvel(n)
          dv_dy = dv_dy + dphi_dy(n)*vvel(n)

          du_dz = du_dz + dphi_dz(n)*uvel(n)
          dv_dz = dv_dz + dphi_dz(n)*vvel(n)

       enddo   ! nNodesPerElement

       ! compute effective strain rate at this quadrature point (PGB 2012, eq. 3 and 9)

       effstrainsq = effstrainsq_min                                       &
                   + ssa_factor * (du_dx**2 + dv_dy**2 + du_dx*dv_dy       &
                                   + 0.25d0 * (dv_dx + du_dy)**2)          &
                   + sia_factor * 0.25d0 * (du_dz**2 + dv_dz**2) * (len0/thk0)**2  !TODO - remove scaling?

       ! compute effective viscosity (PGB 2012, eq. 4)

       efvs = flwafact * effstrainsq**p_effstr
       !TODO - Multiply by scyr/tim0 / tau0?

       if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'Compute efvs, i, j, k, p =', i, j, k, p
          print*, 'du/dx, du/dy, du/dz =', du_dx, du_dy, du_dz
          print*, 'dv/dx, dv/dy, dv/dz =', dv_dx, dv_dy, dv_dz
          print*, 'effstrainsq =', effstrainsq
          print*, 'flwafact =', flwafact
          print*, 'efvs =', efvs
       endif

    case(1)      ! set the effective viscosity to a multiple of the flow factor, 0.5*A^(-1/n)
                 ! (actually, to the flow factor itself)
                 !TODO - Verify that this value is OK

       efvs = flwafact

       if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'Set efvs = flwafact, i, j, k, p =', i, j, k, p
          print*, 'efvs =', efvs       
       endif

!WHL - adding an option to set viscosity to the same value everywhere
!TODO - Check the scaling.  If dividing by tau0, need to make sure rhoi*grav is not on RHS

    case(2)

       efvs = 1.d7    ! Steve recommends 10^6 to 10^7 Pa yr

!WHL - This is the glam-type scaling
!!       efvs = efvs * scyr/tim0 / tau0   ! tau0 = rhoi*grav*thk0

       if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
          print*, ' '
          print*, 'Set efvs = constant, i, j, k, p =', i, j, k, p
          print*, 'efvs =', efvs       
       endif

    end select

  end subroutine compute_effective_viscosity

!****************************************************************************

!WHL - debug - Pass in i, j, k, and p for now

  subroutine increment_element_matrix(nNodesPerElement,                &
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
    ! This subroutine will work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

!TODO - Add optional arguments that reduce the solver to SIA or SSA.

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

    real(dp), intent(in) ::  &
       sia_factor,      & ! = 1. if SIA terms are included, else = 0.
       ssa_factor         ! = 1. if SSA terms are included, else = 0.

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    if (verbose .and. ii==itest .and. jj==jtest .and. k==ktest) then
       print*, ' '
       print*, 'Increment element matrix, p =', p
    endif

    ! Increment the element stiffness matrices for the first-order Blatter-Pattyn equations
    ! The terms in parentheses can be derived from PGB 2012, eq. 13 and 15.
    ! TODO: Check factor of 2 in eq. 15. I think it's been absorbed into parenthetical terms. 
    
    do j = 1, nNodesPerElement      ! columns of K
       do i = 1, nNodesPerElement   ! rows of K

       if (verbose .and. ii==itest .and. jj==jtest .and. k==ktest .and. p==ptest &
                   .and.  i==1     .and.  j==1) then
          print*, 'efvs, wqp, detJ/vol0 =', efvs, wqp, detJ/vol0
          print*, 'dphi_dz(1) =', dphi_dz(1)
          print*, 'dphi_dx(1) =', dphi_dx(1)
          print*, 'Kuu dphi/dz increment(1,1) =', efvs*wqp*detJ/vol0*dphi_dz(1)*dphi_dz(1)
          print*, 'Kuu dphi/dx increment(1,1) =', efvs*wqp*detJ/vol0*4.d0*dphi_dx(1)*dphi_dx(1)
       endif

!WHL - Note volume scaling such that detJ/vol0 is close to unity
!       Compare SIA matrix to total matrix

          Kuu(i,j) = Kuu(i,j) + efvs*wqp*detJ/vol0 *                                               &
                              ( ssa_factor * (4.d0*dphi_dx(j)*dphi_dx(i) + dphi_dy(j)*dphi_dy(i))  &
                              + sia_factor * (dphi_dz(j)*dphi_dz(i)) )

          Kuv(i,j) = Kuv(i,j) + efvs*wqp*detJ/vol0 *                                               &
                                ssa_factor * (2.d0*dphi_dx(j)*dphi_dy(i) + dphi_dy(j)*dphi_dx(i))

          Kvu(i,j) = Kvu(i,j) + efvs*wqp*detJ/vol0 *                                               &
                                ssa_factor * (2.d0*dphi_dy(j)*dphi_dx(i) + dphi_dx(j)*dphi_dy(i))

          Kvv(i,j) = Kvv(i,j) + efvs*wqp*detJ/vol0 *                                               &
                              ( ssa_factor * (4.d0*dphi_dy(j)*dphi_dy(i) + dphi_dx(j)*dphi_dx(i))  &
                              + sia_factor * (dphi_dz(j)*dphi_dz(i)) )

       enddo  ! i (rows)
    enddo     ! j (columns)

  end subroutine increment_element_matrix

!****************************************************************************

!WHL - Pass in i, j, k, p for now

  subroutine increment_load_vector(nNodesPerElement, nid,           &
                                   iNodeIndex,       jNodeIndex,   kNodeIndex,  &
                                   wqp,              detJ,          &
                                   usrf,             phi,           &
                                   dphi_dx,          dphi_dy,       &
                                   bu,               bv,     ii, jj, kk, p)

    integer, intent(in) :: ii, jj, kk, p

    integer, intent(in) :: nNodesPerElement

    integer, dimension(nNodesPerElement), intent(in) ::    &
       nid                ! IDs for the nodes of this element

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

    real(dp), intent(in) ::   &
       wqp,              &! weight for this quadrature point
       detJ               ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
       usrf               ! upper surface elevation at each node

    real(dp), dimension(nNodesPerElement), intent(in) ::   & 
       phi                ! nodal basis functions, evaluated at quad pts

    real(dp), dimension(nNodesPerElement), intent(in) ::   &
       dphi_dx, dphi_dy   ! derivatives of basis functions, evaluated at quad pts

    real(dp), dimension(:,:,:), intent(inout) ::  &
       bu, bv             ! load vector, divided into u and v components

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp) ::    &
       ds_dx, ds_dy       ! horizontal gradient of upper surface elevation

    integer :: i, j, k, n

    ! Evaluate ds_dx and ds_dy at this quadrature point

    ds_dx = 0.d0
    ds_dy = 0.d0

    do n = 1, nNodesPerElement
       ds_dx = ds_dx + dphi_dx(n) * usrf(n)
       ds_dy = ds_dy + dphi_dy(n) * usrf(n)
    enddo

    if (verbose .and. ii==itest .and. jj==jtest .and. kk==ktest) then
       print*, ' '
       print*, 'Increment load vector, p =', p
       print*, 'ds/dx, ds/dy =', ds_dx, ds_dy       
       print*, 'detJ/vol0 =', detJ/vol0
       print*, 'detJ/vol0* (ds/dx, ds/dy) =', detJ/vol0*ds_dx, detJ/vol0*ds_dy
    endif

    ! Increment the load vector with the gravitational contribution from
    ! this quadrature point.

    ! Is the indexing correct?  E.g., no n index for ds_dx?

!       Note: glam code does not include rhoi*grav on RHS.

    do n = 1, nNodesPerElement

       i = iNodeIndex(nid(n))
       j = jNodeIndex(nid(n))
       k = kNodeIndex(nid(n))
         
       bu(k,i,j) = bu(k,i,j) - rhoi*grav * wqp * detJ/vol0 * ds_dx * phi(n)
       bv(k,i,j) = bv(k,i,j) - rhoi*grav * wqp * detJ/vol0 * ds_dy * phi(n)

       if (verbose .and. ii==itest .and. jj==jtest .and. kk==ktest .and. p==ptest) then
!          print*, ' '
          print*, 'n, nid(n), phi(n), delta(bu), delta(bv):', n, nid(n), phi(n), &
                rhoi*grav*wqp*detJ/vol0 * ds_dx * phi(n), &
                rhoi*grav*wqp*detJ/vol0 * ds_dy * phi(n)
       endif

    enddo

  end subroutine increment_load_vector

!****************************************************************************

  subroutine dirichlet_boundary_conditions(nx,  ny,  nz,                      &
                                           active_vertex,   umask_dirichlet,  &
                                           nNodes,          NodeID,           &
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
       nz                       ! number of vertical levels where velocity is computed

    logical, dimension(nx-1,ny-1), intent(in) ::  &
       active_vertex       ! true for active vertices (columns where velocity is computed)

    logical, dimension(nz,nx-1,ny-1), intent(in) ::  &
       umask_dirichlet     ! Dirichlet mask for velocity (if true, u = 0)

    integer, intent(in) ::   &
       nNodes              ! number of nodes where velocity is computed

    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       NodeID              ! ID for each active node

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(inout) ::   &
       Auu, Auv,    &      ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       bu, bv              ! assembled load vector, divided into 2 parts

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------
    
    integer :: i, j, k     ! Cartesian indices of nodes
    integer :: iA, jA, kA  ! i, j, and k offsets of neighboring nodes 

!TODO - don't need row and col?
    integer :: row, col    ! row and column indices of assembled global matrix

    do j = 1, ny-1
       do i = 1, nx-1
          if (active_vertex(i,j)) then
             do k = 1, nz
                if (umask_dirichlet(k,i,j)) then

                   ! determine the row of the assembled global matrix associated with this node
!!                   row = NodeID(k,i,j)       
                   
                   ! loop through potential nonzero matrix values in the row associated with this vertex
                   ! Note: row = NodeID(k,i,j) if we were forming a matrix with one row per node

                   do jA = -1,1
                      do iA = -1,1
                         do kA = -1,1

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

                            Auu( kA,  iA,  jA, k,    i,    j   ) = 0.d0
                            Auu(-kA, -iA, -jA, k+kA, i+iA, j+jA) = 0.d0
                            Avv( kA,  iA,  jA, k,    i,    j   ) = 0.d0
                            Avv(-kA, -iA, -jA, k+kA, i+iA, j+jA) = 0.d0
                            Auv( kA,  iA,  jA, k,    i,    j   ) = 0.d0
                            Auv(-kA, -iA, -jA, k+kA, i+iA, j+jA) = 0.d0
                            Avu( kA,  iA,  jA, k,    i,    j   ) = 0.d0
                            Avu(-kA, -iA, -jA, k+kA, i+iA, j+jA) = 0.d0

                         enddo     ! kA
                      enddo        ! iA
                   enddo           ! jA

                   ! put back 1's on the main diagonal
!                   Auu(0,0,0,row) = 1.d0
!                   Avv(0,0,0,row) = 1.d0
                   Auu(0,0,0,k,i,j) = 1.d0
                   Avv(0,0,0,k,i,j) = 1.d0

                   ! fix u = v = 0 for this node by zeroing the rhs  
!                   bu(row) = 0.d0            
!                   bv(row) = 0.d0
                   bu(k,i,j) = 0.d0            
                   bv(k,i,j) = 0.d0

                endif    ! umask_dirichlet
             enddo       ! k
          endif          ! active_vertex
       enddo             ! i
    enddo                ! j

  end subroutine dirichlet_boundary_conditions

!****************************************************************************

!WHL - Pass in ii, jj, kk for now
  subroutine element_to_global_matrix(nx,           ny,          nz,          &
                                      nNodesPerElement, nNodes,  nid,         &
                                      iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                      Kmat,         Amat, ii, jj, kk)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! So this subroutine is called four times per element.

    integer, intent(in) :: ii, jj, kk

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    integer, intent(in) ::   &
       nNodesPerElement,     &
       nNodes            ! total number of active nodes !TODO - Is this needed?
       
    integer, dimension(nNodesPerElement), intent(in) ::  &
       nid               ! node ID for each node of this element

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: i, j, k
    integer :: iA, jA, kA
    integer :: m, n

    if (verbose .and. ii==itest .and. jj==jtest .and. kk==ktest) then
       print*, 'First row of K:'
       write(6, '(8e12.4)') Kmat(1,:)
    endif

!TODO - Switch loops or switch order of indices in K(m,n)?
!       Inner index m is currently the outer loop

    do m = 1, nNodesPerElement       ! rows of K

       i = iNodeIndex(nid(m))
       j = jNodeIndex(nid(m))
       k = kNodeIndex(nid(m))

       do n = 1, nNodesPerElement    ! columns of K

          kA = kshift(m,n)           ! k index of A into which K(m,n) is summed;
          iA = ishift(m,n)           ! similarly for i and j indices 
          jA = jshift(m,n)           ! these indices can take values -1, 0 and 1

          Amat(kA,iA,jA,k,i,j) = Amat(kA,iA,jA,k,i,j) + Kmat(m,n)

          if (verbose .and. nid(m)==ntest) then
!              print*, ' '
!              print*, 'Contribution to nid', nid(m), 'from element i, j, k =', ii, jj, kk
!              print*, 'm, n, Kmat(m,n) =', m, n, Kmat(m,n)
!              print*, 'kA, iA, jA', kA, iA, jA
!              print*, 'New value of A =', Amat(kA,iA,jA,k,i,j)
          endif

       enddo     ! n
    enddo        ! m

  end subroutine element_to_global_matrix

!****************************************************************************

  subroutine solver_preprocess(nx,           ny,          nz,          &
                               nNodes,       NodeID,                   &
                               iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                               Auu,          Auv,   &
                               Avu,          Avv,   &
                               bu,           bv,    &
                               uvel,         vvel,  &
                               matrix_order, nNonzero,  &
                               matrix,       rhs,   &
                               answer)

    !----------------------------------------------------------------
    ! Using the intermediate matrices (Auu, Auv, Avu, Avv), load vectors (bu, bv),
    ! and velocity components (uvel, vvel), form the matrix and the rhs and answer
    ! vectors in the desired sparse matrix format.
    !
    ! The matrix is formed in ascending row order, so it can easily be transformed
    ! to compressed sparse row (CSR) format without further sorting.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels at which velocity is computed

!TODO - This is no longer needed
    integer, intent(in) ::   &
       nNodes            ! number of active nodes

    integer, dimension(:,:,:), intent(in) ::  &
       NodeID            ! ID for each active node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = node and its 2 neighbors in z direction 
                          ! 2nd dimension = node and its 2 neighbors in x direction 
                          ! 3rd dimension = node and its 2 neighbors in y direction 
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

    integer :: i, j, k, iA, jA, kA, n, ct

    integer :: rowA, colA   ! row and column of A submatrices (order = nNodes)
    integer :: row, col     ! row and column of sparse matrix (order = 2*nNodes) 

    real(dp) :: val         ! value of matrix coefficient
    
    ! Set the nonzero coefficients of the sparse matrix 

    ct = 0

    do rowA = 1, nNodes

       i = iNodeIndex(rowA)
       j = jNodeIndex(rowA)
       k = kNodeIndex(rowA)

       ! Load the nonzero values associated with Auu and Auv
       ! These are assigned a value of matrix%row = 2*rowA - 1

       do jA = -1, 1
       do iA = -1, 1
       do kA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)     &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx)         &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node

            !TODO - Change order so that matrix is filled in ascending row order?

             ! Auu
             val = Auu(kA,iA,jA,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA - 1                    
                matrix%val(ct) = val
             endif

             ! Auv 
             val = Auv(kA,iA,jA,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA - 1
                matrix%col(ct) = 2*colA
                matrix%val(ct) = val
             endif

          endif     ! i+iA, j+jA, and k+kA in bounds

       enddo        ! kA
       enddo        ! iA
       enddo        ! jA

       ! Load the nonzero values associated with Avu and Avv
       ! These are assigned a value of matrix%row = 2*rowA

       do jA = -1, 1
       do iA = -1, 1
       do kA = -1, 1

          if ( (k+kA >= 1 .and. k+kA <= nz)     &
                          .and.                     &
               (i+iA >= 1 .and. i+iA <= nx)         &
                          .and.                     &
               (j+jA >= 1 .and. j+jA <= ny-1) ) then

             colA = NodeID(k+kA, i+iA, j+jA)   ! ID for neighboring node

             ! Avu 
             val = Avu(kA,iA,jA,k,i,j)
             if (val /= 0.d0) then
                ct = ct + 1
                matrix%row(ct) = 2*rowA
                matrix%col(ct) = 2*colA - 1
                matrix%val(ct) = val
             endif

             ! Avv 
             val = Avv(kA,iA,jA,k,i,j)
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

    if (verbose) then
       print*, 'solver_preprocess'
       print*, 'order, nonzeros =', matrix%order, matrix%nonzeros
    endif

    ! Initialize the answer vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodes
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       answer(2*n-1) = uvel(k,i,j)
       answer(2*n)   = vvel(k,i,j)

       if (verbose .and. n==ntest) then
          print*, ' '
          print*, 'n, initial uvel =', n, answer(2*n-1)
          print*, 'n, initial vvel =', n, answer(2*n)
       endif

    enddo

    ! Set the rhs vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodes
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       rhs(2*n-1) = bu(k,i,j)
       rhs(2*n)   = bv(k,i,j)

       if (verbose .and. n==ntest) then
          print*, ' '
          print*, 'n, initial bu =', n, rhs(2*n-1)
          print*, 'n, initial bv =', n, rhs(2*n)
       endif

    enddo

  end subroutine solver_preprocess

!****************************************************************************

  subroutine solver_postprocess(nNodes,       answer,                   &
                                iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                uvel,         vvel)

  ! Extract the velocities from the solution vector.
                                            
!    use parallel

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) :: nNodes        ! number of active nodes

    real(dp), dimension(:), intent(in) ::   &
       answer             ! velocity solution vector

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(:,:,:), intent(inout) ::   &
       uvel, vvel         ! u and v components of velocity

    integer :: i, j, k, n

       if (verbose) then
          print*, ' '
          print*, 'solver postprocess: n, i, j, k, u, v'
       endif

    do n = 1, nNodes
       i = iNodeIndex(n)
       j = jNodeIndex(n)
       k = kNodeIndex(n)

       uvel(k,i,j) = answer(2*n-1)
       vvel(k,i,j) = answer(2*n)

       if (verbose) then
!          print*, n, i, j, k, uvel(k,i,j), vvel(k,i,j)
       endif

    enddo

  end subroutine solver_postprocess

!****************************************************************************

  subroutine compute_residual_vector_structured(nx,          ny,            &
                                                nz,          active_vertex, &
                                                Auu,         Auv,           &
                                                Avu,         Avv,           &
                                                bu,          bv,            &
                                                uvel,        vvel,          &
                                                resid_vec_u, resid_vec_v,   &
                                                L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in structured (x/y/z) format.

    use parallel

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 3 (node and its 2 neighbors in z direction)
                              ! 2nd dimension = 3 (node and its 2 neighbors in x direction)
                              ! 3rd dimension = 3 (node and its 2 neighbors in y direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv              ! assembled load (rhs) vector, divided into 2 parts

   real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       uvel, vvel          ! u and v components of velocity

    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::   &
       resid_vec_u,      & ! residual vector, divided into 2 parts
       resid_vec_v

    real(dp), intent(out) ::    &
       L2_norm             ! L2 norm of residual vector

    integer :: i, j, k
    integer :: iA, jA, kA 

!WHL - debug
    integer :: nid, islap, jslap

    if (verbose) then
       print*, ' '
       i = itest
       j = jtest
       k = ktest
       print*, 'i, j, k, u answer, u rhs:', i, j, k, uvel(k,i,j), bu(k,i,j)
       print*, 'i, j, k, v answer, v rhs:', i, j, k, vvel(k,i,j), bv(k,i,j)
    endif

    resid_vec_u(:,:,:) = 0.d0
    resid_vec_v(:,:,:) = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1

       !TODO - Use indirect addressing to avoid an 'if' here?

       if (active_vertex(i,j)) then

          do k = 1, nz

             do jA = -1,1
             do iA = -1,1
             do kA = -1,1

                if ( (k+kA >= 1 .and. k+kA <= nz)     &
                                .and.                     &
                     (i+iA >= 1 .and. i+iA <= nx-1)         &
                             .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   resid_vec_u(k,i,j) = resid_vec_u(k,i,j)                      & 
                                      + Auu(kA,iA,jA,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                      + Auv(kA,iA,jA,k,i,j)*vvel(k+kA,i+iA,j+jA)

                   resid_vec_v(k,i,j) = resid_vec_v(k,i,j)                      &
                                      + Avu(kA,iA,jA,k,i,j)*uvel(k+kA,i+iA,j+jA)  &
                                      + Avv(kA,iA,jA,k,i,j)*vvel(k+kA,i+iA,j+jA)

                endif   ! in bounds

             enddo   ! kA
             enddo   ! iA
             enddo   ! jA

          enddo   ! k

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j


    if (verbose) then
!       print*, ' '
!       print*, 'Compute residual:'
!       print*, 'k, i, j, Axu, bu:'
    endif

    L2_norm = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1
       if (active_vertex(i,j)) then
          do k = 1, nz

!WHL - debug
!             print*, k, i, j, resid_vec_u(k,i,j), bu(k,i,j)

             resid_vec_u(k,i,j) = resid_vec_u(k,i,j) - bu(k,i,j)
             resid_vec_v(k,i,j) = resid_vec_v(k,i,j) - bv(k,i,j)
             L2_norm = L2_norm + resid_vec_u(k,i,j)*resid_vec_u(k,i,j)  &
                               + resid_vec_v(k,i,j)*resid_vec_v(k,i,j)
          enddo  ! k
       endif     ! active vertex
    enddo        ! i
    enddo        ! j

!TODO - Parallel sum
!    L2_norm = parallel_reduce_sum(L2_norm)

    L2_norm = sqrt(L2_norm)

  end subroutine compute_residual_vector_structured

!****************************************************************************

  subroutine compute_residual_vector(matrix,    answer,   rhs, &
                                     resid_vec, L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.
    ! This subroutine assumes that the matrix is stored in triad (row/col/val) format.

    use parallel

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

!WHL - debug
    integer :: nid, islap, jslap

    if (verbose) then
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

!TODO - Parallel sum
!    L2_norm = parallel_reduce_sum(L2_norm)

    L2_norm = sqrt(L2_norm)

  end subroutine compute_residual_vector

!****************************************************************************

  subroutine compute_residual_velocity(nhalo,  whichresid, &
                                       uvel,   vvel,        &
                                       usav,   vsav,        &
                                       resid_velo)

    use parallel

    integer, intent(in) ::   &
       nhalo,           & ! number of layers of halo cells
       whichresid         ! option for method to use when calculating residual

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

    case(1)   ! max speed difference, excluding the bed

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
       
       resid_velo = parallel_reduce_max(resid_velo)

    case(2)   ! mean relative speed difference

       count = 0

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

       call not_parallel(__FILE__, __LINE__)
       !TODO - parallel reduction for mean

   case default    ! max speed difference, including basal speeds  (case 0 or 3)

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

    ! check that Kuu = (Kuu)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kuu(i,j) - Kuu(j,i)) > eps12) then
             print*, 'Kuu is not symmetric'
             print*, 'i, j, Kuu(i,j), Kuu(j,i):', i, j, Kuu(i,j), Kuu(j,i)
             stop
          endif    
       enddo
    enddo

    ! check that Kvv = (Kvv)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kvv(i,j) - Kvv(j,i)) > eps12) then
             print*, 'Kvv is not symmetric'
             print*, 'i, j, Kvv(i,j), Kvv(j,i):', i, j, Kvv(i,j), Kvv(j,i)
             stop
          endif    
       enddo
    enddo

    ! Check that Kuv = (Kvu)^T

    do j = 1, nNodesPerElement
       do i = 1, nNodesPerElement
          if (abs(Kuv(i,j) - Kvu(j,i)) > eps12) then
             print*, 'Kuv .ne. (Kvu)^T'
             print*, 'i, j, Kuv(i,j), Kvu(j,i):', i, j, Kuv(i,j), Kvu(j,i)
             stop
          endif    
       enddo
    enddo

  end subroutine check_symmetry_element_matrix

!****************************************************************************

  subroutine check_symmetry_assembled_matrix(nx,         ny,         nz,        &
                                             nNodes,     NodeID,    &
                                             iNodeIndex, jNodeIndex, kNodeIndex,  &
                                             Auu, Auv, Avu, Avv)

    !------------------------------------------------------------------
    ! Check that the assembled stiffness matrix is symmetric.
    ! Roughly speaking, this is true provided that (1) Auu = (Auu)^T
    !                                              (2) Avv = (Avv)^T
    !                                              (3) Auv = (Avu)^T
    ! However, the A matrices are assembled in a dense fashion to save storage.
    ! The row of A corresponds to the row of the full matrix (which is not assembled). 
    !  Using the indexing scheme we can determine the column of the full matrix, and then
    !  we can compare A(col,row) to A(row,col) for the full matrix.
    !
    ! This check should not be needed for production runs with a well-tested code,
    !  but is included for now to help with debugging.
    !------------------------------------------------------------------    

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nz                       ! number of vertical levels where velocity is computed

!TODO - not needed?
    integer, intent(in) :: nNodes     ! number of active nodes = number of rows of Auu, etc.

    integer, dimension(:,:,:), intent(in) ::  &
       NodeID             ! ID for each active node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(in) ::   &
             Auu, Auv, Avu, Avv     ! components of assembled stiffness matrix
                                    !
                                    !    Auu  | Auv
                                    !    _____|____          
                                    !    Avu  | Avv
                                    !         |

    integer :: i, j, k, iA, jA, kA, row, col

    real(dp) :: val1, val2          ! values of matrix coefficients

    real(dp) :: maxdiff

    ! check that Auu = (Auu)^T

    maxdiff = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1
    do k = 1, nz

       if (NodeID(k,i,j) > 0) then    ! node is active

          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1

             val1 = Auu(kA, iA, jA, k, i, j)      ! value of Auu(row,col)

             if (val1 /= 0.d0) then

                val2 = Auu(-kA, -iA, -jA, k+kA, i+iA, j+jA)   ! value of Auu(col,row)

                if (abs(val2 - val1) > eps12) then
                   print*, 'Auu is not symmetric: i, j, k, iA, jA, kA =', i, j, k, iA, jA, kA
                   print*, 'row, col, Auu(row,col), Auu(col,row):', &
                            row, col, val1, val2
!!                      stop
                endif

!WHL - debug
                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)

             endif   ! val1 /= 0.
                
          enddo        ! kA
          enddo        ! iA
          enddo        ! jA

       endif        ! NodeID > 0

    enddo        ! k
    enddo        ! i
    enddo        ! j

!WHL - debug
    if (verbose) print*, 'Auu: max difference from symmetry =', maxdiff

    ! check that Avv = (Avv)^T

    !TODO - The following code is basically identical to the above code for Auu.
    !       Combine into one code chunk called twice?

    maxdiff = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1
    do k = 1, nz

       if (NodeID(k,i,j) > 0) then    ! node is active

          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1

             val1 = Avv(kA, iA, jA, k, i, j)      ! value of Avv(row,col)

             if (val1 /= 0.d0) then

                val2 = Avv(-kA, -iA, -jA, k+kA, i+iA, j+jA)   ! value of Avv(col,row)

                if (abs(val2 - val1) > eps12) then
                   print*, 'Avv is not symmetric'
                   print*, 'row, col, Avv(row,col), Avv(col,row):', &
                            row, col, val1, val2
!!                   stop
                endif

                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             endif   ! val1 /= 0.
                
          enddo        ! kA
          enddo        ! iA
          enddo        ! jA

       endif        ! NodeID > 0

    enddo        ! k
    enddo        ! i
    enddo        ! j


!WHL - debug
    if (verbose) print*, 'Avv: max difference from symmetry =', maxdiff

    ! Check that Auv = (Avu)^T

    maxdiff = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1
    do k = 1, nz

       if (NodeID(k,i,j) > 0) then    ! node is active

          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1

             val1 = Auv(kA, iA, jA, k, i, j)      ! value of Auv(row,col)

             if (val1 /= 0.d0) then

                val2 = Avu(-kA, -iA, -jA, k+kA, i+iA, j+jA)   ! value of Avu(col,row)

                if (abs(val2 - val1) > eps12) then
                   print*, 'Auv is not equal to (Avu)^T'
                   print*, 'row, col, Auv(row,col), Avu(col,row):', &
                            row, col, val1, val2
!!                   stop
                endif

                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             endif   ! val1 /= 0.
                
          enddo        ! kA
          enddo        ! iA
          enddo        ! jA

       endif        ! NodeID > 0

    enddo        ! k
    enddo        ! i
    enddo        ! j

!WHL - debug
    if (verbose) print*, 'Auv v. (Avu)^T: max difference from symmetry =', maxdiff

    ! Check that Avu = (Auv)^T
    ! This check may not be necessary given that we've already checked that Auv = (Avu)^T.
    !  Just being extra careful.

    maxdiff = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1
    do k = 1, nz

       if (NodeID(k,i,j) > 0) then    ! node is active

          do jA = -1, 1
          do iA = -1, 1
          do kA = -1, 1

             val1 = Avu(kA, iA, jA, k, i, j)      ! value of Avu(row,col)

             if (val1 /= 0.d0) then

                val2 = Auv(-kA, -iA, -jA, k+kA, i+iA, j+jA)   ! value of Auv(col,row)

                if (abs(val2 - val1) > eps12) then
                   print*, 'Avu is not equal to (Auv)^T'
                   print*, 'row, col, Avu(row,col), Auv(col,row):', &
                            row, col, val1, val2
!!                   stop
                endif

                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             endif   ! val1 /= 0.
                
          enddo        ! kA
          enddo        ! iA
          enddo        ! jA

       endif        ! NodeID > 0

    enddo        ! k
    enddo        ! i
    enddo        ! j

!WHL - debug
    if (verbose) print*, 'Avu v. (Auv)^T: max difference from symmetry =', maxdiff

  end subroutine check_symmetry_assembled_matrix

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
    endif

    ! Solve the linear matrix problem

    if (whichsparse == STANDALONE_PCG_UNSTRUC) then

       if (verbose_test) then
          print*, 'Call standalone cism pcg solver'
       endif

       call cism_sparse_pcg_solve(matrix,  rhs,     answer,   &
                                  err,     niters,  verbose)

    else  ! SLAP solver

       if (verbose_test) then
          print*, 'Call sparse_easy_solve, whichsparse =', whichsparse
       endif

       call sparse_easy_solve(matrix, rhs,    answer,  &
                              err,    niters, whichsparse)

    endif

    if (verbose_test) then
       print*, ' '
       print*, 'answer =', answer
       print*, 'err =', err       
       print*, 'niters =', niters
    endif

    stop

  end subroutine solve_test_matrix

!****************************************************************************


!*************************************************

!WHL - Rename and revise this subroutine
!      Call from a k loop to avoid writing 3D version?

  subroutine calc_eps2_2d (ewn,       nsn,   &
                           dew,       dns,   &
                           uvel,      vvel,  &
                           thckmask,  eps2)

  !---------------------------------------------------------------------------
  ! Calculate the effective strain rate for the 2D shallow shelf approximation
  !
  ! eps^2 = (du/dx)^2 + (dv/dy)^2 + (du/dx + dv/dy)^2 + (1/2)*(du/dy + dv/dx)^2
  !
  !       =       2 * [(du/dx)^2 + (dv/dy)^2 + (du/dx*dv/dy)] + du/dy*dv/dx
  !         + (1/2) * [(dv/dx)^2 + (du/dy)^2]
  !
  ! Each component in this sum is an integrated average over the grid cell,
  ! given a bilinear reconstruction of u and v as a function of the four
  ! nodal values.
  !---------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  ! in-out variables
  !---------------------------------------------------------------------------

  integer, intent(in) ::  &
      ewn, nsn                  ! horizontal grid dimensions
      
  real(dp), intent(in) ::   &
      dew, dns                  ! E-W and N-S grid lengths

  real(dp), dimension(:,:), intent(in) ::    &
      uvel, vvel                ! E-W and N-S velocity components
 
  integer, dimension(:,:), intent(in) ::     &
      thckmask                  ! = 1 where ice is present, else = 0

  real(dp), dimension(:,:), intent(out) ::   &
       eps2                     ! square of effective strain rate for SSA

  !---------------------------------------------------------------------------
  ! local variables
  !---------------------------------------------------------------------------

  integer ::  i, j           ! horizontal indices

  real(dp) ::    &
     fxx, fyy, fxy,         &! geometric factors
     dun, dus, due, duw,    &! difference in uvel between adjacent nodes
     dvn, dvs, dve, dvw,    &! difference in vvel between adjacent nodes
     dudx2, dvdx2,          &! components of eps2
     dudy2, dvdy2,          &
     dudxdvdy, dvdxdudy

  !---------------------------------------------------------------------------
  ! calculate eps2
  !---------------------------------------------------------------------------

!WHL - may be able to speed this up by setting duw(i,j) = due(i-1,j)?

  fxx = 1.d0 / (dew*dew)
  fyy = 1.d0 / (dns*dns)
  fxy = 1.d0 / (dew*dns) 

  do j = 1, nsn
  do i = 1, ewn

     if (thckmask(i,j)==1) then

        dun = uvel(i,j)   - uvel(i-1,j)
        dus = uvel(i,j-1) - uvel(i-1,j-1)
        due = uvel(i,j)   - uvel(i,j-1)
        duw = uvel(i-1,j) - uvel(i-1,j-1)

        dvn = vvel(i,j)   - vvel(i-1,j)
        dvs = vvel(i,j-1) - vvel(i-1,j-1)
        dve = vvel(i,j)   - vvel(i,j-1)
        dvw = vvel(i-1,j) - vvel(i-1,j-1)

        dudx2 = fxx *  ( dun*dun + dus*dus + dun*dus) / 3.d0
        dvdx2 = fxx *  ( dvn*dvn + dvs*dvs + dvn*dvs) / 3.d0

        dudy2 = fyy *  ( due*due + duw*duw + due*duw) / 3.d0
        dvdy2 = fyy *  ( dve*dve + dvw*dvw + dve*dvw) / 3.d0
  
        dudxdvdy = fxy * (due*dvn + due*dvs + duw*dvn + duw*dvs)/4.d0
        dvdxdudy = fxy * (dve*dun + dve*dus + dvw*dun + dvw*dus)/4.d0

        eps2(i,j) = 2.d0 * (dudx2 + dvdy2 + dudxdvdy)  +  dvdxdudy   &
                  + (dvdx2 + dudy2) / 2.d0

     else    ! thckmask = 0

        eps2(i,j) = 0.d0

     endif   ! thckmask

  enddo   ! i
  enddo   ! j

  end subroutine calc_eps2_2d

!****************************************************************************

  end module glissade_velo_higher
!*************************************************************************
