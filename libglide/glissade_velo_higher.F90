! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! +                                                           + 
! +  glissade_velo_higher.F90                                 + 
! +                                                           + 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains routines for computing the ice velocity
!  using a variational finite-element approach.
! See Dukowicz, Price and Lipscomb (J. Glaciology, 2010).
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
    use glimmer_physcon, only: gn, rhoi, rhoo, grav
    use glimmer_log, only: write_log
    use glimmer_sparse_type
    use glimmer_sparse
    
    use parallel, only: main_task

!    use glide_deriv
!    use glide_grids, only: stagvarb

    implicit none

    private
    public :: glissade_velo_higher_init, glissade_velo_higher_solve

!TODO - Declare nhalo in parallel modules (in addition to or instead of lhalo and uhalo).    
    public :: nhalo

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

!TODO - At some point we coud make this code more MPAS-friendly by replacing 
!       the ij loops with indirect addressing.  But keep i and j for now.

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

!whl - Can remove gVerticesOnVertex
!      Might want to replace gVerticesOnCell with gNodesOnElement in FE subroutines
  
!whl - remove this one too?
!    integer, dimension(:,:), allocatable ::    &
!       VertexOnCellID        ! VertexID for each vertex of each cell

!    integer, dimension(:,:), allocatable ::    &
!       gVerticesOnVertex,   &! global index of vertices neighboring a given vertex
!       gVerticesOnCell       ! global index of vertices of a given cell

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
       phi,            &    ! interpolation function, evaluated at quad pts
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

!whl - debug
    logical :: verbose_init = .false.   ! for debug print statements
    logical :: verbose_Jac = .false.
    logical :: verbose = .true.   ! for debug print statements

!whl - debug
    integer, parameter :: &
       itest = 24, jtest = 17, ktest = 1

    integer, parameter :: ntest = 2371  ! nodeID for (24,17,1)

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

!whl - debug
    real(dp) :: sumx, sumy, sumz

    if (verbose_init) then
       print*, ' '
       print*, 'In glissade_velo_higher_init'
    endif

    !----------------------------------------------------------------
    ! Assign indices for the cells and vertices
    ! TODO - Decide whether these indices should be local or global for parallel code,
    !        and rewrite as needed.
    !----------------------------------------------------------------

    ! Define global cell index for cells on this processor.
    ! Note: Cell 1 is at lower left, cell nCells is at upper right

!whl - MPAS has indexToCellID and indexToVertexID.  

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

  subroutine glissade_velo_higher_solve(nx,       ny,       &
                                        nlyr,     &     
                                        sigma,    &
                                        nhalo,    &
                                        thck,     usfc,     &
                                        stagthck, stagusfc, &
                                        thklim,             &
                                        flwa,     &
                                        uvel,     vvel,     &
!                                        beta,               &
!                                        whichbabc,          &
                                        whichefvs,          &
                                        whichresid,             &
                                        whichnonlinear,         &
                                        whichsparse)

    !--------------------------------------------------------
    ! Input-output variables
    !--------------------------------------------------------

!TODO - Compute nx, ny, and nlyr based on size(thck) and size(flwa)?
!Note: nlyr = upn-1; change to nz = upn?

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions
       nlyr,                 &  ! number of vertical layers
       nhalo                    ! number of rows/columns of halo cells
 
    real(dp), dimension(:), intent(in) :: &
       sigma

    real (dp), dimension(:,:), intent(in) ::  &
       thck,                 &  ! ice thickness
       usfc                     ! upper surface elevation

    real (dp), dimension(:,:), intent(in) ::  &
       stagusfc,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

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
       nElements,         &   ! no. of elements (nlyr per active cell)
       nNodes                 ! no. of nodes belonging to these elements

    ! Change dimensions to nCells and nVertices?
    logical, dimension(nx,ny) :: &
       active_cell,       &   ! true for active cells (thck > thklim)
       active_vertex          ! true for vertices of active cells

    integer, dimension(nlyr+1,nx,ny) ::  &
        NodeID             ! ID for each active node

    integer, dimension(nVertices*(nlyr+1)) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nVertices*(nlyr+1)) :: &
       xNode, ynode, znode   ! x, y and z coordinates for each node

    integer, dimension(nNodesPerElement,nlyr,nx,ny) ::  &
       NodeOnElementID     ! node ID for each node of each element

    real(dp), dimension(nlyr,nx,ny) :: &
       visc               ! effective viscosity of each element

    real(dp), dimension(:,:,:), allocatable ::   &
       usav, vsav         ! previous guess for velocity solution

    real(dp), dimension(:,:,:,:), allocatable ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 3 (node and its 2 neighbors in z direction) 
                          ! 2nd dimension = 3 (node and its 2 neighbors in x direction) 
                          ! 3rd dimension = 3 (node and its 2 neighbors in y direction) 
                          ! 4th dimension = nNodes

    real(dp), dimension(:), allocatable ::  &
       bu, bv             ! assembled load vector, divided into 2 parts

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
       iter               ! linear iteration count from sparse_easy_solve

!whl - debug
    real(dp), dimension(:,:), allocatable :: Atest, Auu_test, Auv_test, Avu_test, Avv_test
    integer :: i, j, k, n

    !------------------------------------------------------------------------------
    ! Setup for higher-order solver: Compute nodal geometry, allocate storage, etc.
    !------------------------------------------------------------------------------

    if (verbose) then
       print*, ' '
       print*, 'In glissade_velo_higher_solve'
       print*, 'nx, ny, nlyr:', nx, ny, nlyr
       print*, 'sigma:', sigma
       print*, 'thklim:', thklim
       print*, 'max thck:', maxval(thck(:,:))
       print*, 'max usfc:', maxval(usfc(:,:))
       print*, ' '
    endif

    ! Identify and count the active cells (i.e., cells with thck > thklim)

    call get_nodal_geometry(nx,          ny,       &   
                            nlyr,        nhalo,       sigma,    &
                            thck,        thklim,  &
                            stagusfc,    stagthck, &
                            nElements,   nNodes,      active_cell,    &
                            NodeID,      NodeOnElementID,             &
                            iNodeIndex,  jNodeIndex,  kNodeIndex,   &
                            xNode,       yNode,       zNode)

    if (verbose) then
       print*, ' '
       i = itest; j = jtest; k = ktest
       print*, 'itest, jtest, ktest:', i, j, k
       print*, 'NodeID =', NodeID(k,i,j)
       print*, 'NodeOnElementID =', NodeOnElementID(:,k,i,j)
       print*, ' '
    endif

    ! Allocate space for the intermediate (dense) stiffness matrix and rhs vector
    ! (Not needed if using standalone Trilinos solver?)

    allocate(Auu(3,3,3,nNodes), Auv(3,3,3,nNodes))
    allocate(Avu(3,3,3,nNodes), Avv(3,3,3,nNodes))
    allocate(bu(nNodes), bv(nNodes))

    ! Allocate space for the sparse matrix (A), rhs (b), answer (x), and residual vector (Ax-b).

    matrix_order = 2*nNodes         ! Is this exactly enough?
    nNonzero = matrix_order*27      ! 27 = node plus 26 nearest neighbors in hexahedral lattice   

    if (verbose) then
       print*, 'matrix_order =', matrix_order
       print*, 'nNonzero = ', nNonzero
    endif

!whl - debug
!    matrix_order = 3
!    nNonzero = 9
!    matrix_order = 4
!    nNonzero = 16

!    allocate(Atest(matrix_order,matrix_order))
!    allocate(Auu_test(2,2), Auv_test(2,2), Avu_test(2,2), Avv_test(2,2))

    allocate(matrix%row(nNonzero), matrix%col(nNonzero), matrix%val(nNonzero))
    allocate(rhs(matrix_order), answer(matrix_order), resid_vec(matrix_order))

    ! Allocate space for velocity vectors (saved from previous guess)

    !TODO - Use 'size' function to make sure these are the same size as uvel, vvel.

    allocate(usav(nlyr+1,nx-1,ny-1), vsav(nlyr+1,nx-1,ny-1))  

    !---------------------------------------------------------------
    ! Solve Ax = b
    !---------------------------------------------------------------

    ! Much of the following is copied from glam_strs2.F90

    ! set initial values

    resid_velo = 1.d0
    counter = 1
    L2_norm = 1.0d20

    L2_target = 1.0d-4      !TODO - Is this the best value?  
                            ! Should it stay fixed during the outer loop, or should it evolve?

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

    !TODO - Start the outer loop here

    ! initialize outer loop convergence variables (guarantees at least one loop)
    if (counter == 1) then           
       outer_it_criterion = 1.0
       outer_it_target = 0.0
    endif

    do while (outer_it_criterion >= outer_it_target .and. counter < cmax)

       ! choose outer loop stopping criterion
       if (counter > 1) then
          if (whichresid == 3) then
             outer_it_criterion = L2_norm
             outer_it_target = L2_target
          else
             outer_it_criterion = resid_velo
             outer_it_target = resid_target
          end if
       else
          outer_it_criterion = 1.0d10
          outer_it_target = 1.0d-12
       end if

       if (verbose) then
          print*, 'counter =', counter
          print*, 'whichresid =', whichresid
          print*, 'L2_norm, L2_target =', L2_norm, L2_target
          print*, 'resid_velo, resid_target =', resid_velo, resid_target
       endif

       ! save current velocity
       usav(:,:,:) = uvel(:,:,:)
       vsav(:,:,:) = vvel(:,:,:)

       ! TODO: Compute the effective viscosity for each element 
       !       (ideally by summing over quadrature points)

       ! These are the arguments in glam_strs2
!       call compute_effective_viscosity (ewn,  nsn,  upn,      &
!                                         stagsigma,  counter,  &
!                                         whichefvs,  efvs,     &
!                                         uvel,       vvel,     &
!                                         flwa,       thck,     &
!                                         dusrfdew,   dthckdew, &
!                                         dusrfdns,   dthckdns, &
!                                         umask)


       ! Assemble the stiffness matrix
       !TODO - Remove unneeded arguments.

       if (verbose) then
          print*, ' '
          print*, 'call assemble_stiffness_matrix'
       endif

       call assemble_stiffness_matrix(nx, ny, nlyr,  nhalo,              &
                                      active_cell,                       &
                                      nCells,           nActiveCells,    &
                                      NodeID,           NodeOnElementID, &
                                      iNodeIndex,       jNodeIndex,      &
                                      kNodeIndex,                        &
                                      xNode,            yNode,           &
                                      zNode,                             &
                                      nElements,        nNodes,          &
                                      visc,                              &
                                      Auu,              Auv,             &
                                      Avu,              Avv)
       
       if (verbose) print*, 'Assembled the stiffness matrix'

       stop

       ! TODO: Assemble the load vector

       call assemble_load_vector(nElements,        nNodes,      &
                                 bu,               bv)


       ! Put the nonzero matrix elements into the required sparse format
       ! (For SLAP solver; may not be needed for Trilinos)

!WHL - debug - Test sparse_easy_solve

       matrix%order = matrix_order   ! set above   
       matrix%nonzeros = nNonzero
       matrix%symmetric = .false.

       if (matrix%order == 3) then
          Atest(1,1:3) = (/1.d0,  2.d0,  3.d0 /)
          Atest(2,1:3) = (/4.d0,  5.d0,  6.d0 /)
          Atest(3,1:3) = (/7.d0,  8.d0,  9.d0 /)
          rhs(1:3)   =  (/10.d0, 28.d0, 46.d0 /)   ! answer = (3 2 1)
       else if (matrix%order == 4) then
          Atest(1,1:4) = (/1.d0,  2.d0,  3.d0,  4.d0 /)
          Atest(2,1:4) = (/2.d0,  3.d0,  4.d0,  5.d0 /)
          Atest(3,1:4) = (/3.d0,  4.d0,  5.d0,  6.d0 /)
          Atest(4,1:4) = (/4.d0,  5.d0,  6.d0,  7.d0 /)
          rhs(1:4)    = (/20.d0, 30.d0, 40.d0, 50.d0/)   ! answer = (4 3 2 1)
       endif

       n = 0
       do j = 1, matrix%order
          do i = 1, matrix%order
             n = n + 1
             matrix%row(n) = i
             matrix%col(n) = j
             matrix%val(n) = Atest(i,j)
          enddo
       enddo

!whl - debug - skip solver_preprocess for now

       go to 100

          ! for testing Auu, etc.
          Auu_test(1,1:2) = (/ 1.d0, 2.d0 /)
          Auu_test(2,1:2) = (/ 2.d0, 3.d0 /)

          Auv_test(1,1:2) = (/ 3.d0, 4.d0 /)
          Auv_test(2,1:2) = (/ 4.d0, 5.d0 /)

          Avu_test(1,1:2) = (/ 3.d0, 4.d0 /)
          Avu_test(2,1:2) = (/ 4.d0, 5.d0 /)

          Avv_test(1,1:2) = (/ 5.d0, 6.d0 /)
          Avv_test(2,1:2) = (/ 6.d0, 7.d0 /)

       ! Given the intermediate stiffness matrices (Auu, etc.) and load vector (bu, bv),
       ! form the global matrix (in sparse matrix format) and rhs.

       call solver_preprocess(nNodes,       NodeID,        &
                              iNodeIndex,   jNodeIndex,   kNodeIndex,   &
                              Auu,          Auv,   &
                              Avu,          Avv,   &
                              bu,           bv,    &
                              uvel,         vvel,  &
                              matrix_order, nNonzero,  &
                              matrix,       rhs,   &
                              answer)
                                       

100    print*, 'Calling sparse_easy_solve, whichsparse, order =', whichsparse, matrix%order

       ! Solve the linear matrix problem

       call sparse_easy_solve(matrix, rhs, answer, err, iter, whichsparse)

       print*, 'answer =', answer
       print*, 'err =', err       
       print*, 'iter =', iter

       stop

!whl - end debug

       !TODO - glam_strs2 calls res_vect here--why?
       !!!call res_vect( matrix, vk_1, rhs, size(rhs), g_flag, L2square, whichsparse )
       !! L2norm  = L2square
       !! F(1:pcgsize(1)) = vk_1(:)

       if (whichsparse /= STANDALONE_TRILINOS_SOLVER) then
          call sparse_easy_solve(matrix, rhs, answer, err, iter, whichsparse)
#ifdef TRILINOS
       else
          call solvewithtrilinos(rhs, answer, linearSolveTime)
          totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
          ! write(*,*) 'Total linear solve time so far', totalLinearSolveTime                                           
#endif
       endif

       ! Put the velocity solution back into 3D arrays

       call solver_postprocess(nNodes,       answer,                   &
                               iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                               uvel,         vvel)

       ! Any ghost postprocessing needed, as in glam_strs2?

       ! Halo updates for velocity (to be tested)
       !TODO - Are all four of these calls needed?
!       call staggered_parallel_halo(uvel)
!       call staggered_parallel_halo(vvel)
!       call horiz_bcs_stag_vector_ew(uvel)
!       call horiz_bcs_stag_vector_ns(vvel)


       ! Compute the residual vector and its L2 norm

       call compute_residual_vector(matrix, answer, rhs,  &
                                    resid_vec, L2_norm)


       ! Compute residual quantities based on the velocity solution

       call compute_residual_velocity(nhalo,  whichresid,   &
                                      uvel,   vvel,        &
                                      usav,   vsav,        &
                                      resid_velo)

       ! Write diagnostics (iteration number, max residual, and location of max residual
       ! (send output to the screen or to the log file, per whichever line is commented out) 

       !TODO - Find out if this note from glam_strs2 still applies:
       ! "Can't use main_task flag because main_task is true for all processors in case of parallel_single"
       if (main_task) then
          if (whichresid == 3 )then
             print '(i4,2g20.6)', counter, L2_norm, L2_target
             !write(message,'(i4,3g20.6)') counter, L2_norm, L2_target
             !call write_log (message)
          else
             print '(i4,2g20.6)', counter, resid_velo, resid_target
             !write(message,'(i4,2g20.6)') counter, resid_velo, resid_target
             !call write_log (message)
          end if
       endif

       ! Advance the iteration counter

       counter = counter + 1   

    enddo  ! while (outer_it_criterion >= outer_it_target .and. counter < cmax)

    converged_soln = .true.

    ! Clean up

    deallocate(Auu, Auv, Avu, Avv)
    deallocate(bu, bv)
    deallocate(matrix%row, matrix%col, matrix%val)
    deallocate(rhs, answer, resid_vec)
    deallocate(usav, vsav)

  end subroutine glissade_velo_higher_solve

!****************************************************************************

  subroutine compute_residual_vector(matrix, answer, rhs, &
                                     resid_vec, L2_norm)

    ! Compute the residual vector Ax - b and its L2 norm.

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

    resid_vec(:) = 0.d0

    do n = 1, matrix%nonzeros
       i = matrix%row(n)
       j = matrix%col(n)
       resid_vec(i) = resid_vec(i) + matrix%val(n)*answer(j)
    enddo

    L2_norm = 0.d0
    do n = 1, matrix%order
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


  subroutine solver_preprocess(nNodes,       NodeID,        &
                               iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                               Auu,          Auv,   &
                               Avu,          Avv,   &
                               bu,           bv,    &
                               uvel,         vvel,  &
                               matrix_order, nNonzero,  &
                               matrix,       rhs,   &
                               answer)

    ! Using the intermediate matrices (Auu, Auv, Avu, Avv), load vectors (bu, bv),
    ! and velocity components (uvel, vvel), form the matrix and the rhs and answer
    ! vectors in the desired sparse matrix format.

    !---------------------------------------------------------
    ! Input-output variables
    !---------------------------------------------------------

    integer, intent(in) :: nNodes        ! number of active nodes

    integer, dimension(:,:,:), intent(in) ::  &
        NodeID             ! ID for each active node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(3,3,3,nNodes), intent(in) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 3 (node and its 2 neighbors in z direction) 
                          ! 2nd dimension = 3 (node and its 2 neighbors in x direction) 
                          ! 3rd dimension = 3 (node and its 2 neighbors in y direction) 
                          ! 4th dimension = nNodes

    real(dp), dimension(nNodes), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(:,:,:), intent(in) ::   &
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
    
    ! Set basic matrix parameters

    matrix%order = matrix_order
    matrix%nonzeros = nNonzero
    matrix%symmetric = .false.

    ! Set the nonzero coefficients of the sparse matrix 

    ct = 0

    do rowA = 1, nNodes

       i = iNodeIndex(rowA)
       j = iNodeIndex(rowA)
       k = iNodeIndex(rowA)
       
       do jA = -1, 1
        do iA = -1, 1
         do kA = -1, 1
            colA = NodeID(k+kA, i+iA, j+jA)

            !TODO - Could shorten the following by calling a putval subroutine 4 times

            ! Auu
            val = Auu(kA,iA,jA,rowA)
            if (val /= 0.d0) then   ! Is this adequate?
               ct = ct + 1
               matrix%row(ct) = 2*rowA - 1
               matrix%col(ct) = 2*colA - 1
               matrix%val(ct) = val
            endif

            ! Auv 
            val = Auv(kA,iA,jA,rowA)
            if (val /= 0.d0) then   ! Is this adequate?
               ct = ct + 1
               matrix%row(ct) = 2*rowA - 1
               matrix%col(ct) = 2*colA
               matrix%val(ct) = val
            endif

            ! Avu 
            val = Avu(kA,iA,jA,rowA)
            if (val /= 0.d0) then   ! Is this adequate?
               ct = ct + 1
               matrix%row(ct) = 2*rowA
               matrix%col(ct) = 2*colA - 1
               matrix%val(ct) = val
            endif

            ! Avv 
            val = Avv(kA,iA,jA,rowA)
            if (val /= 0.d0) then
               ct = ct + 1
               matrix%row(ct) = 2*rowA
               matrix%col(ct) = 2*colA
               matrix%val(ct) = val
            endif

         enddo    ! kA
        enddo     ! iA
       enddo      ! jA

    enddo         ! rowA

    ! Set the current guess for the answer vector
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodes
       i = iNodeIndex(n)
       j = iNodeIndex(n)
       k = iNodeIndex(n)

       answer(2*n-1) = uvel(k,i,j)
       answer(2*n)   = vvel(k,i,j)
    enddo

    ! Set the rhs vector (TODO: Is this adequate?) 
    ! For efficiency, put the u and v terms for a given node adjacent in storage.

    do n = 1, nNodes
       rhs(2*n-1) = bu(n)
       rhs(2*n)   = bv(n)
    enddo

  end subroutine solver_preprocess

!****************************************************************************

  subroutine solver_postprocess(nNodes,       answer,                   &
                                iNodeIndex,   jNodeIndex,  kNodeIndex,  &
                                uvel,         vvel)

  ! Extract the velocities from the solution vector.
                                            
!    use parallel

    !---------------------------------------------------------
    ! Input-output variables
    !---------------------------------------------------------

    integer, intent(in) :: nNodes        ! number of active nodes

    real(dp), dimension(:), intent(in) ::   &
       answer             ! velocity solution vector

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(:,:,:), intent(inout) ::   &
       uvel, vvel         ! u and v components of velocity

    integer :: i, j, k, n

    do n = 1, nNodes
       i = iNodeIndex(n)
       j = iNodeIndex(n)
       k = iNodeIndex(n)

       uvel(k,i,j) = answer(2*n-1)
       vvel(k,i,j) = answer(2*n)

    enddo

  end subroutine solver_postprocess

!****************************************************************************

  subroutine get_nodal_geometry(nx,          ny,       &   
                                nlyr,        nhalo,       sigma,    &
                                thck,        thklim,   &
                                stagusfc,    stagthck, &
                                nElements,   nNodes,      active_cell,  &
                                NodeID,      NodeOnElementID,           &
                                iNodeIndex,  jNodeIndex,  kNodeIndex,   &
                                xNode,       yNode,       zNode)
                            
    ! Identify and count the active nodes for the calculations.
    ! All nodes of all elements lying within active cells (i.e., cells with ice present)
    !  are active.
    ! Active nodes may be either free (to be solved for) or constrained (e.g., Dirichlet).

    !---------------------------------------------------------
    ! Input-output variables
    !---------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nlyr,                    &    ! number of vertical layers
       nhalo                    

    real(dp), dimension(nlyr+1), intent(in) :: &
       sigma                  ! sigma vertical coordinate

    real (dp), dimension(nx,ny), intent(in) ::  &
       thck                   ! ice thickness

    real(dp), intent(in) ::   & 
       thklim                 ! minimum ice thickness for active cells

    real (dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusfc,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    integer, intent(out)  ::            &
       nElements,         &   ! no. of elements (nlyr per active cell)
       nNodes                 ! no. of nodes belonging to these elements

    logical, dimension(nx,ny), intent(out) :: &
       active_cell            ! true for active cells (thck > thklim)

    integer, dimension(nlyr+1,nx,ny), intent(out) ::  &
       NodeID             ! ID for each active node

    integer, dimension(nNodesPerElement,nlyr,nx,ny), intent(out) ::  &
       NodeOnElementID     ! node ID for each node of each element

    integer, dimension(nVertices*nlyr), intent(out) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nVertices*nlyr), intent(out) :: &
       xNode, ynode, znode   ! x, y and z coordinates for each node


    !---------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------

    integer ::   &
       nActiveCells,   &      ! number of active cells (thck > thklim)
       nActiveVertices        ! number of active vertices

    logical, dimension(nx,ny) :: &
       active_vertex          ! true for vertices of active cells

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

    ! Count the active elements (nlyr element per active cell)

    nElements = nActiveCells * nlyr

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
          do k = 1, nlyr+1    ! all nodes in column are active
             nNodes = nNodes + 1   
             NodeID(k,i,j) = nNodes   ! unique index for each active node
             iNodeIndex(nNodes) = i
             jNodeIndex(nNodes) = j
             kNodeIndex(nNodes) = k
             xNode(nNodes) = xVertex(i,j)
             yNode(nNodes) = yVertex(i,j)
             zNode(nNodes) = stagusfc(i,j) - sigma(k)*stagthck(i,j)

             if (verbose .and. nNodes==ntest) then
                print*, ' '
                print*, 'i, j, k, n:', i, j, k, nNodes
                print*, 'sigma, stagusfc, stagthck:', sigma(k), stagusfc(i,j), stagthck(i,j)
                print*, 'x, y, z:', xNode(nNodes), yNode(nNodes), zNode(nNodes)
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
          do k = 1, nlyr
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

  subroutine assemble_stiffness_matrix(nx, ny, nlyr,  nhalo,              &
                                       active_cell,                       &
                                       nCells,           nActiveCells,    &
                                       NodeID,           NodeOnElementID, &
                                       iNodeIndex,       jNodeIndex,     kNodeIndex,   &     
                                       xNode,            yNode,          zNode,        &
                                       nElements,        nNodes,       &
                                       visc,                           &
                                       Auu,              Auv,          &
                                       Avu,              Avv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nlyr,                    &    ! number of vertical layers
       nhalo,                   &    ! number of halo layers
       nCells,                  &    ! number of horizontal cells
       nActiveCells                  ! number of cells with ice present

    logical, dimension(nx,ny), intent(in) :: active_cell

    integer, dimension(nVertices*nlyr), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

    real(dp), dimension(nVertices*nlyr), intent(in) ::   &
       xNode, yNode, zNode  ! x, y and z coordinates of nodes

    integer, dimension(:,:,:), intent(in) ::  &
        NodeID             ! ID for each active node

    integer, dimension(nNodesPerElement,nlyr, nx, ny), intent(in) :: &
        NodeOnElementID    ! ID for each node of each element

!TODO - Are these needed?
    integer, intent(in) ::   &
       nElements,    &    ! number of elements
       nNodes             ! number of nodes associated with these elements

    real(dp), dimension(nElements), intent(in) ::  &
       visc               ! effective viscosity

    real(dp), dimension(3,3) ::   &
       Jac, Jinv          ! Jacobian matrix and its inverse

    real(dp) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement) ::   &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis function, evaluated at quad pts

    !----------------------------------------------------------------
    ! Note: Kuu, Kuv, Kvu, and Kvv are 8x8 stiffness matrices for the local element.
    !
    ! Once these matrices are formed, their coefficients are summed into the assembled
    !  matrices Auu, Auv, Avu, Avv.  The A matrices each have as many row as there are
    !  active nodes, but only 27 columns, corresponding to the 27 vertices that belong to
    !  the 8 elements sharing a given node.
    !
    ! At the next step, the terms of the dense A matrices will be put in a sparse matrix
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

    real(dp), dimension(3,3,3,nNodes), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nNodesPerElement) :: x, y, z   ! nodal coordinates

    integer, dimension(nNodesPerElement) :: nid        ! ID for each node of an element

    integer :: i, j, k, n, p


    if (verbose) then
       print*, ' '
       print*, 'In assemble_stiffness_matrix'
    endif

    ! Initialize arrays

    Auu(:,:,:,:) = 0.d0
    Auv(:,:,:,:) = 0.d0
    Avu(:,:,:,:) = 0.d0
    Avv(:,:,:,:) = 0.d0

    ! Sum over elements

!TODO - Could sum over nCells and do not use i and j?

    do j = nhalo+1, ny-nhalo    ! loop over the cells owned by this processor
    do i = nhalo+1, nx-nhalo
       
     if (active_cell(i,j)) then

       do k = 1, nlyr    ! loop over elements in this column
                         ! assume k increases from upper surface to bed

          ! Find node ID and (x,y,z) coordinates of each node for this element

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, 'i, j, k:', i, j, k
             print*, ' '
          endif

          ! initialize element stiffness matrices
          Kuu(:,:) = 0.d0
          Kuv(:,:) = 0.d0
          Kvu(:,:) = 0.d0
          Kvv(:,:) = 0.d0
    
          ! get nodeID and spatial coordinates for each node
          do n = 1, nNodesPerElement
             nid(n) = NodeOnElementID(n,k,i,j)
             x(n) = xNode(nid(n))
             y(n) = yNode(nid(n))
             z(n) = zNode(nid(n))

             if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
                print*, 'n, nid, x, y, z:', n, nid(n), x(n), y(n), z(n)
             endif

          enddo   ! nodes per element

          ! Loop over quadrature points for this element
   
          do p = 1, nQuadPoints

          ! Evaluate the derivatives of the element basis functions
          ! at this quadrature point.

!whl - debug - Pass in i, j, k, and p for now

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

          ! Increment the element stiffness matrix with the contribution from
          ! this quadrature point.

!whl - debug - Pass in i, j, k, and p for now
             call increment_element_matrix(nNodesPerElement,                        &
                                           wqp(p),         detJ,                    &
                                           dphi_dx(:),     dphi_dy(:),  dphi_dz(:), &
                                           Kuu(:,:),       Kuv(:,:),                &
                                           Kvu(:,:),       Kvv(:,:),    i, j, k, p )

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

!whl - Pass in i, j, k for now
          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nid,                         &
                                        Kuu,              Auu, i, j, k)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kuv into Auv'
          endif

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nid,                         &
                                        Kuv,              Auv, i, j, k)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvu into Avu'
          endif

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nid,                         &
                                        Kvu,              Avu, i, j, k)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             print*, ' '
             print*, 'Insert Kvv into Avv'
          endif

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nid,                         &
                                        Kvv,              Avv, i ,j, k)

       ! TODO: Finish assembling the stiffness matrix by incorporating basal sliding?

       enddo   ! nlyr (loop over elements in this column)

     endif   ! active cell

    enddo      ! i
    enddo      ! j

    ! Check symmetry of assembled matrix
    ! This call could be skipped in a well-tested production code.

    call check_symmetry_assembled_matrix(nNodes,      NodeID,      &     
                                         iNodeIndex,  jNodeIndex,  &
                                         kNodeIndex,               &
                                         Auu,         Auv,    &
                                         Avu,         Avv)


  end subroutine assemble_stiffness_matrix

!****************************************************************************

!whl - debug - pass in i, j, k, p for now
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

!whl - debug
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
       !whl - do a proper abort here
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

!whl - debug - Pass in i, j, k, and p for now

  subroutine increment_element_matrix(nNodesPerElement,             &
                                      wqp,      detJ,               &
                                      dphi_dx,  dphi_dy,  dphi_dz,  &
                                      Kuu,      Kuv,                &
                                      Kvu,      Kvv,      ii, jj, k, p)

!TODO - Add optional arguments that reduce the solver to SIA or SSA.

    !------------------------------------------------------------------
    ! Increment the stiffness matrices Kuu, Kuv, Kvu, Kvv with the
    ! contribution from a particular quadrature point, 
    ! based on the Blatter-Pattyn first-order equations.
    !
    ! This subroutine will work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

!whl - debug - Pass in i, j, k, and p for now
    integer, intent(in) :: ii, jj, k, p

    integer, intent(in) :: nNodesPerElement  ! number of nodes per element

    real(dp), intent(in) ::    &
             wqp,        &! weight for this quadrature point
             detJ         ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
             dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis functions,
                                         ! evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(inout) :: &
             Kuu, Kuv, Kvu, Kvv     ! components of element stiffness matrix

    integer :: i, j

    ! Increment the element stiffness matrices for the first-order Blatter-Pattyn equations

    do j = 1, nNodesPerElement      ! columns of K
       do i = 1, nNodesPerElement   ! rows of K

          Kuu(i,j) = Kuu(i,j) + wqp*detJ*(4.d0*dphi_dx(j)*dphi_dx(i)  &
                                         +     dphi_dy(j)*dphi_dy(i)  &
                                         +     dphi_dz(j)*dphi_dz(i) )
          Kuv(i,j) = Kuv(i,j) + wqp*detJ*(2.d0*dphi_dx(j)*dphi_dy(i)  &
                                         +     dphi_dy(j)*dphi_dx(i) )
          Kvu(i,j) = Kvu(i,j) + wqp*detJ*(2.d0*dphi_dy(j)*dphi_dx(i)  &
                                         +     dphi_dx(j)*dphi_dy(i) )
          Kvv(i,j) = Kvv(i,j) + wqp*detJ*(4.d0*dphi_dy(j)*dphi_dy(i)  &
                                         +     dphi_dx(j)*dphi_dx(i)  &
                                         +     dphi_dz(j)*dphi_dz(i) )

       enddo  ! i (rows)
    enddo     ! j (columns)

    if (verbose .and. ii==itest .and. jj==jtest .and. k==ktest) then
       print*, 'Increment element matrix, p =', p
    endif

  end subroutine increment_element_matrix

!****************************************************************************

!whl - Pass in i, j, k for now
  subroutine element_to_global_matrix(nNodesPerElement, nNodes,   &
                                      nid,                        &
                                      Kmat,             Amat, i, j, k)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! So this subroutine is called four times per element.

!whl - Pass in i, j, k for now
    integer, intent(in) :: i, j, k

    integer, intent(in) ::   &
       nNodesPerElement,     &
       nNodes            ! total number of active nodes
       
    integer, dimension(nNodesPerElement), intent(in) ::  &
       nid               ! node ID for each node of this element

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(3,3,3,nNodes), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: m, n, iA, jA, kA

    if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
       print*, ' '
       print*, 'First row of K:'
       write(6, '(8e12.4)') Kmat(1,:)
    endif


    do m = 1, nNodesPerElement       ! rows of K

       do n = 1, nNodesPerElement    ! columns of K

          kA = kshift(m,n)           ! k index of A into which K(m,n) is summed
          iA = ishift(m,n)           ! similarly for i and j indices
          jA = jshift(m,n)         

          !TODO - Add viscosity term
          Amat(kA,iA,jA,nid(m)) = Amat(kA,iA,jA,nid(m)) + Kmat(m,n)

          if (verbose .and. i==itest .and. j==jtest .and. k==ktest) then
             if (m==7) then
!                print*, ' '
!                print*, 'row, col, K(m,n):', m, n, Kmat(m,n)
!                print*, 'nodeID for this row =', nid(m)
!                print*, 'kA, iA, jA', kA, iA, jA
!                print*, 'New value of A =', Amat(kA,iA,jA,nid(m))
             endif
          endif
          
          if (verbose .and. nid(m)==ntest) then
!!              print*, ' '
!!              print*, 'Contribution to nid', nid(m), 'from element i, j, k =', i, j, k
!!              print*, 'm, n, Kmat(m,n) =', m, n, Kmat(m,n)
!!              print*, 'kA, iA, jA', kA, iA, jA
!!              print*, 'New value of A =', Amat(kA,iA,jA,nid(m))
          endif

       enddo     ! n
    enddo        ! m

  end subroutine element_to_global_matrix

!****************************************************************************

!TODO - Make this 3D

  subroutine assemble_load_vector(nElements,        nNodes,      &
                                  bu,               bv)

    integer, intent(in) ::   &
       nElements,    &    ! number of elements on this processor
       nNodes             ! number of nodes associated with these elements

    real(dp), dimension(nNodes) ::  &
       bu, bv             ! assembled load vector, divided into 2 parts

    bu(:) = 0.d0
    bv(:) = 0.d0

!whl - Write this subroutine

!      First the gravitational term

!      Then the pressure term from immersed shelves


  end subroutine assemble_load_vector

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

  subroutine check_symmetry_assembled_matrix(nNodes,     NodeID,    &
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

    integer, intent(in) :: nNodes     ! number of active nodes = number of rows of Auu, etc.

    integer, dimension(:,:,:), intent(in) ::  &
       NodeID             ! ID for each active node

    integer, dimension(:), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

    real(dp), dimension(3,3,3,nNodes), intent(in) ::   &
             Auu, Auv, Avu, Avv     ! component of assembled stiffness matrix
                                    !
                                    !    Auu  | Auv
                                    !    _____|____          
                                    !    Avu  | Avv
                                    !         |

    integer :: i, j, k, iA, jA, kA, row, col

    real(dp) :: val1, val2          ! values of matrix coefficients

    ! check that Auu = (Auu)^T

!whl - debug
    real(dp) :: maxdiff

    maxdiff = 0.d0

    do row = 1, nNodes      ! rows of Auu
       i = iNodeIndex(row)  ! i, j and k values of this node
       j = jNodeIndex(row)
       k = kNodeIndex(row)
!!       print*, ' '
!!       print*, 'row, i, j, k:', row, i, j, k
       do jA = -1, 1
          do iA = -1, 1
             do kA = -1, 1
                col = NodeID(k+kA, i+iA, j+jA)   ! corresponding column of full matrix 
                val1 = Auu(kA,iA,jA,row)         ! value of Auu(row,col)
                val2 = Auu(-kA, -iA, -jA, col)   ! value of Auu(col,row)
!!                print*, ' '
!!                print*, 'iA, jA, kA, col:', iA, jA, kA, col
!!                print*, 'val1, val2:', val1, val2

                if (abs(val2 - val1) > eps12) then
                   print*, 'Auu is not symmetric'
                   print*, 'row, col, Auu(row,col), Auu(col,row):', &
                            row, col, val1, val2
                   stop
                endif

!whl - debug
                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             enddo  ! kA
          enddo     ! iA
       enddo        ! jA
    enddo           ! row of Auu

!whl - debug
    if (verbose) then
       print*, ' '
       print*, 'Auu: max difference from symmetry =', maxdiff
    endif

    ! check that Avv = (Avv)^T
    !TODO - The following code is basically identical to the above code for Auu.
    !       Combine into one code chunk called twice?

    maxdiff = 0.d0

    do row = 1, nNodes      ! rows of Avv
       i = iNodeIndex(row)  ! i, j and k values of this node
       j = jNodeIndex(row)
       k = kNodeIndex(row)
       do jA = -1, 1
          do iA = -1, 1
             do kA = -1, 1
                col = NodeID(k+kA, i+iA, j+jA)   ! corresponding column of full matrix 
                val1 = Avv(kA,iA,jA,row)         ! value of Avv(row,col)
                val2 = Avv(-kA, -iA, -jA, col)   ! value of Avv(col,row)
                if (abs(val2 - val1) > eps12) then
                   print*, 'Avv is not symmetric'
                   print*, 'row, col, Avv(row,col), Avv(col,row):', &
                            row, col, val1, val2
                   stop
                endif

!whl - debug
                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             enddo  ! kA
          enddo     ! iA
       enddo        ! jA
    enddo           ! row of Avv

!whl - debug
    if (verbose) print*, 'Avv: max difference from symmetry =', maxdiff

    ! Check that Auv = (Avu)^T

    maxdiff = 0.d0

    do row = 1, nNodes      ! rows of Auv
       i = iNodeIndex(row)  ! i, j and k values of this node
       j = jNodeIndex(row)
       k = kNodeIndex(row)
       do jA = -1, 1
          do iA = -1, 1
             do kA = -1, 1
                col = NodeID(k+kA, i+iA, j+jA)   ! corresponding column of full matrix 
                val1 = Auv(kA,iA,jA,row)         ! value of Auv(row,col)
                val2 = Avu(-kA, -iA, -jA, col)   ! value of Avu(col,row)
                if (abs(val2 - val1) > eps12) then
                   print*, 'Auv is not equal to (Avu)^T'
                   print*, 'row, col, Auv(row,col), Avu(col,row):', &
                            row, col, val1, val2
                   stop
                endif

!whl - debug
                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             enddo  ! kA
          enddo     ! iA
       enddo        ! jA
    enddo           ! row of Auv

!whl - debug
    if (verbose) print*, 'Auv v. (Avu)^T: max difference from symmetry =', maxdiff

    ! Check that Avu = (Auv)^T
    ! This check may not be necessary given that we've already checked that Auv = (Avu)^T.
    !  Just being extra careful.

    maxdiff = 0.d0

    do row = 1, nNodes      ! rows of Avu
       i = iNodeIndex(row)  ! i, j and k values of this node
       j = jNodeIndex(row)
       k = kNodeIndex(row)
       do jA = -1, 1
          do iA = -1, 1
             do kA = -1, 1
                col = NodeID(k+kA, i+iA, j+jA)   ! corresponding column of full matrix 
                val1 = Avu(kA,iA,jA,row)         ! value of Avu(row,col)
                val2 = Auv(-kA, -iA, -jA, col)   ! value of Auv(col,row)
                if (abs(val2 - val1) > eps12) then
                   print*, 'Auv is not equal to (Avu)^T'
                   print*, 'row, col, Auv(row,col), Avu(col,row):', &
                            row, col, val1, val2
                   stop
                endif

!whl - debug
                if (abs(val2 - val1) > maxdiff) maxdiff = abs(val2 - val1)
                
             enddo  ! kA
          enddo     ! iA
       enddo        ! jA
    enddo           ! row of Avu

!whl - debug
    if (verbose) print*, 'Avu v. (Auv)^T: max difference from symmetry =', maxdiff

    if (verbose) print*, 'The assembled matrix is symmetric.'

  end subroutine check_symmetry_assembled_matrix

!****************************************************************************

!****************************************************************************

!****************************************************************************

!****************************************************************************

!whl - Rename and revise this subroutine
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

!whl - may be able to speed this up by setting duw(i,j) = due(i-1,j)?

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
