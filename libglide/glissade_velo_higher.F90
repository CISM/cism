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
    
!    use glide_deriv
!    use glide_grids, only: stagvarb

    implicit none
    private

    !----------------------------------------------------------------
    ! Here are some definitions:
    !
    ! The horizontal mesh is composed of cells and vertices.
    ! For now, all cells are assumed to be quadrilaterals, but the code can be
    !  generalized later to triangles (e.g., for MPAS mesh).
    ! Each cell can be extruded to form a column with a specified number of layers.
    ! 
    ! An element is a layer of a cell, and a node is a corner of an element.
    ! So elements and nodes live in 3D space, whereas cells and vertices live in
    !  the horizontal plane.
    !
    ! Active elements are those elements lying in cells with ice present
    !  (that is, the ice thickness is above a minimum threshold).
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
    integer, dimension(:,:), allocatable ::    &
       VertexOnCellID        ! VertexID for each vertex of each cell

!    integer, dimension(:,:), allocatable ::    &
!       gVerticesOnVertex,   &! global index of vertices neighboring a given vertex
!       gVerticesOnCell       ! global index of vertices of a given cell

    !----------------------------------------------------------------
    ! Finite element properties
    !
    ! For now, assume 3D hexahedral elements.
    !----------------------------------------------------------------

    integer, parameter ::      &
       nNodesPerElement = 8,   &   ! 8 nodes for hexahedra
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

! storage arrays

    public :: glissade_velo_higher_init, glissade_velo_higher_solve

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
    allocate(xvertex(nx-1,ny-1), yvertex(nx-1,ny-1))

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

!TODO - VertexOnCellID may not be needed (use NodeOnElement instead)
    ! Determine the vertex IDs for each vertex of each cell.
    ! Numbering is counter-clockwise from SW corner
    ! Requires nhalo >= 1

    allocate(VertexOnCellID(nCells, nVerticesOnCell))  ! TODO - Should ncells be on the outside?
    VertexOnCellID(:,:) = 0
    do j = 1+nhalo, ny-nhalo   ! locally owned cells
    do i = 1+nhalo, nx-nhalo
       nc = CellID(i,j)     
       VertexOnCellID(nc,1) = VertexID(i-1,j-1)   ! SW corner
       VertexOnCellID(nc,2) = VertexID(i,j-1)     ! SE corner
       VertexOnCellID(nc,3) = VertexID(i,j)       ! NE corner
       VertexOnCellID(nc,4) = VertexID(i-1,j)     ! NW corner
    enddo
    enddo
    
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

    ! Evaluate basis functions and their derivatives at each quad pt
    ! TODO - Check these carefully.

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

       ! Initialize some matrices that describe how the i, j and k indices of each node
       ! in each element are related to one another.

       ! The kshift matrix describes how the k indices of the 8 nodes are related to one another.
       ! E.g, if kshift (1,5) = -1, this means that the k index of node 5 has a k index
       ! one less than the k index of node 1.  (Assume that k increases downward.)

       kshift(1, 1:8) = (/ 0,  0,  0,  0, -1, -1, -1, -1/)   
       kshift(2, 1:8) = kshift(1, 1:8)
       kshift(3, 1:8) = kshift(1, 1:8)
       kshift(4, 1:8) = kshift(1, 1:8)
       kshift(5, 1:8) = (/ 1,  1,  1,  1,  0,  0,  0,  0/)
       kshift(6, 1:8) = kshift(5, 1:8)
       kshift(7, 1:8) = kshift(5, 1:8)
       kshift(8, 1:8) = kshift(5, 1:8)

       ! The ishift matrix describes how the i indices of the 8 nodes are related to one another.

       ishift(1, 1:8) = (/ 0,  1,  1,  0,  0,  1,  1,  0/)   
       ishift(2, 1:8) = (/-1,  0,  0, -1, -1,  0,  0, -1/)   
       ishift(3, 1:8) = ishift(2,1:8)
       ishift(4, 1:8) = ishift(1,1:8)
       ishift(5, 1:8) = ishift(1,1:8)
       ishift(6, 1:8) = ishift(2,1:8)
       ishift(7, 1:8) = ishift(2,1:8)
       ishift(8, 1:8) = ishift(1,1:8)

       ! And the ishift matrix describes how the j indices of the 8 nodes are related to one another.

       jshift(1, 1:8) = (/ 0,  0,  1,  1,  0,  0,  1,  1/)   
       jshift(2, 1:8) = jshift(1,1:8)
       jshift(3, 1:8) = (/-1, -1,  0,  0, -1, -1,  0, 0/)   
       jshift(4, 1:8) = jshift(3,1:8)
       jshift(5, 1:8) = jshift(1,1:8)
       jshift(6, 1:8) = jshift(1,1:8)
       jshift(7, 1:8) = jshift(3,1:8)
       jshift(8, 1:8) = jshift(3,1:8)

    enddo

!TODO - Write out and debug some of the variables computed above

   go to 100

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
         write(6,*) 'i, j, vertexID:', i, j, VertexID(i,j)
         write(6,*) 'x, y coords:', xVertex(i,j), yVertex(i,j)
      enddo
   enddo

100 return

  end subroutine glissade_velo_higher_init

!****************************************************************************

  subroutine glissade_velo_higher_solve(nx,       ny,       &
                                        nlyr,     sigma,    &
                                        nhalo,    &
                                        thck,     flwa,     &
                                        stagthck, stagusfc, &
                                        thckmin,            &
                                        uvel,     vvel)

    integer, intent(in) ::   &
       nx, ny,             &  ! horizontal grid dimensions
       nlyr,               &  ! number of vertical layers
       nhalo                  ! number of rows/columns of halo cells
 
    real(dp), dimension(nlyr+1), intent(in) :: &
       sigma

    real (dp), dimension(nx,ny), intent(in) ::  &
       thck                   ! ice thickness

    real (dp), dimension(nlyr,nx,ny), intent(in) ::  &
       flwa                   ! flow factor

    real(dp), intent(in) ::   & 
       thckmin                ! minimum ice thickness for active cells

    real (dp), dimension(nx-1,ny-1), intent(in) ::  &
       stagusfc,            & ! upper surface averaged to vertices
       stagthck               ! ice thickness averaged to vertices

    real (dp), dimension(nx-1,ny-1), intent(inout) ::  &
       uvel, vvel             ! velocity components

    integer ::            &
       nActiveCells,      &   ! no. of active cells (thck > thckmin)
       nElements,         &   ! no. of elements (nlyr per active cell)
       nNodes                 ! no. of nodes belonging to these elements

    ! Change dimensions to nCells and nVertices?
    logical, dimension(nx,ny) :: &
       active_cell,       &   ! true for active cells (thck > thckmin)
       active_vertex          ! true for vertices of active cells

!TODO - Not sure if iCellIndex and jCellIndex are needed
    integer, dimension(nCells) ::   &
       iCellIndex, jCellIndex ! i and j indices of active cells

    integer, dimension(nlyr,nx,ny) ::  &
        NodeID             ! ID for each active node

    integer, dimension(nVertices*nlyr) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of active nodes

    real(dp), dimension(nVertices*nlyr) :: &
       xNode, ynode, znode   ! x, y and z coordinates for each node

!TODO - may not be needed
!    logical, dimension(nlyr+1,nx,ny) ::  &
!       active_node            ! true for active nodes

    integer, dimension(nNodesPerElement,nlyr,nx,ny) ::  &
        NodeOnElement     ! node ID for each node of each element

    real(dp), dimension(nlyr,nx,ny) :: &
       visc               ! effective viscosity of each element

    real(dp), dimension(:,:,:,:), allocatable ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 3 (node and its 2 neighbors in z direction) 
                          ! 2nd dimension = 3 (node and its 2 neighbors in x direction) 
                          ! 3rd dimension = 3 (node and its 2 neighbors in y direction) 
                          ! 4th dimension = nNodes

    real(dp), dimension(:), allocatable ::  &
       bu, bv             ! assembled load vector, divided into 2 parts

    integer :: i, j, k, gc, gv, nv, nid

!TODO - Are these needed?
    integer, dimension(nCells) :: gElementIndex   ! global cell index for each element
    integer, dimension(nVertices) :: gNodeIndex   ! global node index for each vertex

!TODO - Put the next few loops in a separate geometry subroutine?

    ! Identify and count the active cells (i.e., cells with thck > thckmin)

    nActiveCells = 0
    active_cell(:,:) = .false.
    do j = 1+nhalo, ny-nhalo
    do i = 1+nhalo, nx-nhalo
       if (thck(i,j) >= thckmin) then
          active_cell(i,j) = .true.
          nActiveCells = nActiveCells + 1
!!          iCellIndex(nActiveCells) = i  ! could use these for indirect addressing
!!          jCellIndex(nActiveCells) = j 
       endif
    enddo
    enddo
        
    ! Count the active elements (nlyr element per active cell)

    nElements = nActiveCells * nlyr

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

!    active_node(:,:,:) = .false.
    nNodes = 0
    iNodeIndex(:) = 0
    jNodeIndex(:) = 0
    kNodeIndex(:) = 0    
    xNode(:) = 0.d0
    yNode(:) = 0.d0
    zNode(:) = 0.d0

    do j = nhalo, ny-nhalo    ! include S edge
    do i = nhalo, nx-nhalo    ! include W edge
       if (active_vertex(i,j)) then
          do k = 1, nlyr+1    ! all nodes in column are active
!             active_node(k,i,j) = .true.
             nNodes = nNodes + 1   
             NodeID(k,i,j) = nNodes   ! unique index for each active node
             iNodeIndex(nNodes) = i
             jNodeIndex(nNodes) = j
             kNodeIndex(nNodes) = k
             xNode(nNodes) = xVertex(i,j)
             yNode(nNodes) = yVertex(i,j)
             zNode(nNodes) = stagusfc(i,j) - sigma(k)*stagthck(i,j)
           enddo   ! k
        endif      ! active cell
    enddo          ! i
    enddo          ! j

    ! Identify the nodes of each element
    ! Assume that k increases from top to bottom 
    ! (Note: The numbering of nodes within each element is from bottom to top--
    !  this can be confusing)

    NodeOnElement(:,:,:,:) = 0  ! indices are (1:8,k,i,j)

    do j = 1+nhalo, ny-nhalo    ! loop over locally owned cells 
    do i = 1+nhalo, nx-nhalo
       if (active_cell(i,j)) then
          do k = 1, nlyr
             NodeOnElement(1,k,i,j) = nodeID(k+1,i-1,j-1)
             NodeOnElement(2,k,i,j) = nodeID(k+1,i,  j-1)
             NodeOnElement(3,k,i,j) = nodeID(k+1,i,  j)
             NodeOnElement(4,k,i,j) = nodeID(k+1,i-1,j)
             NodeOnElement(5,k,i,j) = nodeID(k,  i-1,j-1)
             NodeOnElement(6,k,i,j) = nodeID(k,  i,  j-1)
             NodeOnElement(7,k,i,j) = nodeID(k,  i,  j)
             NodeOnElement(8,k,i,j) = nodeID(k,  i-1,j)
          enddo  ! k
       endif     ! active cell
    enddo        ! i
    enddo        ! j

    ! Allocate space for the stiffness matrix

    allocate(Auu(3,3,3,nNodes))
    allocate(Auv(3,3,3,nNodes))
    allocate(Avu(3,3,3,nNodes))
    allocate(Avv(3,3,3,nNodes))

    ! Allocate space for the row vector
 
    allocate(bu(nNodes))
    allocate(bv(nNodes))

    !TODO - Start the outer loop here

       ! TODO: Compute the effective viscosity for each element
       ! Note dimensions: visc(nElements)
   

       ! Assemble the stiffness matrix

!TODO - Remove unneeded arguments.

       call assemble_stiffness_matrix(nx, ny, nlyr,  nhalo,           &
                                      active_cell,                     &
                                      nCells,           nActiveCells,  &
                                      NodeOnElement,    &
!                                      iCellIndex,       jCellIndex,  &
                                      iNodeIndex,       jNodeIndex,  &
                                      kNodeIndex,                    &
                                      xNode,            yNode,       &
                                      zNode,                         &
                                      nElements,        nNodes,      &
                                      visc,                          &
                                      Auu,              Auv,         &
                                      Avu,              Avv)

       ! Assemble the load vector

       call assemble_load_vector_2d(nElements,        nNodes,      &
                                    gElementIndex,    gNodeIndex,  &
                                    bu,               bv)

       ! Solve Ax = b


    !whl - End the outer loop here 

    ! Clean up

    deallocate(Auu)
    deallocate(Auv)
    deallocate(Avu)
    deallocate(Avv)
    deallocate(bu)
    deallocate(bv)

  end subroutine glissade_velo_higher_solve

!****************************************************************************

  subroutine assemble_stiffness_matrix(nx, ny, nlyr,  nhalo,           &
                                       active_cell,  &
                                       nCells,           nActiveCells, &
                                       NodeOnElement, &
!                                       iCellIndex,       jCellIndex,   &     
                                       iNodeIndex,       jNodeIndex,   &     
                                       kNodeIndex,   &     
                                       xNode,            yNode,       &
                                       zNode,                         &
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

    logical, intent(in), dimension(nx,ny) :: active_cell

!    integer, dimension(nCells, nVerticesOnCell) ::    &
!       VertexOnCellID        ! VertexID for each vertex of each cell

!    integer, dimension(nCells), intent(in) ::   &
!       iCellIndex, jCellIndex      ! i and j indices of active cells

    integer, dimension(nVertices*nlyr), intent(in) ::   &
       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

    real(dp), dimension(nVertices*nlyr), intent(in) ::   &
       xNode, yNode, zNode  ! x, y and z coordinates of nodes

    integer, dimension(nNodesPerElement,nlyr, nx, ny), intent(in) :: NodeOnElement

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

    integer, dimension(nNodesPerElement) :: nodeID     ! ID for each node of an element

    integer :: i, j, k, n, p

    ! Initialize arrays

    Kuu(:,:) = 0.d0
    Kuv(:,:) = 0.d0
    Kvu(:,:) = 0.d0
    Kvv(:,:) = 0.d0
    
    Auu(:,:,:,:) = 0.d0
    Auv(:,:,:,:) = 0.d0
    Avu(:,:,:,:) = 0.d0
    Avv(:,:,:,:) = 0.d0

    ! Sum over elements

!TODO - Could sum over nCells and do not use i and j?

    do j = nhalo+1, ny-nhalo    ! loop over local cells
    do i = nhalo+1, nx-nhalo
       
     if (active_cell(i,j)) then

       do k = 1, nlyr    ! loop over elements in this column
                         ! assume k increases from upper surface to bed

          ! Find node ID and (x,y,z) coordinates of each node for this element

          do n = 1, nNodesPerElement
             nodeID(n) = NodeOnElement(n,k,i,j)
             x(n) = xNode(nodeID(n))
             y(n) = yNode(nodeID(n))
             z(n) = zNode(nodeID(n))
          enddo   ! nodes per element

          ! Loop over quadrature points for this element
   
          do p = 1, nQuadPoints

          ! Evaluate the derivatives of the element basis functions
          ! at this quadrature point.

             call get_basis_function_derivatives(nNodesPerElement,                             &
                                                 x(:),          y(:),          z(:),           &
                                                 dphi_dxr(:,p), dphi_dyr(:,p), dphi_dzr(:,p),  &
                                                 dphi_dx(:),    dphi_dy(:),    dphi_dz(:),     &
                                                 detJ )

          ! Increment the element stiffness matrix with the contribution from
          ! this quadrature point.

             call increment_element_matrix(nNodesPerElement,                        &
                                           wqp(p),         detJ,                    &
                                           dphi_dx(:),     dphi_dy(:),  dphi_dz(:), &
                                           Kuu(:,:),       Kuv(:,:),                &
                                           Kvu(:,:),       Kvv(:,:) )

          enddo   ! nQuadPoints

          call check_symmetry_element_matrix(nNodesPerElement,  &
                                             Kuu, Kuv, Kvu, Kvv)


         ! If solving in parallel with Trilinos, we have the option at this point to 
         ! call a sum_into_global_matrix routine, passing one row at a time of the
         ! element matrix.  Trilinos will handle the rest.
         !
         ! For the serial SLAP solve, we first form a dense intermediate matrix, 
         ! and later put that matrix in sparse format.

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nodeID,                      &
                                        Kuu,              Auu)

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nodeID,                      &
                                        Kuv,              Auv)

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nodeID,                      &
                                        Kvu,              Avu)

          call element_to_global_matrix(nNodesPerElement, nNodes,    &
                                        nodeID,                      &
                                        Kvv,              Avv)

       ! TODO: Finish assembling the stiffness matrix by incorporating basal sliding?

       enddo   ! nlyr (loop over elements in this column)

     endif   ! active cell

    enddo      ! i
    enddo      ! j

  end subroutine assemble_stiffness_matrix

!****************************************************************************

  subroutine get_basis_function_derivatives(nNodesPerElement,       &
                                            xNode,       yNode,     zNode,    &
                                            dphi_dxr,    dphi_dyr,  dphi_dzr, &
                                            dphi_dx,     dphi_dy,   dphi_dz,  &
                                            detJ )

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
                Jac,     &! Jacobian matrix
                Jinv,    &! inverse Jacobian matrix
                cofactor  ! matrix of cofactors

    integer :: i, j, n

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

    do n = 1, nNodesPerElement
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
    cofactor(3,3) =   Jac(1,2)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    detJ = Jac(1,1)*cofactor(1,1) + Jac(1,2)*cofactor(1,2) + Jac(1,3)*cofactor(1,3)

    if (abs(detJ) > 0.d0) then
       do j = 1, 3
          do i = 1, 3
             Jinv(i,j) = cofactor(j,i)
          enddo
       enddo
       Jinv(:,:) = Jinv(:,:) / detJ
    else
       !whl - do a proper abort here
       print*, 'stopping, det J = 0'
       stop
    endif

!whl - TODO - bug check - Verify that J * Jinv = I

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

  end subroutine get_basis_function_derivatives

!****************************************************************************

  subroutine increment_element_matrix(nNodesPerElement,             &
                                      wqp,      detJ,               &
                                      dphi_dx,  dphi_dy,  dphi_dz,  &
                                      Kuu,      Kuv,                &
                                      Kvu,      Kvv )

!TODO - Add optional arguments that reduce the solver to SIA or SSA.

    !------------------------------------------------------------------
    ! Increment the stiffness matrices Kuu, Kuv, Kvu, Kvv with the
    ! contribution from a particular quadrature point, 
    ! based on the Blatter-Pattyn first-order equations.
    !
    ! This subroutine will work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

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

  end subroutine increment_element_matrix

!****************************************************************************

!TODO - Remove this subroutine?  
! Can always implement SSA in 3D by removing terms from HO approximation       

  subroutine increment_element_matrix_ssa(nNodesPerElement,               &
                                          wqp,               detJ,        &
                                          dphi_dx,           dphi_dy,     &
                                          Kuu,               Kuv,         &
                                          Kvu,               Kvv )

    !------------------------------------------------------------------
    ! Increment the stiffness matrices Kuu, Kuv, Kvu, Kvv with the
    ! contribution from a particular quadrature point, 
    ! based on the equations for the shallow-shelf approximation.
    !
    ! This subroutine will work for any 2D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement  ! number of nodes per element

    real(dp), intent(in) ::    &
             wqp,        &! weight for this quadrature point
             detJ         ! determinant of Jacobian for the transformation
                          !  between the reference element and true element

    real(dp), dimension(nNodesPerElement), intent(in) ::  &
             dphi_dx, dphi_dy       ! derivatives of basis functions,
                                    ! evaluated at this quadrature point

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(inout) :: &
             Kuu, Kuv, Kvu, Kvv     ! components of element stiffness matrix

    integer :: i, j

    ! Increment the element stiffness matrices for the SSA

    do j = 1, nNodesPerElement      ! columns of K
       do i = 1, nNodesPerElement   ! rows of K

          Kuu(i,j) = Kuu(i,j) + wqp*detJ*(4.d0*dphi_dx(j)*dphi_dx(i)  &
                                         +     dphi_dy(j)*dphi_dy(i) )
          Kuv(i,j) = Kuv(i,j) + wqp*detJ*(2.d0*dphi_dx(j)*dphi_dy(i)  &
                                         +     dphi_dy(j)*dphi_dx(i) )
          Kvu(i,j) = Kvu(i,j) + wqp*detJ*(2.d0*dphi_dy(j)*dphi_dx(i)  &
                                         +     dphi_dx(j)*dphi_dy(i) )
          Kvv(i,j) = Kvv(i,j) + wqp*detJ*(4.d0*dphi_dy(j)*dphi_dy(i)  &
                                         +     dphi_dx(j)*dphi_dx(i) )

       enddo  ! i (rows)
    enddo     ! j (columns)

  end subroutine increment_element_matrix_ssa

!****************************************************************************

  subroutine check_symmetry_element_matrix(nNodesPerElement,  &
                                           Kuu, Kuv, Kvu, Kvv)

    !------------------------------------------------------------------
    ! Check that the element stiffness matrix is symmetric
    ! This is true provided that (1) Kuu = (Kuu)^T
    !                            (2) Kvv = (Kvv)^T
    !                            (3) Kuv = (Kvu)^T
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

    real(dp), parameter :: eps11 = 1.e-11

    ! check that Kuu = (Kuu)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kuu(i,j) - Kuu(j,i)) > eps11) then
             print*, 'Kuu is not symmetric'
             print*, 'i, j, Kuu(i,j), Kuu(j,i):', i, j, Kuu(i,j), Kuu(j,i)
             stop
          endif    
       enddo
    enddo

    ! check that Kvv = (Kvv)^T

    do j = 1, nNodesPerElement
       do i = j, nNodesPerElement
          if (abs(Kvv(i,j) - Kvv(j,i)) > eps11) then
             print*, 'Kvv is not symmetric'
             print*, 'i, j, Kvv(i,j), Kvv(j,i):', i, j, Kvv(i,j), Kvv(j,i)
             stop
          endif    
       enddo
    enddo

    ! Check that Kuv = (Kvu)^T

    do j = 1, nNodesPerElement
       do i = 1, nNodesPerElement
          if (abs(Kuv(i,j) - Kvu(j,i)) > eps11) then
             print*, 'Kuv .ne. (Kvu)^T'
             print*, 'i, j, Kuv(i,j), Kvu(j,i):', i, j, Kuv(i,j), Kvu(j,i)
             stop
          endif    
       enddo
    enddo

  end subroutine check_symmetry_element_matrix

!****************************************************************************

  subroutine element_to_global_matrix(nNodesPerElement, nNodes,   &
                                      nodeID,                     &
                                      Kmat,             Amat)

    ! Sum terms of element matrix K into dense assembled matrix A
    ! Here we assume that K is partitioned into Kuu, Kuv, Kvu, and Kvv,
    !  and similarly for A.
    ! So this subroutine is called four times.

    integer, intent(in) ::   &
       nNodesPerElement,     &
       nNodes            ! total number of active nodes
       
    integer, dimension(nNodesPerElement), intent(in) ::  &
       nodeID            ! node ID for each node of this element

    real(dp), dimension(nNodesPerElement,nNodesPerElement), intent(in) ::  &
       Kmat              ! element matrix

    real(dp), dimension(3,3,3,nNodes), intent(inout) ::    &
       Amat              ! assembled matrix

    integer :: m, n, nid, iA, jA, kA

    do m = 1, nNodesPerElement       ! rows of K
       nid = nodeID(m)               ! node ID for this row of K

       do n = 1, nNodesPerElement    ! columns of K

          kA = kshift(m,n)           ! k index of A into which K(m,n) is summed
          iA = ishift(m,n)           ! similarly for i and j indices
          jA = jshift(m,n)         

          !TODO - Add viscosity term
          Amat(kA,iA,jA,nid) = Amat(kA,iA,jA,nid) + Kmat(m,n)

       enddo     ! n
    enddo        ! m

  end subroutine element_to_global_matrix

!****************************************************************************

!TODO - Make this 3D

  subroutine assemble_load_vector_2d(nElements,        nNodes,      &
                                     gElementIndex,    gNodeIndex,  &
                                     bu,               bv)

    integer, intent(in) ::   &
       nElements,    &    ! number of elements on this processor
       nNodes             ! number of nodes associated with these elements

    integer, dimension(nElements), intent(in) ::   &
       gElementIndex      ! global cell index for each element

    integer, dimension(nVertices), intent(in) ::   &
       gNodeIndex         ! global node index for each vertex

    real(dp), dimension(nNodes) ::  &
       bu, bv             ! assembled load vector, divided into 2 parts

    bu(:) = 0.d0
    bv(:) = 0.d0

!whl - Write this subroutine

!      First the gravitational term

!      Then the pressure term from immersed shelves


  end subroutine assemble_load_vector_2d

!****************************************************************************

!whl - copied much of this subroutine from glam_velo_solver
!      Not sure how much of this I might need. 

! This is the driver subroutine, called from 'run_ho_diagnostic' in
! 'glide_velo_higher.F90'. In turn, 'run_ho_diagnostic' is called from 'inc_remap_driver' in
! 'glam.F90', and 'inc_remap_driver' is called from 'glide_tstep_p2' in 'glide.F90'.

subroutine glissade_velo_solver_template( ewn,      nsn,    upn,  &
                                 dew,      dns,          &
                                 sigma,    stagsigma,    &
                                 thck,     usrf,         &
                                 lsrf,     topg,         &
                                 dthckdew, dthckdns,     &
                                 dusrfdew, dusrfdns,     &
                                 dlsrfdew, dlsrfdns,     &
                                 stagthck, flwa,         &
                                 mintauf,                &
                                 btraction,              &
                                 umask,                  &
                                 whichbabc,              &
                                 whichefvs,              &
                                 whichresid,             &
                                 whichnonlinear,         &
                                 whichsparse,            &
                                 periodic_ew,periodic_ns,&
                                 beta,                   &
                                 uvel,     vvel,         &
                                 uflx,     vflx,         &
                                 efvs,     tstep )

  implicit none

  integer, intent(in) :: ewn, nsn, upn, tstep  ! JFL to be removed
  integer, dimension(:,:),   intent(inout)  :: umask
  ! NOTE: 'inout' status to 'umask' should be changed to 'in' at some point,
  ! but for now this allows for some minor internal hacks to CISM-defined mask

  real (kind = dp), intent(in) :: dew, dns                                  ! grid dimensions

  real (kind = dp), dimension(:),     intent(in)  :: sigma, stagsigma       ! sigma coords
  real (kind = dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg ! geom vars
  real (kind = dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns     ! thick grads
  real (kind = dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns     ! upper surf grads
  real (kind = dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns     ! basal surf grads
  real (kind = dp), dimension(:,:),   intent(in)  :: stagthck               ! staggered thickness
  real (kind = dp), dimension(:,:),   intent(in)  :: minTauf                ! till yield stress
  real (kind = dp), dimension(:,:,:), intent(inout) :: btraction            ! consistent basal traction array
  real (kind = dp), dimension(:,:,:), intent(in)  :: flwa                   ! flow law rate factor

  real (kind = dp), dimension(:,:),   intent(in)  :: beta                   ! basal friction factor

  integer, intent(in) :: whichbabc    ! options for betasquared field to use
  integer, intent(in) :: whichefvs    ! options for efvs calculation (calculate it or make it uniform)
  integer, intent(in) :: whichresid   ! options for method to use when calculating vel residul
  integer, intent(in) :: whichnonlinear  ! options for which method for doing elliptic solve
  integer, intent(in) :: whichsparse  ! options for which method for doing elliptic solve

  logical, intent(in) :: periodic_ew, periodic_ns  ! options for applying periodic bcs or not
  real (kind = dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)
  real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx  ! horiz fluxs: u_bar*H, v_bar*H
  real (kind = dp), dimension(:,:,:), intent(out) :: efvs        ! effective viscosity

  integer :: ew, ns, up     ! counters for horiz and vert do loops

  real (kind = dp), parameter :: minres = 1.0d-4    ! assume vel fields converged below this resid
  real (kind = dp), parameter :: NL_tol = 1.0d-06   ! to have same criterion
                                                    ! than with JFNK
  real (kind = dp), save, dimension(2) :: resid     ! vector for storing u resid and v resid
  real (kind = dp) :: plastic_resid_norm = 0.0d0    ! norm of residual used in Newton-based plastic bed iteration

  integer, parameter :: cmax = 100                  ! max no. of iterations
  integer :: counter, linit                         ! iteration counter
  character(len=100) :: message                     ! error message

  ! variables used for incorporating generic wrapper to sparse solver
  type(sparse_matrix_type) :: matrix

  end subroutine glissade_velo_solver_template

!****************************************************************************

!whl - This subroutine to be removed

  subroutine fill_matrix_stress_2d (nVertices,  nCells,      &
                                    nVerticesOnVertex,       &
                                    nVerticesOnCell,         &
                                    gVerticesOnCell,         &
                                    efvs,                    &
                                    kuu, kuv, kvu, kvv,      &
                                    auu, auv, avu, avv )

    integer, intent(in) :: nVertices  ! number of locally owned vertices

    integer, intent(in) :: ncells     ! number of cells containing a locally owned vertex

    integer, intent(in) :: nVerticesOnVertex   ! number of vertices bordering each vertex

    integer, intent(in) :: nVerticesOnCell     ! number of vertices of each cell

    integer, dimension(:,:), intent(in) ::    &
       gVerticesOnCell       ! global index of vertices of a given cell

    real(dp), dimension(nCells) ::   &
          efvs                        ! effective viscosity

    real(dp), dimension(nVerticesOnCell, nVerticesonCell), intent(in)  ::   &
          kuu, kuv, kvu, kvv          ! 4 blocks of element stiffness matrix

!whl - Remove the zero dimension?
    real(dp), dimension(0:nVertices, nVerticesonVertex), intent(inout)  ::   &
          auu, auv, avu, avv          ! arrays for assembling nonzero matrix elements
                                      ! Array dimensions start with zero to allow space for 
                                      !  vertices with global index = 0

    integer :: nb, nc, nv1, nv2, nv3, nv4

!whl - this code probably to be removed
       do nc = 1, nCells

          ! Determine indices for the 4 vertices of this cell
          ! Note: will have global index = 0 for halo vertices

!whl - Vertex index does not have to be global here?

          nv1 = gVerticesOnCell(nc,1)
          nv2 = gVerticesOnCell(nc,2)
          nv3 = gVerticesOnCell(nc,3)
          nv4 = gVerticesOnCell(nc,4)

!whl - Rearrange for stride 1?

          ! Accumulate terms in auu (dI/du or dI/dv terms that are multiplied by u or v)

!          auu(nv1,i00) = auu(nv1,i00) + efvs(nc) * kuu(1,1) 
!          auu(nv1,i0e) = auu(nv1,i0e) + efvs(nc) * kuu(1,2) 
!          auu(nv1,ine) = auu(nv1,ine) + efvs(nc) * kuu(1,3) 
!          auu(nv1,in0) = auu(nv1,in0) + efvs(nc) * kuu(1,4) 

!          auu(nv2,i0w) = auu(nv2,i0w) + efvs(nc) * kuu(2,1) 
!          auu(nv2,i00) = auu(nv2,i00) + efvs(nc) * kuu(2,2) 
!          auu(nv2,in0) = auu(nv2,in0) + efvs(nc) * kuu(2,3) 
!          auu(nv2,inw) = auu(nv2,inw) + efvs(nc) * kuu(2,4) 

!          auu(nv3,isw) = auu(nv3,isw) + efvs(nc) * kuu(3,1) 
!          auu(nv3,is0) = auu(nv3,is0) + efvs(nc) * kuu(3,2) 
!          auu(nv3,i00) = auu(nv3,i00) + efvs(nc) * kuu(3,3) 
!          auu(nv3,i0w) = auu(nv3,i0w) + efvs(nc) * kuu(3,4) 

!          auu(nv4,is0) = auu(nv4,is0) + efvs(nc) * kuu(4,1) 
!          auu(nv4,ise) = auu(nv4,ise) + efvs(nc) * kuu(4,2) 
!          auu(nv4,i0e) = auu(nv4,i0e) + efvs(nc) * kuu(4,3) 
!          auu(nv4,i00) = auu(nv4,i00) + efvs(nc) * kuu(4,4) 

!whl - then would need to repeat for auv, avu, avv

       enddo   ! nCells

  end subroutine fill_matrix_stress_2d

!****************************************************************************

!whl - This subroutine to be renamed and rewritten.
!      Need something like this to put terms of Auu in sparse matrix format

  subroutine assemble_matrix (nVertices,          &
                              nVerticesOnVertex,  &
                              gVerticesOnVertex,  &
                              auu)

!whl - Move putpcgc elsewhere?

     use glam_strs2, only: putpcgc

     integer, intent(in) ::  &
         nVertices,          &! number of locally owned vertices
         nVerticesOnVertex    ! number of vertices in neighborhood of each vertex
                              ! 9 in 2d, 27 in 3d

     integer, dimension(nVertices, nVerticesOnVertex), intent(in) ::  &
         gVerticesOnVertex    ! vertex index for the neighbors of a given vertex

     real(dp), dimension(0:nVertices, nVerticesonVertex), intent(in)  ::   &
          auu                         ! array for accumulating nonzero matrix elements
                                      ! Array dimensions start with zero to allow space for 
                                      !  vertices with global index = 0

     integer, dimension (4) ::   &
         row_offset, col_offset

     integer :: nb, nv, i

     integer :: row, col

     real(dp) :: val

     row_offset(1) = 0            ! upper left
     row_offset(2) = nVertices    ! upper right
     row_offset(3) = 0            ! lower left
     row_offset(4) = nVertices    ! lower right

     col_offset(1) = 0            ! upper left
     col_offset(2) = 0            ! upper right
     col_offset(3) = nVertices    ! lower left
     col_offset(4) = nVertices    ! lower right


!whl - This code probably will be removed
     nb = 1
     do i = 1, nVerticesOnVertex
        do nv = 1, nVertices
           row = nv + row_offset(nb)
           col = gVerticesOnVertex(nv,i) + col_offset(nb)
           val = auu(nv,i)
           call putpcgc(val, col, row)
        enddo
     enddo

!whl - would need to add auv, avu, etc.


  end subroutine assemble_matrix

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
