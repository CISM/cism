! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! +                                                           + 
! +  glissade_velo_higher.F90                                 + 
! +                                                           + 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains routines for computing the ice velocity
!  using a variational finite-element approach.
! See Dukowicz, Price and Lipscomb (J. Glac., 2010).
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

!TODO - At some point, we coud make this code MPAS-friendly by replacing 
!       the ij loops with indirect addressing.  But keep i and j for now.

    !----------------------------------------------------------------
    ! basic cell properties (quadrilateral cells)
    !
    !TODO - might want to declare these elsewhere (they would be declared
    !       elsewhere in MPAS code)
    !----------------------------------------------------------------

    integer :: nCells     ! total number of cells on this processor
    integer :: nVertices  ! total number of vertices associated with at least one cell

    integer, parameter ::      &
       nVerticesOnCell = 4,    &   
       nVerticesOnVertex = 9 

!TODO - Use VertexID instead of gVertex Index?
    integer, dimension(:,:), allocatable ::    &
!       gCellIndex,          &! global cell index for each cell (i,j)
!       gVertexIndex,        &!
       CellID,              &! unique local ID for each cell (i,j)
       VertexID              ! unique local ID for each vertex (i,j)

!whl - Can remove gVerticesOnVertex
!      Might want to replace gVerticesOnCell with gNodesOnElement in FE subroutines
  
    integer, dimension(:,:), allocatable ::    &
       VertexOnCellID        ! VertexID for each vertex of each cell

!    integer, dimension(:,:), allocatable ::    &
!       gVerticesOnVertex,   &! global index of vertices neighboring a given vertex
!       gVerticesOnCell       ! global index of vertices of a given cell

    real(dp), dimension(:), allocatable ::     &
       xVertex, yVertex, zvertex    ! x, y and z coordinates of each vertex
                                    ! x and y are fixed, but z will change in time
    !----------------------------------------------------------------
    ! Finite element properties
    !
    ! For now, assume 3D hexahedral elements.
    !----------------------------------------------------------------

    integer, parameter ::      &
       nNodesPerElement = 8,   &   ! 8 nodes for hexahedra
       nQuadPoints = 8             ! number of quadrature points
                                   ! These live at +- 1/sqrt(3) for hexahedra

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

! storage arrays

!TODO - Get rid of these
! Would be nice to get rid of these (assemble Auu, etc. in some other way)

    integer, parameter :: &   ! indices for a vertex and its 8 horizontal neighbors
       isw = 1,    &
       is0 = 2,    &
       ise = 3,    &
       i0w = 4,    &
       i00 = 5,    &
       i0e = 6,    &
       inw = 7,    &
       in0 = 8,    &
       ine = 9

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


    !----------------------------------------------------------------
    ! Here are some definitions:
    !
    ! The horizontal mesh is composed of cells and vertices.
    ! Each cell can be extruded to form a column with a specified number of layers.
    ! 
    ! An element is a layer of a cell, and a node is a corner of an element.
    !
    ! Active elements are those elements lying in cells with ice present
    !  (that is, the ice thickness is above a minimum threshold).
    ! Active nodes are the nodes of these active elements.
    !
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx,  ny,              &     ! number of grid cells in each direction
       nhalo                       ! number of halo layers

    real(dp), intent(in) ::  &
       dx,  dy                     ! grid cell dimensions


    integer :: i, j, n, p, nc, nv

    !----------------------------------------------------------------
    ! Assign indices for the cells and vertices on this mesh.
    !----------------------------------------------------------------

    ! Define global cell index for cells on this processor.
    ! Note: Cell 1 is at lower left, cell nCells is at upper right

!TODO - The following indexing scheme works only in serial.  
!       Must rewrite (or let Trilinos handle the indexing) if in parallel.

!whl - MPAS has indexToCellID and indexToVertexID.  
!      Could these replace gCellIndex and gVertexIndex?

    ! Assign a unique index to each cell on this processor.

!    allocate(gCellIndex(nx,ny))
!    gCellIndex(:,:) = 0

!    nCells = 0
!    do j = 1+nhalo, ny-nhalo
!    do i = 1+nhalo, nx-nhalo     
!       nCells = nCells + 1
!       gCellIndex(i,j) = nCells
!    enddo
!    enddo       

    
    allocate(CellID(nx,ny))
    nCells = 0
    do j = 1, ny   ! OK to count halo cells?  Or should loop be over local cells only?
    do i = 1, nx
       nCells = nCells + 1
       CellID(i,j) = nCells
    enddo
    enddo

    ! Assign a unique index to each vertex on this processor
    ! By convention, vertex (i,j) lies at the NE corner of cell(i,j).

!    allocate(gVertexIndex(nx-1,ny-1))
!    gVertexIndex(:,:) = 0

!    nVertices = 0
!    do j = nhalo, ny-nhalo   ! include south edge (index nhalo)     
!    do i = nhalo, nx-nhalo   ! include west edge  (index nhalo)
!       nVertices = nVertices + 1
!       gVertexIndex(i,j) = nVertices
!    enddo
!    enddo       

    allocate(VertexID(nx-1, ny-1))
    nVertices = 0

    do j = 1, ny-1
    do i = 1, nx-1
       nVertices = nVertices + 1
       VertexID(i,j) = nVertices
    enddo
    enddo

    ! Determine the vertex IDs for each vertex of each cell.
    ! Numbering is counter-clockwise from SW corner
    ! Requires nhalo >= 1

    allocate(VertexOnCellID(nCells, nVerticesOnCell))
    VertexOnCellID(:,:) = 0

    do j = 1+nhalo, ny-nhalo   ! locally owned cells
    do i = 1+nhalo, nx-nhalo
       nc = CellID(i,j)     
       VertexOnCellID(nc,1) = VertexID(i-1,j-1)   ! lower SW corner
       VertexOnCellID(nc,2) = VertexID(i,j-1)     ! lower SE corner
       VertexOnCellID(nc,3) = VertexID(i,j)       ! lower NE corner
       VertexOnCellID(nc,4) = VertexID(i-1,j)     ! lower NW corner
!       VertexOnCellID(nc,5) = VertexID(i-1,j-1)   ! upper SW corner
!       VertexOnCellID(nc,6) = VertexID(i,j-1)     ! upper SE corner
!       VertexOnCellID(nc,7) = VertexID(i,j)       ! upper NE corner
!       VertexOnCellID(nc,8) = VertexID(i-1,j)     ! upper NW corner
    enddo
    enddo
    
    ! Assign x and y coordinates to each vertex.
    ! z values set to zero for now; will be set later based on the ice thickness.

    !TODO - These coordinates are relative to local processor; modify for parallel case?

    allocate(xVertex(nVertices))
    allocate(yVertex(nVertices))
    allocate(zVertex(nVertices))

    xVertex(:) = 0.d0
    yVertex(:) = 0.d0
    zVertex(:) = 0.d0

    do j = 1, ny-1
    do i = 1, nx-1
       xVertex(VertexID(i,j)) = dx * i
       yVertex(VertexID(i,j)) = dy * j
    enddo
    enddo

!TODO: Not sure the following is needed.
    ! Determine the global indices of each vertex and its 8 neighbors
    ! Note: This requires nhalo >= 2.  Must revise if nhalo = 1.

!    allocate(gVerticesOnVertex(nVertices, nVerticesOnVertex))
!    gVerticesOnVertex(:,:) = 0

!    do j = nhalo, ny-nhalo    ! include south edge
!    do i = nhalo, nx-nhalo    ! include west edge
!       nv = gVertexIndex(i,j)
!       gVerticesOnVertex(nv,isw) = gVertexIndex(i-1,j-1) 
!       gVerticesOnVertex(nv,is0) = gVertexIndex(i  ,j-1) 
!       gVerticesOnVertex(nv,ise) = gVertexIndex(i+1,j-1) 
!       gVerticesOnVertex(nv,i0w) = gVertexIndex(i-1,j) 
!       gVerticesOnVertex(nv,i00) = gVertexIndex(i  ,j)
!       gVerticesOnVertex(nv,i0e) = gVertexIndex(i+1,j) 
!       gVerticesOnVertex(nv,inw) = gVertexIndex(i-1,j+1) 
!       gVerticesOnVertex(nv,in0) = gVertexIndex(i  ,j+1) 
!       gVerticesOnVertex(nv,ine) = gVertexIndex(i+1,j+1) 
!    enddo
!    enddo

    !----------------------------------------------------------------
    ! Initialize some time-independent finite element arrays
    !----------------------------------------------------------------

    ! Trilinear basis set for reference hexahedron, x=(-1,1), y=(-1,1), z=(-1,1)             
    ! Indexing is counter-clockwise from SW corner, with 1-4 on lower surface
    !  and 5-8 on upper surface
    !
    ! N1 = (1-x)*(1-y)*(1-z)/8
    ! N2 = (1+x)*(1-y)*(1-z)/8
    ! N3 = (1+x)*(1+y)*(1-z)/8
    ! N4 = (1-x)*(1+y)*(1-z)/8
    ! N5 = (1-x)*(1-y)*(1+z)/8
    ! N6 = (1+x)*(1-y)*(1+z)/8
    ! N7 = (1+x)*(1+y)*(1+z)/8
    ! N8 = (1-x)*(1+y)*(1+z)/8
   
    ! Set coordinates and weights of quadrature points for reference square element
    ! Numbering is counter-clockwise from southwest

    xqp(1) = -rsqrt3
    yqp(1) = -rsqrt3
    zqp(1) = -rsqrt3
    wqp(1) =  1.d0   !TODO - check that weight = 1 and not 1/8

    xqp(2) =  rsqrt3
    yqp(2) = -rsqrt3
    zqp(2) = -rsqrt3
    wqp(2) =  1.d0

    xqp(3) =  rsqrt3
    yqp(3) =  rsqrt3
    zqp(3) = -rsqrt3
    wqp(3) =  1.d0

    xqp(4) = -rsqrt3
    yqp(4) =  rsqrt3
    zqp(4) = -rsqrt3
    wqp(4) =  1.d0

    xqp(5) = -rsqrt3
    yqp(5) = -rsqrt3
    zqp(5) =  rsqrt3
    wqp(5) =  1.d0

    xqp(6) =  rsqrt3
    yqp(6) = -rsqrt3
    zqp(6) =  rsqrt3
    wqp(6) =  1.d0

    xqp(7) =  rsqrt3
    yqp(7) =  rsqrt3
    zqp(7) =  rsqrt3
    wqp(7) =  1.d0

    xqp(8) = -rsqrt3
    yqp(8) =  rsqrt3
    zqp(8) =  rsqrt3
    wqp(8) =  1.d0

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

    enddo

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

    real (dp), dimension(nx,ny,nlyr), intent(in) ::  &
       flwa                   ! flow factor

    integer, intent(in) ::   &
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

!    integer, dimension(nCells) ::   &
!       iCellIndex, jCellIndex ! i and j indices of active cells

!    integer, dimension(nVertices*nlyr) ::   &
!       iNodeIndex, jNodeIndex, kNodeIndex   ! i, j and k indices of nodes

    logical, dimension(nx,ny) :: &
       active_cell,       &   ! true for active cells (thck > thckmin)
       active_vertex          ! true for vertices of active cells

    logical, dimension(nlyr+1,nx,ny) ::  &
       active_node            ! true for active nodes

    real(dp), dimension(:), allocatable :: &
       visc,           &   ! effective viscosity
       vertex_usfc,    &   ! upper surface averaged to vertices
       vertex_thck         ! ice thickness averaged to vertices

    real(dp), dimension(:,:,:,:), allocatable ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv           ! 1st dimension = 3 (node and its 2 neighbors in z direction) 
                          ! 2nd dimension = 3 (node and its 2 neighbors in x direction) 
                          ! 3rd dimension = 3 (node and its 2 neighbors in y direction) 
                          ! 4th dimension = nNodes

    real(dp), dimension(:), allocatable ::  &
       bu, bv             ! assembled load vector, divided into 2 parts

    integer :: i, j, k, gc, gv, nv

!TODO - Are these needed?
    integer, dimension(nCells) :: gElementIndex   ! global cell index for each element
    integer, dimension(nVertices) :: gNodeIndex   ! global node index for each vertex


    ! Identify and count the active cells (i.e., cells with thck > thckmin)

    nActiveCells = 0
    active_cell(:,:) = .false.
!    iCellIndex(:) = 0
!    jCellIndex(:) = 0

!TODO - May not need iCellIndex is using CellID(i,j)?
    do j = 1+nhalo, ny-nhalo
    do i = 1+nhalo, nx-nhalo
       if (thck(i,j) >= thckmin) then
          active_cell(i,j) = .true.
          nActiveCells = nActiveCells + 1
!          iCellIndex(nActiveCells) = i
!          jCellIndex(nActiveCells) = j
       endif
    enddo
    enddo
        
    ! Count the active elements

    nElements = nActiveCells * nlyr

    ! Identify the active vertices (i.e., all vertices of active cells)
    ! Many of these nodes will be set to true multiple times.

    active_vertex(:,:) = .false.

    do j = 1+nhalo, ny-nhalo
    do i = 1+nhalo, nx-nhalo
       if (active_cell(i,j)) then
          active_vertex(i-1:i, j-1:j) = .true.  ! 4 vertices of this cell
       endif
    enddo
    enddo

    ! Identify and count the active nodes

    active_node(:,:,:) = .false.
    nNodes = 0

    do j = nhalo, ny-nhalo    ! include S edge
    do i = nhalo, nx-nhalo    ! include W edge
       if (active_vertex(i,j)) then
          do k = 1, nlyr+1    ! all nodes in column are active
             active_node(k,i,j) = .true.
             nNodes = nNodes + 1   
!             iNodeIndex(nNodes) = i
!             jNodeIndex(nNodes) = j
!             kNodeIndex(nNodes) = k
          enddo
       endif
    enddo
    enddo

!TODO - Is this needed?
    ! Identify the nodes of each element

    do j = 1+nhalo, ny-nhalo
    do i = 1+nhalo, nx-nhalo
        if (active_cell(i,j)) then
           do k = 1, nlyr
!              Node_On_Element(i,j,k,1) = nodeID(i-1,j-1,k-1)
           enddo
        endif
    enddo
    enddo


    ! Allocate space for the stiffness matrix

    allocate(Auu(3,3,3,nNodes))
    allocate(Auv(3,3,3,nNodes))
    allocate(Avu(3,3,3,nNodes))
    allocate(Avv(3,3,3,nNodes))

    ! Allocate space for the row vector
 
    allocate(bu(nNodes))
    allocate(bv(nNodes))

    ! allocate space for the effective viscosity
    allocate(visc(nElements))

    !TODO - Start the outer loop here


       ! Compute the thickness at each vertex
       do j = 1, ny-1
       do i = 1, nx-1
          nv = VertexID(i,j)
          vertex_usfc(nv) = stagusfc(i,j)
          vertex_thck(nv) = stagthck(i,j)
       enddo
       enddo

       ! TODO: Compute the effective viscosity for each element
       ! Note dimensions: visc(nElements)
   

       ! Assemble the stiffness matrix

       call assemble_stiffness_matrix(nx, ny, nlyr,  nhalo,           &
                                      active_cell,                     &
                                      nCells,           nActiveCells,  &
                                      VertexOnCellID,                 &
!                                      iCellIndex,       jCellIndex,  &
!                                      iNodeIndex,       jNodeIndex,  &
!                                      kNodeIndex,                    &
                                      nElements,        nNodes,      &
                                      visc,             sigma,       &
                                      vertex_usfc,      vertex_thck, &
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
    deallocate(visc)

  end subroutine glissade_velo_higher_solve

!****************************************************************************

  subroutine assemble_stiffness_matrix(nx, ny, nlyr,  nhalo,           &
                                       active_cell,  &
                                       nCells,           nActiveCells, &
                                       VertexOnCellID, &
!                                       iCellIndex,       jCellIndex,   &     
!                                       iNodeIndex,       jNodeIndex,   &     
!                                       kNodeIndex,   &     
                                       nElements,        nNodes,       &
                                       visc,             sigma,        &
                                       vertex_usfc,      vertex_thck,  &
                                       Auu,              Auv,          &
                                       Avu,              Avv)

    integer, intent(in) ::      &
       nx, ny,                  &    ! horizontal grid dimensions
       nlyr,                    &    ! number of vertical layers
       nhalo,                   &    ! number of halo layers
       nCells,                  &    ! number of horizontal cells
       nActiveCells                  ! number of cells with ice present

    logical, intent(in), dimension(nx,ny) :: active_cell

    integer, dimension(nCells, nVerticesOnCell) ::    &
       VertexOnCellID        ! VertexID for each vertex of each cell

!    integer, dimension(nCells), intent(in) ::   &
!       iCellIndex, jCellIndex      ! i and j indices of active cells

!    integer, dimension(nVertices*nlyr), intent(in) ::   &
!       iNodeIndex, jNodeIndex, kNodeIndex  ! i, j and k indices of nodes

!TODO - Are these needed?
    integer, intent(in) ::   &
       nElements,    &    ! number of elements
       nNodes             ! number of nodes associated with these elements

    real(dp), dimension(nElements), intent(in) ::  &
       visc               ! effective viscosity

    real(dp), dimension(nlyr+1), intent(in) ::   &
       sigma              ! sigma vertical coordinate

    real(dp), dimension(nVertices), intent(in) ::  &
       vertex_usfc,     & ! upper surface averaged to vertices
       vertex_thck        ! ice thickness averaged to vertices

    real(dp), dimension(3,3,3,nNodes), intent(out) ::  &
       Auu, Auv,    &     ! assembled stiffness matrix, divided into 4 parts
       Avu, Avv                                    

    real(dp), dimension(nNodesPerElement) ::   &
       xNode, yNode, zNode  ! x, y and z coordinates of nodes

    real(dp), dimension(2,2) ::   &
       Jac, Jinv          ! Jacobian matrix and its inverse

    real(dp) ::   &
       detJ               ! determinant of J

    real(dp), dimension(nNodesPerElement) ::   &
       dphi_dx, dphi_dy, dphi_dz   ! derivatives of basis function, evaluated at quad pts

    real(dp), dimension(nNodesPerElement, nNodesPerElement) ::   &   !
       Kuu,          &    ! element stiffness matrix, divided into 4 parts as shown below
       Kuv,          &    !  
       Kvu,          &    !
       Kvv                !    Kuu  | Kuv
                          !    _____|____          
                          !    Kvu  | Kvv
                          !         |
                          ! Kvu may not be needed if matrix is symmetric, but is included for now

    integer :: ii, jj, i, j, k, n, nc, nv, gc, gv, p, r, s, vid

    ! Initialize arrays

    Auu(:,:,:,:) = 0.d0
    Auv(:,:,:,:) = 0.d0
    Avu(:,:,:,:) = 0.d0
    Avv(:,:,:,:) = 0.d0

    Kuu(:,:) = 0.d0
    Kuv(:,:) = 0.d0
    Kvu(:,:) = 0.d0
    Kvv(:,:) = 0.d0
    
    ! Sum over elements

!TODO - Sum over nCells and do not use i and j?

    do jj = nhalo+1, ny-nhalo    ! loop over local cells
    do ii = nhalo+1, nx-nhalo
       
     if (active_cell(ii,jj)) then

       nc = CellID(ii,jj)
          
       do k = 1, nlyr    ! loop over elements in this column
                         ! assume k increases from upper surface to bed

          ! Find x, y and z coordinates of each node for this element
          ! (need to check this)

          do n = 1, nNodesPerElement
             nv = mod(n, nVerticesOnCell)  ! assume consecutive numbering of each layer 
                                           ! e.g., 1-4 on bottom layer and 5-8 on top layer
             vid = VertexOnCellID(nc,nv)
             xNode(n) = xVertex(vid)
             yNode(n) = yVertex(vid)
             if (n <= nVerticesOnCell) then   ! lower surface
                zNode(n) = vertex_usfc(vid) - sigma(k+1)*vertex_thck(vid)
             else
                zNode(n) = vertex_usfc(vid) - sigma(k)*vertex_thck(vid)   
             endif
          enddo   ! nodes per element

       ! Loop over quadrature points
   
       do p = 1, nQuadPoints

!TODO - Modify to include z coordinate

          ! Evaluate the derivatives of the element basis functions
          ! at this quadrature point.

          call get_basis_function_derivatives_2d(nNodesPerElement,               &
                                                 xNode(:),       yNode(:),       &
                                                 dphi_dxr(:,p),  dphi_dyr(:,p),  &
                                                 dphi_dx(:),     dphi_dy(:),     &
                                                 detJ )

          ! Increment the element stiffness matrix with the contribution from
          ! this quadrature point, assuming we are solving the shallow-shelf equations.

          call increment_element_matrix_ssa(nNodesPerElement,               &
                                            wqp(p),            detJ,        &
                                            dphi_dx(:),        dphi_dy(:),  &
                                            Kuu(:,:),          Kuv(:,:),    &
                                            Kvu(:,:),          Kvv(:,:) )

       enddo   ! nQuadPoints

       call check_symmetry_element_matrix(nNodesPerElement,  &
                                          Kuu, Kuv, Kvu, Kvv)

!       call element_to_global_matrix_2d

       ! Insert terms of element stiffness matrices (Kuu, Kuv, Kvu, Kvv) into the 
       ! global stiffness matrices (Auu, Auv, Avu, Avv).
       ! 
       ! The integrated basis functions (computed above) are multiplied by the
       ! viscosity of each element.
       !
       ! The number of rows in the global matrices is nNodes.
       ! The number of columns is 9, corresponding to a vertex and its 8 neighbors. 
       ! In this way we can accumulate and store all the nonzero terms of A in a compact matrix.

       !whl - It might be better if we could store A in the sparse matrix data type.
       !      But I don't know how to accumulate values with this data type; I only know
       !       how to put values that are already accumulated.
       !
       !whl - Change i0e to ie, in0 to in, etc.?

!       i = 1                          ! 1st row of K matrix, associated with vertex 1
!       gv = gVerticesOnCell(gc,i)     ! global vertex index for this node (gc = cell index)

!        gv = 1 ! just for now
!        r = 1  ! just for now
!       r = gNodeIndex(gv)             ! global matrix row corresponding to this node

!       Auu(r,i00) = Auu(r,i00) + Kuu(1,1)*visc(ne)  ! column corresponding to vertex 1
!       Auu(r,i0e) = Auu(r,i0e) + Kuu(1,2)*visc(ne)  ! column corresponding to vertex 2
!       Auu(r,ine) = Auu(r,ine) + Kuu(1,3)*visc(ne)  ! column corresponding to vertex 3
!       Auu(r,in0) = Auu(r,in0) + Kuu(1,4)*visc(ne)  ! column corresponding to vertex 4

!       i = 2                          ! 2nd row of K matrix, associated with vertex 2
!       gv = gVerticesOnCell(gc,i)     ! global vertex index for this node
!       r = gNodeIndex(gv)             ! global matrix row corresponding to this node

!       Auu(r,i0w) = Auu(r,i0w) + Kuu(2,1)*visc(ne)  ! column corresponding to vertex 1
!       Auu(r,i00) = Auu(r,i00) + Kuu(2,2)*visc(ne)  ! column corresponding to vertex 2
!       Auu(r,in0) = Auu(r,in0) + Kuu(2,3)*visc(ne)  ! column corresponding to vertex 3
!       Auu(r,inw) = Auu(r,inw) + Kuu(2,4)*visc(ne)  ! column corresponding to vertex 4

!       i = 3                          ! 3rd row of K matrix, associated with vertex 3
!       gv = gVerticesOnCell(gc,i)     ! global vertex index for this node
!       r = gNodeIndex(gv)             ! global matrix row corresponding to this node

!       Auu(r,isw) = Auu(r,isw) + Kuu(3,1)*visc(ne)  ! column corresponding to vertex 1
!       Auu(r,is0) = Auu(r,is0) + Kuu(3,2)*visc(ne)  ! column corresponding to vertex 2
!       Auu(r,i00) = Auu(r,i00) + Kuu(3,3)*visc(ne)  ! column corresponding to vertex 3
!       Auu(r,i0w) = Auu(r,i0w) + Kuu(3,4)*visc(ne)  ! column corresponding to vertex 4

!       i = 4                          ! 3rd row of K matrix, associated with vertex 3
!       gv = gVerticesOnCell(gc,i)     ! global vertex index for this node
!       r = gNodeIndex(gv)             ! global matrix row corresponding to this node

!       Auu(r,is0) = Auu(r,is0) + Kuu(4,1)*visc(ne)  ! column corresponding to vertex 1
!       Auu(r,ise) = Auu(r,ise) + Kuu(4,2)*visc(ne)  ! column corresponding to vertex 2
!       Auu(r,i0e) = Auu(r,i0e) + Kuu(4,3)*visc(ne)  ! column corresponding to vertex 3
!       Auu(r,i00) = Auu(r,i00) + Kuu(4,4)*visc(ne)  ! column corresponding to vertex 4

!whl - Then need to repeat for Auv, Avu, and Avv
!      How to minimize the number of lines of code?
!      Maybe turn this into a subroutine with inputs Auu/Kuu, Auv/Kuv, etc.?

       ! Here I should finish assembling the stiffness matrix by incorporating basal sliding?


       enddo   ! nlyr

     endif   ! active cell

    enddo      ! i
    enddo      ! j

  end subroutine assemble_stiffness_matrix

!****************************************************************************

  subroutine get_basis_function_derivatives_2d(nNodesPerElement,       &
                                               xNode,       yNode,     &
                                               dphi_dxr,    dphi_dyr,  &
                                               dphi_dx,     dphi_dy,   &
                                               detJ )

    !------------------------------------------------------------------
    ! Evaluate the x and y derivatives of the element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 2D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element
 
    real(dp), dimension(nNodesPerElement), intent(in) :: &
                xNode, yNode,        &! nodal coordinates
                dphi_dxr, dphi_dyr    ! derivatives of basis functions at quad pt
                                      !  wrt x and y in reference element

    real(dp), dimension(nNodesPerElement), intent(out) :: &
                dphi_dx , dphi_dy     ! derivatives of basis functions at quad pt
                                      !  wrt x and y in true coordinates  

    real(dp), intent(out) :: &
                detJ      ! determinant of Jacobian matrix

    real(dp), dimension(2,2) ::  &
                Jac,     &! Jacobian matrix
                Jinv      ! inverse Jacobian matrix

    integer :: n

    !------------------------------------------------------------------
    ! Compute the Jacobian for the transformation from the reference
    ! coordinates to the true coordinates:
    !
    !              |                                                 |
    !              | sum_n{dphi_n/dxr * xn}   sum_n{dphi_n/dxr * yn} |
    !   J(xr,yr) = |                                                 |
    !              | sum_n{dphi_n/dyr * xn}   sum_n{dphi_n/dyr * yn} |
    !              |                                                 |
    !
    ! where (xn,yn) are the true nodal coordinates,
    !       (xr,yr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

    do n = 1, nNodesPerElement
       Jac(1,1) = Jac(1,1) + dphi_dxr(n) * xNode(n)
       Jac(1,2) = Jac(1,2) + dphi_dxr(n) * yNode(n)
       Jac(2,1) = Jac(2,1) + dphi_dyr(n) * xNode(n)
       Jac(2,2) = Jac(2,2) + dphi_dyr(n) * xNode(n)
    enddo

    !------------------------------------------------------------------
    ! Compute the determinant and inverse of J
    !------------------------------------------------------------------

    detJ = Jac(1,1)*Jac(2,2) - Jac(1,2)*Jac(2,1)

    if (abs(detJ) > 0.d0) then
       Jinv(1,1) =  Jac(2,2)
       Jinv(1,2) = -Jac(1,2)
       Jinv(2,1) = -Jac(2,1)
       Jinv(2,2) =  Jac(1,1)
       Jinv(:,:) = Jinv(:,:) / detJ
    else
       !whl - do a proper abort here
       print*, 'stopping, det J = 0'
       stop
    endif

    !------------------------------------------------------------------
    ! Compute the contribution of this quadrature point to dphi/dx and dphi/dy
    ! for each basis function.
    !
    !   | dphi_n/dx |          | dphi_n/dxr |
    !   |           | = Jinv * |            | 
    !   | dphi_n/dy |          | dphi_n/dyr |
    !------------------------------------------------------------------

    dphi_dx(:) = 0.d0
    dphi_dy(:) = 0.d0

    do n = 1, nNodesPerElement
       dphi_dx(n) = dphi_dx(n) + Jinv(1,1)*dphi_dxr(n)  &
                               + Jinv(1,2)*dphi_dyr(n)
       dphi_dy(n) = dphi_dy(n) + Jinv(2,1)*dphi_dxr(n)  &
                               + Jinv(2,2)*dphi_dyr(n)
    enddo

  end subroutine get_basis_function_derivatives_2d

!****************************************************************************

!whl - Write a 3D version of this subroutine too.

  subroutine get_basis_function_derivatives_3d(nNodesPerElement,       &
                                               xNode,     yNode,    znode,     &
                                               dphi_dxr,  dphi_dyr, dphi_dzr,  &
                                               dphi_dx,   dphi_dy,  dphi_dz,   &
                                               detJ )

    !------------------------------------------------------------------
    ! Evaluate the x and y derivatives of the element basis functions
    ! at a particular quadrature point.
    !
    ! Also determine the Jacobian of the transformation between the
    ! reference element and the true element.
    ! 
    ! This subroutine should work for any 3D element with any number of nodes.
    !------------------------------------------------------------------

    integer, intent(in) :: nNodesPerElement   ! number of nodes per element
 
    real(dp), dimension(nNodesPerElement), intent(in) :: &
                xNode, yNode, zNode,         &! nodal coordinates
                dphi_dxr, dphi_dyr, dphi_dzr  ! derivatives of basis functions at quad pt
                                              !  wrt x and y in reference element

    real(dp), dimension(nNodesPerElement), intent(out) :: &
                dphi_dx, dphi_dy, dphi_dz     ! derivatives of basis functions at quad pt
                                              !  wrt x and y in true coordinates  

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
    !                 |                                                                        |
    !                 | sum_n{dphi_n/dxr * xn}  sum_n{dphi_n/dxr * yn}  sum_n{dphi_n/dxr * zn} |
    !   J(xr,yr,zr) = |                                                                        |
    !                 | sum_n{dphi_n/dyr * xn}  sum_n{dphi_n/dyr * yn}  sum_n(dphi_n/dyr * zn  |
    !                 |                                                                        |
    !                 | sum_n{dphi_n/dzr * xn}  sum_n{dphi_n/dzr * yn}  sum_n(dphi_n/dyz * zn  |
    !
    ! where (xn,yn,zn) are the true nodal coordinates,
    !       (xr,yr,zr) are the coordinates of the quad point in the reference element,
    !       and sum_n denotes a sum over nodes.
    !------------------------------------------------------------------

    Jac(:,:) = 0.d0

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

!whl - to do - bug check - Verify that J * Jinv = I

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

  end subroutine get_basis_function_derivatives_3d

!****************************************************************************

!whl - Write a 3D version of this subroutine for the BP equations.

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

  subroutine increment_element_matrix_first_order_3d(nNodesPerElement,             &
                                                     wqp,      detJ,               &
                                                     dphi_dx,  dphi_dy,  dphi_dz,  &
                                                     Kuu,      Kuv,                &
                                                     Kvu,      Kvv )

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

  end subroutine increment_element_matrix_first_order_3d

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

  subroutine glissade_velo_init_3d(nx, ny,   &
!!!                                nz,       &     !whl - not needed?
                                   dx, dy)

!whl - This is an old version.
!      Modify to be consistent with glissade_velo_init_2d
!      Separate out the parts that are independent of the velocity solver (2D or 3D).


  ! Initialize some data structures related to the mesh.
  ! These are similar to data structures in the MPAS framework, but for a grid
  !  that is structured in the horizontal.
  !
  ! Assumptions:
  ! (1) The mesh is logically rectangular in the horizontal with dimensions (nx, ny).
  ! (2) The number of cells is nx * ny.
  ! (3) The number of vertices is (nx-1)*(ny-1).  Hence each vertex is associated with
  !     four cells.  Interior cells are associated with four vertices, but border cells
  !     are associated with only 1 or 2 vertices.  Vertex (i,j) lies at the NE corner
  !     of cell (i,j).
  ! (4) The number of edges is (nx-1)*ny + (ny-1)*nx.  Each edge is associated with
  !     2 cells.  Interior cells are associated with 4 edges, but border cells are
  !     associated with only 2 or 3 edges.
  ! (5) Border cells (i = 1, i = nx, j = 1, j = ny) exist in the data structure, but
  !     they are not computational elements.  Any ice in the border cells is dynamically
  !     inactive.

    integer, intent(in) :: nx, ny       ! number of cells in x and y dimensions
!!!    integer, intent(in) :: nz            ! number of levels in z dimension
    real(dp), intent(in) :: dx, dy      ! grid cell dimensions

!whl - Remove later
    integer, parameter ::      &
       nVerticesOnCell = 8,    &
       nVerticesOnVertex = 27

    integer ::     &
       nCells,              &! number of cells on the mesh
       nVertices,           &! number of vertices
       nEdges                ! number of edges

    integer, dimension(:,:), allocatable ::    &
       gCellIndex,          &! assigns global cell index for each horizontal cell (i,j)
       gVertexIndex          ! assigns global vertex index for each horizontal vertex (i,j)
  
    integer, dimension(:,:), allocatable ::    &
       VerticesOnCell,      &! global indices of 4 vertices of a given cell (CCW from SW corner)
       CellsOnVertex,       &! global indices of 4 cells neighboring a given vertex (CCW from SW
       VerticesOnEdge,      &! global indices of 2 vertices of a given edge (S to N or W to E)
       EdgesOnCell,         &! global indices of 4 edges of a given cell (CCW from S edge)
       CellsOnEdge           ! global indices of 2 cells neighboring a given edge (S to N, or W to E)

    real(dp), dimension(:), allocatable ::   &
       VertexCoordX,        &! x coordinates of vertices
       VertexCoordY          ! x coordinates of vertices

    integer :: i, j, gc, gv 
!!!    integer :: k, n
    integer :: count, ewcount, nscount

!whl - Move or remove later
    ! trilinear basis set for 3D             
    ! N1 = (1-q1)*(1-q2)*(1-q3)             N4----N3
    ! N2 =    q1 *(1-q2)*(1-q3)             |     |    Lower layer
    ! N3 =    q1 * q2   *(1-q3)             |     |
    ! N4 = (1-q1)* q2   *(1-q3)             N1----N2
    !
    ! N5 = (1-q1)*(1-q2)* q3                N5----N8
    ! N6 =    q1 *(1-q2)* q3                |     |    Upper layer
    ! N7 =    q1 * q2   * q3                |     |
    ! N8 = (1-q1)* q2   * q3                N5----N6

    ! Count number of cells, vertices, and edges

    nCells = nx * ny
    nVertices = (nx-1) * (ny-1)
    nEdges = (nx-1)*ny + (ny-1)*nx

    ! Assign global index to each cell 

    allocate(gCellIndex(nx,ny))

    gCellIndex(:,:) = 0
    count = 0   
    do j = 1, ny   ! S to N
    do i = 1, nx   ! W to E
       count = count + 1
       gCellIndex(i,j) = count
    enddo
    enddo

    ! Assign global index to each vertex
    ! Assign (x,y) coordinates to each vertex 

    allocate(gVertexIndex(nx-1,ny-1))
    allocate(VertexCoordX(nVertices))
    allocate(VertexCoordY(nVertices))

    gVertexIndex(:,:) = 0
    count = 0   
    do j = 1, ny-1   ! S to N
    do i = 1, nx-1   ! W to E
       count = count + 1
       gVertexIndex(i,j) = count
       VertexCoordX(count) = (dx*i)
       VertexCoordY(count) = (dy*j)
    enddo
    enddo
    nVertices = count

    ! 4 vertices of each cell
    ! Vertices outside the domain are assigned index 0.

    allocate(VerticesOnCell(nCells,4))

    VerticesOnCell(:,:) = 0
    do j = 1, ny
    do i = 1, nx
       gc = gCellIndex(i,j)
       if (i > 1  .and. j > 1 ) VerticesOnCell(gc,1) = gVertexIndex(i-1,j-1)  ! SW
       if (i < nx .and. j > 1 ) VerticesOnCell(gc,2) = gVertexIndex(i  ,j-1)  ! SE
       if (i < nx .and. j < ny) VerticesOnCell(gc,3) = gVertexIndex(i  ,j  )  ! NE
       if (i > 1  .and. j < ny) VerticesOnCell(gc,4) = gVertexIndex(i-1,j  )  ! NW
    enddo
    enddo

    ! 4 cells neighboring each vertex

    allocate(CellsOnVertex(nVertices,4))

    CellsOnVertex(:,:) = 0
    do j = 1, ny-1
    do i = 1, nx-1
       gv = gVertexIndex(i,j)
       CellsOnVertex(gv,1) = gCellIndex(i  ,j  )  ! SW
       CellsOnVertex(gv,2) = gCellIndex(i+1,j  )  ! SE
       CellsOnVertex(gv,3) = gCellIndex(i+1,j+1)  ! NE
       CellsOnVertex(gv,4) = gCellIndex(i  ,j+1)  ! NW
    enddo
    enddo

    ! Edge data structures:
    ! 2 vertices on each edge, 2 cells on each edge

    allocate(VerticesOnEdge(nEdges, 2))
    allocate(CellsOnEdge(nEdges, 2))
    allocate(EdgesOnCell(nCells, 4))

    VerticesOnEdge(:,:) = 0
    CellsOnEdge(:,:) = 0
    EdgesOnCell(:,:) = 0

    ! Start with edges oriented parallel to E-W axis.
    ! These edges are numbered 1, 3, 5, 7, etc.
    ! Note: Edges along N and S boundaries of domain are not defined.

    ewcount = 0
    do j = 1, ny-1   
    do i = 1, nx
       ewcount = ewcount + 1
       count = 2*ewcount - 1

       gc = gCellIndex(i,j)  ! south of the edge
       CellsOnEdge(count,1) = gc
       EdgesOnCell(gc,3) = count  ! N edge of cell

       gc = gCellIndex(i,j+1)  ! north of the edge
       CellsOnEdge(count,2) = gc 
       EdgesOnCell(gc,1) = count  ! S edge of cell

       if (i > 1)  VerticesOnEdge(count,1) = gVertexIndex(i-1,j)
       if (i < nx) VerticesOnEdge(count,2) = gVertexIndex(i,j)

    enddo
    enddo

    ! Now the edges oriented parallel to N-S axis
    ! These edges are numbered 2, 4, 6, 8, etc.
    ! Note: Edges along E and W boundaries of domain are not defined.

    nscount = 0
    do j = 1, ny   
    do i = 1, nx-1
       nscount = nscount + 1
       count = 2*nscount

       gc = gCellIndex(i,j)  ! west of the edge
       CellsOnEdge(count,1) = gc
       EdgesOnCell(gc,2) = count  ! E edge of cell

       gc = gCellIndex(i+1,j)  ! east of the edge
       CellsOnEdge(count,2) = gc 
       EdgesOnCell(gc,4) = count  ! W edge of cell

       if (j > 1)  VerticesOnEdge(count,1) = gVertexIndex(i,j-1)
       if (j < ny) VerticesOnEdge(count,2) = gVertexIndex(i,j)

    enddo
    enddo

!whl - Write out and debug

   go to 100

   ! Cells
   do j = 1, ny
   do i = 1, nx
      gc = gCellIndex(i,j)
      write(6,*) ' '
      write(6,*) 'i, j, gcell:', i, j, gc
      write(6,*) 'Edges:',    EdgesOnCell(gc,:)
      write(6,*) 'Vertices:', VerticesOnCell(gc,:)
   enddo
   enddo

   ! Vertices
   do j = 1, ny-1
   do i = 1, nx-1
      gv = gVertexIndex(i,j)
      write(6,*) ' '    
      write(6,*) 'i, j, gvertex:', i, j, gv
      write(6,*) 'x, y coords:', VertexCoordX(gv), VertexCoordY(gv)
      write(6,*) 'Cells:', CellsOnVertex(gv,:)
   enddo
   enddo

   ! EW edges
   count = -1
   do j = 1, ny-1
   do i = 1, nx
      count = count + 2
      write(6,*) ' '
      write(6,*) 'nedge =', count
      write(6,*) 'Vertices:', VerticesOnEdge(count,:)
      write(6,*) 'Cells:', CellsOnEdge(count,:)
   enddo
   enddo

   ! NS edges
   count = 0
   do j = 1, ny
   do i = 1, nx-1
      count = count + 2
      write(6,*) ' '
      write(6,*) 'nedge =', count
      write(6,*) 'Vertices:', VerticesOnEdge(count,:)
      write(6,*) 'Cells:', CellsOnEdge(count,:)
   enddo
   enddo

100 return

  end subroutine glissade_velo_init_3d

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

          auu(nv1,i00) = auu(nv1,i00) + efvs(nc) * kuu(1,1) 
          auu(nv1,i0e) = auu(nv1,i0e) + efvs(nc) * kuu(1,2) 
          auu(nv1,ine) = auu(nv1,ine) + efvs(nc) * kuu(1,3) 
          auu(nv1,in0) = auu(nv1,in0) + efvs(nc) * kuu(1,4) 

          auu(nv2,i0w) = auu(nv2,i0w) + efvs(nc) * kuu(2,1) 
          auu(nv2,i00) = auu(nv2,i00) + efvs(nc) * kuu(2,2) 
          auu(nv2,in0) = auu(nv2,in0) + efvs(nc) * kuu(2,3) 
          auu(nv2,inw) = auu(nv2,inw) + efvs(nc) * kuu(2,4) 

          auu(nv3,isw) = auu(nv3,isw) + efvs(nc) * kuu(3,1) 
          auu(nv3,is0) = auu(nv3,is0) + efvs(nc) * kuu(3,2) 
          auu(nv3,i00) = auu(nv3,i00) + efvs(nc) * kuu(3,3) 
          auu(nv3,i0w) = auu(nv3,i0w) + efvs(nc) * kuu(3,4) 

          auu(nv4,is0) = auu(nv4,is0) + efvs(nc) * kuu(4,1) 
          auu(nv4,ise) = auu(nv4,ise) + efvs(nc) * kuu(4,2) 
          auu(nv4,i0e) = auu(nv4,i0e) + efvs(nc) * kuu(4,3) 
          auu(nv4,i00) = auu(nv4,i00) + efvs(nc) * kuu(4,4) 

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
