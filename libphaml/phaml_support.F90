! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  phaml_support.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"
module phaml_support
!    use phaml
!    use phaml_user_mod
!    implicit none

    !------------------------------------------------------------------------------
    !SUBROUTINE: is_ice_edge
    !ARGUMENTS: model (glimmer), ew, ns
    !DESCRIPTION:
    !   ew, ns is a point in the model. The function checks the four adjacent
    !   points on the mask.  if (ew,ns) has ice, but one of the adjacent points
    !   does not have ice, then the function returns true.  This is a way to find 
    !   a pseudo grounding line/ice edge.
    !------------------------------------------------------------------------------
    
    subroutine is_ice_edge(model,ew,ns,ret)
        use glide_types
        use glide_mask
        use phaml_user_mod
        implicit none
        
        type(glide_global_type) :: model
        integer, intent(in) :: ew,ns
        logical, intent(out) :: ret
        logical :: m,m1,m2,m3,m4
        
        !m (ew,ns)
        m = .false.
        if(has_ice(model%geometry%thkmask(ew,ns)) .eqv. .true.) then
            m = .true.
        end if
        
        !m1 (ew,ns+1)
        m1 = .true.
        if(ns .eq. gnsn) then
            m1 = .false.
        else if(has_ice(model%geometry%thkmask(ew,ns+1)) .eqv. .false.) then
            m1 = .false.
        end if
        
        !m2 (ew,ns-1)
        m2 = .true.
        if(ns .eq. 1) then
            m2 = .false.
        else if(has_ice(model%geometry%thkmask(ew,ns-1)) .eqv. .false.) then
            m2 = .false.
        end if
        
        !m3 (ew+1,ns)
        m3 = .true.
        if(ew .eq. gewn) then
            m3 = .false.
        else if(has_ice(model%geometry%thkmask(ew+1,ns)) .eqv. .false.) then
            m3 = .false.
        end if
        
        !m4 (ew-1,ns)        
        m4 = .true.
        if(ew .eq. 1) then
            m4 = .false.
        else if(has_ice(model%geometry%thkmask(ew-1,ns)) .eqv. .false.) then
            m4 = .false.
        end if
        ret = .false.
        if(m .eqv. .true.) then
            if((m1 .eqv. .false.) .or. &
                (m2 .eqv. .false.) .or. &
                (m3 .eqv. .false.) .or. &
                (m4 .eqv. .false.)) then
                ret = .true.    
            end if
        end if
        !is_ice_edge = ret
    end subroutine is_ice_edge
    !------------------------------------------------------------------------------
    !SUBROUTINE: get_bmark
    !ARGUMENTS: model (glimmer), ew1, ns1, ew2, ns2, bmark
    !DESCRIPTION:
    !   ew1, ns1 is a point in the model, and ew2,ns2 is an adjacent point such 
    !   that there is an edge between the two.  This function determines if that
    !   edge is on the grounding line or within the boundaries and sets the bmark
    !   to the appropriate value.
    !------------------------------------------------------------------------------
    subroutine get_bmark(model,ew1,ns1,ew2,ns2,bmark)
        use glide_types
        use glide_mask
        implicit none
        
        type(glide_global_type) :: model
        integer, intent(in) :: ns1,ew1,ns2,ew2
        integer, intent(out) :: bmark
        logical :: ret1, ret2
        ret1 = .false.
        ret2 = .false.
        bmark=0 ! default value
        if ((has_ice(model%geometry%thkmask(ew1,ns1)) .eqv. .true.) .and. &
            (has_ice(model%geometry%thkmask(ew2,ns2)) .eqv. .true.)) then 
            bmark=model%geometry%thkmask(ew1,ns1) !3
        end if
        if ((is_grounding_line(model%geometry%thkmask(ew1,ns1)) .eqv. .true.) .and. &
            (is_grounding_line(model%geometry%thkmask(ew2,ns2)) .eqv. .true.)) then 
            bmark=model%geometry%thkmask(ew1,ns1)!5
        end if
        call is_ice_edge(model,ew1,ns1,ret1)
        call is_ice_edge(model,ew2,ns2,ret2)
        if (ret1 .and. ret2) then  
            bmark=glide_mask_grounding_line !model%geometry%thkmask(ew1,ns1)
        end if
        !get_bmark = bmark
        !return
    end subroutine get_bmark
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: make_poly_file
    !ARGUMENTS: model (glimmer)
    !DESCRIPTION:
    ! This subroutine creates a node file for an nx X ny grid on the rectangular
    ! domain (xmin,xmax) X (ymin,ymax), and runs triangle to create the
    ! required .node, .ele, .neigh and .edge files.
    !
    ! bmark is set to be
    !
    ! Running triangle uses the common extension "call system", which is
    ! compiler dependent.  That call may need to be changed, or may not be
    ! supported at all, on some compilers.
    !------------------------------------------------------------------------------
    subroutine make_poly_file(model)
        use glide
        use glide_mask
        implicit none
        
        integer :: nx, ny, node1ew,node1ns,node2ew,node2ns
        real :: dx,dy   
        integer :: i, j, count,bmark,status,tot_nodes
        logical :: ret
        type(glide_global_type) :: model
        !node_table addressed x,y, stores node number  and is -1 if not included, also bmark
        integer, dimension(:,:,:), allocatable :: node_table
        integer, dimension(:,:), allocatable :: edge_table   
        character string*2
    
        
        !grid variables
        nx = get_ewn(model)
        ny = get_nsn(model)
        dx = get_dew(model)
        dy = get_dns(model)!model%numerics%dns
        
        allocate(node_table(ny,nx,2))
        !declares the max number of edges needed if grid was square 1 is x, 2 is y
        count = max(nx,ny) 
        allocate(edge_table((4*count-2)*(count-1),3))
        
        !write a .poly file with the grid of nodes, boundary edges, and bmark
        
        open(unit=21,file="mesh.poly",status="replace")
        !count the nodes
        tot_nodes = 0
        count = 0
        do i=0,nx-1
            do j=0,ny-1
                bmark=0
                node1ew = j+1
                node1ns = i+1
                if (has_ice(model%geometry%thkmask(node1ew,node1ns)) .eqv. .true.) then
                    bmark=model%geometry%thkmask(node1ew,node1ns)!3
                end if
                if (is_grounding_line(model%geometry%thkmask(node1ew,node1ns)) .eqv. .true.) then 
                    bmark=model%geometry%thkmask(node1ew,node1ns)!5
                end if
                ret = .false.
                call is_ice_edge(model,node1ew,node1ns,ret)
                if (ret .eqv. .true.) then 
                    bmark=glide_mask_grounding_line!model%geometry%thkmask(node1ew,node1ns)!5
                end if
                if (bmark > 0) then
                    count = count + 1
                    node_table(node1ns, node1ew,1) = count
                    node_table(node1ns, node1ew,2) = bmark
                    tot_nodes = tot_nodes + 1
                else
                    node_table(node1ns, node1ew,1) = -1
                end if
                !write(21,*) count, xvals(node1ns), yvals(node1ew), bmark
            end do
        end do
        write(*,*) 'Writing Vertices'
        write(21,*) tot_nodes, 2, 0, 1
        do i=1,nx
            do j=1,ny
                if (node_table(j,i,1) .ge. 0) then
                    write(21,*) node_table(j,i,1), dx*(i-1), dy*(j-1), node_table(j,i,2)
                end if
            end do
        end do
        
        write(*,*) 'Writing Edges'
        count = 0
        !write edges connecting columns (vertical)
        do i=0,nx-1
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 2
                node2ew = i + 1
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                    count = count + 1
                    !write(22,*) count,node_table(node1ns,node1ew,1),node_table(node2ns,node2ew,1),bmark
                    edge_table(count,1) = node_table(node1ns,node1ew,1)
                    edge_table(count,2) = node_table(node2ns,node2ew,1)
                    edge_table(count,3) = bmark
                end if
            end do
        end do
        !write edges connecting rows (horizontal)
        do i=0,nx-2
            do j=0,ny-1
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 1
                node2ew = i + 2
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                    count = count + 1
                    !write(22,*) count,node_table(node1ns,node1ew,1),node_table(node2ns,node2ew,1),bmark
                    edge_table(count,1) = node_table(node1ns,node1ew,1)
                    edge_table(count,2) = node_table(node2ns,node2ew,1)
                    edge_table(count,3) = bmark
                end if
            end do
        end do
        !write forward cross edges    
        do i=0,nx-2
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 1
                node2ns = j + 2
                node2ew = i + 2
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                    count = count + 1
                    !write(22,*) count,node_table(node1ns,node1ew,1),node_table(node2ns,node2ew,1),bmark
                    edge_table(count,1) = node_table(node1ns,node1ew,1)
                    edge_table(count,2) = node_table(node2ns,node2ew,1)
                    edge_table(count,3) = bmark
                end if
            end do
        end do    
        !write backwards cross edges    
        do i=0,nx-2
            do j=0,ny-2
                bmark=0
                node1ns = j + 1
                node1ew = i + 2
                node2ns = j + 2
                node2ew = i + 1
                if ((node_table(node1ns,node1ew,1).ge. 0) .and. &
                    (node_table(node2ns,node2ew,1).ge. 0)) then
                    if ((node_table(j+1,i+1,1).lt. 0) .or. &
                        (node_table(j+2,i+2,1).lt. 0)) then
                        call get_bmark(model,node1ew,node1ns,node2ew,node2ns,bmark)
                        count = count + 1
                        !write(22,*) count,node_table(node1ns,node1ew,1),node_table(node2ns,node2ew,1),bmark
                        edge_table(count,1) = node_table(node1ns,node1ew,1)
                        edge_table(count,2) = node_table(node2ns,node2ew,1)
                        edge_table(count,3) = bmark
                    end if
                end if
            end do
        end do
        write(21,*) count, 1
        do i=1,count
            write(21,*) i, edge_table(i,1), edge_table(i,2), edge_table(i,3)
        end do
        write(21,*) 0
        
        close(unit=21)
        
        !deallocate the dynamic arrays
        deallocate( node_table )
        deallocate( edge_table )
        ! run triangle
        ! NOTE: THIS IS COMPILER DEPENDENT.  You may need to change this statement.
        status = System("triangle -neqj mesh.poly")
    
    end subroutine make_poly_file


!---------------------------------------------------------------------------
end module phaml_support





module phaml_user_mod
    !----------------------------------------------------
    ! This module contains user global data, and functions for formatting
    !
    !----------------------------------------------------
    use global
    use phaml
    implicit none
    
    integer, save :: gnsn, gewn, gdns, gdew, num_arrays
    real(my_real), save, allocatable, dimension(:,:) :: uphaml
    !real(my_real), save, allocatable, dimension(:) :: uphaml_one
contains
    subroutine user_init(ewn,nsn,dew,dns)
        implicit none
        integer, intent(in) :: nsn,ewn
        real(my_real), intent(in) :: dns,dew
        gnsn = nsn  !ny
        gewn = ewn  !nx
        gdns = dns  !dy
        gdew = dew  !dx
        !this defines how many arrays are passed in update_usermod
        num_arrays = 1
        !allocate(uphaml(gewn, gnsn))
        !allocate(uphaml_one(gewn*gnsn))
    end subroutine user_init
    
    subroutine array_init()
        !if (size(uphaml) .le. 0) then
            allocate(uphaml(gewn, gnsn))
        !end if
    end subroutine array_init
    
    subroutine array_deinit()
        if (size(uphaml) > 0) then
            deallocate(uphaml)
            !deallocate(uphaml_one)
        end if
    end subroutine array_deinit
    !-------------------------------------------------------------------------
    subroutine concat_arrays(array1, array2, retarray)
        implicit none
        real(my_real), intent(in),dimension(:) :: array1,array2
        real(my_real), intent(out),dimension(:) :: retarray
        integer :: i,maxind
        
        maxind = size(array1) + size(array2)
        do i=1, maxind
            retarray(i) = array1(i)
            retarray(maxind + i) = array2(i)
        end do
    end subroutine concat_arrays
    !-------------------------------------------------------------------------        
    subroutine split_arrays(array1, array2, retarray, maxind)
        implicit none
        real(my_real), intent(out),dimension(:) :: array1,array2
        real(my_real), intent(in),dimension(:) :: retarray
        integer, intent(in) :: maxind
        integer :: i
        
        do i=1, maxind
            array1(i) = retarray(i)
            array2(i) = retarray(maxind + i)
        end do
    end subroutine split_arrays
    !-------------------------------------------------------------------------    
    subroutine reshape_array_to_one(twod,oned)
        implicit none
        real(my_real), intent(in),dimension(:,:) :: twod
        real(my_real), intent(out),dimension(:) :: oned
        integer :: r,c
        
         do r=0,gnsn-1
            do c=0,gewn-1
               oned(gewn*r+c+1) = twod(c+1,r+1)
            end do
        end do
    end subroutine reshape_array_to_one
    !-------------------------------------------------------------------------    
    subroutine reshape_array_to_two(twod,oned)
        implicit none
        real(my_real), intent(out),dimension(:,:) :: twod
        real(my_real), intent(in),dimension(:) :: oned
        integer :: count,r,c
        
        do count=0,size(oned) - 1
            r = count/gewn + 1
            c = int(Mod(count,gewn)) + 1
            twod(c,r) = real(oned(count+1))        
        end do
    end subroutine reshape_array_to_two
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: get_xyarrays
    !ARGUMENTS: model (glimmer), x, y 
    !DESCRIPTION:
    !------------------------------------------------------------------------------
    subroutine get_xyarrays(x,y)
        implicit none
        real(my_real), intent(inout) :: x(:),y(:)
        integer :: count,r,c

        count = 1
        !populate the arrays with the xy coordinates
        do r=0,gnsn-1
            do c=0,gewn-1
               x(count) = REAL(c*gdew)
               y(count) = REAL(r*gdns)
               count = count + 1
            end do
        end do
    end subroutine get_xyarrays
    
    
    !these two functions get you the array indices to 
    !access the 2d grid given the x,y coordinates in the domain
    integer function getew(x)
        implicit none
        real(my_real), intent(in) :: x
        getew = int((x/gdew)) + 1
    end function getew
    
    integer function getns(y)
        implicit none
        real(my_real), intent(in) :: y
        getns = int((y/gdns)) + 1
    end function getns
    
end module phaml_user_mod


