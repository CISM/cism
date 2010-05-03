! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  phaml_example.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"
module phaml_example
  use glide_types
  use glimmer_global
  use phaml
  use phaml_user_mod
  implicit none
  


!-------------------------------------------------------------------------
    subroutine phaml_init(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call user_init(get_ewn(cism_model),get_nsn(cism_model), &
        get_dew(cism_model),get_dns(cism_model))
    end subroutine phaml_init

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
    subroutine phaml_setup(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use phaml_support
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        
        !initialise phaml and the variables needed for solution
        call phaml_init(cism_model,phaml_solution)
        !generate the mesh file
        call make_poly_file(model)
        !set the mesh and create the phaml solution
        call phaml_create(phaml_solution,nproc=1,update_umod=.true., &
            triangle_files="mesh.1", &
            !spawn_form=DEBUG_SLAVE, &
            draw_grid_who=NO_ONE)
            
        !set the initial conditions array    
        uphaml = cism_model%phaml%uphaml
        
        !set the initial conditions
        !update usermod must be called at least twice before iconds
        !can be set with phaml_solve_pde
        call update_usermod(phaml_solution)
        !this sets the initial conditions
        call phaml_solve_pde(phaml_solution,            &
                     !max_vert=32000,                   &
                     task=SET_INITIAL,                  &
                     refterm=KEEP_NELEM,                &
                     max_eq=3200,                       &
                     error_estimator=INITIAL_CONDITION, &
                     print_header_who=NO_ONE,           &
                     print_trailer_who=NO_ONE,          &
                     degree=3,                          &
                     draw_grid_when=NEVER)
                     
        !necessary here???????????????
        !put the solution in the uphaml variable                 
        call phaml_getsolution(phaml_solution, cism_model%phaml%uphaml)           
    end subroutine phaml_lin_evolve

!-------------------------------------------------------------------------
    

!-------------------------------------------------------------------------
    subroutine phaml_lin_evolve(cism_model,phaml_solution)
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call phaml_copy_soln_to_old(phaml_solution)
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         task = SOLVE_ONLY,             &   !these two options make it
                         max_refsolveloop= 1,           &   !solve just one iteration
                         print_error_who=MASTER,        & 
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
                         
    end subroutine phaml_lin_evolve

!-------------------------------------------------------------------------

!-------------------------------------------------------------------------
    subroutine phaml_nonlin_evolve(cism_model,phaml_solution)
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call phaml_copy_soln_to_old(phaml_solution)
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         task = SOLVE_ONLY,             &   !these two options make it
                         max_refsolveloop= 1,           &   !solve just one iteration
                         print_error_who=MASTER,        & 
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
                         
    end subroutine phaml_nonlin_evolve

!-------------------------------------------------------------------------
!this creates the mesh and evolves, then returns the solution
!-------------------------------------------------------------------------
    subroutine phaml_evolve(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use phaml_support
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        !initialise phaml and the variables needed for solution
        call phaml_init(cism_model,phaml_solution)
        !generate the mesh file
        call make_poly_file(model)
        !set the mesh and create the phaml solution
        call phaml_create(phaml_solution,nproc=1,update_umod=.true., &
            triangle_files="mesh.1", &
            !spawn_form=DEBUG_SLAVE, &
            draw_grid_who=NO_ONE)
            
        !set the initial conditions array    
        uphaml = cism_model%phaml%uphaml
        
        !set the initial conditions
        !update usermod must be called at least twice before iconds
        !can be set with phaml_solve_pde
        call update_usermod(phaml_solution)
        
        !solve the pde
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         !task = SOLVE_ONLY,            &   !these two options make it
                         !max_refsolveloop= 1,          &   !solve just one iteration
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
                         
        !put the solution in the uphaml variable                 
        call phaml_getsolution(phaml_solution, cism_model%phaml%uphaml)
        !close the phaml session
        call phaml_close(phaml_solution)
    end subroutine phaml_evolve

!-------------------------------------------------------------------------
    subroutine phaml_getsolution(phaml_solution, uout)
        use phaml
        use phaml_user_mod
        implicit none
        type(phaml_solution_type) :: phaml_solution
        real(my_real), intent(out), dimension(:,:) :: uout
        real(my_real), allocatable, dimension(:) :: x,y,u
        
        !the globals are set in init and in phaml_user_mod
        allocate(x(gewn*gnsn))
        allocate(y(gewn*gnsn))
        allocate(u(gewn*gnsn))
        !this sets up the locations of our grid points to pass to phaml_evaluate  
        call get_xyarrays(x,y)
        !get the initial solution in u
        call phaml_evaluate(soln,x,y,u)
        !write solution to uphaml in the model
        call reshape_array_to_two(uout,u)
        !free memory
        deallocate(x)
        deallocate(y)
        deallocate(u)
    end subroutine phaml_getsolution(phaml_solution)
!-------------------------------------------------------------------------
    subroutine phaml_close(phaml_solution)
        use phaml
        use phaml_user_mod
        implicit none
        type(phaml_solution_type) :: phaml_solution
        
        !destroy phaml session
        call phaml_destroy(phaml_solution)
        call array_deinit()
    end subroutine phaml_close

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
       
    !----------------------------------------------------
    ! This file contains the user supplied external subroutines that define
    ! the PDE(s) to be solved, and other required external subroutines.
    !   pdecoefs bconds boundary_point boundary_npiece boundary_param iconds trues
    !   truexs trueys update_usermod phaml_integral_kernel
    !----------------------------------------------------
    
    !-------------------------------------------------------------------------
    subroutine pdecoefs(x,y,cxx,cxy,cyy,cx,cy,c,rs)
        !----------------------------------------------------
        ! This subroutine returns the coefficient and right hand side of the PDE
        ! at the point (x,y)
        !
        ! The PDE is
        !
        !    -( cxx(x,y) * u  )  -( cyy(x,y) * u  ) + c(x,y) * u = rs(x,y)
        !                   x  x                y  y
        !
        ! For eigenvalue problems, the right hand side is lambda * u * rs(x,y)
        !
        ! cxy, cx and cy are not currently used and should not be set.  They are
        ! for future expansion.
        !
        ! NOTE: BE CAREFUL TO GET THE SIGNS RIGHT
        ! e.g. cxx=cyy=1 means rs=-(uxx+uyy)
        !
        !----------------------------------------------------
        use phaml
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        real(my_real), intent(out) :: cxx(:,:),cxy(:,:),cyy(:,:),cx(:,:),cy(:,:), &
                                      c(:,:),rs(:)
        !----------------------------------------------------
        ! Begin executable code
        
        cxx(1,1) = 1.0_my_real
        cyy(1,1) = 1.0_my_real
        c(1,1) = 0.0_my_real
        rs(1) = 0.0_my_real
        
        cxy=0; cx=0; cy=0
    
    end subroutine pdecoefs
    
    !-------------------------------------------------------------------------
    subroutine bconds(x,y,bmark,itype,c,rs)
        !----------------------------------------------------
        ! This subroutine returns the boundary conditions at the point (x,y).
        !
        ! Each boundary condition is either
        !
        !    u = rs(x,y) or  u  + c(x,y)*u = rs(x,y)
        !                     n
        !
        ! itype indicates which type of condition applies to each point, using
        ! symbolic constants from module phaml.  It is DIRICHLET for the first
        ! condition, NATURAL for the second condition with c==0, and MIXED for
        ! the second condition with c/=0.
        !
        !----------------------------------------------------
        use phaml
        use phaml_user_mod
        use glide_mask
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: bmark
        integer, intent(out) :: itype(:)
        real(my_real), intent(out) :: c(:,:),rs(:)
        integer :: ew,ns
        !----------------------------------------------------
        ! Non-module procedures used are:
        
        interface
        
           function trues(x,y,comp,eigen) ! real (my_real)
           use phaml
           real (my_real), intent(in) :: x,y
           integer, intent(in) :: comp,eigen
           real (my_real) :: trues
           end function trues
        
        end interface
        
        !----------------------------------------------------
        ! Begin executable code
        ! Dirichlet boundary conditions
        itype = DIRICHLET    
        ew = getew(x)
        ns = getns(y)
        rs(1) = 0.0_my_real
        if (is_grounding_line(bmark) then
           rs(1) = 0.0_my_real
           if(uphaml(ew,ns) .gt. 0) then
               rs(1) = uphaml(ew,ns)
           end if
        !elseif (has_ice(bmark) .eqv. .true.) then
        !   rs(1) = 0.0_my_real
        endif
    
        c = 0.0_my_real
    
    end subroutine bconds
    
    !-------------------------------------------------------------------------
    function iconds(x,y,comp,eigen)
        !----------------------------------------------------
        ! This routine returns the initial condition for a time dependent problem.
        ! It can also be used for the initial guess for a nonlinear problem, or
        ! systems of equations solved by successive substitution.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        ! For problems where there are no initial conditions, it is a dummy.
        !----------------------------------------------------
        use phaml
        use phaml_user_mod
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real(my_real) :: ret_value
        real(my_real) :: iconds
        integer :: ew, ns
        
        !maps an x, y to nearest integer r,c value for array lookup
        ew = getew(x)
        ns = getns(y)
        
        !----------------------------------------------------
        ! Begin executable code
        !if (uphaml(ew,ns) .gt. 0) then
        !    write(*,*) 'iconds:'
        !    write(*,*) uphaml(ew,ns)
        !end if
        ret_value = uphaml(ew,ns)
        iconds = ret_value
    
    end function iconds
    
    !-------------------------------------------------------------------------
    function trues(x,y,comp,eigen) ! real (my_real)
        !----------------------------------------------------
        ! This is the true solution of the differential equation, if known.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        !----------------------------------------------------
        use phaml
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real (my_real) :: trues
        !----------------------------------------------------
        ! Begin executable code
        
        !trues = x**2 + y**2
        trues = huge(0.0_my_real) !return 0 since we don't know
    
    end function trues
    
    !-------------------------------------------------------------------------
    function truexs(x,y,comp,eigen) ! real (my_real)
        !----------------------------------------------------
        ! This is the x derivative of the true solution of the differential
        ! equation, if known.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        !----------------------------------------------------
        use phaml
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real (my_real) :: truexs
        !----------------------------------------------------
        ! Begin executable code
        
        truexs = huge(0.0_my_real)
    
    end function truexs
    
    !-------------------------------------------------------------------------
    function trueys(x,y,comp,eigen) ! real (my_real)
        !----------------------------------------------------
        ! This is the y derivative of the true solution of the differential
        ! equation, if known.
        ! comp,eigen is which solution to use from a coupled system of PDEs or multiple
        ! eigenvectors from an eigenvalue problem, and is ignored in this example.
        !----------------------------------------------------
        use phaml
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(in) :: x,y
        integer, intent(in) :: comp,eigen
        real (my_real) :: trueys
        !----------------------------------------------------
        ! Begin executable code
        
        trueys = huge(0.0_my_real)
    
    end function trueys
    
    !-------------------------------------------------------------------------
    subroutine boundary_point(ipiece,s,x,y)
        !----------------------------------------------------
        ! This routine defines the boundary of the domain.  It returns the point
        ! (x,y) with parameter s on piece ipiece.
        ! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
        ! data files.
        !----------------------------------------------------
        use phaml
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        integer, intent(in) :: ipiece
        real(my_real), intent(in) :: s
        real(my_real), intent(out) :: x,y
        !----------------------------------------------------
        ! Begin executable code
        ! Dummy version
        
        x = 0.0_my_real
        y = 0.0_my_real
    
    end subroutine boundary_point
    
    !-------------------------------------------------------------------------
    function boundary_npiece(hole)
        !----------------------------------------------------
        ! This routine gives the number of pieces in the boundary definition.
        ! If boundary_npiece <= 0 the domain is given by triangle data files.
        !----------------------------------------------------
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        integer, intent(in) :: hole
        integer :: boundary_npiece
        !----------------------------------------------------
        ! Begin executable code
        
        boundary_npiece = 0
    
    end function boundary_npiece
    
    !-------------------------------------------------------------------------
    subroutine boundary_param(start,finish)
        !----------------------------------------------------
        ! This routine gives the range of parameter values for each piece of the
        ! boundary.
        ! If boundary_npiece <= 0 it is a dummy and the domain is given by triangle
        ! data files.
        !----------------------------------------------------
        use phaml
        implicit none
        !----------------------------------------------------
        ! Dummy arguments
        real(my_real), intent(out) :: start(:), finish(:)
        !----------------------------------------------------
        ! Begin executable code
        
        ! Dummy version
        
        start = 0.0_my_real
        finish = 0.0_my_real
    
    end subroutine boundary_param
    
    
    !-------------------------------------------------------------------------
    function phaml_integral_kernel(kernel,x,y)
        use phaml
        integer, intent(in) :: kernel
        real(my_real), intent(in) :: x,y
        real(my_real) :: phaml_integral_kernel
        
        ! Identity function
        
        phaml_integral_kernel = 1.0
    
    end function phaml_integral_kernel
    
    !-------------------------------------------------------------------------
    function regularity(x,y)
        use phaml
        real(my_real), intent(in) :: x(3),y(3)
        real(my_real) :: regularity
        
        ! Dummy version, assume infinitely differentiable everywhere.
        
        regularity = huge(0.0_my_real)
    
    end function regularity


    subroutine update_usermod(phaml_solution) 
        !----------------------------------------------------
        ! This routine updates the module variables on the slave processes by
        ! sending them from the master process
        !----------------------------------------------------
        use phaml
        use phaml_user_mod
        implicit none
        type(phaml_solution_type), intent(in) :: phaml_solution
        integer :: iparam(5)
        real(my_real),allocatable,dimension(:) :: rparam
        logical, save :: first_call = .true.
        
        !the slaves need the globals passed first so that the allocate
        !for rparam has the values for the size.  This means after
        !phaml_create another call to update_usermod must be made
        if (first_call) then
            allocate(rparam(1))
            iparam(1) = gnsn
            iparam(2) = gewn
            iparam(3) = gdns
            iparam(4) = gdew
            iparam(5) = num_arrays
            call master_to_slaves(phaml_solution,iparam,rparam)
            gnsn = iparam(1)
            gewn = iparam(2)
            gdns = iparam(3)
            gdew = iparam(4)
            num_arrays = iparam(5)
            deallocate(rparam)
            call array_init()
            first_call = .false.
        else
            allocate(rparam(gnsn*gewn*num_arrays))
            iparam(1) = gnsn
            iparam(2) = gewn
            iparam(3) = gdns
            iparam(4) = gdew
            iparam(5) = num_arrays
            call reshape_array_to_one(uphaml,rparam)
            !Call the routine that performs the actual exchange.   
            call master_to_slaves(phaml_solution,iparam,rparam)
            gnsn = iparam(1)
            gewn = iparam(2)
            gdns = iparam(3)
            gdew = iparam(4)
            num_arrays = iparam(5)
            if(num_arrays .eq. 1) then
                call reshape_array_to_two(uphaml,rparam)
            end if
            !an else should be added here to handle passing more than one array
            !using the concat_arrays, split_arrays function in combination with the
            !combine and separate functions.
            deallocate(rparam)
        end if
    end subroutine update_usermod


!---------------------------------------------------------------------------
end module phaml_example
