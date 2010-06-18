! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  phaml_example.F90 - part of the Glimmer-CISM ice model      + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


module phaml_example
  use glide_types
  use glimmer_global
  use phaml
  use phaml_user_mod
  implicit none
  
contains
!need to specify which phaml module
!phaml_example is 1
!-------------------------------------------------------------------------
    subroutine phaml_init(cism_model,phaml_solution)
        use phaml
        use phaml_user_mod
        use glide_types
        implicit none
        type(phaml_solution_type) :: phaml_solution
        type(glide_global_type) :: cism_model
        
        call user_init(get_ewn(cism_model),get_nsn(cism_model), &
        get_dew(cism_model),get_dns(cism_model),1)
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
        call make_ice_poly_file(cism_model)
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
    end subroutine phaml_setup

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
        !call phaml_init(cism_model,phaml_solution)
        !generate the mesh file
        call make_ice_poly_file(cism_model)
        !set the mesh and create the phaml solution
        call phaml_create(phaml_solution,nproc=2,update_umod=.true., &
            triangle_files="mesh.1", &
            spawn_form=NORMAL_SPAWN, &!DEBUG_BOTH, &
            draw_grid_who=MASTER)!NO_ONE)
        !set the initial conditions array    
        uphaml = cism_model%phaml%uphaml
        
        !set the initial conditions
        !update usermod must be called at least twice before iconds
        !can be set with phaml_solve_pde
        call update_usermod(phaml_solution)
        write(*,*) 'phaml_solve'
        !solve the pde
        call phaml_solve_pde(phaml_solution,            &
                         draw_grid_when=PHASES, &!NEVER,          &
                         pause_at_start=.false.,        &
                         pause_after_phases=.false.,    &
                         task = SOLVE_ONLY,            &   !these two options make it
                         max_refsolveloop= 1,          &   !solve just one iteration
                         mg_cycles=20,                  &
                         mg_tol=MG_NO_TOL,              &
                         print_header_who=NO_ONE,       &
                         print_trailer_who=NO_ONE)
        write(*,*) 'donequ'                         
        !put the solution in the uphaml variable                 
        call phaml_getsolution(phaml_solution, cism_model%phaml%uphaml)
        !close the phaml session
        !call phaml_close(phaml_solution)
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
        call phaml_evaluate(phaml_solution,x,y,u)
        !write solution to uphaml in the model
        call reshape_array_to_two(uout,u)
        !free memory
        deallocate(x)
        deallocate(y)
        deallocate(u)
    end subroutine phaml_getsolution
!-------------------------------------------------------------------------
    subroutine phaml_close(phaml_solution)
        use phaml
        use phaml_user_mod
        implicit none
        type(phaml_solution_type) :: phaml_solution
        
        !destroy phaml session
        call phaml_destroy(phaml_solution)
        call user_close()
    end subroutine phaml_close

!-------------------------------------------------------------------------


!---------------------------------------------------------------------------
end module phaml_example
