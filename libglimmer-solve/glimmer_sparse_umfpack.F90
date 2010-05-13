module glimmer_sparse_umfpack
    !*FD This module builds on the glimmer_umf module to provide an easy
    !*FD interface to SLAP.  The SLAP interface is intended to be both
    !*FD usable and a guide to implementing other interfaces
    
    use glimmer_sparse_type
    use glimmer_global, only: dp
    implicit none

    type umf_solver_workspace
        !*FD This type contains any working memory needed for the umf solver.
        !*FD It is used to store states between calls to the solver
        !*FD In the SLAP implementation, it is used to store the SLAP workspace
        !*FD This module must have this type, but its contents should be opaque
        !*FD to the user (e.g. client code should only manipulate the
        !*FD umf_solver_workspace as a whole and should never touch its members)

        !UMFPACK status arrays
        real(kind=dp),dimension(20) :: control
        real(kind=dp),dimension(90) :: info

        !UMFPACK pointers
        integer :: numeric
        integer :: symbolic

        logical :: alloc
    end type umf_solver_workspace

    type umf_solver_options
        !*FD This type holds options that are passed to the umf solver, such
        !*FD as preconditioner type, error tolerances, etc.  At a minimum, it
        !*FD must define the tolerance and maxiters field, as these will be
        !*FD common to any iterative umf linear solver.  Other options
        !*FD can be defined as necessary.
        !*FD
        !*FD Design note: the options are seperated from the workspace because
        !*FD one set of options could apply to multiple matrices, and the
        !*FD lifecycles for each could be different (a workspace need only
        !*FD exist as long as the matrix does, the options could persist
        !*FD throughout the entire program)
        real(kind=dp) :: tolerance !*FD Error tolerance
        integer :: maxiters !*FD Max iterations before giving up
        logical :: use_iterative_refinement
    end type umf_solver_options

contains
    subroutine check_umfpack()
#ifndef HAVE_UMFPACK
        write(*,*)"ERROR: Unsymmetric multifrontal method has been requested but CISM"
        write(*,*)"       was not compiled with UMFPACK support."
        stop
#endif
    end subroutine

    subroutine umf_default_options(opt)
        !*FD Populates a umf_solver_options (defined above) with default
        !*FD options.  This is necessary because different solvers may define
        !*FD different options beyond the required fields defined above.
        !*FD Filling them in this function allows client code to pick "good"
        !*FD values in a generic way.
        type(umf_solver_options), intent(out) :: opt
        opt%tolerance  = 1e-5
        opt%maxiters = 2000
        opt%use_iterative_refinement = .true.
    end subroutine umf_default_options

    subroutine umf_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)
        !*FD Allocate solver workspace.  This needs to be done once
        !*FD (when the maximum number of nonzero entries is first known)
        
        !*FD Note that the max_nonzeros argument must be optional, and if
        !*FD the current number of nonzeroes must be used.
        type(sparse_matrix_type) :: matrix
        type(umf_solver_options) :: options
        type(umf_solver_workspace) :: workspace
        integer, optional :: max_nonzeros_arg
    
        !*FD UMFPACK internally allocates all solver memory when we perform
        !*FD factorization.  Thus, the only thing we'll do here is to initialize
        !*FD the workspace
        workspace%alloc = .false.
    
        call check_umfpack()
    
    end subroutine umf_allocate_workspace

    subroutine umf_solver_preprocess(matrix, options, workspace)
        !*FD Performs any preprocessing needed to be performed on the umf
        !*FD matrix.  Workspace must have already been allocated. 
        !*FD This function should be safe to call more than once.
        !*FD 
        !*FD It is an error to call this function on a workspace without already
        !*FD allocated memory.
        !*FD
        !*FD In general umf_allocate_workspace should perform any actions
        !*FD that depend on the *size* of the umf matrix, and
        !*FD sprase_solver_preprocess should perform any actions that depend
        !*FD upon the *contents* of the umf matrix.
        type(sparse_matrix_type) :: matrix
        type(umf_solver_options) :: options
        type(umf_solver_workspace) :: workspace
        call check_umfpack()

#ifdef HAVE_UMFPACK
        !Check if we have previously allocated UMFPACK storage
        if (workspace%alloc) then
            call umf_destroy_workspace(matrix, options, workspace)
        end if

        !Set default control parameters
        call umf4def(workspace%control)

        !Convert to the UMFPACK format.
        !Move from triad format to column format
        call to_column_format(matrix)

        !Sort by row within each column in the column format.
        !to_column_format does not necessarily do this
        call sort_column_format(matrix)

        !Move from 1-based to 0-based indexing
        matrix%row = matrix%row - 1
        matrix%col = matrix%col - 1

        !Pre-order and symbolic analysis
        call umf4sym(matrix%order, matrix%order, matrix%col, matrix%row, matrix%val,  &
                     workspace%symbolic, workspace%control, workspace%info)
        if (workspace%info(1) < 0) then
            write(*,*)"Error occurred in umf4sym: ", workspace%info(1)
            call umf_interpret_error(int(workspace%info(1)))
            write(*,*) matrix%row
            write(*,*) matrix%col
            write(*,*) matrix%val
            stop
        end if

        !Numeric factorization
        call umf4num(matrix%col, matrix%row, matrix%val, &
                     workspace%symbolic, workspace%numeric, workspace%control, workspace%info)
        if (workspace%info(1) < 0) then
            write(*,*)"Error occured in umf4num: ", workspace%info(1)
            call umf_interpret_error(int(workspace%info(1)))
        end if

        workspace%alloc=.true.
#endif
    end subroutine umf_solver_preprocess

    function umf_solve(matrix, rhs, solution, options, workspace,err,niters, verbose)
        !*FD Solves the umf linear system, and reports status information.
        !*FD This function returns an error code that should be zero if the
        !*FD call succeeded and nonzero if it failed.  No additional error codes
        !*FD are defined.  Although this function reports back the final error
        !*FD and the number of iterations needed to converge, these should *not*
        !*FD be relied upon as not every umf linear solver may report them.

        type(sparse_matrix_type), intent(inout) :: matrix 
        !*FD Sparse matrix to solve.  This is inout because the umf solver
        !*FD may have to do some re-arranging of the matrix.
        
        real(kind=dp), dimension(:), intent(in) :: rhs 
        !*FD Right hand side of the solution vector
        
        real(kind=dp), dimension(:), intent(inout) :: solution 
        !*FD Solution vector, containing an initial guess.

        type(umf_solver_options), intent(in) :: options
        !*FD Options such as convergence criteria
        
        type(umf_solver_workspace), intent(inout) :: workspace
        !*FD Internal solver workspace
        
        real(kind=dp), intent(out) :: err
        !*FD Final solution error
        
        integer, intent(out) :: niters
        !*FD Number of iterations required to reach the solution

        logical, optional, intent(in) :: verbose
        !*FD If present and true, this argument may cause diagnostic information
        !*FD to be printed by the solver (not every solver may implement this).
        
        integer :: umf_solve

        integer :: iunit !Unit number to print verbose output to (6=stdout, 0=no output)
        integer :: isym !Whether matrix is symmetric
        integer :: sys

        sys=0
       
        call check_umfpack()
#ifdef HAVE_UMFPACK
        
        call umf4solr(sys, matrix%col, matrix%row, matrix%val, solution, rhs, &
                     workspace%numeric, workspace%control, workspace%info)
#endif
        umf_solve = workspace%info(1) 
    
        !umfpack is a direct method, so the iters and err returns mean nothing
        err = 0
        niters = 0
    end function umf_solve

    subroutine umf_solver_postprocess(matrix, options, workspace)
        type(sparse_matrix_type) :: matrix
        type(umf_solver_options) :: options
        type(umf_solver_workspace) :: workspace
    end subroutine

    subroutine umf_destroy_workspace(matrix, options, workspace)
        !*FD Deallocates all working memory for the umf linear solver.
        !*FD This need *not* be safe to call of an unallocated workspace
        !*FD No umf solver should call this automatically.
        type(sparse_matrix_type) :: matrix
        type(umf_solver_options) :: options
        type(umf_solver_workspace) :: workspace
        !Deallocate all of the working memory
        !Free the Umfpack symbolic analysis
#ifdef HAVE_UMFPACK
        call umf4fsym(workspace%symbolic)
        !Free the Umfpack numeric factorization
        call umf4fnum(workspace%numeric)
#endif
        workspace%alloc = .false.
    end subroutine umf_destroy_workspace

    subroutine umf_interpret_error(error_code, error_string)
        !*FD takes an error code output from umf_solve and interprets it.
        !*FD error_string must be an optional argument.
        !*FD If it is not provided, the error is printed to standard out
        !*FD instead of being put in the string
        integer :: error_code
        character(*), optional, intent(out) :: error_string
        character(256) :: tmp_error_string
        
        select case (error_code)
            case (0)
                tmp_error_string="All went well"
            case (1)
                tmp_error_string="Matrix is singular (This is a warning, but we barf anyway if it happens)"
            case (-1)
                tmp_error_string="Out of memory"
            case (-3)
                tmp_error_string="Invalid Numeric object"
            case (-4)
                tmp_error_string="Invalid Symbolic object" 
            case (-5)
                tmp_error_string="Argument missing"
            case (-6)
                tmp_error_string="n (matrix order given) is nonpositive"
            case (-8)
                tmp_error_string="Invalid umf matrix format"
            case (-11)
                tmp_error_string="Different pattern"
            case(-13)
                tmp_error_string="Invalid system"
            case(-15)
                tmp_error_string="Invalid permutation"
            case(-911)
                tmp_error_string="Internal error, contact the author of UMFPACK"
            case(-17)
                tmp_error_string="File IO error"
            case default
                tmp_error_string="Unrecognized error code"
        end select


        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif
    end subroutine umf_interpret_error
end module glimmer_sparse_umfpack
