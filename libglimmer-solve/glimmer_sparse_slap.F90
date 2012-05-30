module glimmer_sparse_slap
    !*FD This module builds on the glimmer_slap module to provide an easy
    !*FD interface to SLAP.  The SLAP interface is intended to be both
    !*FD usable and a guide to implementing other interfaces
    
    use glimmer_sparse_type
    use glimmer_global, only: dp, size_t
    use glimmer_log
    implicit none

    type slap_solver_workspace
        !*FD This type contains any working memory needed for the slap solver.
        !*FD It is used to store states between calls to the solver
        !*FD In the SLAP implementation, it is used to store the SLAP workspace
        !*FD This module must have this type, but its contents should be opaque
        !*FD to the user (e.g. client code should only manipulate the
        !*FD slap_solver_workspace as a whole and should never touch its members)
        real(kind=dp), dimension(:), pointer :: rwork => NULL()
        integer, dimension(:), pointer :: iwork => NULL()
        integer :: max_nelt !*FD Maximum number of nonzeroes allowed given the allocated workspace
    end type slap_solver_workspace

    type slap_solver_options
        !*FD This type holds options that are passed to the slap solver, such
        !*FD as preconditioner type, error tolerances, etc.  At a minimum, it
        !*FD must define the tolerance and maxiters field, as these will be
        !*FD common to any iterative slap linear solver.  Other options
        !*FD can be defined as necessary.
        !*FD
        !*FD Design note: the options are seperated from the workspace because
        !*FD one set of options could apply to multiple matrices, and the
        !*FD lifecycles for each could be different (a workspace need only
        !*FD exist as long as the matrix does, the options could persist
        !*FD throughout the entire program)
        integer :: itol !*FD Tolerance code, see SLAP documentation
        logical :: use_gmres !*FD Whether to use the GMRES method instead of Biconjugate Gradient
        integer :: gmres_saved_vectors !*FD How many vectors to save while performing GMRES iteration
        type(sparse_solver_options_base), pointer :: base => null() !*FD Pointer to basic options
    end type slap_solver_options

contains
    subroutine slap_default_options(opt, base)
        !*FD Populates a slap_solver_options (defined above) with default
        !*FD options.  This is necessary because different solvers may define
        !*FD different options beyond the required fields defined above.
        !*FD Filling them in this function allows client code to pick "good"
        !*FD values in a generic way.
        type(slap_solver_options), intent(out) :: opt
        type(sparse_solver_options_base), intent(in), target :: base
        opt%itol = 2
        opt%use_gmres = .false.
        opt%gmres_saved_vectors = 20
        opt%base => base
    end subroutine slap_default_options

    subroutine slap_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)
        !*FD Allocate solver workspace.  This needs to be done once
        !*FD (when the maximum number of nonzero entries is first known)
        !*FD This function need not be safe to call on already allocated memory
        !*FD
        !*FD Note that the max_nonzeros argument must be optional, and if
        !*FD it is not supplied the current number of nonzeroes must be used.
        type(sparse_matrix_type) :: matrix
        type(slap_solver_options) :: options
        type(slap_solver_workspace) :: workspace
        integer, optional :: max_nonzeros_arg
        integer :: max_nonzeros
        integer(kind=size_t) :: lenrw
        integer(kind=size_t) :: leniw
        
        if (present(max_nonzeros_arg)) then
            max_nonzeros = max_nonzeros_arg
        else
            max_nonzeros = matrix%nonzeros
        end if
        
        !Only allocate the memory if it hasn't been allocated or it needs
        !to grow
        if (.not. associated(workspace%rwork) .or. workspace%max_nelt < max_nonzeros) then
            !If memory is already allocated get rid of it
            if (associated(workspace%rwork)) then
                deallocate(workspace%rwork)
                deallocate(workspace%iwork)
            end if

            !Figure out how much memory to allocate.  These figures were derived
            !from the SLAP documentation.
            lenrw = 20*max_nonzeros 
            leniw = 20*max_nonzeros

            if (lenrw < 0 .or. leniw < 0) then
                call write_log("The amount of workspace memory that SLAP needs caused a numerical overflow.  " // &
                               "If you are not running on a 64-bit architecture, you will need to decrease" // & 
                               "the size of your data set.  If you are running a 64-bit architecture, try" // & 
                               "modifying size_t in glimmer_global to a larger size and recompiling Glimmer.", GM_FATAL)
            end if

            !write(*,*) "MAX NONZEROS",max_nonzeros
            !write(*,*) "ALLOCATING WORKSPACE",lenrw,leniw 
            allocate(workspace%rwork(lenrw))
            allocate(workspace%iwork(leniw))
            !Recored the number of nonzeros so we know whether to allocate more
            !memory in the future
            workspace%max_nelt = max_nonzeros
        end if
    end subroutine slap_allocate_workspace

    subroutine slap_solver_preprocess(matrix, options, workspace)
        !*FD Performs any preprocessing needed to be performed on the slap
        !*FD matrix.  Workspace must have already been allocated. 
        !*FD This function should be safe to call more than once.
        !*FD
        !*FD It is an error to call this function on a workspace without
        !*FD allocated memory
        !*FD
        !*FD In general slap_allocate_workspace should perform any actions
        !*FD that depend on the *size* of the slap matrix, and
        !*FD sprase_solver_preprocess should perform any actions that depend
        !*FD upon the *contents* of the slap matrix.
        type(sparse_matrix_type) :: matrix
        type(slap_solver_options) :: options
        type(slap_solver_workspace) :: workspace
    end subroutine slap_solver_preprocess

    function slap_solve(matrix, rhs, solution, options, workspace,err,niters, verbose)
        !*FD Solves the slap linear system, and reports status information.
        !*FD This function returns an error code that should be zero if the
        !*FD call succeeded and nonzero if it failed.  No additional error codes
        !*FD are defined.  Although this function reports back the final error
        !*FD and the number of iterations needed to converge, these should *not*
        !*FD be relied upon as not every slap linear solver may report them.
        type(sparse_matrix_type), intent(inout) :: matrix 
        !*FD Sparse matrix to solve.  This is inout because the slap solver
        !*FD may have to do some re-arranging of the matrix.
        
        real(kind=dp), dimension(:), intent(inout) :: rhs 
        !*FD Right hand side of the solution vector
        
        real(kind=dp), dimension(:), intent(inout) :: solution 
        !*FD Solution vector, containing an initial guess.

        type(slap_solver_options), intent(in) :: options
        !*FD Options such as convergence criteria
        
        type(slap_solver_workspace), intent(inout) :: workspace
        !*FD Internal solver workspace
        
        real(kind=dp), intent(out) :: err
        !*FD Final solution error
        
        integer, intent(out) :: niters
        !*FD Number of iterations required to reach the solution

        logical, optional, intent(in) :: verbose
        !*FD If present and true, this argument may cause diagnostic information
        !*FD to be printed by the solver (not every solver may implement this).
        
        integer :: slap_solve

        integer :: ierr !SLAP-provided error code
        integer :: iunit !Unit number to print verbose output to (6=stdout, 0=no output)
        integer :: isym !Whether matrix is symmetric
        
        logical :: allzeros
        integer :: i

        iunit = 0
        if (present(verbose)) then
            if(verbose) then
                iunit=6
                write(*,*) 'Tolerance=',options%base%tolerance
            end if
        end if

        if (matrix%symmetric) then
            isym = 1
        else
            isym = 0
        end if
       
        allzeros = .true.
        !Check if the RHS is zero; if it is, don't iterate!  The biconjugate
        !gradient method doesn't work in this case
        zero_check: do i = 1, size(rhs)
            if (rhs(i) /= 0) then
                allzeros = .false.
                exit zero_check
            end if
        end do zero_check


	!----------------------------------------------
	! RN_20091102: An example of calls to Trilinos solvers
	!#ifdef HAVE_TRILINOS
	!call helloworld()
	!#endif
	!----------------------------------------------


        if (allzeros) then
            err = 0
            ierr = 0
            niters = 0
            solution = 0
            call write_log("RHS of all zeros passed to BCG method; iteration not perfomred.", &
                           GM_WARNING, __FILE__, __LINE__)        
        else
            !Set up SLAP if it hasn't been already
            call slap_solver_preprocess(matrix, options, workspace)

            if (options%use_gmres) then
                call dslugm(matrix%order, rhs, solution, matrix%nonzeros, &
                            matrix%row, matrix%col, matrix%val, &
                            isym, options%gmres_saved_vectors, options%itol, &
                            options%base%tolerance, options%base%maxiters, &
                            niters, err, ierr, iunit, &
                            workspace%rwork, size(workspace%rwork), workspace%iwork, size(workspace%iwork))
            else
                call dslucs(matrix%order, rhs, solution, matrix%nonzeros, &
                            matrix%row, matrix%col, matrix%val, &
                            isym, options%itol, options%base%tolerance, options%base%maxiters,&
                            niters, err, ierr, iunit, &
                            workspace%rwork, size(workspace%rwork), workspace%iwork, size(workspace%iwork))
            end if
        end if
        slap_solve = ierr
    end function slap_solve

    subroutine slap_solver_postprocess(matrix, options, workspace)
        type(sparse_matrix_type) :: matrix
        type(slap_solver_options) :: options
        type(slap_solver_workspace) :: workspace
    end subroutine

    subroutine slap_destroy_workspace(matrix, options, workspace)
        !*FD Deallocates all working memory for the slap linear solver.
        !*FD This need *not* be safe to call of an unallocated workspace
        !*FD No slap solver should call this automatically.
        type(sparse_matrix_type) :: matrix
        type(slap_solver_options) :: options
        type(slap_solver_workspace) :: workspace
        !Deallocate all of the working memory
        deallocate(workspace%rwork)
        deallocate(workspace%iwork)
    end subroutine slap_destroy_workspace

    subroutine slap_interpret_error(error_code, error_string)
        !*FD takes an error code output from slap_solve and interprets it.
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
                tmp_error_string="Insufficient space allocated for WORK or IWORK"
            case (2)
                tmp_error_string="Method failed to converge in ITMAX steps"
            case (3)
                tmp_error_string="Error in user input.  Check input values of N, ITOL."
            case (4)
                tmp_error_string="User error tolerance set too tight." 
            case (5)
                tmp_error_string="Breakdown of the method detected.  (r0,r) approximately 0."
            case (6)
                tmp_error_string="Stagnation of the method detected. (r0, v) approximately 0."
            case (7)
                tmp_error_string="Incomplete factorization broke down and was fudged."
        end select


        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif
    end subroutine slap_interpret_error
end module glimmer_sparse_slap
