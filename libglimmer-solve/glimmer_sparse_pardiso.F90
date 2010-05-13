#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_sparse_pardiso
    !*FD This module builds on the glimmer_sparse module to provide an 'easy'
    !*FD interface to PARDISO.
    
    use glimmer_sparse_type
    use glimmer_global
#ifdef HAVE_ISO_C_BINDING
    use ISO_C_BINDING, only: C_INTPTR_T
#endif
#ifdef _OPENMP
    use OMP_LIB
#endif
    implicit none

    type pardiso_solver_workspace
        !*FD Memory addresses used internally by PARDISO
#ifdef HAVE_ISO_C_BINDING
        integer(kind=C_INTPTR_T), dimension(64) :: pt
#else
        integer(kind=size_t), dimension(64) :: pt
#endif
       !*FD Error messages from pardisoinit: 0 no error, -10 no license, -11
        !*FD licence expired (?), -12 wrong user or hostname.
        integer(kind=4) :: error
        !*FD Work array for multi-recursive solver only
        real(kind=dp) :: dparm(64)
    end type pardiso_solver_workspace

    type pardiso_solver_options
        !*FD Matrix type, we are generally real, non-symmetric (11)
        integer(kind=4) :: mtype
        !*FD Solver, sparse direct is 0. If we ever have symmetric indefinate
        !*FD matrices, we might try 1; 'multi-recusive iterative'
        integer(kind=4) :: solver
        !*FD Storage for passing parameters to and getting info from PARDISO
        integer(kind=4), dimension(64) :: iparm
        ! Tolerance only needed in multirecursive solver
        real(kind=dp) :: tolerance
        ! Tell PARDISO to use the hybrid iterative/direct CGS method
        logical :: iterative
        ! Number of processors
        integer :: processors
    end type pardiso_solver_options

contains
    subroutine pardiso_default_options(opt)
        !*FD Populates a sparse_solver_options (defined above) with default
        !*FD options.  This is necessary because different solvers may define
        !*FD different options beyond the required fields defined above.
        !*FD Filling them in this function allows client code to pick "good"
        !*FD values in a generic way.
        implicit none
        type(pardiso_solver_options), intent(inout) :: opt

        ! TODO: bring these in from configuration
        opt%mtype = 11
        opt%solver = 0
        opt%tolerance = 1.d-9
        opt%iterative = .FALSE.  ! true for iterative/hybrid method
#ifdef _OPENMP
	opt%processors = OMP_GET_MAX_THREADS()
#else
        opt%processors = 1  ! Number of processors
#endif
    end subroutine pardiso_default_options

    subroutine pardiso_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)
        !*FD Allocate solver workspace.  This needs to be done once
        !*FD (when the maximum number of nonzero entries is first known)
        !*FD Note that the max_nonzeros argument must be optional, and if
        !*FD the current number of nonzeroes must be used.
        implicit none
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options),intent(inout) :: options
        type(pardiso_solver_workspace),intent(inout) :: workspace
        integer, optional :: max_nonzeros_arg
#ifdef HAVE_PARDISO
        !...Check license of the solver and initialize the solver
        call pardisoinit(workspace%pt, options%mtype, options%solver,&
                         options%iparm, workspace%dparm, workspace%error)
#endif
        ! These are options passed to PARDISO.
        ! However, they are reset by pardisoinit, so I set them here.
        ! This array is also used to return values. 
        if (options%iterative) then
            options%iparm(1) = 1   ! No defaults (1)
            options%iparm(2) = 2   ! Use Metis reordering (2)
            options%iparm(3) = options%processors   ! Number of processors
            ! LU CGS iteration (1) to 10^-6 (61)
            options%iparm(4) = 0
            !options%iparm(4) = int(-log10(options%tolerance) * 10 + 1)
            options%iparm(5) = 0   ! 0 Do not use permution 
            options%iparm(6) = 0   ! Write solution to seperate vector, not RHS
            options%iparm(8) = 0   ! Max iterative refinement steps
            options%iparm(10) = 13 ! eps pivot, 13 for non-symmetric, 8 symmetric
            options%iparm(11) = 1  ! Use non-symmetric scaling vector (1)
            options%iparm(12) = 0  ! Do not transpose matrix (0)
            options%iparm(13) = 2  ! Greater accuracy for indefinate matrix (2)
            options%iparm(18) = -1 ! Determine the number of non-zeros in LU (-1)
            options%iparm(19) = -1 ! Determine the Mflops for LU fact. (-1)
            options%iparm(21) = 1  ! Pivoting  1x1 and 2x2 Bunch-Kaufman
            options%iparm(24) = 1  ! Parallel numerical factor. (two level=1)
            options%iparm(25) = 1  ! Parallel solve (1)
            options%iparm(26) = 0  ! Splitting of basckward/forward solve (0)
            options%iparm(28) = 1  ! Parallel Metis (1)
            options%iparm(29) = 0  ! 64 bit accuracy (0)
            options%iparm(30) = 0  ! Default supernodes (0)
            options%iparm(31) = 0  ! No partial solutions using permute (0)
            options%iparm(32) = 0  ! Use sparse direct solver (0)
            options%iparm(33) = 0  ! Do not compute determinite (0) 
            options%iparm(34) = 0  ! Do not require identical parallel results(0)
        else
            options%iparm(3) = options%processors   ! Number of processors
        endif
    end subroutine pardiso_allocate_workspace

    subroutine pardiso_solver_preprocess(matrix, options, workspace)
        type(sparse_matrix_type),intent(inout) :: matrix
        type(pardiso_solver_options),intent(in) :: options
        type(pardiso_solver_workspace),intent(in) :: workspace

        ! To compressed row method is not in place, need to create storage
        integer,dimension(matrix%order+1) :: iao
        integer,dimension(matrix%nonzeros) :: jao
        real(dp),dimension(matrix%nonzeros) :: ao

        ! PARDISO is looking for a compressed row format for the matrix.
        ! Covervt triad to compressed row format
        call coocsr(matrix%order,matrix%nonzeros,matrix%val,matrix%row,&
                    matrix%col,ao,jao,iao)

        ! Place new sparse matrix into appropriate fields.
        matrix%row(1:matrix%order+1) = iao(:)
        matrix%col(1:matrix%nonzeros) = jao(:)
        matrix%val(1:matrix%nonzeros) = ao(:)

        ! Sort the column entries, needed by PARDISO
        call sort_row_format(matrix)
    end subroutine pardiso_solver_preprocess

    function pardiso_solve(matrix, rhs, solution,options ,workspace ,err ,niters ,verbose)
        implicit none
        type(sparse_matrix_type), intent(in) :: matrix 
        real(kind=dp), dimension(:), intent(in) :: rhs 
        real(kind=dp), dimension(:), intent(inout) :: solution 

        type(pardiso_solver_options), intent(in) :: options
        type(pardiso_solver_workspace), intent(inout) :: workspace
        integer, intent(out) :: niters
        real(kind=dp),intent(out) :: err
        logical, optional, intent(in) :: verbose
        integer :: pardiso_solve
        integer :: phase,mesglvl

        if (verbose) then 
           mesglvl = 6
        else
           mesglvl = 0
        end if

#ifdef HAVE_PARDISO
        ! This is a hack, needed to check the state of the solver, seems
        ! to work, but the behavior of workspace%pt is undocumented.
        if (workspace%pt(1) == 0) then
            ! Symbolic factorization only done on the first call to solver
            phase = 11
            call pardiso(workspace%pt, 1, 1, options%mtype, phase, matrix%order,&
                         matrix%val(1:matrix%nonzeros),&
                         matrix%row(1:matrix%order+1),&
                         matrix%col(1:matrix%nonzeros),&
                         1, 1, options%iparm, mesglvl, rhs, solution,&
                         pardiso_solve,workspace%dparm)
        endif 
        ! Numeric factorization and solution
        phase = 23
        call pardiso(workspace%pt, 1, 1, options%mtype, phase, matrix%order,&
                     matrix%val(1:matrix%nonzeros),&
                     matrix%row(1:matrix%order+1),&
                     matrix%col(1:matrix%nonzeros),&
                     1, 1, options%iparm, &
                     mesglvl, rhs, solution, pardiso_solve,workspace%dparm)
#endif
        if (options%iterative) then
            niters = options%iparm(20) ! CGS iterations
	else
            niters = options%iparm(7) ! Direct iterations
	endif
        err=1.0
    end function pardiso_solve

    subroutine pardiso_solver_postprocess(matrix, options, workspace)
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options) :: options
        type(pardiso_solver_workspace) :: workspace

    end subroutine

    subroutine pardiso_destroy_workspace(matrix, options, workspace)
        !*FD Deallocates all working memory for the sparse linear solver.
        !*FD This need *not* be safe to call of an unallocated workspace
        !*FD No sparse solver should call this automatically.
        type(sparse_matrix_type) :: matrix
        type(pardiso_solver_options) :: options
        type(pardiso_solver_workspace) :: workspace

        ! I believe that all memory clearing has been accomplised in the solve.

    end subroutine pardiso_destroy_workspace

    subroutine pardiso_interpret_error(error_code, error_string)
        !*FD takes an error code output from sparse_solve and interprets it.
        !*FD error_string must be an optional argument.
        !*FD If it is not provided, the error is printed to standard out
        !*FD instead of being put in the string
        integer :: error_code
        character(*), optional, intent(out) :: error_string
        character(256) :: tmp_error_string
        
        select case (error_code)
            case (0)
                tmp_error_string = "No error"
            case (-1)
                tmp_error_string = "Input inconsistent"
            case (-2)
                tmp_error_string = "Not enough memory"
            case (-3)
                tmp_error_string = "Reordering problem"
            case (-4)
                tmp_error_string = "Zero pivot, numerical fact. or iterative refinement problem"
            case (-5)
                tmp_error_string = "Unclassified (internal) errror"
            case (-6)
                tmp_error_string = "Preordering failed (matrix types 11, 13 only)"
            case (-7)
                tmp_error_string = "Diagonal matrix problem"
            case (-8)
                tmp_error_string = "32-bit integer overflow problem"
            case (-10)
                tmp_error_string = "No license file pardiso.lic found"
            case (-11)
                tmp_error_string = "License expired."
            case (-12)
                tmp_error_string = "Wrong username or hostname."
            case (-100)
                tmp_error_string = "Maximum Krylov iterations reached."
            case (-101)
                tmp_error_string = "Insufficient convergence in Krylov subspace iteration in 25 iterations."
            case (-102)
                tmp_error_string = "Error in Krylov-subspace iteration."
            case (-103)
                tmp_error_string = "Break down in Krylov-subspace iteration."
        end select

        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif
    end subroutine pardiso_interpret_error
end module glimmer_sparse_pardiso
