! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_sparse.f90 - part of the GLIMMER ice model       + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_sparse
  use glimmer_global, only: sp, dp
  use glimmer_sparse_type
  use glimmer_sparse_slap
  use glimmer_sparse_umfpack
  use glimmer_sparse_pardiso
  implicit none

  type sparse_solver_options
        type(sparse_solver_options_base) :: base
        type(slap_solver_options) :: slap
        type(umf_solver_options)  :: umf
        type(pardiso_solver_options)  :: pardiso
  end type

  type sparse_solver_workspace
        type(slap_solver_workspace), pointer :: slap => null()
        type(umf_solver_workspace),  pointer :: umf  => null()
        type(pardiso_solver_workspace),  pointer :: pardiso  => null()
  end type

  integer, parameter :: SPARSE_SOLVER_BICG = 0
  integer, parameter :: SPARSE_SOLVER_GMRES = 1
  integer, parameter :: SPARSE_SOLVER_UMF = 2
  integer, parameter :: SPARSE_SOLVER_PARDISO = 3

!EIB!  !MAKE_RESTART
!EIB!#ifdef RESTARTS
!EIB!#define RST_GLIMMER_SPARSE
!EIB!#include "glimmer_rst_head.inc"
!EIB!#undef RST_GLIMMER_SPARSE
!EIB!#endif

contains

!EIB!#ifdef RESTARTS
!EIB!#define RST_GLIMMER_SPARSE
!EIB!#include "glimmer_rst_body.inc"
!EIB!#undef RST_GLIMMER_SPARSE
!EIB!#endif

    subroutine sparse_solver_default_options(method, opt)
        integer, intent(in) :: method
        type(sparse_solver_options) :: opt

        opt%base%method = method
        opt%base%tolerance  = 1e-11
        opt%base%maxiters = 200

        !Solver specific options
        if (method == SPARSE_SOLVER_BICG) then
            call slap_default_options(opt%slap, opt%base) 

        else if (method == SPARSE_SOLVER_GMRES) then
            call slap_default_options(opt%slap, opt%base)
            opt%slap%use_gmres = .true.

        else if (method == SPARSE_SOLVER_UMF) then
            call umf_default_options(opt%umf)

        else if (method == SPARSE_SOLVER_PARDISO) then
            call pardiso_default_options(opt%pardiso)
        else 
            !call glide_finalise_all(.true.)
            call write_log("Invalid sparse matrix option used.", GM_FATAL)
        end if
    end subroutine

    subroutine sparse_allocate_workspace(matrix, options, workspace, max_nonzeros_arg)
        !*FD Allocate solver workspace.  This needs to be done once
        !*FD (when the maximum number of nonzero entries is first known)
        !*FD This function need not be safe to call on already allocated memory
        !*FD
        !*FD Note that the max_nonzeros argument must be optional, and if
        !*FD it is not supplied the current number of nonzeroes must be used.
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace
        integer, optional :: max_nonzeros_arg
        integer :: max_nonzeros
        
        if (present(max_nonzeros_arg)) then
            max_nonzeros = max_nonzeros_arg
        else
            max_nonzeros = matrix%nonzeros
        end if

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES) then
            allocate(workspace%slap)
            call slap_allocate_workspace(matrix, options%slap, workspace%slap, max_nonzeros)

        else if (options%base%method == SPARSE_SOLVER_UMF) then
            allocate(workspace%umf)
            call umf_allocate_workspace(matrix, options%umf, workspace%umf, max_nonzeros)
        else if (options%base%method == SPARSE_SOLVER_PARDISO) then
            allocate(workspace%pardiso)
            call pardiso_allocate_workspace(matrix, options%pardiso, workspace%pardiso, max_nonzeros)
        end if
    end subroutine sparse_allocate_workspace

    subroutine sparse_solver_preprocess(matrix, options, workspace)
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
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES) then
            call slap_solver_preprocess(matrix, options%slap, workspace%slap)

        else if (options%base%method == SPARSE_SOLVER_UMF) then
            call umf_solver_preprocess(matrix, options%umf, workspace%umf)

        else if (options%base%method == SPARSE_SOLVER_PARDISO) then
            call pardiso_solver_preprocess(matrix, options%pardiso, workspace%pardiso)
        end if

    end subroutine sparse_solver_preprocess

    function sparse_solve(matrix, rhs, solution, options, workspace,err,niters, verbose)
        !*FD Solves the linear system, and reports status information.
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

        type(sparse_solver_options), intent(in) :: options
        !*FD Options such as convergence criteria
        
        type(sparse_solver_workspace), intent(inout) :: workspace
        !*FD Internal solver workspace
        
        real(kind=dp), intent(out) :: err
        !*FD Final solution error
        
        integer, intent(out) :: niters
        !*FD Number of iterations required to reach the solution

        logical, optional, intent(in) :: verbose
        !*FD If present and true, this argument may cause diagnostic information
        !*FD to be printed by the solver (not every solver may implement this).
        
        integer :: sparse_solve

        logical :: verbose_var

        verbose_var = .false.
        if (present(verbose)) then
            verbose_var = verbose
        end if

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES) then
            sparse_solve = slap_solve(matrix, rhs, solution, options%slap, workspace%slap, err, niters, verbose_var)

        else if (options%base%method == SPARSE_SOLVER_UMF) then
            sparse_solve = umf_solve(matrix, rhs, solution, options%umf, workspace%umf, err, niters, verbose_var)

        else if (options%base%method == SPARSE_SOLVER_PARDISO) then
            sparse_solve = pardiso_solve(matrix, rhs, solution, options%pardiso,&
                                         workspace%pardiso, err,niters, verbose_var)
        end if
    end function sparse_solve

    subroutine sparse_solver_postprocess(matrix, options, workspace)
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace

        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES) then
            call slap_solver_postprocess(matrix, options%slap, workspace%slap)

        else if (options%base%method == SPARSE_SOLVER_UMF) then
            call umf_solver_postprocess(matrix, options%umf, workspace%umf)

        else if (options%base%method == SPARSE_SOLVER_PARDISO) then
            call pardiso_solver_postprocess(matrix, options%pardiso, workspace%pardiso)

        end if
    end subroutine

    subroutine sparse_destroy_workspace(matrix, options, workspace)
        !*FD Deallocates all working memory for the slap linear solver.
        !*FD This need *not* be safe to call of an unallocated workspace
        !*FD No slap solver should call this automatically.
        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: options
        type(sparse_solver_workspace) :: workspace
        
        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES) then
            call slap_destroy_workspace(matrix, options%slap, workspace%slap)
            deallocate(workspace%slap)

        else if (options%base%method == SPARSE_SOLVER_UMF) then
            call umf_destroy_workspace(matrix, options%umf, workspace%umf)
            deallocate(workspace%umf)

        else if (options%base%method == SPARSE_SOLVER_PARDISO) then
            call pardiso_destroy_workspace(matrix, options%pardiso, workspace%pardiso)
            deallocate(workspace%pardiso)
        end if

    end subroutine sparse_destroy_workspace

    subroutine sparse_interpret_error(options, error_code, error_string)
        !*FD takes an error code output from slap_solve and interprets it.
        !*FD error_string must be an optional argument.
        !*FD If it is not provided, the error is printed to standard out
        !*FD instead of being put in the string
        type(sparse_solver_options) :: options
        integer :: error_code
        character(*), optional, intent(out) :: error_string
        character(256) :: tmp_error_string
        
        if (options%base%method == SPARSE_SOLVER_BICG .or. &
            options%base%method == SPARSE_SOLVER_GMRES) then
            call slap_interpret_error(error_code, tmp_error_string)

        else if (options%base%method == SPARSE_SOLVER_UMF) then
            call umf_interpret_error(error_code, tmp_error_string)

        else if (options%base%method == SPARSE_SOLVER_PARDISO) then
            call pardiso_interpret_error(error_code, tmp_error_string)
        end if


        if (present(error_string)) then
            error_string = tmp_error_string
        else
            write(*,*) tmp_error_string
        endif
    end subroutine sparse_interpret_error

    subroutine sparse_easy_solve(matrix, rhs, answer, err, iter, method_arg, calling_file, calling_line)
        !This subroutine wraps the basic (though probably the most inefficient)
        !workflow to solve a sparse matrix using the sparse matrix solver
        !framework.  It handles errors gracefully, and reports back the
        !iterations required and the error estimate in the case of an iterative
        !solver.  At the very least it is an encapsulated example of how to
        !use the sparse solver routines, and is easy enough to drop in your
        !code if you don't care about allocating and deallocating workspace
        !every single timestep.

        type(sparse_matrix_type) :: matrix
        real(dp), dimension(:) :: rhs
        real(dp), dimension(:) :: answer
        
        real(dp), intent(out) :: err
        integer, intent(out) :: iter
        
        integer, optional, intent(in) :: method_arg

        character(*), optional :: calling_file
        integer, optional :: calling_line

        type(sparse_solver_options) :: opt
        type(sparse_solver_workspace) :: wk

        integer :: ierr
        integer :: method

        if (present(method_arg)) then
            method = method_arg
        else
            method = SPARSE_SOLVER_BICG
        endif

        call sparse_solver_default_options(method, opt)
        call sparse_allocate_workspace(matrix, opt, wk)
        call sparse_solver_preprocess(matrix, opt, wk)

        ierr = sparse_solve(matrix, rhs, answer, opt, wk, err, iter, .false.)
       
        call sparse_solver_postprocess(matrix, opt, wk)

        if (ierr /= 0) then
            if (present(calling_file) .and. present(calling_line)) then
                call handle_sparse_error(matrix, opt, ierr, calling_file, calling_line)
            else
                call handle_sparse_error(matrix, opt, ierr, __FILE__, __LINE__)
            end if
        end if

        call sparse_destroy_workspace(matrix, opt, wk)

    end subroutine sparse_easy_solve

    subroutine handle_sparse_error(matrix, solver_options, error, error_file, error_line, time)
        !Checks a sparse error flag and, if an error occurred, log it to
        !the GLIMMER log file.  This does not stop Glimmer, it just writes
        !to the log
        !use glide_stop
        use glimmer_log
        use glimmer_filenames
        
        integer :: error
        integer, optional :: error_line
        character(*), optional :: error_file
        real(dp), optional :: time

        type(sparse_matrix_type) :: matrix
        type(sparse_solver_options) :: solver_options
        integer :: isym
        integer :: lunit
        integer :: i

        character(512) :: message
        character(128) :: errfname
        character(256) :: errdesc

        !If no error happened, this routine should be a nop
        if (error == 0 .OR. error == 2 .OR. error == 6) return

        !Aquire a file unit, and open the file
        lunit = get_free_unit()
        errfname = trim(process_path('sparse_dump.txt'))
        open(lunit,file=errfname)

        if (matrix%symmetric) then
            isym = 1
        else
            isym = 0
        end if

        !Output sparse matrix data to the file
        call dcpplt(matrix%order, matrix%nonzeros, matrix%row, matrix%col, matrix%val,&
                    isym, lunit)

        write(lunit,*) '***Sparse matrix structure ends.  Value listing begins'
        do i=1,matrix%nonzeros
            write(lunit,*) matrix%val(i)
        end do

        !Close unit and finish off
        close(lunit)
        
        !Grab the error message from the sparse solver
        call sparse_interpret_error(solver_options, error, errdesc)

        !construct the error message and write it to the log file
        if (present(time)) then
            write(message, *)'Sparse matrix error at time: ', time, &
                             'Error description: ', errdesc, &
                             'Data dumped to ', trim(errfname)
        else
            write(message, *)'Sparse matrix error. Error description: ', errdesc, &
                             'Data dumped to ', trim(errfname)
        end if

        write(*,*)message

        !call glide_finalise_all(.true.)

        if (present(error_file) .and. present(error_line)) then
            call write_log(trim(errdesc), GM_FATAL, error_file, error_line)
        else
            call write_log(trim(errdesc), GM_FATAL, __FILE__, __LINE__)
        end if
    end subroutine handle_sparse_error

end module glimmer_sparse
