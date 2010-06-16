! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_sparse.f90 - part of the Glimmer-CISM ice model  + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glimmer_sparse
  use glimmer_global
  type sparse_matrix_type
     !*FD sparse matrix type
     integer n                                              !*FD number of elements
     integer, dimension(:), pointer :: col => NULL()        !*FD column index
     integer, dimension(:), pointer :: row => NULL()        !*FD row index
     real(kind=dp), dimension(:), pointer :: val => NULL()  !*FD values
  end type sparse_matrix_type

  ! size of sparse matrix 
  integer, parameter, private :: chunksize=1000

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

  subroutine new_sparse_matrix(n,mat)
    !*FD create a new sparse matrix
    implicit none
    integer, intent(in) :: n          !*FD initial size of matrix
    type(sparse_matrix_type) :: mat   !*FD matrix
    
    if (.not.associated(mat%col)) then
       allocate(mat%row(n))
       allocate(mat%col(n))
       allocate(mat%val(n))
    else
       if (size(mat%row).lt.n) then
          call del_sparse_matrix(mat)
          allocate(mat%row(n))
          allocate(mat%col(n))
          allocate(mat%val(n))
       end if
    end if
    mat%n = 0
  end subroutine new_sparse_matrix

  subroutine copy_sparse_matrix(inmat,outmat)
    !*FD copy a sparse matrix
    implicit none
    type(sparse_matrix_type) :: inmat  !*FD matrix to be copied
    type(sparse_matrix_type) :: outmat !*FD result matrix

    call new_sparse_matrix(inmat%n,outmat)
    outmat%row(:) = inmat%row(:)
    outmat%col(:) = inmat%col(:)
    outmat%val(:) = inmat%val(:)
    outmat%n = inmat%n
  end subroutine copy_sparse_matrix

  subroutine grow_sparse_matrix(matrix)
    !*FD grow sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix

    integer, dimension(:), pointer :: newrow,newcol
    real(kind=dp), dimension(:), pointer :: newval
    integer oldsize

    oldsize = size(matrix%val)
    
    allocate(newrow(chunksize+oldsize))
    allocate(newcol(chunksize+oldsize))
    allocate(newval(chunksize+oldsize))

    newcol(1:oldsize) = matrix%col(:)
    newrow(1:oldsize) = matrix%row(:)
    newval(1:oldsize) = matrix%val(:)

    deallocate(matrix%col)
    deallocate(matrix%row)
    deallocate(matrix%val)

    matrix%col => newcol
    matrix%row => newrow
    matrix%val => newval

  end subroutine grow_sparse_matrix

  subroutine del_sparse_matrix(matrix)
    !*FD delete sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix

    deallocate(matrix%col)
    deallocate(matrix%row)
    deallocate(matrix%val)
  end subroutine del_sparse_matrix

  subroutine print_sparse(matrix, unit)
    !*FD print sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix
    integer, intent(in) :: unit        !*FD unit to be printed to

    integer i
    do i = 1, matrix%n
       write(unit,*) matrix%col(i), matrix%row(i), matrix%val(i)
    end do
  end subroutine print_sparse

  subroutine sparse_matrix_vec_prod(matrix, vec, res)
    !*FD sparse matrix vector product
    implicit none
    type(sparse_matrix_type) :: matrix                !*FD matrix
    real(kind=dp), intent(in), dimension(:) :: vec    !*FD input vector
    real(kind=dp), intent(out), dimension(:) :: res   !*FD result vector

    integer i

    res = 0.
    do i=1,matrix%n
       res(matrix%col(i)) = res(matrix%col(i)) + vec(matrix%row(i))*matrix%val(i)
    end do
  end subroutine sparse_matrix_vec_prod

  subroutine sparse_insert_val(matrix, i, j, val)
    !*FD insert value into sparse matrix
    implicit none
    type(sparse_matrix_type) :: matrix !*FD matrix
    integer, intent(in) :: i,j         !*FD column and row
    real(kind=dp), intent(in) :: val   !*FD value

    matrix%n =  matrix%n + 1
    matrix%col(matrix%n) = i
    matrix%row(matrix%n) = j
    matrix%val(matrix%n) = val

    if (matrix%n .eq. size(matrix%val)) then
       call grow_sparse_matrix(matrix)
    end if
  end subroutine sparse_insert_val
end module glimmer_sparse
