!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   isostasy_elastic.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module isostasy_elastic

  !> handle elastic lithosphere

  !> NOTE: The elastic lithosphere calculation is done on a single task.
  !>       For parallel runs, data are gathered onto the main task for the computation,
  !>        then scattered back to local processors.
  !>       This procedure is manageable on a 4 km mesh, but may be too expensive and/or
  !>        memory-intensive at higher resolutions. A more scalable approach is a priority
  !>        for future development.

  use glimmer_global, only : dp
  use glide_types, only: isos_elastic

  implicit none

  real(dp), private, parameter :: r_lr = 6.d0   ! influence of disk load at (0,0) is felt within a radius of r_lr*rbel_r

  private :: init_rbel, rbel_ow, rbel_iw

  logical, parameter :: verbose_isostasy = .false.  ! if true, print diagnostic messages

!-------------------------------------------------------------------------

contains

!-------------------------------------------------------------------------
  
  subroutine init_elastic(rbel, deltax)

    !> initialise elastic lithosphere calculations

    use glimmer_physcon, only : pi
    implicit none

    type(isos_elastic) :: rbel          !> structure holding elastic lithosphere data
    real(dp), intent(in) :: deltax      !> grid spacing

    ! local variables
    real(dp) :: a     ! radius of disk
    real(dp) :: r     ! distance from centre
    integer :: i,j

    ! calculate a so that a circle of radius a is equivalent to a square with size deltax
    a = deltax/sqrt(pi)

    ! initialise w
    call init_rbel(rbel, a)

    ! calculate size of operator
    rbel%wsize = int(r_lr*rbel%lr/deltax)

    ! allocate memory for operator
    allocate(rbel%w(0:rbel%wsize,0:rbel%wsize))

    ! calculating points within disk
    rbel%w(0,0) = rbel_iw(rbel,0.d0)
    r = deltax/rbel%lr
    rbel%w(0,1) = rbel_iw(rbel,r)
    rbel%w(1,0) = rbel%w(0,1)

    ! calculating points outside disk
    do j=0,rbel%wsize
       do i=2,rbel%wsize
          r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do

    do j=2,rbel%wsize
       do i=0,1
          r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
          rbel%w(i,j) = rbel_ow(rbel,r)
       end do
    end do

    i=1
    j=1
    r = deltax * sqrt(real(i)**2 + real(j)**2)/rbel%lr
    rbel%w(i,j) = rbel_ow(rbel,r)

#ifdef DEB_REBOUND
    open(1,file='w.dat',status='UNKNOWN')
    do j=0,rbel%wsize
       do i=0,rbel%wsize
          write(1,*) i,j,rbel%w(i,j)
       end do
    end do
    close(1)
#endif

    !rbel%w=rbel%w/len0

  end subroutine init_elastic

!-------------------------------------------------------------------------

  subroutine calc_elastic(&
       rbel,  load_factors,  load,  &
       idiag, jdiag,                &
       idiag_local, jdiag_local, rdiag_local)

    !> Calculate surface loading effect using elastic lithosphere approximation.
    !> Functionally equivalent to subroutine calc_elastic from Glimmer's original isostasy model.
    !> The main difference is that this subroutine uses a global gather and scatter to compute
    !>  the load for simulations on more than one task.

    use parallel, only : global_ewn, global_nsn, this_rank, main_task, nhalo, &
                         distributed_gather_var, distributed_scatter_var, parallel_halo
    use parallel, only : parallel_reduce_sum  ! diagnostic only

    implicit none

    type(isos_elastic) :: rbel                             !> structure holding elastic litho data
    real(dp), dimension(:,:), intent(in)  :: load_factors  !> load mass divided by mantle density
    real(dp), dimension(:,:), intent(out) :: load          !> loading effect due to load_factors

    ! The following are needed only for diagnostic prints
    integer, intent(in) :: &
         idiag, jdiag                              !> global coordinates of diagnostic point
    integer, intent(in) :: &
         idiag_local, jdiag_local, rdiag_local     !> local coordinates of diagnostic point

    ! local variables

    integer :: ewn, nsn    !> grid dimensions on the local task; includes halo cells
    integer :: i, j, n, m

    real(dp), dimension(:,:), allocatable :: &
         load_global,             & !> global version of the output 'load' array
         load_factors_global        !> global version of the input 'load_factors' array

    real(dp) :: local_sum_load, global_sum_load   !> diagnostic sums

    ! initialize

    ewn = size(load,1)
    nsn = size(load,2)
    load(:,:) = 0.0d0

    if (verbose_isostasy .and. main_task) then
       print*, 'ISOSTASY: calc_elastic'
       print*, 'local ewn/nsn =', ewn, nsn
       print*, 'global_ewn/nsn =', global_ewn, global_nsn
    endif

    ! Gather the local arrays onto the main task
    ! Note: global arrays are allocated in the subroutine
    call distributed_gather_var(load_factors, load_factors_global)
    call distributed_gather_var(load, load_global)

    if (main_task) then
       do j = 1, global_nsn

          if (verbose_isostasy .and. main_task) then
             if (mod(j,100) == 0) print*, 'j =', j   ! to see how fast the calculation is going
          endif
          
          do i = 1, global_ewn

             ! Compute load terms by summing over cells in the radius of influence
             do n = max(1,j-rbel%wsize), min(global_nsn,j+rbel%wsize)
                do m = max(1,i-rbel%wsize), min(global_ewn,i+rbel%wsize)
                   load_global(i,j) = load_global(i,j) + load_factors_global(m,n) * rbel%w(abs(m-i),abs(n-j))
                end do
             end do

          end do  ! i
       end do  ! j
    endif  ! main_task

    ! Scatter the load values back to local arrays
    ! Note: The global array is deallocated in the subroutine
    call distributed_scatter_var(load, load_global)

    ! distributed_scatter_var does not update the halo, so do an update here
    call parallel_halo(load)

    ! Deallocate the other global array (which is intent(in) and does not need to be scattered)
    deallocate(load_factors_global)

    if (verbose_isostasy .and. main_task) then

       ! print value at diagnostic point
       if (this_rank==rdiag_local) then
          i = idiag_local
          j = jdiag_local
          print*, 'ISOSTASY: r, i, j, load:', rdiag_local, i, j, load(i,j)
       endif

    endif  ! verbose_isostasy

  end subroutine calc_elastic

!-------------------------------------------------------------------------

  subroutine finalise_elastic(rbel)
    !> clean-up data structure
    implicit none
    type(isos_elastic) :: rbel     !> structure holding elastic litho data    

    deallocate(rbel%w)
  end subroutine finalise_elastic

!-------------------------------------------------------------------------

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! private subroutines
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine init_rbel(rbel, a)

    !> initialise elastic lithosphere calculations
    use glimmer_paramets, only: len0
    use glimmer_physcon, only: rhom,grav
    use isostasy_kelvin
    implicit none
    type(isos_elastic) :: rbel        !> structure holding elastic litho data
    real(dp), intent(in) :: a             !> radius of disk

    real(dp) :: dummy_a

    call set_kelvin(1.d-10,40)

    rbel%lr = ((rbel%d/(rhom*grav))**0.25d0)/len0
    rbel%a  = a

    dummy_a = rbel%a/rbel%lr
    
    rbel%c1  =  dummy_a * dker(dummy_a)
    rbel%c2  = -dummy_a * dkei(dummy_a)
    rbel%cd3 =  dummy_a * dber(dummy_a)
    rbel%cd4 = -dummy_a * dbei(dummy_a)
    
  end subroutine init_rbel

!-------------------------------------------------------------------------

  function rbel_ow(rbel,r)
    use isostasy_kelvin
    !> calculating deflection outside disk
    implicit none
    real(dp) :: rbel_ow
    real(dp), intent(in) :: r          !> radius, r should be scaled with lr
    type(isos_elastic) :: rbel     !> structure holding elastic litho data
    
    rbel_ow = rbel%cd3*ker(r) + rbel%cd4*kei(r)
  end function rbel_ow

!-------------------------------------------------------------------------

  function rbel_iw(rbel,r)
    use isostasy_kelvin
    !> calculating deflection inside disk
    implicit none
    real(dp) :: rbel_iw
    real(dp), intent(in) :: r          !> radius, r should be scaled with lr
    type(isos_elastic) :: rbel         !> structure holding elastic litho data
    
    rbel_iw = 1.d0 + rbel%c1*ber(r) + rbel%c2*bei(r)
  end function rbel_iw

!-------------------------------------------------------------------------

end module isostasy_elastic

!-------------------------------------------------------------------------
