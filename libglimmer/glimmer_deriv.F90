! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_deriv.f90 - part of the Glimmer-CISM ice model   + 
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

!> This module contains functions for computing derivatives numerically, both
!! for a single value and for an entire field.
!!
!! \author Tim Bocek
!! \date 2008-10-15
!!
!! Note that this module is written with the first index in a matrix corresponding
!! to the x (east-west) coordinate.  If this is not the case (i.e. if the first
!! index corresponds to the y (north-south) coordinate), then transposition
!! will be necessary.  Simply ask for the y-derivative when you mean to ask for
!! the x-derivative, and vice versa.
module glimmer_deriv

  use glimmer_global, only: sp, dp
  !------------------------------------------------------------------
  !First Derivative Estimates, Second Order, 2D
  !------------------------------------------------------------------
contains

  !> Computes derivative with respect to x at a given point.
  !! Applies periodic boundary conditions if needed.
  function dfdx_2d(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: dfdx_2d

    dfdx_2d = (-.5/delta)*f(i-1, j) + (.5/delta)*f(i+1, j)
    !write(*,*), i, j, f(i,j), ip1, im1, delta, dfdx_2d
  end function dfdx_2d

  !> Computes derivative with respect to y at a given point
  function dfdy_2d(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: dfdy_2d

    integer :: jp1, jm1
    jp1 = j + 1
    jm1 = j - 1
    if (jp1 == size(f, 2)+1) jp1 = 2
    if (jm1 == 0) jm1 = size(f, 2)-1

    dfdy_2d = (-.5/delta)*f(i, j-1) + (.5/delta)*f(i, j+1)
  end function dfdy_2d

  !> Computes derivative with respect to x at the equivalent
  !! point on a staggered grid.
  function dfdx_2d_stag(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: dfdx_2d_stag
    dfdx_2d_stag = (f(i+1, j) + f(i+1, j+1) - f(i, j) - f(i, j+1))/(2*delta) 
  end function dfdx_2d_stag

  !> Computes derivative with respect to y at the equivalent
  !! point on a staggered grid.
  function dfdy_2d_stag(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: dfdy_2d_stag
    dfdy_2d_stag = (f(i, j+1) + f(i+1, j+1) - f(i,j) - f(i+1, j))/(2*delta)
  end function dfdy_2d_stag

  !> Computes derivative with respect to x at the given point
  !! using an upwind method (suitable for maximum boundaries)
  function dfdx_2d_upwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: dfdx_2d_upwind
    dfdx_2d_upwind = (.5 * f(i-2,j) - 2 * f(i-1, j) + 1.5 * f(i, j))/delta
  end function dfdx_2d_upwind

  !> Computes derivative with respect to y at the given point
  !! using an upwind method (suitable for maximum boundaries)
  function dfdy_2d_upwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: dfdy_2d_upwind
    dfdy_2d_upwind = (.5 * f(i,j-2) - 2 * f(i, j-1) + 1.5 * f(i, j))/delta
  end function dfdy_2d_upwind

  !> Computes derivative with respect to x at the given point
  !! using a downwind method (suitable for minimum boundaries)
  function dfdx_2d_downwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: dfdx_2d_downwind
    dfdx_2d_downwind = (-1.5 * f(i, j) + 2 * f(i+1, j) - .5 * f(i+2, j))/delta
  end function dfdx_2d_downwind

  !> Computes derivative with respect to y at the given point
  !! using a downwind method (suitable for minimum boundaries)
  function dfdy_2d_downwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: dfdy_2d_downwind
    dfdy_2d_downwind = (-1.5 * f(i, j) + 2 * f(i, j+1) - .5 * f(i, j+2))/delta
  end function dfdy_2d_downwind

  !> Computes derivative fields of the given function.
  subroutine df_field_2d(f, deltax, deltay, out_dfdx, out_dfdy, periodic_x, periodic_y, direction_x, direction_y)
    implicit none
    real(dp), dimension(:, :), intent(in) :: f                      !< field to be derived
    real(dp), intent(in) :: deltax                                  !< grid spacing in x direction
    real(dp), intent(in) :: deltay                                  !< grid spacing in y direction
    real(dp), dimension(:, :), intent(out) :: out_dfdx              !< derivatives in x direction
    real(dp), dimension(:, :), intent(out) :: out_dfdy              !< derivatives in y direction
    real(dp), dimension(:, :), intent(in), optional  :: direction_x !< x-direction used for upwind derivatives
    real(dp), dimension(:, :), intent(in), optional  :: direction_y !< y-direction used for upwind derivatives

    logical :: upwind !Whether or not directions for upwinding were provided

    integer :: grad_x, grad_y !Whether to upwind or downwind at the current point


    integer :: nx, ny, x, y
    logical :: periodic_x, periodic_y

    !Get the size of the field we're working with
    nx = size(f, 1)
    ny = size(f, 2)

    upwind = present(direction_x) .and. present(direction_y)

    !! \todo refactor this, put the if (upwind) test outside the loop

    !For now, we'll use the function calls defined above.
    !Later on we might want to refactor?
    do x = 1, nx
       do y = 1, ny
          grad_x = 0
          grad_y = 0
          if (upwind) then
             if (direction_x(x,y) < 0 .and. x > 2) then !Upstream case
                grad_x = -1
             else if(direction_x(x,y) > 0 .and. x < nx - 1) then !Downstream case
                grad_x = 1
             end if

             if (direction_y(x,y) < 0 .and. y > 2) then !Upstream case
                grad_y = -1
             else if(direction_y(x,y) > 0 .and. y < ny - 1) then !Downstream case
                grad_y = 1
             end if
          end if

          !For each of the variables in x, y, check whether or not
          !we need to use an upwinding or downwinding differentiation
          !scheme.
          if (x == 1 .or. grad_x > 0) then
             out_dfdx(x, y) = dfdx_2d_downwind(f, x, y, deltax)
          else if (x == nx .or. grad_x < 0) then
             out_dfdx(x, y) = dfdx_2d_upwind(f, x, y, deltax)
          else
             out_dfdx(x, y) = dfdx_2d(f, x, y, deltax)
          end if

          if (y == 1 .or. grad_y > 0) then
             out_dfdy(x, y) = dfdy_2d_downwind(f, x, y, deltay)
          elseif (y == ny .or. grad_y < 0) then
             out_dfdy(x, y) = dfdy_2d_upwind(f, x, y, deltay)
          else
             out_dfdy(x, y) = dfdy_2d(f, x, y, deltay)
          end if

       end do
    end do

  end subroutine df_field_2d

  !> Computes derivative fields of the given function.  Places the result
  !! on a staggered grid.  If periodic in one dimension is set, that 
  !! dimension for derivatives must be the same size as the value's dimension.
  !! Otherwise, it should be one less
  !! \todo enable periodic boundary conditions
  subroutine df_field_2d_staggered(f, deltax, deltay, out_dfdx, out_dfdy, periodic_x, periodic_y)
    implicit none
    real(dp), dimension(:, :), intent(in) :: f                      !< field to be derived
    real(dp), intent(in) :: deltax                                  !< grid spacing in x direction
    real(dp), intent(in) :: deltay                                  !< grid spacing in y direction
    real(dp), dimension(:, :), intent(out) :: out_dfdx              !< derivatives in x direction
    real(dp), dimension(:, :), intent(out) :: out_dfdy              !< derivatives in y direction
    logical :: periodic_x                                           !< flag for periodic BC in x
    logical :: periodic_y                                           !< flag for periodic BC in y

    integer :: nx, ny, x, y

    !Get the size of the field we're working with
    nx = size(f, 1)
    ny = size(f, 2)

    !For now, we'll use the function calls defined above.
    !Later on we might want to refactor?    
    do x = 1, nx - 1 !We go to nx - 1 because we're using a staggered grid
       do y = 1, ny - 1
          out_dfdx(x,y) = dfdx_2d_stag(f, x, y, deltax)
          out_dfdy(x,y) = dfdy_2d_stag(f, x, y, deltay)
       end do
    end do

    !               !Deal with periodic boundary conditions.  We will do so by
    !               !providing another set of values at the end of each dimension
    !               !that contains the derivative of the value off the edge of the
    !               !grid.  Because this set of values is at the end, when
    !               !x = nx, x+1 = 1.  This identity has been hard-coded below.
    !               if (periodic_x) then
    !                       do y = 1, ny - 1
    !                               out_dfdx(nx,y) = -(f(1, y) + f(1, y+1) - f(nx, y) - f(nx, y+1))/(2*deltax)
    !                               out_dfdy(nx,y) = dfdy_2d_stag(f, nx, y, deltay)
    !                       end do
    !               end if
    !               
    !               if (periodic_y) then
    !                       do x = 1, nx - 1
    !                           out_dfdx(x,ny) = dfdx_2d_stag(f, x, ny, deltax)
    !                               out_dfdy(x,ny) = -(f(x, 1) + f(x+1, 1) - f(x,ny) - f(x+1, ny))/(2*deltay)
    !                       end do
    !               end if
    !               
    !               !Do the corner that hasn't been done if both dimensions are periodic
    !               if (periodic_x .and. periodic_y) then
    !                       out_dfdx(nx,ny) = (f(1, ny) + f(1, 1) - f(nx, ny) - f(nx, 1))/(2*deltax)
    !                       out_dfdy(nx,ny) = (f(nx, 1) + f(1, 1) - f(nx,ny)  - f(1, ny))/(2*deltay)
    !               end if
    !               
  end subroutine df_field_2d_staggered
  !------------------------------------------------------------------
  !First Derivative Estimates, Second Order, 3D
  !------------------------------------------------------------------

  !> Computes derivative with respect to x at a given point
  function dfdx_3d(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in x-direction
    real(dp) :: dfdx_3d
    dfdx_3d = (-.5/delta)*f(k, i-1, j)  + (.5/delta)*f(k, i+1, j)
  end function dfdx_3d

  !> Computes derivative with respect to y at a given point
  function dfdy_3d(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in y-direction
    real(dp) :: dfdy_3d
    dfdy_3d = (-.5/delta)*f(k, i, j-1) + (.5/delta)*f(k, i, j+1)
  end function dfdy_3d



  !> Computes derivative with respect to z at a given point
  !! where the Z axis uses an irregular grid defined by deltas.
  !! This derivative is given by the formula:
  function dfdz_3d_irregular(f, i, j, k, deltas)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f  !< field to be derived
    integer, intent(in) :: i                     !< first index at which derivative is taken
    integer, intent(in) :: j                     !< second index at which derivative is taken
    integer, intent(in) :: k                     !< third index at which derivative is taken
    real(dp), dimension(:), intent(in) :: deltas !< vertical sigma levels
    real(dp) :: dfdz_3d_irregular

    dfdz_3d_irregular = f(k-1,i,j)*(deltas(k) - deltas(k+1))/((deltas(k) - deltas(k-1))*(deltas(k+1)-deltas(k-1))) + &
         f(k,  i,j)*(deltas(k+1)-2*deltas(k)+deltas(k-1))/((deltas(k)-deltas(k-1))*(deltas(k+1)-deltas(k))) + &
         f(k+1,i,j)*(deltas(k)-deltas(k-1))/((deltas(k+1)-deltas(k))*(deltas(K+1)-deltas(k-1)))
  end function dfdz_3d_irregular

  !> Computes derivative with respect to z at a given point using an upwinding
  !! scheme.  The Z axis uses an irregular grid defined by deltas.
  function dfdz_3d_upwind_irregular(f, i, j, k, deltas)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f  !< field to be derived
    integer, intent(in) :: i                     !< first index at which derivative is taken
    integer, intent(in) :: j                     !< second index at which derivative is taken
    integer, intent(in) :: k                     !< third index at which derivative is taken
    real(dp), dimension(:), intent(in) :: deltas !< vertical sigma levels
    real(dp) :: dfdz_3d_upwind_irregular
    real(dp) :: zkMinusZkm1, zkMinusZkm2, zkm1MinusZkm2
    zkMinusZkm1 = deltas(k) - deltas(k-1)
    zkMinusZkm2 = deltas(k) - deltas(k-2)
    zkm1MinusZkm2 = deltas(k-1) - deltas(k-2)

    dfdz_3d_upwind_irregular = f(k-2, i, j) * zkMinusZkm1 / (zkm1MinusZkm2 * zkMinusZkm2) - &
         f(k-1, i, j) * zkMinusZkm2 / (zkMinusZkm1 * zkm1MinusZkm2) + &
         f(k,   i, j) * (2*deltas(k) - deltas(k-1) - deltas(k-2)) / (zkMinusZkm1 * zkMinusZkm2)
  end function dfdz_3d_upwind_irregular

  !> Computes derivative with respect to z at a given point using a downwinding
  !! scheme.  The Z axis uses an irregular grid defined by \iittext{deltas}.
  function dfdz_3d_downwind_irregular(f, i, j, k, deltas)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f  !< field to be derived
    integer, intent(in) :: i                     !< first index at which derivative is taken
    integer, intent(in) :: j                     !< second index at which derivative is taken
    integer, intent(in) :: k                     !< third index at which derivative is taken
    real(dp), dimension(:), intent(in) :: deltas !< vertical sigma levels
    real(dp) :: dfdz_3d_downwind_irregular
    real(dp) :: zkp1MinusZk, zkp2MinusZk, zkp2MinusZkp1
    zkp1MinusZk = deltas(k+1) - deltas(k)
    zkp2MinusZk = deltas(k+2) - deltas(k)
    zkp2MinusZkp1 = deltas(k+2) - deltas(k+1)

    dfdz_3d_downwind_irregular =f(k,   i, j) * (-zkp1MinusZk - zkp2MinusZk)/(zkp1MinusZk * zkp2MinusZk) + &
         f(k+1, i, j) * zkp2MinusZk / (zkp2MinusZkp1 * zkp1MinusZk) - &
         f(k+2, i, j) * zkp1MinusZk / (zkp2MinusZkp1 * zkp2MinusZk)
  end function dfdz_3d_downwind_irregular

  !> Computes derivative with respect to x at the equivalent
  !! point on a staggered grid.
  function dfdx_3d_stag(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in x-direction
    real(dp) :: dfdx_3d_stag
    dfdx_3d_stag = (f(k, i+1, j) + f(k, i+1, j+1) - f(k, i, j) - f(k, i, j+1))/(2*delta) 
  end function dfdx_3d_stag

  !> Computes derivative with respect to y at the equivalent
  !! point on a staggered grid.
  function dfdy_3d_stag(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in y-direction
    real(dp) :: dfdy_3d_stag
    dfdy_3d_stag = (f(k, i, j+1) + f(k, i+1, j+1) - f(k, i, j) - f(k, i+1, j))/(2*delta)
  end function dfdy_3d_stag

  !> Computes derivative with respect to x at the given point
  !! using an upwind method (suitable for maximum boundaries)
  function dfdx_3d_upwind(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in x-direction
    real(dp) :: dfdx_3d_upwind
    dfdx_3d_upwind = (.5 * f(k, i-2, j) - 2 * f(k, i-1, j) + 1.5 * f(k, i, j))/delta
  end function dfdx_3d_upwind

  !> Computes derivative with respect to y at the given point
  !! using an upwind method (suitable for maximum boundaries)
  function dfdy_3d_upwind(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in y-direction
    real(dp) :: dfdy_3d_upwind
    dfdy_3d_upwind = (.5 * f(k, i, j-2) - 2 * f(k, i, j-1) + 1.5 * f(k, i, j))/delta
  end function dfdy_3d_upwind

  !> Computes derivative with respect to x at the given point
  !! using a downwind method (suitable for minimum boundaries)
  function dfdx_3d_downwind(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in x-direction
    real(dp) :: dfdx_3d_downwind
    dfdx_3d_downwind = (-1.5 * f(k, i, j) + 2 * f(k, i+1, j) - .5 * f(k, i+2, j))/delta
  end function dfdx_3d_downwind

  !> Computes derivative with respect to y at the given point
  !! using a downwind method (suitable for minimum boundaries)
  function dfdy_3d_downwind(f, i, j, k, delta)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta               !< grid spacing in y-direction
    real(dp) :: dfdy_3d_downwind
    dfdy_3d_downwind = (-1.5 * f(k, i, j) + 2 * f(k, i, j+1) - .5 * f(k, i, j+2))/delta
  end function dfdy_3d_downwind

  !> Computes derivative fields of the given function.
  !! The z axis is computed on an irregular grid.
  !!
  !! Field containing the direction that derivatives should be upwinded in.
  !! If 0, centered differences are used.  If negative, then upwinded
  !! derivatives (approaching from the negative side) are used.  If
  !! positive, then downwinded derivatives (approaching from the positive
  !! side) are used.
  subroutine df_field_3d(f, deltax, deltay, deltas, out_dfdx, out_dfdy, out_dfdz, &
       direction_x, direction_y)
    implicit none
    real(dp), dimension(:, :, :), intent(in) :: f          !< field to be derived
    real(dp), intent(in) :: deltax                         !< grid spacing in x-direction
    real(dp), intent(in) :: deltay                         !< grid spacing in y-direction
    real(dp), dimension(:), intent(in) :: deltas           !< vertical sigma levels
    real(dp), dimension(:, :, :), intent(out) :: out_dfdx  !< derivatives in x direction
    real(dp), dimension(:, :, :), intent(out) :: out_dfdy  !< derivatives in y direction
    real(dp), dimension(:, :, :), intent(out) :: out_dfdz  !< derivatives in z direction

    real(dp), dimension(:,:), optional :: direction_x      !< x-direction used for upwind derivatives
    real(dp), dimension(:,:), optional :: direction_y      !< y-direction used for upwind derivatives


    integer :: grad_x, grad_y !Sign of the gradient, used for determining upwinding
    integer :: nx, ny, nz, x, y, z
    logical :: upwind

    !Get the size of the field we're working with
    nx = size(f, 2)
    ny = size(f, 3)
    nz = size(f, 1)

    upwind = present(direction_x) .and. present(direction_y)

    !For now, we'll use the function calls defined above.
    !Later on we might want to refactor?
    do x = 1, nx
       do y = 1, ny
          grad_x = 0
          grad_y = 0
          if (upwind) then
             if (direction_x(x,y) < 0 .and. x > 2) then !Upstream case
                grad_x = -1
             else if(direction_x(x,y) > 0 .and. x < nx - 1) then !Downstream case
                grad_x = 1
             end if

             if (direction_y(x,y) < 0 .and. y > 2) then !Upstream case
                grad_y = -1
             else if(direction_y(x,y) > 0 .and. y < ny - 1) then !Downstream case
                grad_y = 1
             end if
          end if

          do z = 1, nz
             !For each of the variables in x, y, check whether or not
             !we need to use an upwinding or downwinding differentiation
             !scheme.
             if (x == 1 .or. grad_x > 0) then
                out_dfdx(z, x, y) = dfdx_3d_downwind(f, x, y, z, deltax)
                !out_dfdx(x, y, z) = (f(x+1,y,z) - f(x,y,z))/deltax
             else if (x == nx .or. grad_x < 0) then
                out_dfdx(z, x, y) = dfdx_3d_upwind(f, x, y, z, deltax)
                !out_dfdx(x, y, z) = (f(x,y,z) - f(x-1,y,z))/deltax
             else
                out_dfdx(z, x, y) = dfdx_3d(f, x, y, z, deltax)
             end if
             if (y == 1 .or. grad_y > 0) then
                out_dfdy(z, x, y) = dfdy_3d_downwind(f, x, y, z, deltay)
                !out_dfdy(x, y, z) = (f(x,y+1,z) - f(x,y,z))/deltay
             else if (y == ny .or. grad_y < 0) then
                out_dfdy(z, x, y) = dfdy_3d_upwind(f, x, y, z, deltay)
                !out_dfdy(x, y, z) = (f(x,y,z) - f(x,y-1,z))/deltay
             else
                out_dfdy(z, x, y) = dfdy_3d(f, x, y, z, deltay)
             end if
             if (z == 1) then
                out_dfdz(z, x, y) = dfdz_3d_downwind_irregular(f, x, y, z, deltas)
             else if (z == nz) then
                out_dfdz(z, x, y) = dfdz_3d_upwind_irregular(f, x, y, z, deltas)
             else
                out_dfdz(z, x, y) = dfdz_3d_irregular(f, x, y, z, deltas)
             end if
          end do
       end do
    end do

  end subroutine df_field_3d

  !> Computes the derivative fields of the given function.  The X and Y
  !! derivatives are computed on a staggered grid.  The Z derivative
  !! is computed on a nonstaggered but irregular grid.  This means that,
  !! if an array of dimensions (n1, n2, n3), the output arrays should
  !! be of size (n1 - 1, n2 - 1, n3)
  subroutine df_field_3d_stag(f, deltax, deltay, deltas, out_dfdx, out_dfdy, out_dfdz)
    implicit none
    real(dp), dimension(:, :, :), intent(in) :: f          !< field to be derived
    real(dp), intent(in) :: deltax                         !< grid spacing in x-direction
    real(dp), intent(in) :: deltay                         !< grid spacing in y-direction
    real(dp), dimension(:), intent(in) :: deltas           !< vertical sigma levels
    real(dp), dimension(:, :, :), intent(out) :: out_dfdx  !< derivatives in x direction
    real(dp), dimension(:, :, :), intent(out) :: out_dfdy  !< derivatives in y direction
    real(dp), dimension(:, :, :), intent(out) :: out_dfdz  !< derivatives in z direction

    real(dp), dimension(4) :: zDerivs !Temporarily holds derivatives in Z to average
    integer :: nx, ny, nz, x, y, z

    !Get the size of the field we're working with
    nx = size(f, 1)
    ny = size(f, 2)
    nz = size(f, 3)

    do x = 1, nx - 1
       do y = 1, ny - 1
          do z = 1, nz
             !We will never have to compute upstream and downstream
             !derivatives in the horizontal (avoided by the staggered scheme),
             !but we will in the vertical.
             out_dfdx(x,y,z) = dfdx_3d_stag(f, x, y, z, deltax)
             out_dfdy(x,y,z) = dfdy_3d_stag(f, x, y, z, deltay)

             !Even though we are not staggering in the vertical, the points
             !we compute the derivatives at are still staggered in the
             !horizontal.  We'll solve this by computing four
             !derivatives horizontally around the point requested
             !and averaging the results
             if (z == 1) then
                zDerivs(1) = dfdz_3d_downwind_irregular(f, x, y, z, deltas)
                zDerivs(2) = dfdz_3d_downwind_irregular(f, x+1, y, z, deltas)
                zDerivs(3) = dfdz_3d_downwind_irregular(f, x, y+1, z, deltas)
                zDerivs(4) = dfdz_3d_downwind_irregular(f, x+1, y+1, z, deltas)
             else if (z == nz) then
                zDerivs(1) = dfdz_3d_upwind_irregular(f, x, y, z, deltas)
                zDerivs(2) = dfdz_3d_upwind_irregular(f, x+1, y, z, deltas)
                zDerivs(3) = dfdz_3d_upwind_irregular(f, x, y+1, z, deltas)
                zDerivs(4) = dfdz_3d_upwind_irregular(f, x+1, y+1, z, deltas)
             else
                zDerivs(1) = dfdz_3d_irregular(f, x, y, z, deltas)
                zDerivs(2) = dfdz_3d_irregular(f, x+1, y, z, deltas)
                zDerivs(3) = dfdz_3d_irregular(f, x, y+1, z, deltas)
                zDerivs(4) = dfdz_3d_irregular(f, x+1, y+1, z, deltas)
             end if
             out_dfdz(x, y, z) = (zDerivs(1) + zDerivs(2) + zDerivs(3) + zDerivs(4)) / 4
          end do
       end do
    end do

  end subroutine df_field_3d_stag

  !------------------------------------------------------------------
  !Second Derivative Estimates, Second Order
  !------------------------------------------------------------------

  !> Computes 2nd derivative with respect to x at the given point
  function d2fdx2_2d(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: d2fdx2_2d     
    d2fdx2_2d = (f(i+1,j) + f(i-1,j) - 2 * f(i, j))/(delta*delta)
  end function d2fdx2_2d

  !> Computes 2nd derivatives with respect to x at the given point
  !! using an downwind method (suitable for maximum boundaries)
  function d2fdx2_2d_downwind(f,i,j,delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: d2fdx2_2d_downwind   

    d2fdx2_2d_downwind = (3*f(i, j) - 7*f(i+1, j) + 5*f(i+2, j) - f(i+3, j)) / (2*delta**2)

  end function d2fdx2_2d_downwind

  !> Computes 2nd derivatives with respect to x at the given point
  !! using an upwind method (suitable for maximum boundaries)
  function d2fdx2_2d_upwind(f,i,j,delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: d2fdx2_2d_upwind 

    d2fdx2_2d_upwind = (3*f(i, j) - 7*f(i-1, j) + 5*f(i-2, j) - f(i-3, j)) / (2*delta**2)

  end function d2fdx2_2d_upwind

  !> Computes 2nd derivatives with respect to y at the given point
  !! using an downwind method (suitable for maximum boundaries)
  function d2fdy2_2d_downwind(f,i,j,delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: d2fdy2_2d_downwind   

    d2fdy2_2d_downwind = (3*f(i, j) - 7*f(i, j+1) + 5*f(i, j+2) - f(i, j+3)) / (2*delta**2)

  end function d2fdy2_2d_downwind

  !> Computes 2nd derivatives with respect to y at the given point
  !! using an upwind method (suitable for maximum boundaries)
  function d2fdy2_2d_upwind(f,i,j,delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: d2fdy2_2d_upwind 

    d2fdy2_2d_upwind = (3*f(i, j) - 7*f(i, j-1) + 5*f(i, j-2) - f(i, j-3)) / (2*delta**2)

  end function d2fdy2_2d_upwind

  !> Computes second derivative with respect to x at the equivalent
  !! point on a staggered grid.
  function d2fdx2_2d_stag(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: d2fdx2_2d_stag

    !This formula can be derived using two central differences
    !(i to i+2, and i-1 to i+1) to get the derivative at
    !i and i+1, then applying a central difference to that
    !in order to get the 2nd derivative at a staggered point
    d2fdx2_2d_stag = sum(f(i+2, j:j+1) + f(i-1, j:j+1) - f(i+1, j:j+1) - f(i, j:j+1))/(4*delta**2)
  end function d2fdx2_2d_stag

  !> Computes second derivative with respect to x at the given point
  !! using a downwind method (suitable for minimum boundaries)
  function d2fdx2_2d_stag_downwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: d2fdx2_2d_stag_downwind

    d2fdx2_2d_stag_downwind = sum(3*f(i, j:j+1) - 7*f(i+1, j:j+1) + 5*f(i+2, j:j+1) - f(i+3, j:j+1)) / (4*delta**2)
  end function d2fdx2_2d_stag_downwind

  !> Computes second derivative with respect to x at the given point
  !! using a upwind method (suitable for minimum boundaries)
  function d2fdx2_2d_stag_upwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in x-direction
    real(dp) :: d2fdx2_2d_stag_upwind

    d2fdx2_2d_stag_upwind = sum(-3*f(i+1, j:j+1) + 7*f(i, j:j+1) - 5*f(i-1, j:j+1) + f(i-2, j:j+1)) / (4*delta**2)
  end function d2fdx2_2d_stag_upwind

  !> Computes 2nd derivative with respect to y at the given point
  function d2fdy2_2d(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: d2fdy2_2d
    d2fdy2_2d = (f(i, j+1) + f(i, j-1) - 2 * f(i, j))/(delta*delta)
  end function d2fdy2_2d

  !> Computes second derivative with respect to y at the equivalent
  !! point on a staggered grid.
  function d2fdy2_2d_stag(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: d2fdy2_2d_stag

    !This formula can be derived using two central differences
    !(i to i+2, and i-1 to i+1) to get the derivative at
    !i and i+1, then applying a central difference to that
    !in order to get the 2nd derivative at a staggered point
    d2fdy2_2d_stag = sum(f(i:i+1, j+2) + f(i:i+1, j-1) - f(i:i+1, j+1) - f(i:i+1, j))/(4*delta**2)
  end function d2fdy2_2d_stag

  !> Computes second derivative with respect to y at the given point
  !! using a downwind method (suitable for minimum boundaries)
  function d2fdy2_2d_stag_downwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: d2fdy2_2d_stag_downwind

    d2fdy2_2d_stag_downwind = sum(3*f(i:i+1, j) - 7*f(i:i+1, j+1) + 5*f(i:i+1, j+2) - f(i:i+1, j+3)) / (4*delta**2)
  end function d2fdy2_2d_stag_downwind

  !> Computes second derivative with respect to y at the given point
  !! using a upwind method (suitable for minimum boundaries)
  function d2fdy2_2d_stag_upwind(f, i, j, delta)
    implicit none
    real(dp), dimension(:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                  !< first index at which derivative is taken
    integer, intent(in) :: j                  !< second index at which derivative is taken
    real(dp), intent(in) :: delta             !< grid spacing in y-direction
    real(dp) :: d2fdy2_2d_stag_upwind

    d2fdy2_2d_stag_upwind = sum(-3*f(i:i+1, j+1) + 7*f(i:i+1, j) - 5*f(i:i+1, j-1) + f(i:i+1, j-2)) / (4*delta**2)
  end function d2fdy2_2d_stag_upwind


  !> Computes second derivative fields of the given function.
  subroutine d2f_field(f, deltax, deltay, d2fdx2, d2fdy2, direction_x, direction_y)
    implicit none 
    real(dp), dimension(:, :), intent(in) :: f                      !< field to be derived
    real(dp), intent(in) :: deltax                                  !< grid spacing in x direction
    real(dp), intent(in) :: deltay                                  !< grid spacing in y direction
    real(dp), dimension(:, :), intent(out) :: d2fdx2                !< derivatives in x direction
    real(dp), dimension(:, :), intent(out) ::    d2fdy2             !< derivatives in y direction
    real(dp), dimension(:, :), intent(in), optional  :: direction_x !< x-direction used for upwind derivatives
    real(dp), dimension(:, :), intent(in), optional  :: direction_y !< y-direction used for upwind derivatives

    integer :: i,j

    do i = 1,size(f,1)
       do j = 1,size(f,2)
          !non-staggered grid
          if (i == 1) then
             d2fdx2(i,j) = d2fdx2_2d_downwind(f,i,j,deltax)
          else if (i == size(f,1)) then
             d2fdx2(i,j) = d2fdx2_2d_upwind(f,i,j,deltax)
          else
             if (present(direction_x)) then
                if (direction_x(i,j) > 0) then
                   d2fdx2(i,j) = d2fdx2_2d_downwind(f,i,j,deltax)
                else if (direction_x(i,j) < 0) then
                   d2fdx2(i,j) = d2fdx2_2d_upwind(f,i,j,deltax)
                else
                   d2fdx2(i,j) = d2fdx2_2d(f,i,j,deltax)
                end if
             else
                d2fdx2(i,j) = d2fdx2_2d(f,i,j,deltax)
             end if
          end if

          if (j == 1) then
             d2fdy2(i,j) = d2fdy2_2d_downwind(f,i,j,deltax)
          else if (j == size(f,2)) then
             d2fdy2(i,j) = d2fdy2_2d_upwind(f,i,j,deltax)
          else
             if (present(direction_y)) then
                if (direction_y(i,j) > 0) then
                   d2fdy2(i,j) = d2fdy2_2d_downwind(f,i,j,deltax)
                else if (direction_y(i,j) < 0) then
                   d2fdy2(i,j) = d2fdy2_2d_upwind(f,i,j,deltax)
                else
                   d2fdy2(i,j) = d2fdy2_2d(f,i,j,deltax)
                end if
             else
                d2fdy2(i,j) = d2fdy2_2d(f,i,j,deltax)
             end if
          end if
       end do
    end do
  end subroutine d2f_field

  
  !> Computes derivative fields of the given function.  Places the result
  !! on a staggered grid.  If periodic in one dimension is set, that 
  !! dimension for derivatives must be the same size as the value's dimension.
  !! Otherwise, it should be one less
  !! \todo Rewrite this using the existing derivative machinery
  !! \todo enable periodic boundary conditions
  subroutine d2f_field_stag(f, deltax, deltay, d2fdx2, d2fdy2, periodic_x, periodic_y)
    implicit none 

    real(dp), dimension(:, :), intent(in) :: f                      !< field to be derived
    real(dp), intent(in) :: deltax                                  !< grid spacing in x direction
    real(dp), intent(in) :: deltay                                  !< grid spacing in y direction
    real(dp), dimension(:, :), intent(out) :: d2fdx2                !< derivatives in x direction
    real(dp), dimension(:, :), intent(out) ::    d2fdy2             !< derivatives in y direction
    logical :: periodic_x                                           !< flag for periodic BC in x
    logical :: periodic_y                                           !< flag for periodic BC in y

    real(dp) :: dewsq4, dnssq4
    integer :: ew,ns

    integer :: pt(2)
    integer :: nsn
    integer :: ewn

    nsn = size(f,2)
    ewn = size(f,1)

    dewsq4 = 4.0d0 * deltax * deltax
    dnssq4 = 4.0d0 * deltay * deltay

    d2fdx2 = 0
    d2fdy2 = 0

    do ns = 2, nsn-2
       do ew = 2, ewn-2
          d2fdx2(ew,ns) = centerew(ew,ns)
          d2fdy2(ew,ns) = centerns(ew,ns)
       end do
    end do

    ! *** 2nd order boundaries using upwinding

    do ew = 1, ewn-1, ewn-2

       pt = whichway(ew)

       do ns = 2, nsn-2 
          d2fdx2(ew,ns) = boundyew(pt,ns)
          d2fdy2(ew,ns) = centerns(ew,ns)
       end do

    end do

    do ns = 1, nsn-1, nsn-2

       pt = whichway(ns)

       do ew = 2, ewn-2
          d2fdx2(ew,ns) = centerew(ew,ns)
          d2fdy2(ew,ns) = boundyns(pt,ew)
       end do

    end do

    do ns = 1, nsn-1, nsn-2
       do ew = 1, ewn-1, ewn-2
          pt = whichway(ew)
          d2fdx2(ew,ns) = boundyew(pt,ns)
          pt = whichway(ns)
          d2fdy2(ew,ns) = boundyns(pt,ew)
       end do
    end do

  contains

    function centerew(ew,ns)

      implicit none

      real(dp) :: centerew
      integer ns,ew

      centerew = (sum(f(ew+2,ns:ns+1)) + sum(f(ew-1,ns:ns+1)) - &
           sum(f(ew+1,ns:ns+1)) - sum(f(ew,ns:ns+1))) / dewsq4

    end function centerew

    function centerns(ew,ns)

      implicit none

      real(dp) :: centerns
      integer ns,ew

      centerns = (sum(f(ew:ew+1,ns+2)) + sum(f(ew:ew+1,ns-1)) - &
           sum(f(ew:ew+1,ns+1)) - sum(f(ew:ew+1,ns))) / dnssq4

    end function centerns

    function boundyew(pt,ns)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyew
      integer ns

      boundyew = pt(1) * (3.0d0 * sum(f(pt(2),ns:ns+1)) - 7.0d0 * sum(f(pt(2)+pt(1),ns:ns+1)) + &
           5.0d0 * sum(f(pt(2)+2*pt(1),ns:ns+1)) - sum(f(pt(2)+3*pt(1),ns:ns+1))) / dewsq4

    end function boundyew

    function boundyns(pt,ew)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyns
      integer ew

      boundyns = pt(1) * (3.0d0 * sum(f(ew:ew+1,pt(2))) - 7.0d0 * sum(f(ew:ew+1,pt(2)+pt(1))) + &
           5.0d0 * sum(f(ew:ew+1,pt(2)+2*pt(1))) - sum(f(ew:ew+1,pt(2)+3*pt(1)))) / dnssq4

    end function boundyns

    function whichway(i)

      implicit none

      integer, intent(in) :: i
      integer :: whichway(2) 

      if (i == 1) then 
         whichway = (/1,1/)
      else
         whichway = (/-1,i+1/)
      end if

    end function whichway

    !       real(dp), dimension(:,:), intent(in) :: f
    !       real(dp), dimension(:,:), intent(out) :: d2fdx2, d2fdy2
    !       real(dp), intent(in) :: deltax, deltay
    !       logical :: periodic_x, periodic_y
    !       
    !       integer :: nx, x, ny, y
    !       
    !       nx = size(f, 1)
    !       ny = size(f, 2)
    !       
    !       !NOTE: See the field 1st derivative staggered function for
    !       !a discussion of periodic boundary conditions
    !       
    !       !First compute the values that do not fall on any boundaries
    !       !This is the same regardless of whether periodic boundary
    !       !conditions are used
    !       do x = 1, nx-1
    !             do y = 1, ny-1
    !                 if (x == 1) then
    !                     d2fdx2(1,y) = d2fdx2_2d_stag_downwind(f, 1, y, deltax)
    !                 else if (x == nx - 1) then
    !                     d2fdx2(nx-1, y) = d2fdx2_2d_stag_upwind(f, nx-1, y, deltax)
    !                 else
    !                     d2fdx2(x,y) = d2fdx2_2d_stag(f, x, y, deltax)
    !                 end if
    !                 
    !                 if (y == 1) then
    !                     d2fdy2(x,1) = d2fdy2_2d_stag_downwind(f, x, 1, deltay)
    !                 else if (y == ny - 1) then
    !                     d2fdy2(x, ny-1) = d2fdy2_2d_stag_upwind(f, x, ny-1, deltay)
    !                 else
    !                     d2fdy2(x,y) = d2fdy2_2d_stag(f, x, y, deltay)
    !                 end if
    !               end do
    !       end do
    !       
    !       !If we are not using periodic boundary conditions, then we need
    !       !to use an upwinding scheme to get the values when x = 1, y = 1,
    !       !x = nx - 1, or y = ny - 1
    !               !If we are using periodic boundary conditions, then compute the
    !       !boundaries with input from the periodic conditions.  We do not
    !       !upwind or downwind.  Also, because an extra set of values around
    !       !the edges is necessary to correctly maintain periodicity,
    !       !we fill in values where x = nx and where y = ny (whereas we
    !       !do not with nonperiodic boundaries, as the staggered grid
    !       !points fall strictly in the interior of the nonstaggered
    !       !grid)
    !       do y = 1, ny - 2
    !               if (.not.(periodic_x)) then
    !                       d2fdx2(1,y) = d2fdx2_2d_stag_downwind(f, 1, y, deltax)
    ! 
    !                       d2fdx2(nx-1, y) = d2fdx2_2d_stag_upwind(f, nx-1, y, deltax)
    !                       
    !                   else
    !                       !Because of the periodicity, I will simply copy the appropriate values
    !                       !(e.g. u(1) = u(n-2), u(n-1) = u(2)
    !                       d2fdx2(1,y) = d2fdx2(nx-2,y)
    !                       d2fdx2(nx-1,y) = d2fdx2(2,y)
    !               end if
    !               d2fdy2(1,y) = d2fdy2_2d_stag(f, 1, y, deltay)    
    !               d2fdy2(nx-1, y) = d2fdy2_2d_stag(f, nx-1, y, deltay)    
    !       end do
    !       
    !       !See comments for the periodic x boundary case above; the same
    !       !principles apply here.
    !       do x=1, nx-2
    !               if (.not.(periodic_y)) then
    !                       d2fdy2(x,1) = d2fdy2_2d_stag_downwind(f, x, 1, deltay)
    !                       d2fdy2(x, ny-1) = d2fdy2_2d_stag_upwind(f, x, ny-1, deltay)
    !               else
    !                               d2fdy2(x,1) = d2fdy2(x,ny-2)
    !                               d2fdy2(x,nx-1) = d2fdy2(x,2)
    !               end if
    !               d2fdx2(x,1) = d2fdx2_2d_stag(f, x, 1, deltax)
    !               d2fdx2(x, ny-1) = d2fdx2_2d_stag(f, x, ny-1, deltax)
    !       end do
    !       
    !       
    !       !TODO: CHange this to use the scheme above
    !       !We have neglected so far to take care of the four points that occur at the
    !       !maximum two indices in x and y.  If no periodic boundaries are being used,
    !       !we compute the value zt (nx-1, ny-1) using upwinding schemes.
    !       if (.not. periodic_x .and. .not. periodic_y) then
    !               d2fdx2(nx-1, ny-1) = d2fdx2_2d_stag_upwind(f, nx-1, ny-1, deltax)
    !               d2fdy2(nx-1, ny-1) = d2fdy2_2d_stag_upwind(f, nx-1, ny-1, deltay)
    !       else if (.not. periodic_x) then
    !               !y is periodic - this means we need to compute the derivative
    !               !for x=nx-1 and y=ny, ny-1.  We will copy and paste
    !               !y derivatives (for now), as above, and upwind
    !               !the x derivatives
    !               d2fdx2(nx-1, ny-1) = d2fdx2_2d_stag_upwind(f, nx-1, ny-1, deltax)
    !               d2fdy2(nx-1, ny-1) = sum(f(nx-1:nx, 1) + f(nx-1:nx, ny-2) - f(nx-1:nx, ny) - f(nx-1:nx, ny-1))/(4*deltay**2)
    !               
    !               
    !               d2fdx2(nx-1, ny)  = d2fdx2_2d_stag_upwind(f, nx-1, ny, deltax)
    !               d2fdy2(nx-1, ny)  = sum(f(nx-1:nx, 2) + f(nx-1:nx, ny-1) - f(nx-1:nx, 1) - f(nx-1:nx, ny))/(4*deltay**2)
    !               
    !       else if (.not. periodic_y) then
    !               !See comments for the periodic y case above - we are basically using the same
    !               !logic with x and y swapped
    !               d2fdx2(nx-1, ny-1) = sum(f(1, ny-1:ny) + f(nx-2, ny-1:ny) - f(nx, ny-1:ny) - f(nx-1, ny-1:ny))/(4*deltax**2)
    !               d2fdy2(nx-1, ny-1) = d2fdy2_2d_stag_upwind(f, nx-1, ny-1, deltay)
    !               
    !               d2fdx2(nx, ny-1) = sum(f(2, ny-1:ny) + f(nx-1, ny-1:ny) - f(1, ny-1:ny) - f(nx, ny-1:ny))/(4*deltax**2)
    !               d2fdy2(nx, ny-1) = d2fdy2_2d_stag_upwind(f, nx-1, ny-1, deltay)
    !       else
    !               !X and Y are periodic; we will use the periodic forms of the above differences
    !               !Some of these will get very funky because 
    !                       d2fdx2(nx-1, ny-1) = sum(f(1, ny-1:ny) + f(nx-1, ny-1:ny) - f(nx, ny-1:ny) - f(nx-1, ny-1:ny))/(4*deltax**2)
    !                       d2fdy2(nx-1, ny-1) = sum(f(nx-1:nx, 1) + f(nx-1:nx, ny-2) - f(nx-1:nx, ny) - f(nx-1:nx, ny-1))/(4*deltay**2)
    !               
    !               d2fdx2(nx, ny-1) = sum(f(2, ny-1:ny) + f(nx-1, ny-1:ny) - f(1, ny-1:ny) - f(nx, ny-1:ny))/(4*deltax**2)
    !               d2fdy2(nx, ny-1) = ((f(nx, 1) + f(nx, ny-2) - f(nx, ny) - f(nx, ny-1)) + &
    !                                    (f(1, 1) + f(1, ny-2) - f(1, ny) - f(1, ny-1)))/(4*deltay**2)
    !               
    !               d2fdy2(nx-1, ny)  = ((f(1, ny) + f(nx-1, ny) - f(nx, ny) - f(nx-1, ny)) + &
    !                                    (f(1, 1) + f(nx-1, 1) - f(nx, 1) - f(nx-1, 1)))/(4*deltax**2)
    !                       d2fdy2(nx-1, ny)  = sum(f(nx-1:nx, 2) + f(nx-1:nx, ny-1) - f(nx-1:nx, 1) - f(nx-1:nx, ny))/(4*deltay**2)
    !                       
    !                       d2fdx2(nx, ny) = ((f(2, ny) + f(nx-1, ny) - f(1, ny) - f(nx, ny)) + (f(2, 1) + f(nx-1, 1) - f(1, 1) - f(nx, 1))) / (4*deltax**2)
    !                       d2fdy2(nx, ny) = ((f(nx, 2) + f(nx, ny-1) - f(nx, 1) - f(nx, ny)) + (f(1, 2) + f(1, ny-1) - f(1, 1) - f(1, ny)))/(4*deltay**2)
    !               
    !       end if

  end subroutine d2f_field_stag

  !> Computes derivative taken first w.r.t x, then to y at the given point.
  function d2fdxy_3d(f, i, j, k, delta_x, delta_y)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f !< field to be derived
    integer, intent(in) :: i                    !< first index at which derivative is taken
    integer, intent(in) :: j                    !< second index at which derivative is taken
    integer, intent(in) :: k                    !< third index at which derivative is taken
    real(dp), intent(in) :: delta_x             !< grid spacing in x direction
    real(dp), intent(in) :: delta_y             !< grid spacing in y direction

    real(dp) :: d2fdxy_3d

    d2fdxy_3d = (f(k, i-1, j-1) - f(k, i-1, j+1) - f(k, i+1, j-1) + f(k, i+1, j+1))/(4*delta_x*delta_y) 
  end function d2fdxy_3d

  !> Computes derivative taken first w.r.t x, then to z at the given point.
  function d2fdxz_3d(f, i, j, k, delta_x, deltas)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f  !< field to be derived
    integer, intent(in) :: i                     !< first index at which derivative is taken
    integer, intent(in) :: j                     !< second index at which derivative is taken
    integer, intent(in) :: k                     !< third index at which derivative is taken
    real(dp), intent(in) :: delta_x              !< grid spacing in x direction
    real(dp), dimension(:), intent(in) :: deltas !< vertical sigma levels
    real(dp) :: d2fdxz_3d

    d2fdxz_3d = (.5/delta_x) * ( &
         (f(k-1, i+1, j) - f(k-1, i-1, j)) * (deltas(k) - deltas(k+1)) / ( (deltas(k) - deltas(k-1)) * (deltas(k+1) - deltas(k-1)) ) + &
         (f(k,   i+1, j) - f(k,   i-1, j)) * (deltas(k+1) + deltas(k-1) - 2*deltas(k)) / ( (deltas(k) - deltas(k-1)) * (deltas(k+1) - deltas(k)) ) + &
         (f(k+1, i+1, j) - f(k+1, i-1, j)) * (deltas(k) - deltas(k-1)) / ( (deltas(k+1) - deltas(k)) * (deltas(k+1) - deltas(k-1)) ) )
  end function d2fdxz_3d

  !> Computes derivative taken first w.r.t y, then to z at the given point.
  function d2fdyz_3d(f, i, j, k, delta_y, deltas)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f  !< field to be derived
    integer, intent(in) :: i                     !< first index at which derivative is taken
    integer, intent(in) :: j                     !< second index at which derivative is taken
    integer, intent(in) :: k                     !< third index at which derivative is taken
    real(dp), intent(in) :: delta_y              !< grid spacing in y direction
    real(dp), dimension(:), intent(in) :: deltas !< vertical sigma levels
    real(dp) :: d2fdyz_3d

    d2fdyz_3d = (.5/delta_y) * ( &
         (f(k-1, i, j+1) - f(k-1, i, j-1)) * (deltas(k) - deltas(k+1)) / ( (deltas(k) - deltas(k-1)) * (deltas(k+1) - deltas(k-1)) ) + &
         (f(k,   i, j+1) - f(k,   i, j-1)) * (deltas(k+1) + deltas(k-1) - 2*deltas(k)) / ( (deltas(k) - deltas(k-1)) * (deltas(k+1) - deltas(k)) ) + &
         (f(k+1, i, j+1) - f(k+1, i, j-1)) * (deltas(k) - deltas(k-1)) / ( (deltas(k+1) - deltas(k)) * (deltas(k+1) - deltas(k-1)) ) )
  end function d2fdyz_3d

  !> Computes second derivative with respect to z at a given point
  !! where the Z axis uses an irregular grid defined by \ittext{deltas}.
  !! This derivative is given by the formula:
  function d2fdz2_3d_irregular(f, i, j, k, deltas)
    implicit none
    real(dp), dimension(:,:,:), intent(in) :: f  !< field to be derived
    integer, intent(in) :: i                     !< first index at which derivative is taken
    integer, intent(in) :: j                     !< second index at which derivative is taken
    integer, intent(in) :: k                     !< third index at which derivative is taken
    real(dp), dimension(:), intent(in) :: deltas !< vertical sigma levels
    real(dp) :: d2fdz2_3d_irregular
    real(dp) :: zkMinusZkp1, zkMinusZkm1, zkp1MinusZkm1, zkp1MinusZk

    zkMinusZkp1 = deltas(k) - deltas(k+1)
    zkMinusZkm1 = deltas(k) - deltas(k-1)
    zkp1MinusZkm1 = deltas(k+1) - deltas(k-1)
    zkp1MinusZk = -1 * zkMinusZkp1


    d2fdz2_3d_irregular = 2 * f(k-1, i, j) / (zkMinusZkm1 * zkp1MinusZkm1) - &
         2 * f(k,   i, j) / (zkp1MinusZk * zkMinusZkm1) + &
         2 * f(k+1, i, j) / (zkp1Minuszk * zkp1MinusZkm1)    
  end function d2fdz2_3d_irregular

end module glimmer_deriv
