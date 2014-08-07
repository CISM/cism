!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grid_operators.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! Various grid operators for glissade dycore, including routines for computing gradients
! and switching between staggered and unstaggered grids

module glissade_grid_operators

    use glimmer_global, only: dp
    use parallel  ! for debugging only

    implicit none

    private
    public :: glissade_stagger, glissade_unstagger,  &
              glissade_centered_gradient, glissade_upstream_gradient

    logical, parameter :: verbose_gradient = .false.

contains

!----------------------------------------------------------------------------

  subroutine glissade_stagger(nx,           ny,        &
                              var,          stagvar,   &
                              imask,        stag_flag_in)

    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    ! to find values on the staggered grid (dimension nx-1, ny-1).

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! unstaggered field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       imask                    ! = 1 where values are included in the average, else = 0
                                ! Typically imask = 1 where ice is present (or thck > thklim), else = 0

    integer, intent(in), optional ::   &
       stag_flag_in             ! 0 = use all values when interpolating (including zeroes where ice is absent)
                                ! 1 = use only values where ice is present

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stag_flag

    if (present(stag_flag_in)) then
       stag_flag = stag_flag_in
    else
       stag_flag = 1  ! default is to average only over the cells with ice present
    endif

    stagvar(:,:) = 0.d0

    if (stag_flag == 0) then

       ! Average over all four neighboring cells

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          stagvar(i,j) = (var(i,j+1) + var(i+1,j+1) + var(i,j) + var(i+1,j)) / 4.d0	  
       enddo
       enddo  

    elseif (stag_flag==1) then

       ! Average over cells with ice present (imask = 1)

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = imask(i,j+1)*var(i,j+1) + imask(i+1,j+1)*var(i+1,j+1)  &
                 + imask(i,j)  *var(i,j)   + imask(i+1,j)  *var(i+1,j)	  
          summask = real(imask(i,j+1) + imask(i+1,j+1) + imask(i,j) + imask(i+1,j), dp)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

  end subroutine glissade_stagger

!----------------------------------------------------------------------------

  subroutine glissade_unstagger(nx,           ny,          &
                                stagvar,      unstagvar,   &
                                vmask,        stag_flag_in)

    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    ! to find values on the staggered grid (dimension nx-1, ny-1).

    use parallel, only: parallel_halo

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    real(dp), dimension(nx,ny), intent(out) ::    &
       unstagvar                ! unstaggered field, defined at cell centers

    integer, dimension(nx-1,ny-1), intent(in) ::        &
       vmask                    ! = 1 for vertices where the value is used in the average, else = 0
                                ! Note: The user needs to compute this mask in the calling subroutine.
                                !       It will likely be based on the scalar ice mask, but the details are left open.

    integer, intent(in), optional ::   &
       stag_flag_in             ! 0 = use all values when interpolating
                                ! 1 = use only values where vmask = 1

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stag_flag

    if (present(stag_flag_in)) then
       stag_flag = stag_flag_in
    else
       stag_flag = 1  ! default is to average over cells where vmask = 1
    endif

    unstagvar(:,:) = 0.d0

    if (stag_flag == 0) then

       ! Average over all four neighboring cells

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          unstagvar(i,j) = (stagvar(i,j) + stagvar(i-1,j) + stagvar(i,j-1) + stagvar(i-1,j-1)) / 4.d0	  
       enddo
       enddo  

    elseif (stag_flag==1) then

       ! Average over cells with vmask = 1

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          sumvar = vmask(i-1,j)  *stagvar(i-1,j)   + vmask(i,j)  *stagvar(i,j)  &
                 + vmask(i-1,j-1)*stagvar(i-1,j-1) + vmask(i,j-1)*stagvar(i,j-1)  
          summask = real(vmask(i-1,j) + vmask(i,j) + vmask(i-1,j-1) + vmask(i,j-1), dp)
          if (summask > 0.d0) unstagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

    ! Fill in halo values
    call parallel_halo(unstagvar)

  end subroutine glissade_unstagger

!****************************************************************************

  subroutine glissade_centered_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        f,                       &
                                        df_dx,        df_dy,     &
                                        imask,        grad_flag_in)

    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient is evaluated at the four neighboring points and is second-order accurate.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       imask                    ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       grad_flag_in             ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: if one or more values is masked out, set df_dx = df_dy = 0

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    integer :: summask, grad_flag

    !   Gradient at vertex(i,j) is based on f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)


    if (present(grad_flag_in)) then
       grad_flag = grad_flag_in
    else
       grad_flag = 2  ! TODO - Think about the best default
    endif

    if (grad_flag == 0) then

       ! compute gradient as in Glide, using all four scalar values
       ! (provided ice is present in at least one cell)

       do j = 1, ny-1
       do i = 1, nx-1

          summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)

          if (summask >= 1) then
             df_dx(i,j) = (f(i+1,j) + f(i+1,j+1) - f(i,j) - f(i,j+1)) / (2.d0 * dx)
             df_dy(i,j) = (f(i,j+1) + f(i+1,j+1) - f(i,j) - f(i+1,j)) / (2.d0 * dy)
          else
             df_dx(i,j) = 0.d0
             df_dy(i,j) = 0.d0
          endif

       enddo     ! i
       enddo     ! j

    elseif (grad_flag==1) then

       ! set gradient to zero at ice margin

       do j = 1, ny-1
       do i = 1, nx-1

          summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)

          if (summask == 4) then
             df_dx(i,j) = (f(i+1,j) + f(i+1,j+1) - f(i,j) - f(i,j+1)) / (2.d0 * dx)
             df_dy(i,j) = (f(i,j+1) + f(i+1,j+1) - f(i,j) - f(i+1,j)) / (2.d0 * dy)
          else
             df_dx(i,j) = 0.d0
             df_dy(i,j) = 0.d0
          endif

       enddo     ! i
       enddo     ! j

    elseif (grad_flag==2) then

       ! compute gradient using neighbor values that are available

       do j = 1, ny-1
       do i = 1, nx-1

          summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)

          if (summask == 4) then
             df_dx(i,j) = (f(i+1,j) + f(i+1,j+1) - f(i,j) - f(i,j+1)) / (2.d0 * dx)
             df_dy(i,j) = (f(i,j+1) + f(i+1,j+1) - f(i,j) - f(i+1,j)) / (2.d0 * dy)
          else
             ! df_dx
             if (imask(i,j)==1 .and. imask(i+1,j)==1) then
                df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
             elseif (imask(i,j+1)==1 .and. imask(i+1,j+1)==1) then
                df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
             else
                df_dx(i,j) = 0.d0
             endif

             ! df_dy
             if (imask(i,j)==1 .and. imask(i,j+1)==1) then
                df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
             elseif (imask(i+1,j)==1 .and. imask(i+1,j+1)==1) then
                df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
             else
                df_dy(i,j) = 0.d0
             endif

          endif

       enddo    ! i
       enddo    ! j

    endif        ! grad_flag

!WHL - debug
    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Centered gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f6.0)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
          do i = 1, nx-1
             write(6,'(f6.0)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_centered_gradient

!****************************************************************************

  subroutine glissade_upstream_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        f,                       &
                                        df_dx,        df_dy,     &
                                        imask,                   &
                                        accuracy_flag_in,        &
                                        grad_flag1_in,           & 
                                        grad_flag2_in)

    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    !  compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient can be evaluated at two upstream points (for first-order accuracy) 
    !  or at four upstream points (for second-order accuracy).
    ! Note: Upstream is defined by the direction of higher surface elevation
    !  rather than the direction the flow is coming from (though these are
    !  usually the same).
    !
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       f                        ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       imask                    ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       accuracy_flag_in         ! = 1 for 1st order, 2 for 2nd order

    integer, intent(in), optional ::    &
       grad_flag1_in            ! options for 1st order upstream (2 upstream points) 
                                ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: if one or more values is masked out, set df_dx = df_dy = 0
                                ! 2: use values that are available, else set gradient components to 0

    integer, intent(in), optional ::    &
       grad_flag2_in            ! options for 2nd order upstream (4 upstream points)
                                ! 0: 2nd order gradient using all values
                                ! 1: mix 1st order and 2nd order based on relative magnitude of terms
                                ! 2: mix 1st order and 2nd order based on agreement with centered difference

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sum1, sum2
    real(dp) :: diff1, diff2
    integer :: grad_flag1, grad_flag2, accuracy_flag, summask

    real(dp) :: dfdx_c, dfdx_1, dfdx_2
    real(dp) :: dfdy_c, dfdy_1, dfdy_2

    !   Gradient at vertex(i,j) is based on two points out of f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)


    if (present(accuracy_flag_in)) then
       accuracy_flag = accuracy_flag_in
    else
       accuracy_flag = 2   ! default to second-order
    endif

    if (present(grad_flag1_in)) then
       grad_flag1 = grad_flag1_in
    else
       grad_flag1 = 0   ! default is to include all cells in gradient, including ice-free cells
    endif

    if (present(grad_flag2_in)) then
       grad_flag2 = grad_flag2_in
    else
       grad_flag2 = 0   ! default is to use 2nd order gradient everywhere
    endif

    df_dx(:,:) = 0.d0
    df_dy(:,:) = 0.d0

    if (accuracy_flag == 1) then   ! first-order accurate

       !WHL - For the Halfar SIA test, none of the first-order options
       !      are as accurate as the second-order options.
       !      But I've kept the first-order options for completeness.

       if (grad_flag1 == 0) then  ! use all points, including ice-free points

          do j = 1, ny-1
          do i = 1, nx-1

             ! Compute df_dx by taking upstream gradient

             sum1 = f(i+1,j+1) + f(i,j+1)
             sum2 = f(i+1,j) + f(i,j)
             
             if (sum1 > sum2) then
                df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
             else
                df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
             endif

             ! Compute df_dy by taking upstream gradient

             sum1 = f(i+1,j+1) + f(i+1,j)
             sum2 = f(i,j+1) + f(i,j)

             if (sum1 > sum2) then
                df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
             else
                df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
             endif

          enddo     ! i
          enddo     ! j

       elseif (grad_flag1 == 1) then  ! gradient = 0 if any neighbor cell is ice-free

          do j = 1, ny-1
          do i = 1, nx-1

             summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)

             if (summask == 4) then

                ! Compute df_dx by taking upstream gradient

                sum1 = f(i+1,j+1) + f(i,j+1)
                sum2 = f(i+1,j) + f(i,j)

                if (sum1 > sum2) then
                   df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
                else
                   df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
                endif

                ! Compute df_dy by taking upstream gradient
                
                sum1 = f(i+1,j+1) + f(i+1,j)
                sum2 = f(i,j+1) + f(i,j)

                if (sum1 > sum2) then
                   df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
                else
                   df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
                endif

             else

                df_dx(i,j) = 0.d0
                df_dy(i,j) = 0.d0

             endif

          enddo    ! i
          enddo    ! j

       elseif (grad_flag1 == 2) then  ! gradient = 0 if either or both upstream points are ice-free

          do j = 1, ny-1
          do i = 1, nx-1

             ! Compute df_dx by taking upstream gradient

             sum1 = f(i+1,j+1) + f(i,j+1)
             sum2 = f(i+1,j) + f(i,j)

             if (sum1 > sum2 .and. imask(i+1,j+1)==1 .and. imask(i,j+1)==1) then
                df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
             elseif (sum1 <= sum2 .and. imask(i+1,j)==1 .and. imask(i,j)==1) then
                df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
             else
                df_dx(i,j) = 0.d0
             endif

             ! Compute df_dy by taking upstream gradient

             sum1 = f(i+1,j+1) + f(i+1,j)
             sum2 = f(i,j+1) + f(i,j)
             
             if (sum1 > sum2 .and. imask(i+1,j+1)==1 .and. imask(i+1,j)==1) then
                df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
             elseif (sum1 <= sum2 .and. imask(i,j+1)==1 .and. imask(i,j)==1) then
                df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
             else
                df_dy(i,j) = 0.d0
             endif

          enddo     ! i
          enddo     ! j

       endif   ! grad_flag1

    else    ! second-order accurate

       !WHL: I tested all three values of grad_flag2 for the Halfar SIA test
       !     problem.  Results were similar for all three, but I opted for
       !     grad_flag2 = 0 as the default because it is simplest.

       if (grad_flag2 == 0) then  ! use all points, including ice-free points

          do j = 2, ny-2
          do i = 2, nx-2

             summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)
           
             if (summask > 0) then
                ! Compute df_dx by taking upstream gradient

                sum1 = f(i+1,j+1) + f(i,j+1) + f(i+1,j+2) + f(i,j+2)
                sum2 = f(i+1,j) + f(i,j) + f(i+1,j-1) + f(i,j-1)
             
                if (sum1 > sum2) then
                   df_dx(i,j) = (1.5d0 * (f(i+1,j+1) - f(i,j+1))     &
                               - 0.5d0 * (f(i+1,j+2) - f(i,j+2))) / dx
                else
                   df_dx(i,j) = (1.5d0 * (f(i+1,j) - f(i,j))         &
                               - 0.5d0 * (f(i+1,j-1) - f(i,j-1))) / dx
                endif

                ! Compute df_dy by taking upstream gradient

                sum1 = f(i+1,j+1) + f(i+1,j) + f(i+2,j+1) + f(i+2,j)
                sum2 = f(i,j+1) + f(i,j) + f(i-1,j+1) + f(i-1,j)
             
                if (sum1 > sum2) then
                   df_dy(i,j) = (1.5d0 * (f(i+1,j+1) - f(i+1,j))     &
                               - 0.5d0 * (f(i+2,j+1) - f(i+2,j))) / dy
                else
                   df_dy(i,j) = (1.5d0 * (f(i,j+1) - f(i,j))         &
                               - 0.5d0 * (f(i-1,j+1) - f(i-1,j))) / dy
                endif

             endif  ! summask

          enddo     ! i
          enddo     ! j

       elseif (grad_flag2 == 1) then  ! second order in most places, reduces to first order in regions of strong curvature

          do j = 2, ny-2   ! loop does not include all of halo
          do i = 2, nx-2

             summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)

             if (summask > 0) then

                ! Compute df_dx

                sum1 = f(i+1,j+1) + f(i,j+1)
                sum2 = f(i+1,j) + f(i,j)
             
                if (sum1 > sum2) then  ! increasing j is upstream
                   diff1 = f(i+1,j+1) - f(i,j+1)
                   diff2 = f(i+1,j+2) - f(i,j+2)
                   if (abs(diff2) < 3.d0*abs(diff1)) then   ! second-order
                      df_dx(i,j) = (1.5d0 * (f(i+1,j+1) - f(i,j+1))     &
                                  - 0.5d0 * (f(i+1,j+2) - f(i,j+2))) / dx
                   else  ! first-order
                      df_dx(i,j) = (f(i+1,j+1) - f(i,j+1)) / dx
                   endif
                else  ! decreasing j is upstream
                   diff1 = f(i+1,j) - f(i,j)
                   diff2 = f(i+1,j-1) - f(i,j-1)
                   if (abs(diff2) < 3.d0*abs(diff1)) then   ! second-order
                      df_dx(i,j) = (1.5d0 * (f(i+1,j) - f(i,j))         &
                                  - 0.5d0 * (f(i+1,j-1) - f(i,j-1))) / dx
                   else  ! first-order
                      df_dx(i,j) = (f(i+1,j) - f(i,j)) / dx
                   endif
                endif    ! sum1 > sum2

                ! Compute df_dy

                sum1 = f(i+1,j+1) + f(i+1,j)
                sum2 = f(i,j+1) + f(i,j)
             
                if (sum1 > sum2) then  ! increasing i is upstream
                   diff1 = f(i+1,j+1) - f(i+1,j)
                   diff2 = f(i+2,j+1) - f(i+2,j)
                   if (abs(diff2) < 3.d0*abs(diff1)) then   ! second-order
                      df_dy(i,j) = (1.5d0 * (f(i+1,j+1) - f(i+1,j))     &
                                  - 0.5d0 * (f(i+2,j+1) - f(i+2,j))) / dy
                   else  ! first-order
                      df_dy(i,j) = (f(i+1,j+1) - f(i+1,j)) / dy
                   endif
                else  ! decreasing i is upstream
                   diff1 = f(i,j+1) - f(i,j)
                   diff2 = f(i-1,j+1) - f(i-1,j)
                   if (abs(diff2) < 3.d0*abs(diff1)) then   ! second-order
                      df_dy(i,j) = (1.5d0 * (f(i,j+1) - f(i,j))         &
                                  - 0.5d0 * (f(i-1,j+1) - f(i-1,j))) / dy
                   else  ! first-order
                      df_dy(i,j) = (f(i,j+1) - f(i,j)) / dy
                   endif
                endif    ! sum1 > sum2

             endif    ! summask

          enddo     ! i
          enddo     ! j

       elseif (grad_flag2 == 2) then  ! second order, except where first order is closer to (centered) benchmark

          do j = 2, ny-2   ! loop does not include all of halo
          do i = 2, nx-2

             summask = imask(i,j) + imask(i+1,j) + imask(i,j+1) + imask(i+1,j+1)

             if (summask > 0) then

                ! Compute df_dx

                ! centered difference
                dfdx_c = (f(i+1,j+1) - f(i,j+1) + f(i+1,j) - f(i,j)) / (2.d0*dx)

                ! identify upstream direction

                sum1 = f(i+1,j+1) + f(i,j+1)
                sum2 = f(i+1,j) + f(i,j)

                if (sum1 > sum2) then  ! increasing j is upstream

                   ! first-order difference 
                   dfdx_1 = (f(i+1,j+1) - f(i,j+1)) / dx

                   ! second-order difference
                   dfdx_2 = (1.5d0 * (f(i+1,j+1) - f(i,j+1))     &
                                - 0.5d0 * (f(i+1,j+2) - f(i,j+2))) / dx

                else   ! decreasing j is upstream
                   
                   ! first-order difference
                   dfdx_1 = (f(i+1,j) - f(i,j)) / dx

                   ! second-order difference
                   dfdx_2 = (1.5d0 * (f(i+1,j) - f(i,j))         &
                                - 0.5d0 * (f(i+1,j-1) - f(i,j-1))) / dx

                endif   ! sum1 > sum2

                diff1 = abs(dfdx_1 - dfdx_c)
                diff2 = abs(dfdx_2 - dfdx_c)

                if (diff2 < 2.d0*diff1) then  ! use second-order difference
                   df_dx(i,j) = dfdx_2
                else  ! use first-order difference
                   df_dx(i,j) = dfdx_1
                endif

                ! Compute df_dy

                ! centered difference
                dfdy_c = (f(i+1,j+1) - f(i+1,j) + f(i,j+1) - f(i,j)) / (2.d0*dy)

                ! identify upstream direction

                sum1 = f(i+1,j+1) + f(i+1,j)
                sum2 = f(i,j+1) + f(i,j)

                if (sum1 > sum2) then  ! increasing i is upstream

                   ! first-order difference 
                   dfdy_1 = (f(i+1,j+1) - f(i+1,j)) / dy

                   ! second-order difference
                   dfdy_2 = (1.5d0 * (f(i+1,j+1) - f(i+1,j))     &
                           - 0.5d0 * (f(i+2,j+1) - f(i+2,j))) / dy

                else   ! decreasing j is upstream
                   
                   ! first-order difference
                   dfdy_1 = (f(i,j+1) - f(i,j)) / dy

                   ! second-order difference
                   dfdy_2 = (1.5d0 * (f(i,j+1) - f(i,j))         &
                           - 0.5d0 * (f(i-1,j+1) - f(i-1,j))) / dy

                endif   ! sum1 > sum2

                diff1 = abs(dfdy_1 - dfdy_c)
                diff2 = abs(dfdy_2 - dfdy_c)

                if (diff2 < 2.0*diff1) then  ! use second-order difference
                   df_dy(i,j) = dfdy_2
                else  ! use first-order difference
                   df_dy(i,j) = dfdy_1
                endif

             endif  ! summask

          enddo   ! i
          enddo   ! j

       endif   ! grad_flag2

       ! fill in halo values
       call staggered_parallel_halo(df_dx)
       call staggered_parallel_halo(df_dy)

    endif  ! accuracy_flag

!WHL - debug
    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'df_dx:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'df_dy:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo

    endif


  end subroutine glissade_upstream_gradient

!****************************************************************************

  end module glissade_grid_operators

!****************************************************************************
