! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_lithot1d.f90 - part of the Glimmer-CISM ice model  + 
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

! module for 1D temperature calculations in the upper lithosphere

module glide_lithot1d

contains
  subroutine init_lithot1d(model)
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    ! allocate memory for 1D code
    allocate(model%lithot%rhs(model%lithot%nlayer))
    allocate(model%lithot%subd(model%lithot%nlayer))
    allocate(model%lithot%diag(model%lithot%nlayer))
    allocate(model%lithot%supd(model%lithot%nlayer))
    
    ! setup coefficient matrix
    model%lithot%subd(:) =    - model%lithot%zfactors(1,:)
    model%lithot%diag(:) = 1. + model%lithot%zfactors(2,:)
    model%lithot%supd(:) =    - model%lithot%zfactors(3,:)
    ! and the boundary conditions
    ! top face
    ! simply match air temperature where no ice and basal temperature where ice
    model%lithot%subd(1) = 0.
    model%lithot%diag(1) = 1.
    model%lithot%supd(1) = 0.
    ! bottom face
    ! keep constant
    model%lithot%subd(model%lithot%nlayer) = 0.
    model%lithot%diag(model%lithot%nlayer) = 1.
    model%lithot%supd(model%lithot%nlayer) = 0.
  end subroutine init_lithot1d

  subroutine calc_lithot1d(model)
    use glide_types
    use glimmer_utils
    use glide_mask
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    integer i,j,k

    ! loop over grid
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          ! calculate RHS for upper BC
          if (is_ground(model%geometry%thkmask(i,j)) .and. .not. is_thin(model%geometry%thkmask(i,j)) ) then
             model%lithot%rhs(1) = model%temper%temp(model%general%upn,i,j) ! ice basal temperature
             model%lithot%mask(i,j) = .true.
          else
             if (model%lithot%mask(i,j)) then
                if (is_ocean(model%geometry%thkmask(i,j))) then
                   model%lithot%rhs(1) = model%lithot%mart
                else if (is_land(model%geometry%thkmask(i,j))) then
                   model%lithot%rhs(1) = model%climate%artm(i,j) ! air temperature outside ice sheet
                end if
             end if
          end if

          if (model%lithot%mask(i,j)) then
             ! calculate RHS for rest
             do k=2,model%lithot%nlayer-1
                model%lithot%rhs(k) = - model%lithot%subd(k)*model%lithot%temp(i,j,k-1) &
                     + (2.-model%lithot%diag(k))*model%lithot%temp(i,j,k) &
                     - model%lithot%supd(k)*model%lithot%temp(i,j,k+1)
             end do
             model%lithot%rhs(model%lithot%nlayer) = model%lithot%temp(i,j,model%lithot%nlayer)

             ! solve tri-diagonal matrix eqn
             call tridiag(model%lithot%subd(1:), &
                  model%lithot%diag(:), &
                  model%lithot%supd(:model%lithot%nlayer), &
                  model%lithot%temp(i,j,:) ,                 &
                  model%lithot%rhs(:))
          end if
       end do
    end do
  end subroutine calc_lithot1d

  subroutine finalise_lithot1d(model)
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    deallocate(model%lithot%rhs)
    deallocate(model%lithot%subd)
    deallocate(model%lithot%diag)
    deallocate(model%lithot%supd)
  end subroutine finalise_lithot1d

end module glide_lithot1d
