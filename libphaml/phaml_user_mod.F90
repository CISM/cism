! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  phaml_user_mod.F90 - part of the Glimmer-CISM ice model      + 
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



module phaml_user_mod
    !----------------------------------------------------
    ! This module contains user global data, and functions for formatting
    !
    !----------------------------------------------------
    use global
    use phaml
    implicit none


    integer, save :: gnsn, gewn, gdns, gdew, num_arrays, modnum
    real(my_real), save, allocatable, dimension(:,:) :: uphaml
contains        
    subroutine user_init(ewn,nsn,dew,dns,modn)
        implicit none
        integer, intent(in) :: nsn,ewn,modn
        real(my_real), intent(in) :: dns,dew
        gnsn = nsn  !ny
        gewn = ewn  !nx
        gdns = dns  !dy
        gdew = dew  !dx
        modnum = modn
        !this defines how many arrays are passed in update_usermod
        num_arrays = 1
        !allocate(uphaml(gewn, gnsn))
        !allocate(uphaml_one(gewn*gnsn))
    end subroutine user_init
    
    subroutine user_close()
        if (size(uphaml) > 0) then
            deallocate(uphaml)
        end if
    end subroutine user_close
    
    subroutine array_init()
        !if (size(uphaml) .le. 0) then
            allocate(uphaml(gewn, gnsn))
        !end if
    end subroutine array_init
    
    subroutine array_close()
        if (size(uphaml) > 0) then
            deallocate(uphaml)
            !deallocate(uphaml_one)
        end if
    end subroutine array_close
    !-------------------------------------------------------------------------
    subroutine combine_arrays(array1, array2, retarray)
        implicit none
        real(my_real), intent(in),dimension(:) :: array1,array2
        real(my_real), intent(out),dimension(:) :: retarray
        integer :: i,maxind
        
        maxind = size(array1) + size(array2)
        do i=1, maxind
            retarray(i) = array1(i)
            retarray(maxind + i) = array2(i)
        end do
    end subroutine combine_arrays
    !-------------------------------------------------------------------------        
    subroutine split_arrays(array1, array2, retarray, maxind)
        implicit none
        real(my_real), intent(out),dimension(:) :: array1,array2
        real(my_real), intent(in),dimension(:) :: retarray
        integer, intent(in) :: maxind
        integer :: i
        
        do i=1, maxind
            array1(i) = retarray(i)
            array2(i) = retarray(maxind + i)
        end do
    end subroutine split_arrays
    !-------------------------------------------------------------------------    
    subroutine reshape_array_to_one(twod,oned)
        implicit none
        real(my_real), intent(in),dimension(:,:) :: twod
        real(my_real), intent(out),dimension(:) :: oned
        integer :: r,c
        
         do r=0,gnsn-1
            do c=0,gewn-1
               oned(gewn*r+c+1) = twod(c+1,r+1)
            end do
        end do
    end subroutine reshape_array_to_one
    !-------------------------------------------------------------------------    
    subroutine reshape_array_to_two(twod,oned)
        implicit none
        real(my_real), intent(out),dimension(:,:) :: twod
        real(my_real), intent(in),dimension(:) :: oned
        integer :: count,r,c
        
        do count=0,size(oned) - 1
            r = count/gewn + 1
            c = int(Mod(count,gewn)) + 1
            twod(c,r) = real(oned(count+1))        
        end do
    end subroutine reshape_array_to_two
    
    !------------------------------------------------------------------------------
    !SUBROUTINE: get_xyarrays
    !ARGUMENTS: model (glimmer), x, y 
    !DESCRIPTION:
    !------------------------------------------------------------------------------
    subroutine get_xyarrays(x,y)
        implicit none
        real(my_real), intent(inout) :: x(:),y(:)
        integer :: count,r,c

        count = 1
        !populate the arrays with the xy coordinates
        do r=0,gnsn-1
            do c=0,gewn-1
               x(count) = REAL(c*gdew)
               y(count) = REAL(r*gdns)
               count = count + 1
            end do
        end do
    end subroutine get_xyarrays
    
    
    !these two functions get you the array indices to 
    !access the 2d grid given the x,y quadrature coordinates in the domain
    !not safe to use for every x,y
    integer function getew(x)
        implicit none
        real(my_real), intent(in) :: x
        getew = int((x/gdew)) + 1
    end function getew
    
    integer function getns(y)
        implicit none
        real(my_real), intent(in) :: y
        getns = int((y/gdns)) + 1
    end function getns
    
end module phaml_user_mod
