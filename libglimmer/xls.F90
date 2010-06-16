! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  xls.f90 - part of the Glimmer-CISM ice model             + 
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

!Debugging module: contains routines to write out data fields quickly
!that can be read into Matlab
module xls
contains
  subroutine write_xls(fname, data)
    character(len=*) :: fname
    double precision data(:,:)
    integer :: i , j
    open(11, file=fname)

    do i = 1, size(data, 2)
       write(11, *)(data(j,i),j=1,size(data,1))
    end do

    close(11)

  end subroutine write_xls

  subroutine write_xls_int(fname, data)
    character(len=*) :: fname
    integer data(:,:)
    integer :: i , j
    open(11, file=fname)

    do i = 1, size(data, 2)
       write(11, *)(data(j,i),j=1,size(data,1))
    end do

    close(11)

  end subroutine write_xls_int



  subroutine write_xls_3d(fname, data)
    character(len=*) :: fname
    double precision data(:,:,:)
    integer :: i , j, k
    open(11, file=fname)
    do k = 1, size(data, 1)
       do i = 1, size(data, 3)
          write(11, *)(data(k,j,i),j=1,size(data,2))
       end do
    end do
    close(11)
  end subroutine write_xls_3d

  subroutine write_sparse_system(name, ia, ja, a,nelt)
    double precision, dimension(:) :: a
    integer, dimension(:) :: ia, ja
    integer :: nelt, i
    character(len=*) :: name

    open(11, file=name)
    write(11,*)(ia(i),i=1,nelt)
    write(11,*)(ja(i),i=1,nelt)
    write(11,*)(a(i),i=1,nelt)
    close(11)
  end subroutine write_sparse_system
end module xls
