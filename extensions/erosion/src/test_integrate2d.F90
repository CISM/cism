! test_integrate2d.f90
! Magnus Hagdorn, March 2003
!
! testing the 2d integration stuff

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program testintegrate2d
  use erosion_integrate2d
  use glimmer_sparse
  use geometry
  use glimmer_coordinates
  implicit none
  type(geom_poly) :: patch
  real(kind=dp), parameter :: xsize = 200.
  real(kind=dp), parameter :: ysize = 200.
  real(kind=dp) :: deltaxy
  integer :: numi,numj
  integer i
  type(sparse_matrix_type) :: weight
  real(kind=dp), dimension(:,:), allocatable :: vec
  real(kind=dp), dimension(:), allocatable ::res
  real(kind=dp) :: area
  type(coordsystem_type) :: coords

  deltaxy = 10.

  numi = int(xsize/deltaxy)+1
  numj = int(ysize/deltaxy)+1

  coords = coordsystem_new(0.d0,0.d0,deltaxy,deltaxy,numi,numj)

  allocate(res(numi*numj))
  allocate(vec(numi,numj))

  ! initialise grid
  call init_integrate2d

  ! test shape
  patch = poly_new(4)
  call poly_add_vert(patch,60.,15.)
  call poly_add_vert(patch,95.,46.)
  call poly_add_vert(patch,55.,80.)
  call poly_add_vert(patch,10.,42.)

  write(*,*) 'Size :', 0.5*real(poly_area2(patch))


  call new_sparse_matrix(1000,weight)
  call calc_weight(coords,weight, patch, 1)

  vec = 0
  do i=1,weight%n
     vec(mod(weight%row(i),numi)+1,weight%row(i)/numi+1) = weight%val(i)
  end do

  do i=0,8
     deltaxy = 2.**i
     numi = int(xsize/deltaxy) + 1
     numj = int(ysize/deltaxy) + 1
     coords = coordsystem_new(0.d0,0.d0,deltaxy,deltaxy,numi,numj)
     area = 0.
     weight%n = 0
     call calc_weight(coords,weight, patch, 1)
     area = sum(weight%val(1:weight%n))*deltaxy*deltaxy
     write(*,*) deltaxy,weight%n,area
  end do

end program testintegrate2d
