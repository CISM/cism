! test_advect.f90
! Magnus Hagdorn, January 2005
!
! testing the advection code

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program testadvect
  use erosion_advect
  use netcdf
  implicit none

  character(len=30) filename
  integer timeslice
  real(kind=dp) x0,y0

  ! netCDF stuff
  integer status
  integer ncid,uxid,uyid,xid,yid,timeid,dimid
  real(kind=dp), dimension(:,:), allocatable :: ux, uy
  real(kind=dp), dimension(:), allocatable :: time,x,y
  real(kind=dp) :: deltax,deltay,deltat
  integer numx,numy,numt
  type(coordsystem_type) :: coords

  ! output
  integer numtimes,i,t
  real(kind=dp), parameter :: dt = 100.
  real(kind=dp), dimension(:), allocatable :: times, posx,posy

  ! set up
  write(*,*) 'Enter name of netCDF file containing velo field'
  read(*,*) filename
  write(*,*) 'Enter time slice to be read'
  read(*,*) timeslice
  write(*,*) 'Enter starting position'
  read(*,*) x0,y0
  
  ! open netCDF file and get variable sizes
  status = nf90_open(trim(filename),NF90_NOWRITE,ncid)
  call nc_errorhandle(status)
  status = nf90_inq_dimid(ncid, "time", dimid)
  call nc_errorhandle(status)
  status = nf90_inquire_dimension(ncid,dimid,len=numt)
  call nc_errorhandle(status)
  status = nf90_inq_dimid(ncid, "x0", dimid)
  call nc_errorhandle(status)
  status = nf90_inquire_dimension(ncid,dimid,len=numx)
  call nc_errorhandle(status)
  status = nf90_inq_dimid(ncid, "y0", dimid)
  call nc_errorhandle(status)
  status = nf90_inquire_dimension(ncid,dimid,len=numy)
  call nc_errorhandle(status)
  ! allocate arrays
  allocate(ux(numx,numy))
  allocate(uy(numx,numy))
  allocate(x(numx))
  allocate(y(numy))
  allocate(time(numt))
  ! get variable ids and load dimension variables
  status = nf90_inq_varid(ncid,"time",timeid)
  call nc_errorhandle(status)
  status = nf90_get_var(ncid,timeid,time)
  call nc_errorhandle(status)
  status = nf90_inq_varid(ncid,"x0",xid)
  call nc_errorhandle(status)
  status = nf90_get_var(ncid,xid,x)
  call nc_errorhandle(status)
  deltax = x(2)-x(1)
  status = nf90_inq_varid(ncid,"y0",yid)
  call nc_errorhandle(status)
  status = nf90_get_var(ncid,yid,y)
  call nc_errorhandle(status)
  deltay = y(2)-y(1)
  status = nf90_inq_varid(ncid,"uvel",uxid)
  call nc_errorhandle(status)
  status = nf90_inq_varid(ncid,"vvel",uyid)
  call nc_errorhandle(status)
 
  coords = coordsystem_new(x(1),y(1),deltax,deltay,numx,numy)

  call er_advect2d_init(coords)

  open(10,file='test_advect.data',status='unknown')
  write(10,*) time(timeslice), x0, y0
  do t=timeslice,numt-1
     deltat = time(t+1)-time(t)
     write(*,*) 'Processing year ',time(t)
     status = nf90_get_var(ncid,uxid,ux,(/1,1,1,t/),(/numx,numy,1,1/))
     call nc_errorhandle(status)
     status = nf90_get_var(ncid,uyid,uy,(/1,1,1,t/),(/numx,numy,1,1/))
     call nc_errorhandle(status)

     numtimes = int(deltat/dt)
     allocate(times(numtimes), posx(numtimes), posy(numtimes))
     do i=1,numtimes
        times(i) = i*dt
     end do
     
     call set_velos(ux,uy)
     
     call er_advect2d(times,x0,y0,posx,posy)
  
     do i=1,numtimes
        write(10,*) time(t)+times(i),posx(i),posy(i)
     end do
     
     x0 = posx(numtimes)
     y0 = posy(numtimes)
     deallocate(times,posx,posy)
  end do
  close(10)
  status = nf90_close(ncid)
end program testadvect

subroutine nc_errorhandle(status)
  use netcdf
  implicit none
  integer, intent(in) :: status
  if (status.ne.NF90_NOERR) then
     write(*,*) nf90_strerror(status)
     stop
  end if
end subroutine nc_errorhandle
