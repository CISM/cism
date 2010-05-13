
module ppm
!
! V. Lee 2008
!
! ppmthck cannibalises the Lagrangian version of the piecewise parabolic 
! method (ppm) called VIRGINIA HYDRODYNAMICS #1 or VH-1. 
!
!  The Virginia Numerical Bull Session ideal hydrodynamics PPMLR
!  http://wonka.physics.ncsu.edu/pub/VH-1/index.html
! 
! VH-1 is based on 
! Collela & Woodward (1984). The Piecewise Parabolic Method (PPM) for 
! Gas-Dynamical Simulations. J. Comp. Phys. 54, p 174-201.
! 
! Adapted for GLIMMER/CISM by W. Lipscomb, Mar. 2009

  use glimmer_paramets, only : dp, sp

  implicit none
  private
  public :: ppminitial, ppmthck, ppmfinal
    
!   
! PPM zone variables
! 
!  imax  - number of physical zones in the x direction
!  jmax  - number of physical zones in the y direction
!  kmax  - number of physical zones in the z direction
!
  integer :: imax, jmax, kmax

!  zro   - density: zone average
!  zux   - u-velocity: zone edge
!  zux   - v-velocity: zone edge
!
  real (kind = dp), allocatable, dimension(:,:,:) :: zro, zux, zuy

!  zxa   - x grid zone face locations
!  zya   - y grid zone face locations
!  zza   - z grid zone face locations
!  zdx   - xa(i+1) - xa(i)
!  zdy   - ya(i+1) - ya(i)
!  zdz   - za(i+1) - za(i)
!
  real (kind = dp), allocatable, dimension(:) :: zxa, zxc, zdx
  real (kind = dp), allocatable, dimension(:) :: zya, zyc, zdy
  real (kind = dp), allocatable, dimension(:) :: zza, zzc, zdz
  
! PPM global variables
    integer :: ndim                                              ! Number of dimensions
    integer :: ngeomx, ngeomy, ngeomz                            ! Type of coordinates
    integer :: nleftx, nlefty, nleftz, nrightx, nrighty, nrightz ! Boundary condition type
    real (kind = dp) :: thcklx, thckrx, flxinx, flxotx  ! BC: Thickness and flux at left and right hand boundaries
    real (kind = dp) :: thckly, thckry, flxiny, flxoty
    real (kind = dp), dimension(6) :: dinflox, dotflox           ! BC: Density in ghost cells at boundaries
    real (kind = dp), dimension(6) :: dinfloy, dotfloy
    real (kind = dp), dimension(3) :: uinflo, uotflo             ! BC: Velocity in ghost cells
    real (kind = dp), dimension(3) :: vinflo, votflo
  
contains

!----------------------------------------------------------------------

  subroutine ppminitial(ewn,     nsn,  &
                        dew,     dns  )

    implicit none

    integer, intent(in) :: ewn, nsn             ! domain size
    real (kind = dp), intent(in) :: dew, dns    ! grid cell dimensions
    real (kind = dp) :: xmin, xmax, ymin, ymax, zmin, zmax       ! Size of domain 

   call initialisedomain(ewn,    nsn,    &
                         dew,    dns,    &
                         imax,   jmax,   &
                         kmax,           &
                         xmin,   xmax,   &
                         ymin,   ymax,   &
                         zmin,   zmax,   &
                         ndim,   ngeomx, &
                         ngeomy, ngeomz)

   call initialbc( nleftx,nrightx,nlefty,nrighty,nleftz,nrightz, &
       & dinflox,dotflox,dinfloy,dotfloy,uinflo,uotflo,vinflo,votflo, &
       & thcklx,thckrx,thckly,thckry,flxinx,flxotx,flxiny,flxoty )

! Allocate and initialise ppm arrays
    allocate(zro(imax,jmax,kmax)); zro=0.0_dp
    allocate(zux(imax+1,jmax,kmax),zuy(imax,jmax+1,kmax)); zux=0.0_dp; zuy=0.0_dp
    allocate(zxa(imax),zxc(imax),zdx(imax)); zxa=0.0_dp; zxc=0.0_dp; zdx=0.0_dp
    allocate(zya(jmax),zyc(jmax),zdy(jmax)); zya=0.0_dp; zyc=0.0_dp; zdy=0.0_dp
    allocate(zza(kmax),zzc(kmax),zdz(kmax)); zza=0.0_dp; zzc=0.0_dp; zdz=0.0_dp

! PPM zones
!
!     | = zone edge 
!     * = ice thickness location
!     o = uvel location
!     + = extra velocity for ppm (usually set to zero since outside the domain)
!     number of cells = ewn = imax
!
!                x=xmin                              x=xmax       
!         zxa(1)                             zxa(imax)
!                zdx(1)                            zdx(imax)
!                zxc(1)                            zxc(imax)           
!            |           |           |           |           |
!            +     *     o     *     o     *     o     *     +     
!            |           |           |           |           |
!
!
    
    call grid(imax,xmin,xmax,zxa,zxc,zdx)
    call grid(jmax,ymin,ymax,zya,zyc,zdy)
    call grid(kmax,zmin,zmax,zza,zzc,zdz)
  
  end subroutine ppminitial

!----------------------------------------------------------------------
  
  subroutine ppmthck( ewn,         nsn,         &
                      dt,                       &
                      uflx,        vflx,        &
                      thck,        stagthck,    &
                      usrf,        lsrf,        &
                      acab,        leuler)

    implicit none

    integer, intent(in) :: ewn, nsn             ! domain size
    real (kind = dp), intent(in) :: dt          ! time step
    real (kind = dp), intent(in), dimension(:,:) :: uflx, vflx, stagthck,lsrf
    real (kind = sp), intent(in), dimension(:,:) :: acab
    real (kind = dp), intent(inout), dimension(:,:) :: thck
    real (kind = dp), intent(out), dimension(:,:) :: usrf
    logical, intent(in) :: leuler

!
! Local variables
!   
    integer :: ew, ns
    integer :: ii, jj, kk                                        ! Zone indices
    real (kind = dp) :: deltat                                   ! Operator splitting timestep
    real (kind = dp), allocatable, dimension(:,:) :: uvel, vvel  ! Vertically integrated velocity components

! Recover vertically integrated velocity field
    allocate(uvel(ewn-1,nsn-1),vvel(ewn-1,nsn-1))
    uvel=0.0_dp
    vvel=0.0_dp
    do ns=1,nsn-1
       do ew=1,ewn-1
          if (abs(stagthck(ew,ns)) > 1.0e-3) then
             uvel(ew,ns)=uflx(ew,ns)/stagthck(ew,ns)
             vvel(ew,ns)=vflx(ew,ns)/stagthck(ew,ns)
          end if
       end do
    end do

    call assignppmvars (ewn,     nsn,         &
                        imax,    jmax,        &
                        kmax,    zro,         &
                        zux,     zuy,         &
                        thck,                 &
                        uvel,    vvel,        &
                        nleftx,  nrightx,     &
                        nlefty,  nrighty)

    deallocate(uvel,vvel)

    if (ndim > 1) then
       deltat=dt/2.0_dp
    else
       deltat=dt
    end if

! Calculate the new 'density' zro
    if (imax > 1) call sweepx(deltat,imax,jmax,kmax,zro,zux,zxa,zdx,ngeomx, &
         & nleftx,nrightx,dinflox,dotflox,uinflo,uotflo,thcklx,thckrx, &
         & flxinx,flxotx,leuler)
    if (jmax > 1) call sweepy(deltat,imax,jmax,kmax,zro,zuy,zxa,zdx,zya,zdy, &
         & ngeomy,nlefty,nrighty,dinfloy,dotfloy,vinflo,votflo,thckly,thckry, &
         & flxiny,flxoty,leuler)
! Operator splitting: reverse sweep direction for approximate 2nd order
    if (ndim > 1) then
       call sweepy(deltat,imax,jmax,kmax,zro,zuy,zxa,zdx,zya,zdy,ngeomy, &
            & nlefty,nrighty,dinfloy,dotfloy,vinflo,votflo,thckly,thckry, &
            & flxiny,flxoty,leuler)
       call sweepx(deltat,imax,jmax,kmax,zro,zux,zxa,zdx,ngeomx,nleftx,nrightx, &
            & dinflox,dotflox,uinflo,uotflo,thcklx,thckrx,flxinx,flxotx,leuler)
    end if

! Add accumulation
    kk=kmax
    do jj=1,jmax
       ns=jj
       do ii=1,imax
          ew=ii
          zro(ii,jj,kk) = zro(ii,jj,kk)+real(acab(ew,ns),dp)*dt
       end do
    end do
    where(zro < 0.0_dp) zro=0.0_dp

    if (nleftx==4 .and. imax > 1) zro(1,:,:)=thcklx
    if (nrightx==4 .and. imax > 1) zro(imax,:,:)=thckrx
    if (nlefty==4 .and. jmax > 1) zro(:,1,:)=thckly
    if (nrighty==4 .and. jmax > 1) zro(:,jmax,:)=thckry

! Return to thck
    kk=kmax
    do ns=1,nsn
       jj=ns
       do ew=1,ewn
          ii=ew
          thck(ew,ns) = zro(ii,jj,kk)
       end do
    end do

    usrf=thck+lsrf

  end subroutine ppmthck

!----------------------------------------------------------------------

  subroutine ppmfinal
    implicit none

    deallocate(zux, zuy)
    deallocate(zxa,zxc,zdx,zya,zyc,zdy,zza,zzc,zdz)
    deallocate(zro)

  end subroutine ppmfinal

!----------------------------------------------------------------------

  subroutine initialisedomain(ewn,    nsn,    &
                              dew,    dns,    &
                              imax,   jmax,   &
                              kmax,           &
                              xmin,   xmax,   & 
                              ymin,   ymax,   &
                              zmin,   zmax,   &
                              ndim,   ngeomx, &
                              ngeomy, ngeomz)

    implicit none

    integer, intent(in)  :: ewn, nsn
    real (kind = dp), intent(in)  :: dew, dns
    integer, intent(out) :: imax,jmax,kmax, ndim, ngeomx,ngeomy,ngeomz
    real (kind = dp), intent(out) :: xmin,xmax,ymin,ymax,zmin,zmax

! Number of zones
    imax=ewn
    jmax=nsn
    kmax=1

! Size of domain
    xmin=0.0_dp; xmax=dew*real(ewn-1,dp)
    if (imax==1) xmax=1.0_dp
    ymin=0.0_dp; ymax=dns*real(nsn-1,dp)
    if (jmax==1) ymax=1.0_dp
    zmin=0.0_dp; zmax=1.0_dp

! Number of dimensions
    if (imax > 1 .and. jmax > 1) then
       ndim = 2
    else
       ndim = 1
    end if

! Geometry
! ngeom  
!   = 0  :  planar Cartesian
!   = 1  :  cylindrical radial 
!   = 2  :  spherical   radial 
!   = 3  :  cylindrical angle  
!   = 4  :  spherical polar angle (theta)
!   = 5  :  spherical azimu angle (phi)  
    ngeomx = 0
    ngeomy = 0
    ngeomz = 0

  end subroutine initialisedomain

!----------------------------------------------------------------------

  subroutine initialbc( nleftx,nrightx,nlefty,nrighty,nleftz,nrightz, &
       & dinflox,dotflox,dinfloy,dotfloy,uinflo,uotflo,vinflo,votflo, &
       & thcklx,thckrx,thckly,thckry,flxinx,flxotx,flxiny,flxoty )

    integer, intent(out) :: nleftx, nlefty, nleftz, nrightx, nrighty, nrightz ! Boundary condition type
    real (kind = dp), intent(out) :: thcklx, thckrx, flxinx, flxotx   ! BC: Thickness and flux at left and right hand boundaries
    real (kind = dp), intent(out) :: thckly, thckry, flxiny, flxoty
    real (kind = dp), dimension(6), intent(out) :: dinflox, dotflox           ! BC: Density in ghost cells at boundaries
    real (kind = dp), dimension(6), intent(out) :: dinfloy, dotfloy
    real (kind = dp), dimension(3), intent(out) :: uinflo, uotflo             ! BC: Velocity in ghost cells
    real (kind = dp), dimension(3), intent(out) :: vinflo, votflo

! Boundary condition type
! nleft, nright
!   = 0  :  reflecting boundary condition
!   = 1  :  inflow/outflow boundary condition (zero gradient)
!   = 2  :  fixed inflow/outflow boundary condition (set dinflo,dotflo,uinflo,uotflo)
!   = 3  :  periodic
!   = 4  :  fixed thickness (set thckl,thckr,dinflo,dotflo,uinflo,uotflo)
!   = 5  :  fixed flux, not complete (set flxin, flxot, dinflo,dotflo,uinflo,uotflo)
    nleftx = 3; nrightx= 3
    nlefty = 4; nrighty= 4
    nleftz = 0; nrightz= 0

! Assign density in ghost cells for fixed inflow and outflow 
    dinflox(1:6) = 0.0_dp; dotflox(1:6) = 0.0_dp
    dinfloy(1:6) = 0.0_dp; dotfloy(1:6) = 0.0_dp
! Assign velocity in ghost cells
    uinflo(1:3) = 0.0_dp; uotflo(1:3) = 0.0_dp
    vinflo(1:3) = 0.0_dp; votflo(1:3) = 0.0_dp

! Set thickness at boundaries
    thcklx=0.0_dp; thckrx=0.0_dp
    thckly=0.0_dp; thckry=0.0_dp
! Set flux at boundaries
    flxinx=0.0_dp; flxotx=0.0_dp
    flxiny=0.0_dp; flxoty=0.0_dp

  end subroutine initialbc

!----------------------------------------------------------------------

  subroutine grid( nzones, xmin, xmax, xa, xc, dx)
!
! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! VH1 xa(1) is left boundary location - at xmin-dx
!     xa(nzones+1) is right boundary location - at xmax+dx
! Ice sheet xc(1) is left boundary location - at xmin
!     xc(nzones) is right boundary location - at xmax
!
    implicit none

    integer, intent(in) :: nzones
    real (kind = dp), intent(in) :: xmin, xmax
    real (kind = dp), dimension(nzones), intent(out) :: xa, dx, xc
!
! LOCALS
!
    real (kind = dp) :: dxfac
    integer :: n

    if (nzones==1) then
       dxfac = (xmax - xmin) / real(nzones,dp)
       do n = 1, nzones
          xa(n) = xmin + real((n-1),dp)*dxfac
          dx(n) = dxfac
          xc(n) = xa(n) + 0.5_dp*dx(n)
       enddo
    else
       dxfac = (xmax - xmin) / real(nzones-1,dp)
       do n = 1, nzones
          xc(n) = xmin + real(n-1,dp)*dxfac
          dx(n) = dxfac
          xa(n) = xc(n) - 0.5_dp*dx(n)
       enddo
    end if

  end subroutine grid

!----------------------------------------------------------------------

  subroutine assignppmvars(ewn,     nsn,         &
                           imax,    jmax,        &
                           kmax,    zro,         &
                           zux,     zuy,         &
                           thck,                 &
                           uvel,    vvel,        &
                           nleftx,  nrightx,     &
                           nlefty,  nrighty)

    implicit none

    integer, intent(in) :: ewn, nsn
    integer, intent(in) :: imax,jmax,kmax,nleftx,nrightx,nlefty,nrighty
    real (kind = dp), intent(in), dimension(:,:) :: thck, uvel, vvel
    real (kind = dp), intent(inout), dimension(:,:,:) :: zro, zux, zuy

    integer :: ii, jj, kk, ew, ns

! PPM variable arrangement
! Density is zone centred (z*c), velocity is at zone edges (z*a)
!
!     | = zone edge 
!     * = ice thickness location
!     o = uvel location
!     + = extra velocity for ppm (usually set to zero since outside the domain)
!     number of cells = ewn = imax
!
!                x=xmin                              x=xmax       
!         zxa(1)                             zxa(imax)
!                zdx(1)                            zdx(imax)
!                zxc(1)                            zxc(imax)           
!            |           |           |           |           |
!            +     *     o     *     o     *     o     *     +     
!            |           |           |           |           |
!                zro(1)                            zro(imax)
!                thck(1)                           thck(ewn)
!          zux(1)                                         zux(imax+1) 
!                     uvel(1)                uvel(ewn-1)
!
!

    do kk=1,kmax

! Assign ice thickness to zone density
       do jj=1,jmax
          ns=jj
          do ii=1,imax
             ew=ii
             zro(ii,jj,kk)=thck(ew,ns)
          end do
       end do

! Assign u-velocity to zone edge u-velocity, 
!    i.e. ( uvel_{i+1/2,j+1/2} + uvel_{i+1/2,j-1/2} )/2
       if (jmax > 1) then
          do jj=2,jmax-1
             ns=jj
             do ii=1,imax-1
                ew=ii
                zux(ii+1,jj,kk)=uvel(ew,ns)+uvel(ew,ns-1)
             end do
          end do
          zux=0.5_dp*zux
       else
          do jj=1,jmax
             ns=jj
             do ii=1,imax-1
                ew=ii
                zux(ii+1,jj,kk)=uvel(ew,ns)
             end do
          end do
       end if

! Assign v-velocity to zone edge v-velocity,
!    i.e. ( vvel_{i+1/2,j+1/2} + vvel_{i-1/2,j+1/2} )/2    
       if (imax > 1) then
          do jj=1,jmax-1
             ns=jj
             do ii=2,imax-1
                ew=ii
                zuy(ii,jj+1,kk)=vvel(ew,ns)+vvel(ew-1,ns)
             end do
          end do
          zuy=0.5_dp*zuy
       else
          do jj=1,jmax-1
             ns=jj
             do ii=1,imax
                ew=ii
                zuy(ii,jj+1,kk)=vvel(ew,ns)
             end do
          end do
       end if

    end do

! Take advantage of periodic boundary conditions to calculate average velocities
    if (imax > 1 .and. jmax > 1) then
       if (nleftx == 3) then
          do kk=1,kmax
             do jj=1,jmax-1
                ns=jj
                zuy(1,jj+1,kk) = 0.5_dp*(vvel(1,ns)+vvel(ewn-1,ns))
             end do
          end do
       end if
       if (nrightx == 3) then
          do kk=1,kmax
             do jj=1,jmax-1
                ns=jj
                zuy(imax,jj+1,kk) = 0.5_dp*(vvel(1,jj)+vvel(ewn-1,jj))
             end do
          end do
       end if
       if (nlefty == 3) then
          do kk=1,kmax
             do ii=1,imax-1
                ew=ii
                zux(ii+1,1,kk) = 0.5_dp*(uvel(ew,1)+uvel(ew,nsn-1))
             end do
          end do
       end if
       if (nrighty == 3) then
          do kk=1,kmax
             do ii=1,imax-1
                ew=ii
                zux(ii+1,jmax,kk) = 0.5_dp*(uvel(ew,1)+uvel(ew,nsn-1))
             end do
          end do
       end if
    end if

! Periodic boundary conditions: assign extra velocity points 
    if (nleftx == 3 .or. nrightx == 3) then
       do kk=1,kmax
          do jj=1,jmax-1
             if (nleftx == 3) zux(1,jj,kk) = zux(imax,jj,kk)
             if (nrightx == 3) zux(imax+1,jj,kk) = zux(2,jj,kk)
          end do
       end do
    end if
    if (nlefty == 3 .or. nrighty == 3) then
       do kk=1,kmax
          do ii=1,imax-1
             if (nlefty == 3) zuy(ii,1,kk) = zuy(ii,jmax,kk)
             if (nrighty == 3) zuy(ii,jmax+1,kk) = zuy(ii,2,kk)
          end do
       end do
    end if

  end subroutine assignppmvars

!----------------------------------------------------------------------

  subroutine sweepx(dt,imax,jmax,kmax,zro,zux,zxa,zdx,ngeomx,nleftx,nrightx, &
       & dinflox,dotflox,uinflo,uotflo,thcklx,thckrx,flxinx,flxotx,leuler)
!
! This subroutine performs 1D hydro sweeps in the X direction,
! looping over j by k rows.
!
    implicit none

! GLOBALS
    integer, intent(in) :: imax,jmax,kmax, ngeomx, nleftx,nrightx
    real (kind = dp), intent(in) :: dt, thcklx, thckrx, flxinx, flxotx
    real (kind = dp), intent(in), dimension(:) :: dinflox,dotflox, uinflo,uotflo
    real (kind = dp), intent(in), dimension(:) :: zxa, zdx
    real (kind = dp), intent(in), dimension(:,:,:) :: zux
    real (kind = dp), intent(inout), dimension(:,:,:) :: zro
    logical, intent(in) :: leuler

! Sweep variables: data structures used in 1D sweeps, dimensioned maxsweep
!
! parameters:
!  maxsweep - maximum size of any sweep
! variables:
!  nmin  - number of first real zone (=3)
!  nmax  - number of last  real zone (=i/jmax+2)
!  r     - density: zone average
!  xa    - grid zone face locations
!  dx    - xa(j+1) - xa(j)
!  dvol  - zone volume
!  xa0   - grid zone face locations, prior to lagrangian update
!  dx0   - dx, prior to lagrangian update
!  dvol0 - dvol, prior to lagrangian update
!  radius- radius out to current column of angular zones
!
    integer :: maxsweep          ! Dimension of 1D sweeps.  maxsweep must be as long
                                 ! as the longest of the 3D arrays PLUS the ghost 
                                 ! zones:  maxsweep = max(imax,jmax,kmax) + 12
    integer :: nmin, nmax
    real (kind = dp), allocatable, dimension(:) :: r, umid, dvol, dvol0
    real (kind = dp), allocatable, dimension(:) :: xa, xa0, dx, dx0
    real (kind = dp) :: radius

! LOCALS
    integer i, j, k, n


! Add six ghost zones to 1D sweep arrays

    nmin = 7
    nmax = imax + 6
    maxsweep = imax + 12

    allocate(r(maxsweep), umid(maxsweep))
    allocate(xa(maxsweep), xa0(maxsweep), dx(maxsweep), dx0(maxsweep))
    allocate(dvol(maxsweep), dvol0(maxsweep))

! Now Loop over each row...

    do k = 1, kmax
       do j = 1, jmax

! Initialize sweep variables
          umid=0.0_dp
          r=0.0_dp
          xa=0.0_dp
          xa0=0.0_dp
          dx=0.0_dp
          dx0=0.0_dp
          radius=1.0_dp

! Put state variables into 1D arrays, padding with 6 ghost zones
!
!     | = grid location 
!     * = r location
!     o = umid location
!  Number of physical cells = imax
!  Number of numerical cells = imax + 12 ghost cells (6 at each end)        
!  nmin/nmax = location of first/last physical cell
!
! Ghosts    ||               Physical cells                  ||     Ghost
! cells     ||                                               ||     cells
!   6           r(nmin)                             r(nmax)           6
!            |           |           |           |           |
!            o     *     o     *     o     *     o     *     o       
!            |           |           |           |           |
!   3    umid(nmin)                                     umid(nmax+1)  3
!
!

          do i = 1,imax
             n = i + nmin-1
             r  (n) = zro(i,j,k)
             xa0(n) = zxa(i)
             dx0(n) = zdx(i)
             xa (n) = zxa(i)
             dx (n) = zdx(i)
          end do
          do i=1,imax+1
             n = i + nmin-1
             umid(n)=zux(i,j,k)
          end do

! Set boundary conditions, compute volume elements, evolve flow, then remap
          call sweepbc( nleftx,nrightx,nmin,nmax,r,umid,xa,xa0,dx,dx0, &
               & dinflox,dotflox,uinflo,uotflo,thcklx,thckrx )
          call volume ( ngeomx,nmin,nmax,dvol,dvol0,xa,xa0,dx,dx0,radius )
          if (leuler) then
             call advance( dt,ngeomx,nleftx,nrightx,maxsweep,nmin,nmax, &
                  & r,umid,dvol0,xa0,dx0,flxinx,flxotx )
          else
             call evolve( dt,ngeomx,maxsweep,nmin,nmax,r,umid,dvol,xa,dx,radius )
             call remap ( dt,ngeomx,nleftx,nrightx,maxsweep,nmin,nmax,r, &
                  & dvol,dvol0,xa,xa0,dx,dx0,radius,flxinx,flxotx )
          end if

! Put updated values back into 3D arrays, dropping ghost zones

          do i = 1, imax
             n = i + 6
             zro(i,j,k) = r(n)
          enddo

       end do
    end do

    deallocate(r, umid, xa, xa0, dx, dx0, dvol, dvol0)

    return
  end subroutine sweepx

!----------------------------------------------------------------------

  subroutine sweepy(dt,imax,jmax,kmax,zro,zuy,zxa,zdx,zya,zdy,ngeomy,nlefty,nrighty, &
       & dinfloy,dotfloy,vinflo,votflo,thckly,thckry,flxiny,flxoty,leuler)
!
! This subroutine performs sweeps in the Y direction,
! looping over i by k rows.   
!
    implicit none
!
! GLOBALS
    integer, intent(in) :: imax,jmax,kmax, ngeomy, nlefty,nrighty
    real (kind = dp), intent(in) :: dt, thckly, thckry, flxiny, flxoty
    real (kind = dp), intent(in), dimension(:) :: dinfloy,dotfloy, vinflo,votflo
    real (kind = dp), intent(in), dimension(:) :: zxa, zdx, zya, zdy
    real (kind = dp), intent(in), dimension(:,:,:) :: zuy
    real (kind = dp), intent(inout), dimension(:,:,:) :: zro
    logical, intent(in) :: leuler

    integer :: maxsweep          ! Dimension of 1D sweeps.  maxsweep must be as long
                                 ! as the longest of the 3D arrays PLUS the ghost 
                                 ! zones:  maxsweep = max(imax,jmax,kmax) + 12
    integer :: nmin, nmax
    real (kind = dp), allocatable, dimension(:) :: r, umid, dvol, dvol0
    real (kind = dp), allocatable, dimension(:) :: xa, xa0, dx, dx0
    real (kind = dp) :: radius

! LOCALS
    integer i, j, k, n


! Add six ghost zones on each side of 1D arrays

    nmin = 7
    nmax = jmax + 6
    maxsweep = jmax + 12

    allocate(r(maxsweep), umid(maxsweep))
    allocate(xa(maxsweep), xa0(maxsweep), dx(maxsweep), dx0(maxsweep))
    allocate(dvol(maxsweep), dvol0(maxsweep))

! Now Loop over each column...

    do k = 1, kmax
       do i = 1, imax

! Initialize sweep variables
          umid=0.0_dp
          r=0.0_dp
          xa=0.0_dp
          xa0=0.0_dp
          dx=0.0_dp
          dx0=0.0_dp

          radius = zxa(i) + 0.5_dp*zdx(i)

! Put state variables into 1D arrays, padding with 6 ghost zones

          do j = 1, jmax
             n = j + nmin-1
             r  (n) = zro(i,j,k)
             xa (n) = zya(j)
             dx (n) = zdy(j)
             xa0(n) = zya(j)
             dx0(n) = zdy(j)
          end do
          do j=1,jmax+1
             n = j + nmin-1
             umid(n)=zuy(i,j,k)
          end do

! Set boundary conditions, compute volume elements, evolve flow, then remap

          call sweepbc( nlefty, nrighty,nmin,nmax,r,umid,xa,xa0,dx,dx0, &
               & dinfloy,dotfloy,vinflo,votflo,thckly,thckry )
          call volume ( ngeomy,nmin,nmax,dvol,dvol0,xa,xa0,dx,dx0,radius )
          if (leuler) then
             call advance( dt,ngeomy,nlefty,nrighty,maxsweep,nmin,nmax, &
                  & r,umid,dvol0,xa0,dx0,flxiny,flxoty )
          else
             call evolve( dt,ngeomy,maxsweep,nmin,nmax,r,umid,dvol,xa,dx,radius )
             call remap ( dt,ngeomy,nlefty,nrighty,maxsweep,nmin,nmax,r, &
                  & dvol,dvol0,xa,xa0,dx,dx0,radius,flxiny,flxoty )
          end if

! Put updated values into 2D arrays, dropping ghost zones

          do j = 1, jmax
             n = j + 6
             zro(i,j,k) = r(n)
          enddo

       end do
    end do

    deallocate(r, umid, xa, xa0, dx, dx0, dvol, dvol0)

    return
  end subroutine sweepy

!----------------------------------------------------------------------

  subroutine sweepbc(nleft,nright,nmin,nmax,r,umid,xa,xa0,dx,dx0, &
       & dinflo,dotflo,uinflo,uotflo,thckl,thckr)

    implicit none
!
! Impose boundary conditions on ends of 1D sweep arrays
!

! GLOBALS
    integer, intent(in) :: nleft,nright, nmin,nmax
    real (kind = dp), intent(in) :: thckl, thckr
    real (kind = dp), intent(in), dimension(:) :: dinflo,dotflo, uinflo,uotflo
    real (kind = dp), intent(inout), dimension(:) :: r, umid, xa, xa0, dx, dx0

! LOCALS
    integer :: n

!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed inflow (eg, uinflo,pinflo,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))
!    = 4 : fixed density with inflow
!    = 5 : fixed flux with inflow

    select case (nleft)
    case (0)
       do n = 1, 6
          dx(nmin-n) = dx(nmin+n-1)
          xa(nmin-n) = xa(nmin-n+1) - dx(nmin-n)
          dx0(nmin-n)= dx0(nmin+n-1)
          xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
          r (nmin-n) = r (nmin+n-1)
       enddo
       do n = 1, 3
          umid(nmin-n) = umid(nmin+n-1)
       enddo
    case (1)
       do n = 1, 6
          dx(nmin-n) = dx(nmin)
          xa(nmin-n) = xa(nmin-n+1) - dx(nmin-n)
          dx0(nmin-n)= dx0(nmin+n-1)
          xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
          r (nmin-n) = r (nmin)
       enddo
       do n = 1, 3
          umid(nmin-n) = umid(nmin)
       enddo
    case (2,5)
       do n = 1, 6
          dx(nmin-n) = dx(nmin)
          xa(nmin-n) = xa(nmin-n+1) - dx(nmin-n)
          dx0(nmin-n)= dx0(nmin+n-1)
          xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
          r (nmin-n) = dinflo(nmin-n)
       enddo
       do n = 1, 3
          umid(nmin-n) = uinflo(4-n)
       enddo
    case (3)
       do n = 1, 6
          dx(nmin-n) = dx(nmax-n)
          xa(nmin-n) = xa(nmin-n+1) - dx(nmin-n)
          dx0(nmin-n)= dx0(nmax-1)
          xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
          r (nmin-n) = r (nmax-n)
       enddo
       do n = 1, 3
          umid(nmin-n) = umid(nmax-n)
       enddo
    case (4)
       do n = 1, 6
          dx(nmin-n) = dx(nmin)
          xa(nmin-n) = xa(nmin-n+1) - dx(nmin-n)
          dx0(nmin-n)= dx0(nmin+n-1)
          xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
          r (nmin-n) = dinflo(nmin-n)
       enddo
       r(nmin) = thckl
       do n = 1, 3
          umid(nmin-n) = uinflo(4-n)
       enddo
    end select

    select case (nright)
    case (0)
       do n = 1, 6
          dx(nmax+n) = dx(nmax+1-n)
          xa(nmax+n) = xa(nmax+n-1) + dx(nmax+n-1)
          dx0(nmax+n)= dx0(nmax+1-n)
          xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
          r (nmax+n) = r (nmax+1-n)
       enddo
       do n = 1, 3
          umid(nmax+1+n) = umid(nmax+2-n)
       enddo
    case (1)
       do n = 1, 6
          dx(nmax+n) = dx(nmax)
          xa(nmax+n) = xa(nmax+n-1) + dx(nmax+n-1)
          dx0(nmax+n)= dx0(nmax+1-n)
          xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
          r (nmax+n) = r (nmax)
       enddo
       do n = 1, 3
          umid(nmax+1+n) = umid(nmax+1)
       enddo
    case (2,5)
       do n = 1, 6
          dx(nmax+n) = dx (nmax)
          xa(nmax+n) = xa (nmax+n-1) + dx(nmax+n-1)
          dx0(nmax+n)= dx0(nmax+1-n)
          xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
          r (nmax+n) = dotflo(n)
       enddo
       do n = 1, 3
          umid(nmax+n+1) = uotflo(n)
       enddo
    case (3)
       do n = 1, 6
          dx(nmax+n) = dx(nmin+n)
          xa(nmax+n) = xa(nmax+n-1) + dx(nmax+n-1)
          dx0(nmax+n)= dx0(nmin+n)
          xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
          r (nmax+n) = r (nmin+n)
       enddo
       do n = 1, 3
          umid(nmax+n+1) = umid(nmin+n+1)
       enddo
    case (4)
       do n = 1, 6
          dx(nmax+n) = dx (nmax)
          xa(nmax+n) = xa (nmax+n-1) + dx(nmax+n-1)
          dx0(nmax+n)= dx0(nmax+1-n)
          xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
          r (nmax+n) = dotflo(n)
       enddo
       r(nmax) = thckr
       do n = 1, 3
          umid(nmax+n+1) = uotflo(n)
       enddo
    end select

    return
  end subroutine sweepbc

!----------------------------------------------------------------------

  subroutine volume(ngeom,nmin,nmax,dvol,dvol0,xa,xa0,dx,dx0,radius)

    implicit none

! GLOBALS
    integer, intent(in) :: ngeom, nmin,nmax
    real (kind = dp), intent(inout) :: radius
    real (kind = dp), intent(in), dimension(:) :: xa, xa0, dx, dx0
    real (kind = dp), intent(inout), dimension(:) :: dvol,dvol0

! LOCALS
    integer :: l
    real (kind = dp), parameter :: third = 1.0_dp/3.0_dp
!
! Calculate zone volume and average area based on geometry of sweep

    select case (ngeom)
    case(0)
       radius = 1.0_dp
       do l = nmin-3, nmax+4
          dvol (l) = dx (l)
          dvol0(l) = dx0(l)
       enddo
    case(1)
       radius = 1.0_dp
       do l = nmin-3, nmax+4
          dvol (l) = dx (l)*(xa (l)+0.5_dp*dx (l))
          dvol0(l) = dx0(l)*(xa0(l)+0.5_dp*dx0(l))
       enddo
    case(2)
       radius = 1.0_dp
       do l = nmin-3, nmax+4
          dvol (l)=dx (l)*(xa (l)*(xa (l)+dx (l))+dx (l)*dx (l)*third)
          dvol0(l)=dx0(l)*(xa0(l)*(xa0(l)+dx0(l))+dx0(l)*dx0(l)*third)
       enddo
    case(3)
       do l = nmin-3, nmax+4
          dvol (l) = dx (l)*radius
          dvol0(l) = dx0(l)*radius
       enddo
    case(4)
       do l = nmin-3, nmax+4
          dvol (l) = (cos(xa (l))-cos(xa (l+1)))*radius
          dvol0(l) = (cos(xa0(l))-cos(xa0(l+1)))*radius
       enddo
    case(5)
       do l = nmin-3, nmax+4
          dvol (l) = dx (l) * radius
          dvol0(l) = dx0(l) * radius
       enddo
    case default
       write(*,*) 'Geometry', ngeom, ' not implemented.'
       stop
    end select

    return
  end subroutine volume

!----------------------------------------------------------------------

  subroutine advance( dt,ngeom,nleft,nright,maxsweep,nmin,nmax,r,umid,dvol0,xa0,dx0, &
       & flxin,flxot )

    implicit none

! GLOBALS
    integer, intent(in) :: ngeom, nleft,nright, maxsweep,nmin,nmax
    real (kind = dp), intent(in) :: dt, flxin, flxot
    real (kind = dp), intent(in), dimension(:) :: umid, dvol0, xa0, dx0
    real (kind = dp), intent(inout), dimension(:) :: r

! Piecewise-parabolic interpolation coeffs and parameters
    real(kind = dp) :: fractn, fractn2, deltx
    real (kind = dp), parameter :: third = 1.0_dp/3.0_dp
    real (kind = dp), parameter :: fourthd= 4.0_dp * third
    real(kind = dp), allocatable, dimension(:) :: delta, dr, rl, r6
    real(kind = dp), allocatable, dimension(:,:) :: para

    real(kind = dp), allocatable, dimension(:) :: fluxr

! LOCALS
    integer :: n, k
!
!
    allocate(delta(maxsweep), dr(maxsweep), rl(maxsweep), r6(maxsweep))
    allocate(fluxr(maxsweep), para(5,maxsweep))

    delta=0.0_dp
    dr=0.0_dp
    rl=0.0_dp
    r6=0.0_dp
    fluxr=0.0_dp
    para=0.0_dp
!
! Generate interpolation functions
!
    call paraset( para, dx0, xa0, nmin-1, nmax+1, maxsweep )
    call parabola(nmin-1,nmax+1,maxsweep,para,r,dr,r6,rl)

! 
! Calculate change in flux during the timestep
!
    if(ngeom.eq.0) then
       do n = nmin, nmax+1
          delta(n) = umid(n) * dt
       end do
    else
       print *, 'Eulerian advance: Limited geometry available.'
       print *, 'Modification required'
       stop
    endif
!
    do n = nmin, nmax + 1
       deltx = umid(n) * dt

       if(deltx .ge. 0.0_dp) then

          k = n - 1
          fractn  = 0.5_dp*deltx/dx0(k)
          fractn2 = 1.0_dp - fourthd*fractn
          fluxr(n) = delta(n)*( &
               &   rl(k)+dr(k)-fractn*(dr(k)-fractn2*r6(k)))
       else

          fractn   = 0.5_dp*deltx/dx0(n)
          fractn2  = 1.0_dp + fourthd*fractn
          fluxr(n) = delta(n)*( &
               &   rl(n) - fractn*(dr(n)+fractn2*r6(n)))

       endif
    end do
! Applying flux boundary conditions
! CcC Vicky 08/08. Not tested
    if (nleft  == 5) fluxr(nmin)   = dt * flxin
    if (nright == 5) fluxr(nmax+1) = dt * flxot

!
! Advance the continuity equation
!
    do n = nmin, nmax
       r(n) = r(n) + (fluxr(n) - fluxr(n+1))/dx0(n)
       r(n) = max(0.0_dp,r(n))
    end do
!
    deallocate(delta, dr, rl, r6, fluxr, para)

    return
  end subroutine advance

!----------------------------------------------------------------------

  subroutine evolve( dt,ngeom,maxsweep,nmin,nmax,r,umid,dvol,xa,dx,radius )
!
! Use umid to update density
! Physical zones are from nmin to nmax. Zone boundary numbers run from
! nmin to nmax+1
!
    implicit none

! GLOBALS
    integer, intent(in) :: ngeom,maxsweep,nmin,nmax
    real (kind = dp), intent(in) :: dt, radius
    real (kind = dp), intent(in), dimension(:) :: umid
    real (kind = dp), intent(inout), dimension(:) :: r, dvol, xa, dx

    real(kind = dp), allocatable, dimension(:) :: xa1, dvol1  ! Original Lagrangian grid

! LOCALS
    integer :: n
    real (kind = dp), parameter :: third = 1.0_dp/3.0_dp

    allocate(xa1(maxsweep), dvol1(maxsweep))
    xa1=0.0_dp
    dvol1=0.0_dp
!
! grid position evolution
!
    do n = nmin-3, nmax + 4
       xa1  (n) = xa(n)
       dvol1(n) = dvol(n)
       xa   (n) = xa(n) + dt * umid(n) / radius
    end do
    xa1(nmin-4) = xa(nmin-4)
    xa1(nmax+5) = xa(nmax+5)
    do n = nmin-4, nmax+5
       dx (n)   = xa(n+1) - xa(n)
    end do
!
! Calculate dvolume and average area based on geometry of sweep
!
    select case (ngeom)
    case(0)
       do n = nmin-3, nmax+4
          dvol(n) = dx(n)
       end do
    case(1)
       do n = nmin-3, nmax+4
          dvol(n) = dx(n)*(xa(n)+0.5_dp*dx(n))
       end do
    case(2)
       do n = nmin-3, nmax+4
          dvol(n) = dx(n)*(xa(n)*(xa(n)+dx(n))+dx(n)*dx(n)*third)
       end do
    case(3)
       do n = nmin-3, nmax+4
          dvol(n) = dx(n)*radius
       end do
    case(4)
       do n = nmin-3, nmax+4
          dvol(n) = (cos(xa(n))-cos(xa(n+1)))*radius
       end do
    case(5)
       do n = nmin-3, nmax+4
          dvol(n) = dx(n) * radius
       end do
    case default
       write(*,*) 'Geometry', ngeom, ' not implemented.'
       stop
    end select
!     
    do n = nmin-3, nmax+3
!
! density evolution. lagrangian code, so all we have to do is watch the
! change in the geometry.
       r(n) = r(n) * ( dvol1(n) / dvol(n) )
       r(n) = max(r(n),0.0_dp)

    end do
!
    deallocate(xa1 ,dvol1)

    return
  end subroutine evolve

!----------------------------------------------------------------------

  subroutine remap( dt,ngeom,nleft,nright,maxsweep,nmin,nmax,r,dvol,dvol0, &
       & xa,xa0,dx,dx0,radius,flxin,flxot )
!
! Remap mass, momentum, and energy from the updated lagrangian grid
! to the fixed Eulerian grid, using piecewise parabolic functions.
! No flattening is used on remap (pass dummy array: dum to parabola.f).
!
    implicit none

! GLOBALS
    integer, intent(in) :: ngeom, nleft,nright, maxsweep,nmin,nmax
    real (kind = dp), intent(in) :: dt, flxin, flxot,radius
    real (kind = dp), intent(in), dimension(:) :: dvol,dvol0,xa,xa0,dx,dx0
    real (kind = dp), intent(inout), dimension(:) :: r

! Piecewise-parabolic interpolation coeffs and parameters
    real(kind = dp) :: fractn, fractn2, deltx
    real (kind = dp), parameter :: third = 1.0_dp/3.0_dp
    real (kind = dp), parameter :: fourthd= 4.0_dp * third
    real(kind = dp), allocatable, dimension(:) :: delta, dr, rl, r6
    real(kind = dp), allocatable, dimension(:,:) :: para

    real(kind = dp), allocatable, dimension(:) :: dm, dm0, fluxr

! LOCALS
    integer :: n, k
!
!                                       
!
    allocate(delta(maxsweep), dr(maxsweep), rl(maxsweep), r6(maxsweep))
    allocate(dm(maxsweep), dm0(maxsweep), fluxr(maxsweep), para(5,maxsweep))

    delta=0.0_dp
    dr=0.0_dp
    rl=0.0_dp
    r6=0.0_dp
    dm=0.0_dp
    dm0=0.0_dp
    fluxr=0.0_dp
    para=0.0_dp
!
! Generate interpolation functions
!
    call paraset( para, dx, xa, nmin-1, nmax+1, maxsweep )
    call parabola(nmin-1,nmax+1,maxsweep, para,r,dr,r6,rl)

! Calculate the volume of the overlapping subshells (delta)

    if(ngeom.eq.0) then
       do n = nmin, nmax+1
          delta(n) = xa(n) - xa0(n)
       end do
    else if(ngeom.eq.1) then
       do n = nmin, nmax+1
          delta(n) = xa(n) - xa0(n)
          delta(n) = delta(n)*(xa0(n)+0.5_dp*delta(n))
       end do
    else if(ngeom.eq.2) then
       do n = nmin, nmax+1
          delta(n) = xa(n) - xa0(n)
          delta(n) = delta(n)*(xa0(n)*(xa0(n) &
               &               + delta(n)) + delta(n)**2*third)
       end do
    else if(ngeom.eq.3) then
       do n = nmin, nmax+1
          delta(n) = (xa(n) - xa0(n)) * radius
       end do
    else if(ngeom.eq.4) then
       do n = nmin, nmax+1
          delta(n) = (cos(xa0(n)) - cos(xa(n))) * radius
       end do
    else if(ngeom.eq.5) then
       do n = nmin, nmax+1
          delta(n) = (xa(n) - xa0(n)) * radius
       end do
    endif
!
! Calculate the total mass (fluxr)
! in the subshell created by the overlap of the Lagrangian and Eulerican grids.
! If the zone face has moved to the left (deltx > 0), use the integral from the
! left side of zone n (fluxrr).  If the zone face has moved to the right 
! (deltx < 0), use the integral from the right side of zone k=n-1 (fluxrl).
!
    do n = nmin, nmax + 1
       deltx = xa(n) - xa0(n)

       if(deltx .ge. 0.0_dp) then

          k = n - 1
          fractn  = 0.5_dp*deltx/dx(k)
          fractn2 = 1.0_dp - fourthd*fractn
          fluxr(n) = delta(n)*( &
               &   rl(k)+dr(k)-fractn*(dr(k)-fractn2*r6(k)))
       else

          fractn   = 0.5_dp*deltx/dx(n)
          fractn2  = 1.0_dp + fourthd*fractn
          fluxr(n) = delta(n)*( &
               &   rl(n) - fractn*(dr(n)+fractn2*r6(n)))

       endif
    end do
! Applying flux boundary conditions
! CcC Vicky 08/08. Not tested
    if (nleft  == 5) fluxr(nmin)   = dt * flxin
    if (nright == 5) fluxr(nmax+1) = dt * flxot

!
! Advect mass by moving the subshell quantities 
! into the appropriate Eulerian zone. 
!
    do n = nmin, nmax
       dm (n) = r(n) * dvol(n)
       dm0(n) = (dm(n) + fluxr(n) - fluxr(n+1))
       r  (n) = dm0(n)/dvol0(n)
       r  (n) = max(0.0_dp,r(n))
    end do
!
    deallocate(delta, dr, rl, r6, dm, dm0, fluxr, para)

    return
  end subroutine remap

!----------------------------------------------------------------------

  subroutine paraset( para, dx, xa, nmin, nmax, maxsweep )
!
! Colella and Woodward, JCompPhys 54, 174-201 (1984) eq 1.6, 1.7
!
! paraset sets up constants which are re-used each time we want to
! interpolate a parabola on a quantity. First pull out constants
! A, B, and C, and then compute all the equations in terms of those
! quantities. 
!
! the quantities calculated here are stored in a array para,
! to be read by parabola()
!
! nmin/nmax are index range for which one will calculate parabolae

    implicit none

    integer, intent(in) :: nmin, nmax, maxsweep
    real(kind = dp), dimension(:,:), intent(inout) :: para
    real(kind = dp), dimension(:), intent(in) :: dx, xa

    integer :: n
    real(kind = dp), allocatable, dimension(:) :: a, b, c, d, ai, bi, ci

    allocate(a(maxsweep), b(maxsweep), c(maxsweep), d(maxsweep))
    allocate(ai(maxsweep), bi(maxsweep), ci(maxsweep))
    a=0.0_dp
    b=0.0_dp
    c=0.0_dp
    d=0.0_dp
    ai=0.0_dp
    bi=0.0_dp
    ci=0.0_dp
!                                        A =  dX_j +  dX_j+1
!                                        B = 2dX_j +  dX_j+1
!                                        C =  dX_j + 2dX_j+1
!
!                                        ai, bi, and ci are inverse quantities
    do n = nmin-2, nmax+1
       a(n) = dx(n) + dx(n+1)
       ai(n) = 1.0_dp/a(n)
       b(n) = a(n) + dx(n)
       bi(n) = 1.0_dp/b(n)
       c(n) = a(n) + dx(n+1)
       ci(n) = 1.0_dp/c(n)
    enddo

!                                        constants for equation 1.6
!     a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * da(j+1) + C3 * da(j)

    do  n = nmin-1, nmax
       d(n)      = 1.0_dp / (a(n-1) + a(n+1))
       para(1,n) = dx(n) * ai(n) + 2.0_dp * dx(n+1) * dx(n) * &
            &         d(n) * ai(n) * ( a(n-1) * bi(n) - a(n+1) * ci(n) )
       para(2,n) = - d(n) * dx(n)   * a(n-1) * bi(n)
       para(3,n) =   d(n) * dx(n+1) * a(n+1) * ci(n)
    enddo

!                                        constants for equation 1.7
!     da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))

    do n = nmin-1, nmax+1
       d(n) = dx(n) / ( a(n-1) + dx(n+1) )
       para(4,n) = d(n) * b(n-1) * ai(n)
       para(5,n) = d(n) * c(n)   * ai(n-1)
    enddo
!

    deallocate(a, b, c, d, ai, bi, ci)

    return
  end subroutine paraset

!----------------------------------------------------------------------

  subroutine parabola( nmin, nmax, maxsweep, para, a, deltaa, a6, al)
!
! Colella and Woodward, JCompPhys 54, 174-201 (1984) eq 1.5-1.8,1.10
!
! parabola calculates the parabolas themselves. call paraset first
! for a given grid-spacing to set the constants, which can be reused
! each time parabola is called.
!
! flatening coefficients are calculated externally in flatten. 
!
! nmin/nmax are indicies over which the parabolae are calculated

    implicit none

    integer, intent(in) :: nmin, nmax, maxsweep
    real(kind = dp), dimension(:,:), intent(in) :: para
    real(kind = dp), dimension(:), intent(in) :: a
    real(kind = dp), dimension(:), intent(inout) :: deltaa, a6, al

    integer :: n
    real(kind = dp), allocatable, dimension(:) :: da, ar, diffa, scrch1, scrch2

    allocate(da(maxsweep), ar(maxsweep), diffa(maxsweep))
    allocate(scrch1(maxsweep), scrch2(maxsweep))
    da=0.0_dp
    ar=0.0_dp
    diffa=0.0_dp
    scrch1=0.0_dp
    scrch2=0.0_dp
!
!
    do n = nmin-2, nmax+1
       diffa(n) = a(n+1) - a(n)
    enddo

!                                                       Equation 1.7 of C&W
!     da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))
!     zero out da(n) if a(n) is a local max/min

    do n = nmin-1, nmax+1
       if(diffa(n-1)*diffa(n) .lt. 0.0_dp) then
          da(n) = 0.0_dp
       else
          da(n) = para(4,n) * diffa(n) + para(5,n) * diffa(n-1)
          da(n) = sign( min(abs(da(n)), 2.0_dp*abs(diffa(n-1)), &
               &    2.0_dp*abs(diffa(n))), da(n) )
       endif
    enddo

!                                                       Equation 1.6 of C&W
!     a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * dma(j+1) + C3 * dma(j)
! MONOT: Limit ar(n) to the range defined by a(n) and a(n+1)

    do n = nmin-1, nmax
       ar(n) = a(n) + para(1,n)*diffa(n) + para(2,n)*da(n+1) + &
            &    para(3,n)*da(n)
       al(n+1) = ar(n)
    enddo

! MONOTONICITY constraints:

! compute delta_a, a_6
! MONOT: if a is a local max/min, flaten zone structure ar,al -> a.
! MONOT: compute monotonized values using eq. 1.10 of C&W
!        if parabola exceeds al/ar, reset ar/al so that slope -> 0.
! Recalculate delta_a and a_6

    do n = nmin, nmax

       deltaa(n)= ar(n) - al(n)
       a6(n)    = 6.0_dp * (a(n) - 0.5_dp * (al(n) + ar(n)))

       scrch1(n) = (ar(n) - a(n)) * (a(n)-al(n))
       if(scrch1(n) .lt. 0.0_dp) then
          ar(n) = a(n)
          al(n) = a(n)
       endif

       scrch1(n) = deltaa(n) * deltaa(n)
       scrch2(n) = deltaa(n) * a6(n)
       if(scrch1(n) .lt. +scrch2(n)) al(n) = 3.0_dp * a(n) - 2.0_dp * ar(n)
       if(scrch1(n) .lt. -scrch2(n)) ar(n) = 3.0_dp * a(n) - 2.0_dp * al(n)

       deltaa(n)= ar(n) - al(n)
       a6(n)    = 6.0_dp * (a(n) - 0.5_dp * (al(n) + ar(n)))

    enddo

    deallocate(da, ar, diffa, scrch1, scrch2)

    return
  end subroutine parabola

!----------------------------------------------------------------------

end module ppm
