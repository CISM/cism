
! "glam_strs2.F90"
!
! 3d velocity calculation based on Blatter/Pattyn, 1st-order equations, by Tony Payne (Univ.
! of Bristol) and Steve Price (Univ. of Bristol / Los Alamos Nat. Lab.). Boundary conditions
! available include periodic (lateral), free surface, zero slip at bed, specified basal 
! traction at bed, and specified basal yield stress at bed (all three of which are implemented
! through various verions of the specified traction b.c.)
! include macros for glide mask definitions

#include "glide_mask.inc"       

!***********************************************************************

module glam_strs2

!***********************************************************************

use glimmer_paramets, only : dp
use glimmer_physcon,  only : gn, rhoi, rhoo, grav, pi, scyr
use glimmer_paramets, only : thk0, len0, vel0, vis0, vis0_glam, tim0, lambda0, evs0, tau0_glam
use glimmer_log,      only : write_log
use glide_mask
use glimmer_sparse_type
use glimmer_sparse

implicit none

  integer, save :: locplusup
  logical, save :: lateralboundry = .false.
  integer, dimension(6), save :: loc_latbc

  real (kind = dp), allocatable, dimension(:,:,:),     save :: flwafact
  real (kind = dp), allocatable, dimension(:),         save :: dups
  real (kind = dp), allocatable, dimension(:,:,:,:,:), save :: corr
  real (kind = dp), allocatable, dimension(:,:,:,:),   save :: usav
  real (kind = dp), allocatable, dimension(:,:,:),     save :: tvel
  real (kind = dp), allocatable, dimension(:),         save :: dup, dupm

  integer, dimension(:,:), allocatable :: uindx
  integer, dimension(:,:), allocatable :: umask 

  ! regularization constant for eff. strain rate to avoid infinite visc.
  ! NOTE: would be good to explore how small this really needs to be, as 
  ! code converges much better when this value is made larger.
  real (kind = dp), parameter :: effstrminsq = (1.0e-20_dp * tim0)**2

  real (kind = dp) :: p1, p2    ! variants of Glen's "n" (e.g. n, (1-n)/n)
  real (kind = dp) :: dew2, dns2, dew4, dns4

  ! combinations of coeffs. used in momentum balance calcs
  real (kind = dp) :: cdxdy
  real (kind = dp), dimension(2) :: cdxdx
  real (kind = dp), dimension(:),   allocatable :: cdsds, cds
  real (kind = dp), dimension(:), allocatable :: cvert, fvert
  real (kind = dp), dimension(:,:), allocatable :: cdsdx

  real (kind = dp), dimension(:), allocatable :: dsigmadew, dsigmadns
  real (kind = dp), dimension(:), allocatable :: d2sigmadew2, d2sigmadns2, d2sigmadewdns
  real (kind = dp) :: d2sigmadewdsigma, d2sigmadnsdsigma

  ! vectors of coeffs. used for switching symmetric solution subroutines between calc.
  ! of x-comp of vel or y-comp of vel
  real (kind = dp), dimension(2), parameter ::   &
           oneorfour = (/ 1.0_dp, 4.0_dp /),     &
           fourorone = (/ 4.0_dp, 1.0_dp /),     &
           oneortwo  = (/ 1.0_dp, 2.0_dp /),     &
           twoorone  = (/ 2.0_dp, 1.0_dp /)

  ! coeff. for forward differencing template, used for stress bcs at lateral boundaries
  real (kind = dp), dimension(3), parameter ::   &
           onesideddiff = (/ -3.0_dp, 4.0_dp, -1.0_dp /)

  ! geometric 2nd and cross-derivs
  real (kind = dp), dimension(:,:), allocatable :: &
              d2thckdew2, d2usrfdew2, d2thckdns2, d2usrfdns2, d2thckdewdns, d2usrfdewdns

  ! variables for use in sparse matrix calculation
  real (kind = dp), dimension(:), allocatable :: pcgval, rhsd
  integer, dimension(:), allocatable :: pcgcol, pcgrow
  integer, dimension(2) :: pcgsize
  integer :: ct

  integer, parameter :: unin = 90


!***********************************************************************

contains

!***********************************************************************


subroutine glam_velo_fordsiapstr_init( ewn,   nsn,   upn,    &
                                       dew,   dns,           &
                                       sigma, stagsigma )

    ! Allocate arrays and initialize variables.
    implicit none

    integer, intent(in) :: ewn, nsn, upn
    real (kind = dp), intent(in) :: dew, dns

    real (kind = dp), dimension(:), intent(in)  :: sigma
    real (kind = dp), dimension(:), intent(out)  :: stagsigma

    integer :: up

    allocate( dup(upn) )
    allocate( dupm(upn) )
    allocate( cvert(upn) )
    allocate( cdsdx(upn,2) )
    allocate( cdsds(upn) )
    allocate( cds(upn) )
    allocate( fvert(upn) )

    ! NOTE: "dup", the sigma coordinate spacing is defined as a vector to allow it to 
    ! be read in from file for use with non-constant vertical grid spacing. Currently, this
    ! is not working, so the code will not give accurate results if the sigma coordinate is
    ! not regularly spaced. - not working!!) 
    dup = (/ ( (sigma(2)-sigma(1)), up = 1, upn) /) 
    dupm = - 0.25_dp / dup
    stagsigma = (sigma(1:upn-1) + sigma(2:upn)) / 2.0_dp

    ! p1 = -1/n   - used with rate factor in eff. visc. def.
    ! p2 = (1-n)/2n   - used with eff. strain rate in eff. visc. def. 
    p1 = -1.0_dp / real(gn,dp)      
    p2 = (1.0_dp - real(gn,dp)) / (2.0_dp * real(gn,dp))

    dew2 = 2.0_dp * dew; dns2 = 2.0_dp * dns        ! 2x the standard grid spacing
    dew4 = 4.0_dp * dew; dns4 = 4.0_dp * dns        ! 4x the standard grid spacing

    allocate(dsigmadew(upn),  dsigmadns(upn))
    allocate(d2sigmadew2(upn),d2sigmadns2(upn),d2sigmadewdns(upn))

    allocate (d2thckdew2(ewn-1,nsn-1),d2thckdns2(ewn-1,nsn-1),d2thckdewdns(ewn-1,nsn-1), &
              d2usrfdew2(ewn-1,nsn-1),d2usrfdns2(ewn-1,nsn-1),d2usrfdewdns(ewn-1,nsn-1))

    allocate(umask(ewn-1,nsn-1))                        

    allocate(flwafact(1:upn-1,ewn,nsn))  ! NOTE: the vert dim here must agree w/ that of 'efvs'

    allocate(dups(upn)) 

    flwafact = 0.0_dp

     ! define constants used in various FD calculations associated with the 
     ! subroutine 'findcoefst'   
     call calccoeffsinit(upn, dew, dns)

    dups = (/ (sigma(up+1) - sigma(up), up=1,upn-1), 0.0d0 /)

end subroutine glam_velo_fordsiapstr_init


!***********************************************************************

! Note that this is the driver subroutine, called from 'run_ho_diagnostic' in
! 'glide_velo_higher.F90'. In turn, 'run_ho_model' is called from 'inc_remap_driver' in
! 'glam.F90', and 'inc_remap_driver' is called from 'glide_tstep_ps' in 'glide.F90'.

subroutine glam_velo_fordsiapstr(ewn,      nsn,    upn,  &
                                 dew,      dns,          &
                                 sigma,    stagsigma,    &
                                 thck,     usrf,         &
                                 lsrf,     topg,         &
                                 dthckdew, dthckdns,     &
                                 dusrfdew, dusrfdns,     & 
                                 dlsrfdew, dlsrfdns,     &
                                 stagthck, flwa,         & 
                                 mintauf,                & 
                                 umask,                  & 
                                 whichbabc,              &
                                 whichefvs,              &
                                 whichresid,             &
                                 whichsparse,            &
                                 periodic_ew,periodic_ns,&
                                 beta,                   & 
                                 uvel,     vvel,         &
                                 uflx,     vflx,         &
                                 efvs )

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:),   intent(inout)  :: umask  
  ! NOTE: 'inout' status to 'umask' should be changed to 'in' at some point, 
  ! but for now this allows for some minor internal hacks to CISM-defined mask  

  real (kind = dp), intent(in) :: dew, dns

  real (kind = dp), dimension(:),     intent(in)  :: sigma, stagsigma       ! sigma coords
  real (kind = dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg ! geom vars
  real (kind = dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns     ! thick grads
  real (kind = dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns     ! upper surf grads
  real (kind = dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns     ! basal surf grads
  real (kind = dp), dimension(:,:),   intent(in)  :: stagthck               ! staggered thickness
  real (kind = dp), dimension(:,:),   intent(in)  :: minTauf                ! till yield stress
  real (kind = dp), dimension(:,:,:), intent(in)  :: flwa                   ! flow law rate factor

  ! This is the betasquared field from CISM (externally specified), and should eventually
  ! take the place of the subroutine 'calcbetasquared' below. For now, there is simply an option
  ! in the subroutine 'calcbetasquared' (case 9) to use this external, CISM specified value for
  ! the betasquared field as opposed to one of the values calculated internally.
  real (kind = dp), dimension(:,:),   intent(in)  :: beta 

  integer, intent(in) :: whichbabc    ! options for betasquared field to use
  integer, intent(in) :: whichefvs    ! options for efvs calculation (calculate it or make it uniform)
  integer, intent(in) :: whichresid   ! options for method to use when calculating vel residul
  integer, intent(in) :: whichsparse  ! options for which method for doing elliptic solve
  logical, intent(in) :: periodic_ew, periodic_ns  ! options for applying periodic bcs or not

  real (kind = dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)
  real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx  ! horiz fluxs: u_bar*H, v_bar*H
  real (kind = dp), dimension(:,:,:), intent(out) :: efvs        ! effective viscosity

  integer :: ew, ns, up     ! counters for horiz and vert do loops

  real (kind = dp), parameter :: minres = 1.0d-4    ! assume vel fields converged below this resid 
  real (kind = dp), save, dimension(2) :: resid     ! vector for storing u resid and v resid 

  integer, parameter :: cmax = 3000                 ! max no. of iterations
  integer :: counter                                ! iteation counter 
  character(len=100) :: message                     ! error message

  ! variables used for incorporating generic wrapper to sparse solver
  type(sparse_matrix_type) :: matrix
  real (kind = dp), dimension(:), allocatable :: answer
  real (kind = dp) :: err
  integer :: iter


  ! calc geometric 2nd deriv. for generic input variable 'ipvr', returns 'opvr'
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! calc geometric 2nd cross-deriv. for generic input variable 'ipvr', returns 'opvr'
  call geom2derscros(dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(dew, dns, usrf, stagthck, d2usrfdewdns)

  allocate(uindx(ewn-1,nsn-1))

  ! If a point from the 2d array 'mask' is associated with a non-zero ice thickness
  ! assign it a unique number. If not assign a zero.			 
  uindx = indxvelostr(ewn, nsn, upn, umask,pcgsize(1))


!! A hack of the boundary condition mask needed for the Ross Ice Shelf exp.
!! The quick check of whether or not this is the Ross experiment is to look
!! at the domain size.
 if( ewn == 151 .and. nsn == 115 )then
    do ns=1,nsn-1; do ew=1,ewn-1
        if( umask(ew,ns) == 21 .or. umask(ew,ns) == 5 )then
            umask(ew,ns) = 73
        endif
    end do; end do
 end if


  ! allocate space for storing temporary across-flow comp of velocity
  allocate(tvel(upn,ewn-1,nsn-1)) 
  tvel = 0.0_dp
 
  ! allocate space for variables used by 'mindcrash' function (unstable manifold correction)
  allocate(corr(upn,ewn-1,nsn-1,2,2),usav(upn,ewn-1,nsn-1,2))

  ! make an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 20

  ! allocate space matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)),rhsd(pcgsize(1)), &
            pcgval(pcgsize(2)))
  allocate(matrix%row(pcgsize(2)), matrix%col(pcgsize(2)), &
            matrix%val(pcgsize(2)), answer(pcgsize(1)))

  ! set residual and iteration counter to initial values
  resid = 1.0_dp    
  counter = 1

  ! print some info to the screen to update on iteration progress
  print *, ' '
  print *, 'Running Payne/Price higher-order dynamics solver'
  print *, ' '
  print *, 'iter #     uvel resid          vvel resid         target resid'
  print *, ' '

  ! ****************************************************************************************
  ! START of Picard iteration
  ! ****************************************************************************************
 
  ! Picard iteration; continue iterating until resid falls below specified tolerance
  ! or the max no. of iterations is exceeded
  do while ( maxval(resid) > minres .and. counter < cmax)
  !do while ( resid(1) > minres .and. counter < cmax)  ! used for 1d solutions where d*/dy=0 

    ! calc effective viscosity using previously calc vel. field
    call findefvsstr(ewn,  nsn,  upn,      &
                     stagsigma,  counter,    &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     umask)

    ! calculate coeff. for stress balance in y-direction 
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     vvel,        uvel,           &
                     thck,        dusrfdns,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta, counter )

    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, vvel )

    ! solve 'Ax=b' for the y-component of velocity using method "which_sparse"
    call sparse_easy_solve( matrix, rhsd, answer, err, iter, whichsparse )

    ! put vels and coeffs from sparse vector format (soln) back into 3d arrays
    call solver_postprocess( ewn, nsn, upn, uindx, answer, tvel )

    ! NOTE: y-component of velocity that comes out is called "tvel", to differentiate it
    ! from the y-vel solution from the previous iteration, which is maintained as "vvel". 
    ! This is necessary since we have not yet solved for the x-comp of vel, which needs the
    ! old prev. guess as an input (NOT the new guess).


! implement periodic boundary conditions in y (if flagged)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( periodic_ns )then
        tvel(:,:,nsn-1) = tvel(:,:,2)
        tvel(:,:,1) = tvel(:,:,nsn-2)
    end if
    if( periodic_ew )then
        tvel(:,ewn-1,:) = tvel(:,2,:)
        tvel(:,1,:) = tvel(:,ewn-2,:)
    end if
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


    ! calculate coeff. for stress balance calc. in x-direction 
    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        vvel,           &
                     thck,        dusrfdew,       &
                     dusrfdew,    dthckdew,       &
                     d2usrfdew2,  d2thckdew2,     &
                     dusrfdns,    dthckdns,       &
                     d2usrfdns2,  d2thckdns2,     &
                     d2usrfdewdns,d2thckdewdns,   &
                     dlsrfdew,    dlsrfdns,       &
                     stagthck,    whichbabc,      &
                     uindx,       umask,          &
                     lsrf,        topg,           &
                     minTauf,     flwa,           &
                     beta, counter )


    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, uvel )

    ! solve 'Ax=b' for the x-component of velocity using method "which_sparse"
    call sparse_easy_solve( matrix, rhsd, answer, err, iter, whichsparse )

    ! put vels and coeffs from sparse vector format (soln) back into 3d arrays
    call solver_postprocess( ewn, nsn, upn, uindx, answer, uvel )


    ! apply unstable manifold correction to converged velocities
    uvel = mindcrshstr(1,whichresid,uvel,counter,resid(1))
    vvel = mindcrshstr(2,whichresid,tvel,counter,resid(2))


! implement periodic boundary conditions in x (if flagged)
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if( periodic_ns )then
        uvel(:,:,nsn-1) = uvel(:,:,2)
        uvel(:,:,1) = uvel(:,:,nsn-2)
    end if
    if( periodic_ew )then
        uvel(:,ewn-1,:) = uvel(:,2,:)
        uvel(:,1,:) = uvel(:,ewn-2,:)
    end if
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    counter = counter + 1   ! advance the iteration counter

    ! output the iteration status: iteration number, max residual, and location of max residual
    ! (send output to the screen or to the log file, per whichever line is commented out) 
    print '(i4,3g20.6)', counter, resid(1), resid(2), minres
    !write(message,'(" * strs ",i3,3g20.6)') counter, resid(1), resid(2), minres
    !call write_log (message)

  end do

  ! ****************************************************************************************
  ! END of Picard iteration
  ! ****************************************************************************************

  do ns = 1,nsn-1
      do ew = 1,ewn-1 
      ! calc. fluxes from converged vel. fields (needed for input to thickness evolution subroutine)
         if (umask(ew,ns) > 0) then
             uflx(ew,ns) = vertintg(upn, sigma, uvel(:,ew,ns)) * stagthck(ew,ns)
             vflx(ew,ns) = vertintg(upn, sigma, vvel(:,ew,ns)) * stagthck(ew,ns)
         end if
      end do
  end do


  ! de-allocation sparse matrix solution variables 
  deallocate(tvel)
  deallocate(uindx,corr,usav)
  deallocate(pcgval,pcgrow,pcgcol,rhsd)
  deallocate(matrix%row, matrix%col, matrix%val)
  deallocate(answer) 

  return

end subroutine glam_velo_fordsiapstr

!***********************************************************************

function indxvelostr(ewn,  nsn,  upn,  &
                     mask, pointno)
!if a point from the 2d array 'mask' is associated with non-zero ice thickness, 
! (either a boundary or interior point) give it a unique number. If not, give it a zero.

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in), dimension(:,:) :: mask
  integer, intent(out) :: pointno

  integer :: ew, ns
  integer, dimension(size(mask,1),size(mask,2)) :: indxvelostr

  pointno = 1

  do ew = 1,ewn-1
      do ns = 1,nsn-1
        if ( GLIDE_HAS_ICE( mask(ew,ns) ) ) then 
          indxvelostr(ew,ns) = pointno
          pointno = pointno + 1
        else
          indxvelostr(ew,ns) = 0
        end if
      end do
  end do

  ! add two ghost points at upper and lower boundaries (needed for sfc and basal bcs)
  pointno = (pointno - 1) * (upn + 2)

  return

end function indxvelostr

!***********************************************************************

subroutine findefvsstr(ewn,  nsn, upn,       &
                       stagsigma, counter,     &
                       whichefvs, efvs,      &
                       uvel,      vvel,      &
                       flwa,      thck,      &
                       dusrfdew,  dthckdew,  &
                       dusrfdns,  dthckdns,  &
                       mask)

  ! calculate the eff. visc.	
  implicit none 

  integer, intent(in) :: ewn, nsn, upn 
  real (kind = dp), intent(in), dimension(:)     :: stagsigma
  real (kind = dp), intent(in), dimension(:,:,:) :: uvel, vvel, flwa
  real (kind = dp), intent(inout), dimension(:,:,:) :: efvs
  real (kind = dp), intent(in), dimension(:,:) :: thck, dthckdew, dusrfdew, & 
                                                  dusrfdns, dthckdns
  integer, intent(in), dimension(:,:) :: mask
  integer, intent(in) :: whichefvs, counter
       
  integer :: ew, ns, up

  real (kind = dp), dimension(size(efvs,1)) :: effstr, ugradup, vgradup, &
                                               ugradew, ugradns, vgradew, vgradns

  integer, dimension(2) :: mew, mns

  ! This is the factor 1/4(X0/H0)^2 in front of the term ((dv/dz)^2+(du/dz)^2) 
  real (kind = dp), parameter :: f1 = 0.25_dp * (len0 / thk0)**2

  select case(whichefvs)

  case(0)       ! calculate eff. visc. using eff. strain rate

  if (1 == counter) then
    do ns = 2,nsn-1; do ew = 2,ewn-1
    if (thck(ew,ns) > 0.0_dp) then
      ! this is the rate factor term in the expression for the eff. visc: 1/2*A^(-1/n),
      ! which is averaged to midpoints in the vertical (i.e. it lives on a staggered 
      ! grid in the vertical, which is the case for "efvs" as well).
      forall (up = 1:upn-1) flwafact(up,ew,ns) = 0.5_dp * (sum(flwa(up:up+1,ew,ns)) / 2.0_dp)**p1
    end if; end do; end do
  end if

  do ns = 2,nsn-1
      do ew = 2,ewn-1
        if (thck(ew,ns) > 0.0_dp) then

            ugradup = vertideriv(upn, hsum(uvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))
            vgradup = vertideriv(upn, hsum(vvel(:,ew-1:ew,ns-1:ns)), thck(ew,ns))

            ugradew = horizderiv(upn,  stagsigma,        &
                         sum(uvel(:,ew-1:ew,ns-1:ns),3), &
                         dew4, ugradup,                  &             
                         sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                         sum(dthckdew(ew-1:ew,ns-1:ns)))
        
            vgradew = horizderiv(upn,  stagsigma,        &
                         sum(vvel(:,ew-1:ew,ns-1:ns),3), &
                         dew4, vgradup,                  &             
                         sum(dusrfdew(ew-1:ew,ns-1:ns)), &
                         sum(dthckdew(ew-1:ew,ns-1:ns)))

            ugradns = horizderiv(upn,  stagsigma,        &
                         sum(uvel(:,ew-1:ew,ns-1:ns),2), &
                         dns4, ugradup,                  &                              
                         sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                         sum(dthckdns(ew-1:ew,ns-1:ns)))

            vgradns = horizderiv(upn,  stagsigma,        &
                         sum(vvel(:,ew-1:ew,ns-1:ns),2), &
                         dns4, vgradup,                  &                              
                         sum(dusrfdns(ew-1:ew,ns-1:ns)), &
                         sum(dthckdns(ew-1:ew,ns-1:ns)))

            ! "effstr" = eff. strain rate squared
            effstr = ugradew**2 + vgradns**2 + ugradew*vgradns + &
                         0.25_dp * (vgradew + ugradns)**2 + &
!                         f1 * (ugradup**2 + vgradup**2)              ! make line ACTIVE for "capping" version (see note below)   
                         f1 * (ugradup**2 + vgradup**2) + effstrminsq ! make line ACTIVE for new version

    ! -----------------------------------------------------------------------------------
    ! NOTES on capping vs. non-capping version of eff. strain rate calc.
    ! -----------------------------------------------------------------------------------
    !
    ! Set eff. strain rate (squared) to some min value where it falls below some 
    ! threshold value, 'effstrminsq'. Commented out the old version below, which "caps" 
    ! the min eff strain rate (and thus the max eff visc) in favor of a version that 
    ! leads to a "smooth" description of eff strain rate (and eff visc). The change for 
    ! new version is that the value of 'effstrminsq' simply gets added in with the others
    ! (e.g. how it is done in the Pattyn model). The issues w/ the capping approach are 
    ! discussed (w.r.t. sea ice model) in: Lemieux and Tremblay, JGR, VOL. 114, C05009, 
    ! doi:10.1029/2008JC005017, 2009). Long term, the capping version should probably be
    ! available as a config file option or possibly removed altogether.   

    ! Old "capping" scheme       ! these lines must be active to use the "capping" scheme for the efvs calc
!            where (effstr < effstrminsq)
!                   effstr = effstrminsq
!            end where

    ! Note that the vert dims are explicit here, since glide_types defines this 
    ! field as having dims 1:upn. This is something that we'll have to decide on long-term;
    ! should efvs live at cell centroids in the vert (as is assumed in this code)
    ! or should we be doing some one-sided diffs at the sfc/bed boundaries so that it has vert dims 
    ! of upn? For now, we populate ONLY the first 1:upn-1 values of the efvs vector and leave the one
    ! at upn empty (the Pattyn/Bocek/Johnson core would fill all values, 1:upn).

    ! NOTE also that efvs lives on the non-staggered grid in the horizontal. That is, in all of the 
    ! discretizations conducted below, efvs is explicitly averaged from the normal horiz grid onto the 
    ! staggered horiz grid (Thus, in the calculations, efvs is treated as if it lived on the staggered 
    ! horiz grid, even though it does not). 

            ! Below, p2=(1-n)/2n. The 1/2 is from taking the sqr root of the squared eff. strain rate
            efvs(1:upn-1,ew,ns) = flwafact(1:upn-1,ew,ns) * effstr**p2

        else  
           efvs(:,ew,ns) = effstrminsq ! if the point is associated w/ no ice, set to min value
        end if

       end do   ! end ew
   end do       ! end ns

  case(1)       ! set the eff visc to some const value 

    efvs = 1.0_dp

  end select

  return

end subroutine findefvsstr

!***********************************************************************

function vertideriv(upn, varb, thck)

  implicit none 

  integer, intent(in) :: upn
  real (kind = dp), intent(in), dimension(:) :: varb
  real (kind = dp), intent(in) :: thck    

  real (kind = dp), dimension(size(varb)-1) :: vertideriv
  !'dupm' is defined as -1/(2*del_sigma), in which case it seems like 
  ! there should be a '-' in front of this expression ... but note that
  ! the negative sign is implicit in the fact that the vertical index 
  ! increases moving downward in the ice column (up=1 is the sfc, 
  ! up=upn is the bed).

  vertideriv(1:upn-1) = dupm * (varb(2:upn) - varb(1:upn-1)) / thck

  return

end function vertideriv

!***********************************************************************

function horizderiv(upn,     stagsigma,   &
                    varb,    grid,        &
                    dvarbdz, dusrfdx, dthckdx)

  implicit none
  
  integer, intent(in) :: upn
  real (kind = dp), dimension(:), intent(in) :: stagsigma
  real (kind = dp), dimension(:,:), intent(in) :: varb
  real (kind = dp), dimension(:), intent(in) :: dvarbdz
  real (kind = dp), intent(in) :: dusrfdx, dthckdx, grid

  real (kind = dp) :: horizderiv(size(varb,1)-1)
  
  horizderiv = (varb(1:upn-1,2) + varb(2:upn,2) - varb(1:upn-1,1) - varb(2:upn,1)) / grid - &    
                dvarbdz * (dusrfdx - stagsigma * dthckdx) / 4.0_dp

  return

end function horizderiv

!***********************************************************************

function getlocrange(upn, indx)

  implicit none

  integer, intent(in) :: upn 
  integer, intent(in) :: indx
  integer, dimension(2) :: getlocrange

  getlocrange = (indx - 1) * (upn + 2) + 1 + (/ 1, upn /)

  return

end function getlocrange

!***********************************************************************

function getlocationarray(ewn, nsn, upn,  &
                          mask )
    implicit none

    integer, intent(in) :: ewn, nsn, upn 
    integer, dimension(:,:), intent(in) :: mask
    integer, dimension(ewn-1,nsn-1) :: getlocationarray, temparray
    integer :: cumsum, ew, ns

    cumsum = 0

    do ew=1,ewn-1
        do ns=1,nsn-1
        if ( GLIDE_HAS_ICE( mask(ew,ns) ) ) then
            cumsum = cumsum + ( upn + 2 )
            getlocationarray(ew,ns) = cumsum
            temparray(ew,ns) = upn + 2
        else
            getlocationarray(ew,ns) = 0
            temparray(ew,ns) = 1
        end if
        end do
    end do

    getlocationarray = ( getlocationarray + 1 ) - temparray

    return

end function getlocationarray

!***********************************************************************

subroutine solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, vel )

  ! Puts sparse matrix variables in SLAP triad format into "matrix" derived type, 
  ! so that it can be passed to the generic solver wrapper, "sparse_easy_solve". 
  ! Takes place of the old, explicit solver interface to SLAP linear solver.

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  real (kind = dp), dimension(:,:,:), intent(in) :: vel
  integer, dimension(:,:), intent(in) :: uindx
  type(sparse_matrix_type), intent(inout) :: matrix
  real (kind = dp), dimension(:), intent(out) :: answer

  integer :: ew, ns
  integer, dimension(2) :: loc

  pcgsize(2) = ct - 1

  matrix%order = pcgsize(1) 
  matrix%nonzeros = pcgsize(2)
  matrix%symmetric = .false.

  matrix%row = pcgrow
  matrix%col = pcgcol
  matrix%val = pcgval

  ! Initial estimate for vel. field; take from 3d array and put into
  ! the format of a solution vector.
  do ns = 1,nsn-1
    do ew = 1,ewn-1
        if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            answer(loc(1):loc(2)) = vel(:,ew,ns)
            answer(loc(1)-1) = vel(1,ew,ns)
            answer(loc(2)+1) = vel(upn,ew,ns)
        end if
    end do
  end do

end subroutine solver_preprocess

!***********************************************************************

subroutine solver_postprocess( ewn, nsn, upn, uindx, answrapped, ansunwrapped )

  ! Unwrap the vels from the solution vector and place into a 3d array.

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: uindx
  real (kind = dp), dimension(:), intent(in) :: answrapped
  real (kind = dp), dimension(upn,ewn-1,nsn-1), intent(out) :: ansunwrapped 

  integer, dimension(2) :: loc
  integer :: ew, ns

  do ns = 1,nsn-1
      do ew = 1,ewn-1
          if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            ansunwrapped(:,ew,ns) = answrapped(loc(1):loc(2))
          else 
            ansunwrapped(:,ew,ns) = 0.0d0
          end if
      end do
  end do

end subroutine solver_postprocess

!***********************************************************************

function mindcrshstr(pt,whichresid,vel,counter,resid)

  ! Function to perform 'unstable manifold correction' (see Hindmarsch and Payne, 1996,
  ! "Time-step limits for stable solutions of the ice-sheet equation", Annals of 
  ! Glaciology, 23, p.74-85)

  implicit none

  real (kind = dp), intent(in), dimension(:,:,:) :: vel
  integer, intent(in) :: counter, pt, whichresid 

  real (kind = dp), intent(out) :: resid

  real (kind = dp), dimension(size(vel,1),size(vel,2),size(vel,3)) :: mindcrshstr

  real (kind = dp), parameter :: ssthres = 5.0_dp * pi / 6.0_dp, &
                                 critlimit = 10.0_dp / (scyr * vel0), &
                                 small = 1.0e-16_dp

  real (kind = dp), intrinsic :: abs, acos
  
  integer, dimension(2), save :: new = 1, old = 2
  integer :: locat(3)

  integer :: nr
  integer, dimension(size(vel,1),size(vel,2),size(vel,3)) :: vel_ne_0

  if (counter == 1) then
    usav(:,:,:,pt) = 0.0d0
  end if

  corr(:,:,:,new(pt),pt) = vel - usav(:,:,:,pt)           

  if (counter > 1) then

    where (acos((corr(:,:,:,new(pt),pt) * corr(:,:,:,old(pt),pt)) / &
          (abs(corr(:,:,:,new(pt),pt)) * abs(corr(:,:,:,old(pt),pt)) + small)) > &
           ssthres .and. corr(:,:,:,new(pt),pt) - corr(:,:,:,old(pt),pt) /= 0.0_dp )

      mindcrshstr = usav(:,:,:,pt) + &
                    corr(:,:,:,new(pt),pt) * abs(corr(:,:,:,old(pt),pt)) / &
                    abs(corr(:,:,:,new(pt),pt) - corr(:,:,:,old(pt),pt)) 

    elsewhere

      mindcrshstr = vel;
 
    end where
    
  else 

    mindcrshstr = vel;
   
  end if

  if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 1; new(pt) = 2; end if

  select case (whichresid)

  ! options for residual calculation method, as specified in configuration file 
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel ) 
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels 
  ! case(2): use mean of abs( vel_old - vel ) / vel )

   case(0)
    resid = maxval( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel .ne. 0.0_dp)  
    locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel .ne. 0.0_dp)

   case(1)
    nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
                        MASK = vel .ne. 0.0_dp)
    locat = maxloc( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
            MASK = vel .ne. 0.0_dp)

   case(2)
    nr = size( vel, dim=1 )
    vel_ne_0 = 0
    where ( vel .ne. 0.0_dp ) vel_ne_0 = 1

    ! include basal velocities in resid. calculation when using MEAN
    resid = sum( abs((usav(:,:,:,pt) - vel ) / vel ), &
            MASK = vel .ne. 0.0_dp) / sum( vel_ne_0 )

    ! ignore basal velocities in resid. calculation when using MEAN
    ! resid = sum( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),   &
    !           MASK = vel .ne. 0.0_dp) / sum( vel_ne_0(1:nr-1,:,:) )

    ! NOTE that the location of the max residual is somewhat irrelevent here
    !      since we are using the mean resid for convergence testing
    locat = maxloc( abs((usav(:,:,:,pt) - vel ) / vel ), MASK = vel .ne. 0.0_dp)

  end select

    usav(:,:,:,pt) = vel

    ! Additional debugging line, useful when trying to determine if convergence is being consistently 
    ! help up by the residual at one or a few particular locations in the domain.
!    print '("* ",i3,g20.6,3i6,g20.6)', counter, resid, locat, vel(locat(1),locat(2),locat(3))*vel0

  return

end function mindcrshstr

!***********************************************************************

subroutine findcoefstr(ewn,  nsn,   upn,            &
                       dew,  dns,   sigma,          &
                       pt,          efvs,           &
                       thisvel,     othervel,       &
                       thck,        thisdusrfdx,    &
                       dusrfdew,    dthckdew,       &
                       d2usrfdew2,  d2thckdew2,     &
                       dusrfdns,    dthckdns,       &
                       d2usrfdns2,  d2thckdns2,     &
                       d2usrfdewdns,d2thckdewdns,   &
                       dlsrfdew,    dlsrfdns,       &
                       stagthck,    whichbabc,      &
                       uindx,       mask,           &
                       lsrf,        topg,           &
                       minTauf,     flwa,           &    
                       beta, count )

  ! Main subroutine for determining coefficients that go into the LHS matrix A 
  ! in the expression Au = b. Calls numerous other subroutines, including boundary
  ! condition subroutines, which determin "b".

  implicit none

  integer, intent(in) :: ewn, nsn, upn, count 
  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), dimension(:), intent(in) :: sigma

  real (kind = dp), dimension(:,:,:), intent(in) :: efvs, thisvel, &
                                                    othervel
  real (kind = dp), dimension(:,:), intent(in) :: stagthck, thisdusrfdx,     &
                                                  dusrfdew,   dthckdew,      &
                                                  d2usrfdew2, d2thckdew2,    &
                                                  dusrfdns,   dthckdns,      &
                                                  d2usrfdns2, d2thckdns2,    &
                                                  d2usrfdewdns,d2thckdewdns, &
                                                  dlsrfdew,   dlsrfdns,      &
                                                  thck, lsrf, topg 

  real (kind = dp), dimension(:,:), intent(in) :: minTauf

  real (kind = dp), dimension(:,:), intent(in) :: beta

  real (kind = dp), dimension(:,:,:), intent(in) :: flwa

  integer, dimension(:,:), intent(in) :: mask, uindx
  integer, intent(in) :: pt, whichbabc

  real (kind = dp), dimension(ewn-1,nsn-1) :: betasquared
  real (kind = dp), dimension(2,2,2) :: localefvs   
  real (kind = dp), dimension(3,3,3) :: localothervel
  real (kind = dp), dimension(upn) :: boundaryvel
  real (kind = dp) :: flwabar

  integer, dimension(ewn-1,nsn-1) :: loc_array
  integer, dimension(6) :: loc
  integer, dimension(3) :: shift
  integer :: ew, ns, up

  ct = 1        ! index to count the number of non-zero entries in the sparse matrix

  ! calculate/specify the map of 'betasquared', for use in the basal boundary condition. 
  ! Input to the subroutine 'bodyset' (below) ... 
  call calcbetasquared (whichbabc,              &
                        dew,        dns,        &
                        ewn,        nsn,        &
                        lsrf,       topg,       &
                        thck,                   &
                        thisvel(upn,:,:),       &   
                        othervel(upn,:,:),      &
                        minTauf, beta,          &
                        betasquared )

  do ns = 1,nsn-1
    do ew = 1,ewn-1 

     ! Calculate the depth-averaged value of the rate factor, needed below when applying an ice shelf
     ! boundary condition (complicated code so as not to include funny values at boundaries ...
     ! ... kind of a mess and could be redone or made into a function or subroutine).
     flwabar = ( sum( flwa(:,ew,ns), 1, flwa(1,ew,ns)*vis0_glam < 1.0d-10 )/real(upn) + &
               sum( flwa(:,ew,ns+1), 1, flwa(1,ew,ns+1)*vis0_glam < 1.0d-10 )/real(upn)  + &
               sum( flwa(:,ew+1,ns), 1, flwa(1,ew+1,ns)*vis0_glam < 1.0d-10 )/real(upn)  + &
               sum( flwa(:,ew+1,ns+1), 1, flwa(1,ew+1,ns+1)*vis0_glam < 1.0d-10 )/real(upn) ) / &
               ( sum( flwa(:,ew,ns)/flwa(:,ew,ns), 1, flwa(1,ew,ns)*vis0_glam < 1.0d-10 )/real(upn) + &
               sum( flwa(:,ew,ns+1)/flwa(:,ew,ns+1), 1, flwa(1,ew,ns+1)*vis0_glam < 1.0d-10 )/real(upn) + &
               sum( flwa(:,ew+1,ns)/flwa(:,ew+1,ns), 1, flwa(1,ew+1,ns)*vis0 < 1.0d-10 )/real(upn) + &
               sum( flwa(:,ew+1,ns+1)/flwa(:,ew+1,ns+1), 1, flwa(1,ew+1,ns+1)*vis0_glam < 1.0d-10 )/real(upn) )


    if( ns == 1 .and. ew == 1 ) then
           loc_array = getlocationarray(ewn, nsn, upn, mask )
    end if

!  !!!!!!!!! useful for debugging !!!!!!!!!!!!!!
!    print *, 'loc_array = '
!    print *, loc_array
!    pause

    loc(1) = loc_array(ew,ns)

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if ( GLIDE_HAS_ICE(mask(ew,ns)) .and. .not. &
         GLIDE_IS_COMP_DOMAIN_BND(mask(ew,ns)) .and. .not. &
         GLIDE_IS_MARGIN(mask(ew,ns)) .and. .not. &
         GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .and. .not. &
         GLIDE_IS_CALVING(mask(ew,ns) ) ) &
    then
!    print *, 'In main body ... ew, ns = ', ew, ns
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        call calccoeffs( upn,               sigma,              &
                         stagthck(ew,ns),                       &
                         dusrfdew(ew,ns),   dusrfdns(ew,ns),    &
                         dthckdew(ew,ns),   dthckdns(ew,ns),    &
                         d2usrfdew2(ew,ns), d2usrfdns2(ew,ns),  &
                         d2usrfdewdns(ew,ns),                   &
                         d2thckdew2(ew,ns), d2thckdns2(ew,ns),  &
                         d2thckdewdns(ew,ns))

        ! get index of cardinal neighbours
        loc(2) = loc_array(ew+1,ns)
        loc(3) = loc_array(ew-1,ns)
        loc(4) = loc_array(ew,ns+1)
        loc(5) = loc_array(ew,ns-1)

        ! this loop fills coeff. for all vertical layers at index ew,ns (including sfc. and bed bcs)
        do up = 1,upn

            ! Function to adjust indices at sfc and bed so that most correct values of 'efvs' and 'othervel'
            ! are passed to function. Because of the fact that efvs goes from 1:upn-1 rather than 1:upn
            ! we simply use the closest values. This could probably be improved upon at some point
            ! by extrpolating values for efvs at the sfc and bed using one-sided diffs, and it is not clear
            ! how important this simplfication is.
            shift = indshift( 0, ew, ns, up, ewn, nsn, upn, loc_array, stagthck(ew-1:ew+1,ns-1:ns+1) )  

            call bodyset(ew,  ns,  up,        &
                         ewn, nsn, upn,       &
                         dew,      dns,       &
                         pt,       loc_array, &
                         loc,      stagthck,  &
                         thisdusrfdx,         &
                         dusrfdew, dusrfdns,  &
                         dlsrfdew, dlsrfdns,  &
                         efvs(up-1+shift(1):up+shift(1),ew:ew+1,ns:ns+1),  &
                         othervel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         betasquared(ew,ns) ) 

        end do  ! upn

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    elseif ( GLIDE_IS_CALVING( mask(ew,ns) ) .and. .not. &
             GLIDE_IS_COMP_DOMAIN_BND(mask(ew,ns) ) .and. .not. &
             GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) ) &
    then
!    print *, 'At a SHELF boundary ... ew, ns = ', ew, ns
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        call calccoeffs( upn,               sigma,              &
                         stagthck(ew,ns),                       &
                         dusrfdew(ew,ns),   dusrfdns(ew,ns),    &
                         dthckdew(ew,ns),   dthckdns(ew,ns),    &
                         d2usrfdew2(ew,ns), d2usrfdns2(ew,ns),  &
                         d2usrfdewdns(ew,ns),                   &
                         d2thckdew2(ew,ns), d2thckdns2(ew,ns),  &
                         d2thckdewdns(ew,ns))

        do up = 1, upn
           lateralboundry = .true.
           shift = indshift(  1, ew, ns, up,                   &
                               ewn, nsn, upn,                  &
                               loc_array,                      & 
                               stagthck(ew-1:ew+1,ns-1:ns+1) )

            call bodyset(ew,  ns,  up,        &
                         ewn, nsn, upn,       &
                         dew,      dns,       &
                         pt,       loc_array, &
                         loc,      stagthck,  &
                         thisdusrfdx,         &
                         dusrfdew, dusrfdns,  &
                         dlsrfdew, dlsrfdns,  &
                         efvs(up-1+shift(1):up+shift(1),ew:ew+1,ns:ns+1),  &
                         othervel(up-1+shift(1):up+1+shift(1),  &
                         ew-1+shift(2):ew+1+shift(2),  &
                         ns-1+shift(3):ns+1+shift(3)), &
                         betasquared(ew,ns), abar=flwabar, cc=count )        
        end do
        lateralboundry = .false.

! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    elseif ( GLIDE_HAS_ICE(mask(ew,ns)) .and. ( GLIDE_IS_DIRICHLET_BOUNDARY(mask(ew,ns)) .or. &
             GLIDE_IS_COMP_DOMAIN_BND(mask(ew,ns)) ) .or. GLIDE_IS_LAND_MARGIN(mask(ew,ns))) &
    then
!    print *, 'At a NON-SHELF boundary ... ew, ns = ', ew, ns
! >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        ! Put specified value for vel on rhs. NOTE that this is NOT zero by default 
        ! unless the initial guess is zero. It will be set to whatever the initial value 
        ! for the vel at location up,ew,ns is in the initial array!
        locplusup = loc(1)
        call valueset(0.0_dp)
        locplusup = loc(1) + upn + 1
        call valueset(0.0_dp)
        do up = 1,upn
           locplusup = loc(1) + up
           call valueset( thisvel(up,ew,ns) )     ! vel at margin set to initial value 
!           call valueset( 0.0_dp )                ! vel at margin set to 0 
        end do

    end if

    end do;     ! ew 
  end do        ! ns

  return

end subroutine findcoefstr

!***********************************************************************

subroutine bodyset(ew,  ns,  up,           &
                   ewn, nsn, upn,          &
                   dew,      dns,          &
                   pt,       loc_array,    &
                   loc,      stagthck,     &
                   thisdusrfdx,            &
                   dusrfdew, dusrfdns,     &
                   dlsrfdew, dlsrfdns,     &
                   local_efvs,             &
                   local_othervel,         &
                   betasquared,            &
                   local_thisvel,          &
                   abar, cc)

  ! This subroutine does the bulk of the work in calling the appropriate discretiztion routines,
  ! which determine the values for coefficients that will go into the sparse matrix, for points
  ! on and inside of the boundaries.

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in) :: ew, ns, up
  real (kind = dp), intent(in) :: dew, dns
  integer, intent(in) :: pt
  integer, dimension(ewn-1,nsn-1), intent(in) :: loc_array
  integer, dimension(6), intent(in) :: loc

  real (kind = dp), dimension(:,:), intent(in) :: stagthck
  real (kind = dp), dimension(:,:), intent(in) :: dusrfdew, dusrfdns
  real (kind = dp), dimension(:,:), intent(in) :: dlsrfdew, dlsrfdns
  real (kind = dp), dimension(:,:), intent(in) :: thisdusrfdx

  real (kind = dp), dimension(2,2,2), intent(in) :: local_efvs

  ! "local_othervel" is the other vel component (i.e. u when v is being calc and vice versa),
  ! which is taken as a known value (terms involving it are moved to the RHS and treated as sources)
  real (kind = dp), dimension(3,3,3), intent(in) :: local_othervel

  real (kind = dp), intent(in) :: betasquared
  real (kind = dp), intent(in), optional :: local_thisvel
  real (kind = dp), intent(in), optional :: abar
  integer, intent(in), optional :: cc

  ! storage space for coefficients that go w/ the discretization at the local point up, ew, ns
  real (kind = dp), dimension(3,3,3) :: g

  ! source term for the rhs when using ice shelf lateral boundary condition,
  ! e.g. source = rho*g*H/(2*Neff) * ( 1 - rho_i / rho_w ) for ice shelf
  real (kind = dp) :: source

  real (kind = dp) :: slopex, slopey    ! local sfc (or bed) slope terms

  ! lateral boundary normal and vector to indicate use of forward
  ! or bacward one-sided diff. when including specified stress lateral bcs
  real (kind = dp), dimension(2) :: fwdorbwd, normal

  real (kind = dp) :: nz   ! z dir normal vector component at sfc or bed (takes diff value for each)

  integer, dimension(2) :: bcflag  ! indicates choice of sfc and basal bcs ...

  real (kind = dp) :: efvstot   ! both of these vars are used for averaging of eff. vis. near lat boundaries
  integer :: efvscount, i, j, k
  efvstot = 0.0d0; efvscount = 0; i = 0; j = 0; k = 0

  locplusup = loc(1) + up

  if( lateralboundry )then

  ! *********************************************************************************************
  ! lateral boundary conditions 
  
  ! if at sfc or bed, source due to seawater pressure is 0 and bc normal vector
  ! should contain sfc/bed slope components, e.g. (-ds/dx, -ds/dy, 1) or (db/dx, db/dy, -1)
     source = 0.0_dp

     call getlatboundinfo( ew,  ns,  up,                            &
                           ewn, nsn, upn,                           & 
                           stagthck(ew-1:ew+1, ns-1:ns+1),          &
                           loc_array, fwdorbwd, normal, loc_latbc)
                
     if( up == 1 .or. up == upn )then
    
        if( up == 1 )then
           locplusup = loc(1) + up - 1  ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.0_dp
        else
           locplusup = loc(1) + up + 1  ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.0_dp
        end if

        g = normhorizmainbc_lat(dew,           dns,             &
                                slopex,        slopey,          &
                                dsigmadew(up), dsigmadns(up),   &
                                pt,            2,               &
                                dup(up),                        &
                                oneorfour,     fourorone,       &
                                onesideddiff,                   &
                                normal,        fwdorbwd)

        ! add on coeff. associated w/ du/digma  
        g(:,3,3) = g(:,3,3) & 
                 + vertimainbc( stagthck(ew,ns), bcflag,dup(up),     &
                                local_efvs,      betasquared,    nz )

        ! put the coeff. for the b.c. equation in the same place as the prev. equation
        ! (w.r.t. cols), on a new row ...
        call fillsprsebndy( g, locplusup, loc_latbc, up, normal ) 

        ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
        ! which results from moving them from the LHS over to the RHS, has been moved
        ! inside of "croshorizmainbc_lat".
        rhsd(locplusup) = sum( croshorizmainbc_lat(dew,           dns,           &
                                                   slopex,        slopey,        &
                                                   dsigmadew(up), dsigmadns(up), &
                                                   pt,            2,             &
                                                   dup(up),       local_othervel,&
                                                   oneortwo,      twoorone,      &
                                                   onesideddiff,                 &
                                                   normal,fwdorbwd)              &         
                                                 * local_othervel )
                
    end if     ! up = 1 or up = upn (IF at lateral boundary and IF at surface or bed)

    ! If in main body and at ice/ocean boundary, calculate depth-averaged stress
    ! due to sea water, bc normal vector components should be boundary normal 
    locplusup = loc(1) + up

    ! for this bc, the normal vector components are not the sfc/bed slopes but are taken
    ! from a normal to the shelf front in map view (x,y plane); slopex,slopey are simply renamed here 
    slopex = normal(1)          
    slopey = normal(2)

    ! Two options here, (1) use the 1d solution that involves the rate factor (not accurate for 
    !                       3d domains, but generally more stable 
    !                   (2) use the more general solution that involves the eff. visc. and normal
    !                       vector orientation at lateral boundary
    ! ... only one of these options should be active at a time (comment the other lines out)
    ! The default setting is (2), the more general case that should also work in the 1d case.

    ! In some cases, the two options can be used together to improve performance, e.g. for the Ross
    ! ice shelf experiment, a number of early iterations use the more simple bc (option 1) and then
    ! when the solution has converged a bit, we switch to the more realistic implementation (option 2).
    ! That is achieved in the following if construct ...

!    if( cc < 10 )then   ! use this to "pre-condition" the shelf BC w/ the simple, 1d version
!    if( cc >= 0 )then   ! use this to use only the 1d version
    if( cc > 1000000 )then   ! use this to go straight to the full 2d version of the bc

    ! --------------------------------------------------------------------------------------
    ! (1) source term (strain rate at shelf/ocean boundary) from Weertman's analytical solution 
    ! --------------------------------------------------------------------------------------
    ! See eq. 2, Pattyn+, 2006, JGR v.111; eq. 8, Vieli&Payne, 2005, JGR v.110). Note that this 
    ! contains the 1d assumption that ice is not spreading lateraly !(assumes dv/dy = 0 for u along flow)
    ! Note that factor of 2 in front of 'stagthck' is NOT part of the standard bc. Here, it is used to 
    ! correct for the fact that the staggered thickness will be 1/2 of the normal thickness at a boundary 
    ! ... as of summer 2009, this hack has been removed (no more factor of 2) and replaced by a new stagthck
    ! averaging scheme at the margins, which uses only the non-zero values of thickness on the normal grid to
    ! calc. the value of the stag. thickness
    source = abar*vis0_glam * ( 1.0_dp/4.0_dp * rhoi * grav * stagthck(ew,ns)*thk0 * ( 1.0_dp - rhoi/rhoo))**3.0_dp

    ! multiply by 4 so that case where v=0, du/dy = 0, LHS gives: du/dx = du/dx|_shelf 
    ! (i.e. LHS = 4*du/dx, requires 4*du/dx_shelf)
    source = source * 4.0_dp

    ! split source based on the boundary normal orientation and non-dimensinoalize
    ! Note that it is not really appropriate to apply option (1) to 1d flow, since terms other than du/dx in 
    ! eff. strain rate are ignored. For 2d flow, should use option (2) below. 
     source = source * normal(pt)
     source = source * tim0 ! make source term non-dim
    ! --------------------------------------------------------------------------------------

  else

    ! --------------------------------------------------------------------------------------
    ! (2) source term (strain rate at shelf/ocean boundary) from MacAyeal depth-ave solution. 
    ! --------------------------------------------------------------------------------------
    source = (rhoi*grav*stagthck(ew,ns)*thk0) / tau0_glam / 2.0_dp * ( 1.0_dp - rhoi / rhoo )

    ! terms after "/" below count number of non-zero efvs cells ... needed for averaging of the efvs at boundary 
    source = source / ( sum(local_efvs, local_efvs > 1.0d-12) / &
             sum( (local_efvs/local_efvs), local_efvs > 1.0d-12 ) )

    source = source * normal(pt) ! partition according to normal vector at lateral boundary
                                 ! NOTE that source term is already non-dim here 
    ! --------------------------------------------------------------------------------------

  end if

                        
    g = normhorizmainbc_lat(dew,           dns,        &
                            slopex,        slopey,     &
                            dsigmadew(up), dsigmadns(up),  &
                            pt,            1,          &
                            dup(up),                   &
                            oneorfour,     fourorone,  &
                            onesideddiff,              &
                            normal,        fwdorbwd)

    ! put the coeff. for the b.c. equation in the same place as the prev. equation
    ! (w.r.t. cols), on a new row ...
    call fillsprsebndy( g, locplusup, loc_latbc, up, normal )

    ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
    ! which results from moving them from the LHS over to the RHS, has been moved
    ! inside of "croshorizmainbc_lat".
    rhsd(locplusup) = sum( croshorizmainbc_lat(dew,           dns,            &
                                               slopex,        slopey,         &
                                               dsigmadew(up), dsigmadns(up),  &
                                               pt,            1,              &
                                               dup(up),       local_othervel, &
                                               oneortwo,      twoorone,       &
                                               onesideddiff,                  &
                                               normal,        fwdorbwd)       &
                                              * local_othervel ) + source
 
  else   ! NOT at a lateral boundary 

! *********************************************************************************************
! normal discretization for points inside of lateral boundary and inside main body of ice sheet
        
     g = normhorizmain(pt,up,local_efvs)
     g(:,2,2) = g(:,2,2) + vertimain(hsum(local_efvs),up)             
     call fillsprsemain(g,locplusup,loc,up) 
     ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
     ! which results from moving them from the LHS over to the RHS, is explicit and 
     ! hast NOT been moved inside of "croshorizmin" (as is the case for the analogous
     ! boundary condition routines).
     rhsd(locplusup) = thisdusrfdx(ew,ns) - sum(croshorizmain(pt,up,local_efvs) * local_othervel)

  end if

! *********************************************************************************************
! higher-order sfc and bed boundary conditions in main body of ice sheet (NOT at lat. boundry)

  if(  ( up == upn  .or. up == 1 ) .and. .not. lateralboundry) then


     if( up == 1 )then                ! specify necessary variables and flags for free sfc
        bcflag = (/1,0/)
        locplusup = loc(1) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
        slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.0_dp
     else                             ! specify necessary variables and flags for basal bc

        !bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...

                                      ! better to specify very large value for betasquared below
        bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                      ! where betasquared is MacAyeal-type traction parameter

        locplusup = loc(1) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
        slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.0_dp
     end if

     g = normhorizmainbc(dew,           dns,     &
                         slopex,        slopey,  &
                         dsigmadew(up), dsigmadns(up),  &
                         pt,            bcflag,  &
                         dup(up),                &
                         oneorfour,     fourorone)

     ! add on coeff. associated w/ du/digma
     g(:,2,2) = g(:,2,2)   &
              + vertimainbc(stagthck(ew,ns), bcflag,dup(up), local_efvs, betasquared, nz)

     ! put the coeff. for the b.c. equation in the same place as the prev. equation
     ! (w.r.t. cols), on a new row ...
     call fillsprsemain(g,locplusup,loc,up)

     ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
     ! which results from moving them from the LHS over to the RHS, has been moved
     ! inside of "croshorizmainbc".
     rhsd(locplusup) =  sum( croshorizmainbc(dew,           dns,            &
                                            slopex,        slopey,         &
                                            dsigmadew(up), dsigmadns(up),  &
                                            pt,            bcflag,         &
                                            dup(up),       local_othervel, &
                                            oneortwo,      twoorone)       &
                                            * local_othervel )

  end if   ! (up = 1 or up = upn) and lateralboundry = F

! *********************************************************************************************

  return

end subroutine bodyset

!***********************************************************************

subroutine valueset(local_value)

  ! plugs given value into the right location in the rhs vector of matrix equation Ax=rhs

  implicit none

  real (kind = dp), intent(in) :: local_value

  call putpcgc(1.0_dp,locplusup,locplusup)
  rhsd(locplusup) = local_value 

  return

end subroutine valueset

!***********************************************************************

subroutine calccoeffsinit (upn, dew, dns)

  ! determines constants used in various FD calculations associated with 'findcoefst'
  ! In general, the constants contain (1) grid spacing info, (2) numeric constants 
  ! used for averaging of eff. visc. from normal grid in horiz onto stag grid in horiz. 
  implicit none

  integer, intent(in) :: upn
  real (kind = dp), intent(in) :: dew, dns

  ! this coefficient used in finite differences of vertical terms.
  cvert(:) = (len0**2) / (4.0_dp * thk0**2 * dup**2)

  ! these coefficients used in finite differences of horizontal terms
  ! for d/dx(fdu/dx), d/dx(fdu/dy), d/dsigma(fdu/dx), d/dx(fdu/dsigma) and
  ! du/dsigma. 
  cdxdx = (/ 0.25_dp / dew**2, 0.25_dp / dns**2 /)
  cdsdx(:,1) = 0.0625_dp / (dew * dup); cdsdx(:,2) = 0.0625_dp / (dns * dup);
  cdsds = 0.25_dp / (dup * dup)
  cds = 0.0625_dp / dup
  cdxdy = 0.0625_dp / (dew * dns)

  return
     
end subroutine calccoeffsinit

!***********************************************************************

subroutine calccoeffs(upn,        sigma,                    &
                      stagthck,                             &
                      dusrfdew,   dusrfdns,                 &
                      dthckdew,   dthckdns,                 &
                      d2usrfdew2, d2usrfdns2, d2usrfdewdns, &
                      d2thckdew2, d2thckdns2, d2thckdewdns)

  ! Called from 'findcoefst' to find coefficients in stress balance equations
  ! Detemines coeficients needed for finite differencing.
  ! This is a column-based operation. In general these coefficients refer 
  ! to grid transformations and averaging of efvs to half grid points.

  implicit none

  integer, intent(in) :: upn 
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), intent(in) :: stagthck, dusrfdew, dusrfdns, dthckdew, dthckdns, &
                                  d2usrfdew2, d2usrfdns2, d2usrfdewdns, &
                                  d2thckdew2, d2thckdns2, d2thckdewdns
     
  fvert(:) = cvert(:) / stagthck**2     

  dsigmadew = calcdsigmadx(upn, sigma, dusrfdew, dthckdew, stagthck)
  dsigmadns = calcdsigmadx(upn, sigma, dusrfdns, dthckdns, stagthck)

  d2sigmadew2 = calcd2sigmadxdy(upn,        sigma,      &
                                d2usrfdew2, d2thckdew2, &
                                dusrfdew,   dusrfdew,   &
                                dthckdew,   dthckdew,   &
                                stagthck)

  d2sigmadns2 = calcd2sigmadxdy(upn,        sigma,      &
                                d2usrfdns2, d2thckdns2, &
                                dusrfdns,   dusrfdns,   &
                                dthckdns,   dthckdns,   &
                                stagthck)

  d2sigmadewdns = calcd2sigmadxdy(upn,          sigma,         &
                                  d2usrfdewdns, d2thckdewdns,  &
                                  dusrfdew,     dusrfdns,      &
                                  dthckdew,     dthckdns,      &
                                  stagthck)

  d2sigmadewdsigma = calcd2sigmadxdsigma(dthckdew,stagthck)
  d2sigmadnsdsigma = calcd2sigmadxdsigma(dthckdns,stagthck)

  return

end subroutine calccoeffs

!***********************************************************************

function calcdsigmadx(upn,     sigma,    &
                      dusrfdx, dthckdx,  &
                      stagthck)

  implicit none

  integer, intent(in) :: upn  
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), intent(in) :: stagthck, dusrfdx, dthckdx
  real (kind = dp), dimension(upn) :: calcdsigmadx

  calcdsigmadx = (dusrfdx - sigma * dthckdx) / stagthck

  return

end function calcdsigmadx

!***********************************************************************

function calcd2sigmadxdy(upn,        sigma,       &
                         d2usrfdxdy, d2thckdxdy,  &
                         dusrfdx,    dusrfdy,     &
                         dthckdx,    dthckdy,     &
                         stagthck)

  implicit none

  integer, intent(in) :: upn 
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), intent(in) :: d2usrfdxdy, d2thckdxdy, dusrfdx, dusrfdy, &
                                  dthckdx, dthckdy, stagthck 
  real (kind = dp), dimension(upn) :: calcd2sigmadxdy

  calcd2sigmadxdy = (stagthck * d2usrfdxdy - &
                     dusrfdx * dthckdy - dusrfdy * dthckdx + &
                     sigma * (2.0_dp * dthckdx * dthckdy - &
                     stagthck * d2thckdxdy)) / stagthck**2

  return

end function calcd2sigmadxdy

!***********************************************************************

function calcd2sigmadxdsigma(dthckdx,stagthck)

  implicit none

  real (kind = dp), intent(in) :: dthckdx, stagthck 
  real (kind = dp) :: calcd2sigmadxdsigma

  calcd2sigmadxdsigma = - dthckdx / stagthck

  return

end function calcd2sigmadxdsigma

!***********************************************************************

function vertimain(efvs,up)

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs

  real (kind = dp), dimension(3) :: vertimain

  integer, intent(in) :: up

  vertimain(3) = fvert(up) * efvs(2) 
  vertimain(1) = fvert(up) * efvs(1)
  vertimain(2) = - vertimain(3) - vertimain(1)

  return

end function vertimain

!***********************************************************************

function normhorizmain(which,up,efvs)

  ! Called from 'findcoefst' to calculate normal-stress grad terms 
  !      like: d/dx(f(du/dx)), d/dy(f(dv/dy)), etc.  
  ! ... calls FUNCTIONS: horiztermdxdx, horiztermdsdx, horiztermdxds,
  !                      horiztermdsds, horiztermds 
  ! determines coefficients from d/dx(fdu/dx) and d/dy(fdu/dy)

  implicit none

  integer, intent(in) :: which, up
  real (kind = dp), dimension(:,:,:), intent(in) :: efvs

  real (kind = dp), dimension(3,3,3) :: normhorizmain
  real (kind = dp), dimension(3,3,3) :: g, h
  real (kind = dp), dimension(2) :: sumefvsup, sumefvsew, sumefvsns
  real (kind = dp) :: sumefvs

  g = 0.0_dp
  h = 0.0_dp

  sumefvsup = hsum(efvs)
  sumefvsew = sum(sum(efvs,3),1)
  sumefvsns = sum(sum(efvs,2),1)
  sumefvs = sum(efvs)

! for d(f.du/dx)/dx
   
  g(2,:,2) = horiztermdxdx(sumefvsew,cdxdx(1))
  g(:,1:3:2,2) = g(:,1:3:2,2) + horiztermdsdx(dsigmadew(up),sumefvsup,cdsdx(up,1)) 
  g(1:3:2,:,2) = g(1:3:2,:,2) + horiztermdxds(dsigmadew(up),sumefvsew,cdsdx(up,1)) 
  g(:,2,2) = g(:,2,2) + horiztermdsds(dsigmadew(up)**2,sumefvsup,cdsds(up)) 
  g(1:3:2,2,2) = g(1:3:2,2,2) + horiztermds(d2sigmadew2(up)+d2sigmadewdsigma*dsigmadew(up),sumefvs,cds(up))

! for d(f.du/dy)/dy 

  h(2,2,:) = horiztermdxdx(sumefvsns,cdxdx(2))
  h(:,2,1:3:2) = h(:,2,1:3:2) + horiztermdsdx(dsigmadns(up),sumefvsup,cdsdx(up,2)) 
  h(1:3:2,2,:) = h(1:3:2,2,:) + horiztermdxds(dsigmadns(up),sumefvsns,cdsdx(up,2)) 
  h(:,2,2) = h(:,2,2) + horiztermdsds(dsigmadns(up)**2,sumefvsup,cdsds(up))  
  h(1:3:2,2,2) = h(1:3:2,2,2) + horiztermds(d2sigmadns2(up)+d2sigmadnsdsigma*dsigmadns(up),sumefvs,cds(up))

  normhorizmain = g * fourorone(which) + h * oneorfour(which)

  return

end function normhorizmain

!***********************************************************************
   
function croshorizmain(which,up,efvs)

  ! Called from 'findcoefst' to calculate cross-stress grad terms 
  !      like: d/dx(f(du/dy)), d/dy(f(dv/dx)), etc.  
  ! ... calls FUNCTIONS: horiztermdxdy, horiztermdsdx, horiztermdxds, 
  !                      horiztermdsds, horiztermds 
  ! determines coefficients from d/dx(fdu/dy) and d/dy(fdu/dx)

  implicit none

  integer, intent(in) :: which, up
  real (kind = dp), dimension(:,:,:), intent(in) :: efvs

  real (kind = dp), dimension(3,3,3) :: croshorizmain
  real (kind = dp), dimension(3,3,3) :: g = 0.0_dp, h = 0.0_dp
  real (kind = dp), dimension(2) :: sumefvsup, sumefvsew, sumefvsns
  real (kind = dp) :: sumefvs

  g = 0.0_dp
  h = 0.0_dp

  sumefvsup = hsum(efvs)
  sumefvsew = sum(sum(efvs,3),1)
  sumefvsns = sum(sum(efvs,2),1)
  sumefvs = sum(efvs)

! for d(f.du/dy)/dx

  g(2,:,1:3:2) = horiztermdxdy(sumefvsew,cdxdy)
  g(:,2,1:3:2) = g(:,2,1:3:2) + horiztermdsdx(dsigmadew(up),sumefvsup,cdsdx(up,2))
  g(1:3:2,:,2) = g(1:3:2,:,2) + horiztermdxds(dsigmadns(up),sumefvsew,cdsdx(up,1))
  g(:,2,2) = g(:,2,2) + horiztermdsds(dsigmadew(up)*dsigmadns(up),sumefvsup,cdsds(up))   
  g(1:3:2,2,2) = g(1:3:2,2,2) + horiztermds(d2sigmadewdns(up)+d2sigmadnsdsigma*dsigmadew(up),sumefvs,cds(up))

! for d(f.du/dx)/dy 

  h(2,1:3:2,:) = transpose(horiztermdxdy(sumefvsns,cdxdy))
  h(:,1:3:2,2) = h(:,1:3:2,2) + horiztermdsdx(dsigmadns(up),sumefvsup,cdsdx(up,1))
  h(1:3:2,2,:) = h(1:3:2,2,:) + horiztermdxds(dsigmadew(up),sumefvsns,cdsdx(up,2))
  h(:,2,2) = h(:,2,2) + horiztermdsds(dsigmadew(up)*dsigmadns(up),sumefvsup,cdsds(up))
  h(1:3:2,2,2) = h(1:3:2,2,2) + horiztermds(d2sigmadewdns(up)+d2sigmadewdsigma*dsigmadns(up),sumefvs,cds(up))

  croshorizmain = g * twoorone(which) + h * oneortwo(which)

  return

end function croshorizmain

!***********************************************************************

! ***************************************************************************
! start of functions to deal with higher-order boundary conditions at sfc and bed
! ***************************************************************************

function vertimainbc(thck, bcflag, dup, efvs, betasquared, nz)

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma 
   
    implicit none

    real (kind = dp), intent(in) :: dup, thck, betasquared 
    real (kind = dp), intent(in) :: nz                      ! sfc normal vect comp in z-dir
    real (kind = dp), intent(in), dimension(2,2,2) :: efvs
    integer, intent(in), dimension(2) :: bcflag

    real (kind = dp) :: c
    real (kind = dp), dimension(3) :: vertimainbc

    c = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           c = nz / thck / (2*dup) * (len0**2 / thk0**2)   ! value of coefficient

           vertimainbc(:) = 0.0_dp
           vertimainbc(3) = -c 
           vertimainbc(1) = c
           vertimainbc(2) = vertimainbc(3) + vertimainbc(1) ! should = 0


    ! for higher-order BASAL B.C. w/ specified basal traction, add on the necessary source term ...
    if( bcflag(2) == 1 )then

            ! last set of terms is mean visc. of ice nearest to the bed
            vertimainbc(2) = vertimainbc(2)   &
                           + ( betasquared / ( sum( efvs(2,:,:) ) / 4.0_dp ) ) * (len0 / thk0)

    end if


    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this is not often implemented, as it is generally sufficient to implement 
    ! an "almost" no slip BC by just making the coeff. for betasquared very large (and the 
    ! the code converges more quickly/stably in this case than for actual no-slip).
    else if( bcflag(1) == 0 )then

           ! if u,v set to 0, there are no coeff. assoc. with du/digma terms ...
           vertimainbc(:) = 0.0_dp

    end if

    return

end function vertimainbc

!***********************************************************************

function normhorizmainbc(dew,       dns,        &
                         dusrfdew,  dusrfdns,   &
                         dsigmadew, dsigmadns,  &
                         which,     bcflag,     &
                         dup,                   &
                         oneorfour, fourorone)

    ! Determines higher-order surface and basal boundary conditions for LHS of equation.
    ! Gives 3x3x3 coeff. array for either u or v component of velocity, depending on the 
    ! value of the flag 'which'. Example of function call:
    !
    !  g = normhorizmainbc(dusrfew(ew,ns),dusrfnx(ew,ns),dsigmadew(up),dsigmadns(up),which,up,bcflag)	
    !
    ! ... where g is a 3x3x3 array.
    !
    ! 'bcflag' is a 1 x 2 vector to indicate (1) which b.c. is being solved for (surface or bed) and 
    ! (2), if solving for the bed b.c., which type of b.c. to use. For example, bcflag = [ 0, 0 ] 
    ! denotes free sfc bc; bcflag = [ 1, 0 ] denotes basal bc w/ u=v=0, etc. (see also subroutine
    ! "bodyset"). "fourorone" and "oneorfour" are given by vectors: fourorone = [ 4 1 ]; oneorfour = [ 1 4 ].
    ! A single value is chosen from each vector and applied to the calculation of coefficients below.
    ! The "correct" value needed to satisfy the expression is chosen based on the "which" flag, which
    ! takes on a value of 1 for calculations in the x direction and a value of 2 for calculations in 
    ! the y direction. 

    implicit none

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(2) :: oneorfour, fourorone
    real (kind = dp), dimension(3,3,3) :: normhorizmainbc
    real (kind = dp), dimension(3,3,3) :: g
    real (kind = dp) :: c

    integer, intent(in) :: which
    integer, intent(in), dimension(2) :: bcflag

    c = 0.0_dp
    g(:,:,:) = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this handles the case for specified stress at the bed as well, as we 
    ! simply pass in a different value for the normal vector (slope) components (still
    ! called "dusrfdns", "dusrfdew" here, but args passed in are different).
    if( bcflag(1) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated 
           ! with u(1,2,2) and u(3,2,2) ...
           c = ( fourorone(which) * dusrfdew * dsigmadew   &
               + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
           g(3,2,2) = -c 
           g(1,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
           c = fourorone(which) * dusrfdew / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

           c = oneorfour(which) * dusrfdns / (2*dns)
           g(2,2,3) = c
           g(2,2,1) = -c

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    ! note that this requires that rhs(up) be set to 0 as well ...
    else if( bcflag(1) == 0 )then

           g(:,:,:) = 0.0_dp
           g(2,2,2) = 1.0_dp;

    end if

    normhorizmainbc = g

    return

end function normhorizmainbc

!***********************************************************************

function croshorizmainbc(dew,       dns,       &
                         dusrfdew,  dusrfdns,  &
                         dsigmadew, dsigmadns, &
                         which,     bcflag,    &
                         dup,       local_othervel,  &
                         oneortwo,  twoorone)

    ! As described for "normhorizmainbc" above. The vectors "twoorone" and 
    ! "oneortwo" are given by: twoorone = [ 2 1 ]; oneortwo = [ 1 2 ];

    implicit none

    integer, intent(in) :: which
    integer, intent(in), dimension(2) :: bcflag

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in), dimension(2) :: oneortwo, twoorone
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(3,3,3) :: local_othervel
    real (kind = dp), dimension(3,3,3) :: g, croshorizmainbc
    real (kind = dp) :: c

    c = 0.0_dp
    g(:,:,:) = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! NOTE that this handles the case for specified stress at the bed as well, as we 
    ! simply pass in a different value for the normal vector (slope) components (still
    ! called "dusrfdns", "dusrfdew" here, but args passed in are different).
    if( bcflag(1) == 1 )then

           ! first, coeff. that go with du/dsigma, and thus are associated
           ! with u(1,2,2) and u(3,2,2) ...
           c = ( - twoorone(which) * dusrfdew * dsigmadns   &
                 - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
           g(3,2,2) = -c
           g(1,2,2) = c

           ! next, coeff. that go with du/dxhat and du/dyhat terms ...
           c = - oneortwo(which) * dusrfdns / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

           c = - twoorone(which) * dusrfdew / (2*dns)
           g(2,2,3) = c
           g(2,2,1) = -c

    ! for higher-order BASAL B.C. U=V=0, in x ('which'=1) or y ('which'=2) direction ...
    else if( bcflag(1) == 0 )then

       g(:,:,:) = 0.0_dp
       g(2,2,2) = 0.0_dp   ! specify the value of the boundary velocity in the rhs HERE

       ! this forces the multiplication by 'local_otherval' in the main program 
       ! to result in a value of 1, thus leaving the boundary vel. unchanged
       g = g / local_othervel  

    end if

    croshorizmainbc = g

    return

end function croshorizmainbc

!***********************************************************************

function normhorizmainbc_lat(dew,       dns,   &
                             dusrfdew,  dusrfdns,  &
                             dsigmadew, dsigmadns, &
                             which,     what,      &
                             dup,                  &
                             oneorfour, fourorone, &
                             onesideddiff,         &
                             normal,    fwdorbwd)

    ! Analogous to "normhorizmainbc" but for the case of lateral stress (ice shelf)
    ! boundary conditions. Note that the basic form of the equations is the same. 
    ! What changes here is (1) the value of the normal vector that is passed in (at
    ! the sfc and bed we pass in the surface or basal slopes, while at the boundaries
    ! we use the normal vector orientation to the boundary in map view) and (2) we to
    ! to use one sided diffs at the lateral boundaries rather than centerd diffs.

    ! Note that we assume here that du/dz (and thus du/dsigma) is approx. 0 for an ice 
    ! shelf, and also that the sfc/basal slopes of an ice shelf are very flat at/near 
    ! the boundary. Thus, we assume flow is depth independent and we ignore gradients 
    ! in sigma. 

    implicit none

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(2) :: oneorfour, fourorone, normal, fwdorbwd
    real (kind = dp), intent(in), dimension(3) :: onesideddiff

    real (kind = dp), dimension(3,3,3) :: normhorizmainbc_lat
    real (kind = dp), dimension(3,3,3) :: g
    real (kind = dp), dimension(2) :: whichbc
    real (kind = dp) :: c

    integer, intent(in) :: which, what

    c = 0.0_dp; g(:,:,:) = 0.0_dp; whichbc = (/ 0.0_dp, 1.0_dp /)

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    ! (also applies to basal stress bc) 

    ! first, coeff. that go with du/dsigma, and thus are associated with u(1,2,2) and u(3,2,2) ...
    ! ...note that these are stored in an empty column of 'g' (a corner column) so that we don't 
    ! overwrite these values in the case of fwd/bwd horiz. diffs., which require 3 spaces
    c = ( fourorone(which) * dusrfdew * dsigmadew    &
            + oneorfour(which) * dusrfdns * dsigmadns )/(2*dup)
    g(3,3,3) = -c * whichbc(what) 
    g(1,3,3) = c * whichbc(what)

    if( normal(1) .eq. 0.0_dp )then     ! centered in x ...

           c = fourorone(which) * dusrfdew / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) .ne. 0.0_dp )then     ! forward/backward in x ...

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdew / (2*dew)
           g(2,2-int(fwdorbwd(1)),2) = c

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdew / (2*dew)
           g(2,2,2) = c

           c = fourorone(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdew / (2*dew)
           g(2,2+int(fwdorbwd(1)),2) = c

    end if


    if( normal(2) .eq. 0.0_dp ) then   ! centered in y ... 
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = oneorfour(which) * dusrfdns / (2*dns)
           g(1,2,3) = c * whichbc(what)
           g(1,2,1) = -c * whichbc(what)

    elseif( normal(2) .ne. 0.0_dp) then ! forward/backward in y ...

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdns / (2*dns)
           g(1,2,2-int(fwdorbwd(2))) = c

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdns / (2*dns)
           g(1,2,2) = c

           c = oneorfour(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdns / (2*dns)
           g(1,2,2+int(fwdorbwd(2))) = c

    end if

    normhorizmainbc_lat = g

    return

end function normhorizmainbc_lat

!***********************************************************************

function croshorizmainbc_lat (dew,       dns,       &
                              dusrfdew,  dusrfdns,  &
                              dsigmadew, dsigmadns, &
                              which,     what,      &
                              dup,       local_othervel,  &
                              oneortwo,  twoorone,  &
                              onesideddiff,         &
                              normal,    fwdorbwd)

    ! Analagous to "normhorizmainbc_lat" but for cross terms. See notes above.

    implicit none

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in), dimension(2) :: oneortwo, twoorone, fwdorbwd, normal
    real (kind = dp), intent(in), dimension(3) :: onesideddiff
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(3,3,3) :: local_othervel

    integer, intent(in) :: which, what

    real (kind = dp), dimension(3,3,3) :: g, croshorizmainbc_lat
    real (kind = dp), dimension(3) :: gvert
    real (kind = dp), dimension(2) :: whichbc
    real (kind = dp) :: c

    integer, dimension(2) :: inormal

    c = 0.0_dp
    g(:,:,:) = 0.0_dp
    gvert = 0.0_dp
    whichbc = (/ 0.0_dp, 1.0_dp /)
    croshorizmainbc_lat = 0.0_dp

    ! first, coeff. that go with du/dsigma, and thus are associated with u(1,2,2) and u(3,2,2) 
    ! ... note that these are stored in a separate vector (to avoid being overwritten if stored in normal 'g')	
    c = ( - twoorone(which) * dusrfdew * dsigmadns   &
              - oneortwo(which) * dusrfdns * dsigmadew )/(2*dup)
    gvert(3) = -c * whichbc(what)
    gvert(1) = c * whichbc(what)

    if( normal(1) .eq. 0.0_dp )then        ! centered in x ...

           c = -oneortwo(which) * dusrfdns / (2*dew)
           g(2,3,2) = c
           g(2,1,2) = -c

    elseif( normal(1) .ne. 0.0_dp )then    ! forward/backward in x ...
                                           ! (NOTE that x coeff. are stored in g(2,:,:) )

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(1) * dusrfdns / (2*dew)
           g(2,2-int(fwdorbwd(1)),2) = c

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(2) * dusrfdns / (2*dew)
           g(2,2,2) = c

           c = -oneortwo(which) * fwdorbwd(1) * onesideddiff(3) * dusrfdns / (2*dew)
           g(2,2+int(fwdorbwd(1)),2) = c

    end if

    if( normal(2) .eq. 0.0_dp )then    ! centered in y ...
                                       ! (NOTE that y coeff. are stored in g(1,:,:) )

           c = -twoorone(which) * dusrfdew / (2*dns)
           g(1,2,3) = c
           g(1,2,1) = -c

    elseif( normal(2) .ne. 0.0_dp )then ! forward/backward in y ...

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(1) * dusrfdew / (2*dns)
           g(1,2,2-int(fwdorbwd(2))) = c

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(2) * dusrfdew / (2*dns)
           g(1,2,2) = c

           c = -twoorone(which) * fwdorbwd(2) * onesideddiff(3) * dusrfdew / (2*dns)
           g(1,2,2+int(fwdorbwd(2))) = c

    end if

    ! Now rearrange position of coefficients in structure 'g' so that they are multiplied by 
    ! the correct velocity component of 'local_othervel' in 'bodyset' ...
    ! ... this can be done by using the boundary normal vector to shift the indices of the rows/columns
    ! in 'g', in the appropriate direction. First, convert the boundary normal to an integer index ...
    inormal(1) = int( normal(1)/abs(normal(1)) )
    inormal(2) = int( normal(2)/abs(normal(2)) )
    if( abs( inormal(1) ) .ne. 1 )then; inormal(1) = 0; end if
    if( abs( inormal(2) ) .ne. 1 )then; inormal(2) = 0; end if

    croshorizmainbc_lat(2,:,2+inormal(2)) = g(2,:,2)    ! move x-coeffs. appropriate amount
    croshorizmainbc_lat(1,2+inormal(1),:) = g(1,2,:)    ! move y-coeffs. appropriate amount

    ! sum coeffs. that are in same column and flatten so that all coeff. are on level (2,:,:)	
    croshorizmainbc_lat(2,:,:) = croshorizmainbc_lat(2,:,:) + croshorizmainbc_lat(1,:,:)    

    ! set remaining coeff. on this level to to 0 ...
    croshorizmainbc_lat(1,:,:) = 0.0_dp

    ! accounter for vertical terms stored seperately and temporarily in 'gvert'
    croshorizmainbc_lat(1,2+inormal(1),2+inormal(2)) = gvert(1) * whichbc(what)
    croshorizmainbc_lat(3,2+inormal(1),2+inormal(2)) = gvert(3) * whichbc(what)

    return

end function croshorizmainbc_lat

!***********************************************************************

! ---> the following routines are for derivatives in the main body 

function horiztermdxdx(efvs,fact)

  ! this is the d/dx(f.du/dx) and d/dy(f.du/dy) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: fact

  real (kind = dp), dimension(3) :: horiztermdxdx  

  horiztermdxdx(3) = efvs(2) * fact
  horiztermdxdx(1) = efvs(1) * fact
  horiztermdxdx(2) = - horiztermdxdx(3) - horiztermdxdx(1)

  return

end function horiztermdxdx

!***********************************************************************

function horiztermdxdy(efvs,fact)

  ! this is the d/dy(f.du/dx) and d/dx(f.du/dy) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: fact

  real (kind = dp), dimension(3,2) :: horiztermdxdy

  horiztermdxdy(3,2) = efvs(2) * fact 
  horiztermdxdy(2,2) = horiztermdxdy(3,2)
  horiztermdxdy(3,1) = - horiztermdxdy(3,2)
  horiztermdxdy(2,1) = - horiztermdxdy(3,2)

  horiztermdxdy(1,2) = - efvs(1) * fact
  horiztermdxdy(2,2) = horiztermdxdy(2,2) + horiztermdxdy(1,2)
  horiztermdxdy(2,1) = horiztermdxdy(2,1) - horiztermdxdy(1,2)
  horiztermdxdy(1,1) = - horiztermdxdy(1,2)

  return

end function horiztermdxdy

!***********************************************************************

function horiztermdsdx(dsigmadxy,efvs,fact)

  ! this is the d/ds(f.du/dx) and d/ds(f.du/dy) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: dsigmadxy, fact

  real (kind = dp), dimension(3,2) :: horiztermdsdx  

  horiztermdsdx(3,2) = dsigmadxy * efvs(2) * fact
  horiztermdsdx(2,2) = horiztermdsdx(3,2)
  horiztermdsdx(3,1) = - horiztermdsdx(3,2)
  horiztermdsdx(2,1) = - horiztermdsdx(3,2)

  horiztermdsdx(1,2) = - dsigmadxy * efvs(1) * fact
  horiztermdsdx(2,2) = horiztermdsdx(2,2) + horiztermdsdx(1,2)
  horiztermdsdx(2,1) = horiztermdsdx(2,1) - horiztermdsdx(1,2)
  horiztermdsdx(1,1) = - horiztermdsdx(1,2)

  return

end function horiztermdsdx

!***********************************************************************

function horiztermdxds(dsigmadxy,efvs,fact)

  ! this is the d/dx(f.du/ds) and d/dy(f.du/ds) terms

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: dsigmadxy, fact

  real (kind = dp), dimension(2,3) :: horiztermdxds

  horiztermdxds(2,3) = dsigmadxy * efvs(2) * fact
  horiztermdxds(2,2) = horiztermdxds(2,3)
  horiztermdxds(1,3) = - horiztermdxds(2,3)
  horiztermdxds(1,2) = - horiztermdxds(2,3)

  horiztermdxds(2,1) = - dsigmadxy * efvs(1) * fact
  horiztermdxds(2,2) = horiztermdxds(2,2) + horiztermdxds(2,1)
  horiztermdxds(1,2) = horiztermdxds(1,2) - horiztermdxds(2,1)
  horiztermdxds(1,1) = - horiztermdxds(2,1)

  return

end function horiztermdxds

!***********************************************************************

function horiztermdsds(dsigmadxysq,efvs,fact)

  ! this is the d/ds(f.du/ds) term

  implicit none

  real (kind = dp), dimension(2), intent(in) :: efvs
  real (kind = dp), intent(in) :: dsigmadxysq, fact

  real (kind = dp), dimension(3) :: horiztermdsds

  horiztermdsds(3) = dsigmadxysq * efvs(2) * fact
  horiztermdsds(1) = dsigmadxysq * efvs(1) * fact

  horiztermdsds(2) = - horiztermdsds(3) - horiztermdsds(1)

  return

end function horiztermdsds

!***********************************************************************

function horiztermds(d2sigmadxy2etc,efvs,fact)

  ! this is the f.du/ds term

  implicit none

  real (kind = dp), intent(in) :: efvs, d2sigmadxy2etc, fact

  real (kind = dp), dimension(2) :: horiztermds

  horiztermds(2) = d2sigmadxy2etc * efvs * fact
  horiztermds(1) = - horiztermds(2)

  return

end function horiztermds

! ---> end of routines for derivatives in the main body 

!***********************************************************************

subroutine fillsprsemain(inp,locplusup,ptindx,up)

  ! scatter coefficients from 3x3x3 block "g" onto sparse matrix row
  implicit none

  real (kind = dp), dimension(3,3,3), intent(in):: inp
  integer, intent(in) :: locplusup, up
  integer, dimension(6), intent(in) :: ptindx

  ! insert entries to "g" that are on same level
  call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
  call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
  call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
  call putpcgc(inp(2,2,3),ptindx(4)+up,locplusup)  
  call putpcgc(inp(2,2,1),ptindx(5)+up,locplusup)

  ! add points for level above (that is, points in "g"  with a LARGER first index,
  ! which correspond to grid points that are CLOSER TO THE BED than at current level)
  call putpcgc(inp(3,2,2),ptindx(1)+up+1,locplusup)
  call putpcgc(inp(3,3,2),ptindx(2)+up+1,locplusup)
  call putpcgc(inp(3,1,2),ptindx(3)+up+1,locplusup)
  call putpcgc(inp(3,2,3),ptindx(4)+up+1,locplusup)  
  call putpcgc(inp(3,2,1),ptindx(5)+up+1,locplusup)

  ! add points for level below (that is, points in "g" with a SMALLER first index,
  ! which correspond to grid points that are CLOSER TO THE SURFACE than at current level) 
  call putpcgc(inp(1,2,2),ptindx(1)+up-1,locplusup)
  call putpcgc(inp(1,3,2),ptindx(2)+up-1,locplusup)
  call putpcgc(inp(1,1,2),ptindx(3)+up-1,locplusup)
  call putpcgc(inp(1,2,3),ptindx(4)+up-1,locplusup)  
  call putpcgc(inp(1,2,1),ptindx(5)+up-1,locplusup)

  return

end subroutine fillsprsemain

!***********************************************************************

subroutine fillsprsebndy(inp,locplusup,ptindx,up,normal)

  ! scatter coeff. from 3x3x3 block "g" onto sparse matrix row. This subroutine
  ! is specifically for the boundary conditions, which are handled differently
  ! than points in the "main" body of the domain (interior to boundaries).
  implicit none

  integer, intent(in) :: locplusup, up
  integer, dimension(6), intent(in) :: ptindx
  real (kind = dp), dimension(3,3,3), intent(in) :: inp
  real (kind = dp), dimension(2), intent(in) :: normal

  ! at points where mixed centered and one-side diffs. would apply
  if( normal(1) == 0.0_dp )then         ! at boundary normal to y, centered diffs in x 
    if( normal(2) == -1.0_dp )then      ! at boundary w/ normal [0,-1]
           call putpcgc(inp(1,3,3),ptindx(5)+up-1,locplusup)
           call putpcgc( inp(2,3,3)+inp(1,2,1),ptindx(5)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(5)+up+1,locplusup)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup)
    else                                ! at boundary w/ normal [0,1]
           call putpcgc(inp(1,3,3),ptindx(4)+up-1,locplusup)
           call putpcgc(inp(2,3,3)+inp(1,2,3),ptindx(4)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(4)+up+1,locplusup)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup)
    end if
    call putpcgc(inp(1,2,2),ptindx(1)+up,locplusup)
  end if

  if( normal(2) == 0.0_dp )then            ! at boundary normal to x, centered diffs in y 
        if( normal(1) == -1.0_dp )then     ! at boundary w/ normal [-1,0]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup)
           call putpcgc( inp(2,3,3)+inp(2,1,2),ptindx(3)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
        else                                 ! at boundary w/ normal [1,0]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup)
           call putpcgc( inp(2,3,3)+inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
    end if
    call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
  end if

  ! at corners where only one-side diffs. apply
  if( normal(1) .gt. 0.0_dp .and. normal(2) .ne. 0.0_dp )then
    if( normal(2) .gt. 0.0_dp )then      ! corner w/ normal [ 1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(2,3,2)+inp(1,2,3),ptindx(2)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
    else                                 ! corner w/ normal [ 1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(1,2,1)+inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup)
    end if
  end if

  if( normal(1) .lt. 0.0_dp .and. normal(2) .ne. 0.0_dp )then
    if( normal(2) .gt. 0.0_dp )then       ! corner w/ normal [ -1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(1,2,3)+inp(2,1,2),ptindx(3)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup)
    else                                  ! corner w/ normal [ -1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup)
           call putpcgc(inp(2,3,3)+inp(2,1,2)+inp(1,2,1),ptindx(3)+up,locplusup)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup)
    end if
  end if

  return

end subroutine fillsprsebndy

!***********************************************************************

subroutine getlatboundinfo( ew, ns, up, ewn, nsn, upn,  &  
                           thck, loc_array,             &
                           fwdorbwd, normal, loc_latbc)

  ! Calculate map plane normal vector at 45 deg. increments
  ! for regions of floating ice
  implicit none

  integer, intent(in) :: ew, ns, up
  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(ewn-1,nsn-1), intent(in) :: loc_array
  real (kind = dp), dimension(3,3), intent(in) :: thck

  real (kind = dp), dimension(2), intent(out) :: fwdorbwd, normal
  integer, dimension(6), intent(out) :: loc_latbc

  real (kind = dp), dimension(3,3) :: mask, maskcorners
  real (kind = dp), dimension(3,3) :: thckmask
  real (kind = dp), dimension(3) :: testvect
  real (kind = dp) :: phi, deg2rad

  deg2rad = 3.141592654d0 / 180.0d0
  loc_latbc = 0; phi = 0
  mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
  mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
  mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
  maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
  maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

  ! specify new value of 'loc' vector such that fwd/bwd diffs. are set up correctly in sparse matrix
  ! when function 'fillsprsebndy' is called. Also, specify appropriate values for the vectors 'normal'
  ! and 'fwdorbwd', which specify the orientation of the boundary normal and the direction of forward or
  ! backward differencing to be done in the lateral boundary condition functions 'normhorizmainbc_lat'
  ! and 'crosshorizmainbc_lat'

  ! following is algorithm for calculating boundary normal at 45 deg. increments, based on arbitray
  ! boundary shape (based on initial suggestions by Anne LeBrocq)
  where( thck .ne. 0.0d0 )
        thckmask = 0.0_dp
  elsewhere( thck .eq. 0.0d0 )
        thckmask = 1.0d0
  endwhere

  testvect = sum( thckmask * mask, 1 )

    ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 
    ! 90 deg. at 3 O'clock, etc.
    if( sum( sum( thckmask, 1 ) ) .eq. 1.0d0 )then
        phi = sum( sum( thckmask * maskcorners, 1 ) )
    else
        if( any( testvect .eq. 360.0d0 ) )then
            if( sum( testvect ) .eq. 450.0d0 )then
                phi = 45.0d0
            elseif( sum( testvect ) .eq. 630.0d0 )then
                phi = 315.0d0
            else
                phi = 0.0d0
            end if
        elseif( all( testvect .ne. 360 ) )then
            phi = sum( testvect ) / sum( testvect/testvect, testvect .ne. 0.0d0 )
        end if
    end if

    ! define normal vectors and change definition of loc_array based on this angle
    if( phi .eq. 0.0d0 )then
         loc_latbc(1) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         loc_latbc(2) = loc_array(ew+1,ns); loc_latbc(3) = loc_array(ew-1,ns)
         normal = (/ 0.0_dp, 1.0_dp /); fwdorbwd = (/ -1.0_dp, -1.0_dp /)
    elseif( phi .eq. 45.0d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(6) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         normal = (/ 1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ -1.0_dp, -1.0_dp /)
    elseif( phi .eq. 90.0d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(4) = loc_array(ew,ns+1); loc_latbc(5) = loc_array(ew,ns-1)
         normal = (/ 1.0_dp, 0.0_dp /); fwdorbwd = (/ -1.0_dp, -1.0_dp /)
    elseif( phi .eq. 135.0d0 )then
         loc_latbc(1) = loc_array(ew-1,ns); loc_latbc(2) = loc_array(ew,ns); loc_latbc(3) = loc_array(ew-2,ns)
         loc_latbc(6) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns)
         normal = (/ 1.0_dp/sqrt(2.0_dp), -1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ -1.0_dp, 1.0_dp /)
    elseif( phi .eq. 180.0d0 )then
         loc_latbc(1) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns)
         loc_latbc(2) = loc_array(ew+1,ns); loc_latbc(3) = loc_array(ew-1,ns)
         normal = (/ 0.0_dp, -1.0_dp /); fwdorbwd = (/ 1.0_dp, 1.0_dp /)
    elseif( phi .eq. 225.0d0 )then
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(6) = loc_array(ew,ns+1); loc_latbc(4) = loc_array(ew,ns+2); loc_latbc(5) = loc_array(ew,ns);
         normal = (/ -1.0_dp/sqrt(2.0_dp), -1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ 1.0_dp, 1.0_dp /)
    elseif( phi .eq. 270.0d0 )then
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(4) = loc_array(ew,ns+1); loc_latbc(5) = loc_array(ew,ns-1)
         normal = (/ -1.0_dp, 0.0_dp /); fwdorbwd = (/ 1.0_dp, 1.0_dp /)
    else
         loc_latbc(1) = loc_array(ew+1,ns); loc_latbc(2) = loc_array(ew+2,ns); loc_latbc(3) = loc_array(ew,ns)
         loc_latbc(6) = loc_array(ew,ns-1); loc_latbc(4) = loc_array(ew,ns); loc_latbc(5) = loc_array(ew,ns-2)
         normal = (/ -1.0_dp/sqrt(2.0_dp), 1.0_dp/sqrt(2.0_dp) /); fwdorbwd = (/ 1.0_dp, -1.0_dp /)
    end if

  return

end subroutine getlatboundinfo

!***********************************************************************

function indshift( which, ew, ns, up, ewn, nsn, upn, loc_array, thck )  

  ! Subroutine to rearrange indices slightly at sfc,bed, and lateral boundaries,
  ! so that values one index inside of the domain are used for, e.g. eff. visc.
  
  ! Function output is a vector containing necessary index shifts for portions of 'othervel' and 'efvs' 
  ! extracted near domain boundaries. NOTE that this contains duplication of some of the code in the 
  ! subroutine "getlatboundinfo", and the two could be combined at some point.
  implicit none

  integer, intent(in) :: which
  integer, intent(in) :: ew, ns, up, ewn, nsn, upn 
  integer, dimension(ewn-1,nsn-1), intent(in) :: loc_array
  real (kind = dp), dimension(3,3), intent(in) :: thck

  integer, dimension(3) :: indshift
  integer :: upshift = 0, ewshift = 0, nsshift = 0

  real (kind = dp), dimension(3,3) :: mask, maskcorners
  real (kind = dp), dimension(3,3) :: thckmask
  real (kind = dp), dimension(3) :: testvect
  real (kind = dp) :: phi, deg2rad

  deg2rad = 3.141592654d0 / 180.0d0
  mask(:,1) = (/ 0.0d0, 180.0d0, 0.0d0 /)
  mask(:,2) = (/ 270.0d0, 0.0d0, 90.0d0 /)
  mask(:,3) = (/ 0.0d0, 360.0d0, 0.0d0 /)
  maskcorners(:,1) = (/ 225.0d0, 0.0d0, 135.0d0 /)
  maskcorners(:,2) = (/ 0.0d0, 0.0d0, 0.0d0 /)
  maskcorners(:,3) = (/ 315.0d0, 0.0d0, 45.0d0 /)

  if( up == 1 )then   !! first treat bed/sfc, which aren't complicated
      upshift = 1
  elseif( up == upn )then
      upshift = -1
  else
      upshift = 0
  end if

  select case(which)

      case(0)   !! internal to lateral boundaries; no shift to ew,ns indices

          ewshift = 0; nsshift = 0;

      case(1)   !! at lateral boundaries; shift to ew,ns may be non-zero

          where( thck .ne. 0.0d0 ) 
            thckmask = 0.0_dp
          elsewhere( thck .eq. 0.0d0 )
            thckmask = 1.0d0
          endwhere

          testvect = sum( thckmask * mask, 1 )

        ! calculate the angle of the normal in cart. (x,y) system w/ 0 deg. at 12 O'clock, 90 deg. at 3 O'clock, etc.
        if( sum( sum( thckmask, 1 ) ) .eq. 1.0d0 )then
            phi = sum( sum( thckmask * maskcorners, 1 ) )
        else
            if( any( testvect .eq. 360.0d0 ) )then
                if( sum( testvect ) .eq. 450.0d0 )then
                    phi = 45.0d0
                elseif( sum( testvect ) .eq. 630.0d0 )then
                    phi = 315.0d0
                else
                    phi = 0.0d0
                end if
            elseif( all( testvect .ne. 360 ) )then
                phi = sum( testvect ) / sum( testvect/testvect, testvect .ne. 0.0d0 )
            end if
        end if

        ! define shift to indices based on this angle 
        if( phi .eq. 0.0d0 )then
            nsshift = -1; ewshift = 0
        elseif( phi .eq. 45.0d0 )then
            nsshift = -1; ewshift = -1
        elseif( phi .eq. 90.0d0 )then
            nsshift = 0; ewshift = -1
        elseif( phi .eq. 135.0d0 )then
            nsshift = 1; ewshift = -1
        elseif( phi .eq. 180.0d0 )then
            nsshift = 1; ewshift = 0
        elseif( phi .eq. 225.0d0 )then
            nsshift = 1; ewshift = 1
        elseif( phi .eq. 270.0d0 )then
            nsshift = 0; ewshift = 1
        elseif( phi .eq. 315.0d0 )then
            nsshift = -1; ewshift = 1
        end if

  end select

  indshift = (/ upshift, ewshift, nsshift /)

  return

end function indshift

!***********************************************************************

subroutine calcbetasquared (whichbabc,               & 
                            dew,         dns,        &
                            ewn,         nsn,        &
                            lsrf,        topg,       &
                            thck,                    &
                            thisvel,     othervel,   &
                            minTauf, beta,           &
                            betasquared, betafile) 

  ! subroutine to calculate map of betasquared sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").
  implicit none

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), intent(in), dimension(:,:) :: lsrf, topg, thck
  real (kind = dp), intent(in), dimension(:,:) :: thisvel, othervel, minTauf, beta

  real (kind = dp), intent(out), dimension(ewn-1,nsn-1) :: betasquared

  character (len=30), intent(in), optional :: betafile
  real (kind = dp) :: smallnum = 1.0d-4
  real (kind = dp), dimension(ewn) :: grounded
  real (kind = dp) :: alpha, dx, thck_gl, betalow, betahigh, roughness
  integer :: ew, ns

  select case(whichbabc)

    case(0)     ! constant value; useful for debugging and test cases

      betasquared = 1.0d0

    case(1)     ! simple pattern; also useful for debugging and test cases
                ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      betasquared = 1.0d10

      do ew=1, ewn-1; do ns=10, nsn-10
        betasquared(ew,ns) = 1.0d0
      end do; end do

    case(2)     ! take input value for till yield stress and force betasquared to be implemented such
                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).
     
      betasquared = minTauf / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

    case(3)     ! circular ice shelf: set B^2 ~ 0 except for at center, where B^2 >> 0 to enforce u,v=0 there

      betasquared = 1.0d-5
      betasquared( (ewn-1)/2:(ewn-1)/2+1, (nsn-1)/2:(nsn-1)/2+1 ) = 1.0d10

    case(4)    ! frozen (u=v=0) ice-bed interface

      betasquared = 1.0d10

    case(5)    ! use value passed in externally from CISM (NOTE not dimensional when passed in) 

      ! scale CISM input value to dimensional units of (Pa yrs 1/m)
      betasquared = beta * scyr * vel0 * len0 / (thk0**2)    

      ! this is a check for NaNs, which indicate, and are replaced by no slip
      where ( betasquared /= betasquared )  
        betasquared = 1.0d10
      end where    

    case default    ! frozen (u=v=0) ice-bed interface

      betasquared = 1.0d10

  end select
  
  ! convert whatever the specified value is to dimensional units of (Pa s m^-1 ) 
  ! and then non-dimensionalize using PP dyn core specific scaling.
  betasquared = ( betasquared * scyr ) / ( tau0_glam * tim0 / len0 ) 

end subroutine calcbetasquared


!***********************************************************************

function vertintg(upn, sigma, in) 
 
  implicit none 
 
  integer, intent(in) :: upn
  real (kind = dp), dimension(:), intent(in) :: sigma
  real (kind = dp), dimension(:), intent(in) :: in 
  real (kind = dp) :: vertintg 
 
  integer :: up 
 
  vertintg = 0.0d0 
 
  do up = upn-1, 1, -1 
    vertintg = vertintg + sum(in(up:up+1)) * dups(up) 
  end do 
 
  vertintg = vertintg / 2.0d0 
 
  return 
 
end function vertintg 

!***********************************************************************

subroutine geom2derscros(dew,  dns,   &
                         ipvr, stagthck, opvrewns)        

  ! geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'       
 
  implicit none 

  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), intent(out), dimension(:,:) :: opvrewns
  real (kind = dp), intent(in), dimension(:,:) :: ipvr, stagthck
 
  ! consider replacing by a loop over ewn, nsn?
  where (stagthck /= 0.0d0)
    opvrewns = (eoshift(eoshift(ipvr,1,0.0_dp,2),1,0.0_dp,1) + ipvr   &
               - eoshift(ipvr,1,0.0_dp,1) - eoshift(ipvr,1,0.0_dp,2)) / (dew*dns)
  elsewhere
    opvrewns = 0.0d0
  end where
 
  return
 
end subroutine geom2derscros

!***********************************************************************

subroutine geom2ders(ewn,    nsn,  &
                     dew,    dns,  &
                     ipvr,   stagthck,  &
                     opvrew, opvrns)       

  ! geometric 1st deriv. for generic input variable 'ipvr', 
  ! output as 'opvr' (includes 'upwinding' for boundary values)

  implicit none 
 
  integer, intent(in) :: ewn, nsn 
  real (kind = dp), intent(in) :: dew, dns
  real (kind = dp), intent(out), dimension(:,:) :: opvrew, opvrns
  real (kind = dp), intent(in), dimension(:,:) :: ipvr, stagthck
 
  integer :: ew, ns 
  real (kind = dp) :: dewsq4, dnssq4
 
  integer :: pt(2)
 
  dewsq4 = 4.0d0 * dew * dew
  dnssq4 = 4.0d0 * dns * dns

  do ns = 2, nsn-2 
  do ew = 2, ewn-2
    if (stagthck(ew,ns) .gt. 0.0d0) then
      opvrew(ew,ns) = centerew(ew,ns,ipvr,dewsq4)
      opvrns(ew,ns) = centerns(ew,ns,ipvr,dnssq4)
    else
      opvrew(ew,ns) = 0.0d0
      opvrns(ew,ns) = 0.0d0
    end if
  end do
  end do
 
  ! *** 2nd order boundaries using upwinding
 
  do ew = 1, ewn-1, ewn-2
 
    pt = whichway(ew)
 
    do ns = 2, nsn-2 
      if (stagthck(ew,ns) .gt. 0.0d0) then
        opvrew(ew,ns) = boundyew(ns,pt,ipvr,dewsq4)
        opvrns(ew,ns) = centerns(ew,ns,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.0d0
        opvrns(ew,ns) = 0.0d0
      end if
    end do
 
  end do
 
  do ns = 1, nsn-1, nsn-2
 
    pt = whichway(ns)
 
    do ew = 2, ewn-2  
      if (stagthck(ew,ns) .gt. 0.0d0) then
        opvrew(ew,ns) = centerew(ew,ns,ipvr,dewsq4)
        opvrns(ew,ns) = boundyns(ew,pt,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.0d0
        opvrns(ew,ns) = 0.0d0
      end if
    end do
 
  end do
 
  do ns = 1, nsn-1, nsn-2
    do ew = 1, ewn-1, ewn-2
      if (stagthck(ew,ns) .gt. 0.0d0) then
        pt = whichway(ew)
        opvrew(ew,ns) = boundyew(ns,pt,ipvr,dewsq4)
        pt = whichway(ns)
        opvrns(ew,ns) = boundyns(ew,pt,ipvr,dnssq4)
      else
        opvrew(ew,ns) = 0.0d0
        opvrns(ew,ns) = 0.0d0
      end if
    end do
  end do
 
end subroutine geom2ders
 
!***********************************************************************

  function centerew(ew, ns, ipvr, dewsq4)
 
    implicit none

    integer, intent(in) :: ew, ns 
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dewsq4
    real (kind = dp) :: centerew
 
    centerew = (sum(ipvr(ew+2,ns:ns+1)) + sum(ipvr(ew-1,ns:ns+1)) - &
                sum(ipvr(ew+1,ns:ns+1)) - sum(ipvr(ew,ns:ns+1))) / dewsq4
 
    return
    
  end function centerew 
 
!***********************************************************************

  function centerns(ew, ns, ipvr, dnssq4)
 
    implicit none
 
    integer, intent(in) :: ew, ns 
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dnssq4
    real (kind = dp) :: centerns
 
    centerns = (sum(ipvr(ew:ew+1,ns+2)) + sum(ipvr(ew:ew+1,ns-1)) - &
                sum(ipvr(ew:ew+1,ns+1)) - sum(ipvr(ew:ew+1,ns))) / dnssq4
 
    return
    
  end function centerns 
 
!***********************************************************************

  function boundyew(ns,pt,ipvr,dewsq4)
 
    implicit none

    integer, intent(in) :: ns
    integer, intent(in) :: pt(2)
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dewsq4
    real (kind = dp) :: boundyew
 
    boundyew = pt(1) * (3.0d0 * sum(ipvr(pt(2),ns:ns+1)) - 7.0d0 * sum(ipvr(pt(2)+pt(1),ns:ns+1)) + &
               5.0d0 * sum(ipvr(pt(2)+2*pt(1),ns:ns+1)) - sum(ipvr(pt(2)+3*pt(1),ns:ns+1))) / dewsq4
 
    return
 
  end function boundyew
 
!***********************************************************************

  function boundyns(ew,pt,ipvr,dnssq4)
 
    implicit none
 
    integer, intent(in) :: ew
    integer, intent(in) :: pt(2)
    real(kind=dp), intent(in) :: ipvr(:,:)
    real(kind=dp), intent(in) :: dnssq4
    real (kind = dp) :: boundyns
 
    boundyns = pt(1) * (3.0d0 * sum(ipvr(ew:ew+1,pt(2))) - 7.0d0 * sum(ipvr(ew:ew+1,pt(2)+pt(1))) + &
               5.0d0 * sum(ipvr(ew:ew+1,pt(2)+2*pt(1))) - sum(ipvr(ew:ew+1,pt(2)+3*pt(1)))) / dnssq4
 
    return
 
  end function boundyns
 
!***********************************************************************

  function whichway(i)
 
    implicit none
 
    integer, intent(in) :: i
    integer :: whichway(2) 
 
    if (i == 1) then 
      whichway = (/1,1/)
    else
      whichway = (/-1,i+1/)
    end if
 
    return
 
  end function whichway
 

!***********************************************************************

    function hsum(inp) 
 
      implicit none
 
      real (kind = dp), dimension(:,:,:), intent(in) :: inp
      real (kind = dp), dimension(size(inp,dim=1)) :: hsum
 
      hsum = sum(sum(inp(:,:,:),dim=3),dim=2)
 
      return 
 
    end function hsum

!***********************************************************************

subroutine putpcgc(value,col,row) 
 
  implicit none
 
  integer, intent(in) :: row, col 
  real (kind = dp), intent(in) :: value 

      if (value /= 0.0d0) then
        pcgval(ct) = value
        pcgcol(ct) = col
        pcgrow(ct) = row
        ct = ct + 1
      end if

  return
 
end subroutine putpcgc 

!***********************************************************************

end module glam_strs2

!***********************************************************************
