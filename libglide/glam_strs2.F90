! "glam_strs2.F90"
!
! 3d velocity calculation based on Blatter/Pattyn, 1st-order equations, by Tony Payne (Univ.
! of Bristol) and Steve Price (Univ. of Bristol / Los Alamos Nat. Lab.). Boundary conditions
! available include periodic (lateral), free surface, zero slip at bed, specified basal 
! traction at bed, and specified basal yield stress at bed (all three of which are implemented
! through various verions of the specified traction b.c.)
! include macros for glide mask definitions
#include "glide_mask.inc"
#include "config.inc"
!***********************************************************************
module glam_strs2
!***********************************************************************

use glimmer_paramets, only : dp
use glimmer_physcon,  only : gn, rhoi, rhoo, grav, pi, scyr
use glimmer_paramets, only : thk0, len0, vel0, vis0, vis0_glam, tim0, evs0, tau0
use glimmer_log,      only : write_log
use glide_mask
use glimmer_sparse_type
use glimmer_sparse
use glide_types

implicit none

  integer, save :: locplusup
  logical, save :: lateralboundry = .false.
  integer, dimension(6), save :: loc_latbc

  real (kind = dp), allocatable, dimension(:,:,:),     save :: flwafact
  real (kind = dp), allocatable, dimension(:),         save :: dups
  real (kind = dp), allocatable, dimension(:,:,:,:,:), save :: corr
  real (kind = dp), allocatable, dimension(:,:,:,:),   save :: usav
  real (kind = dp), dimension(2),                      save :: usav_avg
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


  real (kind = dp), allocatable, dimension(:,:,:), save  :: ughost 
  real (kind = dp), allocatable, dimension(:,:,:), save  :: vghost

  ! coeff. for forward differencing template, used for stress bcs at lateral boundaries
  real (kind = dp), dimension(3), parameter ::   &
           onesideddiff = (/ -3.0_dp, 4.0_dp, -1.0_dp /)

  ! geometric 2nd and cross-derivs
  real (kind = dp), dimension(:,:), allocatable :: &
              d2thckdew2, d2usrfdew2, d2thckdns2, d2usrfdns2, d2thckdewdns, d2usrfdewdns

  ! variables for plastic-till basal BC iteration using Newton method
  real (kind = dp), dimension(:,:,:), allocatable :: velbcvect, plastic_coeff_lhs, plastic_coeff_rhs, &
                                                     plastic_rhs, plastic_resid
  real (kind = dp), dimension(:,:,:,:), allocatable :: ghostbvel

  ! variables for use in sparse matrix calculation
  real (kind = dp), dimension(:), allocatable :: pcgval, rhsd 
  integer, dimension(:), allocatable :: pcgcol, pcgrow
  integer, dimension(2) :: pcgsize
  ! additional storage needed for off diagonal blocks when using JFNK for nonlinear iteration 
  real (kind = dp), dimension(:), allocatable :: pcgvaluv, pcgvalvu
  integer, dimension(:), allocatable :: pcgcoluv, pcgrowuv, pcgcolvu, pcgrowvu
  integer :: ct, ct2

!RN_20100125: The following are for Trilinos:
  ! This flag switches between: 
  !  0: load Triad format and call Trilinos through sparse_easy_solve
  !  1: load directly into Trilinos format, skipping Triad matrix altogether
  !     Note: only option 1 will work for distributed fortran fill
!  integer :: load_directly_into_trilinos = 1

  !*sfp* NOTE: these redefined here so that they are "in scope" and can avoid being passed as args
  integer :: whatsparse ! needed for putpgcg()
  integer :: nonlinear  ! flag for indicating type of nonlinar iteration (Picard vs. JFNK)
  logical, save :: storeoffdiag = .false. ! true only if using JFNK solver and block, off diag coeffs needed
  logical, save :: calcoffdiag = .false. 

  real (kind = dp) :: linearSolveTime = 0
  real (kind = dp) :: totalLinearSolveTime = 0 ! total linear solve time

!***********************************************************************

contains

!***********************************************************************


subroutine glam_velo_fordsiapstr_init( ewn,   nsn,   upn,    &
                                       dew,   dns,           &
                                       sigma)

    ! Allocate arrays and initialize variables.
    implicit none

    integer, intent(in) :: ewn, nsn, upn
    real (kind = dp), intent(in) :: dew, dns

    real (kind = dp), dimension(:), intent(in)  :: sigma

    integer :: up

    allocate( dup(upn) )
    allocate( dupm(upn) )
    allocate( cvert(upn) )
    allocate( cdsdx(upn,2) )
    allocate( cdsds(upn) )
    allocate( cds(upn) )
    allocate( fvert(upn) )
    allocate(ughost(2,ewn-1,nsn-1))
    allocate(vghost(2,ewn-1,nsn-1))

    ! NOTE: "dup", the sigma coordinate spacing is defined as a vector to allow it to 
    ! be read in from file for use with non-constant vertical grid spacing. Currently, this
    ! is not working, so the code will not give accurate results if the sigma coordinate is
    ! not regularly spaced. 
    dup = (/ ( (sigma(2)-sigma(1)), up = 1, upn) /)
    dupm = - 0.25_dp / dup
!whl - Moved stagsigma calculation to glide_setup module
!!    stagsigma(1:upn-1) = (sigma(1:upn-1) + sigma(2:upn)) / 2.0_dp

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

  ! allocate/initialize variables for plastic-till basal BC iteration using Newton method
    allocate(velbcvect(2,ewn-1,nsn-1),plastic_coeff_rhs(2,ewn-1,nsn-1),plastic_coeff_lhs(2,ewn-1,nsn-1), &
            plastic_rhs(2,ewn-1,nsn-1), plastic_resid(1,ewn-1,nsn-1) )
    allocate(ghostbvel(2,3,ewn-1,nsn-1))        !! for saving the fictious basal vels at the bed !!

    plastic_coeff_rhs(:,:,:) = 0.0d0
    plastic_coeff_lhs(:,:,:) = 0.0d0
    plastic_rhs(:,:,:) = 0.0d0
    plastic_resid(:,:,:) = 0.0d0
    ghostbvel(:,:,:,:) = 0.0d0
    velbcvect(:,:,:) = 0.0d0

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
                                 btraction,              &
                                 umask,                  &
                                 whichbabc,              &
                                 whichefvs,              &
                                 whichresid,             &
                                 whichnonlinear,         &
                                 whichsparse,            &
                                 periodic_ew,periodic_ns,&
                                 beta,                   &
                                 uvel,     vvel,         &
                                 uflx,     vflx,         &
                                 efvs, tstep )


  implicit none

  integer, intent(in) :: ewn, nsn, upn, tstep  ! JFL to be removed
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
  real (kind = dp), dimension(:,:,:), intent(inout) :: btraction            ! consistent basal traction array
  real (kind = dp), dimension(:,:,:), intent(in)  :: flwa                   ! flow law rate factor

  ! This is the betasquared field from CISM (externally specified), and should eventually
  ! take the place of the subroutine 'calcbetasquared' below. For now, there is simply an option
  ! in the subroutine 'calcbetasquared' (case 9) to use this external, CISM specified value for
  ! the betasquared field as opposed to one of the values calculated internally.
  real (kind = dp), dimension(:,:),   intent(in)  :: beta

  integer, intent(in) :: whichbabc    ! options for betasquared field to use
  integer, intent(in) :: whichefvs    ! options for efvs calculation (calculate it or make it uniform)
  integer, intent(in) :: whichresid   ! options for method to use when calculating vel residul
  integer, intent(in) :: whichnonlinear  ! options for which method for doing elliptic solve
  integer, intent(in) :: whichsparse  ! options for which method for doing elliptic solve
  logical, intent(in) :: periodic_ew, periodic_ns  ! options for applying periodic bcs or not

  real (kind = dp), dimension(:,:,:), intent(inout) :: uvel, vvel  ! horiz vel components: u(z), v(z)
  real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx  ! horiz fluxs: u_bar*H, v_bar*H
  real (kind = dp), dimension(:,:,:), intent(out) :: efvs        ! effective viscosity

  integer :: ew, ns, up     ! counters for horiz and vert do loops

  real (kind = dp), parameter :: minres = 1.0d-4    ! assume vel fields converged below this resid 
  real (kind = dp), parameter :: NL_tol = 1.0d-06   ! to have same criterion
                                                    ! than with JFNK
  real (kind = dp), save, dimension(2) :: resid     ! vector for storing u resid and v resid 
  real (kind = dp) :: plastic_resid_norm = 0.0d0    ! norm of residual used in Newton-based plastic bed iteration

  integer, parameter :: cmax = 100                  ! max no. of iterations
  integer :: counter, linit                         ! iteation counter, ???
  character(len=100) :: message                     ! error message

  ! variables used for incorporating generic wrapper to sparse solver
  type(sparse_matrix_type) :: matrix
  real (kind = dp), dimension(:), allocatable :: answer, uk_1, vk_1, F
  real (kind = dp) :: err, L2norm, L2square, NL_target
  integer :: iter, pic
  integer , dimension(:), allocatable :: g_flag ! jfl flag for ghost cells

  ! AGS: partition information for distributed solves
  integer, allocatable, dimension(:) :: myIndices
  integer :: mySize = -1

  ! RN_20100125: assigning value for whatsparse, which is needed for putpcgc()
  whatsparse = whichsparse

  ! assign value for nonlinear iteration flag
  nonlinear = whichnonlinear

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

!!!!!!!!!! Boundary conditions HACKS section !!!!!!!!!!!!!

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

!! hack for basal processes submodel test case, to avoid floatation at downstream
!! end yet still allow for application of a floating ice bc there
!  do ns=1,nsn-1; do ew=1,ewn-1
!      if( umask(ew,ns) == 37 )then
!          umask(ew,ns) = 41
!      endif
!  end do; end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! allocate space for storing temporary across-flow comp of velocity
  allocate(tvel(upn,ewn-1,nsn-1))
  tvel = 0.0_dp

  ! allocate space for variables used by 'mindcrash' function (unstable manifold correction)
  allocate(corr(upn,ewn-1,nsn-1,2,2),usav(upn,ewn-1,nsn-1,2))

  ! make an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 20

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
     ! AGS: Get partition -- later this will be known by distributed glimmer
!      call dopartition(pcgsize(1), mySize) ! JCC - No Trilinos support yet
     allocate(myIndices(mySize))
!      call getpartition(mySize, myIndices) ! JCC - No Trilinos support yet

     ! Now send this partition to Trilinos initialization routines
!      call inittrilinos(20, mySize, myIndices) ! JCC - No Trilinos support yet

     !No Triad matrix needed in this case -- save on memory alloc
     pcgsize(2) = 1

     deallocate(myIndices)
  endif

!==============================================================================
! RN_20100126: End of the block
!==============================================================================

  ! allocate space matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)),rhsd(pcgsize(1)), &
            pcgval(pcgsize(2)))

  allocate(matrix%row(pcgsize(2)), matrix%col(pcgsize(2)), &
            matrix%val(pcgsize(2)), answer(pcgsize(1)))

  allocate( uk_1(pcgsize(1)), vk_1(pcgsize(1)), &
            F(2*pcgsize(1)), g_flag(pcgsize(1)) ) ! jfl for res calc.

  ! set residual and iteration counter to initial values
  resid = 1.0_dp
  counter = 1
  L2norm = 1.0d20
  linit = 0

  ! print some info to the screen to update on iteration progress
  print *, ' '
  print *, 'Running Payne/Price higher-order dynamics solver'
  print *, ' '
!  print *, 'iter #     uvel resid         vvel resid       target resid'
  print *, 'iter #     resid (L2 norm)       target resid'
  print *, ' '

  ! ****************************************************************************************
  ! START of Picard iteration
  ! ****************************************************************************************

  call ghost_preprocess( ewn, nsn, upn, uindx, ughost, vghost, &
                         uk_1, vk_1, uvel, vvel, g_flag) ! jfl_20100430

  ! Picard iteration; continue iterating until resid falls below specified tolerance
  ! or the max no. of iterations is exceeded

  do while ( L2norm .ge. NL_target .and. counter < cmax)    ! use L2 norm for resid calculation
  !do while ( maxval(resid) > minres .and. counter < cmax)   ! standard residual calculation
  !do while ( resid(1) > minres .and. counter < cmax)        ! standard residual (for 1d solutions where d*/dy=0) 

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
                     beta, btraction,             &
                     counter, 0 )

    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, vvel )

!==============================================================================
! jfl 20100412: residual for v comp: Fv= A(u^k-1,v^k-1)v^k-1 - b(u^k-1,v^k-1)  
!==============================================================================

     call res_vect( matrix, vk_1, rhsd, size(rhsd), counter, g_flag, L2square, whichsparse ) ! JCC - No Trilinos support yet

      L2norm  = L2square
      F(1:pcgsize(1)) = vk_1(:)
      
!   call output_res(ewn,nsn,upn,uindx,counter,size(vk_1),vk_1, 2) ! JFL

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

  if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
     call sparse_easy_solve(matrix, rhsd, answer, err, iter, whichsparse)
  else
!      call solvewithtrilinos(rhsd, answer, linearSolveTime) ! JCC - No Trilinos support yet
     totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
!     write(*,*) 'Total linear solve time so far', totalLinearSolveTime
  endif

!==============================================================================
! RN_20100129: End of the block
!==============================================================================

    vk_1 = answer ! jfl for residual calculation

    ! put vels and coeffs from sparse vector format (soln) back into 3d arrays
    call solver_postprocess( ewn, nsn, upn, 2, uindx, answer, tvel, ghostbvel )

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
                     beta, btraction,             &
                     counter, 0 )


    ! put vels and coeffs from 3d arrays into sparse vector format
    call solver_preprocess( ewn, nsn, upn, uindx, matrix, answer, uvel )

!==============================================================================
! jfl 20100412: residual for u comp: Fu= C(u^k-1,v^k-1)u^k-1 - d(u^k-1,v^k-1)  
!==============================================================================

!     call res_vect( matrix, uk_1, rhsd, size(rhsd), counter, g_flag, L2square, whichsparse ) ! JCC - No Trilinos support yet

    L2norm = sqrt(L2norm + L2square)
    F(pcgsize(1)+1:2*pcgsize(1)) = uk_1(:) ! F = [ Fv, Fu ]

!    print *, 'L2 with/without ghost (k)= ', counter, &
!              sqrt(DOT_PRODUCT(F,F)), L2norm
!    if (counter .le. 2) NL_target = NL_tol * L2norm
!    if (counter == 1) NL_target = NL_tol * L2norm
    if (counter == 1) NL_target = 1.0d-4 

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

  if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
     call sparse_easy_solve(matrix, rhsd, answer, err, iter, whichsparse)
  else
!      call solvewithtrilinos(rhsd, answer, linearSolveTime) ! JCC - No Trilinos support yet
     totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
!     write(*,*) 'Total linear solve time so far', totalLinearSolveTime
  endif

!==============================================================================
! RN_20100129: End of the block
!==============================================================================


    uk_1 = answer ! jfl for residual calculation


    ! put vels and coeffs from sparse vector format (soln) back into 3d arrays
    call solver_postprocess( ewn, nsn, upn, 1, uindx, answer, uvel, ghostbvel )

    ! call fraction of assembly routines, passing current vel estimates (w/o manifold
    ! correction!) to calculate consistent basal tractions

    call findcoefstr(ewn,  nsn,   upn,            &
                     dew,  dns,   sigma,          &
                     2,           efvs,           &
                     tvel,        uvel,           &
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
                     beta, btraction,             &
                     counter, 1 )

   call findcoefstr(ewn,  nsn,   upn,             &
                     dew,  dns,   sigma,          &
                     1,           efvs,           &
                     uvel,        tvel,           &
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
                     beta, btraction,             &
                     counter, 1 )

    call plasticbediteration( ewn, nsn, uvel(upn,:,:), tvel(upn,:,:), btraction, minTauf, &
                              plastic_coeff_lhs, plastic_coeff_rhs, plastic_rhs, plastic_resid )

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

!    ! output the iteration status: iteration number, max residual, and location of max residual
!    ! (send output to the screen or to the log file, per whichever line is commented out) 
!    print '(i4,3g20.6)', counter, resid(1), resid(2), minres
!    !write(message,'(" * strs ",i3,3g20.6)') counter, resid(1), resid(2), minres
!    !call write_log (message)

    print '(i4,3g20.6)', counter, L2norm, NL_target    ! Output when using L2norm for convergence

    counter = counter + 1   ! advance the iteration counter

  end do

  ! ****************************************************************************************
  ! END of Picard iteration
  ! ****************************************************************************************

  call ghost_postprocess( ewn, nsn, upn, uindx, uk_1, vk_1, &
                          ughost, vghost )

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
  deallocate(uk_1, vk_1, F, g_flag) ! jfl 

  return

end subroutine glam_velo_fordsiapstr

!***********************************************************************

subroutine JFNK                 (ewn,      nsn,    upn,  &
                                 dew,      dns,          &
                                 sigma,    stagsigma,    &
                                 thck,     usrf,         &
                                 lsrf,     topg,         &
                                 dthckdew, dthckdns,     &
                                 dusrfdew, dusrfdns,     & 
                                 dlsrfdew, dlsrfdns,     &
                                 stagthck, flwa,         & 
                                 mintauf,                & 
                                 btraction,              & 
                                 umask,                  & 
                                 whichbabc,              &
                                 whichefvs,              &
                                 whichresid,             &
                                 whichnonlinear,         &
                                 whichsparse,            &
                                 periodic_ew,periodic_ns,&
                                 beta,                   & 
                                 uvel,     vvel,         &
                                 uflx,     vflx,         &
                                 efvs, tstep )


  implicit none

  integer, intent(in) :: ewn, nsn, upn, tstep
  integer, dimension(:,:),   intent(inout)  :: umask  !*sfp* replaces the prev., internally calc. mask
                                                      ! ... 'inout' status allows for a minor alteration
                                                      ! to cism defined mask, which don't necessarily 
                                                      ! associate all/any boundaries as a unique mask value.
  real (kind = dp), intent(in) :: dew, dns

  real (kind = dp), dimension(:),     intent(in)  :: sigma, stagsigma
  real (kind = dp), dimension(:,:),   intent(in)  :: thck, usrf, lsrf, topg
  real (kind = dp), dimension(:,:),   intent(in)  :: dthckdew, dthckdns
  real (kind = dp), dimension(:,:),   intent(in)  :: dusrfdew, dusrfdns
  real (kind = dp), dimension(:,:),   intent(in)  :: dlsrfdew, dlsrfdns
  real (kind = dp), dimension(:,:),   intent(in)  :: stagthck
  real (kind = dp), dimension(:,:),   intent(in)  :: minTauf
  real (kind = dp), dimension(:,:,:), intent(inout) :: btraction            ! consistent basal traction array
  real (kind = dp), dimension(:,:,:), intent(in)  :: flwa
  
  !*sfp* This is the betasquared field from CISM (externally specified), and should eventually
  ! take the place of the subroutine 'calcbetasquared' below (for now, using this value instead
  ! will simply be included as another option within that subroutine) 
  real (kind = dp), dimension(:,:),   intent(in)  :: beta 

  integer, intent(in) :: whichbabc
  integer, intent(in) :: whichefvs
  integer, intent(in) :: whichresid
  integer, intent(in) :: whichnonlinear
  integer, intent(in) :: whichsparse
  logical, intent(in) :: periodic_ew, periodic_ns

  real (kind = dp), dimension(:,:,:), intent(out) :: uvel, vvel
  real (kind = dp), dimension(:,:),   intent(out) :: uflx, vflx
  real (kind = dp), dimension(:,:,:), intent(out) :: efvs

  integer :: ew, ns, up, nele, k

  real (kind = dp), parameter :: NL_tol = 1.0d-06

  integer, parameter :: kmax = 100, img = 20, img1 = img+1
  character(len=100) :: message

!*sfp* needed to incorporate generic wrapper to solver
  type(sparse_matrix_type) :: matrixA, matrixC, matrixtp, matrixAuv, matrixAvu
  real (kind = dp), dimension(:), allocatable :: answer, uk_1, vk_1
  real (kind = dp), dimension(:), allocatable :: vectp, uk_1_plus, vk_1_plus
  real (kind = dp), dimension(:), allocatable :: dx, F, F_plus
  real (kind = dp), dimension(:), allocatable :: wk1, wk2, rhs
  real (kind = dp), dimension(:,:), allocatable :: vv, wk
  real (kind = dp) :: L2norm, L2norm_wig, tol, gamma_l, epsilon,NL_target
  real (kind = dp) :: crap
  integer :: tot_its, itenb, maxiteGMRES, iout, icode
  integer , dimension(:), allocatable :: g_flag ! jfl flag for ghost cells

  ! AGS: partition information for distributed solves
  integer, allocatable, dimension(:) :: myIndices
  integer :: mySize = -1

  ! RN_20100125: assigning value for whatsparse, which is needed for putpcgc()
  whatsparse = whichsparse
  nonlinear = whichnonlinear

  ! *sfp** geometric 1st deriv. for generic input variable 'ipvr',
  !      output as 'opvr' (includes 'upwinding' for boundary values)
  call geom2ders(ewn, nsn, dew, dns, usrf, stagthck, d2usrfdew2, d2usrfdns2)
  call geom2ders(ewn, nsn, dew, dns, thck, stagthck, d2thckdew2, d2thckdns2)

  ! *sfp** geometric (2nd) cross-deriv. for generic input variable 'ipvr', output as 'opvr'
  call geom2derscros(dew, dns, thck, stagthck, d2thckdewdns)
  call geom2derscros(dew, dns, usrf, stagthck, d2usrfdewdns)

  ! *sfp** make a 2d array identifying if the associated point has zero thickness,
  !      has non-zero thickness and is interior, or has non-zero thickness
  !      and is along a boundary

  !*sfp* This subroutine has been altered from its original form (was a function, still included
  ! below w/ subroutine but commented out) to allow for a tweak to the CISM calculated mask (adds
  ! in an unique number for ANY arbritray boundary, be it land, water, or simply at the edge of
  ! the calculation domain). 
  !
  ! As of late July 2009, call to this function should no longer be necessary, as the mask and 
  ! code here have been altered so that the general mask can be used for flagging the appropriate
  ! boundary conditions.
  ! call maskvelostr(ewn, nsn, thck, stagthck, umask)

  allocate(uindx(ewn-1,nsn-1))

  ! *sfp** if a point from the 2d array 'mask' is associated with non-zero ice thickness,
  !      either a boundary or interior point, give it a unique number. If not, give it a zero			 
  uindx = indxvelostr(ewn, nsn, upn,  &
                      umask,pcgsize(1))

  allocate(tvel(upn,ewn-1,nsn-1)) 
  tvel = 0.0_dp
 
  ! *sfp** allocate space for variables used by 'mindcrash' function
  allocate(corr(upn,ewn-1,nsn-1,2,2),usav(upn,ewn-1,nsn-1,2))

  ! *sfp** an initial guess at the size of the sparse matrix
  pcgsize(2) = pcgsize(1) * 20

!==============================================================================
! RN_20100129: Option to load Trilinos matrix directly bypassing sparse_easy_solve
!==============================================================================

  if (whatsparse == STANDALONE_TRILINOS_SOLVER) then
     ! AGS: Get partition -- later this will be known by distributed glimmer
!      call dopartition(pcgsize(1), mySize) ! JCC - No Trilinos support yet
     allocate(myIndices(mySize))
!      call getpartition(mySize, myIndices) ! JCC - No Trilinos support yet

     ! Now send this partition to Trilinos initialization routines
!      call inittrilinos(25, mySize, myIndices) ! JCC - No Trilinos support yet

     ! Triad sparse matrix not used in this case, so save on memory
     pcgsize(2) = 1

     deallocate(myIndices)
  endif

!==============================================================================
! RN_20100126: End of the block
!==============================================================================

  ! *sfp** allocate space matrix variables
  allocate (pcgrow(pcgsize(2)),pcgcol(pcgsize(2)),rhsd(pcgsize(1)), &
            pcgval(pcgsize(2)))
  allocate(matrixA%row(pcgsize(2)), matrixA%col(pcgsize(2)), &
            matrixA%val(pcgsize(2)), answer(pcgsize(1)))
  allocate(matrixC%row(pcgsize(2)), matrixC%col(pcgsize(2)), &
            matrixC%val(pcgsize(2)))
  allocate(matrixtp%row(pcgsize(2)), matrixtp%col(pcgsize(2)), &
            matrixtp%val(pcgsize(2)))

  !*sfp* allocation for storage of (block matrix) off diagonal terms in coeff sparse matrix
  ! (these terms are usually sent to the RHS and treated as a source term in the operator splitting
  ! done in the standard Picard iteration)
  allocate (pcgrowuv(pcgsize(2)),pcgcoluv(pcgsize(2)),pcgvaluv(pcgsize(2)))
  allocate (pcgrowvu(pcgsize(2)),pcgcolvu(pcgsize(2)),pcgvalvu(pcgsize(2)))
  allocate(matrixAuv%row(pcgsize(2)),matrixAuv%col(pcgsize(2)),matrixAuv%val(pcgsize(2)))
  allocate(matrixAvu%row(pcgsize(2)),matrixAvu%col(pcgsize(2)),matrixAvu%val(pcgsize(2)))

  allocate( uk_1(pcgsize(1)), vk_1(pcgsize(1)),g_flag(pcgsize(1)) )
  allocate( vectp(pcgsize(1)), uk_1_plus(pcgsize(1)), vk_1_plus(pcgsize(1)))
  allocate( F(2*pcgsize(1)), F_plus(2*pcgsize(1)), dx(2*pcgsize(1)))
  allocate( wk1(2*pcgsize(1)), wk2(2*pcgsize(1)), rhs(2*pcgsize(1)))
  allocate( vv(2*pcgsize(1),img1), wk(2*pcgsize(1), img))

  call ghost_preprocess( ewn, nsn, upn, uindx, ughost, vghost, &
                         uk_1, vk_1, uvel, vvel, g_flag) ! jfl_20100430

  print *, ' '
  print *, 'Running Payne/Price higher-order dynamics with JFNK solver' 

!==============================================================================
! Beginning of Newton loop. Solves F(x) = 0 for x where x = [v, u] and
!                                                       F = [Fv(u,v), Fu(u,v)] 
!==============================================================================

  do k = 1, kmax

!==============================================================================
! calculate F(u^k-1,v^k-1)
!==============================================================================

    calcoffdiag = .true.    ! save off diag matrix components

    call calc_F (ewn, nsn, upn, stagsigma, k,                    &
                 whichefvs, efvs, uvel, vvel,                    &
                 dew, dns, sigma, thck,                          &
                 dusrfdew, dthckdew, d2usrfdew2, d2thckdew2,     &
                 dusrfdns, dthckdns, d2usrfdns2, d2thckdns2,     &
                 d2usrfdewdns, d2thckdewdns, dlsrfdew, dlsrfdns, &
                 stagthck, whichbabc,            &
                 uindx, umask,                                   &
                 lsrf, topg, minTauf, flwa, beta, btraction,     &
                 vk_1, uk_1, pcgsize(1), 2*pcgsize(1), g_flag,   &
                 F, L2norm, matrixA, matrixC )

    calcoffdiag = .false.    ! next time calling calc_F, DO NOT save off diag matrix components

    L2norm_wig = sqrt(DOT_PRODUCT(F,F)) ! with ghost

!==============================================================================
! -define nonlinear target (if k=1)
! -check at all k if target is reached
!==============================================================================

!      print *, 'target with L2norm with or without ghost?'

    if (k .eq. 1) NL_target = NL_tol * L2norm

    print *, 'L2 with, without ghost (k)= ', k, L2norm_wig, L2norm

    if (L2norm .lt. NL_target) exit ! nonlinear convergence criterion

!==============================================================================
! solve J(u^k-1,v^k-1)dx = -F(u^k-1,v^k-1) with fgmres, dx = [dv, du]  
!==============================================================================

      rhs = -1d0*F

      dx  = 0d0 ! initial guess

      call forcing_term (k, L2norm_wig, gamma_l)

      tol = gamma_l * L2norm_wig ! setting the tolerance for fgmres

      epsilon = 1d-07 ! for J*vector approximation

      maxiteGMRES = 300
      
      iout   = 0    ! set  higher than 0 to have res(ite)

      icode = 0

 10   CONTINUE
      
!       call fgmres (2*pcgsize(1),img,rhs,dx,itenb,vv,wk,wk1,wk2, &
!                    tol,maxiteGMRES,iout,icode,tot_its) ! JCC - No Trilinos support yet

      IF ( icode == 1 ) THEN   ! precond step: use of Picard linear solver
                               ! wk2 = P^-1*wk1

      call apply_precond ( matrixA, matrixC, pcgsize(1), 2*pcgsize(1), &
                           wk1, wk2, whichsparse ) 

      GOTO 10

      ELSEIF ( icode >= 2 ) THEN  ! matvec step: Jacobian free approach
                                  ! J*wk1 ~ wk2 = (F_plus - F)/epsilon
         
! form  v^k-1_plus = v^k-1 + epsilon*wk1v. We use solver_postprocess to 
! transform vk_1_plus from a vector to a 3D field. (same idea for u^k-1_plus)

         vectp(:) = wk1(1:pcgsize(1)) ! for v
         vk_1_plus = vk_1 + epsilon*vectp

         call solver_postprocess( ewn, nsn, upn, 2, uindx, &
                                  vk_1_plus, vvel, ghostbvel )

         vectp(:) = wk1(pcgsize(1)+1:2*pcgsize(1)) ! for u
         uk_1_plus = uk_1 + epsilon*vectp

         call solver_postprocess( ewn, nsn, upn, 1, uindx, &
                                  uk_1_plus, uvel, ghostbvel )

! form F(x + epsilon*wk1) = F(u^k-1 + epsilon*wk1u, v^k-1 + epsilon*wk1v)

    call calc_F (ewn, nsn, upn, stagsigma, k,                    &
                 whichefvs, efvs, uvel, vvel,                    &
                 dew, dns, sigma, thck,                          &
                 dusrfdew, dthckdew, d2usrfdew2, d2thckdew2,     &
                 dusrfdns, dthckdns, d2usrfdns2, d2thckdns2,     &
                 d2usrfdewdns, d2thckdewdns, dlsrfdew, dlsrfdns, &
                 stagthck, whichbabc,            &
                 uindx, umask,                                   &
                 lsrf, topg, minTauf, flwa, beta, btraction,     &
                 vk_1_plus, uk_1_plus, pcgsize(1), 2*pcgsize(1), g_flag, &
                 F_plus, crap, matrixtp, matrixtp )

! put approximation of J*wk1 in wk2

          wk2 =  ( F_plus - F ) / epsilon

         GOTO 10
      ENDIF

!------------------------------------------------------------------------
! End of FGMRES method    
!------------------------------------------------------------------------

      if (tot_its .eq. maxiteGMRES) then
         print *,'WARNING: FGMRES has not converged'
         stop
      endif

! icode = 0 means that fgmres has finished and sol contains the app. solution

!------------------------------------------------------------------------
! Update solution vectors (x^k = x^k-1 + dx) and 3D fields 
!------------------------------------------------------------------------

      vk_1 = vk_1 + dx(1:pcgsize(1)) 
      uk_1 = uk_1 + dx(pcgsize(1)+1:2*pcgsize(1))
      ! in fact vk and uk but we use vk_1 and uk_1 to save memory

      call solver_postprocess( ewn, nsn, upn, 2, uindx, vk_1, vvel, ghostbvel )
      call solver_postprocess( ewn, nsn, upn, 1, uindx, uk_1, uvel, ghostbvel )

! WATCHOUT FOR PERIODIC BC      

  end do

  call ghost_postprocess( ewn, nsn, upn, uindx, uk_1, vk_1, &
                          ughost, vghost )

  do ns = 1,nsn-1
      do ew = 1,ewn-1 
      ! calc. fluxes from converged vel. fields (needed for input to thickness evolution subroutine)
         if (umask(ew,ns) > 0) then
             uflx(ew,ns) = vertintg(upn, sigma, uvel(:,ew,ns)) * stagthck(ew,ns)
             vflx(ew,ns) = vertintg(upn, sigma, vvel(:,ew,ns)) * stagthck(ew,ns)
         end if
      end do
  end do

  ! *sfp* de-allocation of sparse matrix solution variables 
  deallocate(tvel)
  deallocate(uindx,corr,usav)
  deallocate(pcgval,pcgrow,pcgcol,rhsd)
  deallocate(pcgvaluv,pcgrowuv,pcgcoluv)
  deallocate(pcgvalvu,pcgrowvu,pcgcolvu)
  deallocate(matrixA%row, matrixA%col, matrixA%val)
  deallocate(matrixAuv%row, matrixAuv%col, matrixAuv%val)
  deallocate(matrixAvu%row, matrixAvu%col, matrixAvu%val)
  deallocate(matrixC%row, matrixC%col, matrixC%val)
  deallocate(matrixtp%row, matrixtp%col, matrixtp%val)
  deallocate(uk_1, vk_1, g_flag)
  deallocate(answer, dx, vectp, uk_1_plus, vk_1_plus )
  deallocate(F, F_plus)
  deallocate(wk1, wk2)
  deallocate(vv, wk)

  return

end subroutine JFNK

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
                       stagsigma, counter,   &
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

!whl - If temp and flwa live on the staggered vertical grid (like the effective viscosity),
!      then the size of flwa is (upn-1), and vertical averaging of flwa is not needed here.

     if (size(flwa,1)==upn-1) then   ! temperature and flwa live on staggered vertical grid

        do ns = 2,nsn-1
        do ew = 2,ewn-1
           if (thck(ew,ns) > 0.0_dp) then
              ! This is the rate factor term in the expression for the eff. visc: 1/2*A^(-1/n).
              ! If both temperature and eff. visc. live on a staggered grid in the vertical, then
              !  no vertical averaging is needed.
              flwafact(1:upn-1,ew,ns) = 0.5_dp * flwa(1:upn-1,ew,ns)**p1
           end if
        end do
        end do

     else  ! size(flwa,1)=upn; temperature and flwa live on unstaggered vertical grid

       do ns = 2,nsn-1
       do ew = 2,ewn-1
          if (thck(ew,ns) > 0.0_dp) then
             ! this is the rate factor term in the expression for the eff. visc: 1/2*A^(-1/n),
             ! which is averaged to midpoints in the vertical (i.e. it lives on a staggered 
             ! grid in the vertical, which is the case for "efvs" as well).
             forall (up = 1:upn-1) flwafact(up,ew,ns) = 0.5_dp * (sum(flwa(up:up+1,ew,ns)) / 2.0_dp)**p1
          end if
       end do
       end do

     end if   ! present(flwa_vstag)

  endif       ! counter

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
!                         f1 * (ugradup**2 + vgradup**2)      ! make line ACTIVE for "capping" version (see note below)   
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
!            efvs(:,ew,ns) = flwafact(:,ew,ns) * effstr**p2

        else
           efvs(:,ew,ns) = effstrminsq ! if the point is associated w/ no ice, set to min value
        end if

       end do   ! end ew
   end do       ! end ns

  case(1)       ! set the eff visc to some const value 

!   *sfp* changed default setting for linear viscosity so that the value of the rate
!   factor is taken into account
!    efvs = 1.0_dp
  do ns = 2,nsn-1
      do ew = 2,ewn-1
       if (thck(ew,ns) > 0.0_dp) then
        efvs(1:upn-1,ew,ns) = 0.5_dp * flwa(1:upn-1,ew,ns)**(-1.0_dp)
        else
           efvs(:,ew,ns) = effstrminsq ! if the point is associated w/ no ice, set to min value
       end if
      end do
  end do

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

subroutine solver_postprocess( ewn, nsn, upn, pt, uindx, answrapped, ansunwrapped, ghostbvel )

  ! Unwrap the vels from the solution vector and place into a 3d array.

  implicit none

  integer, intent(in) :: ewn, nsn, upn, pt
  integer, dimension(:,:), intent(in) :: uindx
  real (kind = dp), dimension(:), intent(in) :: answrapped
  real (kind = dp), dimension(upn,ewn-1,nsn-1), intent(out) :: ansunwrapped
  real (kind = dp), dimension(:,:,:,:), intent(inout) :: ghostbvel

  integer, dimension(2) :: loc
  integer :: ew, ns

  do ns = 1,nsn-1
      do ew = 1,ewn-1
          if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            ansunwrapped(:,ew,ns) = answrapped(loc(1):loc(2))
            !! save the fictitious basal velocities for basal traction calculation !!
            ghostbvel(pt,:,ew,ns) = answrapped( loc(2)-1:loc(2)+1 )  
          else
            ansunwrapped(:,ew,ns) = 0.0d0
          end if
      end do
  end do

end subroutine solver_postprocess

!***********************************************************************

subroutine form_matrix( matrix ) ! for JFNK solver

  ! Puts sparse matrix variables in SLAP triad format into "matrix" 
  ! derived type. Similar to solver_preprocess but does not form answer vector

  implicit none

!  integer, intent(in) :: ewn, nsn, upn
  type(sparse_matrix_type), intent(inout) :: matrix

  pcgsize(2) = ct - 1

  matrix%order = pcgsize(1)
  matrix%nonzeros = pcgsize(2)
  matrix%symmetric = .false.

  matrix%row = pcgrow
  matrix%col = pcgcol
  matrix%val = pcgval

end subroutine form_matrix

!***********************************************************************

subroutine forcing_term ( k, L2normk_1, gamma_l ) 

  ! Calculates the forcing term (i.e. the factor that multiplies the initial
  ! L2 norm to determine the tolerance for the linear solve in the JFNK solver)
  ! at iteration k given the L2norm at k-1 and k-2.
  ! jfl, 10 Sept 2010

  ! See eq 2.6 in S.C. Eisenstat, H.F. Walker, Choosing the forcing terms in
  ! an inexact Newton method, SIAM J. Sci. Comput. 17 (1996) 16-32.

  implicit none
      
  integer, intent(in) :: k
  real (kind = dp), intent(in) :: L2normk_1 ! L2 norm at k-1
  real (kind = dp), intent(out):: gamma_l
  real (kind = dp) :: gamma_ini, gamma_min, expo
  real (kind = dp), save :: L2normk_2      ! L2 norm at k-2

      gamma_ini = 0.9d0
      gamma_min = 0.01d0
      expo      = 2d0

      if (k .eq. 1) then
         gamma_l = gamma_ini
      else
         gamma_l = (L2normk_1 / L2normk_2)**expo
      endif

      if (gamma_l .gt. gamma_ini) gamma_l = gamma_ini
      if (gamma_l .lt. gamma_min) gamma_l = gamma_min
      
      L2normk_2 = L2normk_1

end subroutine forcing_term

!***********************************************************************

subroutine apply_precond( matrixA, matrixC, nu1, nu2, wk1, wk2, whichsparse ) 

  ! Apply preconditioner operator for JFNK solver: wk2 = P^-1 *wk1 
  ! The preconditioner operator is in fact taken from the Picard solver
  ! There is a splitting of the v (A matrix) and u (C matrix) equations
  ! Each component is solved to a loose tolerance (as opposed to Picard)

  implicit none

  integer, intent(in) :: nu1, nu2, whichsparse
  integer :: iter
  type(sparse_matrix_type), intent(in) :: matrixA, matrixC
  real (kind = dp), dimension(nu2), intent(in) :: wk1
  real (kind = dp), dimension(nu2), intent(out):: wk2
  real (kind = dp), dimension(nu1) :: answer, vectp
  real (kind = dp) :: err

! precondition v component 
       
      answer = 0d0 ! initial guess
      vectp(:) = wk1(1:nu1) ! rhs for precond v
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
         call sparse_easy_solve(matrixA, vectp, answer, err, iter, whichsparse, nonlinear_solver = nonlinear)
      else
!          call restoretrilinosmatrix(0); ! JCC - No Trilinos support yet
!          call solvewithtrilinos(vectp, answer, linearSolveTime) ! JCC - No Trilinos support yet
         totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
!         write(*,*) 'Total linear solve time so far', totalLinearSolveTime
      endif
      wk2(1:nu1) = answer(:)

! precondition u component 
       
      answer = 0d0 ! initial guess
      vectp(:) = wk1(nu1+1:nu2) ! rhs for precond u
      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
         call sparse_easy_solve(matrixC, vectp, answer, err, iter, whichsparse, nonlinear_solver = nonlinear)
      else
!          call restoretrilinosmatrix(1); ! JCC - No Trilinos support yet
!          call solvewithtrilinos(vectp, answer, linearSolveTime) ! JCC - No Trilinos support yet
         totalLinearSolveTime = totalLinearSolveTime + linearSolveTime
!         write(*,*) 'Total linear solve time so far', totalLinearSolveTime
      endif
      wk2(nu1+1:nu2) = answer(:)

end subroutine apply_precond

!***********************************************************************

subroutine calc_F (ewn, nsn, upn, stagsigma, counter,            &
                 whichefvs, efvs, uvel, vvel,                    &
                 dew, dns, sigma, thck,                          &
                 dusrfdew, dthckdew, d2usrfdew2, d2thckdew2,     &
                 dusrfdns, dthckdns, d2usrfdns2, d2thckdns2,     &
                 d2usrfdewdns, d2thckdewdns, dlsrfdew, dlsrfdns, &
                 stagthck, whichbabc,                            &
                 uindx, umask,                                   &
                 lsrf, topg, minTauf, flwa, beta, btraction,     &
                 vtp, utp, nu1, nu2, g_flag,                     &
                 F, L2norm, matrixA, matrixC)

  ! Calculates either F(x) or F(x+epsilon*vect) for the JFNK method
  ! Recall that x=[v,u]
  ! vtp is either v^k-1 or v^k-1+epsilon*vect_v (same idea for utp)

  implicit none

  integer, intent(in) :: ewn, nsn, upn, counter, whichbabc, whichefvs
  integer, intent(in) :: nu1, nu2
  integer, dimension(nu1), intent(in) :: g_flag ! 0 :reg cell
                                                ! 1 :top ghost, 2 :base ghost
  integer, dimension(:,:), intent(in) :: uindx, umask

  type(sparse_matrix_type), intent(inout) :: matrixA, matrixC

  real (kind = dp) :: L2square
  real (kind = dp), intent(out):: L2norm
  real (kind = dp), dimension(nu1) :: vectp
  real (kind = dp), dimension(nu1), intent(in) :: vtp, utp
  real (kind = dp), dimension(nu2), intent(out) :: F
  real (kind = dp), intent(in) :: dew, dns    
  real (kind = dp), dimension(:), intent(in) :: sigma, stagsigma

  real (kind = dp), dimension(:,:), intent(in) :: thck,               &
                                                  dusrfdew, dthckdew, & 
                                                  d2usrfdew2, d2thckdew2, &
                                                  dusrfdns, dthckdns, &
                                                  d2usrfdns2, d2thckdns2, &
                                                  d2usrfdewdns, d2thckdewdns,&
                                                  dlsrfdew, dlsrfdns, &
                                                  stagthck, lsrf, topg, &
                                                  minTauf, beta

  real (kind = dp), dimension(:,:,:), intent(inout) :: efvs, btraction
  real (kind = dp), dimension(:,:,:), intent(in) :: uvel, vvel, flwa
           
    call findefvsstr(ewn,  nsn,  upn,       &
                     stagsigma,  counter,  &
                     whichefvs,  efvs,     &
                     uvel,       vvel,     &
                     flwa,       thck,     &
                     dusrfdew,   dthckdew, &
                     dusrfdns,   dthckdns, &
                     umask)

!==============================================================================
! jfl 20100412: residual for v comp: Fv= A(utp,vtp)vtp - b(utp,vtp)  
!==============================================================================

    ! *sfp** calculation of coeff. for stress balance calc. 
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
                     beta, btraction,             &
                     counter, 0 )

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
      call form_matrix ( matrixA ) ! to get A(utp,vtp)
    else
!       call savetrilinosmatrix(0); ! JCC - No Trilinos support yet
    end if
    
      vectp = vtp

!     call res_vect(matrixA, vectp, rhsd, nu1, counter, g_flag, L2square, whatsparse)!rhsd = b ! JCC - No Trilinos support yet
    L2norm = L2square

    F(1:nu1) = vectp ! Fv
      
!==============================================================================
! jfl 20100412: residual for u comp: Fu= C(utp,vtp)utp - d(utp,vtp)  
!==============================================================================

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
                     beta, btraction,             &
                     counter, 0 )


    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
      call form_matrix ( matrixC ) ! to get C(utp,vtp)
    else
!       call savetrilinosmatrix(1); ! JCC - No Trilinos support yet
    end if
    
      vectp = utp

!     call res_vect(matrixC, vectp, rhsd, nu1, counter, g_flag, L2square, whatsparse)!rhsd = d ! JCC - No Trilinos support yet
    L2norm = sqrt(L2norm + L2square)

    F(nu1+1:nu2) = vectp ! Fu

end subroutine calc_F

!***********************************************************************

subroutine ghost_preprocess( ewn, nsn, upn, uindx, ughost, vghost, & 
                             uk_1, vk_1, uvel, vvel, g_flag)

! puts vel values in  uk_1, vk_1 (including ghost values) and creates the
! ghost flag vector. uk_1, vk_1 and the ghost flag vector are used for 
! the residual calculation (jfl 20100430)

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: uindx
  integer, dimension(:), intent(out) :: g_flag 
  real (kind = dp), dimension(2,ewn-1,nsn-1), intent(in) ::ughost,vghost 
  real (kind = dp), dimension(:,:,:), intent(in) :: uvel, vvel
  real (kind = dp), dimension(:), intent(out) :: uk_1, vk_1 

  integer :: ew, ns
  integer, dimension(2) :: loc

  g_flag = 0

  do ns = 1,nsn-1
   do ew = 1,ewn-1
        if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            uk_1(loc(1):loc(2)) = uvel(:,ew,ns)
            uk_1(loc(1)-1)      = ughost(1,ew,ns) ! ghost at top
            uk_1(loc(2)+1)      = ughost(2,ew,ns) ! ghost at base

            vk_1(loc(1):loc(2)) = vvel(:,ew,ns)
            vk_1(loc(1)-1)      = vghost(1,ew,ns) ! ghost at top
            vk_1(loc(2)+1)      = vghost(2,ew,ns) ! ghost at base

            g_flag(loc(1)-1) = 1 ! ghost at top
            g_flag(loc(2)+1) = 2 ! ghost at base
        end if
    end do
  end do

end subroutine ghost_preprocess

!***********************************************************************

subroutine ghost_postprocess( ewn, nsn, upn, uindx, uk_1, vk_1, &
                              ughost, vghost )

! puts ghost values (which are now in  uk_1 and vk_1) into ughost and 
! vghost so that they can be used fro the next time step (jfl 20100430)

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, dimension(:,:), intent(in) :: uindx
  real (kind = dp), dimension(:), intent(in) :: uk_1, vk_1
  real (kind = dp), dimension(2,ewn-1,nsn-1), intent(out) :: ughost,vghost

  integer :: ew, ns
  integer, dimension(2) :: loc

  do ns = 1,nsn-1
      do ew = 1,ewn-1
          if (uindx(ew,ns) /= 0) then
            loc = getlocrange(upn, uindx(ew,ns))
            ughost(1,ew,ns) = uk_1(loc(1)-1) ! ghost at top
            ughost(2,ew,ns) = uk_1(loc(2)+1) ! ghost at base
            vghost(1,ew,ns) = vk_1(loc(1)-1) ! ghost at top
            vghost(2,ew,ns) = vk_1(loc(2)+1) ! ghost at base
          else 
            ughost(1,ew,ns) = 0d0
            ughost(2,ew,ns) = 0d0
            vghost(1,ew,ns) = 0d0
            vghost(2,ew,ns) = 0d0
          end if
      end do
  end do

end subroutine ghost_postprocess

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

!      mindcrshstr = vel; ! jfl uncomment this and comment out line above 
!                         ! to avoid the unstable manifold correction

    elsewhere

      mindcrshstr = vel;

    end where

  else

    mindcrshstr = vel;

  end if

  !*sfp* Old version
  if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 1; new(pt) = 2; end if  

  !*sfp* correction from Carl Gladdish
  !if (new(pt) == 1) then; old(pt) = 1; new(pt) = 2; else; old(pt) = 2; new(pt) = 1; end if   

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
    !print '("* ",i3,g20.6,3i6,g20.6)', counter, resid, locat, vel(locat(1),locat(2),locat(3))*vel0

  return

end function mindcrshstr

!***********************************************************************

function mindcrshstr2(pt,whichresid,vel,counter,resid)

  ! Function to perform 'unstable manifold correction' (see Hindmarsch and Payne, 1996,
  ! "Time-step limits for stable solutions of the ice-sheet equation", Annals of 
  ! Glaciology, 23, p.74-85)

  ! Alternate unstable manifold scheme, based on DeSmedt, Pattyn, and De Goen, J. Glaciology 2010
  ! Written by Carl Gladdish

  implicit none

  real (kind = dp), intent(in), dimension(:,:,:) :: vel
  integer, intent(in) :: counter, pt, whichresid 
  real (kind = dp), intent(out) :: resid

  real (kind = dp), dimension(size(vel,1),size(vel,2),size(vel,3)) :: mindcrshstr2

  integer, parameter :: start_umc = 3
  real (kind=dp), parameter :: cvg_accel = 2.0_dp
  real (kind = dp), parameter :: small = 1.0e-16_dp

  real (kind=dp) in_prod, len_new, len_old, mean_rel_diff, sig_rel_diff
  real (kind=dp) :: theta
  real (kind = dp), intrinsic :: abs, acos
  
  integer, dimension(2), save :: new = 1, old = 2
  integer :: locat(3)
  
  integer :: nr
  integer,      dimension(size(vel,1),size(vel,2),size(vel,3)) :: vel_ne_0
  real(kind=dp),dimension(size(vel,1),size(vel,2),size(vel,3)) :: rel_diff

  if (counter == 1) then
    usav(:,:,:,pt) = 0.0d0
    corr(:,:,:,:,:) = 0.0d0
  end if

  corr(:,:,:,new(pt),pt) = vel - usav(:,:,:,pt)           

  if (counter >= start_umc) then
  
  in_prod = sum( corr(:,:,:,new(pt),pt) * corr(:,:,:,old(pt),pt) )
  len_new = sqrt(sum( corr(:,:,:,new(pt),pt) * corr(:,:,:,new(pt),pt) ))
  len_old = sqrt(sum( corr(:,:,:,old(pt),pt) * corr(:,:,:,old(pt),pt) ))
 
  theta = acos( in_prod / (len_new * len_old + small) )
    
   if (theta  < (1.0/8.0)*pi) then
        mindcrshstr2 = usav(:,:,:,pt) + cvg_accel * corr(:,:,:,new(pt),pt)
!        print *, theta/pi, 'increased correction'
   else if(theta < (19.0/20.0)*pi) then
        mindcrshstr2 = vel
!        print *, theta/pi, 'standard correction'
   else 
        mindcrshstr2 = usav(:,:,:,pt) + (1.0/cvg_accel) * corr(:,:,:,new(pt),pt)
!        print *, theta/pi, 'decreasing correction'
   end if

  else 

    mindcrshstr2 = vel;
 !   print *, 'Not attempting adjustment to correction'  
   
  end if


  ! now swap slots for storing the previous correction
  if (new(pt) == 1) then
      old(pt) = 1; new(pt) = 2
  else
      old(pt) = 2; new(pt) = 1
  end if

  if (counter == 1) then
        usav_avg = 1.0_dp
  else
        usav_avg(1) = sum( abs(usav(:,:,:,1)) ) / size(vel)  ! a x-dir transport velocity scale
        usav_avg(2) = sum( abs(usav(:,:,:,2)) ) / size(vel)  ! a y-dir transport velocity scale
  end if

!  print *, 'usav_avg(1)',usav_avg(1),'usav_avg(2)',usav_avg(2)

  select case (whichresid)

  ! options for residual calculation method, as specified in configuration file 
  ! (see additional notes in "higher-order options" section of documentation)
  ! case(0): use max of abs( vel_old - vel ) / vel ) 
  ! case(1): use max of abs( vel_old - vel ) / vel ) but ignore basal vels 
  ! case(2): use mean of abs( vel_old - vel ) / vel )

   case(0)
    rel_diff = 0.0_dp
    vel_ne_0 = 0
    where ( mindcrshstr2 .ne. 0.0_dp )
        vel_ne_0 = 1
        rel_diff = abs((usav(:,:,:,pt) - mindcrshstr2) / mindcrshstr2) & 
                           * usav_avg(pt)/sqrt(sum(usav_avg ** 2.0))
    end where

    resid = maxval( rel_diff, MASK = mindcrshstr2 .ne. 0.0_dp )
    locat = maxloc( rel_diff, MASK = mindcrshstr2 .ne. 0.0_dp )

!    mean_rel_diff = sum(rel_diff) / sum(vel_ne_0)
!    sig_rel_diff = sqrt( sum((rel_diff - mean_rel_diff) ** 2.0 )/ sum(vel_ne_0) )
!    print *, 'mean', mean_rel_diff, 'sig', sig_rel_diff

    !write(*,*) 'locat', locat
    !call write_xls('resid1.txt',abs((usav(1,:,:,pt) - mindcrshstr2(1,:,:)) / (mindcrshstr2(1,:,:) + 1e-20)))

   case(1)
    !**cvg*** should replace vel by mindcrshstr2 in the following lines, I belive
    nr = size( vel, dim=1 ) ! number of grid points in vertical ...
    resid = maxval( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
                        MASK = vel .ne. 0.0_dp)
    locat = maxloc( abs((usav(1:nr-1,:,:,pt) - vel(1:nr-1,:,:) ) / vel(1:nr-1,:,:) ),  &
            MASK = vel .ne. 0.0_dp)

   case(2)
    !**cvg*** should replace vel by mindcrshstr2 in the following lines, I belive
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

  usav(:,:,:,pt) = mindcrshstr2

    ! Additional debugging line, useful when trying to determine if convergence is being consistently 
    ! held up by the residual at one or a few particular locations in the domain.
!    print '("* ",i3,g20.6,3i6,g20.6)', counter, resid, locat, vel(locat(1),locat(2),locat(3))*vel0

  return

end function mindcrshstr2

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
                       beta, btraction,             &
                       count, assembly )

  ! Main subroutine for determining coefficients that go into the LHS matrix A 
  ! in the expression Au = b. Calls numerous other subroutines, including boundary
  ! condition subroutines, which determin "b".

  implicit none

  integer, intent(in) :: ewn, nsn, upn, count, assembly
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
  real (kind = dp), dimension(:,:,:), intent(inout) :: btraction

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
  integer :: ew, ns, up, up_start

  ct = 1        ! index to count the number of non-zero entries in the sparse matrix
  ct2 = 1

  if( assembly == 1 )then   ! for normal assembly (assembly=0), start vert index at sfc and go to bed
      up_start = upn        ! for boundary traction calc (assembly=1), do matrix assembly on for equations at bed
  else
      up_start = 1
  end if

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
        do up = up_start, upn

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
                         betasquared(ew,ns),           &
                         btraction,                    &
                         whichbabc, assembly )

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

        do up = up_start, upn
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
                         betasquared(ew,ns),           &
                         btraction,                    &
                         whichbabc, assembly,          &              
                         abar=flwabar, cc=count )
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
        do up = up_start, upn
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
                   btraction,              &
                   whichbabc, assembly,    &
                   local_thisvel,          &
                   abar, cc)

  ! This subroutine does the bulk of the work in calling the appropriate discretiztion routines,
  ! which determine the values for coefficients that will go into the sparse matrix, for points
  ! on and inside of the boundaries.

  implicit none

  integer, intent(in) :: ewn, nsn, upn
  integer, intent(in) :: ew, ns, up
  real (kind = dp), intent(in) :: dew, dns
  integer, intent(in) :: pt, whichbabc, assembly
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
  real (kind = dp), dimension(:,:,:), intent(inout) :: btraction
  real (kind = dp), intent(in), optional :: local_thisvel
  real (kind = dp), intent(in), optional :: abar
  integer, intent(in), optional :: cc

  ! storage space for coefficients that go w/ the discretization at the local point up, ew, ns.
  ! Note that terms other than 'g' are used for storing particular parts needed for calculation
  ! of the basal traction vector.
  real (kind = dp), dimension(3,3,3) :: g, h, g_cros, g_vert, g_norm, g_vel_lhs, g_vel_rhs

  ! source term for the rhs when using ice shelf lateral boundary condition,
  ! e.g. source = rho*g*H/(2*Neff) * ( 1 - rho_i / rho_w ) for ice shelf
  real (kind = dp) :: source

  real (kind = dp) :: slopex, slopey    ! local sfc (or bed) slope terms

  ! lateral boundary normal and vector to indicate use of forward
  ! or bacward one-sided diff. when including specified stress lateral bcs
  real (kind = dp), dimension(2) :: fwdorbwd, normal

  real (kind = dp) :: nz   ! z dir normal vector component at sfc or bed (takes diff value for each)

  integer, dimension(2) :: bcflag  ! indicates choice of sfc and basal bcs ...

  real (kind = dp) :: scalebabc

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

        if( up == 1 )then                ! specify necessary variables and flags for free sfc
           bcflag = (/1,0/)
           locplusup = loc(1) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.0_dp
        else                             ! specify necessary variables and flags for basal bc
   
           if( whichbabc == 6 )then
                bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...
                                             ! better to specify very large value for betasquared below
           elseif( whichbabc >=0 .and. whichbabc <= 5 )then
                bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                              ! where betasquared is MacAyeal-type traction parameter
           elseif( whichbabc == 7 )then
                bcflag = (/1,2/)              ! plastic bed iteration using Newton implementation
           end if   

           locplusup = loc(1) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
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

        ! add on coeff. associated w/ du/dsigma  
        g(:,3,3) = g(:,3,3) &
                 + vertimainbc( stagthck(ew,ns), bcflag, dup(up),     &
                                local_efvs,      betasquared,   g_vert,    nz, &
                                plastic_coeff=plastic_coeff_lhs(pt,ew,ns)  )

        !! scale basal bc coeffs when using JFNK solver 
        scalebabc = scalebasalbc( g, bcflag, lateralboundry, betasquared, local_efvs )
        g = g / scalebabc

        ! put the coeff. for the b.c. equation in the same place as the prev. equation
        ! (w.r.t. cols), on a new row ...
        call fillsprsebndy( g, locplusup, loc_latbc, up, normal, pt )

        ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
        ! which results from moving them from the LHS over to the RHS, has been moved
        ! inside of "croshorizmainbc_lat".
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
                                                 * local_othervel ) / scalebabc

        if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
            storeoffdiag = .true.
            h = -croshorizmainbc_lat(dew,           dns,           & 
                                slopex,        slopey,        &
                                dsigmadew(up), dsigmadns(up), & 
                                pt,            2,             & 
                                dup(up),       local_othervel,& 
                                oneortwo,      twoorone,      & 
                                onesideddiff,                 &
                                 normal,fwdorbwd) / scalebabc 
            call fillsprsebndy( h, locplusup, loc_latbc, up, normal, pt ) 
            storeoffdiag = .false.
        end if     

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

    if( cc < 5 )then   ! use this to "pre-condition" the shelf BC w/ the simple, 1d version
!    if( cc >= 0 )then   ! use this to use only the 1d version
!    if( cc > 1000000 )then   ! use this to go straight to the full 2d version of the bc

    ! --------------------------------------------------------------------------------------
    ! (1) source term (strain rate at shelf/ocean boundary) from Weertman's analytical solution 
    ! --------------------------------------------------------------------------------------
    ! See eq. 2, Pattyn+, 2006, JGR v.111; eq. 8, Vieli&Payne, 2005, JGR v.110). Note that this 
    ! contains the 1d assumption that ice is not spreading lateraly !(assumes dv/dy = 0 for u along flow)
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
    source = (rhoi*grav*stagthck(ew,ns)*thk0) / tau0 / 2.0_dp * ( 1.0_dp - rhoi / rhoo )

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
    call fillsprsebndy( g, locplusup, loc_latbc, up, normal, pt )

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

     if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
         storeoffdiag = .true.
         h = -croshorizmainbc_lat(dew,           dns,            &
                                 slopex,        slopey,         &
                                 dsigmadew(up), dsigmadns(up),  &
                                 pt,            1,              &
                                 dup(up),       local_othervel, &
                                 oneortwo,      twoorone,       &
                                 onesideddiff,                  &
                                 normal,        fwdorbwd)
         call fillsprsebndy( h, locplusup, loc_latbc, up, normal, pt )
         storeoffdiag = .false.
     end if

  else   ! NOT at a lateral boundary 

! *********************************************************************************************
! normal discretization for points inside of lateral boundary and inside main body of ice sheet

     g = normhorizmain(pt,up,local_efvs)
     g(:,2,2) = g(:,2,2) + vertimain(hsum(local_efvs),up)
     call fillsprsemain(g,locplusup,loc,up,pt)
     ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
     ! which results from moving them from the LHS over to the RHS, is explicit and 
     ! hast NOT been moved inside of "croshorizmin" (as is the case for the analogous
     ! boundary condition routines).

     ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
     ! which results from moving them from the LHS over to the RHS, is explicit and 
     ! hast NOT been moved inside of "croshorizmin" (as is the case for the analogous
     ! boundary condition routines).
     rhsd(locplusup) = thisdusrfdx(ew,ns) - sum(croshorizmain(pt,up,local_efvs) * local_othervel)

     if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
         storeoffdiag = .true.
         h = croshorizmain(pt,up,local_efvs)   
         call fillsprsemain(h,locplusup,loc,up,pt)
         storeoffdiag = .false.
     end if     

  end if

! *********************************************************************************************
! higher-order sfc and bed boundary conditions in main body of ice sheet (NOT at lat. boundry)

  if(  ( up == upn  .or. up == 1 ) .and. .not. lateralboundry) then

        if( up == 1 )then                ! specify necessary variables and flags for free sfc
           bcflag = (/1,0/)
           locplusup = loc(1) + up - 1   ! reverse the sparse matrix / rhs vector row index by 1 ...
           slopex = -dusrfdew(ew,ns); slopey = -dusrfdns(ew,ns); nz = 1.0_dp
        else                             ! specify necessary variables and flags for basal bc

           if( whichbabc == 6 )then
                bcflag = (/0,0/)             ! flag for u=v=0 at bed; doesn't work well so commented out here...
                                             ! better to specify very large value for betasquared below
           elseif( whichbabc >=0 .and. whichbabc <= 5 )then
                bcflag = (/1,1/)              ! flag for specififed stress at bed: Tau_zx = betasquared * u_bed,
                                              ! where betasquared is MacAyeal-type traction parameter
           elseif( whichbabc == 7 )then
                bcflag = (/1,2/)              ! plastic bed iteration (new)
           end if
           
           locplusup = loc(1) + up + 1   ! advance the sparse matrix / rhs row vector index by 1 ...
           slopex = dlsrfdew(ew,ns); slopey = dlsrfdns(ew,ns); nz = -1.0_dp
        
        end if

     g = normhorizmainbc(dew,           dns,     &
                         slopex,        slopey,  &
                         dsigmadew(up), dsigmadns(up),  &
                         pt,            bcflag,  &
                         dup(up),                &
                         oneorfour,     fourorone)

     g_norm = g              ! save for basal traction calculation

     ! add on coeff. associated w/ du/dsigma
     g(:,2,2) = g(:,2,2)   &
              + vertimainbc( stagthck(ew,ns),bcflag,dup(up),local_efvs,betasquared, &
                            g_vert, nz, plastic_coeff=plastic_coeff_lhs(pt,ew,ns) )


     !! scale basal bc coeffs when using JFNK solver 
     scalebabc = scalebasalbc( g, bcflag, lateralboundry, betasquared, local_efvs )
     g = g / scalebabc

     ! put the coeff. for the b.c. equation in the same place as the prev. equation
     ! (w.r.t. cols), on a new row ...
     call fillsprsemain(g,locplusup,loc,up,pt)

     ! NOTE that in the following expression, the "-" sign on the crosshoriz terms, 
     ! which results from moving them from the LHS over to the RHS, has been moved
     ! inside of "croshorizmainbc".

     if( bcflag(2) == 2 )then

          rhsd(locplusup) = sum( croshorizmainbc(dew,           dns,            &
                                             slopex,        slopey,         &
                                             dsigmadew(up), dsigmadns(up),  &
                                             pt,            bcflag,         &
                                             dup(up),       local_othervel, &
                                             local_efvs,                    &
                                             oneortwo,      twoorone,       &
                                             g_cros,                       &
                                             velbc=velbcvect(pt,ew,ns),     &
                                             plastic_coeff=plastic_coeff_rhs(pt,ew,ns) ) &
                                              * local_othervel )            &
                                             + plastic_rhs(pt,ew,ns) / (sum(local_efvs(2,:,:))/4.0d0) &
                                             *(len0/thk0) / scalebabc

         if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag )then
             storeoffdiag = .true.
             h = -croshorizmainbc(dew,           dns,            &
                                 slopex,        slopey,         &
                                 dsigmadew(up), dsigmadns(up),  &
                                 pt,            bcflag,         &
                                 dup(up),       local_othervel, &
                                 local_efvs,                    &
                                 oneortwo, twoorone, g_cros ) / scalebabc
             call fillsprsemain(h,locplusup,loc,up,pt)
             storeoffdiag = .false.
         end if

     else if( bcflag(2) /= 2 )then

          rhsd(locplusup) = sum( croshorizmainbc(dew,           dns,            &
                                             slopex,        slopey,         &
                                             dsigmadew(up), dsigmadns(up),  &
                                             pt,            bcflag,         &
                                             dup(up),       local_othervel, &
                                             local_efvs,                    &
                                             oneortwo, twoorone, g_cros )  &
                                              * local_othervel ) / scalebabc 

         if( nonlinear == HO_NONLIN_JFNK .and. calcoffdiag)then
             storeoffdiag = .true.
             h = -croshorizmainbc(dew,           dns,            &
                                 slopex,        slopey,         &
                                 dsigmadew(up), dsigmadns(up),  &
                                 pt,            bcflag,         &
                                 dup(up),       local_othervel, &
                                 local_efvs,                    &
                                 oneortwo, twoorone, g_cros ) / scalebabc
             call fillsprsemain(h,locplusup,loc,up,pt)
             storeoffdiag = .false.
         end if

      end if

      ! The following calculates the basal traction AFTER an updated solution is obtain by passing the new
      ! values of uvel, vvel back to the matrix assembly routines, and thus obtaining updated values of the 
      ! relevant coefficients. The if construct allows the assembly routines to be called for only the vert
      ! layers that are needed to cacluate the basal traction (as opposed to all vert levels 1:upn).
      if( assembly == 1 )then

      select case( pt )
         case(1)
           g_vel_lhs(:,:,:) = ghostbvel(1,:,ew-1:ew+1,ns-1:ns+1)
           g_vel_rhs(:,:,:) = ghostbvel(2,:,ew-1:ew+1,ns-1:ns+1)
         case(2)
           g_vel_lhs(:,:,:) = ghostbvel(2,:,ew-1:ew+1,ns-1:ns+1)
           g_vel_rhs(:,:,:) = ghostbvel(1,:,ew-1:ew+1,ns-1:ns+1)
       end select

       btraction(pt,ew,ns) = sum( (g_norm+g_vert)*g_vel_lhs*thk0/len0*sum(local_efvs(2,:,:))/4.0d0 ) &
                           - sum( g_cros*g_vel_rhs*thk0/len0*sum(local_efvs(2,:,:))/4.0d0 )

     end if



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

function vertimainbc(thck, bcflag, dup, efvs, betasquared, g_vert, nz, plastic_coeff)

! altered form of 'vertimain' that calculates coefficients for higher-order
! b.c. that go with the 'normhorizmain' term: -(X/H)^2 * dsigma/dzhat * du/dsigma 

    implicit none

    real (kind = dp), intent(in) :: dup, thck, betasquared 
    real (kind = dp), intent(in) :: nz                      ! sfc normal vect comp in z-dir
    real (kind = dp), intent(in), dimension(2,2,2) :: efvs
    real (kind = dp), intent(out), dimension(3,3,3) :: g_vert
    real (kind = dp), optional, intent(in) :: plastic_coeff
    integer, intent(in), dimension(2) :: bcflag

    real (kind = dp) :: c
    real (kind = dp), dimension(3) :: vertimainbc

    c = 0.0_dp
    g_vert = 0.0_dp

    ! for higher-order FREE SURFACE B.C. for x ('which'=1) or y ('which'=2) direction ...
    if( bcflag(1) == 1 )then

           c = nz / thck / (2*dup) * (len0**2 / thk0**2)   ! value of coefficient

           vertimainbc(:) = 0.0_dp
           vertimainbc(3) = -c
           vertimainbc(1) = c
           vertimainbc(2) = vertimainbc(3) + vertimainbc(1) ! should = 0

           ! this is the part of the vertimain coeff. block that we want to keep for calc
           ! of boundary tractions (note that it DOES NOT include terms from boundary forcing)
           g_vert(:,2,2) = vertimainbc

   ! for higher-order BASAL B.C. w/ specified basal traction, add on the necessary source term ...
    if( bcflag(2) == 1 )then

            ! last set of terms is mean visc. of ice nearest to the bed
            vertimainbc(2) = vertimainbc(2)   &
                           + ( betasquared / ( sum( efvs(2,:,:) ) / 4.0_dp ) ) * (len0 / thk0)
    end if

    ! for higher-order BASAL B.C. w/ plastic yield stress iteration ...
    if( bcflag(2) == 2 )then

             ! last set of terms is mean visc. of ice nearest to the bed
            vertimainbc(2) = vertimainbc(2)   &
                           + ( plastic_coeff / ( sum( efvs(2,:,:) ) / 4.0_dp ) ) * (len0 / thk0)
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
                         efvs,                       &
                         oneortwo,  twoorone,        &
                         g_cros, velbc, plastic_coeff )

    ! As described for "normhorizmainbc" above. The vectors "twoorone" and 
    ! "oneortwo" are given by: twoorone = [ 2 1 ]; oneortwo = [ 1 2 ];

    implicit none

    integer, intent(in) :: which
    integer, intent(in), dimension(:) :: bcflag

    real (kind = dp), intent(in) :: dew, dns
    real (kind = dp), intent(in), dimension(:) :: oneortwo, twoorone
    real (kind = dp), intent(in) :: dusrfdew, dusrfdns, dsigmadew, dsigmadns, dup
    real (kind = dp), intent(in), dimension(:,:,:) :: local_othervel
    real (kind = dp), intent(in), dimension(:,:,:) :: efvs 
    real (kind = dp), intent(in), optional :: velbc, plastic_coeff
    real (kind = dp), intent(out),dimension(:,:,:) :: g_cros


    real (kind = dp), dimension(3,3,3) :: g, croshorizmainbc
    real (kind = dp) :: c
    integer :: nz

    c = 0.0_dp
    g(:,:,:) = 0.0_dp
    g_cros = g
    nz = 0

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
    ! This forces the multiplication by 'local_otherval' in the main program 
    ! to result in a value of 1, thus leaving the boundary vel. unchanged
    ! ... conditional makes sure there is no div by zero if the bc value IS also zero
    else if( bcflag(1) == 0 )then

        g(:,:,:) = 0.0_dp

        where( local_othervel /= 0.0d0 )
            g = 1
        elsewhere
            g = 0.0d0
        endwhere

        nz = sum( g )
        g(:,:,:) = 0.0_dp

        where( local_othervel /= 0.0d0 )
            g = ( velbc / nz ) / local_othervel
        elsewhere
            g = 0.0d0
        endwhere

     end if

     if( bcflag(2) == 2 )then        ! add on coeff. associated w/ plastic bed iteration

         ! NOTE: here we define 'g_cros' FIRST, because we want the value w/o the plastic
         ! bed coeff. included (needed for estimate of basal traction in plastic bed iteration)
         g_cros = g
         g(2,2,2) = g(2,2,2) + plastic_coeff / ( sum( efvs(2,:,:) ) / 4.0_dp ) * (len0 / thk0)

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
           g(2,3,2) = c * whichbc(what)
           g(2,1,2) = -c * whichbc(what)

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
 
subroutine fillsprsemain(inp,locplusup,ptindx,up,pt)

  ! scatter coefficients from 3x3x3 block "g" onto sparse matrix row
  ! scatter coefficients from 3x3x3 block "g" onto sparse matrix row
  implicit none

  real (kind = dp), dimension(3,3,3), intent(in):: inp
  integer, intent(in) :: locplusup, up, pt
  integer, dimension(6), intent(in) :: ptindx

  ! insert entries to "g" that are on same level
  call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
  call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
  call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
  call putpcgc(inp(2,2,3),ptindx(4)+up,locplusup,pt)
  call putpcgc(inp(2,2,1),ptindx(5)+up,locplusup,pt)

  ! add points for level above (that is, points in "g"  with a LARGER first index,
  ! which correspond to grid points that are CLOSER TO THE BED than at current level)
  call putpcgc(inp(3,2,2),ptindx(1)+up+1,locplusup,pt)
  call putpcgc(inp(3,3,2),ptindx(2)+up+1,locplusup,pt)
  call putpcgc(inp(3,1,2),ptindx(3)+up+1,locplusup,pt)
  call putpcgc(inp(3,2,3),ptindx(4)+up+1,locplusup,pt)
  call putpcgc(inp(3,2,1),ptindx(5)+up+1,locplusup,pt)

  ! add points for level below (that is, points in "g" with a SMALLER first index,
  ! which correspond to grid points that are CLOSER TO THE SURFACE than at current level) 
  call putpcgc(inp(1,2,2),ptindx(1)+up-1,locplusup,pt)
  call putpcgc(inp(1,3,2),ptindx(2)+up-1,locplusup,pt)
  call putpcgc(inp(1,1,2),ptindx(3)+up-1,locplusup,pt)
  call putpcgc(inp(1,2,3),ptindx(4)+up-1,locplusup,pt)
  call putpcgc(inp(1,2,1),ptindx(5)+up-1,locplusup,pt)

  return

end subroutine fillsprsemain

!***********************************************************************

subroutine fillsprsebndy(inp,locplusup,ptindx,up,normal,pt)

  ! scatter coeff. from 3x3x3 block "g" onto sparse matrix row. This subroutine
  ! is specifically for the boundary conditions, which are handled differently
  ! than points in the "main" body of the domain (interior to boundaries).
  implicit none

  integer, intent(in) :: locplusup, up, pt
  integer, dimension(6), intent(in) :: ptindx
  real (kind = dp), dimension(3,3,3), intent(in) :: inp
  real (kind = dp), dimension(2), intent(in) :: normal

  ! at points where mixed centered and one-side diffs. would apply
  if( normal(1) == 0.0_dp )then         ! at boundary normal to y, centered diffs in x 
    if( normal(2) == -1.0_dp )then      ! at boundary w/ normal [0,-1]
           call putpcgc(inp(1,3,3),ptindx(5)+up-1,locplusup,pt)
           call putpcgc( inp(2,3,3)+inp(1,2,1),ptindx(5)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(5)+up+1,locplusup,pt)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup,pt)
    else                                ! at boundary w/ normal [0,1]
           call putpcgc(inp(1,3,3),ptindx(4)+up-1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(1,2,3),ptindx(4)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(4)+up+1,locplusup,pt)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup,pt)
    end if
    call putpcgc(inp(1,2,2),ptindx(1)+up,locplusup,pt)
  end if

  if( normal(2) == 0.0_dp )then            ! at boundary normal to x, centered diffs in y 
        if( normal(1) == -1.0_dp )then     ! at boundary w/ normal [-1,0]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup,pt)
           call putpcgc( inp(2,3,3)+inp(2,1,2),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
        else                                 ! at boundary w/ normal [1,0]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup,pt)
           call putpcgc( inp(2,3,3)+inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup,pt)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
    end if
    call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
  end if

  ! at corners where only one-side diffs. apply
  if( normal(1) .gt. 0.0_dp .and. normal(2) .ne. 0.0_dp )then
    if( normal(2) .gt. 0.0_dp )then      ! corner w/ normal [ 1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(2,3,2)+inp(1,2,3),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup,pt)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
    else                                 ! corner w/ normal [ 1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(2)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(2)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(1,2,1)+inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(2,1,2),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup,pt)
    end if
  end if

  if( normal(1) .lt. 0.0_dp .and. normal(2) .ne. 0.0_dp )then
    if( normal(2) .gt. 0.0_dp )then       ! corner w/ normal [ -1/sqrt(2), 1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(1,2,3)+inp(2,1,2),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(1,2,1),ptindx(5)+up,locplusup,pt)
    else                                  ! corner w/ normal [ -1/sqrt(2), -1/sqrt(2) ]
           call putpcgc(inp(1,3,3),ptindx(3)+up-1,locplusup,pt)
           call putpcgc(inp(3,3,3),ptindx(3)+up+1,locplusup,pt)
           call putpcgc(inp(2,3,3)+inp(2,1,2)+inp(1,2,1),ptindx(3)+up,locplusup,pt)
           call putpcgc(inp(2,2,2),ptindx(1)+up,locplusup,pt)
           call putpcgc(inp(1,2,2),ptindx(6)+up,locplusup,pt)
           call putpcgc(inp(2,3,2),ptindx(2)+up,locplusup,pt)
           call putpcgc(inp(1,2,3),ptindx(4)+up,locplusup,pt)
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
  real (kind = dp) :: smallnum = 1.0d-2
  real (kind = dp), dimension(ewn) :: grounded
  real (kind = dp) :: alpha, dx, thck_gl, betalow, betahigh, roughness
  integer :: ew, ns

  select case(whichbabc)

    case(0)     ! constant value; useful for debugging and test cases

      betasquared = 10.0d0

    case(1)     ! simple pattern; also useful for debugging and test cases
                ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      betasquared = 1.0d4

      do ew=1, ewn-1; do ns=5, nsn-5
        betasquared(ew,ns) = 10.0d1 
      end do; end do

!      betasquared = betasquared / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

    case(2)     ! take input value for till yield stress and force betasquared to be implemented such
                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

!      betasquared = minTauf*tau0 / dsqrt( (thisvel*vel0*scyr)**2 + (othervel*vel0*scyr)**2 + (smallnum)**2 )

!! *sfp* This hack for use w/ the basal processes test case, although as of Oct. 2010, it is working well
!! enough to just apply the yield stress basal bc everywhere in the domain.
!      betasquared = minTauf(1:ewn-9,:)*tau0 / dsqrt( (thisvel(1:ewn-9,:)*vel0*scyr)**2 + &
!                    (othervel(1:ewn-9,:)*vel0*scyr)**2 + (smallnum)**2 )
!      betasquared(ewn-8:ewn-1,:) = 5.0d0

!! *sfp* similar version to above, but holds a constant yield stress for the downstream portion of the domain 
!! rather than holding/applying a constant B^2. This necessary so that till module dynamics do not control the sliding 
!! at the downstream end of the domain, which we want to remain constant (we want it to remain slipper, as if it were
!! an active ice plain, regardless of the till dynamics in the ice stream proper)
      betasquared(1:ewn-9,:) = minTauf(1:ewn-9,:)*tau0 / dsqrt( (thisvel(1:ewn-9,:)*vel0*scyr)**2 + &
                    (othervel(1:ewn-9,:)*vel0*scyr)**2 + (smallnum)**2 )
      betasquared(ewn-8:ewn-1,:) = 1.0d3 / dsqrt( (thisvel(ewn-8:ewn-1,:)*vel0*scyr)**2 + &
                    (othervel(ewn-8:ewn-1,:)*vel0*scyr)**2 + (smallnum)**2 )

    case(3)     ! circular ice shelf: set B^2 ~ 0 except for at center, where B^2 >> 0 to enforce u,v=0 there

      betasquared = 1.0d-5
      betasquared( (ewn-1)/2:(ewn-1)/2+1, (nsn-1)/2:(nsn-1)/2+1 ) = 1.0d10

    case(4)    ! frozen (u=v=0) ice-bed interface

      betasquared = 1.0d8

    case(5)    ! use value passed in externally from CISM (NOTE not dimensional when passed in) 

      ! scale CISM input value to dimensional units of (Pa yrs 1/m)
      betasquared = beta * scyr * vel0 * len0 / (thk0**2)

      ! this is a check for NaNs, which indicate, and are replaced by no slip
      where ( betasquared /= betasquared )
        betasquared = 1.0d10
      end where

  end select

  ! convert whatever the specified value is to dimensional units of (Pa s m^-1 ) 
  ! and then non-dimensionalize using PP dyn core specific scaling.
  betasquared = ( betasquared * scyr ) / ( tau0 * tim0 / len0 )

end subroutine calcbetasquared

!***********************************************************************

subroutine plasticbediteration( ewn, nsn, uvel0, vvel0, btraction, minTauf, &
                            plastic_coeff_lhs, plastic_coeff_rhs, plastic_rhs, plastic_resid )

    implicit none

    integer, intent(in) :: ewn, nsn

    real (kind = dp), intent(in), dimension(:,:) :: uvel0, vvel0

    real (kind = dp), intent(in), dimension(:,:,:) :: btraction

    real (kind = dp), intent(in), dimension(ewn-1,nsn-1) :: minTauf      

    real (kind = dp), intent(out), dimension(2,ewn-1,nsn-1) :: plastic_coeff_lhs
    real (kind = dp), intent(out), dimension(2,ewn-1,nsn-1) :: plastic_coeff_rhs
    real (kind = dp), intent(out), dimension(2,ewn-1,nsn-1) :: plastic_rhs
    real (kind = dp), intent(inout), dimension(1,ewn-1,nsn-1) :: plastic_resid


    real (kind = dp), dimension(ewn-1,nsn-1) :: gamma, beta, tau

    real (kind = dp) :: c, big_C, relax
    real (kind = dp) :: F11, F12, F21, F22, M11, M12, M21, M22, h1, h2, IM11, IM12, IM21, IM22, &
                        L11, L12, L21, L22, K11, K12, K21, K22, det, rhs1, rhs2, alpha, &
                        theta, maxx, plastic_residx, plastic_residy, sigma_x0, sigma_y0

    integer :: ew, ns, countstuck, countslip

    countstuck = 0; countslip = 0; sigma_x0 = 0; sigma_y0 = 0

    tau = minTauf

    !! *sfp* the following is a hack to allow the new plastic bed implementation to be used with 
    !! the basal processes test case
    !tau(ewn-8:ewn-1,:) = 1.0d3 / tau0 

!    print *, 'inside plastic bed iteration subroutine (near top)'

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! initialize some values

    relax = 1.0d0       ! r.elaxation factor

!    c = 5.0d3 / (1.0d3/scyr)                ! "c" has unitis Pa * a/m (like B^2); choose value so that c*u is 
                                            ! approx equal to expected yield stress.
                                            ! here: 5000 Pa / (1000 m/a) = 5 Pa * a / m = 1.57e8 Pa * s / m
!    c = c * vel0 / tau0                     ! non-dimensionalize c                                        

!    c = 1.0d-10
!    c = 1.0d-5
!    c = 1.0d-3 
!    c = 1.0d-2      ! c < 1 needed for n=1 cases

    c = 1.0d0      ! c >= 1 needed for n=3 cases

!    c = 1.0d1
!    c = 1.0d3
!    c = 1.0d5
!    c = 1.0d8 

    big_C = 1.0d8 / tau0           ! regularization constant to allow small motion on boundary; dim value in Pa / tau0 

    gamma = 1.0d0 / big_C

     !! simple test case

!    tau = 100.0d3 / tau0      !! yield stress in Pa (NOTE that this eventually needs to be passed in) !!
!    tau(3:ewn-3,3:nsn-3) = 7.5d3 / tau0
!    tau(5:ewn-5,5:nsn-5) = 5.0d3 / tau0
!    tau(8:ewn-8,8:nsn-8) = 2.0d3 / tau0

!      tau = 1.0d10 / tau0    !! ice stream analytical soln test case
!      tau(3:ewn-3,3:nsn-3) = 6.5d3 / tau0

!!      do ew=1, ewn-1; do ns=4, nsn-4
!      do ew=1, ewn-1; do ns=8-1, nsn-8+1
!!      do ew=1, ewn-1; do ns=16-3, nsn-16+3
!!        tau(ew,ns) = 5.0d3 / tau0
!        tau(ew,ns) = 6.50d3 / tau0
!      end do; end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! loop through and calc. vel. bc based on conditional
    do ns=1, nsn-1
        do ew=1, ewn-1

            ! calculate sigma^t vector 
            sigma_x0 = -btraction(1,ew,ns)
            sigma_y0 = -btraction(2,ew,ns)

            beta(ew,ns) = sqrt( (sigma_x0 + c*uvel0(ew,ns))**2 + (sigma_y0 + c*vvel0(ew,ns))**2 )

            plastic_residx = max( tau(ew,ns), beta(ew,ns) ) * sigma_x0 - tau(ew,ns) * (sigma_x0 + &
                             c*uvel0(ew,ns))

            plastic_residy = max( tau(ew,ns), beta(ew,ns) ) * sigma_y0 - tau(ew,ns) * (sigma_y0 + &
                             c*vvel0(ew,ns))

            plastic_resid(1,ew,ns) = sqrt( plastic_residx**2 + plastic_residy**2 )

            ! if at or below yield stress, set boundary vels to 0 
            if( beta(ew,ns) <= ( tau(ew,ns)*(1.0d0 + c*gamma(ew,ns)) ) )then

!                print *, 'inside of plastic bed iteration loop (stuck nodes)'
                countstuck = countstuck + 1

                plastic_coeff_lhs(1,ew,ns) = big_C
                plastic_coeff_lhs(2,ew,ns) = big_C
                plastic_coeff_rhs(1,ew,ns) = 0.0d0
                plastic_coeff_rhs(2,ew,ns) = 0.0d0
                plastic_rhs(1,ew,ns) = 0.0d0
                plastic_rhs(2,ew,ns) = 0.0d0

            ! if above yield stress, solve for boundary vels
            elseif( beta(ew,ns) > ( tau(ew,ns)*(1.0d0 + c*gamma(ew,ns)) ) )then

!                print *, 'inside of plastic bed iteration loop (slipping nodes)'
                countslip = countslip + 1

                ! prefactor to matrix "M"
                theta = tau(ew,ns) / beta(ew,ns)

                ! max vel test in denom of several terms 
                maxx = maxval( (/ tau(ew,ns), sqrt( sigma_x0**2 + sigma_y0**2 ) /) )

                ! prefactor to matrix "F"
                alpha = 1.0d0 / ( beta(ew,ns) * maxx )

                ! entries to 2x2 matrix "F"
                F11 = alpha * sigma_x0 * ( sigma_x0 + c*uvel0(ew,ns) )
                F12 = alpha * sigma_x0 * ( sigma_y0 + c*vvel0(ew,ns) )
                F21 = alpha * sigma_y0 * ( sigma_x0 + c*uvel0(ew,ns) )
                F22 = alpha * sigma_y0 * ( sigma_y0 + c*vvel0(ew,ns) )

                ! entries to 2x2 matrix "M"
                M11 = theta * ( 1.0d0 - F11 )
                M12 = - theta * F12
                M21 = - theta * F21
                M22 = theta * ( 1.0d0 - F22 )


                ! entries to 2x2 matrix "I - M"
                IM11 = 1.0d0 - M11
                IM12 = -M12
                IM21 = -M21
                IM22 = 1.0d0 - M22

                ! entries to vector "h"
                h1 = tau(ew,ns) / maxx * sigma_x0
                h2 = tau(ew,ns) / maxx * sigma_y0

                ! entries to inverse matrix, K = (I-M)^(-1)
                det = IM11*IM22 - IM12*IM21
                K11 = IM22 / det
                K12 = -IM12 / det
                K21 = -IM21 / det
                K22 = IM11 / det

                ! entries to inverse matrix, L = ((I-M)^(-1)-I)
                L11 = c * (K11 - 1.0d0)
                L12 = c * K12
                L21 = c * K21
                L22 = c * (K22 - 1.0d0)

                ! calculate rhs vector components (NEW values, to be added on to existing RHS for u, v solns
                ! and later divided by N_eff 
                rhs1 = ( K11*h1 + K12*h2 )
                rhs2 = ( K21*h1 + K22*h2 )

                ! compile coeff. to send out of subroutine for new solution, which are the coeff. that go with u(2,2,2) 
                ! (in coeff. matrix on the LHS) (-L11), and v(2,2,2) (in summation RHS) (-L22) when solving for u 
                ! (and vice versa when solving for v). In each case, these get divided by the effective viscosity (???)

                plastic_coeff_lhs(1,ew,ns) = L11
                plastic_coeff_lhs(2,ew,ns) = L22
                plastic_coeff_rhs(1,ew,ns) = -L12
                plastic_coeff_rhs(2,ew,ns) = -L21
                plastic_rhs(1,ew,ns) = -rhs1
                plastic_rhs(2,ew,ns) = -rhs2

            end if

        end do
    end do

!    print *, '# stuck nodes = ', countstuck
!    print *, '# slip nodes = ', countslip 

end subroutine plasticbediteration

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

subroutine putpcgc(value,col,row,pt) 
 
  implicit none
 
  integer, intent(in) :: row, col
  integer, intent(in), optional :: pt
  real (kind = dp), intent(in) :: value 

   !*sfp* For now, ignoring the possibility of using JFNK w/ Trilinos ...
   if( nonlinear == HO_NONLIN_PICARD )then

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then
        ! Option to load entry into Triad sparse matrix format
        if (value /= 0.0d0) then
          pcgval(ct) = value
          pcgcol(ct) = col
          pcgrow(ct) = row
          ct = ct + 1
        end if
    else
        ! Option to load entry directly into Trilinos sparse matrix 
        if (value /= 0.0d0) then
           !AGS: If we find that sparsity changes inside a time step,
           !     consider adding entry even for value==0.
!            call putintotrilinosmatrix(row, col, value) ! JCC - No Trilinos support yet
        end if
    end if
 
 
   !*sfp* if using JFNK, store the main block diagonal coeffs and off diag coeffs 
   elseif ( nonlinear == HO_NONLIN_JFNK )then

    if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then      ! if using Triad format to store matrix entries

      if ( .not. storeoffdiag ) then ! block diag coeffs in normal storage space
          ! load entry into Triad sparse matrix format
          if (value /= 0.0d0) then
            pcgval(ct) = value
            pcgcol(ct) = col
            pcgrow(ct) = row
            ct = ct + 1
          end if
      else if ( storeoffdiag ) then ! off-block diag coeffs in other storage
          ! load entry into Triad sparse matrix format
          if( pt == 1 )then ! store uv coeffs 
              if (value /= 0.0d0) then
                pcgvaluv(ct2) = value
                pcgcoluv(ct2) = col
                pcgrowuv(ct2) = row
                ct2 = ct2 + 1
              end if
          else if( pt == 2 )then ! store vu coeffs
              if (value /= 0.0d0) then
                pcgvalvu(ct2) = value
                pcgcolvu(ct2) = col
                pcgrowvu(ct2) = row
                ct2 = ct2 + 1
              end if
          end if
      end if

    else    ! if storing matrix entires directly in Trilinos sparse format

      ! *sfp* NOTE: that this option does not allow for storing offidag terms at present
      ! (I assume that Andy can grab these when we know they are correct)
      if (.not. storeoffdiag) then
        if (value /= 0.0d0) then
           !AGS: If we find that sparsity changes inside a time step,
           !     consider adding entry even for value==0.
!          call putintotrilinosmatrix(row, col, value) ! JCC - No Trilinos support yet
        end if
      end if

    end if  ! end of "if using Triad or Trilinos storage format" construct
 
   end if   ! end of "if using Picard or JFNK for nonlinear solve" construct
   
  return
 
end subroutine putpcgc 

!***********************************************************************

function scalebasalbc( coeffblock, bcflag, lateralboundry, beta, efvs )

  ! *sfp* This function is used to scale the matrix coeffs and rhs vector coeff
  ! of the basal boundary condition when using JFNK for the nonlinear iteration
  ! (iteration on viscosity). 
  implicit none

  integer, dimension(2), intent(in) :: bcflag         
  logical :: lateralboundry
  real (kind = dp), dimension(:,:,:), intent(in) :: coeffblock 
  real (kind = dp), dimension(:,:,:), intent(in) :: efvs       
  real (kind = dp), intent(in) :: beta

  real (kind = dp) :: scale, scalebasalbc 

    if( nonlinear == 1 )then
        if( bcflag(1) == 1 )then

           ! use the dominant terms in the coeff associated with the velocity under consideration
           !scale = beta / ( sum( efvs(2,:,:) ) / 4.0_dp ) * (len0 / thk0)

           ! Use the magnitude of the coeff associated with the vert stress gradients. 
           ! NOTE that relevant coeffs are stored in diff parts of block depending 
           ! on type of boudnary     
           if( lateralboundry )then
               scale = abs( coeffblock(3,3,3) );  
           else
               scale = abs( coeffblock(3,2,2) );     
           end if                

           if( scale .le. 0.0d0 )then
            scale = 1.0d0
           end if

        else
            scale = 1.0d0
        end if

    else
        scale = 1.0d0
    end if

    scalebasalbc = scale

  return

end function scalebasalbc   

!  Subroutine for calculation of residual vect using L2norm. 
!  uvec is either u^k-1 or v^k-1 on input and Av-b or Cu-d on output
!  Written by J.F. Lemeiux
subroutine res_vect ( matrix, uvec, bvec, nu, counter, g_flag, L2square, whatsparse)


use glimmer_paramets, only : dp
use glimmer_sparse_type
use glimmer_sparse
use glide_mask

implicit none

integer :: i, j, counter, nu, nele, whatsparse ! nu: size of uvec and bvec
integer, dimension(nu), intent(in) :: g_flag ! 0 :reg cell
                                             ! 1 :top ghost, 2 :base ghost

type(sparse_matrix_type),  intent(in) :: matrix

real (kind = dp), dimension(nu), intent(in) :: bvec
real (kind = dp), dimension(nu), intent(inout) :: uvec
real (kind = dp), dimension(nu) :: Au_b_wig
real (kind = dp), intent(out) :: L2square
! 
real (kind = dp) :: scale_ghosts = 0.0d0

! calculate residual vector of the u OR v component

      Au_b_wig = 0d0 ! regular+ghost cells

!      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then

        do nele = 1, matrix%nonzeros

           i = matrix%row(nele)
           j = matrix%col(nele)
           Au_b_wig(i) = Au_b_wig(i) + matrix%val(nele) * uvec(j)

        enddo

!      else
!        call matvecwithtrilinos(uvec, Au_b_wig);
!      endif

      do i = 1, nu
         Au_b_wig(i) = Au_b_wig(i) - bvec(i)
      enddo

      uvec = Au_b_wig

! AGS: Residual norm includes scaling to decrease importance of ghost values
! By calling it a redefinition of an inner product, it is kosher.
      L2square = 0.0
      do i = 1, nu
         if (g_flag(i) .eq. 0) then
            L2square = L2square + Au_b_wig(i) * Au_b_wig(i)
         else
            L2square = L2square + scale_ghosts * Au_b_wig(i) * Au_b_wig(i)
         endif
      end do

      return

end subroutine res_vect

subroutine output_res( ewn, nsn, upn, uindx, counter, nu, reswrapped, u_or_v)
! This routine can be used for debugging and evalution of the nonlinear 
! numerical scheme (Picard iteration...and later Newton). It unwraps the 
! residual on the grid, either for the u comp (u_or_v = 1) or the 
! v comp (u_or_v = 2).
! Written by J.F. Lemeiux

use glimmer_paramets, only : dp
use netcdf

  implicit none

  integer, intent(in) :: ewn, nsn, upn, counter, nu, u_or_v
  integer, dimension(ewn-1,nsn-1), intent(in) :: uindx
  real (kind = dp), dimension(nu), intent(in) :: reswrapped
  real (kind = dp), dimension(ewn-1,nsn-1, upn) :: resgrid

  integer, dimension(2) :: loc
  integer :: ew, ns, up, upl, nc_id, ns_id, ew_id, up_id, res_id, status
  integer :: dimids(3)

  character filename*30

      do ns = 1,nsn-1 !had to flip indices in matrix because of netcdf
         do ew = 1,ewn-1
            if (uindx(ew,ns) /= 0) then
               loc(:) = (uindx(ew,ns) - 1) * (upn + 2) + 1 + (/ 1, upn /)
! loc is the same as given by the function getlocrange in glam_str2.F90
               resgrid(ew,ns,:) = abs(reswrapped(loc(1):loc(2)))
            else
!               resgrid(ew,ns,:) = -999d0 ! no ice
               resgrid(ew,ns,:) = -1d0 ! no ice
            end if
         end do
      end do

      if ( u_or_v .eq. 1) then
         write (filename, '("res_ucomp_k",i2.2,".nc")') counter
      elseif ( u_or_v .eq. 2) then
         write (filename, '("res_vcomp_k",i2.2,".nc")') counter
      endif

      call check(nf90_create(filename, NF90_CLOBBER, nc_id))
      call check(nf90_def_dim(nc_id,'n-s',nsn-1,ns_id))
      call check(nf90_def_dim(nc_id,'e-w',ewn-1,ew_id))
      call check(nf90_def_dim(nc_id,'up',upn,up_id))

      dimids(3) = up_id
      dimids(1) = ew_id
      dimids(2) = ns_id

      call check(nf90_def_var(nc_id,'residual',NF90_DOUBLE, dimids, res_id))
      call check( nf90_enddef(nc_id) )
      call check(nf90_put_var(nc_id,res_id,resgrid))
      call check(nf90_close(nc_id))

end subroutine output_res

subroutine check(status_code)
      use netcdf
      integer,intent(in) :: status_code

      if (status_code /= 0) then
         write(*,*) 'fatal netcdf error:',nf90_strerror(status_code)
         stop
      end if
end subroutine check


!***********************************************************************

end module glam_strs2

!***********************************************************************
