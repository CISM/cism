
!TODO - Add header comments here.

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module cism_sparse_pcg

  use glimmer_global, only: dp
  use glimmer_sparse_type, only: sparse_matrix_type

  implicit none

!WHL - debug
    logical :: verbose_pcg

!WHL - debug
    integer, parameter :: &
        ndiagmax = 10      ! number of values to print out for debugging

!WHL - debug

   integer, parameter :: &
       itest = 24, jtest = 17, ktest = 1

   integer, parameter :: ntest = 2371  ! nodeID for (24,17,1)

contains

!****************************************************************************

  subroutine pcg_solver_structured(nx,        ny,            &
                                   nz,        active_vertex, &
                                   Auu,       Auv,           &
                                   Avu,       Avv,           &
                                   bu,        bv,            &
                                   xu,        xv,            &
                                   err,       niters,        &
                                   NodeID,    verbose)   ! NodeID is temporary for debugging

    !---------------------------------------------------------------
    !  This subroutine uses a preconditioned conjugate-gradient solver to
    !  solve the equation $Ax=b$.
    !  Convergence is checked every {\em ncheck} steps.
    !
    !  It is based on the barotropic solver in the POP ocean model 
    !  (author Phil Jones, LANL).  Input and output arrays are located
    !  on a structured (i,j,k) grid as defined in the glissade_velo_higher
    !  module.  The global matrix is sparse, but its nonzero elements
    !  are stored in four dense matrices called Auu, Avv, Auv, and Avu.
    !  Each matrix has 3x3x3 = 27 potential nonzero elements per
    !  node (i,j,k).
    !
    !  The current preconditioning options are
    !  (0) no preconditioning
    !  (1) diagonal preconditioning
    !  (2) preconditioning using a physics-based SIA solver
    ! 
    !  For the dome test case with higher-order dynamics, option (2) is best. 
    !
    !  This subroutine is less flexible than the generic PCG solver below,
    !  which works with grid-independent matrices in triad format.
    !  However, this subroutine is easier to parallelize because
    !  it can use existing halo updates and global reductions.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! horizontal grid dimensions (for scalars)
                                ! velocity grid has dimensions (nx-1,ny-1)
       nz                       ! number of vertical levels where velocity is computed

    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F
 
    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 3 (node and its 2 neighbors in z direction)
                              ! 2nd dimension = 3 (node and its 2 neighbors in x direction)
                              ! 3rd dimension = 3 (node and its 2 neighbors in y direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::  &
       bu, bv             ! assembled load (rhs) vector, divided into 2 parts

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       xu, xv             ! u and v components of solution (i.e., uvel and vvel)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

!WHL - Temporary for debugging
    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       NodeID

    logical, intent(in), optional :: &
       verbose                           ! if true, print diagnostic output
    
    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  i, j, k      ! grid indices
    integer ::  iA, jA, kA   ! grid offsets ranging from -1 to 1
    integer ::  m            ! iteration counter

    real(dp) ::           &
       eta0, eta1, rr      ! scalar inner product results

    ! vectors (each of these is split into u and v components)
    real(dp), dimension(nz,nx-1,ny-1) ::  &
       Adiagu, Adiagv,    &! diagonal terms of matrices Auu and Avv
       ru, rv,            &! residual vector (b-Ax)
       du, dv,            &! conjugate direction vector
       yu, yv,            &! result of a matvec multiply
       work0u, work0v,    &! cg intermediate results
       work1u, work1v      ! cg intermediate results

    real(dp) ::  &
       L2_resid,          &! L2 norm of residual vector Ax-b
       L2_rhs              ! L2 norm of rhs vector b
                           ! solver converges when L2_resid/L2_rhs < tolerance

    real(dp), dimension(-1:1,nz,nx-1,ny-1) ::  &
       Muu, Mvv            ! simplified SIA matrices for preconditioning

    !---------------------------------------------------------------
    ! Solver parameters
    ! TODO: Pass these in as arguments?
    !---------------------------------------------------------------

    real(dp), parameter ::   &
       tolerance = 1.d-11    ! tolerance for linear solver

    integer, parameter ::    &
       maxiters = 200         ! max number of linear iterations before quitting

    integer, parameter :: &
        solv_ncheck = 1    ! check for convergence every solv_ncheck iterations
                           ! for now, hardwire for debugging
                           ! later, let this be an argument

    integer, parameter :: &
        precond = 2        ! = 0 for no preconditioning
                           ! = 1 for diagonal preconditioning (does not work well for SIA-dominated flow)
                           ! = 2 for preconditioning with SIA solver (works well for SIA-dominated flow)

!WHL - debug
    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .true.   ! for debugging
    endif

    if (verbose_pcg) then
       print*, ' '
       print*, 'In structured pcg solver'
       print*, 'tolerance, maxiters =', tolerance, maxiters
       print*, ' '
    endif

    ! Set up preconditioning

    if (precond == 1) then    ! form diagonal matrix for preconditioning

       Adiagu(:,:,:) = Auu(0,0,0,:,:,:)
       Adiagv(:,:,:) = Avv(0,0,0,:,:,:)

!TODO - Get rid of NodeID later
       if (verbose_pcg) then
          print*, ' '
          print*, 'Using diagonal solver for preconditioning'
          print*, ' '
          print*, 'i, j, k, n, diagonal entries, initial guess, residual:'
          m = 0
          do j = 1, ny-1
          do i = 1, nx-1
          do k = 1, nz
             if (Adiagu(k,i,j) > 0.d0 .and. m < ndiagmax) then
                m = m + 2
                print*, i, j, k, NodeID(k,i,j), Adiagu(k,i,j), xu(k,i,j), ru(k,i,j)
                print*, i, j, k, NodeID(k,i,j), Adiagv(k,i,j), xv(k,i,j), rv(k,i,j)
             else
                exit
             endif
          enddo
          enddo
          enddo
       endif  ! verbose_pcg

    elseif (precond == 2) then  ! form SIA matrices Muu and Mvv with vertical coupling only

       Muu(:,:,:,:) = 0.d0
       Mvv(:,:,:,:) = 0.d0

       do j = 1, ny-1
       do i = 1, nx-1
       do k = 1, nz
          do kA = -1,1
             Muu(kA,k,i,j) = Auu(kA,0,0,k,i,j)   ! remove horizontal coupling
             Mvv(kA,k,i,j) = Avv(kA,0,0,k,i,j)
          enddo
       enddo
       enddo
       enddo

!WHL - debug
       if (verbose_pcg) then
          print*, ' '
          print*, 'Using easy SIA solver for preconditioning'
          i = itest
          j = jtest
          print*, ' '
          print*, 'i, j =', i, j
          print*, ' '
          print*, 'k, Muu(-1:1):'
          do k = 1, nz
             print*, k, Muu(-1:1,k,i,j)
          enddo
          print*, ' '
          print*, 'k, Mvv(-1:1):'
          do k = 1, nz
             print*, k, Mvv(-1:1,k,i,j)
          enddo
       endif

    endif      ! precond

    !  Compute initial residual and initialize the direction vector d

    !TODO - Could use active_vertex array to set up indirect addressing
    !       and avoid an 'if' in the matvec subroutine

    call matvec_multiply_structured(nx,        ny,            &
                                    nz,        active_vertex, &
                                    Auu,       Auv,           &
                                    Avu,       Avv,           &
                                    xu,        xv,            &
                                    yu,        yv, NodeID)

    ru(:,:,:) = bu(:,:,:) - yu(:,:,:)
    rv(:,:,:) = bv(:,:,:) - yv(:,:,:)

    du(:,:,:) = 0.d0
    dv(:,:,:) = 0.d0

    ! POP does a halo update of the residual here.
    ! (I think this is needed to make sure r is correct in halo cells.
    !  The r computed above will be wrong in halo cells if y is not correct in halo cells.)
    ! call halo_update(r)

    ! Initialize fields and scalars

    eta0 = 1.d0 
    niters = maxiters 

    ! Compute L2 norm of rhs vectors
    ! (Goal is to obtain L2_resid/L2_rhs < tolerance)

    work0u(:,:,:) = bu(:,:,:)*bu(:,:,:)
    work0v(:,:,:) = bv(:,:,:)*bv(:,:,:)

    ! POP does a halo update here before taking the global sum
    ! (I don't see why this is needed)

    L2_rhs = sum(work0u) + sum(work0v)       ! (b, b) = squared L2 norm
    ! In parallel, this would need to be a global sum: L2_rhs = global_sum(work0)

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of rhs

    ! Iterate to solution

    iter_loop: do m = 1, maxiters

       ! Compute (PC)r

       if (precond == 0) then      ! no preconditioning

           work1u(:,:,:) = ru(:,:,:)         ! PC(r) = r 
           work1v(:,:,:) = rv(:,:,:)         ! PC(r) = r 

       elseif (precond == 1 ) then  ! diagonal preconditioning

          do j = 1, ny-1
          do i = 1, nx-1
          do k = 1, nz
             if (Adiagu(k,i,j) /= 0.d0) then
                work1u(k,i,j) = ru(k,i,j) / Adiagu(k,i,j)   ! PC(r), where PC is formed from diagonal elements of A
             else                                        
                work1u(k,i,j) = 0.d0
             endif
             if (Adiagv(k,i,j) /= 0.d0) then
                work1v(k,i,j) = rv(k,i,j) / Adiagv(k,i,j)  
             else                                        
                work1v(k,i,j) = 0.d0
             endif
          enddo    ! k
          enddo    ! i
          enddo    ! j

!WHL - debug
          if (verbose_pcg) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'Diagonal preconditioning, i, j =', i, j
             print*, ' '
             print*, 'k, ru, PC(ru):'
             do k = 1, nz
                print*, k, ru(k,i,j), work1u(k,i,j)
             enddo
             print*, ' '
             print*, 'k, rv, PC(rv):'
             do k = 1, nz
                print*, k, rv(k,i,j), work1v(k,i,j)
             enddo
           endif

       elseif (precond == 2) then   ! local vertical shallow-ice solver for preconditioning

          call easy_sia_solver(nx,   ny,    nz,       &
                               active_vertex,         &
                               Muu,  ru,   work1u)      ! solve Muu*work1u = ru 

          call easy_sia_solver(nx,   ny,    nz,       &
                               active_vertex,         &
                               Mvv,  rv,   work1v)      ! solve Mvv*work1v = rv 

!WHL - debug
          if (verbose_pcg) then
             i = itest
             j = jtest
             print*, ' '
             print*, 'SIA preconditioning, i, j =', i, j
             print*, ' '
             print*, 'k, ru, PC(ru):'
             do k = 1, nz
                print*, k, ru(k,i,j), work1u(k,i,j)
             enddo
             print*, ' '
             print*, 'k, rv, PC(rv):'
             do k = 1, nz
                print*, k, rv(k,i,j), work1v(k,i,j)
             enddo
           endif

       endif    ! precond

       work0u(:,:,:) = ru(:,:,:)*work1u(:,:,:)    ! terms of dot product (r, PC(r))
       work0v(:,:,:) = rv(:,:,:)*work1v(:,:,:)    ! terms of dot product (r, PC(r))

       ! POP has a halo update here for parallel solve with external preconditioner (lprecond = T)
       ! if (lprecond) call halo_update(work1)
       ! Not sure why this is done after computing work0 = r*work1

       ! Compute (r, PC(r))

       ! In parallel, this would need to be a global sum: eta1 = global_sum(work0u) + global_sum(work0v)
       eta1 = sum(work0u) + sum(work0v)

!       if (verbose_pcg) print*,'eta0, eta1 =', eta0, eta1

       ! Update conjugate direction vector d

       du(:,:,:) = work1u(:,:,:) + du(:,:,:)*(eta1/eta0)   ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
       dv(:,:,:) = work1v(:,:,:) + dv(:,:,:)*(eta1/eta0)   !
                                                           !                    (r_(i+1), PC(r_(i+1)))
                                                           ! where beta_(i+1) = --------------------  
                                                           !                        (r_i, PC(r_i)) 
                                                           ! Initially eta0 = 1  
                                                           ! For m >=2, eta0 = old eta1

       ! Compute y = A*d

       call matvec_multiply_structured(nx,        ny,            &
                                       nz,        active_vertex, &
                                       Auu,       Auv,           &
                                       Avu,       Avv,           &
                                       du,        dv,            &
                                       yu,        yv, NodeID)

       work0u(:,:,:) = yu(:,:,:) * du(:,:,:)       ! terms of dot product (d, Ad)
       work0v(:,:,:) = yv(:,:,:) * dv(:,:,:)

       ! Compute next solution and residual

       ! Here POP does a halo update for y = A*d
       ! call halo_update(y)

        eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 
                                  

        ! These sums need to be global sums in parallel.
                                                     !          (r, PC(r))
        eta1 = eta0 / (sum(work0u) + sum(work0v))    ! alpha = ----------
                                                     !          (d, A*d)


        if (verbose_pcg) then
!           print*, '(r, PC(r))=', eta0
!           print*, 'alpha =', eta1
        endif

        xu(:,:,:) = xu(:,:,:) + eta1 * du(:,:,:)    ! new solution, x_(i+1) = x_i + alpha*d
        xv(:,:,:) = xv(:,:,:) + eta1 * dv(:,:,:)

        ru(:,:,:) = ru(:,:,:) - eta1 * yu(:,:,:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)
        rv(:,:,:) = rv(:,:,:) - eta1 * yv(:,:,:)

        if (verbose_pcg) then
!           print*, ' '
!           print*, 'New solution, residual:'
!           do i = 1, min(ndiagmax,order)
!              print*, x(i), r(i)
!           enddo
        endif

!TODO - Can we evaluate the L2 norm of the residual each iteration without doing another matvec multiply?
!       Just have to use the residual computed above instead of computing Ax - b

        if (mod(m,solv_ncheck) == 0) then

           ! Compute y = Ax
           
          call matvec_multiply_structured(nx,        ny,            &
                                          nz,        active_vertex, &
                                          Auu,       Auv,           &
                                          Avu,       Avv,           &
                                          xu,        xv,            &
                                          yu,        yv, NodeID)

           ! Compute residual and terms of its L2 norm
           ! TODO - How does this differ from the residual computed above?

           ru(:,:,:) = bu(:,:,:) - yu(:,:,:)
           rv(:,:,:) = bv(:,:,:) - yv(:,:,:)

           work0u(:,:,:) = ru(:,:,:)*ru(:,:,:)
           work0v(:,:,:) = rv(:,:,:)*rv(:,:,:)

           ! POP does a halo update here (Why? Should it be done before computing work0?)
           ! call halo_update(r)
           
           ! In parallel, these would need to be global sums
           L2_resid = sum(work0u) + sum(work0v)        ! (r, r) = squared L2 norm

           L2_resid = sqrt(L2_resid)                   ! L2 norm of residual

           err = L2_resid/L2_rhs

           if (verbose_pcg) print*, 'iter, L2_resid, error =', m, L2_resid, err

           if (err < tolerance) then
              niters = m
              exit iter_loop
           endif            

        endif    ! solv_ncheck

     enddo iter_loop

!WHL - Without good preconditioning, convergence is slow, but the solution after maxiters may be good enough.
 
     if (niters == maxiters) then
        print*, 'Glissade PCG solver not converged'
        print*, 'niters, err, tolerance:', niters, err, tolerance
!!!        stop   ! TODO - Abort cleanly
     endif

  end subroutine pcg_solver_structured

!****************************************************************************

!WHL - NodeID is temporary

  subroutine matvec_multiply_structured(nx,        ny,            &
                                        nz,        active_vertex, &
                                        Auu,       Auv,           &
                                        Avu,       Avv,           &
                                        xu,        xv,            &
                                        yu,        yv, NodeID)

    !---------------------------------------------------------------
    ! Compute the matrix-vector product $y = Ax$.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

!WHL - Temporary for debugging
    integer, dimension(nz,nx-1,ny-1), intent(in) ::  &
       NodeID

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz                     ! number of vertical levels

    !TODO - Replace with indirect addressing?
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(-1:1,-1:1,-1:1,nz,nx-1,ny-1), intent(in) ::   &
       Auu, Auv, Avu, Avv     ! four components of assembled matrix
                              ! 1st dimension = 3 (node and its 2 neighbors in z direction)
                              ! 2nd dimension = 3 (node and its 2 neighbors in x direction)
                              ! 3rd dimension = 3 (node and its 2 neighbors in y direction)
                              ! other dimensions = (z,x,y) indices
                              !
                              !    Auu  | Auv
                              !    _____|____
                              !    Avu  | Avv
                              !         |

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       xu, xv             ! current guess for solution


    real(dp), dimension(nz,nx-1,ny-1), intent(out) ::  &
       yu, yv             ! y = Ax

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, k
    integer :: iA, jA, kA
    
!WHL= debug
    real(dp) :: maxyu, maxyv
    integer :: imaxu, jmaxu, kmaxu
    integer :: imaxv, jmaxv, kmaxv
!!    integer :: nid

    if (verbose_pcg) then
!       print*, ' '
!       print*, 'In structured matvec'
!       print*, 'max(Auu), max(Auv), max(Avu), max(Avv) =', &
!                maxval(Auu), maxval(Auv), maxval(Avu), maxval(Avv)
!       print*, 'sum(xu), sum(xv), sumtot =', sum(xu), sum(xv), sum(xu) + sum(xv)
!       print*, 'max(xu), max(xv) =', maxval(xu), maxval(xv)
    endif

    ! Initialize the result vector.

    yu(:,:,:) = 0.d0
    yv(:,:,:) = 0.d0

    ! Compute y = Ax

    !TODO - Change loop bounds so we don't go out of bounds?
    !       For now the full bounds should work as long as there is no
    !        ice at the global domain boundaries.
    do j = 1, ny-1
    do i = 1, nx-1

       !TODO - Use indirect addressing to avoid an 'if' here?
       if (active_vertex(i,j)) then

!WHL - debug
!!          print*, ' '
!!          print*, 'Active: i, j =', i, j

          do k = 1, nz

!!             nid = NodeID(k,i,j)
!!             if (nid == ntest) print*, 'ntest, xu, xv:', nid, xu(k,i,j), xv(k,i,j)

             !TODO - Could replace these three short loops with long multadds for better GPU efficiency
             do jA = -1,1
             do iA = -1,1
             do kA = -1,1

!TODO - Can we somehow get rid of this 'if' statement and still keep the xu/xv indices in bounds?
                if ( (k+kA >= 1 .and. k+kA <= nz)     &
      	                        .and.                     &
                     (i+iA >= 1 .and. i+iA <= nx-1)         &
                             .and.                     &
                     (j+jA >= 1 .and. j+jA <= ny-1) ) then

                   yu(k,i,j) = yu(k,i,j)   &
                             + Auu(kA,iA,jA,k,i,j)*xu(k+kA,i+iA,j+jA)  &
                             + Auv(kA,iA,jA,k,i,j)*xv(k+kA,i+iA,j+jA)

                   yv(k,i,j) = yv(k,i,j)   &
                             + Avu(kA,iA,jA,k,i,j)*xu(k+kA,i+iA,j+jA)  &
                             + Avv(kA,iA,jA,k,i,j)*xv(k+kA,i+iA,j+jA)
 
!WHL - debug
!!                   if (NodeID(k,i,j) == ntest) then
!!                      print*, 'nid, iA, jA, kA, Auu, xu, u, Auu*xu:',  &
!!                               nid, iA, jA, kA, Auu(kA,iA,jA,k,i,j), xu(k+kA,i+iA,j+jA), Auu(kA,iA,jA,k,i,j)*xu(k+kA,i+iA,j+jA)
!!                      print*, 'nid, iA, jA, kA, Auv, xv, v, Auv*xv:',  &
!!                               nid, iA, jA, kA, Auv(kA,iA,jA,k,i,j), xv(k+kA,i+iA,j+jA), Auv(kA,iA,jA,k,i,j)*xv(k+kA,i+iA,j+jA)
!!                      print*, 'new yu:', yu(k,i,j)
!!                   endif
                   
                endif   ! k+kA, i+iA, j+jA in bounds

             enddo   ! kA
             enddo   ! iA
             enddo   ! jA

          enddo   ! k

       endif   ! active_vertex

    enddo   ! i
    enddo   ! j

    if (verbose_pcg) then
!       print*, ' '
!       print*, 'Find max yu, yv'
!       maxyu = -9999.d0
!       maxyv = -9999.d0
!       do j = 1, ny-1
!       do i = 1, nx-1
!       do k = 1, nz
!          if (yu(k,i,j) > maxyu) then
!             maxyu = yu(k,i,j)
!             imaxu = i
!             jmaxu = j
!             kmaxu = k
!          endif
!          if (yv(k,i,j) > maxyv) then
!             maxyv = yv(k,i,j)
!             imaxv = i
!             jmaxv = j
!             kmaxv = k
!          endif
!       enddo
!       enddo
!       enddo

!       print*, 'i, j, k, NodeID, maxyu:', imaxu, jmaxu, kmaxu, NodeID(kmaxu,imaxu,jmaxu), maxyu
!       print*, 'i, j, k, NodeID, maxyv:', imaxv, jmaxv, kmaxv, NodeID(kmaxv,imaxv,jmaxv), maxyv

    endif  ! verbose_pcg

  end subroutine matvec_multiply_structured
 
!****************************************************************************

  subroutine easy_sia_solver(nx,   ny,   nz,    &
                             active_vertex,     & 
                             A,    b,    x)

    !---------------------------------------------------------------
    ! Solve the problem Ax = y where A is a local shallow-ice matrix,
    !  with coupling in the vertical but not the horizontal.
    ! We simply solve a tridiagonal matrix for each column.
    !---------------------------------------------------------------
   
    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       nx, ny,             &  ! horizontal grid dimensions (for scalars)
       nz                     ! number of vertical levels
           
    logical, dimension(nx-1,ny-1), intent(in) ::   &
       active_vertex          ! T for columns (i,j) where velocity is computed, else F

    real(dp), dimension(-1:1,nz,nx-1,ny-1), intent(in) ::   &
       A                      ! matrix with vertical coupling only
                              ! 1st dimension = node and its upper and lower neighbors

    real(dp), dimension(nz,nx-1,ny-1), intent(in) ::   &
       b                      ! right-hand side

    real(dp), dimension(nz,nx-1,ny-1), intent(inout) ::   &
       x                      ! solution

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    real(dp), dimension(nz) ::  &
       sbdiag,         &      ! subdiagonal matrix entries
       diag,           &      ! diagonal matrix entries
       spdiag,         &      ! superdiagonal matrix entries
       rhs,            &      ! right-hand side
       soln                   ! tridiagonal solution

    integer :: i, j, k

!TODO - Can remove this temporary test code

!WHL - debug
    logical, parameter :: test_tridiag = .false.

    if (test_tridiag) then

       print*, ' '
       print*, 'Testing tridiag solver'

       diag(:) = 0.d0
       sbdiag(:) = 0.d0
       spdiag(:) = 0.d0
       soln(:) = 0.d0

       do k = 1, nz
          diag(k) = 2.d0
          if (k > 1)  sbdiag(k) = -1.d0
          if (k < nz) spdiag(k) = -1.d0
       enddo

       rhs(1) = 1.d0
       rhs(nz) = 1.d0
       rhs(2:nz-1) = 0.d0              ! answer = (1 1 1 ... 1 1 1)

       call tridiag_solver(nz,    sbdiag,   &
                           diag,  spdiag,   &
                           rhs,   soln)
 
       print*, ' '
       print*, 'Solution =', soln(:)

       stop

    endif   ! test_tridiag

!WHL - Real code starts here

    do j = 1, ny-1
    do i = 1, nx-1

       if (active_vertex(i,j)) then

          ! initialize rhs and solution 

          rhs(:)  = b(:,i,j)
          soln(:) = x(:,i,j)

          ! top layer

          k = 1
          sbdiag(k) = 0.d0
          diag(k)   = A(0,k,i,j)
          spdiag(k) = A(1,k,i,j)

          ! intermediate layers

          do k = 2, nz-1
             sbdiag(k) = A(-1,k,i,j)
             diag(k)   = A( 0,k,i,j)
             spdiag(k) = A( 1,k,i,j)
          enddo

          ! bottom layer

          k = nz 
          sbdiag(k) = A(-1,k,i,j)
          diag(k)   = A( 0,k,i,j)
          spdiag(k) = 0.d0           
         
          ! solve

          call tridiag_solver(nz,    sbdiag,   &
                              diag,  spdiag,   &
                              rhs,   soln)

          x(:,i,j) = soln(:)

       endif  ! active_vertex

    enddo     ! i
    enddo     ! j

  end subroutine easy_sia_solver
  
!****************************************************************************

  subroutine tridiag_solver(order,    sbdiag,   &
                            diag,     spdiag,   &
                            rhs,      soln)

    !---------------------------------------------------------------
    ! Solve a 1D tridiagonal matrix problem.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       order            ! matrix dimension

    real(dp), dimension(order), intent(in) :: &
       sbdiag,        & ! sub-diagonal matrix elements
       diag,          & ! diagonal matrix elements
       spdiag,        & ! super-diagonal matrix elements
       rhs              ! right-hand side

    real(dp), dimension(order), intent(inout) :: &
       soln             ! solution vector

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: &
       k               ! row counter

    real(dp) :: &
       beta            ! temporary matrix variable

    real(dp), dimension(order) :: &
       gamma           ! temporary matrix variable


    ! Solve

    beta = diag(1)
    soln(1) = rhs(1) / beta

    do k = 2, order
       gamma(k) = spdiag(k-1) / beta
       beta = diag(k) - sbdiag(k)*gamma(k)
       soln(k) = (rhs(k) - sbdiag(k)*soln(k-1)) / beta
    enddo

    do k = order-1, 1, -1
       soln(k) = soln(k) - gamma(k+1)*soln(k+1)
    enddo

  end subroutine tridiag_solver

!****************************************************************************

!TODO - Move the rest of this module to a different module?
!       The remaining subroutines are associated with a SLAP-type PCG 
!        solver that uses diagonal preconditioning.
!       They reproduce SLAP functionality without having to call SLAP itself,
!        but otherwise they may not be very useful because diagonal
!        preconditioning isn't generally very efficient.

!****************************************************************************

  subroutine cism_sparse_pcg_solve(matrix,  rhs,     solution, &
                                   err,     niters,  verbose)

    !---------------------------------------------------------------
    ! Preconditioned conjugate gradient solver based on the PCG solver 
    ! in the POP ocean model (author Phil Jones, LANL).
    ! Also includes some Fortran90 versions of SLAP subroutines.
    !
    ! This solver is functionally equivalent to the SLAP PCG solver with 
    ! diagonal preconditioning, but a new Fortran90 version is included
    ! here so that we can test and modify it more easily than the SLAP solver.
    ! Currently serial only, but could be made parallel with some more work.
    !
    ! The incoming matrix is assumed be in triad (row, col, val) format.
    ! It is converted here to compressed sparse column (CSC) or compressed 
    ! sparse row (CRS) format for efficient matrix-vector operations. 
    ! See details below.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    type(sparse_matrix_type), intent(in) ::  &
       matrix       ! matrix (assumed to be symmetric positive-definite for PCG)

    real(dp), dimension(:), intent(in) :: &
       rhs          ! right-hand side vector

    real(dp), dimension(:), intent(inout) :: &
       solution     ! solution vector

    real(dp), intent(out) ::  &
       err          ! error in final solution
        
    integer, intent(out) ::  &
       niters       ! number of iterations required to reach the solution

    logical, optional, intent(in) ::  &
       verbose      ! if present and true, then print diagnostics
       
    !---------------------------------------------------------------
    ! Local variables    
    !---------------------------------------------------------------

    integer, dimension(matrix%nonzeros) ::  &
       iA             ! initally refers to the row of each nonzero entry,
                      ! but for CSR storage is converted to an array that holds offsets
                      ! into the A and jA arrays for the start of each row

    integer, dimension(matrix%nonzeros) ::  &
       jA             ! initially refers to the column of each nonzero entry,
                      ! but for CSC storage is converted to an array that holds offsets
                      ! into the A and iA arrays for the start of each row

    real(dp), dimension(matrix%nonzeros) ::   &
       A              ! nonzero matrix entries

    logical ::   &
       sorted         ! true if triad matrix entries are in ascending row/column order

    integer :: i, j, irow, icol, newcol, newrow

    !---------------------------------------------------------------
    ! Solver parameters
    ! TODO: Pass these in as arguments?
    !---------------------------------------------------------------

    real(dp), parameter ::   &
       tolerance = 1.d-11    ! tolerance for linear solver

    integer, parameter ::    &
       maxiters = 200         ! max number of linear iterations before quitting

    integer, parameter ::   &
       storage = 2            ! = 1 for compressed sparse column (CSC)
                              ! = 2 for compressed sparse row (CSR)
                              ! (See comments below for details)

    integer :: n

!WHL - debug
    integer, dimension(matrix%nonzeros) :: itemp

    if (present(verbose)) then
       verbose_pcg = verbose
    else
       verbose_pcg = .true.   ! for debugging
    endif
    
!WHL - debug
        if (verbose_pcg) then
           print*, ' '
           print*, 'Using standalone PCG solver'
           print*, 'matrix%order =', matrix%order
           print*, 'matrix%nonzeros =', matrix%nonzeros
           print*, 'size(rhs) =', size(rhs)
           print*, 'size(solution) =', size(solution)
           print*, 'size(row) =', size(matrix%row)
           print*, 'size(col) =', size(matrix%col)
           print*, 'size(val) =', size(matrix%val)
           print*, 'storage (1 = CSC, 2 = CSR) =', storage
        endif

    ! Make a local copy of the nonzero matrix entries.
    ! This is done so that the original matrix (in triad format) 
    !   can be passed as intent(in) and restored to the solver unchanged.
    ! Note: The arrays matrix%row/col/val may have size > matrix%nonzeros,
    !       with the extra entries all equal to zero.  Here we copy
    !       only the nonzero entries.

    do n = 1, matrix%nonzeros
       iA(n) = matrix%row(n)
       jA(n) = matrix%col(n)
       A(n)  = matrix%val(n)
    enddo

    !---------------------------------------------------------------
    ! Some notes on matrix storage: 
    !
    ! Here we convert from the triad (row, col, val) format to either the 
    ! compressed sparse column (CSC) or compressed sparse row (CSR) format.
    !
    ! In either case, the array A contains the nonzero matrix entries.
    !
    ! For CSC, the array iA contains the row for each entry, whereas jA holds
    ! the offsets into the iA and A arrays for the beginning of each column.
    ! That is, matrix entry jA(icol) is the first entry in column icol,
    ! and entry jA(icol+1)-1 is the last entry in column icol.
    ! We set ja(order+1) = nonzeros + 1.
    !
    ! For CSR, the array jA contains the column for each entry, whereas iA holds
    ! the offsets into the jA and A arrays for the beginning of each row.
    ! That is, matrix entry iA(irow) is the first entry in row irow,
    ! and entry iA(irow+1)-1 is the last entry in row irow.
    ! We set ia(order+1) = nonzeros + 1.  
    !
    ! For a given storage format, we first check to see if the triad entries
    ! are ordered correctly (in ascending column order for CSC, or in ascending
    ! row order for CSR).  If so, we simply redefine jA (for CSC) or iA (for CSR).
    ! If not, we use some SLAP code (subroutines DS2Y and QS2I1D, rewritten below 
    ! in Fortran90) to reorder the triad entries and convert to CSC or CSR.
    !
    ! The SLAP subroutines place the diagonal entries as the first entries 
    ! in each column, but this is not required by the PCG solver.
    !---------------------------------------------------------------

!WHL - debug - Uncomment to switch row and column indices for testing
!       itemp(:) = iA(:)
!       iA(:) = jA(:)
!       jA(:) = itemp(:)

    if (verbose_pcg) then
       print*, ' '
       print*, 'Triad index, iA, jA, A:'
       do i = 1, min(ndiagmax, matrix%nonzeros)
          print*, i, iA(i), jA(i), A(i)
       enddo
    endif

!TODO - Put the following code in a subroutine

    if (storage == 1) then    ! convert triad to compressed sparse column

       if (verbose_pcg) print*, 'Convert triad to CSC:'

       call convert_triad_to_compressed_column(matrix%order, matrix%nonzeros,  &
                                               iA,    jA,    A)      

       if (verbose_pcg) then
          print*, ' '
          print*, 'CSC format, i, iA, jA, A:'
          do i = 1, min(ndiagmax, matrix%nonzeros)
             print*, i, iA(i), jA(i), A(i)
          enddo
       endif
  
    else          ! convert triad to compressed sparse row

       if (verbose_pcg) print*, 'Convert triad to CSR:'

       call convert_triad_to_compressed_row(matrix%order, matrix%nonzeros,  &
                                            iA,    jA,    A)      

       if (verbose_pcg) then
          print*, ' '
          print*, 'CSR format, i, jA, iA, A:'
          do i = 1, min(ndiagmax, matrix%nonzeros)
             print*, i, jA(i), iA(i), A(i)
          enddo
       endif

    endif         ! storage (CSC or CSR)
 
    if (verbose_pcg) print*, 'Call pcg solver'

    call pcg_solver(matrix%order,      matrix%nonzeros,           &
                    iA,                jA,             A,         & 
                    rhs,               solution,                  & 
                    tolerance,         maxiters,       storage,   &
                    err,               niters)

  end subroutine cism_sparse_pcg_solve

!****************************************************************************

  subroutine pcg_solver(order,     nonzeros,                &
                        iA,        jA,         A,           &
                        b,         x,                       &
                        tolerance, maxiters,   storage,     &
                        err,       niters)

    !---------------------------------------------------------------
    !  This subroutine uses a preconditioned conjugate-gradient solver to
    !  solve the equation $Ax=b$, using a diagonal preconditioner.
    !  Convergence is checked every {\em ncheck} steps.
    !
    !  It is based on the barotropic solver in the POP ocean model 
    !  (author Phil Jones, LANL).  POP, however, uses a 9-point 
    !  stencil, whereas here the matrix $A$ is stored explicitly
    !  in compressed row or compressed column format.
    !---------------------------------------------------------------

    integer, intent(in) :: order         ! order (number of rows) of the matrix

    integer, intent(in) :: nonzeros      ! number of nonzeros in the matrix

    integer, dimension(nonzeros), intent(in) ::  &
       iA          ! for CSC: array holding the row for each nonzero value
                   ! for CSR: array holding offsets into jA and A for start of each row

    integer, dimension(nonzeros), intent(in) ::  &
       jA          ! for CSC: array holding offsets into iA and A for start of each column
                   ! for CSR: array holding the column for each nonzero value

    real(dp), dimension(nonzeros), intent(in) ::  &
       A           ! array of nonzero values

    real(dp), intent(in), dimension(order)  ::  &
       b                                 ! right-hand-side vector

    real(dp), intent(inout), dimension(order)  ::  &
       x                                 ! solution vector (initial guess on input,
                                         ! final solution on output)

    real(dp), intent(in) ::   &
       tolerance                         ! error tolerance for solution

    integer, intent(in) ::    &
       maxiters                          ! max number of iterations before quitting

    integer, intent(in) ::  &
       storage          ! = 1 for compressed sparse column (CSC)
                        ! = 2 for compressed sparse row (CSR)

    real(dp), intent(out) ::  &
       err                               ! error (L2 norm of residual) in final solution

    integer, intent(out) ::   &
       niters                            ! iterations needed to solution

    !---------------------------------------------------------------
    ! Local variables and parameters   
    !---------------------------------------------------------------

    integer ::  &
       i, m,              &!
       irow, icol,        &! row and column indices
       ibgn, iend          ! first and last entries in row or column

    real(dp) ::           &
       eta0, eta1, rr      ! scalar inner product results

    real(dp), dimension(order) ::  &
       Adiag,             &! diagonal terms of matrix A
       r,                 &! residual vector (b-Ax)
       d,                 &! conjugate direction vector
       y,                 &! result of a matvec multiply
       work0, work1        ! various cg intermediate results

    real(dp) ::  &
       L2_resid,          &! L2 norm of residual vector Ax-b
       L2_rhs              ! L2 norm of rhs vector b
                           ! solver converges when L2_resid/L2_target < tolerance

!    character (char_len) :: & 
!       noconvrg           ! error message for no convergence

    integer, parameter :: &
        solv_ncheck = 1    ! check for convergence every solv_ncheck iterations
                           ! for now, hardwire for debugging
                           ! later, let this be an argument

    integer, parameter :: &
        precond = 1        ! = 0 for no preconditioning
                           ! = 1 for diagonal preconditioning
                           ! = 2 for something better?

    ! Check for valid parameter values

    if (storage /= 1 .and. storage /= 2) then 
       print*, 'Invalid storage option; must be 1 (CSC) or 2 (CSR)'
       stop    !TODO - Abort cleanly
    endif

    if (verbose_pcg) then
       print*, ' '
       print*, 'In unstructured standalone pcg solver'
       print*, 'tolerance, maxiters =', tolerance, maxiters
       print*, 'L2_rhs =', L2_rhs
!       print*, ' '
!       print*, 'sum(x):', sum(x)
!       print*, 'max(x):', maxval(x)
!       print*, 'sum(b):', sum(b)
!       print*, 'max(b)',  maxval(b)
!       print*, ' '
!       print*, 'Call matvec, y = Ax'
    endif

    !  Compute initial residual and initialize the direction vector d

    call matvec_multiply(order,  nonzeros,  &
                         x,      y,         &
                         iA,     jA,    A,  &
                         storage)

    r(:) = b(:) - y(:)
    d(:) = 0.d0

    ! Initialize fields and scalars

    ! POP does a halo update here.
    ! To parallelize, we would need a halo data structure that keeps
    !  track of which rows of the matrix are owned by this processor,
    !  and which by other processors.
    ! Else we need to make sure that there is no duplication of
    !  rows across processors.

    ! call halo_update(r)

    eta0 = 1.d0 
    niters = maxiters 

    ! Compute L2 norm of rhs vector b

    work0(:) = b(:)*b(:)

    ! POP does a halo update here before taking the global sum
    ! (should this be done before computing b*b?)

    L2_rhs = sum(work0)         ! (b, b) = squared L2 norm
    ! In parallel, this would need to be a global sum: L2_rhs = global_sum(work0)

    L2_rhs = sqrt(L2_rhs)       ! L2 norm of rhs
                                !TODO - check that L2_rhs > 0

    if (precond == 1) then    ! Form diagonal matrix for preconditioning

       Adiag(:) = 0.d0

       if (storage == 1) then   ! CSC

          do icol = 1, order
             ibgn = jA(icol)
             iend = jA(icol+1) - 1
             do i = ibgn, iend
                if (iA(i) == icol) then  ! row = column
                   Adiag(icol) = A(i)
                   exit
                endif
            enddo 
          enddo

       elseif (storage == 2) then   ! CSR

          do irow = 1, order
             ibgn = iA(irow)
             iend = iA(irow+1) - 1
             do i = ibgn, iend
                if (jA(i) == irow) then  ! row = column
                   Adiag(irow) = A(i)
                   exit
                endif
             enddo 
          enddo

       endif   ! storage

    endif      ! precond = 1

    if (verbose_pcg) then
!       print*, ' '
!       print*, 'sum(Adiag) =', sum(Adiag)
       print*, ' '
       print*, 'diagonal entries, initial guess, residual:'
       do i = 1, min(ndiagmax,order)
          print*, i, Adiag(i), x(i), r(i)
       enddo
    endif

    ! Iterate to solution

    iter_loop: do m = 1, maxiters

       if (verbose_pcg) then
!          print*, ' '
!          print*, 'iter =', m
       endif

       ! Compute (PC)r using a diagonal preconditioner
       ! TODO: Consider using a different preconditioner if the
       !       diagonal version is inefficient
       !
       ! Could we use the SIA matrix as a preconditioner?
       ! I.e. solve Py = r for y, thus y = work1 = P^(-1)r
       ! Need an SIA solver that (generically) solves Ax = b.
       ! Might be easy to write since A is tridiagonal in every column.

       if (precond == 0) then      ! no preconditioning

           work1(:) = r(:)         ! PC(r) = r 

       elseif (precond ==1 ) then  ! diagonal preconditioning

          do icol = 1, order
             if (Adiag(icol) /= 0.d0) then
                work1(icol) = r(icol)/Adiag(icol)  ! This is PC(r), where PC is formed from diagonal elements of A
             else                                  ! All diagonal elements should be nonzero (TODO: check this)                                       
                work1(icol) = 0.d0
             endif
          enddo    ! icol

       endif    ! precond

       work0(:) = r(:)*work1(:)    ! terms of dot product (r, PC(r))

       ! POP has a halo update here for parallel solve with external preconditioner (lprecond = T)
       ! if (lprecond) call halo_update(work1)

       ! Compute (r, PC(r))

       ! In parallel, this would need to be a global sum: eta1 = global_sum(work0)
       eta1 = sum(work0)

       if (verbose_pcg) then
!          print*, 'eta0, eta1 =', eta0, eta1
       endif

       ! Update conjugate direction vector d

       d(:) = work1(:) + d(:)*(eta1/eta0)   ! d_(i+1) = PC(r_(i+1)) + beta_(i+1)*d_i
                                            !
                                            !                    (r_(i+1), PC(r_(i+1)))
                                            ! where beta_(i+1) = --------------------  
                                            !                        (r_i, PC(r_i)) 
                                            ! Initially eta0 = 1  
                                            ! For m >=2, eta0 = the old eta1

!       print*, 'work0 =', work0
!       print*, 'd =', d(:)

       if (verbose_pcg) then
!          print*, 'sum(work1) =', sum(work1)
!          print*, 'sum(d) =', sum(d)
!          print*, ' '
!          print*, 'Call matvec, y = Ad'
       endif

       ! Compute y = A*d

       call matvec_multiply(order,  nonzeros,  &
                            d,      y,         &
                            iA,     jA,    A,  &
                            storage)

       work0 = y(:) * d(:)       ! terms of dot product (d, Ad)

!       print*, 'Ad =', y(:)
!       print*, '(d, Ad) =', work0(:)

       ! Compute next solution and residual

       ! Here POP does a halo update for y = A*d
       ! call halo_update(y)

        eta0 = eta1               ! (r_(i+1), PC(r_(i+1))) --> (r_i, PC(r_i)) 

                                  !          (r, PC(r))
        eta1 = eta0/sum(work0)    ! alpha = ----------
                                  !          (d, A*d)


       if (verbose_pcg) then
!          print*, 'sum(y=Ad) =', sum(y)
!          print*, 'max(y=Ad) =', maxval(y)
!          print*, ' '
!          print*, 'sum(work0 = d*Ad) =', sum(work0)
!          print*, 'new eta0, eta1 =', eta0, eta1
       endif

        ! In parallel, this would need to be a global sum: eta1 = eta0/global_sum(work0)

!        print*, '(r, PC(r))=', eta0
!        print*, 'alpha =', eta1

        x(:) = x(:) + eta1 * d(:)    ! new solution, x_(i+1) = x_i + alpha*d
        r(:) = r(:) - eta1 * y(:)    ! new residual, r_(i+1) = r_i - alpha*(Ad)

        if (verbose_pcg) then
!           print*,'new sum(x) =', sum(x)
!           print*,'new max(x) =', maxval(x)
!           print*,'new sum(r) =', sum(r)
!           print*,'new max(r) =', maxval(r)
        endif

        if (verbose_pcg) then
!           print*, ' '
!           print*, 'New solution, residual:'
!           do i = 1, min(ndiagmax,order)
!              print*, x(i), r(i)
!           enddo
        endif

        if (mod(m,solv_ncheck) == 0) then

           if (verbose_pcg) then
!              print*, ' '
!              print*, 'Call matvec, y = Ax'
           endif

           ! Compute y = Ax
           
           call matvec_multiply(order,  nonzeros,  &
                                x,      y,         &
                                iA,     jA,    A,  &
                                storage)

           ! Compute residual and terms of its L2 norm
           ! TODO - How does this differ from the residual computed above?

           r(:) = b(:) - y(:)
           work0(:) = r(:)*r(:)

           ! POP does a halo update here (Why? Should it be done before computing work0?)
           ! call halo_update(r)
           
           L2_resid = sum(work0)         ! (r, r) = squared L2 norm
           ! In parallel, this would need to be a global sum: L2_resid = global_sum(work0)

           L2_resid = sqrt(L2_resid)     ! L2 norm of residual

           err = L2_resid/L2_rhs

           if (verbose_pcg) print*, 'iter, L2_resid, error =', m, L2_resid, err

           if (err < tolerance) then
              niters = m
              exit iter_loop
           endif            

        endif    ! solv_ncheck

     enddo iter_loop

!WHL - Without good preconditioning, convergence is slow, but the solution after maxiters may be good enough.
!      Diagonal preconditioning is not good for the cases I've tested so far (e.g., dome).
 
     if (niters == maxiters) then
        print*, 'Glissade PCG solver not converged'
        print*, 'niters, err, tolerance:', niters, err, tolerance
!!!        stop   ! TODO - Abort cleanly

     endif

  end subroutine pcg_solver

!****************************************************************************

  subroutine matvec_multiply(order, nonzeros,    &
                             x,     y,           &
                             iA,    jA,   A,     &
                             matrix_storage)

    !---------------------------------------------------------------
    ! Compute the matrix-vector product $y = Ax$.
    !
    ! A large fraction of the computational time is spent in this subroutine,
    ! so efficiency is important.
    ! 
    ! The input matrix can be in one of two formats: compressed sparse
    !  column (CSC) or compressed sparse row (CSR).
    ! In serial, these should be about equally efficient.
    ! CSR is better for fine-grained parallelism because each row
    !  of y can be computed independently.
    !---------------------------------------------------------------

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::  &
       order,     &! order of the matrix (number of rows)
       nonzeros    ! number of nonzero values

    real(dp), dimension(order), intent(in) ::  &
       x           ! vector to be multiplied by the matrix

    real(dp), dimension(order), intent(out) ::  &
       y           ! output vector A*x

    integer, dimension(nonzeros), intent(in) ::  &
       iA          ! for CSC: array holding the row for each nonzero value
                   ! for CSR: array holding offsets into jA and A for start of each row

    integer, dimension(nonzeros), intent(in) ::  &
       jA          ! for CSC: array holding offsets into iA and A for start of each column
                   ! for CSR: array holding the column for each nonzero value

    real(dp), dimension(nonzeros), intent(in) ::  &
       A           ! array of nonzero values

    integer, intent(in), optional ::  &
       matrix_storage     ! = 1 for compressed sparse column storage (CSC, as in SLAP)
                          ! = 2 for compressed sparse row storage (CSR, as in Trilinos)

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: icol         ! column index
    integer :: irow         ! row index

    integer :: ibgn, iend   ! starting and ending indices for a given row or column

    integer :: i            ! intermediate indices from ibgn to iend

    integer :: storage      ! = 1 for CSC, = 2 for CSR


!WHL - debug
!!    integer :: rowtest
!!    rowtest = 2*ntest - 1
    
    if (verbose_pcg) then
!       print*, ' '
!       print*, 'In unstructured matvec'
!       print*, 'max(A) =', maxval(A)
!       print*, 'sum(x):', sum(x)
!       print*, 'max(x):', maxval(x)
!!       print*, 'rowtest, x(rowtest):', rowtest, x(rowtest) 
    endif

    if (present(matrix_storage)) then
       storage = matrix_storage   ! 1 for CSC, 2 for CSR
    else
       storage = 1                ! CSC by default
    endif


    ! Initialize the result vector.

    y(:) = 0.d0

    ! Compute y = Ax

    if (storage == 2) then   ! compressed sparse row

       do irow = 1, order
          ibgn = iA(irow)       ! index of first nonzero in irow
          iend = iA(irow+1)-1   ! index of last nonzero in irow
          do i = ibgn, iend
             y(irow) = y(irow) + A(i)*x(jA(i))

!WHL- debug
!!             if (irow == rowtest) then
!!                print*, 'row, A, x, Ax, ynew:', irow, A(i), x(jA(i)), A(i)*x(jA(i)), y(irow)
!!             endif

          enddo
       enddo

    else                     ! compressed sparse column

       do icol = 1, order
          ibgn = jA(icol)       ! index of first nonzero in icol 
          iend = jA(icol+1)-1   ! index of last nonzero in icol 
          do i = ibgn, iend
             y(iA(i)) = y(iA(i)) + A(i)*x(icol)
          enddo
       enddo

    endif

  end subroutine matvec_multiply
 
!****************************************************************************

  subroutine convert_triad_to_compressed_column(order, nonzeros,  &
                                                iA,    jA,    A)      

    !---------------------------------------------------------------
    ! input-output variables    
    !---------------------------------------------------------------

    integer, intent(in) :: &
       order,       & ! order of the matrix (number of rows)
       nonzeros       ! number of nonzero entries

    integer, dimension(nonzeros), intent(inout) ::  &
       iA             ! row of each nonzero entry

    integer, dimension(nonzeros), intent(inout) ::  &
       jA             ! initially refers to the column of each nonzero entry,
                      ! but for CSC storage is converted to an array that holds offsets
                      ! into the A and iA arrays for the start of each row

    real(dp), dimension(nonzeros), intent(inout) ::   &
       A              ! nonzero matrix entries

    !---------------------------------------------------------------
    ! local variables    
    !---------------------------------------------------------------

    logical ::   &
       sorted         ! true if triad matrix entries are in ascending column order

    integer :: i, j, icol, newcol

    ! check whether triad data are already in ascending column order

    sorted = .true.
    icol = jA(1)
    do i = 2, nonzeros      ! loop through all nonzero entries
       newcol = jA(i)           
       if (newcol < icol) then
          sorted = .false.
          exit
       elseif (newcol > icol) then 
          icol = newcol
       endif
    enddo

    if (sorted) then      ! convert to CSC

       if (verbose_pcg) print*, 'Already in ascending column order; redefine jA'
       
       jA(1) = 1
       do icol = 1, order-1
          do j = jA(icol)+1, nonzeros
             if (jA(j) /= icol) then
                jA(icol+1) = j
                exit
             endif
          enddo
       enddo
       jA(order+1) = nonzeros + 1

    else      ! use SLAP subroutines to reorder and convert to CSR

       if (verbose_pcg) print*, 'Not sorted; convert to SLAP column'

       call slap_convert_triad(order, nonzeros,    &
                               iA,    jA,    A)      

    endif     ! sorted


  end subroutine convert_triad_to_compressed_column

!****************************************************************************

  subroutine convert_triad_to_compressed_row(order, nonzeros,  &
                                             iA,    jA,    A)      

    !---------------------------------------------------------------
    ! input-output variables    
    !---------------------------------------------------------------

    integer, intent(in) :: &
       order,       & ! order of the matrix (number of rows)
       nonzeros       ! number of nonzero entries

    integer, dimension(nonzeros), intent(inout) ::  &
       iA             ! initally refers to the row of each nonzero entry,
                      ! but for CSR storage is converted to an array that holds offsets
                      ! into the A and jA arrays for the start of each row

    integer, dimension(nonzeros), intent(inout) ::  &
       jA             ! column of each nonzero entry

    real(dp), dimension(nonzeros), intent(inout) ::   &
       A              ! nonzero matrix entries

    !---------------------------------------------------------------
    ! local variables    
    !---------------------------------------------------------------

    logical ::   &
       sorted         ! true if triad matrix entries are in ascending row order

    integer :: i, j, irow, newrow

       ! check whether triad data are already in ascending row order

       sorted = .true.
       irow = iA(1)
       do i = 2, nonzeros   ! loop through all nonzero entries
          newrow = iA(i)           
          if (newrow < irow) then
             sorted = .false.
             exit
          elseif (newrow > irow) then 
             irow = newrow
          endif
       enddo

       if (sorted) then      ! convert to CSR

          if (verbose_pcg) print*, 'Already in ascending row order; redefine iA'

          iA(1) = 1
          do irow = 1, order-1
             do i = iA(irow)+1, nonzeros
                if (iA(i) /= irow) then
                   iA(irow+1) = i
                   exit
                endif
            enddo
          enddo
          iA(order+1) = nonzeros + 1

       else      ! use SLAP subroutines to reorder and convert to CSR

          if (verbose_pcg) print*, 'Not sorted; convert to SLAP row'

          ! Note: jA is passed in first slot and iA in second slot
          !       for conversion to CSR

          call slap_convert_triad(order, nonzeros,  &
                                  jA,    iA,    A)

       endif     ! sorted


  end subroutine convert_triad_to_compressed_row

!****************************************************************************

!TODO - Clean up this subroutine

  subroutine slap_convert_triad(order, nonzeros,   &
                                iA,    jA,    A)         

    !---------------------------------------------------------------
    ! Convert sparse matrix storage format from triad (row, col, val) to
    ! compressed sparse column format for more efficient matrix operations.
    !
    ! By passing the jA array in the iA slot and passing the iA array
    ! in the jA slot, this subroutine can be used to convert from 
    ! triad to compressed sparse row format.
    !
    ! Based on the SLAP subroutine DS2Y (author Mark K. Seager, LLNL)
    !---------------------------------------------------------------
 
    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) :: &
       order,          & ! order of the matrix (number of rows)
       nonzeros          ! number of nonzero entries

    integer, dimension(nonzeros), intent(inout) ::  &
       iA                ! row of each nonzero entry

    integer, dimension(nonzeros), intent(inout) ::  &
       jA                ! initially refers to the column of each nonzero entry,
                         ! but converted to an array that holds offsets into
                         !  the A and iA arrays for the start of each column

    real(dp), dimension(nonzeros) ::   &
       A                 ! nonzero matrix entries

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer :: i, j, icol, ibgn, iend, itemp

    real(dp) :: temp
  
    if (verbose_pcg) then
       print*, ' '
       print*, 'Convert triad' 
       print*, 'Initial i, iA, jA, A:'
       do i = 1, min(ndiagmax, nonzeros)
          print*, i, iA(i), jA(i), A(i)
        enddo
    endif

    ! Check to see if the (iA,jA,A) arrays are in SLAP column 
    ! format.  If not, then transform from SLAP Triad.

! SLAP code:
!      IF( JA(N+1).EQ.NELT+1 ) RETURN

! New code:
    if (jA(order+1) == nonzeros+1) return

    ! Sort into ascending order by column using a quicksort algorithm.

! SLAP code:
!      CALL QS2I1D( JA, IA, A, NELT, 1 )

    if (verbose_pcg) print*, 'Call triad quicksort'

    ! Note that jA appears before iA in the argument list
    ! Switching jA with iA would result in a row quicksort

     call triad_quicksort(nonzeros, jA, iA, A)    ! column quicksort (jA first)

    ! Use this for a row quicksort
!!     call triad_quicksort(nonzeros, iA, jA, A)    ! row quicksort (iA first)


     if (verbose_pcg) then
        print*, ' '
        print*, 'After triad quicksort:'
        print*, 'index, iA, jA, A:'
        do i = 1, min(ndiagmax, nonzeros)
           print*, i, iA(i), jA(i), A(i)
        enddo
    endif

    ! Loop over each column to see where the column indices change 
    ! in the column index array jA.  This marks the beginning of the
    ! next column.
         
! SLAP code:
!      CVD$R NOVECTOR
!      JA(1) = 1
!      DO 20 ICOL = 1, N-1
!         DO 10 J = JA(ICOL)+1, NELT
!            IF( JA(J).NE.ICOL ) THEN
!               JA(ICOL+1) = J
!               GOTO 20
!            ENDIF
! 10      CONTINUE
! 20   CONTINUE
!      JA(N+1) = NELT+1

! New code:
    jA(1) = 1
    do icol = 1, order-1
       do j = jA(icol)+1, nonzeros
          if (jA(j) /= icol) then
             jA(icol+1) = j
             exit  ! make sure this only exits the inner loop
          endif
       enddo
    enddo
    jA(order+1) = nonzeros + 1
    

!WHL - debug
    if (verbose_pcg) then
       print*, 'icol, jA(icol):'
       do icol = 1, min(ndiagmax,order+1)
         print*, icol, jA(icol)
       enddo
    endif

!   Mark the order+2 element so that future calls to a subroutine 
!   using the SLAP Column storage format will be able to tell.
!TODO - Not necessary?

! SLAP code:        
!      JA(N+2) = 0

! New code:
    jA(order+2) = 0

!   Loop thru the iA(i) array making sure that the diagonal
!   matrix element appears first in the column.  Then sort the
!   rest of the column in ascending order.

! SLAP code:
!      DO 70 ICOL = 1, N
!         IBGN = JA(ICOL)
!         IEND = JA(ICOL+1)-1
!         DO 30 I = IBGN, IEND
!            IF( IA(I).EQ.ICOL ) THEN
!C         Swap the diag element with the first element in the column.
!               ITEMP = IA(I)
!               IA(I) = IA(IBGN)
!               IA(IBGN) = ITEMP
!               TEMP = A(I)
!               A(I) = A(IBGN)
!               A(IBGN) = TEMP
!               GOTO 40
!            ENDIF
! 30      CONTINUE
! 40      IBGN = IBGN + 1
!         IF( IBGN.LT.IEND ) THEN
!            DO 60 I = IBGN, IEND
!               DO 50 J = I+1, IEND
!                  IF( IA(I).GT.IA(J) ) THEN
!                     ITEMP = IA(I)
!                     IA(I) = IA(J)
!                     IA(J) = ITEMP
!                     TEMP = A(I)
!                     A(I) = A(J)
!                     A(J) = TEMP
!                  ENDIF
! 50            CONTINUE
! 60         CONTINUE
!         ENDIF
! 70   CONTINUE
!      RETURN

! New code:
    do icol = 1, order
       ibgn = jA(icol)
       iend = jA(icol+1) - 1
       do i = ibgn, iend
          if (iA(i) == icol) then  ! swap diagonal element with first element in column
             itemp = iA(i)
             iA(i) = iA(ibgn)
             iA(ibgn) = itemp
             temp = A(i)
             A(i) = A(ibgn)
             A(ibgn) = temp
             exit  ! TODO: make sure we exit only the inner loop
          endif
       enddo
       ibgn = ibgn + 1
       if (ibgn < iend) then
          do i = ibgn, iend
             do j = i+1, iend
                if (iA(i) > iA(j)) then
                   itemp = iA(i)
                   iA(i) = iA(j)
                   iA(j) = itemp
                   temp = A(i)
                   A(i) = A(j)
                   A(j) = temp
                endif
             enddo   ! j
          enddo      ! i
       endif        ! ibgn < iend
    enddo           ! icol

    if (verbose_pcg) then
       print*, ' '
       print*, 'After diagonal swap:'
       print*, 'index, iA, jA, A:'
       do i = 1, min(ndiagmax, nonzeros)
          print*, i, iA(i), jA(i), A(i)
       enddo
    endif

  end subroutine slap_convert_triad

!****************************************************************************

  subroutine triad_quicksort(nonzeros, iA, jA, A, kflag)

    !---------------------------------------------------------------
    ! Sort the iA array such that matrix entries are listed 
    ! in either ascending or descending order.  
    ! Make the same interchanges in the jA and A arrays.
    !
    ! Note: If the incoming iA array is the row array, matrix entries
    !       will be sorted by rows.  If the incoming iA array is the
    !       column array, matrix entries will be sorted by columns.
    !
    ! This is a Fortran90 rewrite of SLAP subroutine QS2I1D. 
    ! Authors: R.E. Jones (SNLA), D.K. Kahaner (NBS), M.K. Seager (LLNL),
    !          J.A. Wisniewski (SNLA)
    ! Based on this reference:
    !     Singleton, R. C., Algorithm 347, "An Efficient Algorithm for 
    !     Sorting with Minimal Storage", cacm, Vol. 12, No. 3, 1969, 185-187.
    !---------------------------------------------------------------

    ! TODO: Clean up this subroutine.
    ! This subroutine has many "go to" statements that make it hard to follow the logic.
    !  I'd like to remove them, but this would take some work, so I'm leaving them for now.

    !---------------------------------------------------------------
    ! input-output arguments
    !---------------------------------------------------------------

    integer, intent(in) ::  &
       nonzeros    ! number of values to be sorted

    integer, dimension(nonzeros), intent(inout) ::  &
       iA,      & ! integer array to be sorted
       jA         ! integer array to be carried along

    real(dp), dimension(nonzeros), intent(inout) ::  &
       A          ! real(dp) array to be carried along

    integer, intent(in), optional ::   &
       kflag      ! =  1 for sort in ascending order
                  ! = -1 for sort in descending order

    !---------------------------------------------------------------
    ! local variables
    !---------------------------------------------------------------

    integer, dimension(21) :: il, iu   ! why 21?

    integer :: nn      ! number of nonzeros
    integer :: i, k    ! indices for first part of array
    integer :: j, l    ! indices for second part of array
    integer :: ij      ! index starting somewhere in middle of array
    integer :: kk, m
    integer :: it, jt, iit, jjt  ! temporary values of iA, jA

!TODO - Change it to itmp, etc.?

    real(dp) :: tA, ttA          ! temporary values of A

    real(dp) :: r      ! magic number for quick sorting

      ! Check for valid values of nonzeros and kflag

! SLAP code:
!
!***FIRST EXECUTABLE STATEMENT  QS2I1D
!      NN = N
!      IF (NN.LT.1) THEN
!         CALL XERROR ( 'QS2I1D- the number of values to be sorted was no
!     $T POSITIVE.',59,1,1)
!         RETURN
!      ENDIF
!      IF( N.EQ.1 ) RETURN
!      KK = IABS(KFLAG)
!      IF ( KK.NE.1 ) THEN
!         CALL XERROR ( 'QS2I1D- the sort control parameter, k, was not 1
!     $ OR -1.',55,2,1)
!         RETURN
!      ENDIF

!      IF( KFLAG.LT.1 ) THEN
!         DO 20 I=1,NN
!            IA(I) = -IA(I)
! 20      CONTINUE
!      ENDIF

! New code:

      nn = nonzeros

      print*, 'Start of triad_quicksort'
      print*, 'i, iA, jA, A:'
      do i = 1, min(ndiagmax,nonzeros)
         print*, i, iA(i), jA(i), A(i)
      enddo


      if (nn < 1) then
         print*, 'Triad quicksort: must have nonzeros >= 1'
         stop    ! TODO - clean abort
      elseif (nn == 1 ) then
         return
      endif

      if (present(kflag)) then

         kk = iabs(kflag)
         if (kk /= 1) then
            print*, 'Triad quicksort: must have kflag = 1 or -1'
            stop    !TODO - clean abort
         endif

         ! Flip sign of rowA to get decreasing order if requested.
         if (kflag < 1) then
            do i = 1, nn
               iA(i) = -iA(i)
            enddo
         endif

      endif  ! present(kflag)

!     Sort iA, carrying jA and A along.

!     First a little black magic...
!     WHL: I don't know where these constants come from.
!          Maybe they nudge the algorithm toward optimal N*logN scaling?

! SLAP code:

!      M = 1
!      I = 1
!      J = NN
!      R = .375
! 210  IF( R.LE.0.5898437 ) THEN
!         R = R + 3.90625E-2
!      ELSE
!         R = R-.21875
!      ENDIF
! 225  K = I

! New code:

      m = 1
      i = 1
      j = nn  

      r = 0.375d0

 210  if (r <= 0.5898437d0) then
         r = r + 3.90625d-2
      else
         r = r - 0.21875d0
      endif

 225  k = i

!      print*, ' '
!      print*, 'r, i, j, k, m =', r, i, j, k, m

!     Select a central element of the array and save it in location 
!     it, jt, at.

! SLAP code:
!      IJ = I + IDINT( DBLE(J-I)*R )
!      IT = IA(IJ)
!      JT = JA(IJ)
!      TA = A(IJ)

! New code:
      ij = i + int((j-i)*r)   !TODO - Check that int corresponds to idint
      it = iA(ij)
      jt = jA(ij)
      tA =  A(ij)

!      print*, 'Choose central element, ij, i, j, A =', ij, it, jt, tA

      ! If first element of array is greater than it, interchange with it.

! SLAP code:
!      IF( IA(I).GT.IT ) THEN
!         IA(IJ) = IA(I)
!         IA(I)  = IT
!         IT     = IA(IJ)
!         JA(IJ) = JA(I)
!         JA(I)  = JT
!         JT     = JA(IJ)
!         A(IJ)  = A(I)
!         A(I)   = TA
!         TA     = A(IJ)
!      ENDIF
!      L=J

! New code:

      if (iA(i) > it) then
         iA(ij) = iA(i)
         iA(i)  = it
         it     = iA(ij)
         jA(ij) = jA(i)
         jA(i)  = jt
         jt     = jA(ij)
         A(ij)  = A(i)
         A(i)   = tA
         tA     = A(ij)
      endif
      l = j
                           
!     If last element of array is less than it, swap with it.

!SLAP code:
!      IF( IA(J).LT.IT ) THEN
!         IA(IJ) = IA(J)
!         IA(J)  = IT
1         IT     = IA(IJ)
!         JA(IJ) = JA(J)
!         JA(J)  = JT
!         JT     = JA(IJ)
!         A(IJ)  = A(J)
!         A(J)   = TA
!         TA     = A(IJ)

! New code:
      if (iA(j) < it) then
         iA(ij) = iA(j)
         iA(j)  = it
         it     = iA(ij)
         jA(ij) = jA(j)
         jA(j)  = jt          
         jt     = jA(ij)
         A(ij)  = A(j)
         A(j)   = tA
         tA     = A(ij)

         ! If first element of array is greater than it, swap with it.

! SLAP code:
!         IF ( IA(I).GT.IT ) THEN
!            IA(IJ) = IA(I)
!            IA(I)  = IT
!            IT     = IA(IJ)
!            JA(IJ) = JA(I)
!            JA(I)  = JT
!            JT     = JA(IJ)
!            A(IJ)  = A(I)
!            A(I)   = TA
!            TA     = A(IJ)
!         ENDIF
!      ENDIF

! New code:
         if (iA(i) > it) then
            iA(ij) = iA(i)
            iA(i)  = it
            it     = iA(ij)
            jA(ij) = jA(i)
            jA(i)  = jt
            jt     = jA(ij)
            A(ij)  = A(i)
            A(i)   = tA
            tA     = A(ij)
         endif   ! iA(i) > it

      endif      ! iA(j) < it
 

!     Find an element in the second half of the array which is 
!     smaller than it.

! SLAP code:
!  240 L=L-1
!      IF( IA(L).GT.IT ) GO TO 240

! New code:
 240  l = l-1    ! TODO - Get rid of goto's
      if (ia(l) > it) go to 240

!      Find an element in the first half of the array which is 
!      greater than it.

! SLAP code:
!  245 K=K+1
!      IF( IA(K).LT.IT ) GO TO 245

! New code:
 245  k = k+1
      if (ia(k) < it) go to 245

!      Interchange these elements.

! SLAP code:
!      IF( K.LE.L ) THEN
!         IIT   = IA(L)
!         IA(L) = IA(K)
!         IA(K) = IIT
!         JJT   = JA(L)
!         JA(L) = JA(K)
!         JA(K) = JJT
!         TTA   = A(L)
!         A(L)  = A(K)
!         A(K)  = TTA
!         GOTO 240
!      ENDIF

! New code:
       if (k <= l) then   ! lower-case "L", not 1

          iit = iA(l)
          iA(l) = iA(k)
          ia(k) = iit
          jjt = jA(l)
          jA(l) = jA(k)
          jA(k) = jjt
          ttA = A(l)
          A(l) = A(k)
          A(k) = ttA
          go to 240   ! Find element in 2nd half smaller than it, and in 1st half greater than it
       endif

!      Save upper and lower subscripts of the array yet to be sorted.

! SLAP code:
!      IF( L-I.GT.J-K ) THEN
!         IL(M) = I
!         IU(M) = L
!         I = K
!         M = M+1
!      ELSE
!         IL(M) = K
!         IU(M) = J
!         J = L
!         M = M+1
!      ENDIF
!      GO TO 260

! New code:
      if (l-i > j-k) then
         il(m) = i
         iu(m) = l
         i = k
         m = m + 1
      else
         il(m) = k
         iu(m) = j
         j = l
         m = m + 1
      endif
      go to 260

!      Begin again on another portion of the unsorted array.
                                  
! SLAP code:        
!  255 M = M-1
!      IF( M.EQ.0 ) GO TO 300
!      I = IL(M)
!      J = IU(M)
!  260 IF( J-I.GE.1 ) GO TO 225
!      IF( I.EQ.J ) GO TO 255
!      IF( I.EQ.1 ) GO TO 210
!      I = I-1
!  265 I = I+1
!      IF( I.EQ.J ) GO TO 255
!      IT = IA(I+1)
!      JT = JA(I+1)
!      TA =  A(I+1)
!      IF( IA(I).LE.IT ) GO TO 265
!      K=I
!  270 IA(K+1) = IA(K)
!      JA(K+1) = JA(K)
!      A(K+1)  =  A(K)
!      K = K-1
!      IF( IT.LT.IA(K) ) GO TO 270
!      IA(K+1) = IT
!      JA(K+1) = JT
!      A(K+1)  = TA
!      GO TO 265
 
! New code:
 255  m = m - 1

!      print*, 'Begin again, m =', m

      if (m == 0) go to 300    ! we are done  (how do we know this?)
      i = il(m)
      j = iu(m)
 260  if (j-i >= 1) go to 225  ! set k = i
      if (i == j) go to 255    ! set m = m+1
      if (i == 1) go to 210    ! change r
      i = i - 1       
 265  i = i + 1
      if (i == j) go to 255    ! set m = m+1
      it = iA(i+1)
      jt = jA(i+1)
      tA =  A(i+1)
      if (iA(i) <= it) go to 265   ! set i = i+1
      k = i
 270  iA(k+1) = iA(k)
      jA(k+1) = jA(k)
      A(k+1)  =  A(k)
      k = k - 1
      if (it < iA(k)) go to 270
      iA(k+1) = it
      jA(k+1) = jt
      A(k+1)  = tA

!      print*, 'go to 265, i =', i

      go to 265

      ! Clean up, if necessary.

! SLAP code:
!  300 IF( KFLAG.LT.1 ) THEN
!         DO 310 I=1,NN
!            IA(I) = -IA(I)
! 310     CONTINUE
!      ENDIF


 300  if (present(kflag)) then
         if (kflag < 1) then
            do i = 1, nn
               iA(i) = -iA(i)
            enddo
         endif
      endif   ! present(kflag)

     print*, 'End of triad_quicksort'
     print*, 'i, iA, jA, A:'
     do i = 1, min(ndiagmax, NN)
        print*, i, iA(i), jA(i), A(i)
     enddo


  end subroutine triad_quicksort

!****************************************************************************

end module cism_sparse_pcg
