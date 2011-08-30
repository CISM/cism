!  uvec is either u^k-1 or v^k-1 on input and Av-b or Cu-d on output

subroutine res_vect ( matrix, uvec, bvec, nu, g_flag, L2square, whatsparse)

use parallel

use glimmer_paramets, only : dp
use glimmer_sparse_type
use glimmer_sparse
use glide_mask

implicit none

integer :: i, j, nu, nele, whatsparse ! nu: size of uvec and bvec
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

      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then

        do nele = 1, matrix%nonzeros 

           i = matrix%row(nele)
           j = matrix%col(nele)
           Au_b_wig(i) = Au_b_wig(i) + matrix%val(nele) * uvec(j)

        enddo

      else 
        call matvecwithtrilinos(uvec, Au_b_wig);
      endif 

      do i = 1, nu
         Au_b_wig(i) = Au_b_wig(i) - bvec(i)
      enddo

      uvec = Au_b_wig
!      print*, "matrix%val,uvec",maxval(matrix%val), maxval(uvec)

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

      !JEFF Sum L2square across nodes
      L2square = parallel_reduce_sum(L2square)

      return

end subroutine res_vect

subroutine res_vect_jfnk ( matrixA, matrixC, uvec, bvec, nu1, nu2, g_flag, L2square, whatsparse)

! similar to res_vect, but state vector uvec and rhs vector bvec are now both velocities (kje 101005)
! as an intermediate step, right now A and C matrices are separate, but eventually they will be combined

use glimmer_paramets, only : dp
use glimmer_sparse_type
use glimmer_sparse
use glide_mask

implicit none

integer :: i, j, nu1, nu2, nele, whatsparse ! nu2: size of uvec and bvec, size of u, v within
!integer, dimension(:), intent(in) :: g_flag ! 0 :reg cell
                                             ! 1 :top ghost, 2 :base ghost

type(sparse_matrix_type),  intent(in) :: matrixA, matrixC

integer, dimension(nu2) :: g_flag 
real (kind = dp), dimension(nu2), intent(in) :: bvec
real (kind = dp), dimension(nu2), intent(inout) :: uvec
real (kind = dp), dimension(nu1) :: Au_b_wig, Cv_d_wig
real (kind = dp), intent(out) :: L2square
! 
real (kind = dp) :: scale_ghosts = 0.0d0

! calculate residual vector of the u and v component


      Au_b_wig = 0d0 ! regular+ghost cells
      Cv_d_wig = 0d0 ! regular+ghost cells

      if (whatsparse /= STANDALONE_TRILINOS_SOLVER) then

        do nele = 1, matrixA%nonzeros

           i = matrixA%row(nele)
           j = matrixA%col(nele)
           Au_b_wig(i) = Au_b_wig(i) + matrixA%val(nele) * uvec(j)

        enddo

        do nele = 1, matrixC%nonzeros

           i = matrixC%row(nele)
           j = matrixC%col(nele)
           Cv_d_wig(i) = Cv_d_wig(i) + matrixC%val(nele) * uvec(nu1+j)

        enddo

      else

        call matvecwithtrilinos(uvec(1:nu1), Au_b_wig);
        call matvecwithtrilinos(uvec(nu1+1:nu2), Cv_d_wig);
      endif

      do i = 1, nu1

         Au_b_wig(i) = Au_b_wig(i) - bvec(i)
         Cv_d_wig(i) = Cv_d_wig(i) - bvec(nu1+i)

      enddo

! to do: combine A and C

      do i = 1, nu1

         uvec(i)    = Au_b_wig(i)
         uvec(nu1+i) = Cv_d_wig(i)

      enddo
! AGS: Residual norm includes scaling to decrease importance of ghost values
! By calling it a redefinition of an inner product, it is kosher.
!      L2square = 0.0
!      do i = 1, nu1
!         if (g_flag(i) .eq. 0) then
!            L2square = L2square + Au_b_wig(i) * Au_b_wig(i)
!         else
!            L2square = L2square + scale_ghosts * Au_b_wig(i) * Au_b_wig(i)
!         endif
!      end do
!
!      do i = 1, nu1
!         if (g_flag(nu1+i) .eq. 0) then
!            L2square = L2square + Cv_d_wig(i) * Cv_d_wig(i)
!         else
!            L2square = L2square + scale_ghosts * Cv_d_wig(i) * Cv_d_wig(i)
!         endif
!      end do
! when the combined version is used, convergence wrong
      do i = 1, nu2
         if (g_flag(i) .eq. 0) then
            L2square = L2square + uvec(i) * uvec(i)
         else
            L2square = L2square + scale_ghosts * uvec(i) * uvec(i)
         endif
      end do


      return

end subroutine res_vect_jfnk
