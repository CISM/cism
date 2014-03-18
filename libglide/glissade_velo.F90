!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!
!TODO - Are all these includes needed?
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

!TODO - What is this?
#define shapedbg(x) write(*,*) "x", shape(x)

!WHL - This module used to contain calls to both the glam and glissade velo solvers. 
!      I moved the glam part to a separate module (glam_velo.F90), leaving just the
!      glissade part in this module.

module glissade_velo

    use parallel

    ! Driver for glam and glissade higher-order velocity solvers

    implicit none
    
contains
        
    subroutine glissade_velo_driver(model)

      ! Glissade higher-order velocity driver

      use glimmer_global, only : dp
      use glimmer_physcon, only: gn, scyr
      use glimmer_paramets, only: thk0, len0, vel0, vis0, tau0, evs0
      use glimmer_log
      use glide_types
      use glissade_velo_higher, only: glissade_velo_higher_solve
      use glissade_basal_traction, only: calcbeta
      use glide_mask

      type(glide_global_type),intent(inout) :: model


      !WHL - temporary velocity arrays to remove scaling
      real(dp), dimension(model%general%upn, model%general%ewn-1, model%general%nsn-1) ::  &
           uvel, vvel    ! uvel and vvel with scale factor (vel0) removed

      !TODO - Replace with a different mask?
      !       This is used for now as an input to subroutine calcbeta.
      integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag

!WHL - debug
      integer :: i, j        

      real(dp) :: maxbeta, minbeta
      !-------------------------------------------------------------------
      ! Compute masks
      ! These masks will not change even if we iterate on the velocity
      !  multiple times while solving the thickness transport equation.
      !-------------------------------------------------------------------

      !TODO - This subroutine is used for now to compute geom_mask_stag,
      !       an input to calcbeta.  Might replace later.

      call glide_set_mask(model%numerics,                                     &
                          model%geomderv%stagthck, model%geomderv%stagtopg,   &
                          model%general%ewn-1,     model%general%nsn-1,       &
                          model%climate%eus,       geom_mask_stag)

      !-------------------------------------------------------------------
      ! Compute or prescribe the basal traction field 'beta'.                                                              
      !
      ! Note: The initial value of model%velocity%beta can change depending on
      !       the value of model%options%which_ho_babc.     
      !-------------------------------------------------------------------

      call calcbeta (model%options%which_ho_babc,                  & 
                     model%numerics%dew,      model%numerics%dns,  &
                     model%general%ewn,       model%general%nsn,   &
                     model%velocity%uvel(model%general%upn,:,:),   &
                     model%velocity%vvel(model%general%upn,:,:),   &
                     model%temper%bwat,                            &
                     model%paramets%ho_beta_const,                 &
                     model%basalproc%mintauf,                      &
                     geom_mask_stag,                               &
                     model%velocity%beta)

      call staggered_parallel_halo(model%velocity%beta)

      maxbeta = maxval(model%velocity%beta)
      maxbeta = parallel_reduce_max(maxbeta)
      minbeta = minval(model%velocity%beta)
      minbeta = -parallel_reduce_max(-minbeta)

      !WHL - debug
      if (this_rank==model%numerics%rdiag_local) then
!!         print*, ' '
!!         print*, 'task, which_ho_babc =', this_rank, model%options%which_ho_babc
!!         print*, 'After calcbeta: max, min of beta (Pa/(m/yr)) =', &
!!              tau0/vel0/scyr * maxbeta, tau0/vel0/scyr * minbeta
      endif

      !----------------------------------------------------------------
      ! Note: The glissade solver uses SI units.
      ! Thus we have grid cell dimensions and ice thickness in meters,
      !  velocity in m/s, and the rate factor in Pa^(-n) s(-1).
      !----------------------------------------------------------------

      if (model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

         ! Note: The geometry fields (thck, topg, and usrf) must be updated in halos
         !        before calling glissade_velo_higher-solver.
         !       These updates are done in subroutine glissade_diagnostic_variable_solve
         !        in module glissade.F90.

            
!WHL - debug
         if (this_rank==model%numerics%rdiag_local) then
!!            print*, ' '
!!            print*, 'Call glissade_velo_higher_solve, task', this_rank
!!            print*, 'nx, ny =', model%general%ewn, model%general%nsn
!!            print*, 'size(kinbcmask) =', size(model%velocity%kinbcmask,1), size(model%velocity%kinbcmask,2)
!!            print*, ' '
!!            print*, 'kinbcmask before halo update:'
            do j = model%general%nsn-1, 1, -1
               do i = 1, model%general%ewn-1
!!                  write(6,'(46i3)',advance='no') model%velocity%kinbcmask(i,j)
               enddo
!!               write(6,*) ' '
            enddo
         endif    ! main_task

         !WHL - Instead of assuming that kinbcmask is periodic, extrapolate
         !       the kinbcmask into the global halo region
         !       (and also into the north and east rows of the global domain,
         !       which are not included on the global staggered grid).

         call staggered_parallel_halo_extrapolate (model%velocity%kinbcmask)  ! = 1 for Dirichlet BCs

!WHL - debug
         if (this_rank==model%numerics%rdiag_local) then
            print*, ' '
!!            print*, 'kinbcmask after halo update:'
            do j = model%general%nsn-1, 1, -1
               do i = 1, model%general%ewn-1
!!                  write(6,'(46i3)',advance='no') model%velocity%kinbcmask(i,j)
               enddo
!!               write(6,*) ' '
            enddo
         endif

         call glissade_velo_higher_solve(model,                                             &
                                         model%general%ewn,      model%general%nsn,         &
                                         model%general%upn)

!WHL - debug
         if (this_rank==model%numerics%rdiag_local) then
!            print*, ' '
!            print*, 'efvs (Pa yr, k = 1):'
            do j = model%general%nsn, 1, -1
               do i = 1, model%general%ewn
!                  write(6,'(8e10.2)',advance='no') model%stress%efvs(1,i,j) * evs0/scyr
               enddo
!               write(6,*) ' '
            enddo
         endif

      else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK

         !TODO - Create a JFNK solver that can work with an arbitrary calcF routine
         !       (i.e., glissade as well as glam)
         ! noxsolve could eventually go here 

         call write_log('JFNK not supported for glissade velocity solver',GM_FATAL)

      else   
         call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
      end if

    end subroutine glissade_velo_driver

end module glissade_velo
