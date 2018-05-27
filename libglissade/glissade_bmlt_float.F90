!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_bmlt_float.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

module glissade_bmlt_float

  use glimmer_global, only: dp
  use glimmer_physcon, only: rhoo, rhow, grav, lhci, scyr, pi
  use glimmer_log
  use glide_types
  use parallel

  implicit none
  
  private
  public :: glissade_basal_melting_float

  logical :: verbose_velo = .true.
  logical :: verbose_continuity = .true.
  logical :: verbose_melt = .true.

contains

!****************************************************

  subroutine glissade_basal_melting_float(whichbmlt_float,              &
                                          ewn,         nsn,             &
                                          dew,         dns,             &
                                          itest,       jtest,    rtest, &
                                          x1,                           &
                                          thck,        lsrf,            &
                                          topg,        eus,             &
                                          basal_melt)

    use glissade_masks, only: glissade_get_masks
    use glimmer_paramets, only: tim0, thk0

    ! Compute the rate of basal melting for floating ice by one of several methods.

    !-----------------------------------------------------------------
    ! Input/output arguments
    !-----------------------------------------------------------------

    integer, intent(in) :: whichbmlt_float            ! method for computing melt rate of floating ice

    integer, intent(in) ::  &
         ewn, nsn,             & ! grid dimensions
         itest, jtest, rtest     ! coordinates of diagnostic point

    real(dp), intent(in) :: &
         dew, dns                ! grid spacing in x, y (m)

    real(dp), dimension(:), intent(in) :: &
         x1                      ! x1 grid coordinates (m), ice grid
                                 ! used with bmlt_float_xlim for MISMIP+ Ice2r

    real(dp), dimension(:,:), intent(in) :: &
         lsrf,                 & ! elevation of lower ice surface (m)
         thck                    ! ice thickness (m)

    real(dp), dimension(:,:), intent(in) :: &
         topg                    ! elevation of bed topography (m)

    real(dp), intent(in) :: &
         eus                     ! eustatic sea level (m), = 0. by default

    type(glide_basal_melt), intent(inout) :: &
         basal_melt              ! derived type with fields and parameters related to basal melting

    !-----------------------------------------------------------------
    ! Note: The basal_melt derived type includes the 2D output field bmlt_float,
    !        along with a number of prescribed parameters for MISMIP+:
    !
    !       MISMIP+ Ice1
    !       - bmlt_float_omega     ! time scale for basal melting (s-1), default = 0.2/yr
    !       - bmlt_float_h0        ! scale for sub-shelf cavity thickness (m), default = 75 m
    !       - bmlt_float_z0        ! scale for ice draft (m), default = -100 m
    !
    !       MISMIP+ Ice2
    !       - bmlt_float_const     ! constant melt rate (m/s), default = 100 m/yr
    !       - bmlt_float_xlim      ! melt rate = 0 for abs(x) < bmlt_float_xlim (m), default = 480000 m
    !
    ! Note: Basal melt rates are > 0 for melting, < 0 for freeze-on
    !-----------------------------------------------------------------
 
    !----------------------------------------------------------------
    ! Local variables and pointers set to components of basal_melt derived type
    !----------------------------------------------------------------      

    real(dp), dimension(:,:), pointer :: &
         bmlt_float             ! basal melt rate for floating ice (m/s) (> 0 for melt, < 0 for freeze-on)

    !-----------------------------------------------------------------
    ! Local variables
    !-----------------------------------------------------------------

    integer, dimension(ewn,nsn) ::  &
         ice_mask,            & ! = 1 where ice temperature is computed (thck > thklim), else = 0
         floating_mask,       & ! = 1 where ice is present and floating, else = 0
         ocean_mask             ! = 1 where topg is below sea level and ice is absent

    integer :: i, j
    real(dp) :: h_cavity        ! depth of ice cavity beneath floating ice (m)
    real(dp) :: z_draft         ! elevation of lower surface of floating ice (m)

    real(dp) :: frz_ramp_factor    ! multiplying factor for linear ramp at depths with basal freezing
    real(dp) :: melt_ramp_factor   ! multiplying factor for linear ramp at depths with basal melting

    logical, parameter :: verbose_bmlt = .false.

    !-----------------------------------------------------------------
    ! Compute the basal melt rate for floating ice
    !-----------------------------------------------------------------

    if (main_task .and. verbose_bmlt) print*, 'Computing bmlt_float, whichbmlt_float =', whichbmlt_float

    ! Set bmlt_float pointer and initialize
    bmlt_float  => basal_melt%bmlt_float
    bmlt_float(:,:) = 0.0d0

    ! Compute masks:
    ! - ice_mask = 1 where thck > 0
    ! - floating_mask = 1 where thck > 0 and ice is floating;
    ! - ocean_mask = 1 where topg is below sea level and ice is absent
    !Note: The '0.0d0' argument is thklim. Here, any ice with thck > 0 gets ice_mask = 1.

    call glissade_get_masks(ewn,           nsn,            &
                            thck,          topg,           &
                            eus,           0.0d0,          &
                            ice_mask,                      &
                            floating_mask = floating_mask, &
                            ocean_mask = ocean_mask)

    if (whichbmlt_float == BMLT_FLOAT_CONSTANT) then

       ! Set melt rate to a constant value for floating ice.
       ! This includes ice-free ocean cells, in case ice is advected to those cells by the transport scheme.

       do j = 1, nsn
          do i = 1, ewn

             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                ! Note: For MISMIP+ experiment Ice2r, melting is masked out where x < 480 km

                if (abs(x1(i)) >= basal_melt%bmlt_float_xlim) then   ! melting is allowed
                   bmlt_float(i,j) = basal_melt%bmlt_float_const
                endif

                !WHL - debug
                if (j==jtest .and. this_rank==rtest) then
!!                if (i==itest .and. j==jtest .and. this_rank==rtest) then
!!                   print*, 'rank, i, j, bmlt_float:', this_rank, i, j, bmlt_float(i,j)
                endif
                   
             endif   ! ice is present and floating

          enddo
       enddo

    elseif (whichbmlt_float == BMLT_FLOAT_MISMIP) then   ! MISMIP+

       ! compute melt rate based on bed depth and cavity thickness
       !
       ! The MISMIP+ formula is as follows:
       !
       ! bmlt_float = omega * tanh(H_c/H_0) * max(z_0 - z_d, 0)
       !
       ! where H_c = lsrf - topg is the cavity thickness
       !       z_d = lsrf - eus is the ice draft
       !       omega = a time scale = 0.2 yr^{-1} by default
       !       H_0 = 75 m by default
       !       z_0 = 100 m by default

       do j = 1, nsn
          do i = 1, ewn

             ! compute basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme

             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                h_cavity = lsrf(i,j) - topg(i,j)
                z_draft = lsrf(i,j) - eus
                bmlt_float(i,j) = basal_melt%bmlt_float_omega * tanh(h_cavity/basal_melt%bmlt_float_h0) &
                                * max(basal_melt%bmlt_float_z0 - z_draft, 0.0d0)

                !debug
!                if (j == jtest .and. verbose_bmlt) then
!                   print*, 'cavity, tanh, thck, draft, melt rate (m/yr):', i, j, h_cavity, &
!                         tanh(h_cavity/basal_melt%bmlt_float_h0), thck(i,j), z_draft, bmlt_float(i,j)*scyr
!                endif

             endif   ! ice is present and floating

          enddo
       enddo

       !WHL - debug
       if (verbose_bmlt .and. this_rank == rtest) then
          print*, 'itest, jtest, rtest =', itest, jtest, rtest
          print*, ' '
          print*, 'thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'bmlt_float (m/yr), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') bmlt_float(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
       endif

    elseif (whichbmlt_float == BMLT_FLOAT_DEPTH) then

       ! Compute melt rates as a piecewise linear function of depth, generally with greater melting at depth.
       ! This scheme is similar to the MISMIP scheme, with the additional option of near-surface freezing.
       ! The maximum melting and freezing rates are set independently, with melting usually of greater magnitude.
       ! The melting/freezing rates fall linearly from their max values to zero over ranges defined by
       !  zmeltmax, zmelt0 and zfrzmax.
       ! The melt rate is set to a maximum value where z_draft <= zmeltmax,
       !  then decreases linearly to 0 as z_draft increases from zmeltmax to zmelt0.
       ! The freezing rate is set to a maximum value where z_draft >= zfrzmax,
       !  then decreases linearly to 0 as z_draft decreases from zfrzmax to zmelt0.
       ! (Here, z_draft < 0 by definition.)

       ! Compute ramp factors
       ! These factors are set to avoid divzero whe zmelt0 = zmeltmax, or zmelt0 = zfrzmax

       if (basal_melt%bmlt_float_depth_zfrzmax > basal_melt%bmlt_float_depth_zmelt0) then
          frz_ramp_factor = 1.0d0 / (basal_melt%bmlt_float_depth_zfrzmax - basal_melt%bmlt_float_depth_zmelt0)
       else
          frz_ramp_factor = 0.0d0
       endif

       if (basal_melt%bmlt_float_depth_zmelt0 > basal_melt%bmlt_float_depth_zmeltmax) then
          melt_ramp_factor = 1.0d0 / (basal_melt%bmlt_float_depth_zmelt0 - basal_melt%bmlt_float_depth_zmeltmax)
       else
          melt_ramp_factor = 0.0d0
       endif

       do j = 1, nsn
          do i = 1, ewn

             ! compute basal melt in ice-free ocean cells, in case ice is advected to those cells by the transport scheme
             if (floating_mask(i,j) == 1 .or. ocean_mask(i,j) == 1) then   ! ice is present and floating, or ice-free ocean

                z_draft = lsrf(i,j) - eus

                if (z_draft > basal_melt%bmlt_float_depth_zfrzmax) then
                   ! max freezing
                   bmlt_float(i,j) = -basal_melt%bmlt_float_depth_frzmax   ! frzmax >=0 by definition
                elseif (z_draft > basal_melt%bmlt_float_depth_zmelt0) then
                   ! freezing with a linear taper from frzmax to zero
                   bmlt_float(i,j) = -basal_melt%bmlt_float_depth_frzmax * &
                        frz_ramp_factor * (z_draft - basal_melt%bmlt_float_depth_zmelt0)
                elseif (z_draft > basal_melt%bmlt_float_depth_zmeltmax) then
                   ! melting with a linear taper from meltmax to zero
                   bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmax * &
                        melt_ramp_factor * (basal_melt%bmlt_float_depth_zmelt0 - z_draft)
                elseif (z_draft <= basal_melt%bmlt_float_depth_meltmax) then
                   ! max melting
                   bmlt_float(i,j) = basal_melt%bmlt_float_depth_meltmax
                endif

                ! As with the MISMIP scheme, reduce the melting as the cavity thickness approaches zero.
                ! A small value of bmlt_float_h0 allows more melting in very thin cavities.
                h_cavity = max(lsrf(i,j) - topg(i,j), 0.0d0)
                bmlt_float(i,j) = bmlt_float(i,j) * tanh(h_cavity/basal_melt%bmlt_float_h0)

             endif   ! ice is present and floating

          enddo
       enddo

       !debug
       if (verbose_bmlt .and. this_rank == rtest) then
          print*, 'itest, jtest, rtest =', itest, jtest, rtest
          print*, ' '
          print*, 'topg (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') topg(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'thck (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') thck(i,j)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'z_draft (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') min(lsrf(i,j) - eus, 0.0d0)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'h_cavity (m):'
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') max(lsrf(i,j) - topg(i,j), 0.0d0)
             enddo
             write(6,*) ' '
          enddo
          print*, ' '
          print*, 'bmlt_float (m/yr), rank =', rtest
          do j = jtest+3, jtest-3, -1
             do i = itest-3, itest+3
                write(6,'(f10.2)',advance='no') bmlt_float(i,j)*scyr
             enddo
             write(6,*) ' '
          enddo
       endif

    endif   ! whichbmlt_float

  end subroutine glissade_basal_melting_float

!****************************************************

end module glissade_bmlt_float

!****************************************************
