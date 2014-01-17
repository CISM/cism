!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_basal_traction.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!!#include "glide_mask.inc"
#include "config.inc"

  module glissade_basal_traction

  !-----------------------------------------------------------------------------
  ! Compute or prescribe the basal traction coefficient 'beta' as required by
  ! the higher-order velocity solver.
  ! 
  ! Note that beta is assumed to be a positive constant.  In earlier versions of
  ! the code it was called 'betasquared'.
  !
  ! The units are Pa/(m/yr) if we assume a linear sliding law of the form
  !    taub_x = -beta * u, taub_y = -beta * v
  !
  ! However, the units are Pa if beta is treated as a till yield stress.
  !
  ! Current options are as follows:
  ! 
  ! [0] constant value of 10 Pa/(m/yr) (useful for debugging)
  ! [1] simple hard-coded pattern (useful for debugging)
  ! [2] treat beta value as a till yield stress (in Pa) using Picard iteration 
  ! [3] linear (inverse) function of basal water depth (bwat) 
  ! [4] very large value for beta to enforce no slip everywhere 
  ! [5] beta field passed in from .nc input file as part of standard i/o
  ! [6] no slip everywhere (using Dirichlet BC rather than large beta)
  ! [7] treat beta value as till yield stress (in Pa) using Newton-type iteration (in devel.)
  !
  ! TODO - Renumber so that, for example, the no-slip options have small numbers. 
  !-----------------------------------------------------------------------------

  use glimmer_paramets, only : dp
  use glimmer_physcon,  only : scyr
  use glimmer_paramets, only : vel0, tau0
  use glimmer_log,      only : write_log
  use glide_types
    
  !use glide_mask

  implicit none

!***********************************************************************

contains

!***********************************************************************

  subroutine calcbeta (whichbabc,                  &
                       dew,         dns,           &
                       ewn,         nsn,           &
                       thisvel,     othervel,      &
                       bwat,        beta_const,    &
                       mask,        beta,          &
                       betafile)

  ! subroutine to calculate map of beta sliding parameter, based on 
  ! user input ("whichbabc" flag, from config file as "which_ho_babc").

  use glimmer_paramets, only: len0
  use parallel, only: nhalo

  implicit none

  ! Input/output arguments

  integer, intent(in) :: whichbabc
  integer, intent(in) :: ewn, nsn

  real(dp), intent(in) :: dew, dns
  real(dp), intent(in), dimension(:,:) :: thisvel, othervel

  real(dp), intent(in), dimension(:,:) :: bwat  ! basal water depth

  real(dp), intent(in) :: beta_const  ! spatially uniform beta (Pa yr/m)

  integer, intent(in), dimension(:,:) :: mask 

  real(dp), intent(inout), dimension(ewn-1,nsn-1) :: beta

  character (len=30), intent(in), optional :: betafile

  ! Local variables

  real(dp) :: smallnum = 1.0d-2

  integer :: ew, ns

  ! SFP added for making beta a function of basal water flux 
  real(dp), dimension(:,:), allocatable :: unstagbeta
  real(dp) :: C, m

  real(dp) :: Ldomain   ! size of full domain
  real(dp) :: omega     ! frequency of beta field
  real(dp) :: dx, dy
  integer :: ilo, ihi, jlo, jhi  ! limits of beta field for ISHOM C case
  integer :: i, j

  ! Note that the dimensional scale (tau0 / vel0 / scyr ) is used here for making the basal traction coeff.
  ! beta dimensional, within the subroutine (mainly for debugging purposes), and then non-dimensional 
  ! again before being sent back out for use in the code. This scale is the same as "scale_beta" defined in 
  ! libglimmer/glimmer_scales.F90. See additional notes where that scale is defined. In general, below, it is
  ! assumed that values for beta being accessed from inside the code are already dimensionless and 
  ! any hardwired values have units of Pa yr/m. 

!TODO - Remove scaling here?


  select case(whichbabc)

    case(HO_BABC_CONSTANT)  ! spatially uniform value; useful for debugging and test cases

!      beta(:,:) = 10.d0       ! This is the default value (Pa yr/m)
      beta(:,:) = beta_const   ! Pa yr/m

    case(HO_BABC_SIMPLE)    ! simple pattern; also useful for debugging and test cases
                            ! (here, a strip of weak bed surrounded by stronger bed to simulate an ice stream)

      beta(:,:) = 1.d4        ! Pa yr/m

!TODO - Should these 5's be hardwired?  This will give strange results in parallel.
!       Could fix by applying small value of beta on global domain and scattering to local.  
!TODO - Is 10.d1 correct?  (Change to 100.d0?)
      do ns=5, nsn-5
      do ew=1, ewn-1
        beta(ew,ns) = 10.d1      ! Pa yr/m
      end do
      end do

    case(HO_BABC_YIELD_PICARD)  ! take input value for till yield stress and force beta to be implemented such
                                ! that plastic-till sliding behavior is enforced (see additional notes in documentation).

      !!! NOTE: Eventually, this option will provide the till yield stress as calculate from the basal processes
      !!! submodel. Currently, to enable sliding over plastic till, simple specify the value of "beta" as 
      !!! if it were the till yield stress (in units of Pascals).

      beta(:,:) = ( beta(:,:) * ( tau0 / vel0 / scyr ) ) &     ! Pa yr/m
                         / dsqrt( (thisvel(:,:)*vel0*scyr)**2 + (othervel(:,:)*vel0*scyr)**2 + (smallnum)**2 )

    case(HO_BABC_BETA_BWAT)  ! set value of beta as proportional to value of bwat                                         

      C = 10.d0
      m = 1.d0

      allocate(unstagbeta(ewn,nsn))

      unstagbeta(:,:) = 200.d0       

      where ( bwat > 0.d0 .and. unstagbeta > 200.d0 )
          unstagbeta = C / ( bwat**m )
      endwhere

      ! average beta from unstag grid onto stag grid
      beta = 0.5d0 * ( unstagbeta(1:ewn-1,:) + unstagbeta(2:ewn,:) )
      beta = 0.5d0 * ( unstagbeta(:,1:nsn-1) + unstagbeta(:,2:nsn) )
   
      deallocate(unstagbeta) 

    case(HO_BABC_LARGE_BETA)      ! frozen (u=v=0) ice-bed interface

      beta(:,:) = 1.d10       ! Pa yr/m

    case(HO_BABC_ISHOMC)          ! prescribe according to ISMIP-HOM test C

       !Note: Ideally, beta would be read in from an external netCDF file.
       !      However, this is not possible given that the global velocity grid is smaller
       !       than the ice grid and hence not able to fit the full beta field.
       !      The following code sets beta on the full grid as prescribed by Pattyn et al. (2008).

       !TODO - Abort if domain is not square
       Ldomain = (ewn-2*nhalo) * dew*len0   ! size of full domain (must be square)
       omega = 2.d0*pi / Ldomain

       ilo = nhalo
       ihi = ewn-nhalo
       jlo = nhalo
       jhi = nsn-nhalo
       
       print*, 'ilo, ihi =', ilo, ihi
       print*, 'jlo, jhi =', jlo, jhi

       ! Prescribe beta as in Pattyn et al., The Cryosphere, 2008
       beta(:,:) = 0.d0
       do j = jlo, jhi
          do i = ilo, ihi
             dx = dew*len0 * (i-ilo)
             dy = dns*len0 * (j-jlo)
             beta(i,j) = 1000.d0 + 1000.d0 * sin(omega*dx) * sin(omega*dy)
          enddo
       enddo

       print*, ' '
       print*, 'In calcbeta: beta, Pa/(m/yr):'
       do j = nsn-1, 1, -1
          do i = 1, ewn-1
             write(6,'(f8.0)',advance='no') beta(i,j)
          enddo
          print*, ' '
       enddo


    case(HO_BABC_EXTERNAL_BETA)   ! use value passed in externally from CISM (NOTE not dimensional when passed in) 

      ! scale CISM input value to dimensional units of (Pa yr/m)

      beta(:,:) = beta(:,:) * ( tau0 / vel0 / scyr )

      ! this is a check for NaNs, which indicate, and are replaced by no slip
      !TODO: Not sure I follow the logic of this ... keep/omit? Added by the UMT crew at some point

      do ns=1, nsn-1
      do ew=1, ewn-1 
        if( beta(ew,ns) /= beta(ew,ns) )then
          beta(ew,ns) = 1.d10     ! Pa yr/m
        endif 
      end do
      end do

      ! check for areas where ice is floating (or grounding line?) and make sure beta in these regions is 0  

      !TODO: Ideally, these mask values should not be hardwired, but keeping it this way for now until
      ! we decide which mask values to keep/remove

      do ns=1, nsn-1
      do ew=1, ewn-1 
        !if( ( mask(ew,ns) >= 21 .and. mask(ew,ns) <= 23 ) .or. ( mask(ew,ns) >= 41 .and. mask(ew,ns) <= 57 ) &
        !! less agressive than apply beta = 0 at g.l., which will make some test cases fail (e.g. circ. shelf)
        !! because of lack of fully grounded area.
        if( ( mask(ew,ns) >= 41 .and. mask(ew,ns) <= 43 ) &     
             .or. mask(ew,ns) == 9 .or. mask(ew,ns) == 11 )then
           beta(ew,ns) = 0.d0
        endif
      end do
      end do

    ! NOTE: cases (HO_BABC_NO_SLIP) and (HO_BABC_YIELD_NEWTON) are handled external to this subroutine

  end select

  ! convert the dimensional value of beta to non-dimensional units by dividing by scale factor.
  beta(:,:) = beta(:,:) / ( tau0 / vel0 / scyr )    !! scale in parentheses is: Pa * sec/m * yr/sec = Pa yr/m

end subroutine calcbeta

!***********************************************************************

end module glissade_basal_traction

!***********************************************************************
