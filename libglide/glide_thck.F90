! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_thck.f90 - part of the Glimmer-CISM ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glide_thck

  use glide_types

  private
  public :: init_thck, thck_nonlin_evolve, thck_lin_evolve, stagvarb, timeders, stagleapthck

#ifdef DEBUG_PICARD
  ! debugging Picard iteration
  integer, private, parameter :: picard_unit=101
  real, private, parameter    :: picard_interval=500.
  integer, private            :: picard_max=0
#endif

contains

  subroutine init_thck(model)
    !*FD initialise work data for ice thickness evolution
    use glimmer_log
    implicit none
    type(glide_global_type) :: model

    
    model%pcgdwk%fc2 = (/ model%numerics%alpha * model%numerics%dt / (2.0d0 * model%numerics%dew * model%numerics%dew), &
         model%numerics%dt, (1.0d0-model%numerics%alpha) / model%numerics%alpha, &
         1.0d0 / model%numerics%alpha, model%numerics%alpha * model%numerics%dt / &
         (2.0d0 * model%numerics%dns * model%numerics%dns), 0.0d0 /) 

#ifdef DEBUG_PICARD
    call write_log('Logging Picard iterations')
    open(picard_unit,name='picard_info.data',status='unknown')
    write(picard_unit,*) '#time    max_iter'
#endif

    ! allocate memory for ADI scheme
    if (model%options%whichevol.eq.1) then
       allocate(model%thckwk%alpha(max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%beta (max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%gamma(max(model%general%ewn, model%general%nsn)))
       allocate(model%thckwk%delta(max(model%general%ewn, model%general%nsn)))
    end if
  end subroutine init_thck

!---------------------------------------------------------------------------------

  subroutine thck_lin_evolve(model,newtemps)

    !*FD this subroutine solves the linearised ice thickness equation by computing the
    !*FD diffusivity from quantities of the previous time step

    use glide_velo
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if
       call slipvelo(model,                &
            2,                             &
            model%velocity% btrc,          &
            model%velocity% ubas,          &
            model%velocity% vbas)

       ! calculate diffusivity
       call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%velocity%diffu)

       ! get new thicknesses
       call thck_evolve(model,.true.,model%geometry%thck,model%geometry%thck)

       ! calculate horizontal velocity field
       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
            model%velocity%ubas,           &
            model%velocity%vbas)
       call velo_calc_velo(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if
  end subroutine thck_lin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_nonlin_evolve(model,newtemps)

    !*FD this subroutine solves the ice thickness equation by doing an outer, 
    !*FD non-linear iteration to update the diffusivities and in inner, linear
    !*FD iteration to calculate the new ice thickness distrib

    use glimmer_global, only : dp
    use glide_velo
    use glide_setup
    use glimmer_deriv, only : df_field_2d_staggered 
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    ! local variables
    integer, parameter :: pmax=50                       !*FD maximum Picard iterations
    real(kind=dp), parameter :: tol=1.0d-6
    real(kind=dp) :: residual
    integer p
    logical first_p


    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
#ifdef DEBUG
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if

       first_p = .true.
       model%thckwk%oldthck = model%geometry%thck
       ! do Picard iteration
       do p=1,pmax
          model%thckwk%oldthck2 = model%geometry%thck

          call stagvarb(model%geometry% thck, &
               model%geomderv% stagthck,&
               model%general%  ewn, &
               model%general%  nsn)

          call df_field_2d_staggered(model%geometry%usrf, &
               model%numerics%dew, model%numerics%dns, &
               model%geomderv%dusrfdew, & 
               model%geomderv%dusrfdns, &
               .false., .false.)

          call df_field_2d_staggered(model%geometry%thck, &
               model%numerics%dew, model%numerics%dns, &
               model%geomderv%dthckdew, & 
               model%geomderv%dthckdns, &
               .false., .false.)

          call slipvelo(model,                &
               2,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)

          ! calculate diffusivity
          call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
               model%geomderv%dusrfdns,model%velocity%diffu)

          ! get new thicknesses
          call thck_evolve(model,first_p,model%thckwk%oldthck,model%geometry%thck)

          first_p = .false.
          residual = maxval(abs(model%geometry%thck-model%thckwk%oldthck2))
          if (residual.le.tol) then
             exit
          end if
          
       end do
#ifdef DEBUG_PICARD
       picard_max=max(picard_max,p)
       if (model%numerics%tinc > mod(model%numerics%time,picard_interval)) then
          write(picard_unit,*) model%numerics%time,p
          picard_max = 0
       end if
#endif

       ! calculate horizontal velocity field
       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
            model%velocity%ubas,           &
            model%velocity%vbas)
       call velo_calc_velo(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if
  end subroutine thck_nonlin_evolve

!---------------------------------------------------------------------------------

  subroutine thck_evolve(model,calc_rhs,old_thck,new_thck)

    !*FD set up sparse matrix and solve matrix equation to find new ice thickness distribution
    !*FD this routine does not override the old thickness distribution

    use glide_setup, only: glide_calclsrf
    use glimmer_global, only : dp

    implicit none

    ! subroutine arguments -------------------------------------------------------------

    type(glide_global_type) :: model
    logical,intent(in) :: calc_rhs                      !*FD set to true when rhs should be calculated 
                                                        !*FD i.e. when doing lin solution or first picard iteration
    real(dp), intent(in), dimension(:,:) :: old_thck    !*FD contains ice thicknesses from previous time step
    real(dp), intent(inout), dimension(:,:) :: new_thck !*FD on entry contains first guess for new ice thicknesses
                                                        !*FD on exit contains ice thicknesses of new time step

    ! local variables ------------------------------------------------------------------

    real(dp), dimension(5) :: sumd 
    real(dp) :: err
    integer :: linit
    integer :: ew,ns

    ! the number of grid points
    model%pcgdwk%pcgsize(1) = model%geometry%totpts

    ! Zero the arrays holding the sparse matrix
    model%pcgdwk%pcgval = 0.0
    model%pcgdwk%pcgcol = 0 
    model%pcgdwk%pcgrow = 0
    model%pcgdwk%ct = 1

    ! Boundary Conditions ---------------------------------------------------------------
    ! lower and upper BC
    do ew = 1,model%general%ewn
       ns=1
       if (model%geometry%mask(ew,ns) /= 0) then
          call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
          if (calc_rhs) then
             model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
          end if
          model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
       end if
       ns=model%general%nsn
       if (model%geometry%mask(ew,ns) /= 0) then
          call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
          if (calc_rhs) then
             model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
          end if
          model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
       end if
    end do

    !left and right BC
    if (model%options%periodic_ew.eq.1) then
       do ns=2,model%general%nsn-1
          ew = 1
          if (model%geometry%mask(ew,ns) /= 0) then
             call findsums(model%general%ewn-2,model%general%ewn-1,ns-1,ns)
             call generate_row(model%general%ewn-2,ew,ew+1,ns-1,ns,ns+1)
          end if
          ew=model%general%ewn
          if (model%geometry%mask(ew,ns) /= 0) then
             call findsums(1,2,ns-1,ns)
             call generate_row(ew-1,ew,3,ns-1,ns,ns+1)
          end if
       end do
    else
       do ns=2,model%general%nsn-1
          ew=1
          if (model%geometry%mask(ew,ns) /= 0) then
             call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
             if (calc_rhs) then
                model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
             end if
             model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
          end if
          ew=model%general%ewn
          if (model%geometry%mask(ew,ns) /= 0) then
             call putpcgc(model%pcgdwk,1.0d0, model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))
             if (calc_rhs) then
                model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) = old_thck(ew,ns) 
             end if
             model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)
          end if
       end do
    end if

    ! ice body -------------------------------------------------------------------------

    do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1

          if (model%geometry%mask(ew,ns) /= 0) then
                
             call findsums(ew-1,ew,ns-1,ns)
             call generate_row(ew-1,ew,ew+1,ns-1,ns,ns+1)

          end if
       end do
    end do

    ! Calculate the total number of points
    model%pcgdwk%pcgsize(2) = model%pcgdwk%ct - 1 

    ! Solve the system using SLAP
    call slapsolv(model,linit,err)   

    ! Rejig the solution onto a 2D array
    do ns = 1,model%general%nsn
       do ew = 1,model%general%ewn 

          if (model%geometry%mask(ew,ns) /= 0) then
             new_thck(ew,ns) = model%pcgdwk%answ(model%geometry%mask(ew,ns))
          end if

       end do
    end do

    new_thck = max(0.0d0, new_thck)

#ifdef DEBUG
    print *, "* thck ", model%numerics%time, linit, model%geometry%totpts, &
         real(thk0*new_thck(model%general%ewn/2+1,model%general%nsn/2+1)), &
         real(vel0*maxval(abs(model%velocity%ubas))), real(vel0*maxval(abs(model%velocity%vbas))) 
#endif

    ! calculate upper and lower surface
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  contains

    subroutine generate_row(ewm,ew,ewp,nsm,ns,nsp)
      ! calculate row of sparse matrix equation
      implicit none
      integer, intent(in) :: ewm,ew,ewp  ! ew index to left, central, right node
      integer, intent(in) :: nsm,ns,nsp  ! ns index to lower, central, upper node

      ! fill sparse matrix
      call putpcgc(model%pcgdwk,sumd(1), model%geometry%mask(ewm,ns), model%geometry%mask(ew,ns))       ! point (ew-1,ns)
      call putpcgc(model%pcgdwk,sumd(2), model%geometry%mask(ewp,ns), model%geometry%mask(ew,ns))       ! point (ew+1,ns)
      call putpcgc(model%pcgdwk,sumd(3), model%geometry%mask(ew,nsm), model%geometry%mask(ew,ns))       ! point (ew,ns-1)
      call putpcgc(model%pcgdwk,sumd(4), model%geometry%mask(ew,nsp), model%geometry%mask(ew,ns))       ! point (ew,ns+1)
      call putpcgc(model%pcgdwk,1.0d0 + sumd(5), model%geometry%mask(ew,ns), model%geometry%mask(ew,ns))! point (ew,ns)

      ! calculate RHS
      if (calc_rhs) then
         model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) =                    &
              old_thck(ew,ns) * (1.0d0 - model%pcgdwk%fc2(3) * sumd(5))     &
            - model%pcgdwk%fc2(3) * (old_thck(ewm,ns) * sumd(1)             &
                                   + old_thck(ewp,ns) * sumd(2)             &
                                   + old_thck(ew,nsm) * sumd(3)             &
                                   + old_thck(ew,nsp) * sumd(4))            &
            - model%pcgdwk%fc2(4) * (model%geometry%lsrf(ew,ns)  * sumd(5)  &
                                   + model%geometry%lsrf(ewm,ns) * sumd(1)  &
                                   + model%geometry%lsrf(ewp,ns) * sumd(2)  &
                                   + model%geometry%lsrf(ew,nsm) * sumd(3)  &
                                   + model%geometry%lsrf(ew,nsp) * sumd(4)) &
            + model%climate%acab(ew,ns) * model%pcgdwk%fc2(2)
         if(model%options%basal_mbal==1) then
            model%pcgdwk%rhsd(model%geometry%mask(ew,ns)) =                    &
                 model%pcgdwk%rhsd(model%geometry%mask(ew,ns))                 &
                 - model%temper%bmlt(ew,ns) * model%pcgdwk%fc2(2) ! basal melt is +ve for mass loss
         end if
      end if

      model%pcgdwk%answ(model%geometry%mask(ew,ns)) = new_thck(ew,ns)      

    end subroutine generate_row

    subroutine findsums(ewm,ew,nsm,ns)
      ! calculate diffusivities
      implicit none
      integer, intent(in) :: ewm,ew  ! ew index to left, right
      integer, intent(in) :: nsm,ns  ! ns index to lower, upper

      ! calculate sparse matrix elements
      sumd(1) = model%pcgdwk%fc2(1) * (&
           (model%velocity%diffu(ewm,nsm) + model%velocity%diffu(ewm,ns)) + &
           (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ewm,ns)))
      sumd(2) = model%pcgdwk%fc2(1) * (&
           (model%velocity%diffu(ew,nsm) + model%velocity%diffu(ew,ns)) + &
           (model%velocity%ubas (ew,nsm) + model%velocity%ubas (ew,ns)))
      sumd(3) = model%pcgdwk%fc2(5) * (&
           (model%velocity%diffu(ewm,nsm) + model%velocity%diffu(ew,nsm)) + &
           (model%velocity%ubas (ewm,nsm) + model%velocity%ubas (ew,nsm)))
      sumd(4) = model%pcgdwk%fc2(5) * (&
           (model%velocity%diffu(ewm,ns) + model%velocity%diffu(ew,ns)) + &
           (model%velocity%ubas (ewm,ns) + model%velocity%ubas (ew,ns)))
      sumd(5) = - (sumd(1) + sumd(2) + sumd(3) + sumd(4))
    end subroutine findsums
  end subroutine thck_evolve

!---------------------------------------------------------------------------------

  subroutine slapsolv(model,iter,err)

    use glimmer_global, only : dp 
    use glide_stop
    use glimmer_log
    use glimmer_filenames

    implicit none

    type(glide_global_type) :: model  !*FD Model data type
    integer,  intent(out)   :: iter   !*FD Number of iterations
    real(dp), intent(out)   :: err    !*FD Error estimate of result

    ! For call to dslucs
    real(dp), parameter :: tol   = 1.0d-12
    integer,  parameter :: isym  = 0
    integer,  parameter :: itol  = 2
    integer,  parameter :: itmax = 101
    real(dp), dimension(:),allocatable :: rwork ! Real work array
    integer,  dimension(:),allocatable :: iwork ! Integer work array
    integer :: ierr, mxnelt

    ! Variables for error handling
    character(200) :: message
    character(100) :: errfname
    integer :: lunit

    mxnelt = 20 * model%pcgdwk%pcgsize(1)
    allocate(rwork(mxnelt),iwork(mxnelt))

    ! solve the problem using the SLAP package routines     

    call dslucs(model%pcgdwk%pcgsize(1), &  ! n  ... order of matrix a (in)
                model%pcgdwk%rhsd,       &  ! b  ... right hand side vector (in)
                model%pcgdwk%answ,       &  ! x  ... initial quess/final solution vector (in/out)
                model%pcgdwk%pcgsize(2), &  ! nelt ... number of non-zeroes in A (in)
                model%pcgdwk%pcgrow,     &  ! ia  ... sparse matrix format of A (in)
                model%pcgdwk%pcgcol,     &  ! ja  ... sparse matrix format of A (in)
                model%pcgdwk%pcgval,     &  ! a   ... matrix (in)
                isym,                    &  ! isym ... storage method (0 is complete) (in)
                itol,                    &  ! itol ... convergence criteria (2 recommended) (in)
                tol,                     &  ! tol  ... criteria for convergence (in)
                itmax,                   &  ! itmax ... maximum number of iterations (in)
                iter,                    &  ! iter  ... returned number of iterations (out)
                err,                     &  ! err   ... error estimate of solution (out)
                ierr,                    &  ! ierr  ... returned error message (0 is ok) (out)
                0,                       &  ! iunit ... unit for error writes during iteration (0 no write) (in)
                rwork,                   &  ! rwork ... workspace for SLAP routines (in)
                mxnelt,                  &  ! lenw
                iwork,                   &  ! iwork ... workspace for SLAP routines (in)
                mxnelt)                     ! leniw

    ! Handle errors gracefully
    if (ierr /= 0) then

       ! Acquire a file unit, and open the file
       lunit = get_free_unit()
       errfname = trim(process_path('slap_dump.txt'))
       open(lunit,file=errfname,status='unknown')

       ! Output data to file
       call dcpplt(model%pcgdwk%pcgsize(1), &
                   model%pcgdwk%pcgsize(2), &
                   model%pcgdwk%pcgrow,     &
                   model%pcgdwk%pcgcol,     &
                   model%pcgdwk%pcgval,     &
                   isym,                    &
                   lunit)
       write(lunit,*) '***SLAP data ends. PCGVAL follows'
       write(lunit,*) model%pcgdwk%pcgval

       ! Close unit and finish off
       close(lunit)
       call glide_finalise(model,.true.)
       write(message,*)'SLAP solution error at time: ',model%numerics%time, &
            '. Data dumped to ',trim(errfname)
       call write_log(trim(message),GM_FATAL,__FILE__,__LINE__)
    end if

    deallocate(rwork,iwork)

  end subroutine slapsolv 

!---------------------------------------------------------------------------------

  subroutine putpcgc(pcgdwk,value,col,row)

    use glimmer_global, only : dp

    implicit none

    type(glide_pcgdwk) :: pcgdwk
    integer, intent(in) :: row, col
    real(dp), intent(in) :: value

    if (value /= 0.0d0) then
      pcgdwk%pcgval(pcgdwk%ct) = value
      pcgdwk%pcgcol(pcgdwk%ct) = col
      pcgdwk%pcgrow(pcgdwk%ct) = row
      pcgdwk%ct = pcgdwk%ct + 1
    end if

  end subroutine putpcgc

!---------------------------------------------------------------------------------

  subroutine stagvarb(ipvr,opvr,ewn,nsn)

    use glimmer_global, only : dp ! ewn, nsn
 
    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    integer :: ewn,nsn

    opvr(1:ewn-1,1:nsn-1) = (ipvr(2:ewn,1:nsn-1) + ipvr(1:ewn-1,2:nsn) + &
                             ipvr(2:ewn,2:nsn)   + ipvr(1:ewn-1,1:nsn-1)) / 4.0d0

  end subroutine stagvarb

!---------------------------------------------------------------------------------

  subroutine timeders(thckwk,ipvr,opvr,mask,time,which)

    !*FD Calculates the time-derivative of a field. This subroutine is used by 
    !*FD the temperature solver only.

    use glimmer_global, only : dp, sp
    use glimmer_paramets, only : conv

    implicit none 

    type(glide_thckwk) :: thckwk    !*FD Derived-type containing work data
    real(dp), intent(out), dimension(:,:) :: opvr  !*FD Input field
    real(dp), intent(in),  dimension(:,:) :: ipvr  !*FD Output (derivative) field
    real(sp), intent(in)                  :: time  !*FD current time
    integer,  intent(in),  dimension(:,:) :: mask  !*FD mask for calculation
    integer,  intent(in)                  :: which !*FD selector for stored field

    real(sp) :: factor

    factor = (time - thckwk%oldtime)
    if (factor .eq.0) then
       opvr = 0.0d0
    else
       factor = 1./factor
       where (mask /= 0)
          opvr = conv * (ipvr - thckwk%olds(:,:,which)) * factor
       elsewhere
          opvr = 0.0d0
       end where
    end if

    thckwk%olds(:,:,which) = ipvr

    if (which == thckwk%nwhich) then
      thckwk%oldtime = time
    end if

  end subroutine timeders

!---------------------------------------------------------------------------------

  subroutine filterthck(thck,ewn,nsn)

    use glimmer_global, only : dp ! ew, ewn, ns, nsn

    implicit none

    real(dp), dimension(:,:), intent(inout) :: thck
    real(dp), dimension(:,:), allocatable :: smth
    integer :: ewn,nsn

    real(dp), parameter :: f = 0.1d0 / 16.0d0
    integer :: count
    integer :: ns,ew

    allocate(smth(ewn,nsn))
    count = 1

    do ns = 3,nsn-2
      do ew = 3,ewn-2

        if (all((thck(ew-2:ew+2,ns) > 0.0d0)) .and. all((thck(ew,ns-2:ns+2) > 0.0d0))) then
          smth(ew,ns) =  thck(ew,ns) + f * &
                        (thck(ew-2,ns) - 4.0d0 * thck(ew-1,ns) + 12.0d0 * thck(ew,ns) - &
                         4.0d0 * thck(ew+1,ns) + thck(ew+2,ns) + &
                         thck(ew,ns-2) - 4.0d0 * thck(ew,ns-1) - &
                         4.0d0 * thck(ew,ns+1) + thck(ew,ns+2))
          count = count + 1
        else
          smth(ew,ns) = thck(ew,ns)
        end if

      end do
    end do

    thck(3:ewn-2,3:nsn-2) = smth(3:ewn-2,3:nsn-2)
    print *, count

    deallocate(smth)            

  end subroutine filterthck

!----------------------------------------------------------------------

  subroutine swapbndh(bc,a,b,c,d)

    use glimmer_global, only : dp

    implicit none

    real(dp), intent(out), dimension(:) :: a, c
    real(dp), intent(in), dimension(:) :: b, d
    integer, intent(in) :: bc

    if (bc == 0) then
      a = b
      c = d
    end if

  end subroutine swapbndh

  !-----------------------------------------------------------------------------
  ! ADI routines
  !-----------------------------------------------------------------------------

  subroutine stagleapthck(model,newtemps)
    
    !*FD this subroutine solves the ice sheet thickness equation using the ADI scheme
    !*FD diffusivities are updated for each half time step

    use glide_setup, only: glide_calclsrf
    use glide_velo
    use glimmer_utils
    implicit none
    ! subroutine arguments
    type(glide_global_type) :: model
    logical, intent(in) :: newtemps                     !*FD true when we should recalculate Glen's A

    ! local variables
    integer ew,ns, n

    if (model%geometry%empty) then

       model%geometry%thck = dmax1(0.0d0,model%geometry%thck + model%climate%acab * model%pcgdwk%fc2(2))
#ifdef DEBUG       
       print *, "* thck empty - net accumulation added", model%numerics%time
#endif
    else

       ! calculate basal velos
       if (newtemps) then
          call slipvelo(model,                &
               1,                             &
               model%velocity% btrc,          &
               model%velocity% ubas,          &
               model%velocity% vbas)
          ! calculate Glen's A if necessary
          call velo_integrate_flwa(model%velowk,model%geomderv%stagthck,model%temper%flwa)
       end if
       call slipvelo(model,                &
            2,                             &
            model%velocity% btrc,          &
            model%velocity% ubas,          &
            model%velocity% vbas)

       ! calculate diffusivity
       call velo_calc_diffu(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%velocity%diffu)

       model%velocity%total_diffu(:,:) = model%velocity%diffu(:,:) + model%velocity%ubas(:,:)

       ! first ADI step, solve thickness equation along rows j
       n = model%general%ewn
       do ns=2,model%general%nsn-1
          call adi_tri ( model%thckwk%alpha,                 &
                         model%thckwk%beta,                  &
                         model%thckwk%gamma,                 &
                         model%thckwk%delta,                 &
                         model%geometry%thck(:,ns),          &
                         model%geometry%lsrf(:,ns),          &
                         model%climate%acab(:,ns)-real(model%options%basal_mbal)*real(model%temper%bmlt(:,ns),sp),           &
                         model%velocity%vflx(:,ns),          &
                         model%velocity%vflx(:,ns-1),        &
                         model%velocity%total_diffu(:,ns),   &
                         model%velocity%total_diffu(:,ns-1), &
                         model%numerics%dt,                  &
                         model%numerics%dew,                 &
                         model%numerics%dns )

          call tridiag(model%thckwk%alpha(1:n),    &
                       model%thckwk%beta(1:n),     &
                       model%thckwk%gamma(1:n),    &
                       model%thckwk%oldthck(:,ns), &
                       model%thckwk%delta(1:n))
       end do

       model%thckwk%oldthck(:,:) = max(model%thckwk%oldthck(:,:), 0.d0)

       ! second ADI step, solve thickness equation along columns i
       n = model%general%nsn
       do ew=2,model%general%ewn-1
          call adi_tri ( model%thckwk%alpha,                 &
                         model%thckwk%beta,                  &
                         model%thckwk%gamma,                 &
                         model%thckwk%delta,                 &
                         model%thckwk%oldthck(ew,:),         &
                         model%geometry%lsrf(ew, :),         &
                         model%climate%acab(ew, :)-real(model%options%basal_mbal)*real(model%temper%bmlt(ew, :),sp),          &
                         model%velocity%uflx(ew,:),          &
                         model%velocity%uflx(ew-1,:),        &
                         model%velocity%total_diffu(ew,:),   &
                         model%velocity%total_diffu(ew-1,:), &
                         model%numerics%dt,                  &
                         model%numerics%dns,                 &
                         model%numerics%dew )

          call tridiag(model%thckwk%alpha(1:n),    &
                       model%thckwk%beta(1:n),     &
                       model%thckwk%gamma(1:n),    &
                       model%geometry%thck(ew, :), &
                       model%thckwk%delta(1:n))
       end do

       model%geometry%thck(:,:) = max(model%geometry%thck(:,:), 0.d0)

       ! Apply boundary conditions
       model%geometry%thck(1,:) = 0.0
       model%geometry%thck(model%general%ewn,:) = 0.0
       model%geometry%thck(:,1) = 0.0
       model%geometry%thck(:,model%general%nsn) = 0.0

       ! calculate horizontal velocity field
       call slipvelo(model,                &
            3,                             &
            model%velocity%btrc,           &
            model%velocity%ubas,           &
            model%velocity%vbas)
       call velo_calc_velo(model%velowk,model%geomderv%stagthck,model%geomderv%dusrfdew, &
            model%geomderv%dusrfdns,model%temper%flwa,model%velocity%diffu,model%velocity%ubas, &
            model%velocity%vbas,model%velocity%uvel,model%velocity%vvel,model%velocity%uflx,model%velocity%vflx)
    end if

    !------------------------------------------------------------
    ! calculate upper and lower surface
    !------------------------------------------------------------
    call glide_calclsrf(model%geometry%thck, model%geometry%topg, model%climate%eus, model%geometry%lsrf)
    model%geometry%usrf = max(0.d0,model%geometry%thck + model%geometry%lsrf)

  end subroutine stagleapthck

!---------------------------------------------------------------------------------

  subroutine adi_tri(a,b,c,d,thk,tpg,mb,flx_p,flx_m,dif_p,dif_m,dt,ds1, ds2)
    !*FD construct tri-diagonal matrix system for a column/row
    use glimmer_global, only : dp, sp
    implicit none
    
    real(dp), dimension(:), intent(out) :: a !*FD alpha (subdiagonal)
    real(dp), dimension(:), intent(out) :: b !*FD alpha (diagonal)
    real(dp), dimension(:), intent(out) :: c !*FD alpha (superdiagonal)
    real(dp), dimension(:), intent(out) :: d !*FD right-hand side
    
    real(dp), dimension(:), intent(in) :: thk   !*FD ice thickness
    real(dp), dimension(:), intent(in) :: tpg   !*FD lower surface of ice
    real(sp), dimension(:), intent(in) :: mb    !*FD mass balance
    real(dp), dimension(:), intent(in) :: flx_p !*FD flux +1/2
    real(dp), dimension(:), intent(in) :: flx_m !*FD flux -1/2
    real(dp), dimension(:), intent(in) :: dif_p !*FD diffusivity +1/2
    real(dp), dimension(:), intent(in) :: dif_m !*FD diffusivity -1/2
    
    real(dp), intent(in) :: dt !*FD time step
    real(dp), intent(in) :: ds1, ds2 !*FD spatial steps inline and transversal

    ! local variables
    real(dp) :: f1, f2, f3
    integer :: i,n
    
    n = size(thk)

    f1 = dt/(4*ds1*ds1)
    f2 = dt/(4*ds2)
    f3 = dt/2.

    a(:) = 0.
    b(:) = 0.
    c(:) = 0.
    d(:) = 0.

    a(1) = 0.
    do i=2,n
       a(i) = f1*(dif_m(i-1)+dif_p(i-1))
    end do
    do i=1,n-1
       c(i) = f1*(dif_m(i)+dif_p(i))
    end do
    c(n) = 0.
    b(:) = -(a(:)+c(:))

    ! calculate RHS
    do i=2,n-1
       d(i) = thk(i) - &
            f2 * (flx_p(i-1) + flx_p(i) - flx_m(i-1) - flx_m(i)) + &
            f3 * mb(i) - &
            a(i)*tpg(i-1) - b(i)*tpg(i) - c(i)*tpg(i+1)
    end do

    b(:) = 1.+b(:)

  end subroutine adi_tri

end module glide_thck

