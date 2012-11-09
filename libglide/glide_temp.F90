!CLEANUP - glide_temp.F90 
! Moved the contents of glide_temp_utils.F90 back to this module.
! Decided it was easier to support two different versions of these routines
!  rather than make one version work for both SIA and HO

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_temp.f90 - part of the GLIMMER ice model           +
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009, 2010
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

#include "glide_mask.inc"

! some macros used to disable parts of the temperature equation
! vertical diffusion
#ifdef NO_VERTICAL_DIFFUSION
#define VERT_DIFF 0.
#else
#define VERT_DIFF 1.
#endif

! horizontal advection
#ifdef NO_HORIZONTAL_ADVECTION
#define HORIZ_ADV 0.
#else
#define HORIZ_ADV 1.
#endif

! vertical advection
#ifdef NO_VERTICAL_ADVECTION
#define VERT_ADV 0.
#else
#define VERT_ADV 1.
#endif

! strain heating
#ifdef NO_STRAIN_HEAT
#define STRAIN_HEAT 0.
#else
#define STRAIN_HEAT 1.
#endif

module glide_temp

  use glide_types

  implicit none

  private
  public :: glide_init_temp, glide_temp_driver, glide_calcbmlt

contains

!------------------------------------------------------------------------------------

!TODO - Replace model derived type with explicit arguments?

  subroutine glide_init_temp(model)

    !*FD initialise temperature module
    use glimmer_physcon, only : rhoi, shci, coni, scyr, grav, gn, lhci, rhow
    use glimmer_paramets, only : tim0, thk0, acc0, len0, vis0, vel0
    use glimmer_global, only : dp 
    use glimmer_log
    use glide_bwater, only : find_dt_wat

    type(glide_global_type), intent(inout) :: model       ! model instance

    integer, parameter :: p1 = gn + 1  
    integer up, ns, ew
    real(dp) :: estimate

    if (VERT_DIFF==0.)   call write_log('Vertical diffusion is switched off')
    if (HORIZ_ADV==0.)   call write_log('Horizontal advection is switched off')
    if (VERT_ADV==0.)    call write_log('Vertical advection is switched off')
    if (STRAIN_HEAT==0.) call write_log('Strain heating is switched off')

    ! horizontal advection stuff
    allocate(model%tempwk%hadv_u(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%hadv_v(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%initadvt(model%general%upn,model%general%ewn,model%general%nsn))

    allocate(model%tempwk%inittemp(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%dissip(model%general%upn,model%general%ewn,model%general%nsn))
    allocate(model%tempwk%compheat(model%general%upn,model%general%ewn,model%general%nsn))
    model%tempwk%compheat = 0.0d0
    allocate(model%tempwk%dups(model%general%upn,3))

    allocate(model%tempwk%c1(model%general%upn))

    allocate(model%tempwk%dupa(model%general%upn),model%tempwk%dupb(model%general%upn))
    allocate(model%tempwk%dupc(model%general%upn))

    allocate(model%tempwk%smth(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%wphi(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatu(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bwatv(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%fluxew(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%fluxns(model%general%ewn,model%general%nsn))
    allocate(model%tempwk%bint(model%general%ewn-1,model%general%nsn-1))

    model%tempwk%advconst(1) = HORIZ_ADV*model%numerics%dttem / (16.0d0 * model%numerics%dew)
    model%tempwk%advconst(2) = HORIZ_ADV*model%numerics%dttem / (16.0d0 * model%numerics%dns)

    model%tempwk%dups = 0.0d0

    do up = 2, model%general%upn-1
       model%tempwk%dups(up,1) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) * &
            (model%numerics%sigma(up)   - model%numerics%sigma(up-1)))
       model%tempwk%dups(up,2) = 1.d0/((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) *  &
            (model%numerics%sigma(up+1) - model%numerics%sigma(up)))
       model%tempwk%dups(up,3) = 1.d0/(model%numerics%sigma(up+1)  - model%numerics%sigma(up-1))
    end do

    model%tempwk%zbed = 1.0d0 / thk0
    model%tempwk%dupn = model%numerics%sigma(model%general%upn) - model%numerics%sigma(model%general%upn-1)

! In dimensional units, wmax = thk0 / (tim0/scyr) = 2000 m / 400 yr = 5 m/yr
! In nondimensional units, wmax = 5 m/yr / (thk0*scyr/tim0) = 1.0
! If we remove scaling, then tim0 = thk0 = 1, and wmax = 5 m/yr / scyr.  
! The following expression is correct if scaling is removed.

    model%tempwk%wmax = 5.0d0 * tim0 / (scyr * thk0)

!TODO - Delete line with tau0
    model%tempwk%cons = (/ 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2), &
         model%numerics%dttem / 2.0d0, &
         VERT_DIFF*2.0d0 * tim0 * model%numerics%dttem / (thk0 * rhoi * shci), &
         VERT_ADV*tim0 * acc0 * model%numerics%dttem / coni, &
!!         ( tau0 * vel0 / len0 ) / ( rhoi * shci ) * ( model%numerics%dttem * tim0 ) /)  
         0.d0 /)   !whl - last term no longer needed
         !*sfp* added last term to vector above for use in HO & SSA dissip. cacl
         

    model%tempwk%c1 = STRAIN_HEAT *(model%numerics%sigma * rhoi * grav * thk0**2 / len0)**p1 * &
         2.0d0 * vis0 * model%numerics%dttem * tim0 / (16.0d0 * rhoi * shci)

    model%tempwk%dupc = (/ (model%numerics%sigma(2) - model%numerics%sigma(1)) / 2.0d0, &
         ((model%numerics%sigma(up+1) - model%numerics%sigma(up-1)) / 2.0d0, &
         up=2,model%general%upn-1), (model%numerics%sigma(model%general%upn) - &
         model%numerics%sigma(model%general%upn-1)) / 2.0d0  /)
    model%tempwk%dupa = (/ 0.0d0, 0.0d0, &
         ((model%numerics%sigma(up) - model%numerics%sigma(up-1)) / &
         ((model%numerics%sigma(up-2) - model%numerics%sigma(up-1)) * &
         (model%numerics%sigma(up-2) - model%numerics%sigma(up))), &
         up=3,model%general%upn)  /)
    model%tempwk%dupb = (/ 0.0d0, 0.0d0, &
         ((model%numerics%sigma(up) - model%numerics%sigma(up-2)) / &
         ((model%numerics%sigma(up-1) - model%numerics%sigma(up-2)) * &
         (model%numerics%sigma(up-1) - model%numerics%sigma(up))), &
         up=3,model%general%upn)  /)
    
!TODO - Delete line with tau0
    model%tempwk%f = (/ tim0 * coni / (thk0**2 * lhci * rhoi), &
         tim0 / (thk0 * lhci * rhoi), &
         tim0 * thk0 * rhoi * shci /  (thk0 * tim0 * model%numerics%dttem * lhci * rhoi), &
         tim0 * thk0**2 * vel0 * grav * rhoi / (4.0d0 * thk0 * len0 * rhoi * lhci), &
!!         tim0 * vel0 * tau0 / (4.0d0 * thk0 * rhoi * lhci) /)      
         0.d0 /)   !whl - last term no longer needed
         !*sfp* added the last term in the vect above for HO and SSA dissip. calc. 

    ! setting up some factors for sliding contrib to basal heat flux
    model%tempwk%slide_f = (/ VERT_DIFF * grav * thk0 * model%numerics%dttem/ shci, & ! vert diffusion
         VERT_ADV * rhoi*grav*acc0*thk0*thk0*model%numerics%dttem/coni /)             ! vert advection

    select case(model%options%whichbwat)
       case(0)
          model%paramets%hydtim = tim0 / (model%paramets%hydtim * scyr)
          estimate = 0.2d0 / model%paramets%hydtim
          !EIB! following not in lanl glide_temp
          call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 
          
          model%tempwk%c = (/ model%tempwk%dt_wat, 1.0d0 - 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, &
               1.0d0 + 0.5d0 * model%tempwk%dt_wat * model%paramets%hydtim, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /) 
       case(1)
          !EIB! not in lanl
          !model%tempwk%watvel = model%paramets%hydtim * tim0 / (scyr * len0)
          !estimate = (0.2d0 * model%tempwk%watvel) / min(model%numerics%dew,model%numerics%dns)
          !call find_dt_wat(model%numerics%dttem,estimate,model%tempwk%dt_wat,model%tempwk%nwat) 
          
          !print *, model%numerics%dttem*tim0/scyr, model%tempwk%dt_wat*tim0/scyr, model%tempwk%nwat

          !model%tempwk%c = (/ rhow * grav, rhoi * grav, 2.0d0 * model%numerics%dew, 2.0d0 * model%numerics%dns, &
          !     0.25d0 * model%tempwk%dt_wat / model%numerics%dew, 0.25d0 * model%tempwk%dt_wat / model%numerics%dns, &
          !     0.5d0 * model%tempwk%dt_wat / model%numerics%dew, 0.5d0 * model%tempwk%dt_wat / model%numerics%dns /)
          
       case(3)

          write(*,*)"ERROR: Current release does not support use of the basal processes module."
          write(*,*)"       Re-run code with alternate option selected for 'whichbwat'."
          stop

       end select
       ! JCC - was in LANL but not in parallel
       !      model%temper%temp = -10.0


       !MJH: Initialize ice temperature.
       !This block of code is identical to that in glissade_init_temp
!TODO - Remove hardwired constant (-273.15)
       if (model%temper%temp(1,1,1) < -273.15) then
           call write_log("No initial ice temperature supplied - setting temp to artm.")
           ! temp array still has initialized values - no values have been read in. 
           ! Initialize ice temperature to air temperature (for each column). 
           do ns = 1,model%general%nsn
              do ew = 1,model%general%ewn
                 model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns)))
              end do
           end do
       else
           ! Values have been read in - do nothing
       endif

       ! MJH: Calculate initial value of flwa
      ! If flwa is loaded (e.g. hotstart), use the flwa field in the input file instead
      ! Note: Implementing flwa initialization in this way, I don't think hotstart=1 does anything. 
!       if (model%options%hotstart  /=  1) then
       if (model%temper%flwa(1,1,1) < 0.d0) then
         call write_log("No initial flwa supplied - calculating initial flwa.")

!TODO - Check spelling of 'Glen', make sure it's consistent throughout code
         ! Calculate Glen's A --------------------------------------------------------   
         call calcflwa(model%numerics%sigma,        &
                       model%numerics%thklim,       &
                       model%temper%flwa,           &
                       model%temper%temp(:,1:model%general%ewn,1:model%general%nsn), &
                       model%geometry%thck,         &
                       model%paramets%flow_factor,  &
                       model%paramets%default_flwa, &
                       model%options%whichflwa) 
       endif
!       endif

  end subroutine glide_init_temp

!****************************************************    

!TODO - Replace model derived type with explicit arguments?
!       Could be a lot of work.

  subroutine glide_temp_driver(model,whichtemp,which_ho_diagnostic)

    !*FD Calculates the ice temperature, according to one
    !*FD of several alternative methods.

    use glimmer_utils, only: tridiag
    use glimmer_global, only : dp
    use glimmer_paramets, only : thk0, GLC_DEBUG
    use glide_velo
    use glide_thck
    use glide_grids
    use glide_bwater

    use parallel

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    type(glide_global_type),intent(inout) :: model                  ! model instance
    integer,                intent(in)    :: whichtemp              ! flag to choose method.
    integer,                intent(in)    :: which_ho_diagnostic 

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp),dimension(size(model%numerics%sigma)) :: subd, diag, supd, rhsd
    real(dp),dimension(size(model%numerics%sigma)) :: prevtemp, iteradvt, diagadvt
    real(dp) :: tempresid

    integer :: iter
    integer :: ew,ns

    real(dp),parameter :: tempthres = 0.001d0, floatlim = 10.0d0 / thk0
    integer, parameter :: mxit = 100
    integer, parameter :: ewbc = 1, nsbc = 1 

    real(dp), dimension(size(model%numerics%sigma)) :: weff

!WHLdebug
!!    integer :: i, j
!!    integer, parameter :: itest = 10, jtest = 10    

    !------------------------------------------------------------------------------------
    ! ewbc/nsbc set the type of boundary condition aplied at the end of
    ! the domain. a value of 0 implies zero gradient.
    !------------------------------------------------------------------------------------

    select case(whichtemp)

    !whl - In distributed code, case 0 is now handled by glissade_temp_driver.
    !      But leaving this option in glide_temp_driver in case it's ever needed.
 
    case(0) ! Set column to surface air temperature -------------------------------------

       ! JEFF - Ok for distributed since using air temperature at grid point to initialize.

       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             model%temper%temp(:,ew,ns) = dmin1(0.0d0,dble(model%climate%artm(ew,ns)))
          end do
       end do

    case(1) ! Do full temperature solution as in Glimmer---------------------------------

!TODO - Pretty sure these calls are not needed, since they are now at the end of tstep_p3
!       (so that wgrd can be written to the hotstart file and used for restart).
                         
       ! Calculate time-derivatives of thickness and upper surface elevation ------------

       call timeders(model%thckwk,   &
                     model%geometry%thck,     &
                     model%geomderv%dthckdtm, &
                     model%geometry%mask,     &
                     model%numerics%time,     &
                     1)

       call timeders(model%thckwk,   &
                     model%geometry%usrf,     &
                     model%geomderv%dusrfdtm, &
                     model%geometry%mask,     &
                     model%numerics%time,     &
                     2)

       ! Calculate the vertical velocity of the grid ------------------------------------
    
       ! Calculate the actual vertical velocity; method depends on whichwvel ------------

       ! *sfp* Added the if clause here so that this calc only gets done if has NOT
       ! *sfp* already been done at the end of the HO velocity calculation. If using HO
       ! *sfp* velocity calc, then this has already been done at the end of call to
       ! *sfp* 'run_ho_diagnostic' in glide_velo_higher. 
       if( which_ho_diagnostic == 0 )then

          call gridwvel(model%numerics%sigma,  &
                        model%numerics%thklim, &
                        model%velocity%uvel,   &
                        model%velocity%vvel,   &
                        model%geomderv,        &
                        model%geometry%thck,   &
                        model%velocity%wgrd)
    
          select case(model%options%whichwvel)

!TODO - Here and below, replace derived types (geomderv, numerics, etc.) with explicit arguments

          case(0) 
              ! Usual vertical integration
              call wvelintg(model%velocity%uvel,                        &
                            model%velocity%vvel,                        &
                            model%geomderv,                             &
                            model%numerics,                             &
                            model%velowk,                               &
                            model%velocity%wgrd(model%general%upn,:,:), &
                            model%geometry%thck,                        &
                            model%temper%bmlt,                          &
                            model%velocity%wvel)
    
          case(1)
              ! Vertical integration constrained so kinematic upper BC obeyed.
              call wvelintg(model%velocity%uvel,                        &
                            model%velocity%vvel,                        &
                            model%geomderv,                             &
                            model%numerics,                             &
                            model%velowk,                               &
                            model%velocity%wgrd(model%general%upn,:,:), &
                            model%geometry%thck,                        &
                            model%temper%  bmlt,                        &
                            model%velocity%wvel)

              call chckwvel(model%numerics,                             &
                            model%geomderv,                             &
                            model%velocity%uvel(1,:,:),                 &
                            model%velocity%vvel(1,:,:),                 &
                            model%velocity%wvel,                        &
                            model%geometry%thck,                        &
                            model%climate% acab)
          end select
       end if
       ! apply periodic ew BC
       if (model%options%periodic_ew) then
          call wvel_ew(model)
       end if

       model%tempwk%inittemp = 0.0d0
       model%tempwk%initadvt = 0.0d0
       !*MH model%tempwk%dissip   = 0.0d0  is also set to zero in finddisp

       ! ----------------------------------------------------------------------------------
!TODO - I think efvs is not needed here
       call glide_finddisp(model,          &
                           model%geometry%thck,     &
                           model%options%which_disp,&
                           model%stress%efvs, &
                           model%geomderv%stagthck, &
                           model%geomderv%dusrfdew, &
                           model%geomderv%dusrfdns, &
                           model%temper%flwa)

!WHLdebug
!!       i = itest
!!       j = jtest
!!       write(6,*) ' '
!!       write(6,*) 'i, j =', i, j
!!       write(6,*) 'which_disp =', model%options%which_disp
!!       write(6,*) ' '
!!       write(6,*) 'efvs =',   model%stress%efvs(:,i,j)
!!       write(6,*) ' '
!!       write(6,*) 'flwa =',   model%temper%flwa(:,i,j)
!!       write(6,*) ' '
!!       write(6,*) 'dissip =', model%tempwk%dissip(:,i,j)
!!       write(6,*) ' '
    
       !JCC - Don't use ho velocity fields unless we're using the ho model
       if (model%options%which_ho_diagnostic == 0 ) then
          ! translate velo field
          do ns = 2,model%general%nsn-1
              do ew = 2,model%general%ewn-1
                model%tempwk%hadv_u(:,ew,ns) = model%tempwk%advconst(1) * ( model%velocity%uvel(:,ew-1,ns-1) &
                    + model%velocity%uvel(:,ew-1,ns) + model%velocity%uvel(:,ew,ns-1) + model%velocity%uvel(:,ew,ns) )
                model%tempwk%hadv_v(:,ew,ns) = model%tempwk%advconst(2) * ( model%velocity%vvel(:,ew-1,ns-1) &
                    + model%velocity%vvel(:,ew-1,ns) + model%velocity%vvel(:,ew,ns-1) + model%velocity%vvel(:,ew,ns) )
              end do
          end do
      else ! Using ho physics
          ! translate velo field
          do ns = 2,model%general%nsn-1
              do ew = 2,model%general%ewn-1
                model%tempwk%hadv_u(:,ew,ns) = model%tempwk%advconst(1) * ( model%velocity%uvel(:,ew-1,ns-1) &
                    + model%velocity%uvel(:,ew-1,ns) + model%velocity%uvel(:,ew,ns-1) + model%velocity%uvel(:,ew,ns) )
                model%tempwk%hadv_v(:,ew,ns) = model%tempwk%advconst(2) * ( model%velocity%vvel(:,ew-1,ns-1) &
                    + model%velocity%vvel(:,ew-1,ns) + model%velocity%vvel(:,ew,ns-1) + model%velocity%vvel(:,ew,ns) )
              end do
          end do
      end if ! model%options%which_ho_diagnostic == 0

       call hadvall(model, &
            model%temper%temp, &
            model%geometry%thck)

       ! zeroth iteration
       iter = 0
       tempresid = 0.0d0

       do ns = 2,model%general%nsn-1
          do ew = 2,model%general%ewn-1
             if(model%geometry%thck(ew,ns) > model%numerics%thklim) then

                weff = model%velocity%wvel(:,ew,ns) - model%velocity%wgrd(:,ew,ns)

!TODO - It seems odd to zero out weff when it's big.  Why not set to wmax?
                if (maxval(abs(weff)) > model%tempwk%wmax) then
                   weff = 0.0d0
                end if

                call hadvpnt(iteradvt,                       &
                             diagadvt,                               &
                             model%temper%temp(:,ew-2:ew+2,ns),      &
                             model%temper%temp(:,ew,ns-2:ns+2),      &
                             model%tempwk%hadv_u(:,ew,ns), &
                             model%tempwk%hadv_v(:,ew,ns))
               
                call findvtri(model,ew,ns,subd,diag,supd,diagadvt, &
                     weff, &
                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                call findvtri_init(model,ew,ns,subd,diag,supd,weff,model%temper%temp(:,ew,ns), &
                     model%geometry%thck(ew,ns),GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                call findvtri_rhs(model,ew,ns,model%climate%artm(ew,ns),iteradvt,rhsd, &
                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                prevtemp(:) = model%temper%temp(:,ew,ns)

                call tridiag(subd(1:model%general%upn), &
                             diag(1:model%general%upn), &
                             supd(1:model%general%upn), &
                             model%temper%temp(1:model%general%upn,ew,ns), &
                             rhsd(1:model%general%upn))

                call corrpmpt(model%temper%temp(:,ew,ns),     &
                              model%geometry%thck(ew,ns),     &
                              model%temper%bwat(ew,ns),       &
                              model%numerics%sigma,           &
                              model%general%upn)

                tempresid = max(tempresid,maxval(abs(model%temper%temp(:,ew,ns)-prevtemp(:))))

             endif
          end do
       end do

       do while (tempresid > tempthres .and. iter <= mxit)
          tempresid = 0.0d0

          do ns = 2,model%general%nsn-1
             do ew = 2,model%general%ewn-1
                if(model%geometry%thck(ew,ns) > model%numerics%thklim) then

                   weff = model%velocity%wvel(:,ew,ns) - model%velocity%wgrd(:,ew,ns)
                   if (maxval(abs(weff)) > model%tempwk%wmax) then
                      weff = 0.0d0
                   end if

                   call hadvpnt(iteradvt,                       &
                                diagadvt,                               &
                                model%temper%temp(:,ew-2:ew+2,ns),      &
                                model%temper%temp(:,ew,ns-2:ns+2),      &
                                model%tempwk%hadv_u(:,ew,ns), &
                                model%tempwk%hadv_v(:,ew,ns))

                   call findvtri(model,ew,ns,subd,diag,supd,diagadvt, &
                                 weff, &
                                 GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                   call findvtri_rhs(model,ew,ns,model%climate%artm(ew,ns),iteradvt,rhsd, &
                                     GLIDE_IS_FLOAT(model%geometry%thkmask(ew,ns)))

                   prevtemp(:) = model%temper%temp(:,ew,ns)

                   call tridiag(subd(1:model%general%upn), &
                                diag(1:model%general%upn), &
                                supd(1:model%general%upn), &
                                model%temper%temp(1:model%general%upn,ew,ns), &
                                rhsd(1:model%general%upn))

                   call corrpmpt(model%temper%temp(:,ew,ns),     &
                                 model%geometry%thck(ew,ns),     &
                                 model%temper%bwat(ew,ns),       &
                                 model%numerics%sigma,           &
                                 model%general%upn)


                   tempresid = max(tempresid, maxval(abs(model%temper%temp(:,ew,ns)-prevtemp(:))))

                endif
             end do
          end do

          iter = iter + 1
       end do

       model%temper%niter = max(model%temper%niter, iter )
       
       ! set temperature of thin ice to the air temperature and set ice free nodes to zero
       do ns = 1,model%general%nsn
          do ew = 1,model%general%ewn
             if (GLIDE_IS_THIN(model%geometry%thkmask(ew,ns))) then
                model%temper%temp(:,ew,ns) = min(0.0d0,dble(model%climate%artm(ew,ns)))
             else if (GLIDE_NO_ICE(model%geometry%thkmask(ew,ns))) then
                !JCC - one of these
                !model%temper%temp(:,ew,ns) = min(0.0d0,dble(model%climate%artm(ew,ns)))
                model%temper%temp(:,ew,ns) = 0.0d0
             end if
          end do
       end do

       ! apply periodic ew BC
       if (model%options%periodic_ew) then
          model%temper%temp(:,0,:) = model%temper%temp(:,model%general%ewn-2,:)
          model%temper%temp(:,1,:) = model%temper%temp(:,model%general%ewn-1,:)
          model%temper%temp(:,model%general%ewn,:) = model%temper%temp(:,2,:)
          model%temper%temp(:,model%general%ewn+1,:) = model%temper%temp(:,3,:)
       end if

       ! Calculate basal melt rate --------------------------------------------------

       call glide_calcbmlt(model, &
                           model%options%which_bmelt, &
                           model%temper%temp, &
                           model%geometry%thck, &
                           model%geomderv%stagthck, &
                           model%geomderv%dusrfdew, &
                           model%geomderv%dusrfdns, &
                           model%velocity%ubas, &
                           model%velocity%vbas, &
                           model%temper%bmlt, &
                           GLIDE_IS_FLOAT(model%geometry%thkmask))

       ! Calculate basal water depth ------------------------------------------------

       call calcbwat(model, &
                     model%options%whichbwat, &
                     model%temper%bmlt, &
                     model%temper%bwat, &
                     model%temper%bwatflx, &
                     model%geometry%thck, &
                     model%geometry%topg, &
                     model%temper%temp(model%general%upn,:,:), &
                     GLIDE_IS_FLOAT(model%geometry%thkmask), &
                     model%tempwk%wphi)

!HALO - If bwat is computed only on local cells, do we need a halo update here?

    case(2) ! *sfp* stealing this un-used option ... 

        ! DO NOTHING. That is, hold T const. at initially assigned value
        !TODO whl - Should the do-nothing option have a different case number, such as 0 or -1? 

    end select   ! whichtemp

    ! Calculate Glen's A --------------------------------------------------------

    call calcflwa(model%numerics%sigma,        &
                  model%numerics%thklim,       &
                  model%temper%flwa,           &
                  model%temper%temp(:,1:model%general%ewn,1:model%general%nsn), &
                  model%geometry%thck,         &
                  model%paramets%flow_factor,  &
                  model%paramets%default_flwa, &
                  model%options%whichflwa) 

    ! Output some information ----------------------------------------------------

    if (GLC_DEBUG) then
       print *, "* temp ", model%numerics%time, iter, model%temper%niter, &
            real(model%temper%temp(model%general%upn,model%general%ewn/2+1,model%general%nsn/2+1))
    end if

  end subroutine glide_temp_driver

  !-------------------------------------------------------------------------

  subroutine hadvpnt(iteradvt,diagadvt,tempx,tempy,u,v)

    use glimmer_global, only : dp

    real(dp), dimension(:), intent(in) :: u,v
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    if (u(1) > 0.0d0) then
       iteradvt = u * (- 4.0d0*tempx(:,2) + tempx(:,1))
       diagadvt = u * 3.0d0
    else if (u(1) < 0.0d0) then
       iteradvt = u * (4.0d0*tempx(:,4) - tempx(:,5))
       diagadvt = - u * 3.0d0
    end if

    if (v(1) > 0.0d0) then
       iteradvt = iteradvt + v * (- 4.0d0*tempy(:,2) + tempy(:,1))
       diagadvt = diagadvt + v * 3.0d0
    else if (v(1) < 0.0d0) then
       iteradvt = iteradvt + v * (4.0d0*tempy(:,4) - tempy(:,5))
       diagadvt = diagadvt - v * 3.0d0
    end if

  end subroutine hadvpnt

  !-------------------------------------------------------------------------

  subroutine fohadvpnt(tempwk,iteradvt,diagadvt,tempx,tempy,uvel,vvel)

    use glimmer_global, only : dp
    use glimmer_utils, only: hsum

    type(glide_tempwk) :: tempwk
    real(dp), dimension(:,:,:), intent(in) :: uvel, vvel
    real(dp), dimension(:,:), intent(in) :: tempx, tempy
    real(dp), dimension(:), intent(out) :: iteradvt, diagadvt

    real(dp), dimension(size(iteradvt)) :: u, v

    iteradvt = 0.0d0
    diagadvt = 0.0d0

    u = tempwk%advconst(1) * hsum(uvel(:,:,:))
    v = tempwk%advconst(2) * hsum(vvel(:,:,:))

    if (u(1) > 0.0d0) then
       iteradvt = - u * 2.0d0 * tempx(:,1)
       diagadvt = 2.0d0 * u 
    else if (u(1) < 0.0d0) then
       iteradvt = u * 2.0d0 * tempx(:,3)
       diagadvt = - 2.0d0 * u 
    end if

    if (v(1) > 0.0d0) then
       iteradvt = iteradvt - v * 2.0d0 * tempy(:,1) 
       diagadvt = diagadvt + 2.0d0 * v 
    else if (v(1) < 0.0d0) then
       iteradvt = iteradvt + v * 2.0d0 * tempy(:,3)
       diagadvt = diagadvt - 2.0d0 * v 
    end if

  end subroutine fohadvpnt

  !-------------------------------------------------------------------------

  subroutine hadvall(model,temp,thck)

    use glimmer_global, only : dp 

    type(glide_global_type) :: model
    real(dp), dimension(:,0:,0:), intent(in) :: temp
    real(dp), dimension(:,:), intent(in) :: thck

    real(dp), dimension(size(temp,dim=1)) :: diagadvt

    integer :: ew,ns

    model%tempwk%initadvt = 0.0d0

    do ns = 2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then

             call hadvpnt(model%tempwk%initadvt(:,ew,ns), &
                  diagadvt,                       &
                  temp(:,ew-2:ew+2,ns),           &
                  temp(:,ew,ns-2:ns+2),           &
                  model%tempwk%hadv_u(:,ew,ns), &
                  model%tempwk%hadv_v(:,ew,ns))
          end if
       end do
    end do

  end subroutine hadvall

  !-------------------------------------------------------------------------

  subroutine findvtri(model,ew,ns,subd,diag,supd,diagadvt,weff,float)

    use glimmer_global, only : dp

    type(glide_global_type) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(in) :: weff,  diagadvt
    real(dp), dimension(:), intent(out) :: subd, diag, supd
    logical, intent(in) :: float

    real(dp) :: fact(3)

! These constants are precomputed:
! model%tempwk%cons(1) = 2.0d0 * tim0 * model%numerics%dttem * coni / (2.0d0 * rhoi * shci * thk0**2)
! model%tempwk%cons(2) = model%numerics%dttem / 2.0d0 

    fact(1) = VERT_DIFF*model%tempwk%cons(1) / model%geometry%thck(ew,ns)**2
    fact(2) = VERT_ADV*model%tempwk%cons(2) / model%geometry%thck(ew,ns)    
    
    subd(2:model%general%upn-1) = fact(2) * weff(2:model%general%upn-1) * &
         model%tempwk%dups(2:model%general%upn-1,3)

    supd(2:model%general%upn-1) = - subd(2:model%general%upn-1) - fact(1) * &
         model%tempwk%dups(2:model%general%upn-1,2)

    subd(2:model%general%upn-1) = subd(2:model%general%upn-1) - fact(1) * &
         model%tempwk%dups(2:model%general%upn-1,1)

    diag(2:model%general%upn-1) = 1.0d0 - subd(2:model%general%upn-1) &
         - supd(2:model%general%upn-1) &
         + diagadvt(2:model%general%upn-1)

    supd(1) = 0.0d0
    subd(1) = 0.0d0
    diag(1) = 1.0d0

    ! now do the basal boundary
    ! for grounded ice, a heat flux is applied
    ! for floating ice, temperature held constant

     if (float) then

       supd(model%general%upn) = 0.0d0
       subd(model%general%upn) = 0.0d0
       diag(model%general%upn) = 1.0d0

     else 
 
        supd(model%general%upn) = 0.0d0 
        subd(model%general%upn) = -0.5*fact(1)/(model%tempwk%dupn**2)
        diag(model%general%upn) = 1.0d0 - subd(model%general%upn) + diagadvt(model%general%upn)
 
     end if

  end subroutine findvtri

  !-------------------------------------------------------------------------

  subroutine findvtri_init(model,ew,ns,subd,diag,supd,weff,temp,thck,float)
    !*FD called during first iteration to set inittemp
    use glimmer_global, only : dp
    use glimmer_paramets, only: vel0, vel_scale

    type(glide_global_type) :: model
    integer, intent(in) :: ew, ns
    real(dp), dimension(:), intent(in) :: temp,diag,subd,supd,weff
    real(dp), intent(in) :: thck
    logical, intent(in) :: float    

    ! local variables
    real(dp) :: slterm
    integer ewp,nsp
    integer slide_count

    model%tempwk%inittemp(2:model%general%upn-1,ew,ns) = temp(2:model%general%upn-1) * &
         (2.0d0 - diag(2:model%general%upn-1)) &
         - temp(1:model%general%upn-2) * subd(2:model%general%upn-1) &
         - temp(3:model%general%upn) * supd(2:model%general%upn-1) & 
         - model%tempwk%initadvt(2:model%general%upn-1,ew,ns) &
         + model%tempwk%dissip(2:model%general%upn-1,ew,ns)
    
    if (float) then
       model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) 
       !EIB old!model%tempwk%inittemp(model%general%upn,ew,ns) = pmpt(thck)
    else 
       ! sliding contribution to basal heat flux
       slterm = 0.
       slide_count = 0

       !whl - BUG! - The following expression for taub*ubas is valid only for the SIA
       !             Need a different expression for HO dynamics

       ! only include sliding contrib if temperature node is surrounded by sliding velo nodes
       do nsp = ns-1,ns
          do ewp = ew-1,ew

!SCALING - WHL: Multiply ubas by vel0/vel_scale so we get the same result in these two cases:
!           (1) Old Glimmer with scaling:         vel0 = vel_scale = 500/scyr, and ubas is non-dimensional
!           (2) New Glimmer-CISM without scaling: vel0 = 1, vel_scale = 500/scyr, and ubas is in m/s.

!!!             if ( abs(model%velocity%ubas(ewp,nsp)) > 0.000001 .or. &
!!!                  abs(model%velocity%vbas(ewp,nsp)) > 0.000001 ) then
             if ( abs(model%velocity%ubas(ewp,nsp))*(vel0/vel_scale) > 1.d-6 .or. &
                  abs(model%velocity%vbas(ewp,nsp))*(vel0/vel_scale) > 1.d-6 ) then

                slide_count = slide_count + 1
                slterm = slterm + (&
                     model%geomderv%dusrfdew(ewp,nsp) * model%velocity%ubas(ewp,nsp) + &
                     model%geomderv%dusrfdns(ewp,nsp) * model%velocity%vbas(ewp,nsp))
             end if
          end do
       end do
       if (slide_count >= 4) then
          slterm = 0.25*slterm
       else
          slterm = 0.
       end if
       model%tempwk%inittemp(model%general%upn,ew,ns) = temp(model%general%upn) * &
            (2.0d0 - diag(model%general%upn)) &
            - temp(model%general%upn-1) * subd(model%general%upn) &
            - 0.5*model%tempwk%cons(3) * model%temper%bheatflx(ew,ns) / (thck * model%tempwk%dupn) & ! geothermal heat flux (diff)
            - model%tempwk%slide_f(1)*slterm/ model%tempwk%dupn &                                    ! sliding heat flux    (diff)
            - model%tempwk%cons(4) * model%temper%bheatflx(ew,ns) * weff(model%general%upn) &        ! geothermal heat flux (adv)
            - model%tempwk%slide_f(2)*thck*slterm* weff(model%general%upn) &                         ! sliding heat flux    (adv)
            - model%tempwk%initadvt(model%general%upn,ew,ns)  &
            + model%tempwk%dissip(model%general%upn,ew,ns)
    end if

  end subroutine findvtri_init

  !-----------------------------------------------------------------------

  subroutine findvtri_rhs(model,ew,ns,artm,iteradvt,rhsd,float)

    !*FD RHS of temperature tri-diag system
    use glimmer_global, only : dp, sp 

    type(glide_global_type) :: model
    integer, intent(in) :: ew, ns
    real(sp), intent(in) :: artm 
    real(dp), dimension(:), intent(in) :: iteradvt
    real(dp), dimension(:), intent(out) :: rhsd
    logical, intent(in) :: float    

    ! upper boundary condition
    rhsd(1) = artm
    if (float) then
       rhsd(model%general%upn) = model%tempwk%inittemp(model%general%upn,ew,ns)    
    else
       rhsd(model%general%upn) = model%tempwk%inittemp(model%general%upn,ew,ns) - iteradvt(model%general%upn)
    end if
    rhsd(2:model%general%upn-1) = model%tempwk%inittemp(2:model%general%upn-1,ew,ns) - iteradvt(2:model%general%upn-1)

  end subroutine findvtri_rhs

  !-----------------------------------------------------------------------

  subroutine glide_calcbmlt(model,whichbmelt,temp,thck,stagthck,dusrfdew,dusrfdns,ubas,vbas,bmlt,floater)

    use glimmer_global, only : dp 

    type(glide_global_type) :: model
    real(dp), dimension(:,0:,0:), intent(in) :: temp
    real(dp), dimension(:,:), intent(in) :: thck,  stagthck, dusrfdew, dusrfdns, ubas, vbas  
    real(dp), dimension(:,:), intent(inout) :: bmlt
    logical, dimension(:,:), intent(in) :: floater
    integer, intent(in) :: whichbmelt

    real(dp), dimension(size(model%numerics%sigma)) :: pmptemp
    real(dp) :: slterm, newmlt
 
    integer :: ewp, nsp, up, ew, ns

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim .and. .not. floater(ew,ns)) then

             call calcpmpt(pmptemp,thck(ew,ns),model%numerics%sigma)

             if (abs(temp(model%general%upn,ew,ns)-pmptemp(model%general%upn)) < 0.001) then

                slterm = 0.0d0

                select case( whichbmelt )    !*sfp* added for calculating differently based on model physics

                case( SIA_BMELT )                   ! 0-order SIA approx. --> Tau_d = Tau_b                                     

                    do nsp = ns-1,ns
                        do ewp = ew-1,ew
                            slterm = slterm - stagthck(ewp,nsp) * &
                            (dusrfdew(ewp,nsp) * ubas(ewp,nsp) + dusrfdns(ewp,nsp) * vbas(ewp,nsp))
                        end do
                    end do

                    !*sfp* NOTE that multiplication by this term has been moved up from below
                    slterm = model%tempwk%f(4) * slterm 

!TODO - Remove this case.  Commenting out for now.

!!                case( FIRSTORDER_BMELT )                   ! 1st-order SIA approx. (HO model)
!!                    do nsp = ns-1,ns 
!!                        do ewp = ew-1,ew
!!                        !*sfp* Note that vel and stress arrays have diff vert dims (upn for vel, upn-1 for stress)
!!                        ! so that for now, basal vel at upn is multiplied by basal stress at upn-1 to get frictional
!!                        ! heating term. This may not be entirely correct ... 
!!                             slterm = slterm - &
                                
!!                                !! NEW version: uses consistent basal tractions                
!!                                ( -model%velocity%btraction(1,ewp,nsp) * &
!!                                   model%velocity%uvel(model%general%upn,ewp,nsp)  &
!!                                  -model%velocity%btraction(2,ewp,nsp) * &
!!                                   model%velocity%vvel(model%general%upn,ewp,nsp) )
                                
!!                                !!!! OLD version: uses HO basal shear stress calc. from FD                
!!                                !( model%stress%tau%xz(model%general%upn-1,ewp,nsp) * &
!!                                !  model%stress%uvel(model%general%upn,ewp,nsp) + &
!!                                !  model%stress%tau%yz(model%general%upn-1,ewp,nsp) * &
!!                                !  model%stress%vvel(model%general%upn,ewp,nsp) )

!!                        end do
!!                    end do

!!                    slterm = model%tempwk%f(5) * slterm

!                case( SSA_BMELT )                  ! 1st-order, depth-integrated approx. (SSA) 
!                                            ! NOTE: need to pass 2d basal shear stress arrays from SSA model
!                   do nsp = ns-1,ns 
!                         do ewp = ew-1,ew
!                             slterm = slterm - &
!                                 ( taubxs(upn,ewp,nsp) * ubas(ewp,nsp) + taubys(upn,ewp,nsp) * vbas(ewp,nsp))
!                         end do
!                   end do
!
!                   slterm = model%tempwk%f(5) * slterm

                end select

                bmlt(ew,ns) = 0.0d0

                !*sfp* changed this so that 'slterm' is multiplied by f(4) const. above ONLY for the 0-order SIA case,
                ! since for the HO and SSA cases a diff. const. needs to be used

                ! OLD version
!                newmlt = model%tempwk%f(4) * slterm - model%tempwk%f(2)*model%temper%bheatflx(ew,ns) + model%tempwk%f(3) * &
!                     model%tempwk%dupc(model%general%upn) * &
!                     thck(ew,ns) * model%tempwk%dissip(model%general%upn,ew,ns)

                ! NEW version (sfp)
                newmlt = slterm - model%tempwk%f(2)*model%temper%bheatflx(ew,ns)   &
                        + model%tempwk%f(3) * model%tempwk%dupc(model%general%upn) * &
                          thck(ew,ns) * model%tempwk%dissip(model%general%upn,ew,ns)

                up = model%general%upn - 1

!TODO - Change reals to DP
                do while (abs(temp(up,ew,ns)-pmptemp(up)) < 0.001 .and. up >= 3)
                   bmlt(ew,ns) = bmlt(ew,ns) + newmlt
                   newmlt = model%tempwk%f(3) * model%tempwk%dupc(up) * thck(ew,ns) * model%tempwk%dissip(up,ew,ns)
                   up = up - 1
                end do

                up = up + 1

                if (up == model%general%upn) then
                   bmlt(ew,ns) = newmlt - &
                        model%tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * model%tempwk%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * model%tempwk%dupb(up) ) / thck(ew,ns) 
                else
                   bmlt(ew,ns) = bmlt(ew,ns) + max(0.0d0, newmlt - &
                        model%tempwk%f(1) * ( (temp(up-2,ew,ns) - pmptemp(up-2)) * model%tempwk%dupa(up) &
                        + (temp(up-1,ew,ns) - pmptemp(up-1)) * model%tempwk%dupb(up) ) / thck(ew,ns)) 
                end if

             else

                bmlt(ew,ns) = 0.0d0

             end if

          !EIB! else if (model%options%use_plume == 1) then

          ! do nothing because the plume model will have written the bmlt field
          else

              bmlt(ew,ns) = 0.0d0

          end if
       end do
    end do

    ! apply periodic BC

    if (model%options%periodic_ew) then
       do ns = 2,model%general%nsn-1
          bmlt(1,ns) = bmlt(model%general%ewn-1,ns)
          bmlt(model%general%ewn,ns) = bmlt(2,ns)
       end do
    end if
  end subroutine glide_calcbmlt

!-------------------------------------------------------------------

  subroutine glide_finddisp(model,thck,whichdisp,efvs,stagthck,dusrfdew,dusrfdns,flwa)

    ! Compute the dissipation source term associated with strain heating.
    ! Note that the dissipation is computed in the same way on either a staggered or an
    !  unstaggered vertical grid.  
    ! Note also that dissip and flwa must have the same vertical dimension 
    !  (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
    
    use glimmer_global, only : dp
    use glimmer_physcon, only : gn

    type(glide_global_type) :: model
    real(dp), dimension(:,:), intent(in) :: thck, stagthck, dusrfdew, dusrfdns
    real(dp), dimension(:,:,:), intent(in) :: flwa, efvs
    integer, intent(in) :: whichdisp

    integer, parameter :: p1 = gn + 1  
    integer :: ew, ns
    integer :: iew, ins    !*sfp* for HO and SSA dissip. calc.

    real(dp) :: c2 

    !*sfp* The next 2 declarations needed for HO and SSA dissip. calc. ... only needed 
    ! for internal work, so not clear if it is necessary to declare/allocate them elsewhere 
    real(dp) :: c4                         
    real(dp), dimension(model%general%upn) :: c5     
    
    select case( whichdisp ) 

    case( SIA_DISP )
    !*sfp* 0-order SIA case only 
    ! two methods of doing this. 
    ! 1. find dissipation at u-pts and then average
    ! 2. find dissipation at H-pts by averaging quantities from u-pts
    ! 2. works best for eismint divide (symmetry) but 1 likely to be better for full expts

    model%tempwk%dissip(:,:,:) = 0.0d0

    do ns = 2, model%general%nsn-1
       do ew = 2, model%general%ewn-1
          if (thck(ew,ns) > model%numerics%thklim) then
             
             c2 = (0.25*sum(stagthck(ew-1:ew,ns-1:ns)) * dsqrt((0.25*sum(dusrfdew(ew-1:ew,ns-1:ns)))**2 &
                  + (0.25*sum(dusrfdns(ew-1:ew,ns-1:ns)))**2))**p1
             
             model%tempwk%dissip(:,ew,ns) = c2 * model%tempwk%c1(:) * ( &
                  flwa(:,ew-1,ns-1) + flwa(:,ew-1,ns+1) + flwa(:,ew+1,ns+1) + flwa(:,ew+1,ns-1) + &
                  2*(flwa(:,ew-1,ns)+flwa(:,ew+1,ns)+flwa(:,ew,ns-1)+flwa(:,ew,ns+1)) + &
                  4*flwa(:,ew,ns)) 

          end if
       end do
    end do

    !the compensatory heating(compheat) is initialized to zero and allows
    !for modifying the calculated dissip.  This is needed for exact verification tests,
    !model%tempwk%dissip = model%tempwk%dissip + model%tempwk%compheat 

    case( FIRSTORDER_DISP )

!TODO - Print an error message and exit gracefully

!   case( SSA_DISP )     !!! Waiting for an SSA solver !!!
!    !*sfp* 1st-order, depth-integrated case only (SSA model) 
!    ! NOTE: this needs taus and efvss (2d arrays), which are depth-integrated and averaged 
!    ! effective stress and effective viscosity fields calculated from and/or consistent
!    ! with the SSA model.
!
!    model%tempwk%dissip = 0.0d0
!    do ns = 2, model%general%nsn-1
!       do ew = 2, model%general%ewn-1
!          if (thck(ew,ns) > model%numerics%thklim) then
!            c4 = 0.0d0
!            do ins = ns-1,ns; do iew = ew-1,ew; 
!                if (efvss(iew,ins)  /=  0.0d0) then                     
!                    c4 = c4 + taus(iew,ins)**2 / efvss(iew,ins)
!                end if; 
!            end do; end do
!            model%tempwk%dissip(:,ew,ns) = c4 * model%tempwk%cons(5)
!          end if
!       end do
!    end do

    end select

  end subroutine glide_finddisp

  !-----------------------------------------------------------------------------------

  subroutine corrpmpt(temp,thck,bwat,sigma,upn)

    use glimmer_global, only : dp

    real(dp), dimension(:), intent(inout) :: temp
    real(dp), intent(in) :: thck, bwat
    integer,intent(in) :: upn
    real(dp),dimension(:),intent(in) :: sigma

    real(dp), dimension(:) :: pmptemp(size(temp))

    ! corrects a temperature column for melting point effects
    ! 1. if temperature at any point in column is above pressure melting point then 
    !    set temperature to pressure melting point 
    ! 2. if bed is wet set basal temperature to pressure melting point 

    call calcpmpt(pmptemp,thck,sigma)

    temp = dmin1(temp,pmptemp)

    if (bwat > 0.0d0) temp(upn) = pmptemp(upn)

  end subroutine corrpmpt

  !-------------------------------------------------------------------

  subroutine calcpmpt(pmptemp,thck,sigma)

    use glimmer_global, only : dp !, upn
    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), dimension(:), intent(out) :: pmptemp
    real(dp), intent(in) :: thck
    real(dp),intent(in),dimension(:) :: sigma

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp(:) = fact * thck * sigma(:)

  end subroutine calcpmpt

  !-------------------------------------------------------------------

  subroutine calcpmptb(pmptemp,thck)

    use glimmer_global, only : dp
    use glimmer_physcon, only : rhoi, grav, pmlt 
    use glimmer_paramets, only : thk0

    real(dp), intent(out) :: pmptemp
    real(dp), intent(in) :: thck

    real(dp), parameter :: fact = - grav * rhoi * pmlt * thk0

    pmptemp = fact * thck 

  end subroutine calcpmptb

  !-------------------------------------------------------------------

  subroutine calcflwa(sigma, thklim, flwa, temp, thck, flow_factor, default_flwa_arg, flag)

    !*FD Calculates Glen's $A$ over the three-dimensional domain,
    !*FD using one of three possible methods.

    use glimmer_physcon
    use glimmer_paramets, only : thk0, vis0

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

!      Note: dummy argument sigma can be either model%numerics%sigma or model%numerics%stagsigma.
!      The flwa, temp, and sigma arrays must have the same vertical dimension
!       (1:upn on an unstaggered vertical grid, or 1:upn-1 on a staggered vertical grid).
!
    real(dp),dimension(:),      intent(in)    :: sigma     !*FD Vertical coordinate
    real(dp),                   intent(in)    :: thklim    !*FD thickness threshold
    real(dp),dimension(:,:,:),  intent(out)   :: flwa      !*FD The calculated values of $A$
!whl - Changed this so that the temp array is assumed to start with horizontal index 1 instead of 0.
!!    real(dp),dimension(:,0:,0:),intent(in)    :: temp      !*FD The 3D temperature field
    real(dp),dimension(:,:,:),intent(in)      :: temp      !*FD The 3D temperature field
    real(dp),dimension(:,:),    intent(in)    :: thck      !*FD The ice thickness
    real(dp)                                  :: flow_factor !*FD Fudge factor in arrhenius relationship
    real(dp),                   intent(in)    :: default_flwa_arg !*FD Glen's A to use in isothermal case 
    integer,                    intent(in)    :: flag      !*FD Flag to select the method
                                                           !*FD of calculation:
    !*FD \begin{description}
    !*FD \item[0] {\em Paterson and Budd} relationship.
    !*FD \item[1] {\em Paterson and Budd} relationship, with temperature set to
    !*FD -5$^{\circ}$C.
    !*FD \item[2] Set constant, {\em but not sure how this works at the moment\ldots}
    !*FD \end{description}

    !------------------------------------------------------------------------------------
    ! Internal variables
    !------------------------------------------------------------------------------------

    real(dp), parameter :: fact = grav * rhoi * pmlt * thk0
    real(dp), parameter :: contemp = -5.0d0
    real(dp) :: default_flwa
    real(dp),dimension(4) :: arrfact
    integer :: ew,ns,up,ewn,nsn,upn

!    real(dp), dimension(size(sigma)) :: tempcor
!   KJE give a hard number to satisfy gnu complier 
    real(dp), dimension(1000) :: tempcor

    !------------------------------------------------------------------------------------ 
   
    default_flwa = flow_factor * default_flwa_arg / (vis0*scyr) 

    upn=size(flwa,1) ; ewn=size(flwa,2) ; nsn=size(flwa,3)

    !write(*,*)"Default flwa = ",default_flwa

    arrfact = (/ flow_factor * arrmlh / vis0, &   ! Value of a when T* is above -263K
                 flow_factor * arrmll / vis0, &   ! Value of a when T* is below -263K
                 -actenh / gascon,        &       ! Value of -Q/R when T* is above -263K
                 -actenl / gascon/)               ! Value of -Q/R when T* is below -263K
 
    select case(flag)
    case(FLWA_PATERSON_BUDD)

      ! This is the Paterson and Budd relationship

      do ns = 1,nsn
        do ew = 1,ewn
          if (thck(ew,ns) > thklim) then
            
            ! Calculate the corrected temperature

            tempcor(:) = min(0.0d0, temp(:,ew,ns) + thck(ew,ns) * fact * sigma(:))
            tempcor(:) = max(-50.0d0, tempcor(:))

            ! Calculate Glen's A

!TODO - This call and the one below could be inlined.

            call patebudd(tempcor(:), flwa(:,ew,ns), arrfact) 

          else
            flwa(:,ew,ns) = default_flwa
          end if
        end do
      end do

    case(FLWA_PATERSON_BUDD_CONST_TEMP)

      ! This is the Paterson and Budd relationship, but with the temperature held constant
      ! at -5 deg C

      do ns = 1,nsn
        do ew = 1,ewn
          if (thck(ew,ns) > thklim) then

            ! Calculate Glen's A with a fixed temperature.

            call patebudd((/(contemp, up=1,upn)/),flwa(:,ew,ns),arrfact) 

          else
            flwa(:,ew,ns) = default_flwa
          end if
        end do
      end do

    case(FLWA_CONST_FLWA) 

      flwa(:,:,:) = default_flwa
  
    end select

  end subroutine calcflwa 

!------------------------------------------------------------------------------------------

  subroutine patebudd(tempcor,calcga,fact)

    !*FD Calculates the value of Glen's $A$ for the temperature values in a one-dimensional
    !*FD array. The input array is usually a vertical temperature profile. The equation used
    !*FD is from \emph{Paterson and Budd} [1982]:
    !*FD \[
    !*FD A(T^{*})=a \exp \left(\frac{-Q}{RT^{*}}\right)
    !*FD \]
    !*FD This is equation 9 in {\em Payne and Dongelmans}. $a$ is a constant of proportionality,
    !*FD $Q$ is the activation energy for for ice creep, and $R$ is the universal gas constant.
    !*FD The pressure-corrected temperature, $T^{*}$ is given by:
    !*FD \[
    !*FD T^{*}=T-T_{\mathrm{pmp}}+T_0
    !*FD \] 
    !*FD \[
    !*FD T_{\mathrm{pmp}}=T_0-\sigma \rho g H \Phi
    !*FD \]
    !*FD $T$ is the ice temperature, $T_{\mathrm{pmp}}$ is the pressure melting point 
    !*FD temperature, $T_0$ is the triple point of water, $\rho$ is the ice density, and 
    !*FD $\Phi$ is the (constant) rate of change of melting point temperature with pressure.

    use glimmer_physcon, only : trpt

    !------------------------------------------------------------------------------------
    ! Subroutine arguments
    !------------------------------------------------------------------------------------

    real(dp),dimension(:), intent(in)    :: tempcor  !*FD Input temperature profile. This is 
                                                     !*FD {\em not} $T^{*}$, as it has $T_0$
                                                     !*FD added to it later on; rather it is
                                                     !*FD $T-T_{\mathrm{pmp}}$.
    real(dp),dimension(:), intent(out)   :: calcga   !*FD The output values of Glen's $A$.
    real(dp),dimension(4), intent(in)    :: fact     !*FD Constants for the calculation. These
                                                     !*FD are set when the velo module is initialised

    !------------------------------------------------------------------------------------

    ! Actual calculation is done here - constants depend on temperature -----------------

    where (tempcor >= -10.0d0)         
      calcga = fact(1) * exp(fact(3) / (tempcor + trpt))
    elsewhere
      calcga = fact(2) * exp(fact(4) / (tempcor + trpt))
    end where

  end subroutine patebudd

end module glide_temp

