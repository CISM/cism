!CLEANUP - New module based on glissade_velo.F90.
! For now, it is simply a wrapper to the Payne-Price dycore.
! Later, it will also be a wrapper to the variational dycore.
! Should probably make the JFNK solver independent of the underlying
!  velocity solver (e.g., Payne-Price v. variational)
! Note: glam_velo_fordsiapstr is now called glam_velo_solver
!
!TODO - Are all these includes needed?
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

!TODO - What is this?
#define shapedbg(x) write(*,*) "x", shape(x)

module glissade_velo

    ! Glissade higher-order velocity solver

    use glam_strs2, only: glam_velo_solver, JFNK_velo_solver

    !globals
    use glimmer_global, only : dp
    use glimmer_physcon, only: gn

    !Other modules that this needs to call out to
    use glide_types
    use glide_grids, only: stagvarb
    use glide_mask
    use glide_grids

    implicit none
    
contains
        
!TODO - Pass explicit arguments instead of model derived type?
    subroutine glissade_velo_driver(model)

        ! Glissade higher-order velocity driver

        use glide_mask
        use stress_hom, only : glide_calcstrsstr

!TODO - Don't think we need glide_thckmask
        use glide_thckmask
        
        type(glide_global_type),intent(inout) :: model

        !For HO masking
        logical :: empty
        integer :: totpts
        real(sp), dimension(model%general%ewn-1, model%general%nsn-1) :: stagmassb

        integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag
        real(dp), dimension(model%general%ewn-1, model%general%nsn-1) :: latbc_norms_stag

        !-------------------------------------------------------------------
        ! Velocity prep that is independent of the solver
        !-------------------------------------------------------------------

!TODO - Would this be a good place to compute geometry derivatives?

!TODO - Verify that glide_set_mask works correctly when the input field is on the velo grid.
!       Would be safer to call a set_mask_staggered subroutine?

        !Compute the "geometry mask" (type of square) for the staggered grid

        call glide_set_mask(model%numerics,                                     &
                            model%geomderv%stagthck, model%geomderv%stagtopg,   &
                            model%general%ewn-1,     model%general%nsn-1,       &
                            model%climate%eus,       geom_mask_stag) 

!TODO - What exactly does this do?  Is it solver-dependent?
        !Augment masks with kinematic boundary condition info
        call augment_kinbc_mask(model%geometry%thkmask, model%velocity%kinbcmask)
        call augment_kinbc_mask(geom_mask_stag, model%velocity%kinbcmask)

!TODO - Remove this call?  Don't think it is ever used.
        !Compute the normal vectors to the marine margin for the staggered grid
!!        call glide_marine_margin_normal(model%geomderv%stagthck, geom_mask_stag, latbc_norms_stag)

        ! save the final mask to 'dynbcmask' for exporting to netCDF output file
        model%velocity%dynbcmask = geom_mask_stag

!TODO - HO_DIAG_PP is the only supported option for now.
!       Later we will add an option for the variational dycore
 
        !-------------------------------------------------------------------
        ! Compute the velocity field
        !-------------------------------------------------------------------

        if (model%options%which_ho_diagnostic == HO_DIAG_PP) then

           if ( model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

!       Are all these options still supported?  Probably can remove periodic_ew/ns.

             call t_startf('glam_velo_solver')
              call glam_velo_solver( model%general%ewn,       model%general%nsn,                 &
                                     model%general%upn,                                          &
                                     model%numerics%dew,      model%numerics%dns,                &
                                     model%numerics%sigma,    model%numerics%stagsigma,          &
                                     model%geometry%thck,     model%geometry%usrf,               &
                                     model%geometry%lsrf,     model%geometry%topg,               &
                                     model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                     model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                     model%geomderv%dlsrfdew, model%geomderv%dlsrfdns,           & 
                                     model%geomderv%stagthck, model%temper%flwa,                 &
                                     model%basalproc%minTauf,                                    & 
                                     model%velocity%btraction,                                   & 
                                     geom_mask_stag,                                             &
                                     model%options%which_ho_babc,                                &
                                     model%options%which_ho_efvs,                                &
                                     model%options%which_ho_resid,                               &
                                     model%options%which_ho_nonlinear,                           &
                                     model%options%which_ho_sparse,                              &
                                     model%options%periodic_ew,                                  &
                                     model%options%periodic_ns,                                  &
                                     model%velocity%beta,                                        & 
                                     model%velocity%uvel, model%velocity%vvel,                   &
                                     model%velocity%uflx, model%velocity%vflx,                   &
                                     model%stress%efvs )
             call t_stopf('glam_velo_solver')

           else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK (solver in development...)

!TODO - Create a JFNK solver that can work with an arbitrary calcF routine
!       (e.g., variational as well as Payne-Price)
! noxsolve could eventually go here 

             call t_startf('JFNK_velo_solver')
              call JFNK_velo_solver (model, geom_mask_stag) 
             call t_stopf('JFNK_velo_solver')

           else
              call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
           end if

        end if   ! which_ho_diagnostic

        !-------------------------------------------------------------------
        ! Velocity-related computations that are independent of the solver
        !-------------------------------------------------------------------

        ! compute the velocity norm

!TODO - Since velnorm is strictly diagnostic, it probably could be computed only for I/O.
        model%velocity%velnorm = sqrt(model%velocity%uvel**2 + model%velocity%vvel**2)
        model%velocity%is_velocity_valid = .true.
        
    end subroutine glissade_velo_driver

end module glissade_velo
