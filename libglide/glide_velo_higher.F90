
#ifdef HAVE_CONFIG_H 
#include "config.inc" 
#endif
#include "glide_nan.inc"
#include "glide_mask.inc"

#define shapedbg(x) write(*,*) "x", shape(x)

module glide_velo_higher
    !Includes for the higher-order velocity computations that this calls out to
    use glam_strs2, only: glam_velo_fordsiapstr, JFNK, umask

    !globals
    use glimmer_global, only : dp
    use glimmer_paramets, only : vis0, vis0_glam 
    use glimmer_physcon, only: gn

    !Other modules that this needs to call out to
    use glide_types
    use glide_grids, only: stagvarb
    use glide_mask
    use glide_grids
    implicit none
    
contains
        
    !This is a temporary wrapper function to get all the HO setup and calling
    !code out of glide_thck so it keeps the clutter down.  I'm just passing the
    !model right now because the analysis on what needs to be passed has yet to
    !be done.
    subroutine run_ho_diagnostic(model)
        use glide_thckmask
        use glide_mask
        
        type(glide_global_type),intent(inout) :: model
        !For HO masking
        logical :: empty
        integer :: totpts
        real(sp), dimension(model%general%ewn-1, model%general%nsn-1) :: stagmassb

        !TEMPORARY arrays, these should at some point be placed in Model probably
        integer, dimension(model%general%ewn-1, model%general%nsn-1)  :: geom_mask_stag
        real(dp), dimension(model%general%ewn-1, model%general%nsn-1) :: latbc_norms_stag

        !Compute the "geometry mask" (type of square) for the staggered grid
        call glide_set_mask(model%numerics, model%geomderv%stagthck, model%geomderv%stagtopg, &
                                model%general%ewn-1, model%general%nsn-1, model%climate%eus, &
                                geom_mask_stag) 

        !Augment masks with kinematic boundary condition info
        call augment_kinbc_mask(model%geometry%thkmask, model%velocity_hom%kinbcmask)
        call augment_kinbc_mask(geom_mask_stag, model%velocity_hom%kinbcmask)

        !Compute the normal vectors to the marine margin for the staggered grid
        call glide_marine_margin_normal(model%geomderv%stagthck, geom_mask_stag, latbc_norms_stag)

        ! save the final mask to 'dynbcmask' for exporting to netCDF output file
        model%velocity_hom%dynbcmask = geom_mask_stag


        if (model%options%which_ho_diagnostic == HO_DIAG_PP) then

          if ( model%options%which_ho_nonlinear == HO_NONLIN_PICARD ) then ! Picard (standard solver)

            call glam_velo_fordsiapstr( model%general%ewn,       model%general%nsn,                 &
                                        model%general%upn,                                          &
                                        model%numerics%dew,      model%numerics%dns,                &
                                        model%numerics%sigma,    model%numerics%stagsigma,          &
                                        model%geometry%thck,     model%geometry%usrf,               &
                                        model%geometry%lsrf,     model%geometry%topg,               &
                                        model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                        model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                        model%geomderv%dlsrfdew, model%geomderv%dlsrfdns,           & 
                                        model%geomderv%stagthck, model%temper%flwa*vis0/vis0_glam,  &
                                        model%basalproc%minTauf,                                    & 
                                        model%velocity_hom%btraction,                               & 
                                        geom_mask_stag,                                             &
                                        model%options%which_ho_babc,                                &
                                        model%options%which_ho_efvs,                                &
                                        model%options%which_ho_resid,                               &
                                        model%options%which_ho_nonlinear,                           &
                                        model%options%which_ho_sparse,                              &
                                        model%options%periodic_ew,                                  &
                                        model%options%periodic_ns,                                  &
                                        model%velocity_hom%beta,                                    & 
                                        model%velocity_hom%uvel, model%velocity_hom%vvel,           &
                                        model%velocity_hom%uflx, model%velocity_hom%vflx,           &
                                        model%velocity_hom%efvs )

          else if ( model%options%which_ho_nonlinear == HO_NONLIN_JFNK ) then ! JFNK (solver in development...)

            call JFNK                  ( model%general%ewn,       model%general%nsn,                 &
                                        model%general%upn,                                          &
                                        model%numerics%dew,      model%numerics%dns,                &
                                        model%numerics%sigma,    model%numerics%stagsigma,          &
                                        model%geometry%thck,     model%geometry%usrf,               &
                                        model%geometry%lsrf,     model%geometry%topg,               &
                                        model%geomderv%dthckdew, model%geomderv%dthckdns,           &
                                        model%geomderv%dusrfdew, model%geomderv%dusrfdns,           &
                                        model%geomderv%dlsrfdew, model%geomderv%dlsrfdns,           & 
                                        model%geomderv%stagthck, model%temper%flwa*vis0/vis0_glam,  &
                                        model%basalproc%minTauf,                                    & 
                                        model%velocity_hom%btraction,                               & 
                                        geom_mask_stag,                                             &
                                        model%options%which_ho_babc,                                &
                                        model%options%which_ho_efvs,                                &
                                        model%options%which_ho_resid,                               &
                                        model%options%which_ho_nonlinear,                           &
                                        model%options%which_ho_sparse,                              &
                                        model%options%periodic_ew,                                  &
                                        model%options%periodic_ns,                                  &
                                        model%velocity_hom%beta,                                    & 
                                        model%velocity_hom%uvel, model%velocity_hom%vvel,           &
                                        model%velocity_hom%uflx, model%velocity_hom%vflx,           &
                                        model%velocity_hom%efvs )
           else
              call write_log('Invalid which_ho_nonlinear option.',GM_FATAL)
           end if

        end if
        !Compute the velocity norm - this is independant of the methods used to compute the u and v components so
        !we put it out here
        model%velocity_hom%velnorm = sqrt(model%velocity_hom%uvel**2 + model%velocity_hom%vvel**2)
        model%velocity_hom%is_velocity_valid = .true.
        
    end subroutine

end module glide_velo_higher
