
!***********************************************************************
module glam         
!***********************************************************************

    ! 1st-order ice sheet dynamics from Payne/Price OR Pattyn/Bocek/Johonson solver 
    ! (in config file, "diagnostic_scheme" = 3 (PP, B-grid) or 1 (PB&J, A grid) or 2 (PB&J, B grid)
    ! and thickness evolution using LANL incremental remapping (see "remap_advection.F90" for 
    ! documentation)

    use glide_types
    use glimmer_paramets, only : vis0, vis0_glam
    use glimmer_physcon, only :
    use glide_mask
    use remap_advection, only: horizontal_remap
    use remap_glamutils
    use glide_velo_higher
    use glide_thck

    implicit none
    private

    public :: inc_remap_driver

    ! NOTE: Relevant initializtion routines are in the init section of "glide.F90" 

    contains

    ! This driver is called from "glide.F90", in subroutine "glide_tstep_p2"
    subroutine inc_remap_driver( model )

        type(glide_global_type), intent(inout) :: model

        integer :: k    ! layer index

! JCC - Periodic support disabled
!         integer ::         & 
!           ntrace_ir       ,&! number of tracers to be remapped
!           nghost_ir         ! number of ghost cells used in remapping scheme
        ! Compute the new geometry derivatives for this time step

        call geometry_derivs(model)
        call geometry_derivs_unstag(model)

        ! Compute higher-order ice velocities

        print *, ' '
        print *, '(dH/dt using incremental remapping)'
        print *, 'time = ', model%numerics%time

        call run_ho_diagnostic(model)   ! in glide_velo_higher.F90

        ! Use incremental remapping algorithm to evolve the ice thickness
        ! (and optionally, temperature and other tracers)

        if (model%options%whichtemp == TEMP_REMAP_ADV) then   ! Use IR to advect temperature
                                                              ! (and other tracers, if present)

           ! Put relevant model variables into a format that inc. remapping code wants.
           ! (This subroutine lives in module remap_glamutils)
           ! Assume that the remapping grid has horizontal dimensions (ewn-1, nsn-1), so that
           !  scalar and velo grids are the same size.
           ! Alse assume that temperature points are staggered in the vertical relative
           !  to velocity points, with temp(0) at the top surface and temp(upn) at the
           !  bottom surface. 
           ! Do not advect temp(0), since fixed at artm
           ! At least for now, do not advect temp(upn) either.

           call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                     model%geometry%thck(1:model%general%ewn-1,1:model%general%nsn-1),  &
                                     model%velocity%uflx, model%velocity%vflx,               &
                                     model%geomderv%stagthck, model%numerics%thklim,                 &
                                     model%options%periodic_ew, model%options%periodic_ns,           &
                                     model%velocity%uvel, model%velocity%uvel,               &
                                     model%temper%temp  (1:model%general%upn-1,                      &
                                                         1:model%general%ewn-1,1:model%general%nsn-1) &
                                     )

           ! Remap temperature and fractional thickness for each layer
           ! (This subroutine lives in module remap_advection)

           do k = 1, model%general%upn-1
              
              call horizontal_remap( model%remap_wk%dt_ir,                                         &
                                     model%general%ewn-1,        model%general%nsn-1,              &
                                     ntrace_ir,                  nghost_ir,                        &
                                     model%remap_wk%uvel_ir (:,:,k),                               &
                                     model%remap_wk%vvel_ir (:,:,k),                               &
                                     model%remap_wk%thck_ir (:,:,k),                               &
                                     model%remap_wk%trace_ir(:,:,:,k),                             &
                                     model%remap_wk%dew_ir,      model%remap_wk%dns_ir,            &
                                     model%remap_wk%dewt_ir,     model%remap_wk%dnst_ir,           &
                                     model%remap_wk%dewu_ir,     model%remap_wk%dnsu_ir,           &
                                     model%remap_wk%hm_ir,       model%remap_wk%tarear_ir)

           enddo

           ! Interpolate tracers back to sigma coordinates

           call vertical_remap( model%general%ewn-1,     model%general%nsn-1,               &
                                model%general%upn,       ntrace_ir,                         &
                                model%numerics%sigma,    model%remap_wk%thck_ir(:,:,:),     & 
                                model%remap_wk%trace_ir)

           ! put output from inc. remapping code back into format that model wants
           ! (this subroutine lives in module remap_glamutils)

           call horizontal_remap_out (model%remap_wk,     model%geometry%thck,            &
                                      model%climate%acab, model%numerics%dt,              &
                                      model%temper%temp(1:model%general%upn-1,:,:) )

        else  ! Use IR to transport thickness only

           call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                     model%geometry%thck(1:model%general%ewn-1,1:model%general%nsn-1),    &
                                     model%velocity%uflx, model%velocity%vflx,               &
                                     model%geomderv%stagthck, model%numerics%thklim,                 &
                                     model%options%periodic_ew,             model%options%periodic_ns)

        ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)
        ! (this subroutine lives in module remap_advection)

           call horizontal_remap( model%remap_wk%dt_ir,                                         &
                                  model%general%ewn-1,        model%general%nsn-1,              &
                                  ntrace_ir,                  nghost_ir,                        &
                                  model%remap_wk%uvel_ir,     model%remap_wk%vvel_ir,           &
                                  model%remap_wk%thck_ir,     model%remap_wk%trace_ir,          &
                                  model%remap_wk%dew_ir,      model%remap_wk%dns_ir,            &
                                  model%remap_wk%dewt_ir,     model%remap_wk%dnst_ir,           &
                                  model%remap_wk%dewu_ir,     model%remap_wk%dnsu_ir,           &
                                  model%remap_wk%hm_ir,       model%remap_wk%tarear_ir)

        ! Put output from inc. remapping code back into format that model wants.
        ! (This subroutine lives in module remap_glamutils)

           call horizontal_remap_out (model%remap_wk,     model%geometry%thck,                 &
                                      model%climate%acab, model%numerics%dt )

        endif   ! whichtemp

    end subroutine inc_remap_driver

        ! NOTE finalization routine, to be written for PP HO core needs to be written (e.g.
        ! glam_velo_fordsiapstr_final( ) ), added to glam_strs2.F90, and called from glide_stop.F90

!***********************************************************************
end module glam
!***********************************************************************
