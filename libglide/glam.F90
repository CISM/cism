
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

    use remap_advection
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

        integer :: ewn, nsn

        integer ::         & 
          ntrace_ir       ,&! number of tracers to be remapped
          nghost_ir         ! number of ghost cells used in remapping scheme

        ! This argument is not already included in the model derived type, so add it here
        ! (needs to be added at some point in the future)
        real (kind=dp), dimension(model%general%ewn-1,model%general%nsn-1) :: minTauf
        minTauf = 0.0d0

        ! Compute the new geometry derivatives for this time step
        call geometry_derivs(model)
        call geometry_derivs_unstag(model)  

        print *, ' '
        print *, 'time = ', model%numerics%time

        ! This driver is called from "glide_velo_higher.F90"
        call run_ho_diagnostic(model)

        ! put relevant model variables into a format that inc. remapping code wants
        ! (this subroutine lives in "remap_glamutils.F90")
        call horizontal_remap_in(model%remap_wk, model%numerics%dt, model%geometry%thck(1:ewn-1,1:nsn-1),  &
                                  ntrace_ir,               nghost_ir,                             &
                                  model%numerics%dew,      model%numerics%dns,                    &
                                  model%velocity_hom%uflx, model%velocity_hom%vflx,               &
                                  model%geomderv%stagthck, model%numerics%thklim,                 &
                                  model%options%periodic_ew,             model%options%periodic_ns)

        ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)
        ! (this subroutine lives in "remap_advection.F90")
        call horizontal_remap( model%remap_wk%dt_ir,                                         & 
                               model%remap_wk%ewn_ir,      model%remap_wk%nsn_ir,            &
                               ntrace_ir,   nghost_ir,                                       &
                               model%remap_wk%ubar_ir,     model%remap_wk%vbar_ir,           &
                               model%remap_wk%thck_ir,     model%remap_wk%trace_ir,          &
                               model%remap_wk%dew_ir,      model%remap_wk%dns_ir,            &
                               model%remap_wk%dewt_ir,     model%remap_wk%dnst_ir,           &
                               model%remap_wk%dewu_ir,     model%remap_wk%dnsu_ir,           &
                               model%remap_wk%hm_ir,       model%remap_wk%tarea_ir)

        ! put output from inc. remapping code back into format that model wants
        ! (this subroutine lives in "remap_glamutils.F90")
        call horizontal_remap_out (model%remap_wk, model%geometry%thck,                 &
                                   model%climate%acab, model%numerics%dt,               &
                                   model%options%periodic_ew, model%options%periodic_ns)


    end subroutine inc_remap_driver 


!***********************************************************************
end module glam
!***********************************************************************
