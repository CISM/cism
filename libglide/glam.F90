!CLEANUP - glam.F90
! This module is now decremented.
! For now, it is still possible to call subroutine inc_remap_driver
!  from glissade_tstep by setting call_inc_remap_driver = .true.
! But the default option is to call the velocity solver and transport
!  routines from glissade_tstep, where I inlined code from glam.F90.
!  I inlined the new remapping scheme only.

! The TODOs in this module can be ignored.
!***********************************************************************
module glam         
!***********************************************************************

!TODO - Remove PBJ references in this module
    ! 1st-order ice sheet dynamics from Payne/Price OR Pattyn/Bocek/Johonson solver 
    ! (in config file, "diagnostic_scheme" = 3 (PP, B-grid) or 1 (PB&J, A grid) or 2 (PB&J, B grid)
    ! and thickness evolution using LANL incremental remapping (see "remap_advection.F90" for 
    ! documentation)

    use parallel
    use glide_types
!TODO - Remove scaling (vis0 and vis0_glam not used in this module anyway)
    use glimmer_paramets, only : vis0, vis0_glam
    use glide_mask

!TODO - Remove this one?
    use glimmer_physcon, only :

    use glide_velo_higher

    use glissade_transport, only: glissade_transport_driver,  &
                                  nghost_transport, ntracer_transport

!TODO - When the new remapping routines (glissade_transport and glissade_remap)
!        are known to be working, we can remove remap_advection and remap_glamutils
    use remap_advection, only: horizontal_remap
    use remap_glamutils


!TODO - May not need to use glide_thck
    use glide_thck

    implicit none
    private

!TODO - Remove old_remapping
    public :: inc_remap_driver, old_remapping

    ! NOTE: Relevant initialization routines are in the init section of "glide.F90" 

!TODO - temporary; remove later when the new scheme is the default
    logical, parameter :: old_remapping = .false.  ! if false, then use new remapping scheme
                                                  ! if true, revert to older remapping scheme
    logical, parameter :: write_verbose = .false. ! if true, write state variable fields to log file

    contains

    ! This driver is called from "glide.F90", in subroutine "glide_tstep_p2"
    subroutine inc_remap_driver( model )
        use parallel
        implicit none

        type(glide_global_type), intent(inout) :: model

        integer :: k, sc    ! layer index

! JCC - Periodic support disabled
!         integer ::         & 
!           ntrace_ir       ,&! number of tracers to be remapped
!           nghost_ir         ! number of ghost cells used in remapping scheme

!TODO - Remove debug code.
!debug
        integer :: i, j
        integer, parameter :: idiag=10, jdiag=15

        ! Compute the new geometry derivatives for this time step

!HALO - Would be better to have one derivative call per field.
!       Also could consider computing derivatives in glam_strs2.F90, so as
!        not to have to pass them through the interface.
!       Do halo updates as needed (e.g., thck) before computing derivatives.
!       If we need a derivative on the staggered grid (for all locally owned cells),
!        then we need one layer of halo scalars before calling the derivative routine.

        call geometry_derivs(model)

!HALO - Pretty sure this is not needed
        call geometry_derivs_unstag(model)

        ! Compute higher-order ice velocities

        if (main_task) then
           print *, ' '
           print *, 'Compute higher-order ice velocities, time =', model%numerics%time
        endif

!TODO - Rename run_ho_diagnostic?
       call t_startf('run_ho_diagnostic')
        call run_ho_diagnostic(model)   ! in glide_velo_higher.F90
       call t_stopf('run_ho_diagnostic')

        if (main_task) then
           print *, ' '
           print *, 'Compute dH/dt'
        endif

!TODO - Remove the old remapping scheme.

!TODO - I'm not entirely clear on which options are to be supported.
!       Will we continue to support serial runs that use remapping for advection while
!        using the old Glimmer temperature scheme?

!whl   Introduced a choice here between old and new remapping schemes.
!      The old scheme requires gathering data to the main processor.
!      The new scheme is better designed for distributed parallelism.
!      Old remapping scheme to be removed after the new scheme has been implemented in parallel and tested.

      if( model%numerics%tend > model%numerics%tstart) then
        if (old_remapping) then

!HALO - These are for serial advection with parallel velocity solver. 
!       Can remove when we switch to fully distributed with parallel advection.

           ! Glue code to gather the distributed variables back to main_task processor.
           ! These are outputs from run_ho_diagnostic and are gathered presuming they will be used
          call t_startf('old_remap_gathers')
           call distributed_gather_var(model%stress%efvs, gathered_efvs)
           call distributed_gather_var(model%velocity%uvel, gathered_uvel)
           call distributed_gather_var(model%velocity%vvel, gathered_vvel)
           call distributed_gather_var(model%velocity%uflx, gathered_uflx)
           call distributed_gather_var(model%velocity%vflx, gathered_vflx)
           call distributed_gather_var(model%velocity%velnorm, gathered_velnorm)

           !Gather vars required for following remap routines.
           call distributed_gather_var(model%geometry%thck, gathered_thck)
           call distributed_gather_var(model%geomderv%stagthck, gathered_stagthck)
           call distributed_gather_var(model%climate%acab, gathered_acab)
           call distributed_gather_var(model%temper%temp, gathered_temp, &
                                       lbound(model%temper%temp,1), ubound(model%temper%temp,1))

           !After gathering, then update nsn and ewn to full values (and zero halos?)
           model%general%ewn = global_ewn
           model%general%nsn = global_nsn
          call t_stopf('old_remap_gathers')

           if (main_task) then

            ! Use older incremental remapping algorithm to evolve the ice thickness
            ! (and optionally, temperature and other tracers)

            if (model%options%whichtemp == TEMP_REMAP_ADV) then   ! Use IR to advect temperature
                                                                  ! (and other tracers, if present)

               ! Put relevant model variables into a format that inc. remapping code wants.
               ! Assume that the remapping grid has horizontal dimensions (ewn-1, nsn-1), so that
               !  scalar and velo grids are the same size.
               ! Also assume that temperature points are staggered in the vertical relative
               !  to velocity points, with temp(0) at the top surface and temp(upn) at the
               !  bottom surface.
               ! Do not advect temp(0), since fixed at artm
               ! At least for now, do not advect temp(upn) either.

              call t_startf('horizontal_remap_in')
               call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                         gathered_thck(1:model%general%ewn-1,1:model%general%nsn-1),  &
                                         gathered_uflx, gathered_vflx,               &
                                         gathered_stagthck, model%numerics%thklim,                 &
                                         model%options%periodic_ew,             model%options%periodic_ns, &
                                         gathered_uvel, gathered_vvel,               &
                                         gathered_temp  (1:model%general%upn-1,                        &
                                                         1:model%general%ewn-1,1:model%general%nsn-1))
              call t_stopf('horizontal_remap_in')

               ! Remap temperature and fractional thickness for each layer

              call t_startf('horizontal_remap')
               model%remap_wk%dt_ir = model%remap_wk%dt_ir / model%numerics%subcyc
               do sc = 1 , model%numerics%subcyc
                  write(*,*) 'Processing subcycling step: ',sc
                  do k = 1, model%general%upn-1
                     !all variables going into horizontal_remap are relative to remap_wk, so no gathering required.
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
               enddo
               model%remap_wk%dt_ir = model%remap_wk%dt_ir * model%numerics%subcyc
              call t_stopf('horizontal_remap')

               ! Interpolate tracers back to sigma coordinates
               ! sigma is "Sigma values for vertical spacing of model levels" which is the same on all nodes.
               ! Rest of parameters vertical_remap are relative to remap_wk, so no gathering required.
              call t_startf('vertical_remap')
               call vertical_remap( model%general%ewn-1,     model%general%nsn-1,               &
                                    model%general%upn,       ntrace_ir,                         &
                                    model%numerics%sigma,    model%remap_wk%thck_ir(:,:,:),     &
                                    model%remap_wk%trace_ir)
              call t_stopf('vertical_remap')

                ! gathered_thck and gathered_temp updated in this procedure.  
                ! put output from inc. remapping code back into format that model wants
              call t_startf('horizontal_remap_out')
               call horizontal_remap_out(model%remap_wk, gathered_thck,                 &
                                         gathered_acab, model%numerics%dt,              &
                                         gathered_temp(1:model%general%upn-1,:,:) )
              call t_stopf('horizontal_remap_out')

            else  ! Use IR to transport thickness only
                ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)

                ! put relevant model variables into a format that inc. remapping code wants

              call t_startf('horizontal_remap_in')
               call horizontal_remap_in (model%remap_wk, model%numerics%dt,                             &
                                         gathered_thck(1:model%general%ewn-1,1:model%general%nsn-1),    &
                                         gathered_uflx, gathered_vflx,               &
                                         gathered_stagthck, model%numerics%thklim,   &
                                         model%options%periodic_ew, model%options%periodic_ns)
              call t_stopf('horizontal_remap_in')

                ! All variables going into horizontal_remap are relative to remap_wk, so no gathering required.
                ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)
              call t_startf('horizontal_remap')
                model%remap_wk%dt_ir = model%remap_wk%dt_ir / model%numerics%subcyc
                do sc = 1 , model%numerics%subcyc
                   write(*,*) 'Processing subcycling step: ',sc
                   call horizontal_remap( model%remap_wk%dt_ir,                                         &
                                          model%general%ewn-1,        model%general%nsn-1,              &
                                          ntrace_ir,                  nghost_ir,                        &
                                          model%remap_wk%uvel_ir,     model%remap_wk%vvel_ir,           &
                                          model%remap_wk%thck_ir,     model%remap_wk%trace_ir,          &
                                          model%remap_wk%dew_ir,      model%remap_wk%dns_ir,            &
                                          model%remap_wk%dewt_ir,     model%remap_wk%dnst_ir,           &
                                          model%remap_wk%dewu_ir,     model%remap_wk%dnsu_ir,           &
                                          model%remap_wk%hm_ir,       model%remap_wk%tarear_ir)
                enddo
                model%remap_wk%dt_ir = model%remap_wk%dt_ir * model%numerics%subcyc
              call t_stopf('horizontal_remap')

                ! gathered_thck is updated in this procedure
                ! put output from inc. remapping code back into format that model wants
              call t_startf('horizontal_remap_out')
               call horizontal_remap_out( model%remap_wk, gathered_thck,   &
                                          gathered_acab, model%numerics%dt )
              call t_stopf('horizontal_remap_out')

            endif   ! whichtemp

           endif    ! main_task

          call t_startf('old_remap_scatters')
           !scatter thck (and temp) back to other processes
           call distributed_scatter_var(model%geometry%thck, gathered_thck)
           if (model%options%whichtemp == TEMP_REMAP_ADV) then
             call distributed_scatter_var(model%temper%temp, gathered_temp)
             !If advecting other tracers, add scatter here
           endif
          call t_stopf('old_remap_scatters')

           !After scattering, reset nsn and ewn to distributed values
           model%general%ewn = local_ewn
           model%general%nsn = local_nsn

        else   ! new remapping scheme

!HALO - Need halo updates here for thck, temp (and any other advected tracers), uvel and vvel.
!       If nhalo >= 2, then no halo updates should be needed inside glissade_transport_driver.

!PW FOLLOWING NECESSARY?
!HALO - These halo updates could be moved up a level to the new glissade driver.

           ! Halo updates for velocities, thickness and tracers
          call t_startf('new_remap_halo_upds')
           call staggered_parallel_halo(model%velocity%uvel)
           call staggered_parallel_halo(model%velocity%vvel)
           call parallel_halo(model%geometry%thck)
           if (model%options%whichtemp == TEMP_REMAP_ADV) then
              !If advecting other tracers, add parallel_halo update here
              call parallel_halo(model%temper%temp)
           endif
          call t_stopf('new_remap_halo_upds')

          call t_startf('glissade_transport_driver')
          model%numerics%dt = model%numerics%dt / model%numerics%subcyc
          do sc = 1 , model%numerics%subcyc
           if (model%options%whichtemp == TEMP_REMAP_ADV) then  ! Use IR to transport thickness, temperature
                                                                ! (and other tracers, if present)

              call glissade_transport_driver(model%numerics%dt * tim0,                             &  
                                             model%numerics%dew * len0, model%numerics%dns * len0, &
                                             model%general%ewn,         model%general%nsn,         &
                                             model%general%upn-1,       model%numerics%sigma,      &
                                             nghost_transport,          ntracer_transport,         &
                                             model%velocity%uvel(:,:,:) * vel0,                    &
                                             model%velocity%vvel(:,:,:) * vel0,                    &
                                             model%geometry%thck(:,:),                             &
                                             model%temper%temp(1:model%general%upn-1,:,:) )

!TODO - Will we continue to support this option?

           else  ! Use IR to transport thickness only
                 ! Note: In glissade_transport_driver, the ice thickness is transported layer by layer,
                 !       which is inefficient if no tracers are being transported.  (It would be more
                 !       efficient to transport thickness in one layer only, using a vertically
                 !       averaged velocity.)  But this option probably will not be used in practice;
                 !       it is left in the code just to ensure backward compatibility with an
                 !       older remapping scheme for transporting thickness only.

              call glissade_transport_driver(model%numerics%dt * tim0,                             &  
                                             model%numerics%dew * len0, model%numerics%dns * len0, &
                                             model%general%ewn,         model%general%nsn,         &
                                             model%general%upn-1,       model%numerics%sigma,      &
                                             nghost_transport,          1,                         &           
                                             model%velocity%uvel(:,:,:) * vel0,                &
                                             model%velocity%vvel(:,:,:) * vel0,                &
                                             model%geometry%thck(:,:))

           endif  ! whichtemp
          enddo
          model%numerics%dt = model%numerics%dt * model%numerics%subcyc
          call t_stopf('glissade_transport_driver')

!TODO - Remove these gathers later
           if (write_verbose) then
              call distributed_gather_var(model%geometry%thck, gathered_thck)
              call distributed_gather_var(model%temper%temp, gathered_temp, &
                                          lbound(model%temper%temp,1), ubound(model%temper%temp,1))
           endif

        endif  ! old v. new remapping

!TODO - Would this be the appropriate place to add/remove ice at the upper and lower surfaces?  
!       Should probably do this before vertical remapping.
!       (Vertical remapping currently lives in glissade_transport_driver)
!       Note that vertical remapping is needed to return to standard sigma levels,
!        as assumed by both the temperature and velocity solvers.

        !Update halos of modified fields
       call t_startf('after_remap_haloupds')

!HALO - Move these updates to the new glissade driver.

        call parallel_halo(model%geometry%thck)

        if (model%options%whichtemp /= TEMP_GLIMMER) then   ! Glimmer temperature arrays have 
                                                            ! extra horizontal dimensions; not supported
                                                            ! for halo updates
           call parallel_halo(model%temper%temp)
        endif
       call t_stopf('after_remap_haloupds')

!TODO - optional diagnostics - remove later
!whl - write gathered_thck and temp

        if ((write_verbose) .and. (main_task)) then
           write(50,*) ' '
           write(50,*) 'After transport:'
           write(50,*) ' '
           write(50,*) 'Gathered thickness field (gathered_thck):'
           write(50,*) ' '
           do j = size(gathered_thck,2), 1, -1
              write(50,*) 'J=',j
              write(50,'(16e12.4)') gathered_thck(1:16,j)
           enddo
           write(50,*) ' '
           write(50,*) 'Gathered temperature field (gathered_temp):'
!!          do k = 0, model%general%upn
           k = model%general%upn-1
              write(50,*) ' '
              write(50,*) 'temp, k =', k
              write(50,*) ' '
              do j = size(gathered_temp,3), 1, -1
                 write(50,*) 'J=',j
                 write(50,'(16e12.4)') gathered_temp(k,1:16,j)
              enddo
!!          enddo

           i = idiag
           j = jdiag
           write(50,*) ' '
           write(50,*), 'gathered_thck:', i-2, j-2, gathered_thck(i-2,j-2)
           write(50,*) ' '
           write(50,*), 'column temps:'
           do k = 1, model%general%upn
             write(50,*) k, gathered_temp(k,i-2,j-2)
           enddo
         endif   ! write_verbose

      endif     ! if tend > tstart

    end subroutine inc_remap_driver

!TODO - Has this been done?
        ! NOTE finalization routine, to be written for PP HO core needs to be written (e.g.
        ! glam_velo_fordsiapstr_final( ) ), added to glam_strs2.F90, and called from glide_stop.F90

!***********************************************************************
end module glam
!***********************************************************************
