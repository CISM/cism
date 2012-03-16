

!***********************************************************************
module glam         
!***********************************************************************

    ! 1st-order ice sheet dynamics from Payne/Price OR Pattyn/Bocek/Johonson solver 
    ! (in config file, "diagnostic_scheme" = 3 (PP, B-grid) or 1 (PB&J, A grid) or 2 (PB&J, B grid)
    ! and thickness evolution using LANL incremental remapping (see "remap_advection.F90" for 
    ! documentation)

    use parallel
    use glide_types
    use glimmer_paramets, only : vis0, vis0_glam
    use glimmer_physcon, only :
    use glide_mask

!whl - This is the new transport driver (for upwind or remapping)
    use glissade_transport, only: glissade_transport_driver,  &
                                  nghost_transport, ntracer_transport

!whl - to do - When the new remapping routines (glissade_transport and glissade_remap)
!              are known to be working, we can remove remap_advection and remap_glamutils
    use remap_advection, only: horizontal_remap
    use remap_glamutils

    use glide_velo_higher
    use glide_thck

    implicit none
    private

    public :: inc_remap_driver

    ! NOTE: Relevant initializtion routines are in the init section of "glide.F90" 

!whl - temporary; remove later when the new scheme is the default
    logical, parameter :: old_remapping = .true.  ! if false, then use new remapping scheme
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

!whl - debug
        integer :: i, j
        integer, parameter :: idiag=10, jdiag=15

        ! Compute the new geometry derivatives for this time step

        call geometry_derivs(model)
        call geometry_derivs_unstag(model)

        ! Compute higher-order ice velocities

        if (main_task) then
           print *, ' '
           print *, 'Compute higher-order ice velocities, time =', model%numerics%time
        endif

        call run_ho_diagnostic(model)   ! in glide_velo_higher.F90

        if (main_task) then
           print *, ' '
           print *, 'Compute dH/dt'
        endif

!whl - diagnostics for debugging. Remove later.
        if (write_verbose) then
           write(50,*) ' '
           write(50,*) 'After run_ho_diagnostic:'
           write(50,*) 'Sfc velocity field (model%velocity%uvel):'
!!           do k = 1, model%general%upn
           do k = 1, 1
              write(50,*) ' '
              write(50,*) 'uvel, k =', k
              write(50,*) ' '
              do j = size(model%velocity%uvel,3), 1, -1
                 write(50,'(17e12.4)') model%velocity%uvel(k,1:17,j)
              enddo
              write(50,*) ' '
              write(50,*) 'vvel, k =', k
              write(50,*) ' '
              do j = size(model%velocity%vvel,3), 1, -1
                 write(50,'(17e12.4)') model%velocity%vvel(k,1:17,j)
              enddo
           enddo
        endif

#ifdef JEFFTEST
        ! JEFF Test code to confirm that distributed_gather and distributed_scatter are indeed inverses.
        call distributed_gather_var(model%velocity_hom%efvs, gathered_efvs)
        call distributed_scatter_var(model%velocity_hom%efvs, gathered_efvs, .false.)
        call distributed_gather_var(model%velocity_hom%efvs, gathered_efvs2)

        ! ANY True if any value is true (LOGICAL)
        if (main_task) then
            if (ANY((gathered_efvs(:,:,:) - gathered_efvs2(:,:,:)) /= 0.0)) then
               write(*,*) "Something isn't right.  Gather/Scatter are not inverses."
               call parallel_stop(__FILE__, __LINE__)
            endif
        endif

        deallocate(gathered_efvs)
        deallocate(gathered_efvs2)
#endif

        call parallel_barrier  ! Other procs hang out here waiting for test to complete on main_task.

!whl - When the transport is fully parallel, these gathers will not be needed.
!      But we do need to update the thickness, tracers, and velocity in halo cells.
 
        ! JEFF Glue Code to serialize
        ! Glue code to gather the distributed variables back to main_task processor.
        ! These are outputs from run_ho_diagnostic and are gathered presuming they will be used
        call distributed_gather_var(model%stress%efvs, gathered_efvs)
        call distributed_gather_var(model%velocity%uvel, gathered_uvel)
        call distributed_gather_var(model%velocity%vvel, gathered_vvel)
        call distributed_gather_var(model%velocity%uflx, gathered_uflx)
        call distributed_gather_var(model%velocity%vflx, gathered_vflx)
        call distributed_gather_var(model%velocity%velnorm, gathered_velnorm)

        !Verifying that distributed_gather is working.
        ! call distributed_print("uvel_distributed", model%velocity_hom%uvel)
        ! call parallel_print("uvel_gathered", gathered_uvel)
        !After gathering, then update nsn and ewn to full values (and zero halos?)
        model%general%ewn = global_ewn
        model%general%nsn = global_nsn

        !Gather vars required for following remap routines.
        call distributed_gather_var(model%geometry%thck, gathered_thck)
        call distributed_gather_var(model%geomderv%stagthck, gathered_stagthck)
        call distributed_gather_var(model%climate%acab, gathered_acab)
        call distributed_gather_var(model%temper%temp, gathered_temp, &
                                    lbound(model%temper%temp,1), ubound(model%temper%temp,1))


!whl - Note to Jeff - Introduced a choice here between old and new remapping schemes.
!      The old scheme requires gathering data to the main processor.
!      The new scheme is better designed for distributed parallelism.
!      Old remapping scheme to be removed after the new scheme has been implemented in parallel and tested.

      if( model%numerics%tend > model%numerics%tstart) then
        if (old_remapping) then

           if (main_task) then

            ! Use older incremental remapping algorithm to evolve the ice thickness
            ! (and optionally, temperature and other tracers)

            if (model%options%whichtemp == TEMP_REMAP_ADV) then   ! Use IR to advect temperature
                                                                  ! (and other tracers, if present)

               ! Put relevant model variables into a format that inc. remapping code wants.
               ! (This subroutine lives in module remap_glamutils)
               ! Assume that the remapping grid has horizontal dimensions (ewn-1, nsn-1), so that
               !  scalar and velo grids are the same size.
               ! Also assume that temperature points are staggered in the vertical relative
               !  to velocity points, with temp(0) at the top surface and temp(upn) at the
               !  bottom surface.
               ! Do not advect temp(0), since fixed at artm
               ! At least for now, do not advect temp(upn) either.

#ifdef JEFFORIG
               call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                         model%geometry%thck(1:model%general%ewn-1,1:model%general%nsn-1),  &
                                         model%velocity_hom%uflx, model%velocity_hom%vflx,               &
                                         model%geomderv%stagthck, model%numerics%thklim,                 &
                                         model%velocity_hom%uvel, model%velocity_hom%vvel,               &
                                         model%temper%temp  (1:model%general%upn-1,                        &
                                                             1:model%general%ewn-1,1:model%general%nsn-1))
#endif
               call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                         gathered_thck(1:model%general%ewn-1,1:model%general%nsn-1),  &
                                         gathered_uflx, gathered_vflx,               &
                                         gathered_stagthck, model%numerics%thklim,                 &
                                         model%options%periodic_ew,             model%options%periodic_ns, &
                                         gathered_uvel, gathered_vvel,               &
                                         gathered_temp  (1:model%general%upn-1,                        &
                                                         1:model%general%ewn-1,1:model%general%nsn-1))

               ! Remap temperature and fractional thickness for each layer
               ! (This subroutine lives in module remap_advection)

               model%remap_wk%dt_ir = model%remap_wk%dt_ir / model%numerics%subcyc
               do sc = 1 , model%numerics%subcyc
                  write(*,*) 'Processing subcycling step: ',sc
                  do k = 1, model%general%upn-1
                     !JEFF all variables going into horizontal_remap are relative to remap_wk, so no gathering required.
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

               ! Interpolate tracers back to sigma coordinates
               !JEFF sigma is "Sigma values for vertical spacing of model levels" which is the same on all nodes.
               !JEFF Rest of parameters vertical_remap are relative to remap_wk, so no gathering required.
               call vertical_remap( model%general%ewn-1,     model%general%nsn-1,               &
                                    model%general%upn,       ntrace_ir,                         &
                                    model%numerics%sigma,    model%remap_wk%thck_ir(:,:,:),     &
                                    model%remap_wk%trace_ir)

#ifdef JEFFORIG
               ! put output from inc. remapping code back into format that model wants
               ! (this subroutine lives in module remap_glamutils)
               call horizontal_remap_out (model%remap_wk,     model%geometry%thck,            &
                                          model%climate%acab, model%numerics%dt,              &
                                          model%temper%temp(1:model%general%upn-1,:,:) )
#endif
                !JEFF gathered_thck is updated in this procedure.  I don't know about gathered_temp yet.
                ! put output from inc. remapping code back into format that model wants
                ! (this subroutine lives in "remap_glamutils.F90")
                call horizontal_remap_out(model%remap_wk, gathered_thck,                 &
                                          gathered_acab, model%numerics%dt,              &
                                          gathered_temp(1:model%general%upn-1,:,:) )

            else  ! Use IR to transport thickness only
                ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)
                ! (this subroutine lives in module remap_advection)

                ! put relevant model variables into a format that inc. remapping code wants
                ! (this subroutine lives in "remap_glamutils.F90")
#ifdef JEFFORIG
               Updated call call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                         model%geometry%thck(1:model%general%ewn-1,1:model%general%nsn-1),    &
                                         model%velocity_hom%uflx, model%velocity_hom%vflx,               &
                                         model%geomderv%stagthck, model%numerics%thklim)
#endif
               call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                         gathered_thck(1:model%general%ewn-1,1:model%general%nsn-1),    &
                                         gathered_uflx, gathered_vflx,               &
                                         gathered_stagthck, model%numerics%thklim,   &
                                         model%options%periodic_ew,             model%options%periodic_ns)

                !JEFF all variables going into horizontal_remap are relative to remap_wk, so no gathering required.
                ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)
                ! (this subroutine lives in module remap_advection)
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

#ifdef JEFFORIG
                ! put output from inc. remapping code back into format that model wants
                ! (this subroutine lives in module remap_glamutils.F90)
               call horizontal_remap_out (model%remap_wk,     model%geometry%thck,                 &
                                          model%climate%acab, model%numerics%dt )
#endif
                ! gathered_thck is updated in this procedure
                ! put output from inc. remapping code back into format that model wants
                ! (this subroutine lives in "remap_glamutils.F90")
                call horizontal_remap_out(model%remap_wk, gathered_thck,                 &
                                           gathered_acab, model%numerics%dt )

!whl - debug - Gather temperature for optional diagnostics.
                if (write_verbose) then
                   write(50,*) 'Gather temperature'
                   call distributed_gather_var(model%temper%temp, gathered_temp,  &
                                               lbound(model%temper%temp,1), ubound(model%temper%temp,1))
                endif

            endif   ! whichtemp

           endif    ! main_task

           call parallel_barrier   ! Other tasks hold here until main_task completes

!whl - Note to Jeff - Here is the code used to call the new remapping scheme
!      (To be modified for distributed parallelism)

        else   ! new remapping scheme

!whl - Note to Jeff - Need halo updates here for thck, temp (and any other advected tracers), uvel and vvel.
!                     If nhalo >= 2, then no halo updates should be needed inside glissade_transport_driver.

           ! Halo updates for thickness and tracers.

!           call ice_HaloUpdate (thck,             halo_info,     &
!                                field_loc_center, field_type_scalar)

!           call ice_HaloUpdate (temp,             halo_info,     &
!                                field_loc_center, field_type_scalar)

!           call ice_HaloUpdate (uvel,               halo_info,     &
!                                field_loc_NEcorner, field_type_vector)

!           call ice_HaloUpdate (vvel,               halo_info,     &
!                                field_loc_NEcorner, field_type_vector)

!whl - The new remapping routines assume that the thickness array has horizontal dimensions 
!      equal to model%general%ewn and nsn (as defined at initialization).
!      Above, Jeff redefined model%general%ewn = global_ewn, model%general%nsn = global_nsn.
!      This is *not* correct for the new remapping, so I redefined the dimensions here.
!      When the code is fully parallel, the following two lines should be unnecessary.

           model%general%ewn = size(model%geometry%thck,1)  !whl - just for the call to glissade_transport_driver
           model%general%nsn = size(model%geometry%thck,2)

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

!whl - Gather new thck and temp.
           write(50,*) 'Gather thickness and temperature'
           call distributed_gather_var(model%geometry%thck, gathered_thck)
           call distributed_gather_var(model%temper%temp, gathered_temp,  &
                                       lbound(model%temper%temp,1), ubound(model%temper%temp,1))

!whl -     Reset model%general%ewn and nsn to the global values.  
!          These were temporarily set to the local values before the new remapping was called.
!          The following two lines should be removed when the code is fully distributed.

           model%general%ewn = global_ewn
           model%general%nsn = global_nsn

!whl - Note to Jeff - Think about whether we should do any halo updates after the transport is done.

        endif  ! old v. new remapping

!whl - optional diagnostics - remove later
!whl - write gathered_thck and temp, since the old remapping scheme does not update the local arrays

        if (write_verbose) then
           write(50,*) ' '
           write(50,*) 'After transport:'
           write(50,*) ' '
           write(50,*) 'Gathered thickness field (gathered_thck):'
           write(50,*) ' '
           do j = size(gathered_thck,2), 1, -1
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

        ! NOTE finalization routine, to be written for PP HO core needs to be written (e.g.
        ! glam_velo_fordsiapstr_final( ) ), added to glam_strs2.F90, and called from glide_stop.F90

!***********************************************************************
end module glam
!***********************************************************************
