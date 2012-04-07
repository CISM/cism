

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

    public :: inc_remap_driver, old_remapping

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
        if ((write_verbose) .and. (main_task)) then
           write(50,*) ' '
           write(50,*) 'After run_ho_diagnostic:'
           write(50,*) 'Sfc velocity field (model%velocity%uvel):'
!!           do k = 1, model%general%upn
           do k = 1, 1
              write(50,*) ' '
              write(50,*) 'uvel, k =', k
              write(50,*) ' '
              do j = size(model%velocity%uvel,3), 1, -1
              do i = size(model%velocity%uvel,2), 1, -1
                 write(50,'(i4,1x,i4,1x,e12.4)') i, j, model%velocity%uvel(k,i,j)
              enddo
              enddo
              write(50,*) ' '
              write(50,*) 'vvel, k =', k
              write(50,*) ' '
              do j = size(model%velocity%vvel,3), 1, -1
              do i = size(model%velocity%vvel,2), 1, -1
                 write(50,'(i4,1x,i4,1x,e12.4)') i, j, model%velocity%vvel(k,i,j)
              enddo
              enddo
           enddo
        endif

!whl   Introduced a choice here between old and new remapping schemes.
!      The old scheme requires gathering data to the main processor.
!      The new scheme is better designed for distributed parallelism.
!      Old remapping scheme to be removed after the new scheme has been implemented in parallel and tested.

      if( model%numerics%tend > model%numerics%tstart) then
        if (old_remapping) then

           ! Glue code to gather the distributed variables back to main_task processor.
           ! These are outputs from run_ho_diagnostic and are gathered presuming they will be used
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

               call horizontal_remap_in (model%remap_wk,          model%numerics%dt,                     &
                                         gathered_thck(1:model%general%ewn-1,1:model%general%nsn-1),  &
                                         gathered_uflx, gathered_vflx,               &
                                         gathered_stagthck, model%numerics%thklim,                 &
                                         model%options%periodic_ew,             model%options%periodic_ns, &
                                         gathered_uvel, gathered_vvel,               &
                                         gathered_temp  (1:model%general%upn-1,                        &
                                                         1:model%general%ewn-1,1:model%general%nsn-1))

               ! Remap temperature and fractional thickness for each layer

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

               ! Interpolate tracers back to sigma coordinates
               ! sigma is "Sigma values for vertical spacing of model levels" which is the same on all nodes.
               ! Rest of parameters vertical_remap are relative to remap_wk, so no gathering required.
               call vertical_remap( model%general%ewn-1,     model%general%nsn-1,               &
                                    model%general%upn,       ntrace_ir,                         &
                                    model%numerics%sigma,    model%remap_wk%thck_ir(:,:,:),     &
                                    model%remap_wk%trace_ir)

                ! gathered_thck and gathered_temp updated in this procedure.  
                ! put output from inc. remapping code back into format that model wants
                call horizontal_remap_out(model%remap_wk, gathered_thck,                 &
                                          gathered_acab, model%numerics%dt,              &
                                          gathered_temp(1:model%general%upn-1,:,:) )

            else  ! Use IR to transport thickness only
                ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)

                ! put relevant model variables into a format that inc. remapping code wants

               call horizontal_remap_in (model%remap_wk, model%numerics%dt,                             &
                                         gathered_thck(1:model%general%ewn-1,1:model%general%nsn-1),    &
                                         gathered_uflx, gathered_vflx,               &
                                         gathered_stagthck, model%numerics%thklim,   &
                                         model%options%periodic_ew, model%options%periodic_ns)

                ! All variables going into horizontal_remap are relative to remap_wk, so no gathering required.
                ! call inc. remapping code for thickness advection (i.e. dH/dt calcualtion)
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

                ! gathered_thck is updated in this procedure
                ! put output from inc. remapping code back into format that model wants
                call horizontal_remap_out( model%remap_wk, gathered_thck,   &
                                           gathered_acab, model%numerics%dt )

            endif   ! whichtemp

           endif    ! main_task

           !scatter thck (and temp) back to other processes
           call distributed_scatter_var(model%geometry%thck, gathered_thck)
           if (model%options%whichtemp == TEMP_REMAP_ADV) then
             call distributed_scatter_var(model%temper%temp, gathered_temp)
             !If advecting other tracers, add scatter here
           endif

           !After scattering, reset nsn and ewn to distributed values
           model%general%ewn = local_ewn
           model%general%nsn = local_nsn

        else   ! new remapping scheme

!whl Need halo updates here for thck, temp (and any other advected tracers), uvel and vvel.
!    If nhalo >= 2, then no halo updates should be needed inside glissade_transport_driver.

!PW FOLLOWING NECESSARY?
           ! Halo updates for velocities, thickness and tracers
           call staggered_parallel_halo(model%velocity%uvel)
           call staggered_parallel_halo(model%velocity%vvel)
           call parallel_halo(model%geometry%thck)

           if (model%options%whichtemp == TEMP_REMAP_ADV) then  ! Use IR to transport thickness, temperature
                                                                ! (and other tracers, if present)

              !If advecting other tracers, add parallel_halo update here
              call parallel_halo(model%temper%temp)
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

           if (write_verbose) then
              call distributed_gather_var(model%geometry%thck, gathered_thck)
              call distributed_gather_var(model%temper%temp, gathered_temp, &
                                          lbound(model%temper%temp,1), ubound(model%temper%temp,1))
           endif
        endif  ! old v. new remapping

        !Update halos of modified fields
        call parallel_halo(model%geometry%thck)
        if (model%options%whichtemp == TEMP_REMAP_ADV) then
           call parallel_halo(model%temper%temp)
           !If advecting other tracers, add parallel_halo update here
        endif

!whl - optional diagnostics - remove later
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

        ! NOTE finalization routine, to be written for PP HO core needs to be written (e.g.
        ! glam_velo_fordsiapstr_final( ) ), added to glam_strs2.F90, and called from glide_stop.F90

!***********************************************************************
end module glam
!***********************************************************************
