
module remap_glamutils      

  ! Contains utilities needed for using LANL incremental remapping code
  ! for thickness evolution with higher-order velocity solvers

    use glimmer_paramets, only: sp, dp, len0, thk0, tim0, vel0
    use glide_grids,      only: periodic_boundaries, periodic_boundaries_3d
    use xls

    type remap_glamutils_workspace
        ! arrays needed to pass GLAM variables to/from inc. remapping solver
        real (kind = dp), pointer, dimension(:,:,:) ::   &
            thck_ir,            &
            dew_ir,   dns_ir,   &
            dewt_ir,  dnst_ir,  &
            dewu_ir,  dnsu_ir,  &
            hm_ir,    tarea_ir, &
            ubar_ir,  vbar_ir
        real (kind = dp), pointer, dimension(:,:,:,:) :: trace_ir
        real (kind = dp) :: dt_ir

        ! *sfp* mask to apply for domains where the initial ice thickness limits are equivalent
        ! to the domain edge (e.g. the locations where bcs are applied). Application of this
        ! mask at the end of the advection calc essentially throws away any material that was
        ! transported beyond the boundaries of the initial domain (i.e. anywhere the ice thick-
        ! was initially zero will be forced back to zero, regardless of if material was 
        ! advected into it or not). This is mainly a hack to deal with problems on simplified
        ! domains.
        real (kind = dp), pointer, dimension(:,:) :: mask_ir

        integer :: ewn_ir, nsn_ir
    end type remap_glamutils_workspace

    contains

!----------------------------------------------------------------------

    subroutine horizontal_remap_init (wk, ewn, nsn, periodic_ew, periodic_ns)
    ! initialize variables for use in inc. remapping code   

      implicit none

      type(remap_glamutils_workspace) :: wk

      ! horizontal dimensions
      integer, intent(in) :: ewn, nsn                   

      ! flags from config file for applying periodic boundary conditions
      logical, intent(in) :: periodic_ew, periodic_ns   

      ! NOTE: number of tracers to be mapped (must be >1, even though for 
      ! now no tracers are being mapped)
      integer, parameter  :: ntrace = 1    

      wk%dt_ir = 0.0_dp     ! time step

      if (periodic_ew) then
        wk%ewn_ir = ewn + 2
     else
        wk%ewn_ir = ewn + 4
      end if

      if (periodic_ns) then
        wk%nsn_ir = nsn + 2
      else
        wk%nsn_ir = nsn + 4
      end if

      ! allocate arrays/vars 
      allocate( wk%thck_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%thck_ir = 0.0_dp
      allocate( wk%dew_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );          wk%dew_ir = 0.0_dp
      allocate( wk%dns_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );          wk%dns_ir = 0.0_dp
      allocate( wk%dewt_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%dewt_ir = 0.0_dp
      allocate( wk%dnst_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%dnst_ir = 0.0_dp
      allocate( wk%dewu_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%dewu_ir = 0.0_dp
      allocate( wk%dnsu_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%dnsu_ir = 0.0_dp
      allocate( wk%hm_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );           wk%hm_ir = 0.0_dp
      allocate( wk%tarea_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );        wk%tarea_ir = 0.0_dp
      allocate( wk%ubar_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%ubar_ir = 0.0_dp
      allocate( wk%vbar_ir(1:wk%ewn_ir,1:wk%nsn_ir,1) );         wk%vbar_ir = 0.0_dp
      allocate( wk%trace_ir(1:wk%ewn_ir,1:wk%nsn_ir,ntrace,1) ); wk%trace_ir = 0.0_dp
      allocate( wk%mask_ir(1:wk%ewn_ir,1:wk%nsn_ir) );           wk%mask_ir = 0.0_dp

    end subroutine horizontal_remap_init

!----------------------------------------------------------------------

    subroutine horizontal_remap_final( wk )
    ! deallocate variables for use in inc. remapping code   
    
      implicit none

      type(remap_glamutils_workspace) :: wk

      deallocate( wk%thck_ir )
      deallocate( wk%dew_ir );  deallocate( wk%dns_ir )
      deallocate( wk%dewt_ir ); deallocate( wk%dnst_ir )
      deallocate( wk%dewu_ir ); deallocate( wk%dnsu_ir )
      deallocate( wk%hm_ir )
      deallocate( wk%tarea_ir )
      deallocate( wk%ubar_ir ); deallocate( wk%vbar_ir )
      deallocate( wk%trace_ir )
      deallocate( wk%mask_ir )

    end subroutine horizontal_remap_final

!----------------------------------------------------------------------

    subroutine horizontal_remap_in( wk, dt,       thck,     &
                                    ntrace,   nghost,       &
                                    dew,      dns,          &
                                    uflx,     vflx,         &
                                    stagthck, thklim,       &
                                    periodic_ew, periodic_ns)
    ! put variables in from model into format that remapping code wants 

    implicit none

    type(remap_glamutils_workspace) :: wk


    integer, intent(out) :: ntrace, &   ! no. of tracers to be remapped   
                            nghost      ! no. of ghost cells

    real (kind = dp), dimension(:,:), intent(in) :: thck, uflx, vflx, stagthck
    real (kind = dp), intent(in) :: dew, dns, dt, thklim
    real (kind = dp) :: dt_cfl

    logical, intent(in) :: periodic_ew, periodic_ns

    integer  :: ewn,  nsn

    ! Number of *extra* ghost cells in ew, ns directions 
    ! (depending on whether periodic bcs are enabled) 
    integer  :: ngew, ngns 


    ! NOTE: number of tracers to be mapped (must be >1, even though for 
    ! now no tracers are being mapped)
    ntrace = 1

    ! No. of ghost cells only comes in if code is parallel. For remapping with serial
    ! code nghost=2 is ideal because then no boundary updates are needed
    nghost = 2

    ewn = size(thck, 1)
    nsn = size(thck, 2)

    ngew = (wk%ewn_ir - ewn)/2
    ngns = (wk%nsn_ir - nsn)/2
   
    where( thck*thk0 > thklim )
        wk%thck_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = thck(:,:)*thk0
    elsewhere
        wk%thck_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0d0
    end where
 
    wk%thck_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = thck(:,:)*thk0
    wk%dew_ir(:,:,1)  = dew*len0; wk%dns_ir(:,:,1) = dns*len0
    wk%dewt_ir(:,:,1) = dew*len0; wk%dnst_ir(:,:,1) = dns*len0
    wk%dewu_ir(:,:,1) = dew*len0; wk%dnsu_ir(:,:,1) = dns*len0
    wk%hm_ir(:,:,1) = 1.0_dp
    wk%tarea_ir = 1.0_dp / ( wk%dew_ir * wk%dns_ir )

    call write_xls("uflx.txt", uflx)
    call write_xls("vflx.txt", vflx)

    !Copy fluxes over
    !*tjb** - Some rejiggering is needed here, because, while IR and Glide both place
    !velocities on a B-grid, the IR B-grid is the same size as the A-grid
    !(that is, there is an extra row of cells on the extreme right and bottom
    !of the domain).
    !If periodic BCs are used, then the fluxes already have a row of ghost cells
    !on the left and right.  This means that, like thickness, a second row
    !is needed on the left.  However, unlike thickness, *two* extra rows are needed
    !on the right, to account for the extra B-grid row.  Same goes for top and bottom.
    where( stagthck > 0.0_dp )
        wk%ubar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = uflx/stagthck*vel0;
        wk%vbar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = vflx/stagthck*vel0;
    elsewhere
        wk%ubar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0_dp
        wk%vbar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0_dp
    endwhere

    call periodic_boundaries(wk%thck_ir(:,:,1), periodic_ew, periodic_ns, 2)
    call periodic_boundaries(wk%ubar_ir(:wk%ewn_ir-1,:wk%nsn_ir-1,1), periodic_ew, periodic_ns, 2)
    call periodic_boundaries(wk%vbar_ir(:wk%ewn_ir-1,:wk%nsn_ir-1,1), periodic_ew, periodic_ns, 2)

    !*tjb* Copy the extra set of ghost cells over
    !Hard coded 5 as the source for these ghost cells,
    !because we go in two rows for the low-end ghost cells,
    !then go in two more for the source of the first two
    !high-end ghost cells.
    if (periodic_ew) then
        wk%ubar_ir(wk%ewn_ir,:,:) = wk%ubar_ir(5,:,:)
        wk%vbar_ir(wk%ewn_ir,:,:) = wk%ubar_ir(5,:,:)
    end if

    if (periodic_ns) then
        wk%ubar_ir(:,wk%nsn_ir,:) = wk%ubar_ir(:,5,:)
        wk%vbar_ir(:,wk%nsn_ir,:) = wk%ubar_ir(:,5,:)
    end if

    ! eventually tracer arrays will be filled with temperature and other tracers
    wk%trace_ir(:,:,:,1) = 1.0_dp;
    wk%dt_ir = dt * tim0

    where( wk%thck_ir(:,:,1) > 0.0_dp )
        wk%mask_ir  = 1.0_dp
    end where


    !*tjb** Checks for IR's CFL-like condition.  These only print a warning for now.
    !Use the conservative, "highly divergent" criterion for now.
    dt_cfl = .5 * max(maxval(wk%dew_ir/wk%ubar_ir), maxval(wk%dns_ir/wk%vbar_ir))
    if (dt_cfl < wk%dt_ir) then
        write(*,*) "WARNING: CFL violation in incremental remapping scheme.  Time step should be <= ", dt_cfl
    end if

    end subroutine horizontal_remap_in

!----------------------------------------------------------------------

    subroutine horizontal_remap_out( wk, thck, acab, dt, periodic_ew, periodic_ns )
    ! take output from inc. remapping and put back into format that model likes

    implicit none

    type(remap_glamutils_workspace) :: wk


    real (kind = dp), intent(in) :: dt
    real (kind = sp), intent(in), dimension(:,:) :: acab
    real (kind = dp), dimension(:,:), intent(inout) :: thck
    logical, intent(in) :: periodic_ew, periodic_ns

    integer :: ewn, nsn, ngew, ngns

    ewn = size(thck, 1)
    nsn = size(thck, 2)

    ngew = (wk%ewn_ir - ewn)/2
    ngns = (wk%nsn_ir - nsn)/2

    call periodic_boundaries(wk%thck_ir(:,:,1), periodic_ew, periodic_ns, 2)

    !Map from IR thickness field back to Glide thickness field
    thck = wk%thck_ir(1+ngew:ngew+ewn, 1+ngns:ngns+nsn,1) / thk0
    
    !Apply accumulation
    thck = thck + acab*dt
    
    ! *sfp* Remove thickness from previously ice free locations.
    ! NOTE: this really only applys to a domain where we don't want the
    ! ice to advance via remapping. This is ok for now, as we don't have a 
    ! good scheme for calculating a marginal velocity (except in the case of 
    ! floating ice) in which case the flux at a lateral margin into a neighboring
    ! (previously) ice free grid cell would be suspect anyway. Nevertheless, this
    ! should probably be added as a config file option (apply the mask or not) if
    ! only to remind us that we are making that choice at present. 
    thck = thck * wk%mask_ir(1+ngew:ngew+ewn, 1+ngns:ngns+nsn)
    
    end subroutine horizontal_remap_out

!----------------------------------------------------------------------

end module remap_glamutils
