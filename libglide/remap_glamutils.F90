!whl - to do - Add standard Glimmer-CISM preamble here.

module remap_glamutils      

  ! Contains utilities needed for using LANL incremental remapping code
  ! for thickness evolution with higher-order velocity solvers

!whl - Modified Dec. 2010 to allow remapping of temperature and tracers
!      as well as ice thickness

    use glimmer_paramets, only: sp, dp, len0, thk0, tim0, vel0
    use xls

  !whl - to do - Hardwired ntrace and nghost for now.  Better to pass them in at initialization?
  !whl - The current remapping scheme requires nghost = 2.

    integer, parameter  :: ntrace_ir = 1    ! number of tracers (e.g., temperature, ice age)
    integer, parameter  :: nghost_ir = 2    ! number of ghost cell (halo) layers
    
    ! *sfp** arrays needed to pass GLAM variables to/from inc. remapping solver
    type remap_glamutils_workspace
        integer :: ewn_ir = 0      ! remapping grid dimensions
        integer :: nsn_ir = 0
        integer :: upn_ir = 1
        real (dp), pointer, dimension(:,:,:) ::   &
            thck_ir,             &
            dew_ir,   dns_ir,    &
            dewt_ir,  dnst_ir,   &
            dewu_ir,  dnsu_ir,   &
            hm_ir,    tarear_ir, &
            uvel_ir,  vvel_ir
        real (dp), pointer, dimension(:,:,:,:) :: trace_ir
        real (dp), pointer, dimension(:) :: dsigma_ir
        real (dp) :: dt_ir

        ! *sfp* mask to apply for domains where the initial ice thickness limits are equivalent
        ! to the domain edge (e.g. the locations where bcs are applied). Application of this
        ! mask at the end of the advection calc essentially throws away any material that was
        ! transported beyond the boundaries of the initial domain (i.e. anywhere the ice thick-
        ! was initially zero will be forced back to zero, regardless of if material was 
        ! advected into it or not). This is mainly a hack to deal with problems on simplified
        ! domains.
        real (kind = dp), pointer, dimension(:,:) :: mask_ir

    end type remap_glamutils_workspace

    contains

!----------------------------------------------------------------------

    subroutine horizontal_remap_init (wk,          &
                                      dew,   dns,  &
                                      ewn,   nsn,  &
                                      periodic_ew, periodic_ns, &
                                      upn,   sigma)
    ! initialize variables for use in inc. remapping code   

      implicit none

      type(remap_glamutils_workspace) :: wk

      real(dp), intent(in) :: dew, dns          ! horizontal grid cell size (m)
      integer,  intent(in) :: ewn, nsn          ! horizontal dimensions
      integer,  intent(in), optional :: upn     ! vertical dimensions
                                                ! (number of velocity levels)
      real(dp), dimension(:), intent(in), optional :: sigma  ! vertical coordinate

      ! flags from config file for applying periodic boundary conditions
      ! Currently not supported
      logical, intent(in) :: periodic_ew, periodic_ns

      integer :: k
      wk%dt_ir = 0.0_dp                         ! time step

      ! set remapping array sizes
      ! For IR, the scalar and velocity arrays are assumed to be the same size.

      wk%ewn_ir = ewn - 1   ! to make staggered and unstaggered grids the same size
      wk%nsn_ir = nsn - 1   ! 

!whl - If the optional input argument upn > 1, then temperature and layer thickness
!       will be remapped for each of (upn-1) layers.
!      In this case temperature and velocity are staggered in the vertical,
!       with velocities at sigma levels and temperatures in between.

      if (present(upn)) then
         wk%upn_ir = upn - 1   ! number of layers to be remapped                   
      else
         wk%upn_ir = 1
      endif

      ! allocate arrays/vars 
      !whl - IR assumes that the outer rows and columns of the domain are ghost cells.
      !      Given these array sizes, IR will not work correctly if the outer rows
      !       and columns are actually part of the physical domain with nonzero thickness.
      !      This problem will be fixed in the distributed parallel version of the code.

      !whl - to do - For second-order-accurate vertical remapping, include the upper and lower
      !              surface temperatures in the temperature array. 
      allocate( wk%dew_ir   (wk%ewn_ir,wk%nsn_ir,1) );            wk%dew_ir = 0.0_dp
      allocate( wk%dns_ir   (wk%ewn_ir,wk%nsn_ir,1) );            wk%dns_ir = 0.0_dp
      allocate( wk%dewt_ir  (wk%ewn_ir,wk%nsn_ir,1) );           wk%dewt_ir = 0.0_dp
      allocate( wk%dnst_ir  (wk%ewn_ir,wk%nsn_ir,1) );           wk%dnst_ir = 0.0_dp
      allocate( wk%dewu_ir  (wk%ewn_ir,wk%nsn_ir,1) );           wk%dewu_ir = 0.0_dp
      allocate( wk%dnsu_ir  (wk%ewn_ir,wk%nsn_ir,1) );           wk%dnsu_ir = 0.0_dp
      allocate( wk%hm_ir    (wk%ewn_ir,wk%nsn_ir,1) );             wk%hm_ir = 0.0_dp
      allocate( wk%tarear_ir(wk%ewn_ir,wk%nsn_ir,1) );         wk%tarear_ir = 0.0_dp
      allocate( wk%mask_ir  (ewn,   nsn) );                      wk%mask_ir = 0.0_dp

      allocate( wk%thck_ir (wk%ewn_ir,wk%nsn_ir,wk%upn_ir) );           wk%thck_ir = 0.0_dp
      allocate( wk%uvel_ir (wk%ewn_ir,wk%nsn_ir,wk%upn_ir) );           wk%uvel_ir = 0.0_dp
      allocate( wk%vvel_ir (wk%ewn_ir,wk%nsn_ir,wk%upn_ir) );           wk%vvel_ir = 0.0_dp
      allocate( wk%trace_ir(wk%ewn_ir,wk%nsn_ir,ntrace_ir,wk%upn_ir)); wk%trace_ir = 0.0_dp
      allocate( wk%dsigma_ir(wk%upn_ir));                             wk%dsigma_ir = 0.0_dp

      if (present(sigma)) then
         do k = 1, wk%upn_ir
            wk%dsigma_ir(k) = sigma(k+1) - sigma(k)
         enddo
      else
            wk%dsigma_ir(1) = 1.0_dp
      endif

    ! Set grid quantities

      wk%dew_ir (:,:,1) = dew*len0; wk%dns_ir (:,:,1) = dns*len0
      wk%dewt_ir(:,:,1) = dew*len0; wk%dnst_ir(:,:,1) = dns*len0
      wk%dewu_ir(:,:,1) = dew*len0; wk%dnsu_ir(:,:,1) = dns*len0
      wk%hm_ir  (:,:,1) = 1.0_dp
      wk%tarear_ir = 1.0_dp / ( wk%dew_ir * wk%dns_ir )

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
      deallocate( wk%tarear_ir )
      deallocate( wk%uvel_ir ); deallocate( wk%vvel_ir )
      deallocate( wk%trace_ir )
      deallocate( wk%mask_ir )
      deallocate( wk%dsigma_ir)

    end subroutine horizontal_remap_final

!----------------------------------------------------------------------

    subroutine horizontal_remap_in( wk,       dt,           &
                                    thck,                   &
                                    uflx,     vflx,         &
                                    stagthck, thklim,       &
                                    periodic_ew, periodic_ns, &
                                    uvel,     vvel,         &
                                    temp,     age           ) ! Periodic not currently supported
    ! put variables in from model into format that remapping code wants 

    ! *sfp** get GLAM variables in order for use in inc. remapping code   
    !whl - Optionally, can transport temperature and other tracers such as ice age.
    !      If desired, can add more 3D tracer fields.

    implicit none

    type(remap_glamutils_workspace) :: wk

! JCC - ntrace, nghost used for periodic, disabled
!     integer, intent(out) :: ntrace, &   ! no. of tracers to be remapped
!                             nghost      ! no. of ghost cells

    real (kind = dp), dimension(:,:), intent(in) :: thck, uflx, vflx, stagthck
!     real (kind = dp), intent(in) :: dew, dns, dt, thklim
    real (kind = dp), intent(in) :: dt, thklim
    real (kind = dp) :: dt_cfl
    real (kind = dp), dimension(:,:,:), intent(in), optional :: uvel, vvel
    real (kind = dp), dimension(:,:,:), intent(in), optional :: temp, age

    ! flags from config file for applying periodic boundary conditions
    ! Currently not supported
    logical, intent(in) :: periodic_ew, periodic_ns

    integer :: i, j, k                    ! indices

!     integer  :: ewn,  nsn ! JCC - Used for periodic, disabled
    ! Load thickness and velocities
    ! Note that the IR arrays have horizontal dimensions (ewn-1, nsn-1).
    ! IR assumes that k is the outer index (not the inner indes as in Glide).

    ! Number of *extra* ghost cells in ew, ns directions 
    ! (depending on whether periodic bcs are enabled) 
    ! integer  :: ngew, ngns

    if (wk%upn_ir > 1) then   ! Use IR to transport thickness and tracers in multiple layers,
                              !  assuming that temperature and velocity are staggered in the vertical.
                              ! Average the velocity field from layer interfaces to midpoints
       do k = 1, wk%upn_ir
          wk%thck_ir(:,:,k) = thck(:,:)*thk0 * wk%dsigma_ir(k)
          if (present(uvel)) wk%uvel_ir(:,:,k) = 0.5d0*(uvel(k,:,:) + uvel(k+1,:,:)) * vel0
          if (present(vvel)) wk%vvel_ir(:,:,k) = 0.5d0*(vvel(k,:,:) + vvel(k+1,:,:)) * vel0
       enddo

    else                ! one layer only

    ! NOTE: number of tracers to be mapped (must be >1, even though for 
    ! now no tracers are being mapped)
!     ntrace = 1  ! JCC - Used for periodic, disabled
!       where( thck > thklim )
!           wk%thck_ir(:,:,1) = thck(:,:)*thk0
!       elsewhere
!           wk%thck_ir(:,:,1) = 0.0d0
!       end where

    ! No. of ghost cells only comes in if code is parallel. For remapping with serial
    ! code nghost=2 is ideal because then no boundary updates are needed
!     nghost = 2 ! JCC - Used for periodic, disabled
       wk%thck_ir(:,:,1) = thck(:,:)*thk0

       where( stagthck > 0.0_dp )
           wk%uvel_ir(:,:,1) = uflx/stagthck * vel0
           wk%vvel_ir(:,:,1) = vflx/stagthck * vel0
       elsewhere
           wk%uvel_ir(:,:,1) = 0.0_dp
           wk%vvel_ir(:,:,1) = 0.0_dp
       endwhere

! JCC - Following used for periodic, disabled.
!     ngew = (wk%ewn_ir - ewn)/2
!     ngns = (wk%nsn_ir - nsn)/2
! 
!     where( thck > thklim )
!         wk%thck_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = thck(:,:)*thk0
!     elsewhere
!         wk%thck_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0d0
!     end where
!     
!     wk%dew_ir(:,:,1)  = dew*len0; wk%dns_ir(:,:,1) = dns*len0
!     wk%dewt_ir(:,:,1) = dew*len0; wk%dnst_ir(:,:,1) = dns*len0
!     wk%dewu_ir(:,:,1) = dew*len0; wk%dnsu_ir(:,:,1) = dns*len0
!     wk%hm_ir(:,:,1) = 1.0_dp
!     wk%tarea_ir = 1.0_dp / ( wk%dew_ir * wk%dns_ir )
! JCC - Above used for periodic, disabled.
      call write_xls("uflx.txt", uflx)
      call write_xls("vflx.txt", vflx)

    endif

    ! Load tracers (including temperature), if present
    ! Note: IR puts the vertical index on the outside

    if (present(temp)) then
       do k = 1, wk%upn_ir
          wk%trace_ir(:,:,1,k) = temp(k,:,:)
       enddo
    else
          wk%trace_ir(:,:,1,1) = 1.0_dp
    endif

! JCC - Following used for periodic, disabled.
    !Copy fluxes over
    !*tjb** - Some rejiggering is needed here, because, while IR and Glide both place
    !velocities on a B-grid, the IR B-grid is the same size as the A-grid
    !(that is, there is an extra row of cells on the extreme right and bottom
    !of the domain).
    !If periodic BCs are used, then the fluxes already have a row of ghost cells
    !on the left and right.  This means that, like thickness, a second row
    !is needed on the left.  However, unlike thickness, *two* extra rows are needed
    !on the right, to account for the extra B-grid row.  Same goes for top and bottom.
!     where( stagthck > 0.0_dp )
!         wk%ubar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = uflx/stagthck*vel0;
!         wk%vbar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = vflx/stagthck*vel0;
!     elsewhere
!         wk%ubar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0_dp
!         wk%vbar_ir(1+ngew:ngew+ewn,1+ngns:ngns+nsn,1) = 0.0_dp
!     endwhere
! 
!     call periodic_boundaries(wk%thck_ir(:,:,1), periodic_ew, periodic_ns, 2)
!     call periodic_boundaries(wk%ubar_ir(:wk%ewn_ir-1,:wk%nsn_ir-1,1), periodic_ew, periodic_ns, 2)
!     call periodic_boundaries(wk%vbar_ir(:wk%ewn_ir-1,:wk%nsn_ir-1,1), periodic_ew, periodic_ns, 2)
! 
    !Copy the extra set of ghost cells over
    !Hard coded 5 as the source for these ghost cells,
    !because we go in two rows for the low-end ghost cells,
    !then go in two more for the source of the first two
    !high-end ghost cells.
!     if (periodic_ew) then
!         wk%ubar_ir(wk%ewn_ir,:,:) = wk%ubar_ir(5,:,:)
!         wk%vbar_ir(wk%ewn_ir,:,:) = wk%ubar_ir(5,:,:)
!     end if
! 
!     if (periodic_ns) then
!         wk%ubar_ir(:,wk%nsn_ir,:) = wk%ubar_ir(:,5,:)
!         wk%vbar_ir(:,wk%nsn_ir,:) = wk%ubar_ir(:,5,:)
!     end if
! JCC - Above used for periodic, disabled.
    if (present(age) .and. ntrace_ir >= 2) then
       do k = 1, wk%upn_ir
          wk%trace_ir(:,:,2,k) = age(k,:,:)
       enddo
    endif

    !whl - Add other tracer fields here, if desired

    wk%dt_ir = dt * tim0

    ! Note: thck has dims (ewn-1,nsn-1), but mask has dims (ewn,nsn)

    wk%mask_ir(:,:)  = 0.0_dp
    do j = 1, wk%nsn_ir
    do i = 1, wk%ewn_ir
       if (thck(i,j) > 0.0_dp) then
          wk%mask_ir(i,j) = 1.0_dp
       endif
    enddo
    enddo

    !*tjb** Checks for IR's CFL-like condition.  These only print a warning for now.
    !Use the conservative, "highly divergent" criterion for now.

    dt_cfl = .5 * max(maxval(wk%dew_ir/wk%uvel_ir), maxval(wk%dns_ir/wk%vvel_ir))
    if (dt_cfl < wk%dt_ir) then
        write(*,*) "WARNING: CFL violation in incremental remapping scheme.  Time step should be <= ", dt_cfl
    end if

    end subroutine horizontal_remap_in

!----------------------------------------------------------------------

    subroutine horizontal_remap_out( wk,     thck,   &
                                     acab,   dt,     &
                                     temp,   age)

    ! take output from inc. remapping and put back into format that model likes

    implicit none

    type(remap_glamutils_workspace) :: wk


    real (kind = dp), intent(in) :: dt
    real (kind = sp), intent(in), dimension(:,:) :: acab
    real (kind = dp), dimension(:,:), intent(inout) :: thck
    real (kind = dp), dimension(:,:,:), intent(inout), optional :: temp, age

!    integer :: ewn, nsn, ngew, ngns

    integer :: i, j, k                 ! indices

    ! Map from IR thickness field back to Glide thickness field

    if (wk%upn_ir > 1) then

       thck(:,:) = 0.d0
       do k = 1, wk%upn_ir
          do j = 1, wk%nsn_ir
          do i = 1, wk%ewn_ir
             thck(i,j) = thck(i,j) + wk%thck_ir(i,j,k)
          enddo
          enddo
       enddo
       thck(:,:) = thck(:,:) / thk0

    else   ! upn_ir = 1

!!    thck(::) = wk%thck_ir(:,:,1) / thk0
       do j = 1, wk%nsn_ir
       do i = 1, wk%ewn_ir
          thck(i,j) = wk%thck_ir(i,j,1) / thk0
       enddo
       enddo

    endif   ! upn_ir > 1

    ! Map from IR tracer fields back to Glide fields

    if (present(temp)) then   
       do j = 1, wk%nsn_ir
       do i = 1, wk%ewn_ir
       do k = 1, wk%upn_ir
          temp(k,i,j) = wk%trace_ir(i,j,1,k)
       enddo
       enddo
       enddo
    endif   ! present(temp)

    if (present(age)) then   
       do j = 1, wk%nsn_ir   ! Note: (ewn,nsn) may be greater than (ewn_ir,nsn_ir)
       do i = 1, wk%ewn_ir
       do k = 1, wk%upn_ir
          age(k,i,j) = wk%trace_ir(i,j,2,k)
       enddo
       enddo
       enddo
    endif   ! present(age)

    !whl - Add other tracers here, if present

    !Apply accumulation
    thck = thck + acab*dt

!whl - to do - Allow ice to expand into places where it did not originally exist.
 
    ! *sfp* Remove thickness from previously ice free locations.
    ! NOTE: this really only applys to a domain where we don't want the
    ! ice to advance via remapping. This is ok for now, as we don't have a 
    ! good scheme for calculating a marginal velocity (except in the case of 
    ! floating ice) in which case the flux at a lateral margin into a neighboring
    ! (previously) ice free grid cell would be suspect anyway. Nevertheless, this
    ! should probably be added as a config file option (apply the mask or not) if
    ! only to remind us that we are making that choice at present. 

    thck = thck * wk%mask_ir
    
    !Do the same for temperature and tracers.

    if (present(temp)) then   
       do j = 1, wk%nsn_ir   ! Note: (ewn,nsn) may be greater than (ewn_ir,nsn_ir)
       do i = 1, wk%ewn_ir
       do k = 1, wk%upn_ir
          temp(k,i,j) = temp(k,i,j) * wk%mask_ir(i,j)
       enddo
       enddo
       enddo
    endif

    if (present(age)) then   
       do j = 1, wk%nsn_ir   ! Note: (ewn,nsn) may be greater than (ewn_ir,nsn_ir)
       do i = 1, wk%ewn_ir
       do k = 1, wk%upn_ir
          age(k,i,j) = age(k,i,j) * wk%mask_ir(i,j)
       enddo
       enddo
       enddo
    endif
    
    end subroutine horizontal_remap_out

!----------------------------------------------------------------------

    subroutine vertical_remap(ewn,      nsn,       &
                              upn,      ntrace,    &
                              sigma,    hlyr,      &
                              trcr)
 
    ! Conservative remapping of tracer fields from one set of vertical 
    ! coordinates to another.  The remapping is first-order accurate.
    !
    !whl - to do - Add a 2nd-order accurate vertical remapping scheme.
    !
    ! Author: William Lipscomb, LANL

    implicit none
 
    ! in-out arguments
 
    integer, intent(in) ::  &
         ewn, nsn,   &! number of cells in EW and NS directions
         upn,        &! number of vertical interfaces
         ntrace       ! number of tracer fields

    real(dp), dimension (ewn, nsn, upn-1), intent(inout) ::  &
         hlyr         ! layer thickness

    real(dp), dimension (upn), intent(in) ::  &
         sigma        ! sigma vertical coordinate

    real(dp), dimension (ewn, nsn, ntrace, upn-1), intent(inout) ::   &
         trcr         ! tracer field to be remapped
                      ! tracer(k) = value at midpoint of layer k
 
    ! local variables
 
    integer :: i, j, k, k1, k2, nt
 
    real(dp), dimension (ewn, nsn, upn) ::  &
         zi1,        &! layer interfaces in old coordinate system
                      ! zi1(1) = 0. = value at top surface
                      ! zi1(k) = value at top of layer k
                      ! zi1(nlyr+1) = value at bottom surface (= 1 in sigma coordinates)
         zi2          ! layer interfaces in new coordinate system
                      ! Note: zi1(1) = zi2(1) = 0
                      !       zi1(nlyr+1) = zi2(nlyr+1) = 1

    real(dp), dimension(ewn, nsn) ::        &
         thck,       &! total thickness
         rthck        ! reciprocal of total thickness
 
    real(dp), dimension(ewn,nsn,ntrace,upn-1) ::       &
         htsum        ! sum of thickness*tracer in a layer         

    real(dp) :: zlo, zhi, hovlp

       !-----------------------------------------------------------------
       ! Compute total thickness and reciprocal thickness
       !-----------------------------------------------------------------

       thck(:,:) = 0.d0
       do k = 1, upn-1
          thck(:,:) = thck(:,:) + hlyr(:,:,k)
       enddo

       do j = 1, nsn
          do i = 1, ewn
             if (thck(i,j) > 0.d0) then
                rthck(i,j) = 1.d0/thck(i,j)
             else
                rthck(i,j) = 0.d0
             endif
          enddo
       enddo

       !-----------------------------------------------------------------
       ! Determine vertical coordinate zi1, given input layer thicknesses.
       ! These are the coordinates from which we start.
       !-----------------------------------------------------------------

       zi1(:,:,1) = 0.d0
       do k = 2, upn-1
          zi1(:,:,k) = zi1(:,:,k-1) +  hlyr(:,:,k-1)*rthck(:,:) 
       enddo
       zi1(:,:,upn) = 1.d0
                       
       !-----------------------------------------------------------------
       ! Set vertical coordinate zi2, given sigma.
       ! These are the coordinates to which we remap in the vertical.
       !-----------------------------------------------------------------

       zi1(:,:,1) = 0.d0
       do k = 2, upn-1
          zi2(:,:,k) = sigma(k)
       enddo
       zi2(:,:,upn) = 1.d0

       !-----------------------------------------------------------------
       ! Compute new layer thicknesses (zi2 coordinates)
       !-----------------------------------------------------------------

       do k = 1, upn-1
          hlyr(:,:,k) = (zi2(:,:,k+1) - zi2(:,:,k)) * thck(:,:)
       enddo

       !-----------------------------------------------------------------
       ! Compute sum of h*T for each new layer (k2) by integrating
       ! over the regions of overlap with old layers (k1).
       ! Note: It might be worth trying a more efficient
       !       search algorithm if the number of layers is large.
       !       This algorithm scales as upn^2.
       !       Also, may want to rearrange loop order if there are many tracers.
       !-----------------------------------------------------------------

       do k2 = 1, upn-1
          htsum(:,:,:,k2) = 0.d0 
          do k1 = 1, upn-1
             do nt = 1, ntrace
                do j = 1, nsn
                   do i = 1, ewn
                      zhi = min (zi1(i,j,k1+1), zi2(i,j,k2+1)) 
                      zlo = max (zi1(i,j,k1), zi2(i,j,k2))
                      hovlp = max (zhi-zlo, 0.d0) * thck(i,j)
                      htsum(i,j,nt,k2) = htsum(i,j,nt,k2)    &
                                       +  trcr(i,j,nt,k1) * hovlp
                   enddo   ! i
                enddo      ! j
             enddo         ! nt
          enddo            ! k1
       enddo               ! k2
 
       !-----------------------------------------------------------------
       ! Compute tracer values in new layers
       !-----------------------------------------------------------------
 
       do k = 1, upn-1
          do nt = 1, ntrace
             do j = 1, nsn
                do i = 1, ewn
                   if (hlyr(i,j,k) > 0.d0) then
                      trcr(i,j,nt,k) = htsum(i,j,nt,k) / hlyr(i,j,k)
                   else
                      trcr(i,j,nt,k) = 0.d0
                   endif
                enddo   ! i
             enddo      ! j
          enddo         ! nt
       enddo            ! k

    end subroutine vertical_remap

!----------------------------------------------------------------------

end module remap_glamutils
