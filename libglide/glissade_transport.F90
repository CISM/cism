!=======================================================================
!BOP
!
! !MODULE: glissade_transport - driver for transport of mass and tracers
!
! !DESCRIPTION:
!
! Drivers for remapping and upwind ice transport
!
! author: William H. Lipscomb, LANL 
!
! This version was created from ice_transport_driver in CICE,
!  revision 313, 6 Jan. 2011.
! The repository is here: http://oceans11.lanl.gov/svn/CICE

! !INTERFACE:

      module glissade_transport
!
! !USES:
      use glimmer_global, only: dp
      use glimmer_log
      use glissade_remap, only: glissade_horizontal_remap, make_remap_mask, puny
!whl - Note to Jeff - no halo updates, but parallel code needed for global sums
!!      use parallel 

!whl - remove later
      use glimmer_paramets, only: thk0
!
!EOP
!
      implicit none
      save
      private
      public :: glissade_transport_driver, nghost_transport, ntracer_transport

      logical, parameter ::  &
         prescribed_area = .false.  ! if true, prescribe the area fluxed across each edge

!whl - to do - Should these be declared elsewhere?
      integer, parameter :: nghost_transport = 2
      integer, parameter :: ntracer_transport = 1

!=======================================================================

      contains

!=======================================================================
!BOP
!
! !IROUTINE: glissade_transport_driver - driver for mass and tracer transport
!
! !INTERFACE:
!
      subroutine glissade_transport_driver(dt,                   &
                                           dx,       dy,         &
                                           nx,       ny,         &
                                           nlyr,     sigma,      &
                                           nghost,   ntracer,    &
                                           uvel,     vvel,       &
                                           thck,     temp,       &
                                           age,                  &
                                           upwind_transport_in)
!
! !DESCRIPTION:
!
! This subroutine solves the transport equations for one timestep
! using the conservative remapping scheme developed by John Dukowicz
! and John Baumgardner (DB) and modified for sea ice by William
! Lipscomb and Elizabeth Hunke.
!
! This scheme preserves monotonicity of ice area and tracers.  That is,
! it does not produce new extrema.  It is second-order accurate in space,
! except where gradients are limited to preserve monotonicity. 
!
! Optionally, the remapping scheme can be replaced with a simple
! first-order upwind scheme.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
!
! !INPUT/OUTPUT PARAMETERS:
!
      real(dp), intent(in) ::     &
         dt,                   &! time step (s)
         dx, dy                 ! gridcell dimensions (m)
                                ! (cells assumed to be rectangular)

      integer, intent(in) ::   &
         nx, ny,               &! horizontal array size
         nlyr,                 &! number of vertical layers
         ntracer,              &! number of tracers (should be >= 1)
         nghost                 ! number of rows of ghost cells

      real(dp), intent(in), dimension(nlyr+1) :: &
         sigma                  ! layer interfaces in sigma coordinates
                                ! top sfc = 0, bottom sfc = 1

      !Note dimensions of uvel and vvel
      real(dp), intent(in), dimension(nlyr+1,nx-1,ny-1) :: &
         uvel, vvel             ! horizontal velocity components (m/s)
                                ! (defined at horiz cell corners, vertical interfaces)

      real(dp), intent(inout), dimension(nx,ny) :: &
         thck                   ! ice thickness (m)
                                ! (defined at horiz cell centers)

      real(dp), intent(inout), dimension(nlyr,nx,ny), optional :: &
         temp                   ! ice temperature
                                ! (defined at horiz cell centers, vertical layer midpts)

      real(dp), intent(inout), dimension(nlyr,nx,ny), optional :: &
         age                    ! ice age

      logical, intent(in), optional ::  &
         upwind_transport_in    ! if true, do first-order upwind transport
! 
!EOP
!
      ! local variables

      integer ::     &
         i, j, k         ,&! cell indices
         ilo,ihi,jlo,jhi ,&! beginning and end of physical domain
         nt                ! tracer index

      real (kind=dp), dimension (nx,ny) ::     &
         thck_mask         ! = 1. if ice is present, = 0. otherwise

      real (kind=dp), dimension (nx-1,ny-1) ::     &
         uvel_layer      ,&! uvel averaged to layer midpoint (m/s)
         vvel_layer        ! vvel averaged to layer midpoint (m/s)

      real (kind=dp), dimension (nx,ny,nlyr) ::     &
         thck_layer        ! ice layer thickness

      real (kind=dp), dimension (nx,ny,ntracer,nlyr) ::     &
         tracer            ! tracer values

      integer ::     &
         icells            ! number of cells with ice

      integer, dimension(nx*ny) ::     &
         indxi, indxj      ! compressed i/j indices

    !-------------------------------------------------------------------
    ! If prescribed_area is true, the area of each departure region is
    !  computed in advance (e.g., by taking the divergence of the 
    !  velocity field and passed to locate_triangles.  The departure 
    !  regions are adjusted to obtain the desired area.
    ! If false, edgearea is computed in locate_triangles and passed out.
    !-------------------------------------------------------------------

      real (kind=dp), dimension(nx,ny) ::   &
         edgearea_e     ,&! area of departure regions for east edges
         edgearea_n       ! area of departure regions for north edges

      logical, parameter ::     &
         conservation_check = .true. ! if true, check global conservation

      real (kind=dp) ::     &
         msum_init      ,&! initial global ice mass
         msum_final       ! final global ice mass

      real (kind=dp), dimension(ntracer) ::     &
         mtsum_init     ,&! initial global ice mass*tracer
         mtsum_final      ! final global ice mass*tracer

      logical ::     &
         l_stop            ! if true, abort the model

      character(len=100) :: message

      real (kind=dp), dimension (:,:,:), allocatable :: &
         worku            ! work array

      real (kind=dp), dimension(nx,ny) ::      &
         uee, vnn            ! cell edge velocities for upwind transport

      logical ::     &
         upwind_transport    ! if true, do first-order upwind transport

      if (present(upwind_transport_in)) then
         upwind_transport = upwind_transport_in
      else
         upwind_transport = .false.
      endif

!whl - convert thickness to m
!Note: This is not really necessary, but do this for now to minimize roundoff differences
!      in comparison to old remapping scheme

      thck(:,:) = thck(:,:) * thk0

!---!-------------------------------------------------------------------
!---! Prepare for remapping.
!---! Initialize, update ghost cells, fill tracer arrays.
!---!-------------------------------------------------------------------

!whl - Note to Jeff: pass these in as arguments or compute appropriately
!(ilo,ihi) and (jlo,jhi) are the lower and upper bounds of the local domain
! (i.e., grid cells owned by this processor)

      ilo = 1 + nghost
      ihi = nx - nghost
      jlo = 1 + nghost
      jhi = ny - nghost

      l_stop = .false.

    !-------------------------------------------------------------------
    ! Fill thickness and tracer arrays.
    ! Assume that temperature (if present) is tracer 1, and age (if present)
    !  is tracer 2.  Add more tracers as needed.
    ! If no tracers are present, then only the ice thickness is transported.
    !  In this case we define a dummy tracer array, since glissade_horizontal_remap
    !  requires that a tracer array be passed in.
    !-------------------------------------------------------------------

      do k = 1, nlyr

         thck_layer(:,:,k) = thck(:,:) * (sigma(k+1) - sigma(k))

         if (present(temp)) then
            tracer(:,:,1,k) = temp(k,:,:)
         else
            tracer(:,:,1,k) = 1.d0    ! dummy array
         endif

         if (present(age) .and. ntracer >= 2) tracer(:,:,2,k) = age(k,:,:)

      enddo

    !-------------------------------------------------------------------
    ! Compute initial values of globally conserved quantities (optional)
    !-------------------------------------------------------------------

      if (conservation_check) then

!whl - Note to Jeff: Replace with appropriate global sums.

         ! Compute initial values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         msum_init = sum (thck(:,:))  ! local for now
!         msum_init = global_sum (thck(:,:))  

         do nt = 1, ntracer
            mtsum_init(nt) =  sum (tracer(:,:,nt,:) * thck_layer(:,:,:))  ! local for now
!            mtsum_init(nt) =  global_sum (tracer(:,:,nt,:) * thck_layer(:,:,:))
         enddo

      endif                     ! conservation_check
      
      if (upwind_transport) then

         allocate (worku(nx,ny,0:ntracer))

         do k = 1, nlyr

            ! Average corner velocities at interfaces to edge velocities at layer midpoints.
      
            do j = jlo-1, jhi
            do i = ilo-1, ihi
               uee(i,j) = 0.5d0 * (uvel(k,i,j)   + uvel(k,i,j-1)    &
                                 + uvel(k+1,i,j) + uvel(k+1,i,j-1))
               vnn(i,j) = 0.5d0 * (vvel(k,i,j)   + vvel(k,i-1,j)    &
                                 + vvel(k+1,i,j) + vvel(k+1,i-1,j))
            enddo
            enddo

            ! Fill work array for transport

            worku(:,:,0) = thck_layer(:,:,k)
            do nt = 1, ntracer
               worku(:,:,nt) = thck_layer(:,:,k) * tracer(:,:,nt,k)
            enddo

            !-----------------------------------------------------------------
            ! Upwind transport
            !-----------------------------------------------------------------

            do nt = 0, ntracer
               call upwind_field (nx,             ny,                  &
                                  ilo, ihi,       jlo, jhi,            &
                                  dx,             dy,                  &
                                  dt,             worku(:,:,nt),       &
                                  uee(:,:),       vnn    (:,:))

            enddo   ! ntracer

            ! Recompute thickness and tracers

            thck_layer(:,:,k) = worku(:,:,0)
            do nt = 1, ntracer
               do j = jlo, jhi
               do i = ilo, ihi
                  if (thck_layer(i,j,k) > puny) then
                     tracer(i,j,nt,k) = worku(i,j,nt) / thck_layer(i,j,k)
                  else
                     tracer(i,j,nt,k) = 0.d0
                  endif
               enddo   ! i
               enddo   ! j
            enddo      ! ntracer

         enddo         ! nlyr

         deallocate (worku)

      else    ! remapping transport

      !-------------------------------------------------------------------
      ! Define a mask: = 1 where ice is present, = 0 otherwise         
      ! The mask is used to prevent tracer values in cells without ice from
      !  being used to compute tracer gradients.
      ! Compute here to avoid recomputing for each layer.
      !-------------------------------------------------------------------

         call make_remap_mask (nx,           ny,                 &
                               ilo, ihi,     jlo, jhi,           &
                               nghost,       icells,             &
                               indxi(:),     indxj(:),           &
                               thck(:,:),    thck_mask(:,:))

    !-------------------------------------------------------------------
    ! Remap ice thickness and tracers; loop over layers
    !-------------------------------------------------------------------

         do k = 1, nlyr

            ! Average velocities to the midpoint of the layer

            uvel_layer(:,:) = 0.5d0 * (uvel(k,:,:) + uvel(k+1,:,:))
            vvel_layer(:,:) = 0.5d0 * (vvel(k,:,:) + vvel(k+1,:,:))

            edgearea_e(:,:) = 0.d0
            edgearea_n(:,:) = 0.d0

            ! If prescribed_area is true, compute edgearea by taking the divergence
            !  of the velocity field.

            if (prescribed_area) then

               do j = jlo, jhi
               do i = ilo-1, ihi
                  edgearea_e(i,j) = (uvel_layer(i,j) + uvel_layer(i,j-1))     &
                                    * 0.5d0 * dy * dt
               enddo
               enddo
 
               do j = jlo-1, jhi
               do i = ilo, ihi
                  edgearea_n(i,j) = (vvel_layer(i,j) + vvel_layer(i-1,j))     &
                                   * 0.5d0 * dx * dt
               enddo
               enddo

            endif

         !-------------------------------------------------------------------
         ! Main remapping routine: Step ice thickness and tracers forward in time.
         !-------------------------------------------------------------------

            call glissade_horizontal_remap (dt,                                  &
                                            dx,                dy,               &
                                            nx,                ny,               &
                                            ntracer,           nghost,           &
                                            thck_mask(:,:),    icells,           &
                                            indxi(:),          indxj(:),         &
                                            uvel_layer(:,:),   vvel_layer(:,:),  &
                                            thck_layer(:,:,k), tracer(:,:,:,k),  &
                                            edgearea_e(:,:),   edgearea_n(:,:))

         enddo    ! nlyr

      endif       ! remapping v. upwind transport

      !-------------------------------------------------------------------
      ! Interpolate tracers back to sigma coordinates 
      !-------------------------------------------------------------------

      call glissade_vertical_remap(nx,                ny,      &
                                   nlyr,              ntracer, &
                                   sigma(:),                   &
                                   thck_layer(:,:,:),          &
                                   tracer(:,:,:,:) )

    !-------------------------------------------------------------------
    ! Recompute thickness, temperature and other tracers.
    !-------------------------------------------------------------------

      thck(:,:) = 0.d0

      do k = 1, nlyr

         thck(:,:) = thck(:,:) + thck_layer(:,:,k)

         if (present(temp)) temp(k,:,:) = tracer(:,:,1,k)
         if (present(age) .and. ntracer >= 2) age(k,:,:) = tracer(:,:,2,k)

      enddo

    !-------------------------------------------------------------------
    ! Compute final values of globally conserved quantities.
    ! Check global conservation of mass and mass*tracers.  (Optional)
    !-------------------------------------------------------------------

      if (conservation_check) then

         ! Compute final values of globally conserved quantities.
         ! Assume gridcells of equal area, ice of uniform density.

         msum_final = sum (thck(:,:))  ! local for now
!         msum_final = global_sum (thck(:,:))  

         do nt = 1, ntracer
            mtsum_final(nt) =  sum (tracer(:,:,nt,:) * thck_layer(:,:,:))  ! local for now
!            mtsum_final(nt) =  global_sum (tracer(:,:,nt,k) * thck_layer(:,:,k))
         enddo

!whl - Note to Jeff: Call this routine from master proc only
!         if (my_task == master_task) then
            call global_conservation (msum_init,     msum_final,      &
                                      l_stop,        ntracer,         &
                                      mtsum_init(:), mtsum_final(:))
            if (l_stop) then
               write(message,*) 'Aborting'
               call write_log(message,GM_FATAL)
            endif
!         endif

      endif                     ! conservation_check

!whl - convert thickness back to scaled units
!Note: This is not really necessary, but do this for now to minimize roundoff differences
!      in comparison to old remapping scheme

      thck(:,:) = thck(:,:) / thk0

      end subroutine glissade_transport_driver

!=======================================================================
!
!BOP
!
! !IROUTINE: global_conservation - check for changes in conserved quantities
!
! !INTERFACE:
!
      subroutine global_conservation (msum_init,  msum_final,     &
                                      l_stop,     ntracer,        &
                                      mtsum_init, mtsum_final)
!
! !DESCRIPTION:
!
! Check whether values of conserved quantities have changed.
! An error probably means that ghost cells are treated incorrectly.
!
! !REVISION HISTORY:
!
! author William H. Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      real (kind=dp), intent(in) ::     &
         msum_init   ,&! initial global ice mass 
         msum_final    ! final global ice mass

      logical, intent(inout) ::     &
         l_stop    ! if true, abort on return

      integer, intent(in) ::  &
         ntracer       ! number of tracers

      real (kind=dp), dimension(ntracer), intent(in), optional :: &
         mtsum_init  ,&! initial global ice mass*tracer
         mtsum_final   ! final global ice mass*tracer

      character(len=100) :: message
!
!EOP
!
      integer ::     &
           nt            ! tracer index

      real (kind=dp) ::     &
           diff          ! difference between initial and final values

      if (msum_init > puny) then
         diff = msum_final - msum_init
         if (abs(diff/msum_init) > puny) then
            l_stop = .true.
            write (message,*) 'glissade_transport: ice mass conservation error'
            call write_log(message)
            write (message,*) 'Initial global mass =', msum_init
            call write_log(message)
            write (message,*) 'Final global mass =', msum_final
            call write_log(message)
            write (message,*) 'Fractional error =', abs(diff)/msum_init
            call write_log(message)
         endif
      endif

      if (present(mtsum_init)) then
       do nt = 1, ntracer
         if (abs(mtsum_init(nt)) > puny) then
            diff = mtsum_final(nt) - mtsum_init(nt)
            if (abs(diff/mtsum_init(nt)) > puny) then
               l_stop = .true.
               write (message,*) 'glissade_transport: mass*tracer conservation error'
               call write_log(message)
               write (message,*) 'tracer index =', nt
               call write_log(message)
               write (message,*) 'Initial global mass*tracer =', mtsum_init(nt)
               call write_log(message)
               write (message,*) 'Final global mass*tracer =', mtsum_final(nt)
               call write_log(message)
               write (message,*) 'Fractional error =', abs(diff)/mtsum_init(nt)
               call write_log(message)
            endif
         endif
       enddo
      endif                     ! present(mtsum_init)

      end subroutine global_conservation

!----------------------------------------------------------------------

    subroutine glissade_vertical_remap(nx,       ny,        &
                                       nlyr,     ntracer,   &
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
         nx, ny,     &! number of cells in EW and NS directions
         nlyr,       &! number of vertical layers
         ntracer      ! number of tracer fields

    real(dp), dimension (nx, ny, nlyr), intent(inout) ::  &
         hlyr         ! layer thickness

    real(dp), dimension (nlyr+1), intent(in) ::  &
         sigma        ! sigma vertical coordinate (at layer interfaces)

    real(dp), dimension (nx, ny, ntracer, nlyr), intent(inout) ::   &
         trcr         ! tracer field to be remapped
                      ! tracer(k) = value at midpoint of layer k
 
    ! local variables
 
    integer :: i, j, k, k1, k2, nt
 
    real(dp), dimension (nx, ny, nlyr+1) ::  &
         zi1,        &! layer interfaces in old coordinate system
                      ! zi1(1) = 0. = value at top surface
                      ! zi1(k) = value at top of layer k
                      ! zi1(nlyr+1) = value at bottom surface (= 1 in sigma coordinates)
         zi2          ! layer interfaces in new coordinate system
                      ! Note: zi1(1) = zi2(1) = 0
                      !       zi1(nlyr+1) = zi2(nlyr+1) = 1

    real(dp), dimension(nx, ny) ::        &
         thck,       &! total thickness
         rthck        ! reciprocal of total thickness
 
    real(dp), dimension(nx,ny,ntracer,nlyr) ::       &
         htsum        ! sum of thickness*tracer in a layer         

    real(dp) :: zlo, zhi, hovlp

       !-----------------------------------------------------------------
       ! Compute total thickness and reciprocal thickness
       !-----------------------------------------------------------------

       thck(:,:) = 0.d0
       do k = 1, nlyr
          thck(:,:) = thck(:,:) + hlyr(:,:,k)
       enddo

       do j = 1, ny
          do i = 1, nx
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
       do k = 2, nlyr
          zi1(:,:,k) = zi1(:,:,k-1) +  hlyr(:,:,k-1)*rthck(:,:) 
       enddo
       zi1(:,:,nlyr+1) = 1.d0
                       
       !-----------------------------------------------------------------
       ! Set vertical coordinate zi2, given sigma.
       ! These are the coordinates to which we remap in the vertical.
       !-----------------------------------------------------------------

       zi2(:,:,1) = 0.d0
       do k = 2, nlyr
          zi2(:,:,k) = sigma(k)
       enddo
       zi2(:,:,nlyr+1) = 1.d0

       !-----------------------------------------------------------------
       ! Compute new layer thicknesses (zi2 coordinates)
       !-----------------------------------------------------------------

       do k = 1, nlyr
          hlyr(:,:,k) = (zi2(:,:,k+1) - zi2(:,:,k)) * thck(:,:)
       enddo

       !-----------------------------------------------------------------
       ! Compute sum of h*T for each new layer (k2) by integrating
       ! over the regions of overlap with old layers (k1).
       ! Note: It might be worth trying a more efficient
       !       search algorithm if the number of layers is large.
       !       This algorithm scales as nlyr^2.
       !       Also, may want to rearrange loop order if there are many tracers.
       !-----------------------------------------------------------------

       do k2 = 1, nlyr
          htsum(:,:,:,k2) = 0.d0 
          do k1 = 1, nlyr
             do nt = 1, ntracer
                do j = 1, ny
                   do i = 1, nx
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
 
       do k = 1, nlyr
          do nt = 1, ntracer
             do j = 1, ny
                do i = 1, nx
                   if (hlyr(i,j,k) > 0.d0) then
                      trcr(i,j,nt,k) = htsum(i,j,nt,k) / hlyr(i,j,k)
                   else
                      trcr(i,j,nt,k) = 0.d0
                   endif
                enddo   ! i
             enddo      ! j
          enddo         ! nt
       enddo            ! k

    end subroutine glissade_vertical_remap

!=======================================================================
!BOP
!
! !IROUTINE: upwind_field - first-order upwind transport
!
! !INTERFACE:
!
      subroutine upwind_field (nx,       ny,         &
                               ilo, ihi, jlo, jhi,   &
                               dx,       dy,         &
                               dt,       phi,        &
                               uee,      vnn)
!
! !DESCRIPTION:
!
! upwind transport algorithm
!
! !REVISION HISTORY:
!
! Authors: Elizabeth Hunke and William Lipscomb, LANL
!
! !USES:
!
! !INPUT/OUTPUT PARAMETERS:
!
      integer, intent (in) ::     &
         nx, ny             ,&! block dimensions
         ilo,ihi,jlo,jhi      ! beginning and end of physical domain

      real (kind=dp), intent(in) ::         &
         dx, dy             ,&! x and y gridcell dimensions
         dt                   ! time step

      real (kind=dp), dimension(nx,ny), &
         intent(inout) ::                       &
         phi                  ! scalar field

      real (kind=dp), dimension(nx,ny),         &
         intent(in)::     &
         uee, vnn             ! cell edge velocities

!
!EOP
!
      integer ::     &
         i, j                   ! standard indices

      real (kind=dp) ::        &
         upwind, y1, y2, a, h   ! function

      real (kind=dp), dimension(nx,ny) ::  &
         worka, workb

    !-------------------------------------------------------------------
    ! Define upwind function
    !-------------------------------------------------------------------

      upwind(y1,y2,a,h) = 0.5d0*dt*h*((a+abs(a))*y1+(a-abs(a))*y2)

    !-------------------------------------------------------------------
    ! upwind transport
    !-------------------------------------------------------------------

      do j = jlo-1, jhi
      do i = ilo-1, ihi
         worka(i,j)= upwind(phi(i,j),phi(i+1,j),uee(i,j),dy)
         workb(i,j)= upwind(phi(i,j),phi(i,j+1),vnn(i,j),dx)
      enddo
      enddo

      do j = jlo, jhi
      do i = ilo, ihi
         phi(i,j) = phi(i,j) - ( worka(i,j)-worka(i-1,j)      &
                               + workb(i,j)-workb(i,j-1) )    &
                               / (dx*dy)
      enddo
      enddo

      end subroutine upwind_field

!=======================================================================

      end module glissade_transport

!=======================================================================
