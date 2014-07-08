!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_velo_sia.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!
! This module contains routines for computing the ice velocity
!  using the shallow-ice approximation.
!
! It is roughly based on module glissade_velo_higher but is much simpler, 
!  computing only the SIA velocities. 
! It is called with whichdycore = DYCORE_GLISSADE, which_ho_approx = SIMPLE_APPROX_SIA. 
!
! The calculation is similar to the one in Glide, except that it
!  computes only the SIA velocity profile and does not evolve thickness.  
! Thickness must be evolved separately using the glissade_transport module.
!
! Unlike the Glide SIA calculation, this calculation is fully parallel.
!
! The function of the module is to allow a parallel SIA calculation using an 
!  explicit transport scheme, instead of an implicit diffusion calculation as in Glide.
! This can also be done using glissade_velo_higher, setting which_ho_approx = HO_APPROX_SIA.
!  This uses the same numerical techniques as the higher-order solve but with only the
!  SIA matrix elements.  However, HO_APPROX_SIA is more expensive than SIMPLE_APPROX_SIA, 
!  and is somewhat less accurate for the Halfar test problem.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

  module glissade_velo_sia

    use glimmer_global, only: dp
    use glimmer_physcon, only: gn, rhoi, grav, scyr
    use glimmer_paramets, only: thk0, len0, vel0, vis0, tau0
!    use glimmer_log, only: write_log

    use glide_types

    use parallel

    implicit none

    private
    public :: glissade_velo_sia_solve

!WHL - debug
    logical, parameter :: verbose = .false.
    logical, parameter :: verbose_geom = .false.
    logical, parameter :: verbose_bed = .false.
    logical, parameter :: verbose_interior = .false.

  contains

!****************************************************************************

  subroutine glissade_velo_sia_solve(model,                &
                                     nx,     ny,     nz)

    use glissade_grid_operators, only: glissade_stagger, glissade_centered_gradient
    use glissade_temp, only: glissade_calcbpmp
    
!TODO - Remove nx, ny, nz from argument list?
!       Would then have to allocate some local arrays.

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    type(glide_global_type), intent(inout) :: model   ! derived type holding ice-sheet info

    !----------------------------------------------------------------
    ! Note that the glissade solver uses SI units.
    ! Thus we have grid cell dimensions and ice thickness in meters,
    !  velocity in m/s, and the rate factor in Pa^(-n) s(-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Note: nx and ny are the horizontal dimensions of scalar arrays (e.g., thck and temp).
    !       The velocity arrays have horizontal dimensions (nx-1, ny-1).
    !       nz is the number of levels at which uvel and vvel are computed.
    !       The scalar variables generally live at layer midpoints and have
    !         vertical dimension nz-1.
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,               &  ! number of grid cells in each direction
       nz                       ! number of vertical levels where velocity is computed
                                ! (same as model%general%upn)

    !----------------------------------------------------------------
    ! Local variables and pointers set to components of model derived type 
    !----------------------------------------------------------------

    real(dp) ::  &
       dx,  dy                  ! grid cell length and width (m)
                                ! assumed to have the same value for each grid cell

    real(dp), dimension(:), pointer :: &
       sigma                    ! vertical sigma coordinate, [0,1]

    real(dp)  ::   & 
       thklim, &                ! minimum ice thickness for active cells (m)
       btrc_const               ! constant basal traction ((m/yr)/Pa) for whichbtrc options

    integer :: &
       whichbtrc                ! basal traction option for SIA
                                ! Note: Several, but not all, of the Glide options are supported

    real(dp), dimension(:,:), pointer ::  &
       thck,                 &  ! ice thickness (m)
       usrf,                 &  ! upper surface elevation (m)
       bwat,                 &  ! basal water depth (m)
       btrc                     ! basal traction parameter (m/yr)/Pa), = 1/beta

    real(dp), dimension(:,:,:), pointer ::  &
       uvel, vvel,  &           ! velocity components (m/yr)
       temp,        &           ! temperature (deg C)
       flwa                     ! flow factor in units of Pa^(-n) yr^(-1)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    real(dp), dimension(nx,ny) ::  &
       bpmp                   ! basal pressure melting point temperature (deg C)

    real(dp), dimension(nx-1,ny-1) ::   &
       stagthck,            & ! ice thickness averaged to vertices (m)
       dusrf_dx, dusrf_dy,  & ! gradient of upper surface elevation (m/m)
       stagbwat,            & ! basal water depth averaged to vertices (m)
       stagbtemp,           & ! basal ice temperature averaged to vertices (deg C)
       stagbpmp,            & ! bpmp averaged to vertices (deg C)
       ubas,     vbas         ! basal velocity components (m/yr)

    real(dp), dimension(nz-1,nx-1,ny-1) ::  &
       stagflwa               ! flwa averaged to staggered grid, Pa^(-n) yr^(-1)

    integer, dimension(nx,ny) ::     &
       imask                  ! = 1 where ice is present, else = 0

    integer :: k

!WHL - debug
    integer :: i, j

    !--------------------------------------------------------
    ! Assign local pointers and variables to derived type components
    !--------------------------------------------------------

!    nx = model%general%ewn
!    ny = model%general%nsn
!    nz = model%general%upn

     dx = model%numerics%dew
     dy = model%numerics%dns

     thklim = model%numerics%thklim
     btrc_const = model%velowk%btrac_const
     whichbtrc = model%options%whichbtrc

     sigma    => model%numerics%sigma(:)
     thck     => model%geometry%thck(:,:)
     usrf     => model%geometry%usrf(:,:)
     bwat     => model%temper%bwat(:,:)
     btrc     => model%velocity%btrc(:,:)
     temp     => model%temper%temp(:,:,:)
     flwa     => model%temper%flwa(:,:,:)
     uvel     => model%velocity%uvel(:,:,:)
     vvel     => model%velocity%vvel(:,:,:)

    !--------------------------------------------------------
    ! Convert input variables to appropriate units for this solver.
    ! (Mainly SI, except that time units in flwa, velocities,
    !  and btrc are years instead of seconds)
    !--------------------------------------------------------

    call glissade_velo_sia_scale_input(dx,     dy,            &
                                       thck,   usrf,          &
                                       thklim, flwa,          &
                                       bwat,   btrc_const,    &
                                       uvel,   vvel)


    ! WHL: Here we follow Glide in computing stagthck and stagflwa at cell vertices
    !      (where Glide computes the diffusivity).  It would be possible to stagger
    !      these quantities to cell edges instead.  I do not know if this would
    !      make much difference in terms of accuracy and stability.
    !       

    where (model%geometry%thck > model%numerics%thklim)
       imask = 1
    elsewhere
       imask = 0
    endwhere


    ! Compute staggered variables
    !
    ! stag_flag_in = 0 gives Glide-style averaging
    !  (ice-free cells are included in the average)
    ! stag_flag_in = 1 omits ice-free cells from the average
    !
    ! Setting stag_flag = 1 for stagthck gives slightly more accurate results 
    !  for the Halfar SIA test than does stag_flag = 0.
    !
    ! Note: Glide in effect has stag_flag_in = 0 for all staggered variables.
    ! This seems wrong for temp, bwat, bpmp and flwa.

    call glissade_stagger(nx,      ny,       &
                          thck,    stagthck, &
                          imask,   stag_flag_in = 1)

    do k = 1, nz-1
       call glissade_stagger(nx,          ny,               &
                             flwa(k,:,:), stagflwa(k,:,:),  &
                             imask,       stag_flag_in = 1)
    enddo

    if (whichbtrc == BTRC_CONSTANT_BWAT) then

       ! stag_flag_in = 1 omits empty cells from the average

       call glissade_stagger(nx,           ny,              &
                             bwat(:,:),    stagbwat(:,:),   &
                             imask,        stag_flag_in = 1)

    elseif (whichbtrc == BTRC_CONSTANT_TPMP) then

       call glissade_stagger(nx,           ny,         &
                             temp(nz,:,:), stagbtemp,  &
                             imask,        stag_flag_in = 1)
       
       ! Note: glissade_calcbpmp expects dimensionless thickness

       call glissade_calcbpmp(nx,             ny,    &
                              thck(:,:)/thk0, bpmp(:,:))

       call glissade_stagger(nx,           ny,           &
                             bpmp(:,:),    stagbpmp(:,:), &
                             imask,        stag_flag_in = 1)

    endif    ! whichbtrc

    ! Compute surface elevation gradient
    !
    ! grad_flag_in = 0 gives Glide-style gradient
    !  (ice-free cells included in the gradient)

    call glissade_centered_gradient(nx,       ny,         &
                                    dx,       dy,         &
                                    usrf,                 &
                                    dusrf_dx, dusrf_dy,   &
                                    imask,    grad_flag_in = 0)

!WHL - debug
    if (verbose .and. main_task) then
       print*, ' '
       print*, 'In glissade_velo_sia_solve'
    endif

    if (verbose_geom .and. main_task) then

       print*, ' '
       print*, 'stagthck (m):'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.2)',advance='no') stagthck(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'dusrf_dx:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') dusrf_dx(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'dusrf_dy:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') dusrf_dy(i,j)
          enddo
          print*, ' '
       enddo

    endif

    ! Compute velocity at the bed (ubas, vbas)

    call glissade_velo_sia_bed(nx,         ny,             &
                               stagthck,   thklim,         &
                               dusrf_dx,   dusrf_dy,       &
                               whichbtrc,  stagbwat,       &
                               stagbtemp,  stagbpmp,       &
                               btrc,       btrc_const,     &
                               ubas,       vbas)

    if (verbose_bed .and. main_task) then

       print*, ' '
       print*, 'whichbtrc, btrc_const =', whichbtrc, btrc_const

       print*, ' '
       print*, 'btrc:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') btrc(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'ubas:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.2)',advance='no') ubas(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'vbas:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.2)',advance='no') vbas(i,j)
          enddo
          print*, ' '
       enddo

    endif  ! verbose_bed

    ! Compute velocity in the ice interior

    call glissade_velo_sia_interior(nx,       ny,       nz,  &
                                    dx,       dy,            &
                                    sigma,    thklim,        &
                                    usrf,     stagthck,     &
                                    dusrf_dx, dusrf_dy,     &
                                    stagflwa,               &
                                    ubas,     vbas,          &
                                    uvel,     vvel)

    if (verbose_interior .and. main_task) then

       k = 1

       print*, ' '
       print*, 'uvel, k = 1:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.2)',advance='no') uvel(k,i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'vvel, k = 1:'
       do i = 1, nx-1
          write(6,'(i7)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i4)',advance='no') j
          do i = 1, nx-1
             write(6,'(f7.2)',advance='no') vvel(k,i,j)
          enddo
          print*, ' '
       enddo

    endif   ! verbose_interior

    ! Convert back to dimensionless units before returning

    call glissade_velo_sia_scale_output(thck,    usrf,       &
                                        flwa,                &
                                        bwat,    btrc,       &
                                        uvel,    vvel)

  end subroutine glissade_velo_sia_solve

!*********************************************************************

  subroutine glissade_velo_sia_scale_input(dx,     dy,            &
                                           thck,   usrf,          &
                                           thklim, flwa,          &
                                           bwat,   btrc_const,    &
                                           uvel,   vvel)

    !--------------------------------------------------------
    ! Convert input variables (generally dimensionless)
    ! to appropriate units for the glissade_velo_sia solver.
    !--------------------------------------------------------

    real(dp), intent(inout) ::   &
       dx, dy                  ! grid cell length and width 

    real(dp), dimension(:,:), intent(inout) ::   &
       thck,                &  ! ice thickness
       usrf                    ! upper surface elevation

    real(dp), intent(inout) ::   &
       thklim,  &              ! minimum ice thickness for active cells
       btrc_const              ! constant basal traction ((m/yr)/Pa) for whichbtrc options

    real(dp), dimension(:,:,:), intent(inout) ::  &
       flwa                    ! flow factor in units of Pa^(-n) yr^(-1)

    real(dp), dimension(:,:), intent(inout)  ::  &
       bwat                    ! basal water depth (m)

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel              ! velocity components (m/yr)

    ! grid cell dimensions: rescale from dimensionless to m
    dx = dx * len0
    dy = dy * len0

    ! ice geometry: rescale from dimensionless to m
    thck = thck * thk0
    usrf = usrf * thk0
    thklim = thklim * thk0

    ! rate factor: rescale from dimensionless to Pa^(-n) yr^(-1)
    flwa = flwa * (vis0*scyr)

    ! bwat: rescale from dimensionless to m
    bwat = bwat * thk0

    ! btrc: rescale from dimensionless to (m/yr)/Pa
!    btrc_const = btrc_const * (vel0*scyr) / tau0
    btrc_const = btrc_const * (vel0*scyr) * len0 / thk0**2

    ! ice velocity: rescale from dimensionless to m/yr
    uvel = uvel * (vel0*scyr)
    vvel = vvel * (vel0*scyr)

    end subroutine glissade_velo_sia_scale_input

!*********************************************************************

    subroutine glissade_velo_sia_scale_output(thck,    usrf,     &
                                              flwa,              &
                                              bwat,    btrc,     &
                                              uvel,    vvel)

    !--------------------------------------------------------
    ! Convert output variables to appropriate Glimmer-CISM units
    ! (generally dimensionless)
    !--------------------------------------------------------

    real(dp), dimension(:,:), intent(inout) ::  &
       thck,                 &  ! ice thickness
       usrf                     ! upper surface elevation

    real(dp), dimension(:,:,:), intent(inout) ::  &
       flwa                     ! flow factor in units of Pa^(-n) yr^(-1)

    real(dp), dimension(:,:), intent(inout)  ::  &
       bwat,  &                 ! basal water depth (m)
       btrc                     ! basal traction parameter ((m/yr)/Pa)

    real(dp), dimension(:,:,:), intent(inout) ::  &
       uvel, vvel               ! velocity components (m/yr)

    ! Convert geometry variables from m to dimensionless units
    thck = thck / thk0
    usrf = usrf / thk0

    ! Convert flow factor from Pa^(-n) yr^(-1) to dimensionless units
    flwa = flwa / (vis0*scyr)

    ! Convert bwat from m to dimensionless units
    bwat = bwat / thk0

    ! Convert btrc from (m/yr)/Pa to dimensionless units
    btrc = btrc / ((vel0*scyr)/tau0)

    ! Convert velocity from m/yr to dimensionless units
    uvel = uvel / (vel0*scyr)
    vvel = vvel / (vel0*scyr)

    end subroutine glissade_velo_sia_scale_output

!*********************************************************************

  subroutine glissade_velo_sia_bed(nx,         ny,             &
                                   stagthck,   thklim,         &
                                   dusrf_dx,   dusrf_dy,       &
                                   whichbtrc,  stagbwat,       &
                                   stagbtemp,  stagbpmp,       &
                                   btrc,       btrc_const,     &
                                   ubas,       vbas)

    !----------------------------------------------------------------
    ! Compute the basal traction coefficient (btrc = 1/beta) and
    ! the resulting basal velocities, assuming these velocities are
    ! a linear function of the gravitational driving stress.
    !
    ! Note: Not all the Glide whichbtrc options are supported, but
    ! we support the ones needed for the EISMINT tests.
    !----------------------------------------------------------------
 
    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! number of grid cells in each direction

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       stagthck,              & ! ice thickness averaged to vertices (m)
       dusrf_dx, dusrf_dy,    & ! gradient of upper surface elevation (m/m)
       stagbwat,              & ! basal water depth averaged to vertices (m)
       stagbtemp,             & ! basal temperature averaged to vertices (deg C)
       stagbpmp                 ! basal pressure melting point temperature averaged to vertices (deg C)

    real(dp), intent(in) ::   &
       thklim                   ! minimum ice thickness for active cells

    integer, intent(in) :: &
       whichbtrc                ! basal traction option for SIA
                                ! Note: Several, but not all, of the Glide options are supported
                                !       We support the ones used for EISMINT

    real(dp), intent(inout) :: & 
       btrc_const               ! constant basal traction ((m/yr)/Pa) for whichbtrc options

    real(dp), dimension(nx-1,ny-1), intent(out) ::   &
       btrc,                  & ! basal traction parameter ((m/yr)/Pa), = 1/beta
       ubas,     vbas           ! basal velocity components (m/yr)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    ! Compute ubas and vbas at cell vertices, as in Glide.
    ! TODO - Test computing ubas and vbas at cell edges?

    do j = 1, ny-1
       do i = 1, nx-1

          if (stagthck(i,j) > thklim) then

             ! Compute basal traction coefficient, btrc = 1/beta

             select case(whichbtrc)

             case(BTRC_CONSTANT)

                btrc(i,j) = btrc_const

             case(BTRC_CONSTANT_BWAT)

                ! btrc is constant where basal melt water is present, else no slip
                ! This option can be used for EISMINT-2 experiment H, provided that 
                ! basal water is present where T = Tpmp (e.g., BWATER_LOCAL)

                if (stagbwat(i,j) > 0.d0) then
                   btrc(i,j) = btrc_const
                else
                   btrc(i,j) = 0.d0
                end if

             case(BTRC_CONSTANT_TPMP)

                ! constant where basal temperature equal to pressure melting point, else = 0
                ! This is the actual condition for EISMINT-2 experiment H, which may not be 
                ! the same as case BTRC_CONSTANT_BWAT above, depending on the hydrology

                if (abs(stagbpmp(i,j) - stagbtemp(i,j)) < 1.d-3) then
                   btrc(i,j) = btrc_const
                else
                   btrc(i,j) = 0.d0
                end if
                
             case default  ! includes BTRC_ZERO

                ! no sliding
                ! This is used for EISMINT-2 experiments A to F

                btrc(i,j) = 0.d0

             end select

             ! Compute basal velocity as a linear function of gravitational driving stress

             ubas(i,j) = -btrc(i,j) * rhoi * grav * stagthck(i,j) * dusrf_dx(i,j)
             vbas(i,j) = -btrc(i,j) * rhoi * grav * stagthck(i,j) * dusrf_dy(i,j)

          else   ! stagthck <= thklim

             btrc(i,j) = 0.d0
             ubas(i,j) = 0.d0
             vbas(i,j) = 0.d0

          endif

       end do
    end do

  end subroutine glissade_velo_sia_bed

!*********************************************************************

  subroutine glissade_velo_sia_interior(nx,       ny,      nz,  &
                                        dx,       dy,           &
                                        sigma,    thklim,       &
                                        usrf,     stagthck,     &
                                        dusrf_dx, dusrf_dy,     &
                                        stagflwa,               &
                                        ubas,     vbas,         &
                                        uvel,     vvel)

    use parallel

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,                & ! number of grid cells in each direction
       nz                       ! number of vertical levels where velocity is computed

    real(dp), intent(in) ::   &
       dx, dy,                & ! grid cell length and width (m)
       thklim                   ! minimum ice thickness for active cells

    real(dp), dimension(nz) ::   &
       sigma                    ! vertical sigma coordinate, [0,1]

    real(dp), dimension(nx,ny), intent(in) ::   &
       usrf                     ! upper surface elevation (m)

    real(dp), dimension(nx-1,ny-1), intent(in) ::   &
       stagthck,              & ! ice thickness averaged to vertices (m)
       dusrf_dx, dusrf_dy,    & ! gradient of upper surface elevation at vertices (m/m)
       ubas, vbas               ! basal velocity components at vertices (m/yr)

    real(dp), dimension(nz-1, nx-1,ny-1), intent(in) ::   &
       stagflwa                 ! flwa averaged to vertices (Pa^(-n) yr^(-1))

    real(dp), dimension(nz, nx-1,ny-1), intent(out) ::   &
       uvel, vvel               ! velocity components at vertices (m/yr)

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j, k

    real(dp) ::   &
       siafact,               & ! SIA premultiply factor
       dsdx_edge, dsdy_edge     ! surface elevation gradient at cell edges

    real(dp), dimension(nz,nx-1,ny-1) ::   &
       vintfact,              & ! vertically integrated SIA factor
       uedge, vedge             ! velocity components at edges (m/yr)
                                ! u on E edge, v on N edge (C grid)

!WHL - debug
    real(dp), dimension(nx-1,ny-1) :: diffu, uflx, vflx

    ! Compute vertically integrated factor for velocity calculation.
    ! As in Glide, this factor is located at cell vertices.
    ! TODO: Try putting this factor on cell edges?

    do j = 1, ny-1
       do i = 1, nx-1
          
          if (stagthck(i,j) > thklim) then

             siafact = 2.d0 * (rhoi*grav)**gn * stagthck(i,j)**(gn+1)       &
                            * (dusrf_dx(i,j)**2 + dusrf_dy(i,j)**2) ** ((gn-1)/2)

             vintfact(nz,i,j) = 0.d0

             do k = nz-1, 1, -1

                vintfact(k,i,j) = vintfact(k+1,i,j) -                              &
                                  siafact * stagflwa(k,i,j)                       &
                                          * ((sigma(k) + sigma(k+1))/2.d0) ** gn  &
                                          * (sigma(k+1) - sigma(k))

             enddo   ! k

          else   ! stagthck <= thklim

             vintfact(:,i,j) = 0.d0

          endif

       enddo      ! i
    enddo         ! j

    ! Compute ice velocity components at cell edges
    ! (u at E edge, v at N edge)

    do j = 1, ny-1
       do i = 1, nx-1

          ! east edges

          if (stagthck(i,j) > thklim .and. stagthck(i,j-1) > thklim) then

             uedge(nz,i,j) = (ubas(i,j) + ubas(i,j-1)) / 2.d0

             dsdx_edge = (usrf(i+1,j) - usrf(i,j)) / dx

             do k = nz-1, 1, -1

                uedge(k,i,j) = uedge(nz,i,j) +                                  &
                                (vintfact(k,i,j) + vintfact(k,i,j-1))/2.d0  * dsdx_edge

             enddo   ! k

          else  ! stagthck < thklim

             uedge(:,i,j) = 0.d0

          endif

          ! north edges

          if (stagthck(i,j) > thklim .and. stagthck(i-1,j) > thklim) then

             vedge(nz,i,j) = (vbas(i,j) + vbas(i-1,j)) / 2.d0

             dsdy_edge = (usrf(i,j+1) - usrf(i,j)) / dy

             do k = nz-1, 1, -1

                vedge(k,i,j) = vedge(nz,i,j) +                                  &
                                (vintfact(k,i,j) + vintfact(k,i-1,j))/2.d0 * dsdy_edge

             enddo   ! k

          else  ! stagthck < thklim

             vedge(:,i,j) = 0.d0

          endif

          if (verbose_interior .and. main_task) then
             !WHL - debug - Compute diffusivitity (as defined by Glide) at vertex(i,j)
             diffu(i,j) = 0.d0
             do k = 1, nz-1
                diffu(i,j) = diffu(i,j) - (vintfact(k,i,j) + vintfact(k+1,i,j))/2.d0 * (sigma(k+1) - sigma(k)) * stagthck(i,j)
             enddo

             !WHL - debug - Compute uflx and vflx that would be associated with uedge and vedge
             uflx(i,j) = 0.d0 !uedge(nz,i,j)*stagthck(i,j)
             vflx(i,j) = 0.d0 !uedge(nz,i,j)*stagthck(i,j)
             do k = 1, nz-1
                uflx(i,j) = uflx(i,j) + uedge(k,i,j) * (sigma(k+1) - sigma(k)) * stagthck(i,j)
                vflx(i,j) = vflx(i,j) + vedge(k,i,j) * (sigma(k+1) - sigma(k)) * stagthck(i,j)
             enddo
          endif

       enddo     ! i
    enddo        ! j

    ! Average edge velocities to vertices
    ! Do this for locally owned vertices only, then do halo update

    do j = nhalo+1, ny-nhalo
       do i = nhalo+1, nx-nhalo

          uvel(:,i,j) = (uedge(:,i,j) + uedge(:,i,j+1)) / 2.d0 
          vvel(:,i,j) = (vedge(:,i,j) + vedge(:,i+1,j)) / 2.d0 
 
       enddo
    enddo

    call staggered_parallel_halo(uvel)
    call staggered_parallel_halo(vvel)

!WHL - debug
    if (verbose_interior .and. main_task) then
       print*, ' '
       print*, 'diffu (m^2/yr):'
       do i = 1, nx-1
          write(6,'(i8)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i3)',advance='no') j
          do i = 1, nx-1
             write(6,'(f8.0)',advance='no') diffu(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'uflx (m^2/yr):'
       do i = 1, nx-1
          write(6,'(i8)',advance='no') i
       enddo
       print*, ' '
       do j = ny-1, 1, -1
          write(6,'(i3)',advance='no') j
          do i = 1, nx-1
             write(6,'(f8.0)',advance='no') uflx(i,j)
          enddo
          print*, ' '
       enddo
    endif   ! verbose_interior

  end subroutine glissade_velo_sia_interior

!*********************************************************************

  end module glissade_velo_sia

!*********************************************************************
