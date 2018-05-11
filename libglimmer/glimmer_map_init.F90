!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_map_init.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2018
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

!> initialise map projection routines
module glimmer_map_init

  use glimmer_map_types

  implicit none

contains

  !> read projection configuration from file
  subroutine glimmap_readconfig(proj,config,dx,dy)
    use glimmer_config
    use glimmer_log
    use glimmer_global, only: dp
    implicit none
    type(glimmap_proj),intent(inout) :: proj !< The projection parameters to be initialised
    type(ConfigSection), pointer :: config   !< structure holding sections of configuration file 
    real(dp),intent(in) :: dx                !< grid resolution in x
    real(dp),intent(in) :: dy                !< grid resolution in y

    ! local variables
    type(ConfigSection), pointer :: section
    real(dp) :: lonc,latc,efalse,nfalse,stdp1,stdp2,scale_factor,cpx,cpy
    real(dp),dimension(:),pointer :: std_par
    character(10) :: ptype
    logical :: stdp,scfac
    integer :: ptval,ptold
    logical :: compute_area_factor

    ptype   = ''
    lonc    = 0.d0 ; latc   = 0.d0
    efalse  = 0.d0 ; nfalse = 0.d0
    compute_area_factor = .false.
    std_par => null()
    scale_factor = 0.d0
    stdp1 = 0.d0
    stdp2 = 0.d0

    call GetSection(config,section,'projection')
    if (associated(section)) then
       call GetValue(section,'type',ptype)
       call GetValue(section,'centre_longitude',lonc)
       call GetValue(section,'centre_latitude',latc)
       call GetValue(section,'false_easting',efalse)
       call GetValue(section,'false_northing',nfalse)
       call GetValue(section,'standard_parallel',std_par)
       call GetValue(section,'scale_factor',scale_factor)
       call GetValue(section,'compute_area_factor', compute_area_factor)

       ! Parse the projection type
       if (index(ptype,'LAEA')/=0 .or. index(ptype,'laea')/=0) then
          ptval = GMAP_LAEA
       else if (index(ptype,'AEA')/=0 .or. index(ptype,'aea')/=0) then
          ptval = GMAP_AEA
       else if (index(ptype,'LCC')/=0 .or. index(ptype,'lcc')/=0) then
          ptval = GMAP_LCC
       else if (index(ptype,'STERE')/=0 .or. index(ptype,'stere')/=0) then
          ptval = GMAP_STERE
       else
          call write_log('Unrecognised type in [projection]', &
               GM_FATAL,__FILE__,__LINE__)
       end if
 
       ! Deal with presence or not of standard parallel(s)
       if (associated(std_par)) then
          stdp = .true.
          select case (size(std_par))
          case(1)
             stdp1 = std_par(1) ; stdp2 = std_par(1)
          case(2)
             stdp1 = std_par(1) ; stdp2 = std_par(2)
          case(0)
             stdp=.false.
          case default
             call write_log('More than two Standard parallels given', &
                  GM_FATAL,__FILE__,__LINE__)
          end select
       else
          stdp = .false.
       end if

       ! Deal with scale factor
       if (scale_factor /= 0.d0) then
          scfac = .true.
       else
          scfac = .false.
       end if

    else
       call GetSection(config,section,'GLINT projection')
       if(.not.associated(section)) return
       call write_log('Using [GLINT projection] config section',GM_WARNING)
       call write_log('This config option has been deprecated, and will be removed at some point.',GM_WARNING)
       call write_log('Use [projection] instead',GM_WARNING)
       call GetValue(section,'projection',ptold)
       call GetValue(section,'lonc',lonc)
       call GetValue(section,'latc',latc)
       call GetValue(section,'cpx',cpx)
       call GetValue(section,'cpy',cpy)
       call GetValue(section,'std_parallel',stdp1)
       select case(ptold)
       case(1)
          ptval = GMAP_LAEA
       case(2:4)
          ptval = GMAP_STERE
       case default
          call write_log('Unsupported projection in [GLINT projection] config section',GM_FATAL)
       end select
       efalse = dx*(cpx-1)
       nfalse = dy*(cpy-1)
       if (stdp1 /= 0.d0) then
          stdp2 = stdp1
          stdp = .true.
       else
          stdp = .false.
       end if
       scfac=.false.
    end if


    ! Check for conflict

    if (stdp.and.scfac) then
       call write_log('You cannot specify both a standard parallel and a scale factor.', &
            GM_FATAL,__FILE__,__LINE__)
    end if

    ! Initialise the projection

    if (stdp) then
       call glimmap_proj_define(proj,ptval, &
       lonc,latc,efalse,nfalse, &
       compute_area_factor = compute_area_factor, &
       standard_parallel = stdp1, &
       standard_parallel_2 = stdp2)
    else if (scfac) then
       call glimmap_proj_define(proj,ptval, &
       lonc,latc,efalse,nfalse, &
       compute_area_factor = compute_area_factor, &
       scale_factor_at_proj_origin = scale_factor)
    else
       call glimmap_proj_define(proj,ptval, &
       lonc,latc,efalse,nfalse, &
       compute_area_factor = compute_area_factor)
    end if

  end subroutine glimmap_readconfig

  !-------------------------------------------------------------------------

  !> print projection info to log
  subroutine glimmap_printproj(proj)
    use glimmer_log
    use glimmer_global, only : msg_length

    type(glimmap_proj),intent(in) :: proj !< the projection

    character(len=msg_length) :: message

    call write_log('Projection')
    call write_log('----------')
    if (.not.proj%found) then
       call write_log('No projection found')
       return
    end if

    if (associated(proj%laea)) then

       call write_log('Type: Lambert Azimuthal Equal Area')
       write(message,*)'Longitude of central meridian: ',proj%laea%longitude_of_central_meridian
       call write_log(message)
       write(message,*)'Latitude of projection origin: ',proj%laea%latitude_of_projection_origin
       call write_log(message)
       write(message,*)'False easting:  ',proj%laea%false_easting
       call write_log(message)
       write(message,*)'False northing: ',proj%laea%false_northing
       call write_log(message)

    else if (associated(proj%aea)) then

       call write_log('Type: Albers Equal Area Conic')
       write(message,*)'Longitude of central meridian: ',proj%aea%longitude_of_central_meridian
       call write_log(message)
       write(message,*)'Latitude of projection origin: ',proj%aea%latitude_of_projection_origin
       call write_log(message)
       write(message,*)'False easting:  ',proj%aea%false_easting
       call write_log(message)
       write(message,*)'False northing: ',proj%aea%false_northing
       call write_log(message)
       write(message,*)'Standard parallels: ', &
            proj%aea%standard_parallel(1),proj%aea%standard_parallel(2)
       call write_log(message)

    else if (associated(proj%lcc)) then

       call write_log('Type: Lambert Conformal Conic')
       write(message,*)'Longitude of central meridian: ',proj%lcc%longitude_of_central_meridian
       call write_log(message)
       write(message,*)'Latitude of projection origin: ',proj%lcc%latitude_of_projection_origin
       call write_log(message)
       write(message,*)'False easting:  ',proj%lcc%false_easting
       call write_log(message)
       write(message,*)'False northing: ',proj%lcc%false_northing
       call write_log(message)
       write(message,*)'Standard parallels: ', &
            proj%lcc%standard_parallel(1),proj%lcc%standard_parallel(2)
       call write_log(message)

    else if (associated(proj%stere)) then

       call write_log('Type: Stereographic')
       write(message,*)'Longitude of central meridian: ',proj%stere%longitude_of_central_meridian
       call write_log(message)
       write(message,*)'Latitude of projection origin: ',proj%stere%latitude_of_projection_origin
       call write_log(message)
       write(message,*)'False easting:  ',proj%stere%false_easting
       call write_log(message)
       write(message,*)'False northing: ',proj%stere%false_northing
       call write_log(message)
       write(message,*)'Standard parallel: ',proj%stere%standard_parallel
       call write_log(message)
       write(message,*)'Scale factor: ',proj%stere%scale_factor_at_proj_origin
       call write_log(message)
       write(message,*)'compute_area_factor:',proj%stere%compute_area_factor
    end if

  end subroutine glimmap_printproj

  !-------------------------------------------------------------------------

  !> Defines a projection from scratch, and initialises 
  !! the other elements appropriately.
  subroutine glimmap_proj_define(cfp,ptype, &
       longitude_of_central_meridian, &
       latitude_of_projection_origin, &
       false_easting, &
       false_northing, &
       compute_area_factor, &
       scale_factor_at_proj_origin, &
       standard_parallel, &
       standard_parallel_2)

    use glimmer_log

    type(glimmap_proj),intent(inout) :: cfp                      !< the projection data type
    integer,intent(in) :: ptype                                  !< the projection ID
    real(dp),intent(in) :: longitude_of_central_meridian         !< the longitude of the central meridian
    real(dp),intent(in) :: latitude_of_projection_origin         !< the latitude of the projection origin
    real(dp),intent(in) :: false_easting                         !< false easting
    real(dp),intent(in) :: false_northing                        !< false northing
    logical, optional,intent(in) :: compute_area_factor          !< if true, compute area distortion factors
    real(dp),optional,intent(in) :: scale_factor_at_proj_origin  !< scale factor
    real(dp),optional,intent(in) :: standard_parallel            !< standard parallel 1
    real(dp),optional,intent(in) :: standard_parallel_2          !< standard parallel 2

    if (associated(cfp%laea))  deallocate(cfp%laea)
    if (associated(cfp%aea))   deallocate(cfp%aea)
    if (associated(cfp%lcc))   deallocate(cfp%lcc)
    if (associated(cfp%stere)) deallocate(cfp%stere)

    cfp%found = .true.
    select case(ptype)
    case(GMAP_LAEA)
       allocate(cfp%laea)
       cfp%laea%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%laea%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%laea%false_easting  = false_easting
       cfp%laea%false_northing = false_northing
       call glimmap_laea_init(cfp%laea)
    case(GMAP_AEA)
       allocate(cfp%aea)
       cfp%aea%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%aea%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%aea%false_easting  = false_easting
       cfp%aea%false_northing = false_northing
       if (present(standard_parallel).and.present(standard_parallel_2)) then
          cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel_2 /)
       else if (present(standard_parallel).and..not.present(standard_parallel_2)) then
          cfp%aea%standard_parallel = (/ standard_parallel,standard_parallel /)
       else
          call write_log('Albers Equal Area: you must supply at least one standard parallel',&
               GM_FATAL,__FILE__,__LINE__)
       end if
       call glimmap_aea_init(cfp%aea)
    case(GMAP_LCC)
       allocate(cfp%lcc)
       cfp%lcc%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%lcc%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%lcc%false_easting  = false_easting
       cfp%lcc%false_northing = false_northing
       if (present(standard_parallel).and.present(standard_parallel_2)) then
          cfp%lcc%standard_parallel = (/ standard_parallel,standard_parallel_2 /)
       else if (present(standard_parallel).and..not.present(standard_parallel_2)) then
          cfp%lcc%standard_parallel = (/ standard_parallel,standard_parallel /)
       else
          call write_log('Lambert Conformal Conic: you must supply at least one standard parallel',&
               GM_FATAL,__FILE__,__LINE__)
       end if
       call glimmap_lcc_init(cfp%lcc)
    case(GMAP_STERE)
       allocate(cfp%stere)
       cfp%stere%longitude_of_central_meridian = longitude_of_central_meridian
       cfp%stere%latitude_of_projection_origin = latitude_of_projection_origin
       cfp%stere%false_easting  = false_easting
       cfp%stere%false_northing = false_northing
       if(present(scale_factor_at_proj_origin) .and. present(standard_parallel)) then
          if (scale_factor_at_proj_origin/=0.d0 .and. standard_parallel/=0.d0) &
               call write_log('Both standard parallel and scale factor specified', &
               GM_FATAL,__FILE__,__LINE__)
       end if
       if(present(scale_factor_at_proj_origin)) &
            cfp%stere%scale_factor_at_proj_origin = scale_factor_at_proj_origin
       if(present(standard_parallel)) &
            cfp%stere%standard_parallel = standard_parallel
       if(present(compute_area_factor)) &
            cfp%stere%compute_area_factor = compute_area_factor
       call glimmap_stere_init(cfp%stere)
    case default
       call write_log('Unrecognised projection type', &
            GM_FATAL,__FILE__,__LINE__)
    end select

  end subroutine glimmap_proj_define

  !> initialise Lambert azimuthal equal area projection
  subroutine glimmap_laea_init(params)

    type(proj_laea),intent(inout) :: params

    params%sinp=sin(params%latitude_of_projection_origin*D2R)
    params%cosp=cos(params%latitude_of_projection_origin*D2R)

    ! Check whether polar

    if (abs(params%latitude_of_projection_origin-90.d0)<CONV_LIMIT) then
       params%pole = 1
    else if (abs(params%latitude_of_projection_origin+90.d0)<CONV_LIMIT) then
       params%pole = -1
    else
       params%pole = 0
    end if

  end subroutine glimmap_laea_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> initialise Lambert azimuthal equal area projection
  subroutine glimmap_aea_init(params)

    type(proj_aea),intent(inout) :: params

    params%n = 0.5d0*(sin(params%standard_parallel(1)*D2R) &
         + sin(params%standard_parallel(2)*D2R))
    params%i_n = 1.d0/params%n
    params%c = cos(params%standard_parallel(1)*D2R)**2.d0 &
         + 2.d0*params%n*sin(params%standard_parallel(1)*D2R)
    params%rho0_R = params%i_n * sqrt(params%c - &
           2.d0*params%n*sin(params%latitude_of_projection_origin*D2R))
    params%rho0 = params%rho0_R * EQ_RAD

  end subroutine glimmap_aea_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> initialise Lambert conformal conic projection
  subroutine glimmap_lcc_init(params)

    type(proj_lcc),intent(inout) :: params

    if (abs(params%standard_parallel(1)-params%standard_parallel(2))<CONV_LIMIT) then
       params%n = sin(params%standard_parallel(1)*D2R)
    else
       params%n = log(cos(params%standard_parallel(1)*D2R)/cos(params%standard_parallel(2)*D2R))/ &
            log(tan(M_PI_4+params%standard_parallel(2)*D2R/2.d0)/tan(M_PI_4+params%standard_parallel(1)*D2R/2.d0))
    end if

    params%i_n = 1.d0/params%n

    params%f = params%i_n*cos(params%standard_parallel(1)*D2R)* &
         (tan(M_PI_4+params%standard_parallel(1)*D2R/2.d0))**params%n
    params%rho0 = EQ_RAD * params%f/(tan(M_PI_4+params%latitude_of_projection_origin*D2R/2.d0))**params%n

  end subroutine glimmap_lcc_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> initialise stereographic projection
  subroutine glimmap_stere_init(params)

    use glimmer_log

    type(proj_stere),intent(inout) :: params

    ! Determine polar/equatorial, etc.

    if (abs(params%latitude_of_projection_origin-90.d0) < CONV_LIMIT) then
       params%pole = 1
    else if (abs(params%latitude_of_projection_origin+90.d0) < CONV_LIMIT) then
       params%pole = -1
    else
       params%pole = 0
       if (abs(params%latitude_of_projection_origin) < CONV_LIMIT) then
          params%equatorial = .true.
       else
          params%equatorial = .false.
       end if
    end if

    ! Set up constants accordingly

    if (params%pole==1 .or. params%pole==-1) then
       if (params%standard_parallel /= 0.d0) then
          if (params%pole==1)  params%k0 = EQ_RAD * (1.d0 + sin(D2R*params%standard_parallel))/2.d0
          if (params%pole==-1) params%k0 = EQ_RAD * (1.d0 - sin(D2R*params%standard_parallel))/2.d0
       else if (params%scale_factor_at_proj_origin /= 0.d0) then
          params%k0 = EQ_RAD * params%scale_factor_at_proj_origin
       else
          params%k0 = EQ_RAD
       end if
    else
       if (params%scale_factor_at_proj_origin /= 0.d0) then
          params%k0 = EQ_RAD * params%scale_factor_at_proj_origin
       else
          params%k0 = EQ_RAD
       end if
       if (params%standard_parallel /= 0.d0) &
            call write_log('Stereographic projection not polar: ignoring standard parallel',GM_WARNING)
       params%sinp = sin(D2R * params%latitude_of_projection_origin)
       params%cosp = cos(D2R * params%latitude_of_projection_origin)
    end if

    params%ik0 = 1.d0/params%k0

  end subroutine glimmap_stere_init

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> compute local area scale factors for stereographic projection
  subroutine glimmap_stere_area_factor(params, ewn, nsn, dx, dy)

    ! Compute area scale factors for each grid cell.
    ! These scale factors describe the distortion of areas in a stereographic projection.
    !
    ! This code is adapted a Matlab script provided by Heiko Goelzer, based on this reference:
    ! J. P. Snyder (1987): Map Projections--A Working Manual, US Geological Survey Professional Paper 1395.
    !
    ! Note: This subroutine should not be called until the input file has been read in,
    !       and we have the relevant grid info (ewn, nsn, dx, dy).

    use glimmer_log
    use glimmer_physcon, only: pi, rearth
    use parallel, only: parallel_globalindex, parallel_reduce_max, parallel_reduce_min

    type(proj_stere),intent(inout) :: params

    integer, intent(in) :: ewn, nsn  ! grid dimensions in x and y
    real(dp), intent(in) :: dx, dy   ! grid cell size in x and y (m)

    ! Local variables

    real(dp) :: lat_c, lon_c        ! latitude and longitude of projection origin (degrees)
    real(dp) :: delta_x, delta_y    ! displacement of lower left corner relative to projection origin (m)
    real(dp) :: deg2rad             ! conversion factor for degrees to radians
    real(dp) :: f                   ! flattening of ellipse
    real(dp) :: ecc                 ! eccentricity of ellipse

    real(dp) :: lambda_0            ! longitude of projection origin (radians)
    real(dp) :: phi_c               ! latitude of projection origin (radians)
    real(dp) :: x, y                ! x and y relative to projection origin (m)
    real(dp) :: rho                 ! distance from projection origin (m)
    real(dp) :: m_c, t_c, t         ! coefficients in Snyder formulas
    real(dp) :: phi, xi, m          ! more coefficients
    real(dp) :: max_area_factor, min_area_factor  ! diagnostic info

    integer :: i, j                 ! local horizontal grid indices
    integer :: iglobal, jglobal     ! global horizontal grid indices

    character(len=100) :: message

    ! allocate the area_factor array and set it to a sensible default
    ! NOTE: If model%projection%stere is not allocated, then area_factor will not be allocated.
    !       In this case, area_factor should not be listed as an output variable in the config file; 
    !        else the code will segfault in glide_io_create.

    allocate(params%area_factor(ewn,nsn))

    if (params%compute_area_factor) then  ! compute area factors

       ! latitude and longitude of projection origin
       lat_c = params%standard_parallel
       lon_c = params%longitude_of_central_meridian

       ! lower left corner of grid is located at (delta_x,delta_y) relative to projection origin
       delta_x = -params%false_easting
       delta_y = -params%false_northing

       ! compute some factors needed for the area_factor calculation

       deg2rad = pi/180.d0   ! degrees to radians

       f = 1.0d0 / 298.257d0        ! flatness parameter (Snyder p. 12)
       ecc = sqrt(2.d0*f - f**2)    ! eccentricity

!!    lambda_0 = deg2rad * lon_c   !WHL - lambda0 is not needed to compute area_factor
       phi_c = deg2rad * lat_c

       m_c = cos(phi_c) / (1.d0 - ecc**2 * (sin(phi_c))**2)**0.5d0   ! Eq. 14-15
       t_c = tan(pi/4.d0 - phi_c/2.d0) / ( ((1.d0 - ecc*sin(phi_c)) / (1.d0 + ecc*sin(phi_c)))**(ecc/2.d0) )  ! Eq. 15-9

       ! compute area_factor for each grid cell
       do j = 1, nsn
          do i = 1, ewn

             call parallel_globalindex(i, j, iglobal, jglobal)

             ! compute x and y at cell center, relative to the projection origin
             x = delta_x + (real(iglobal,dp)-0.5d0)*dx
             y = delta_y + (real(jglobal,dp)-0.5d0)*dy

             ! compute other coefficients from Snyder
             rho = sqrt(x**2 + y**2)       ! Eq. 20-18

             if (abs(lat_c - 90.d0) < 1.d-6) then
                write(message,*) 'WARNING: Area factors may not be correct for lon_c = 90 deg; testing needed'
                call write_log(trim(message))
                ! Note: Snyder defines k0 ~ 1 and includes a factor or rearth in the denominator
                !       Glimmer defines k0 ~ rearth
                t = rho * sqrt((1.d0+ecc)**(1.d0+ecc) * (1.d0+ecc)**(1.d0-ecc)) / (2.d0*params%k0)  ! Eq. 21-39
             else  ! lat_c /= 90.d0) then
                t = rho * t_c/(rearth*m_c)    ! Eq. 21-40
             endif

             xi = pi/2.d0 - 2.d0*atan(t)   ! Eq. 7-13

             phi = (xi + (ecc**2/2.d0 + 5.d0*ecc**4/24.d0 + ecc**6/12.d0 + 13.d0*ecc**8/360.d0) * sin(2.d0*xi) + &
                  (7.d0*ecc**4/48.d0 + 29.d0*ecc**6/120.d0 + 811.d0*ecc**8/11520.d0) * sin(4.d0*xi) + &
                  (7.d0*ecc**6/120.d0 + 81.d0*ecc**8/1120.d0) * sin(6.d0*xi) + &
                  (4279.d0*ecc**8/161280.d0) * sin(8.d0*xi))     ! Eq. 3-5

             m = cos(phi) / (1.d0 - (ecc**2)*(sin(phi))**2)**0.5d0   ! Eq. 14-15

             params%area_factor(i,j) = rho/(rearth*m)      ! Eq. 21-32 (k = area_factor)

          enddo
       enddo

       call write_log ('Computed area scale factors for polar stereographic projection')
       max_area_factor = parallel_reduce_max(maxval(params%area_factor))
       min_area_factor = parallel_reduce_min(minval(params%area_factor))
       write(message,*) 'max, min area_factor:', max_area_factor, min_area_factor
       call write_log(trim(message))

    else  ! compute_area_factor = F

       ! set to a sensible default; assume no area distortion
       params%area_factor(:,:) = 1.0d0
       call write_log ('Set area scale factor = 1 for polar stereographic projection')

    endif  ! compute_area_factor


  end subroutine glimmap_stere_area_factor

end module glimmer_map_init
