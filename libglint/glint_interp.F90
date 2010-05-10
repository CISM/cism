! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glint_interp.f90 - part of the GLIMMER ice model         + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004 GLIMMER contributors - see COPYRIGHT file 
! for list of contributors.
!
! This program is free software; you can redistribute it and/or 
! modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation; either version 2 of 
! the License, or (at your option) any later version.
!
! This program is distributed in the hope that it will be useful, 
! but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License 
! along with this program; if not, write to the Free Software 
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 
! 02111-1307 USA
!
! GLIMMER is maintained by:
!
! Ian Rutt
! School of Geographical Sciences
! University of Bristol
! University Road
! Bristol
! BS8 1SS
! UK
!
! email: <i.c.rutt@bristol.ac.uk> or <ian.rutt@physics.org>
!
! GLIMMER is hosted on berliOS.de:
!
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_interp

  !*FD Downscaling and upscaling routines for use in GLIMMER.

  use glimmer_global
  use glimmer_map_types
  use glint_mpinterp

#ifdef GLC_DEBUG
  use glimmer_paramets, only: stdout, itest, jtest, jjtest, itest_local, jtest_local
#endif

  implicit none

  type downscale

     !*FD Derived type containing indexing 
     !*FD information for downscaling. This type was 
     !*FD included for speed. Four of the arrays contained in it
     !*FD are arrays of the indices of the corners 
     !*FD of the global grid-boxes within which the given
     !*FD local grid point lies.

     real(rk),dimension(:,:),pointer :: llats => null() !*FD The latitude of each point in x-y space.
     real(rk),dimension(:,:),pointer :: llons => null() !*FD The longitude of each point in x-y space.

     integer, dimension(:,:,:),pointer :: xloc => null() !*FD The x-locations of the corner points of the
                                                         !*FD interpolation domain.
     integer, dimension(:,:,:),pointer :: yloc => null() !*FD The y-locations of the corner points of the
                                                         !*FD interpolation domain.
     integer, dimension(:,:), pointer :: xin => null() !*FD x-locations of global cell the point is in
     integer, dimension(:,:), pointer :: yin => null() !*FD y-locations of global cell the point is in

     real(rk),dimension(:,:),  pointer :: xfrac => null()
     real(rk),dimension(:,:),  pointer :: yfrac => null()
     real(rk),dimension(:,:),pointer :: sintheta => NULL()  !*FD sines of grid angle relative to north.
     real(rk),dimension(:,:),pointer :: costheta => NULL()  !*FD coses of grid angle relative to north.
     type(mpinterp) :: mpint !*FD Parameters for mean-preserving interpolation
     logical :: use_mpint = .false. !*FD set true if we're using mean-preserving interpolation
     integer,dimension(:,:),pointer :: lmask => null()  !*FD mask = 1 where downscaling is valid 
                                                        !*FD mask = 0 elsewhere

  end type downscale

  type upscale

     !*FD Derived type containing indexing information
     !*FD for upscaling by areal averaging.

     integer, dimension(:,:),pointer :: gboxx => null() !*FD $x$-indicies of global grid-box 
                                                        !*FD containing given local grid-box.
     integer, dimension(:,:),pointer :: gboxy => null() !*FD $y$-indicies of global grid-box 
                                                        !*FD containing given local grid-box.
     integer, dimension(:,:),pointer :: gboxn => null() !*FD Number of local grid-boxes 
                                                        !*FD contained in each global box.
     logical                         :: set = .false.   !*FD Set if the type has been initialised.
  end type upscale

  interface mean_to_global
     module procedure mean_to_global_sp,mean_to_global_dp
  end interface

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLINT_INTERP
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLINT_INTERP
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLINT_INTERP
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLINT_INTERP
!MH!#endif

  subroutine new_downscale(downs,proj,ggrid,lgrid,mpint)

    use glint_global_grid
    use glimmer_map_trans
    use glimmer_map_types
    use glimmer_coordinates

    !*FD Initialises a downscale variable,
    !*FD according to given projected and global grids

    ! Arguments

    type(downscale),intent(out)      :: downs   !*FD Downscaling variable to be set
    type(glimmap_proj),intent(in)    :: proj    !*FD Projection to use
    type(global_grid),intent(in)     :: ggrid   !*FD Global grid to use
    type(coordsystem_type),intent(in) :: lgrid  !*FD Local (ice) grid 
    logical,optional :: mpint !*FD Set true if we're using mean-preserving interp

    ! Internal variables

    real(rk) :: llat,llon
    integer :: i,j
    type(upscale) :: ups
    integer,dimension(:,:),pointer :: upsm

    upsm => null()
    ! Allocate arrays

    allocate(downs%xloc(lgrid%size%pt(1),lgrid%size%pt(2),4))
    allocate(downs%yloc(lgrid%size%pt(1),lgrid%size%pt(2),4))
    call coordsystem_allocate(lgrid,downs%xfrac)
    call coordsystem_allocate(lgrid,downs%yfrac)
    call coordsystem_allocate(lgrid,downs%llons)
    call coordsystem_allocate(lgrid,downs%llats)
    call coordsystem_allocate(lgrid,downs%sintheta)
    call coordsystem_allocate(lgrid,downs%costheta)
    call coordsystem_allocate(lgrid,downs%xin)
    call coordsystem_allocate(lgrid,downs%yin)
    call coordsystem_allocate(lgrid,upsm)
    call coordsystem_allocate(lgrid,downs%lmask)

    ! index local boxes

    call index_local_boxes(downs%xloc,  downs%yloc,   &
                           downs%xfrac, downs%yfrac,  &
                           ggrid, proj, lgrid,        &
                           downs%lmask )

    ! Calculate grid angle

    call calc_grid_angle(downs,proj,lgrid)

    ! Find lats and lons

    do i=1,lgrid%size%pt(1)
       do j=1,lgrid%size%pt(2)
          call glimmap_xy_to_ll(llon,llat,real(i,rk),real(j,rk),proj,lgrid)
          downs%llons(i,j)=llon
          downs%llats(i,j)=llat
       end do
    end do

    ! Initialise mean-preserving interpolation if necessary
    if (present(mpint)) then
       if (mpint) then
          call new_mpinterp(downs%mpint,ggrid)
          downs%use_mpint = .true.
       end if
    end if

    call new_upscale(ups,ggrid,proj,upsm,lgrid)
    downs%xin = ups%gboxx
    downs%yin = ups%gboxy
    deallocate(upsm)

  end subroutine new_downscale

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_wind_to_local(lgrid,zonwind,merwind,downs,xwind,ywind)

    !*FD Interpolates a global wind field 
    !*FD (or any vector field) onto a given projected grid.

    use glimmer_utils
    use glimmer_coordinates

    ! Argument declarations

    type(coordsystem_type), intent(in)     :: lgrid            !*FD Target grid
    real(rk),dimension(:,:),intent(in)     :: zonwind          !*FD Zonal component (input)
    real(rk),dimension(:,:),intent(in)     :: merwind          !*FD Meridional components (input)
    type(downscale),        intent(inout)  :: downs            !*FD Downscaling parameters
    real(rk),dimension(:,:),intent(out)    :: xwind,ywind      !*FD x and y components on the projected grid (output)

    ! Declare two temporary arrays to hold the interpolated zonal and meridional winds

    real(dp),dimension(size(xwind,1),size(xwind,2)) :: tempzw,tempmw

    ! Check input arrays are conformal to one another

    call check_conformal(zonwind,merwind,'interp_wind 1')
    call check_conformal(xwind,ywind,'interp_wind 2')

    ! Interpolate onto the projected grid

    call interp_to_local(lgrid,zonwind,downs,localdp=tempzw)
    call interp_to_local(lgrid,merwind,downs,localdp=tempmw)

    ! Apply rotation

    xwind=tempzw*downs%costheta-tempmw*downs%sintheta
    ywind=tempzw*downs%sintheta+tempmw*downs%costheta

  end subroutine interp_wind_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine interp_to_local (lgrid,     global,      downs,   &
                              localsp,   localdp,     localrk, &
                              global_fn, z_constrain,          &
                              gmask,     maskval)

    !*FD Interpolate a global scalar field onto a projected grid. 
    !*FD 
    !*FD This uses a simple bilinear interpolation, which assumes
    !*FD that the global grid boxes are rectangular - i.e. it works
    !*FD in lat-lon space.
    !*FD
    !*FD Either localsp or localdp must be present (or both), depending
    !*FD which precision output is required.

! Cell indexing for (xloc,yloc) is as follows:
!
!       4---------3
!       |         |
!       |         |
!       |         |
!       1---------2
!

    use glimmer_utils
    use glimmer_coordinates
    use glimmer_log

    ! Argument declarations

    type(coordsystem_type),  intent(in)           :: lgrid     !*FD Local grid
    real(rk), dimension(:,:),intent(in)           :: global    !*FD Global field (input)
    type(downscale),         intent(inout)        :: downs     !*FD Downscaling parameters
    real(sp),dimension(:,:), intent(out),optional :: localsp   !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:), intent(out),optional :: localdp   !*FD Local field on projected grid (output) dp
    real(rk),dimension(:,:), intent(out),optional :: localrk   !*FD Local field on projected grid (output) rk
    real(sp),optional,external                    :: global_fn !*FD Function returning values in global field. This  
                                                               !*FD may be used as an alternative to passing the
                                                               !*FD whole array in \texttt{global} if, for instance the
                                                               !*FD data-set is in a large file, being accessed point by point.
                                                               !*FD In these circumstances, \texttt{global}
                                                               !*FD may be of any size, and its contents are irrelevant.
    logical,optional :: z_constrain
    integer, dimension(:,:), intent(in),optional  :: gmask     !*FD = 1 where global data are valid, else = 0
    real(dp), intent(in), optional                :: maskval   !*FD Value to write for masked-out cells 

    ! Local variable declarations

    integer  :: i,j                          ! Counter variables for main loop
    real(rk),dimension(4) :: f               ! Temporary array holding the four points in the 
                                             ! interpolation domain.
    real(rk), dimension(size(global,1),size(global,2)) :: g_loc
    logical,  dimension(size(global,1),size(global,2)) :: zeros
    logical :: zc

    integer :: x1, x2, x3, x4
    integer :: y1, y2, y3, y4

#ifdef GLC_DEBUG
    integer :: n
#endif

    if (present(z_constrain)) then
       zc=z_constrain
    else
       zc=.false.
    end if

    ! check we have one output at least...

    if (.not.(present(localsp).or.present(localdp).or.present(localrk))) then
       call write_log('Interp_to_local has no output',GM_WARNING,__FILE__,__LINE__)
    endif

    ! Do stuff for mean-preserving interpolation

    if (downs%use_mpint) then
       call mean_preserve_interp(downs%mpint,global,g_loc,zeros)
    end if

    ! Main interpolation loop

    do i=1,lgrid%size%pt(1)
       do j=1,lgrid%size%pt(2)

#ifdef GLC_DEBUG
!          if (i==itest_local .and. j==jtest_local) then
!              write(stdout,*) ' '
!              write(stdout,*) 'Interpolating, i, j, lmask = ', i, j, downs%lmask(i,j)
!              write(stdout,*) 'xloc, yloc:'
!              do n = 1, 4
!                 write(stdout,*) downs%xloc(i,j,n), downs%yloc(i,j,n)
!              enddo
!          endif
#endif

          ! Compile the temporary array f from adjacent points 

!lipscomb - to do - This could be handled more efficiently by precomputing arrays that
!  specify which neighbor gridcell supplies values in each masked-out global gridcell.
 
          if (present(gmask) .and. present(maskval)) then

             if (downs%lmask(i,j) == 0) then

                f(1) = maskval
                f(2) = maskval
                f(3) = maskval
                f(4) = maskval

             else

                x1 = downs%xloc(i,j,1); y1 = downs%yloc(i,j,1)
                x2 = downs%xloc(i,j,2); y2 = downs%yloc(i,j,2)
                x3 = downs%xloc(i,j,3); y3 = downs%yloc(i,j,3)
                x4 = downs%xloc(i,j,4); y4 = downs%yloc(i,j,4)

                ! if a gridcell is masked out, try to assign a value from a
                ! neighbor that is not masked out

                if (gmask(x1,y1) /= 0) then
                   f(1) = global(x1,y1)
                elseif (gmask(x2,y2) /= 0) then
                   f(1) = global(x2,y2)
                elseif (gmask(x4,y4) /= 0) then
                   f(1) = global(x4,y4)
                elseif  (gmask(x3,y3) /= 0) then
                   f(1) = global(x3,y3)
                else
                   f(1) = maskval
                endif

                if (gmask(x2,y2) /= 0) then
                   f(2) = global(x2,y2)
                elseif (gmask(x1,y1) /= 0) then
                   f(2) = global(x1,y1)
                elseif (gmask(x3,y3) /= 0) then
                   f(2) = global(x3,y3)
                elseif  (gmask(x4,y4) /= 0) then
                   f(2) = global(x4,y4)
                else
                   f(2) = maskval
                endif

                if (gmask(x3,y3) /= 0) then
                   f(3) = global(x3,y3)
                elseif (gmask(x4,y4) /= 0) then
                   f(3) = global(x4,y4)
                elseif (gmask(x2,y2) /= 0) then
                   f(3) = global(x2,y2)
                elseif  (gmask(x1,y1) /= 0) then
                   f(3) = global(x1,y1)
                else
                   f(3) = maskval
                endif

                if (gmask(x4,y4) /= 0) then
                   f(4) = global(x4,y4)
                elseif (gmask(x3,y3) /= 0) then
                   f(4) = global(x3,y3)
                elseif (gmask(x1,y1) /= 0) then
                   f(4) = global(x1,y1)
                elseif  (gmask(x2,y2) /= 0) then
                   f(4) = global(x2,y2)
                else
                   f(4) = maskval
                endif

             endif    ! lmask = 0

          else        ! gmask and maskval not present

             if (present(global_fn)) then
                f(1)=global_fn(downs%xloc(i,j,1),downs%yloc(i,j,1))
                f(2)=global_fn(downs%xloc(i,j,2),downs%yloc(i,j,2))
                f(3)=global_fn(downs%xloc(i,j,3),downs%yloc(i,j,3))
                f(4)=global_fn(downs%xloc(i,j,4),downs%yloc(i,j,4))
             else
                if (downs%use_mpint) then
                   f(1)=g_loc(downs%xloc(i,j,1),downs%yloc(i,j,1))
                   f(2)=g_loc(downs%xloc(i,j,2),downs%yloc(i,j,2))
                   f(3)=g_loc(downs%xloc(i,j,3),downs%yloc(i,j,3))
                   f(4)=g_loc(downs%xloc(i,j,4),downs%yloc(i,j,4))
                else
                   f(1)=global(downs%xloc(i,j,1),downs%yloc(i,j,1))
                   f(2)=global(downs%xloc(i,j,2),downs%yloc(i,j,2))
                   f(3)=global(downs%xloc(i,j,3),downs%yloc(i,j,3))
                   f(4)=global(downs%xloc(i,j,4),downs%yloc(i,j,4))
                end if
             end if

          endif  ! gmask and maskval present

          ! Apply the bilinear interpolation

          if (zc.and.zeros(downs%xin(i,j),downs%yin(i,j)).and.downs%use_mpint) then
             if (present(localsp)) localsp(i,j)=0.0_sp
             if (present(localdp)) localdp(i,j)=0.0_dp
             if (present(localrk)) localrk(i,j)=0.0_rk
          else
             if (present(localsp)) localsp(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)
             if (present(localdp)) localdp(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)
             if (present(localrk)) localrk(i,j)=bilinear_interp(downs%xfrac(i,j),downs%yfrac(i,j),f)
          end if

       enddo
    enddo

  end subroutine interp_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_local(proj,lgrid,ggrid,global,localsp,localdp,global_fn)

    !*FD Average a high-resolution global field onto the projected grid
    !*FD This assumes that the global field is sufficiently high-resolution 
    !*FD compared with the local grid - it just averages the points contained 
    !*FD in each local grid-box.

    use glimmer_map_types
    use glimmer_map_trans
    use glimmer_coordinates
    use glimmer_utils
    use glimmer_log
    use glint_global_grid

    ! Argument declarations

    type(glimmap_proj),              intent(in)  :: proj      !*FD Target map projection
    type(coordsystem_type),          intent(in)  :: lgrid     !*FD Local grid information
    type(global_grid),               intent(in)  :: ggrid     !*FD Global grid information
    real(rk),dimension(:,:),         intent(in)  :: global    !*FD Global field (input)
    real(sp),dimension(:,:),optional,intent(out) :: localsp   !*FD Local field on projected grid (output) sp
    real(dp),dimension(:,:),optional,intent(out) :: localdp   !*FD Local field on projected grid (output) dp
    real(sp),optional, external                  :: global_fn !*FD Function returning values in global field. This  
                                                              !*FD may be used as an alternative to passing the
                                                              !*FD whole array in \texttt{global} if, for instance the
                                                              !*FD data-set is in a large file, being accessed point by point.
                                                              !*FD In these circumstances, \texttt{global}
                                                              !*FD may be of any size, and its contents are irrelevant.

    integer :: i,j,xbox,ybox
    real(rk) :: lat,lon,x,y
    real(dp),dimension(lgrid%size%pt(1),lgrid%size%pt(2)) :: temp_out
    real(rk),dimension(lgrid%size%pt(1),lgrid%size%pt(2)) :: mean_count

    if (.not.present(global_fn)) then 
       if ((lgrid%size%pt(1)/=size(ggrid%lons)).or.(lgrid%size%pt(2)/=size(ggrid%lats))) then
          call write_log('Size mismatch in interp_to_local',GM_FATAL,__FILE__,__LINE__)
       end if
    end if

    ! check we have one output at least...

    if (.not.(present(localsp).or.present(localdp))) then
       call write_log('mean_to_local has no output',GM_WARNING,__FILE__,__LINE__)
    endif

    ! Zero some things

    mean_count=0
    temp_out=0.0

    ! Loop over all global points

    do i=1,lgrid%size%pt(1)

       lon=ggrid%lons(i)

       do j=1,lgrid%size%pt(2)

          ! Find location in local coordinates

          lat=ggrid%lats(j)  ! (Have already found lat above)
          call glimmap_ll_to_xy(lon,lat,x,y,proj,lgrid)
          xbox=nint(x)
          ybox=nint(y)

          ! Add to appropriate location and update count

          if (xbox.ge.1.and.xbox.le.lgrid%size%pt(1).and. &
               ybox.ge.1.and.ybox.le.lgrid%size%pt(2)) then
             if (present(global_fn)) then
                temp_out(xbox,ybox)=temp_out(xbox,ybox)+global_fn(i,j)*ggrid%box_areas(xbox,ybox)
             else
                temp_out(xbox,ybox)=temp_out(xbox,ybox)+global(i,j)*ggrid%box_areas(xbox,ybox)
             end if
             mean_count(xbox,ybox)=mean_count(xbox,ybox)+ggrid%box_areas(xbox,ybox)
          end if

       end do
    end do

    ! Divide by number of contributing points and copy to output

    if (present(localsp)) localsp=temp_out/real(mean_count,dp)
    if (present(localdp)) localdp=temp_out/real(mean_count,dp)

  end subroutine mean_to_local

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine pointwise_to_global(proj,lgrid,local,lons,lats,global)

    !*FD Upscale to global domain by
    !*FD pointwise sampling.
    !*FD
    !*FD Note that this is the mathematically inverse process of the 
    !*FD \texttt{interp\_to\_local} routine.

    use glimmer_coordinates
    use glimmer_map_trans

    ! Arguments

    type(glimmap_proj),     intent(in)  :: proj      !*FD Projection to use
    type(coordsystem_type), intent(in)  :: lgrid     !*FD Local grid
    real(rk),dimension(:,:),intent(in)  :: local     !*FD Local field (input)
    real(rk),dimension(:,:),intent(out) :: global    !*FD Global field (output)
    real(rk),dimension(:),  intent(in)  :: lats      !*FD Latitudes of grid-points (degrees)
    real(rk),dimension(:),  intent(in)  :: lons      !*FD Longitudes of grid-points (degrees)

    ! Internal variables

    real(rk),dimension(2,2) :: f
    integer :: nxg,nyg,nxl,nyl,i,j,xx,yy
    real(rk) :: x,y
    real(rk),dimension(size(local,1),size(local,2)) :: tempmask

    nxg=size(global,1) ; nyg=size(global,2)
    nxl=size(local,1)  ; nyl=size(local,2)

    do i=1,nxg
       do j=1,nyg
          call glimmap_ll_to_xy(lons(i),lats(j),x,y,proj,lgrid)
          xx=int(x) ; yy=int(y)
          if (nint(x)<=1.or.nint(x)>nxl-1.or.nint(y)<=1.or.nint(y)>nyl-1) then
             global(i,j)=0.0
          else
             f=local(xx:xx+1,yy:yy+1)*tempmask(xx:xx+1,yy:yy+1)
             global(i,j)=bilinear_interp((x-real(xx))/real(1.0,rk),(y-real(yy))/real(1.0,rk),f)
          endif
       enddo
    enddo

  end subroutine pointwise_to_global

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global_sp(ups,local,global,mask)

    !*FD Upscale to global domain by
    !*FD areal averaging.
    !*FD
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item \texttt{gboxx} and \texttt{gboxy} are the same size as \texttt{local}
    !*FD \item \texttt{gboxn} is the same size as \texttt{global}
    !*FD \item This method is \emph{not} the mathematical inverse of the
    !*FD \texttt{interp\_to\_local} routine.
    !*FD \end{itemize}

    ! Arguments

    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(sp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(rk),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer, dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl,nyl,i,j
    real(rk),dimension(size(local,1),size(local,2)) :: tempmask

    ! Beginning of code

    nxl=size(local,1) ; nyl=size(local,2)

    global=0.0

    if (present(mask)) then
       tempmask=mask
    else
       tempmask=1
    endif

    do i=1,nxl
       do j=1,nyl
          global(ups%gboxx(i,j),ups%gboxy(i,j))= &
               global(ups%gboxx(i,j),ups%gboxy(i,j))+local(i,j)*tempmask(i,j)
       enddo
    enddo

    where (ups%gboxn.ne.0)
       global=global/ups%gboxn
    elsewhere
       global=0.0
    endwhere

  end subroutine mean_to_global_sp

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global_dp(ups,local,global,mask)

    !*FD Upscale to global domain by
    !*FD areal averaging.
    !*FD
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item \texttt{gboxx} and \texttt{gboxy} are the same size as \texttt{local}
    !*FD \item \texttt{gboxn} is the same size as \texttt{global}
    !*FD \item This method is \emph{not} the mathematical inverse of the
    !*FD \texttt{interp\_to\_local} routine.
    !*FD \end{itemize}

    ! Arguments

    type(upscale),          intent(in)  :: ups    !*FD Upscaling indexing data.
    real(dp),dimension(:,:),intent(in)  :: local  !*FD Data on projected grid (input).
    real(rk),dimension(:,:),intent(out) :: global !*FD Data on global grid (output).
    integer, dimension(:,:),intent(in),optional :: mask !*FD Mask for upscaling

    ! Internal variables

    integer :: nxl,nyl,i,j
    real(rk),dimension(size(local,1),size(local,2)) :: tempmask

    ! Beginning of code

    nxl=size(local,1) ; nyl=size(local,2)

    global=0.0

    if (present(mask)) then
       tempmask=mask
    else
       tempmask=1
    endif

    do i=1,nxl
       do j=1,nyl
          global(ups%gboxx(i,j),ups%gboxy(i,j))= &
               global(ups%gboxx(i,j),ups%gboxy(i,j))+local(i,j)*tempmask(i,j)
       enddo
    enddo

    where (ups%gboxn.ne.0)
       global=global/ups%gboxn
    elsewhere
       global=0.0
    endwhere

  end subroutine mean_to_global_dp

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine mean_to_global_mec(ups,                &
                                nxl,      nyl,      &
                                nxg,      nyg,      &
                                nec,      topomax,  &
                                local,    global,   &
                                ltopo,    mask)
 
    ! Upscale from the local domain to a global domain with multiple elevation classes
    ! by areal averaging.
    !
    ! This subroutine is adapted from subroutine mean_to_global in GLIMMER.
    ! The difference is that local topography is upscaled to multiple elevation classes
    !  in each global grid cell.
    !
    ! Note: This method is not the inverse of the interp_to_local routine.
    ! Also note that each local grid cell is weighted equally.
    ! In the future we will probably want to use the CCSM coupler for upscaling.
 
    use glimmer_log

    ! Arguments
 
    type(upscale),            intent(in)    :: ups     ! upscaling indexing data
    integer,                  intent(in)    :: nxl,nyl ! local grid dimensions 
    integer,                  intent(in)    :: nxg,nyg ! global grid dimensions 
    integer,                  intent(in)    :: nec     ! number of elevation classes 
    real(dp),dimension(0:nec),intent(in)    :: topomax ! max elevation in each class 
    real(dp),dimension(nxl,nyl),  intent(in)      :: local   ! data on local grid
    real(dp),dimension(nxg,nyg,nec),intent(out)   :: global  ! data on global grid
    real(dp),dimension(nxl,nyl),  intent(in)      :: ltopo   ! surface elevation on local grid (m)
    integer, dimension(nxl,nyl),intent(in),optional :: mask ! mask for upscaling

    ! Internal variables
 
    integer ::  &
       i, j, n,    &! indices
       ig, jg       ! indices

    integer, dimension(nxl,nyl) ::  &
        tempmask,    &! temporary mask
        gboxec        ! elevation class into which local data is upscaled

    integer, dimension(nxg,nyg,nec) ::  &
        gnumloc       ! no. of local cells within each global cell in each elevation class


    integer :: il, jl
    real(dp) :: lsum, gsum
 
    if (present(mask)) then
       tempmask(:,:) = mask(:,:)
    else
       tempmask(:,:) = 1
    endif
 
    ! Compute global elevation class for each local grid cell
    ! Also compute number of local cells within each global cell in each elevation class

    gboxec(:,:) = 0
    gnumloc(:,:,:) = 0

    do n = 1, nec
       do j = 1, nyl
       do i = 1, nxl
          if (ltopo(i,j) >= topomax(n-1) .and. ltopo(i,j) < topomax(n)) then
             gboxec(i,j) = n
             if (tempmask(i,j)==1) then
                ig = ups%gboxx(i,j)
                jg = ups%gboxy(i,j)
                gnumloc(ig,jg,n) = gnumloc(ig,jg,n) + 1
             endif
          endif
       enddo
       enddo
    enddo

    global(:,:,:) = 0._dp


    do j = 1, nyl
    do i = 1, nxl
       ig = ups%gboxx(i,j)
       jg = ups%gboxy(i,j)
       n = gboxec(i,j)
       if (n==0) then
#ifdef GLC_DEBUG
          write(stdout,*) 'Upscaling error: local topography out of bounds'
          write(stdout,*) 'i, j, topo:', i, j, ltopo(i,j)
          write(stdout,*) 'topomax(0) =', topomax(0)
#endif
          call write_log('Upscaling error: local topography out of bounds', &
               GM_FATAL,__FILE__,__LINE__)
       endif

#ifdef GLC_DEBUG
       if (i==itest_local .and. j==jtest_local) then
          write(stdout,*) ' '
          write(stdout,*) 'il, jl =', i, j
          write(stdout,*) 'ig, jg, n =', ig, jg, n
          write(stdout,*) 'Old global val =', global(ig,jg,n)
          write(stdout,*) 'local, mask =', local(i,j), tempmask(i,j)
       endif
#endif
       global(ig,jg,n) = global(ig,jg,n) + local(i,j)*tempmask(i,j)

#ifdef GLC_DEBUG
       if (i==itest_local .and. j==jtest_local) then
          write(stdout,*) 'New global val =', global(ig,jg,n)
       endif
#endif

    enddo
    enddo
 
    do n = 1, nec
       do j = 1, nyg
       do i = 1, nxg
          if (gnumloc(i,j,n) /= 0) then
             global(i,j,n) = global(i,j,n) / gnumloc(i,j,n)
          else
             global(i,j,n) = 0._dp
          endif
       enddo
       enddo
    enddo

    ! conservation check

    lsum = 0._dp
    do j = 1, nyl
    do i = 1, nxl
       lsum = lsum + local(i,j)*tempmask(i,j)
    enddo    
    enddo

    gsum = 0._dp
    do n = 1, nec
    do j = 1, nyg
    do i = 1, nxg
       gsum = gsum + global(i,j,n)*gnumloc(i,j,n)
    enddo
    enddo
    enddo

    if (abs(lsum) > 1.e-10_dp) then
       if (abs(gsum-lsum)/abs(lsum) > 1.e-10_dp) then 
#ifdef GLC_DEBUG
          write(stdout,*) 'local and global sums disagree'
          write (stdout,*) 'lsum, gsum =', lsum, gsum 
#endif
          call write_log('Upscaling error: local and glocal sums disagree', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    else  ! lsum is close to zero
       if (abs(gsum-lsum) > 1.e-10_dp) then
#ifdef GLC_DEBUG
          write(stdout,*) 'local and global sums disagree'
          write (stdout,*) 'lsum, gsum =', lsum, gsum 
#endif
          call write_log('Upscaling error: local and glocal sums disagree', &
               GM_FATAL,__FILE__,__LINE__)
       endif
    endif

  end subroutine mean_to_global_mec
  
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function bilinear_interp(xp,yp,f)

    !*FD Performs bilinear interpolation 
    !*FD in a rectangular domain. Note that the bilinear interpolation formula is:
    !*FD  \[f_{\mathtt{x},\mathtt{y}}=(1-X')(1-Y')f_{1}+X'(1-Y')f_{2}+X'Y'f_{3}+(1-X')Y'f_{4}\]
    !*FD where $X'$ and $Y'$ are the fractional displacements of the target point within the domain.
    !*RV The value of \texttt{f} at \texttt{x,y}

    ! Argument declarations

    real(rk),             intent(in) :: xp    !*FD The fractional $x$-displacement of the target.
    real(rk),             intent(in) :: yp    !*FD The fractional $y$-displacement of the target.
    real(rk),dimension(4),intent(in) :: f     !*FD The interpolation domain;
                                              !*FD i.e. the four points surrounding the
                                              !*FD target, presented anticlockwise from bottom-
                                              !*FD left 
    ! Apply bilinear interpolation formula

    bilinear_interp=(1-xp)*(1-yp)*f(1)+xp*(1-yp)*f(2)+xp*yp*f(3)+(1-xp)*yp*f(4)

  end function bilinear_interp

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine find_ll_index(il,jl,lon,lat,lons,lats)

    !*FD Find the global gridpoint at the first corner of the box surrounding
    !*FD a given location in lat-lon space.

    use glimmer_utils

    ! Arguments

    real(rk),             intent(in)  :: lon    !*FD Longitude of location to be indexed (input)
    real(rk),             intent(in)  :: lat    !*FD Latitude of location to be indexed (input)
    real(rk),dimension(:),intent(in)  :: lats   !*FD Latitudes of global grid points 
    real(rk),dimension(:),intent(in)  :: lons   !*FD Longitudes of global grid points 
    integer,              intent(out) :: il     !*FD $x$-gridpoint index (output)
    integer,              intent(out) :: jl     !*FD $y$-gridpoint index (output)

    ! Internal variables

    integer :: nx,ny,i

    nx=size(lons) ; ny=size(lats)

    il=nx

    do i=1,nx-1
       if (lon_between(lons(i),lons(i+1),lon)) then
          il=i
          exit
       endif
    enddo

    if ((lat<lats(ny)).and.(lat>-90.0)) then
       jl=ny
       return
    endif

    if ((lat>lats(1)).and.(lat<90.0)) then
       jl=1
       return
    endif

    jl=1
    do 
       if (lat>lats(jl)) exit
       jl=jl+1
    enddo

  end subroutine find_ll_index

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine index_local_boxes (xloc, yloc, xfrac, yfrac, ggrid, proj, lgrid, lmask)

    !*FD Indexes the corners of the
    !*FD global grid box in which each local grid box sits.

    use glimmer_utils
    use glint_global_grid
    use glimmer_coordinates
    use glimmer_map_trans

    ! Arguments

    integer, dimension(:,:,:),intent(out) :: xloc,yloc   !*FD Array of indicies (see \texttt{downscale} type)
    real(rk),dimension(:,:),  intent(out) :: xfrac,yfrac !*FD Fractional off-sets of grid points
    type(global_grid),        intent(in)  :: ggrid       !*FD Global grid to be used
    type(glimmap_proj),       intent(in)  :: proj        !*FD Projection to be used
    type(coordsystem_type),   intent(in)  :: lgrid       !*FD Local grid
    integer, dimension(:,:),  intent(out) :: lmask  !*FD Mask of local cells for which interpolation is valid

    ! Internal variables

    integer :: i,j,il,jl,temp
    real(rk) :: ilon,jlat,xa,ya,xb,yb,xc,yc,xd,yd

#ifdef GLC_DEBUG
    integer :: nx, ny, nxg, nyg, n
    nx = lgrid%size%pt(1)
    ny = lgrid%size%pt(2)
    nxg = size(ggrid%mask,1)
    nyg = size(ggrid%mask,2)

    write(stdout,*) ' '
    write(stdout,*) 'nx,  ny =', nx, ny
    write(stdout,*) 'nxg, nyg =', nxg, nyg
    write(stdout,*) 'Indexing local boxes'
#endif

    do i=1,lgrid%size%pt(1)
       do j=1,lgrid%size%pt(2)

          ! Find out where point i,j is in lat-lon space

          call glimmap_xy_to_ll(ilon,jlat,real(i,rk),real(j,rk),proj,lgrid)

          ! Index that location onto the global grid

          call find_ll_index(il,jl,ilon,jlat,ggrid%lons,ggrid%lats)

          xloc(i,j,1)=il  ! This is the starting point - we now need to find
          yloc(i,j,1)=jl  ! three other points that enclose the interpolation target

          if (jlat>ggrid%lats(ggrid%ny)) then

             ! For all points except on the bottom row

             xloc(i,j,2)=il+1
             yloc(i,j,2)=jl

             xloc(i,j,3)=il+1
             yloc(i,j,3)=jl-1

             xloc(i,j,4)=il
             yloc(i,j,4)=jl-1

             call fix_bcs2d(xloc(i,j,2) ,yloc(i,j,2),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,3) ,yloc(i,j,3),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,4) ,yloc(i,j,4),ggrid%nx,ggrid%ny)

             if (jl==1) then
                temp=xloc(i,j,3)
                xloc(i,j,3)=xloc(i,j,4)
                xloc(i,j,4)=temp
             endif

          else

             ! The bottom row

             xloc(i,j,2)=il-1
             yloc(i,j,2)=jl

             xloc(i,j,3)=il-1
             yloc(i,j,3)=jl+1

             xloc(i,j,4)=il
             yloc(i,j,4)=jl+1

             call fix_bcs2d(xloc(i,j,2) ,yloc(i,j,2),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,3) ,yloc(i,j,3),ggrid%nx,ggrid%ny)
             call fix_bcs2d(xloc(i,j,4) ,yloc(i,j,4),ggrid%nx,ggrid%ny)

             temp=xloc(i,j,3)
             xloc(i,j,3)=xloc(i,j,4)
             xloc(i,j,4)=temp

          endif

#ifdef GLC_DEBUG
          if (i==itest_local .and. j==jtest_local) then
             write(stdout,*) ' '
             write(stdout,*) 'Index local boxes, i, j =', i, j
             write(stdout,*) 'xloc, yloc'
             do n = 1, 4
                write(stdout,*) xloc(i,j,n), yloc(i,j,n)
             enddo
          endif
#endif

          ! Now, find out where each of those points is on the projected
          ! grid, and calculate fractional displacements accordingly

          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,1)),ggrid%lats(yloc(i,j,1)),xa,ya,proj,lgrid)
          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,2)),ggrid%lats(yloc(i,j,2)),xb,yb,proj,lgrid)
          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,3)),ggrid%lats(yloc(i,j,3)),xc,yc,proj,lgrid)
          call glimmap_ll_to_xy(ggrid%lons(xloc(i,j,4)),ggrid%lats(yloc(i,j,4)),xd,yd,proj,lgrid)

          call calc_fractional(xfrac(i,j),yfrac(i,j),real(i,rk),real(j,rk), &
               xa,ya,xb,yb,xc,yc,xd,yd)

          ! If all four global gridcells surrounding this local gridcell
          ! are masked out, then mask out the local gridcell

          if ( (ggrid%mask(xloc(i,j,1),yloc(i,j,1)) == 0) .and.  &
               (ggrid%mask(xloc(i,j,2),yloc(i,j,2)) == 0) .and.  &
               (ggrid%mask(xloc(i,j,3),yloc(i,j,3)) == 0) .and.  &
               (ggrid%mask(xloc(i,j,4),yloc(i,j,4)) == 0) ) then
             lmask(i,j) = 0
          else
             lmask(i,j) = 1
          endif

       enddo
    enddo

#ifdef GLC_DEBUG
          write(stdout,*) ' '
          write(stdout,*) 'Mask in neighborhood of i, j = ', itest_local, jtest_local
          do j = jtest_local-1, jtest_local+1
             write(stdout,*) lmask(itest_local-1:itest_local+1,j)
          enddo

          write(stdout,*) ' '
          write(stdout,*) 'Global mask near Greenland'
          do j = 1, 20
             write(stdout,150) ggrid%mask(nxg-29:nxg,j)
          enddo

          write(stdout,*) ' '
          write(stdout,*) 'Local mask'
          do j = ny, 1, -1
             write(stdout,200) lmask(1:nx,j)
          enddo

  100     format(144i2)
  150     format(30i2)
  200     format(76i2)
#endif

  end subroutine index_local_boxes

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_grid_angle(downs,proj,lgrid)

    !*FD Calculates the angle the projected 
    !*FD grid makes with north at each point and stores the cos 
    !*FD and sin of that angle in the relevant arrays in \texttt{proj}.

    use glimmer_coordinates
    use glimmer_map_trans

    type(downscale),intent(inout) :: downs !*FD The projection to be used
    type(glimmap_proj),intent(in) :: proj
    type(coordsystem_type),intent(in) :: lgrid

    integer :: i,j
    real(rk) :: latn,lonn,lats,lons,lat,lon,dlat,dlon,temp

    do i=1,lgrid%size%pt(1)

       ! Main, central block

       do j=2,lgrid%size%pt(2)-1
          call glimmap_xy_to_ll(lonn,latn,real(i,rk),real(j+1,rk),proj,lgrid)
          call glimmap_xy_to_ll(lon,lat,real(i,rk),real(j,rk),proj,lgrid)
          call glimmap_xy_to_ll(lons,lats,real(i,rk),real(j-1,rk),proj,lgrid)
          dlat=latn-lats
          dlon=lonn-lons
          if (dlon<-90) dlon=dlon+360
          temp=atan(dlon/dlat)
          downs%sintheta(i,j)=sin(temp)
          downs%costheta(i,j)=cos(temp)
       enddo

       ! bottom row

       call glimmap_xy_to_ll(lonn,latn,real(i,rk),real(2,rk),proj,lgrid)
       call glimmap_xy_to_ll(lon,lat,real(i,rk),real(1,rk),proj,lgrid)
       dlat=latn-lat
       dlon=lonn-lon
       if (dlon<-90) dlon=dlon+360
       temp=atan(dlon/dlat)
       downs%sintheta(i,1)=sin(temp)
       downs%costheta(i,1)=cos(temp)

       ! top row

       call glimmap_xy_to_ll(lon,lat,real(i,rk),real(lgrid%size%pt(2),rk),proj,lgrid)
       call glimmap_xy_to_ll(lons,lats,real(i,rk),real(lgrid%size%pt(2)-1,rk),proj,lgrid)
       dlat=lat-lats
       dlon=lon-lons
       if (dlon<-90) dlon=dlon+360
       temp=atan(dlon/dlat)
       downs%sintheta(i,lgrid%size%pt(2))=sin(temp)
       downs%costheta(i,lgrid%size%pt(2))=cos(temp)

    enddo

  end subroutine calc_grid_angle

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine new_upscale(ups,ggrid,proj,mask,lgrid)

    use glint_global_grid
    use glimmer_log
    use glimmer_map_trans
    use glimmer_coordinates

    !*FD Compiles an index of which global grid box contains a given
    !*FD grid box on the projected grid, and sets derived type \texttt{ups}
    !*FD accordingly.

    ! Arguments

    type(upscale),         intent(out) :: ups        !*FD Upscaling type to be set
    type(global_grid),     intent(in)  :: ggrid      !*FD Global grid to be used
    type(glimmap_proj),    intent(in)  :: proj       !*FD Projection being used
    integer,dimension(:,:),intent(in)  :: mask       !*FD Upscaling mask to be used
    type(coordsystem_type),intent(in)  :: lgrid      !*FD local grid

    ! Internal variables

    integer  :: i,j,ii,jj,nx,ny,gnx,gny
    real(rk) :: plon,plat

    ! Beginning of code

    if (associated(ups%gboxx)) deallocate(ups%gboxx)
    if (associated(ups%gboxy)) deallocate(ups%gboxy)
    if (associated(ups%gboxn)) deallocate(ups%gboxn)

    allocate(ups%gboxx(lgrid%size%pt(1),lgrid%size%pt(2)))
    allocate(ups%gboxy(lgrid%size%pt(1),lgrid%size%pt(2)))     
    allocate(ups%gboxn(ggrid%nx,ggrid%ny))

    gnx=ggrid%nx ; gny=ggrid%ny
    nx =lgrid%size%pt(1) ; ny =lgrid%size%pt(2)

    ups%gboxx=0 ; ups%gboxy=0

    do i=1,nx
       do j=1,ny
          call glimmap_xy_to_ll(plon,plat,real(i,rk),real(j,rk),proj,lgrid)
          ii=1 ; jj=1
          do
             ups%gboxx(i,j)=ii
             if (ii>gnx) then
                call write_log('global index failure',GM_FATAL,__FILE__,__LINE__)
             endif
             if (lon_between(ggrid%lon_bound(ii),ggrid%lon_bound(ii+1),plon)) exit
             ii=ii+1
          enddo

          jj=1

          do
             ups%gboxy(i,j)=jj
             if (jj>gny) then
                call write_log('global index failure',GM_FATAL,__FILE__,__LINE__)
             endif
             if ((ggrid%lat_bound(jj)>=plat).and.(plat>ggrid%lat_bound(jj+1))) exit
             jj=jj+1
          enddo

       enddo
    enddo

    ups%gboxn=0

    do i=1,nx
       do j=1,ny
          ups%gboxn(ups%gboxx(i,j),ups%gboxy(i,j))=ups%gboxn(ups%gboxx(i,j),ups%gboxy(i,j))+mask(i,j)
       enddo
    enddo

    ups%set=.true.

  end subroutine new_upscale

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine copy_upscale(in,out)

    use glimmer_log

    type(upscale),intent(in)  :: in
    type(upscale),intent(out) :: out

    if (.not.in%set) then
       call write_log('Attempt to copy un-initialised upscale type',GM_FATAL,&
            __FILE__,__LINE__)
    endif

    if (associated(out%gboxx)) deallocate(out%gboxx)
    if (associated(out%gboxy)) deallocate(out%gboxy)
    if (associated(out%gboxn)) deallocate(out%gboxn)

    allocate(out%gboxx(size(in%gboxx,1),size(in%gboxx,2)))
    allocate(out%gboxy(size(in%gboxy,1),size(in%gboxy,2)))
    allocate(out%gboxn(size(in%gboxn,1),size(in%gboxn,2)))

    out%gboxx=in%gboxx
    out%gboxy=in%gboxy
    out%gboxn=in%gboxn

    out%set=.true.

  end subroutine copy_upscale

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  logical function lon_between(a,b,x)

    !*FD Checks to see whether a 
    !*FD longitudinal coordinate is between two bounds,
    !*FD taking into account the periodic boundary conditions.
    !*RV Returns \texttt{.true.} if $\mathtt{x}\geq \mathtt{a}$ and $\mathtt{x}<\mathtt{b}$.

    ! Arguments

    real(rk),intent(in) :: a  !*FD Lower bound on interval for checking
    real(rk),intent(in) :: b  !*FD Upper bound on interval for checking
    real(rk),intent(in) :: x  !*FD Test value (degrees)

    ! Internal variables

    real(rk) :: ta,tb

    ! Beginning of code

    if (a<b) then
       lon_between=((x>=a).and.(x<b))
    else
       if (x<a) then
          ta=a-360.0
          tb=b
       else 
          ta=a
          tb=b+360.0
       endif
       lon_between=((x>=ta).and.(x<tb))
    endif

  end function lon_between

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine calc_fractional(x,y,xp,yp,xa,ya,xb,yb,xc,yc,xd,yd)

    !*FD Performs a coordinate transformation to locate the point
    !*FD $(X',Y')$ fractionally within an arbitrary quadrilateral, 
    !*FD defined by the points $(x_A,y_A)$, $(x_B,y_B)$, 
    !*FD $(x_C,y_C)$ and $(x_D,y_D)$, which are ordered 
    !*FD anticlockwise.

    real(rk),intent(out) :: x !*FD The fractional $x$ location.
    real(rk),intent(out) :: y !*FD The fractional $y$ location.
    real(rk),intent(in)  :: xp,yp,xa,ya,xb,yb,xc,yc,xd,yd

    real(dp) :: a,b,c
    real(dp),parameter :: small=1d-8

    a=(yb-ya)*(xc-xd)-(yc-yd)*(xb-xa)

    b=xp*(yc-yd)-yp*(xc-xd) &
         +xd*(yb-ya)-yd*(xb-xa) &
         -xp*(yb-ya)+yp*(xb-xa) &
         -xa*(yc-yd)+ya*(xc-xd) 

    c=xp*(yd-ya)+yp*(xa-xd)+ya*xd-xa*yd

    if (abs(a)>small) then
       x=(-b-sqrt(b**2-4*a*c))/(2*a)
    else
       x=-c/b
    endif

    y=(yp-ya-x*(yb-ya))/(yd+x*(yc-yd-yb+ya)-ya)

#ifdef GLC_DEBUG
! Could use the following code if points are degenerate (a=b, c=d, etc.)
!    if (abs(a) > small) then
!       x=(-b-sqrt(b**2-4*a*c))/(2*a)
!    elseif (abs(b) > small) then
!       x=-c/b
!    else
!       x=0._rk
!    endif
!
!    if (abs(yd+x*(yc-yd-yb+ya)-ya) > small) then
!       y=(yp-ya-x*(yb-ya))/(yd+x*(yc-yd-yb+ya)-ya)
!    else
!       y=0._rk
!    endif
#endif

  end subroutine calc_fractional

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


end module glint_interp
