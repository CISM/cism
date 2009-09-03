! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_ncdf.f90 - part of the GLIMMER ice model         + 
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

#define NCO outfile%nc
#define NCI infile%nc

module glimmer_ncdf  
  !*FD netCDF type definitions and functions for managing linked lists
  !*FD written by Magnus Hagdorn, 2004

  use glimmer_global, only: fname_length
  use netcdf

  integer, parameter :: glimmer_nc_meta_len = 100
  !*FD maximum length for meta data
  character(len=*), parameter :: glimmer_nc_mapvarname = 'mapping'
  !*FD name of the grid mapping variable
  real, parameter :: glimmer_nc_max_time=1.e10
  !*FD maximum time that can be written

  type glimmer_nc_stat
     !*FD Data structure holding netCDF file description

     logical :: define_mode = .TRUE.
     !*FD set to .TRUE. when we are in define mode
     logical :: just_processed = .FALSE.
     !*FD set to .TRUE. if the file was used during the last time step
     real :: processsed_time = 0.0
     !*FD the time when the file was last processed
     character(len=fname_length) :: filename = " "
     !*FD name of netCDF file
     integer id
     !*FD id of netCDF file

     integer :: nlevel = 0
     !*FD size of vertical coordinate

     integer timedim
     !*FD id of time dimension
     integer timevar
     !*FD id of time variable 
     character(len=310) vars
     !*FD string containing variables to be processed
     logical :: hotstart = .false.
     !*FD Set to true if we're writing a hotstart file
     character(len=310) vars_copy
     !*FD string containing variables to be processed (retained copy)
  end type glimmer_nc_stat     

  type glimmer_nc_meta
     !*FD Data structure holding netCDF meta data, see CF user guide
     
     character(len=glimmer_nc_meta_len) :: title = ''
     !*FD title of netCDF file
     character(len=glimmer_nc_meta_len) :: institution = ''
     !*FD where the data was produced
     character(len=glimmer_nc_meta_len) :: references = ''
     !*FD list of references
     character(len=glimmer_nc_meta_len) :: source = ''
     !*FD this string will hold the GLIMMER version
     character(len=glimmer_nc_meta_len) :: history = ''
     !*FD netCDF file history string
     character(len=glimmer_nc_meta_len) :: comment = ''
     !*FD some comments
     character(len=10000) :: config = ''
     !*FD the contents of the glide config file
  end type glimmer_nc_meta

  type glimmer_nc_output
     !*FD element of linked list describing netCDF output file
     !NO_RESTART previous

     type(glimmer_nc_stat) :: nc
     !*FD structure containg file info
     real :: freq=1000
     !*FD frequency at which data is written to file
     real :: next_write=0
     !*FD next time step at which data is dumped
     real :: end_write=glimmer_nc_max_time
     !*FD stop writing after this year
     integer :: timecounter=1
     !*FD time counter
     real :: total_time
     !*FD accumulate time steps (used for taking time averages)

     logical :: do_averages = .False. !*FD set to .True. if we need to handle averages
     
     type(glimmer_nc_meta) :: metadata
     !*FD structure holding metadata

     type(glimmer_nc_output), pointer :: next=>NULL()
     !*FD next element in list
     type(glimmer_nc_output), pointer :: previous=>NULL()
     !*FD previous element in list
     logical :: append = .false.
     !*FD Set to true if we are appending onto an existing file.
  end type glimmer_nc_output

  type glimmer_nc_input
     !*FD element of linked list describing netCDF input file
     !NO_RESTART previous
     type(glimmer_nc_stat) :: nc
     !*FD structure containg file info
     integer, pointer, dimension(:) :: times => NULL()     
     !*FD pointer to array holding times
     integer                        :: nt, current_time=1
     !*FDnumber of elements in times and current time index
     integer                        :: get_time_slice = 1     
     !*FD -1 if all times should be loaded, > 0 to load particular slice and then close file

     type(glimmer_nc_input), pointer :: next=>NULL()
     !*FD next element in list
     type(glimmer_nc_input), pointer :: previous=>NULL()
     !*FD previous element in list
  end type glimmer_nc_input

  interface delete
     module procedure delete_output, delete_input
  end interface

  interface add
     module procedure add_output, add_input
  end interface

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_NCDF
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_NCDF
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_NCDF
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_NCDF
!MH!#endif

  function delete_output(oc, cf)
    !*FD remove element from linked list
    use glimmer_log
    implicit none
    type(glimmer_nc_output), pointer :: delete_output
    type(glimmer_nc_output), pointer :: oc
    logical, intent(in), optional :: cf
    ! local variables
    logical closefile
    integer status

    if (present(cf)) then
       closefile = cf
    else
       closefile = .true.
    end if

    if (associated(oc)) then
       if (associated(oc%previous)) then
          oc%previous%next => oc%next
       end if
       if (associated(oc%next)) then
          oc%next%previous => oc%previous
          delete_output => oc%next
       else
          delete_output => NULL()
       end if
       if (closefile) then
          status = nf90_close(oc%nc%id)
          call write_log_div
          call write_log('Closing output file '//trim(oc%nc%filename))
       end if
       deallocate(oc)
    end if
  end function delete_output
  
  function delete_input(ic,cf)
    !*FD remove element from linked list
    use glimmer_log
    implicit none
    type(glimmer_nc_input), pointer :: delete_input
    type(glimmer_nc_input), pointer :: ic
    logical, intent(in), optional :: cf

    ! local variables
    logical closefile
    integer status

    if (present(cf)) then
       closefile = cf
    else
       closefile = .true.
    end if

    if (associated(ic)) then
       if (associated(ic%previous)) then
          ic%previous%next => ic%next
       end if
       if (associated(ic%next)) then
          ic%next%previous => ic%previous
          delete_input => ic%next
       else
          delete_input => NULL()
       end if
       if (closefile) then
          status = nf90_close(ic%nc%id)
          call write_log_div
          call write_log('Closing input file '//trim(ic%nc%filename))
       end if
       deallocate(ic%times)
       deallocate(ic)
    end if
  end function delete_input

  function add_output(oc)
    !*FD add new element to linked list
    implicit none
    type(glimmer_nc_output), pointer :: add_output
    type(glimmer_nc_output), pointer :: oc

    allocate(add_output)

    if (associated(oc)) then
       add_output%previous => oc
       if (associated(oc%next)) then
          add_output%next => oc%next
          oc%next%previous => add_output
       end if
       oc%next => add_output
    end if
  end function add_output

  function add_input(ic)
    !*FD add new element to linked list
    implicit none
    type(glimmer_nc_input), pointer :: add_input
    type(glimmer_nc_input), pointer :: ic

    allocate(add_input)

    if (associated(ic)) then
       add_input%previous => ic
       if (associated(ic%next)) then
          add_input%next => ic%next
          ic%next%previous => add_input
       end if
       ic%next => add_input
    end if
  end function add_input

  recursive subroutine nc_print_output(output)

    !*FD For debugging

    type(glimmer_nc_output),pointer :: output

    if (.not.associated(output)) then
       Print*,'*** Output section not associated'
       return
    end if

    call nc_print_stat(output%nc)
    print*,'freq:       ',output%freq
    print*,'next_write: ',output%next_write
    print*,'timecounter:',output%timecounter
    ! call nc_print_meta(output%metadata)
    if (associated(output%next)) call nc_print_output(output%next)
    
  end subroutine nc_print_output

  subroutine nc_print_stat(stat)

    type(glimmer_nc_stat) :: stat

    print*,'define_mode:    ',stat%define_mode
    print*,'just_processed: ',stat%just_processed
    print*,'processsed_time:',stat%processsed_time
    print*,'filename:       ',stat%filename
    print*,'id:             ',stat%id
    print*,'nlevel:         ',stat%nlevel
    print*,'timedim:        ',stat%timedim
    print*,'timevar:        ',stat%timevar
    print*,'vars:           ',trim(stat%vars)

  end subroutine nc_print_stat

  subroutine nc_repair_outpoint(output)

    implicit none

    !*FD Sets up previous points in the linked list correctly
    !*FD This is needed after a restart, as trying to save both
    !*FD next and previous pointers would cause problems
    !*FD Also resets some other internal components

    type(glimmer_nc_output),pointer :: output
    type(glimmer_nc_output),pointer :: most_recent
    type(glimmer_nc_output),pointer :: tmp

    most_recent => null()
    if (.not.associated(output)) return
    tmp => output

    do
       if (associated(most_recent)) tmp%previous => most_recent
       tmp%nc%vars=tmp%nc%vars_copy
       if (.not.associated(tmp%next)) exit
       most_recent => tmp
       tmp => tmp%next
    end do

  end subroutine nc_repair_outpoint

  subroutine nc_repair_inpoint(input)

    implicit none

    !*FD Sets up previous points in the linked list correctly
    !*FD This is needed after a restart, as trying to save both
    !*FD next and previous pointers would cause problems

    type(glimmer_nc_input),pointer :: input
    type(glimmer_nc_input),pointer :: most_recent
    type(glimmer_nc_input),pointer :: tmp

    most_recent => null()
    if (.not.associated(input)) return
    tmp => input

    do
       if (associated(most_recent)) tmp%previous => most_recent
       if (.not.associated(tmp%next)) exit
       most_recent => tmp
       tmp => tmp%next
    end do

  end subroutine nc_repair_inpoint

  subroutine nc_prefix_outfiles(output,prefix)

    !*FD Adds a prefix to all the filenames stored in the linked list.
    !*FD Used for restarts.

    type(glimmer_nc_output),pointer :: output
    character(*) :: prefix

    type(glimmer_nc_output),pointer :: tmp

    tmp => output
    do
       tmp%nc%filename=trim(prefix)//trim(tmp%nc%filename)
       if (.not.associated(tmp%next)) exit
       tmp => tmp%next
    end do

  end subroutine nc_prefix_outfiles

end module glimmer_ncdf

module glimmer_scales
  !*FD this module holds scales for various fields

  use glimmer_global, only : dp

  real(dp) :: scale2d_f1, scale2d_f2, scale2d_f3, scale2d_f4, scale2d_f5, scale2d_f6, scale2d_f7, scale2d_f8, scale2d_f9
  real(dp) :: scale3d_f1, scale3d_f2, scale3d_f3, scale3d_f4, scale3d_f5, scale3d_f6, scale3d_f7, scale3d_f8

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_SCALES
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_SCALES
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_SCALES
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_SCALES
!MH!#endif

  subroutine glimmer_init_scales
    !*FD calculate scale factors (can't have non-integer powers)
    use physcon, only : scyr, gn
    use paramets, only : thk0, tim0, vel0, vis0, len0, tau0, acc0
    implicit none

    scale2d_f1 = scyr * thk0 / tim0
    scale2d_f2 = scyr * vel0 * thk0
    scale2d_f3 = vel0 / (vis0 * len0)
    scale2d_f4 = vel0 * scyr * len0
    scale2d_f5 = scyr * vel0
    scale2d_f6 = scyr * vel0 * len0 / (thk0**2)
    scale2d_f7 = tau0
    scale2d_f8 = tau0 * len0 / (scyr * vel0)
    scale2d_f9 = scyr * acc0

    scale3d_f1 = scyr * vel0
    scale3d_f2 = vis0 * (vel0/len0)**(gn - 1)
    scale3d_f3 = scyr * thk0
    scale3d_f4 = vel0/(vis0*len0)
    scale3d_f5 = 1.0d0/scale3d_f2**(1.0/gn)
    scale3d_f6 = scale3d_f4**(1.0/gn)
    scale3d_f7 = scyr * thk0/tim0
    scale3d_f8 = vis0*scyr
  end subroutine glimmer_init_scales
end module glimmer_scales

subroutine nc_errorhandle(file,line,status)
  !*FD handle netCDF error
  use netcdf
  use glimmer_log
  implicit none
  character(len=*), intent(in) :: file
  !*FD name of f90 file error occured in
  integer, intent(in) :: line
  !*FD line number error occured at
  integer, intent(in) :: status
  !*FD netCDF return value
  
  if (status.ne.NF90_NOERR) then
     call write_log(nf90_strerror(status),type=GM_FATAL,file=file,line=line)
  end if
end subroutine nc_errorhandle
