! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_ncdf.f90 - part of the Glimmer-CISM ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2004, 2005, 2006, 2007, 2008, 2009, 2010
! Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
! This file is part of Glimmer-CISM.
!
! Glimmer-CISM is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 2 of the License, or (at
! your option) any later version.
!
! Glimmer-CISM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Glimmer-CISM.  If not, see <http://www.gnu.org/licenses/>.
!
! Glimmer-CISM is hosted on BerliOS.de:
! https://developer.berlios.de/projects/glimmer-cism/
!
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#define NCO outfile%nc
#define NCI infile%nc

!> netCDF type definitions and functions for managing linked lists
!!
!! \author Magnus Hagdorn
!! \date 2004
module glimmer_ncdf  

  use glimmer_global, only: fname_length
  use netcdf

  integer, parameter :: glimmer_nc_meta_len = 100 !< maximum length for meta data
  character(len=*), parameter :: glimmer_nc_mapvarname = 'mapping' !< name of the grid mapping variable
  real, parameter :: glimmer_nc_max_time=1.e10 !< maximum time that can be written

  !> Data structure holding netCDF file description
  type glimmer_nc_stat

     logical :: define_mode = .TRUE.               !< set to .TRUE. when we are in define mode
     logical :: just_processed = .FALSE.           !< set to .TRUE. if the file was used during the last time step
     real :: processsed_time = 0.0                 !< the time when the file was last processed
     character(len=fname_length) :: filename = " " !< name of netCDF file
     integer id                                    !< id of netCDF file

     integer :: nlevel = 0                         !< size of vertical coordinate

     integer timedim                               !< id of time dimension
     integer timevar                               !< id of time variable 
     character(len=310) vars                       !< string containing variables to be processed
     logical :: hotstart = .false.                 !< Set to true if we're writing a hotstart file
     character(len=310) vars_copy                  !< string containing variables to be processed (retained copy)
  end type glimmer_nc_stat     

  !> Data structure holding netCDF meta data, see CF user guide
  type glimmer_nc_meta
     
     character(len=glimmer_nc_meta_len) :: title = ''       !< title of netCDF file
     character(len=glimmer_nc_meta_len) :: institution = '' !< where the data was produced
     character(len=glimmer_nc_meta_len) :: references = ''  !< list of references
     character(len=glimmer_nc_meta_len) :: source = ''      !< this string will hold the GLIMMER version
     character(len=glimmer_nc_meta_len) :: history = ''     !< netCDF file history string
     character(len=glimmer_nc_meta_len) :: comment = ''     !< some comments
     character(len=10000) :: config = ''                    !< the contents of the glide config file
  end type glimmer_nc_meta

  !> element of linked list describing netCDF output file
  type glimmer_nc_output
     !NO_RESTART previous

     type(glimmer_nc_stat) :: nc                          !< structure containg file info
     real :: freq=1000                                    !< frequency at which data is written to file
     real :: next_write=0                                 !< next time step at which data is dumped
     real :: end_write=glimmer_nc_max_time                !< stop writing after this year
     integer :: timecounter=1                             !< time counter
     real :: total_time = 0.0                             !< accumulate time steps (used for taking time averages)

     integer :: default_xtype = NF90_REAL                 !< the default external type for storing flointing point values

     logical :: do_averages = .False.                     !< set to .True. if we need to handle averages
     
     type(glimmer_nc_meta) :: metadata                    !< structure holding metadata

     type(glimmer_nc_output), pointer :: next=>NULL()     !< next element in list
     type(glimmer_nc_output), pointer :: previous=>NULL() !< previous element in list
     logical :: append = .false.                          !< Set to true if we are appending onto an existing file.
  end type glimmer_nc_output

  !> element of linked list describing netCDF input file
  type glimmer_nc_input
     !NO_RESTART previous
     type(glimmer_nc_stat) :: nc                          !< structure containg file info
     integer, pointer, dimension(:) :: times => NULL()    !< pointer to array holding times
     integer                        :: nt, current_time=1 !< number of elements in times and current time index
     integer                        :: get_time_slice = 1 !< -1 if all times should be loaded, > 0 to load particular slice and then close file

     type(glimmer_nc_input), pointer :: next=>NULL()      !< next element in list
     type(glimmer_nc_input), pointer :: previous=>NULL()  !< previous element in list
  end type glimmer_nc_input

  interface delete
     module procedure delete_output, delete_input
  end interface

  interface add
     module procedure add_output, add_input
  end interface

contains

  !> remove output file from linked list
  !!
  !! \return pointer to the next element in the linked list or NULL()
  function delete_output(oc, cf)
    use glimmer_log
    implicit none
    type(glimmer_nc_output), pointer :: delete_output
    type(glimmer_nc_output), pointer :: oc !< the output file to be removed
    logical, intent(in), optional :: cf    !< set to .True. if file should be closed
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
  
  !> remove input file from linked list
  !!
  !! \return the next input file or NULL()
  function delete_input(ic,cf)
    !*FD remove element from linked list
    use glimmer_log
    implicit none
    type(glimmer_nc_input), pointer :: delete_input
    type(glimmer_nc_input), pointer :: ic !< the input file to be removed
    logical, intent(in), optional :: cf   !< set to .True. if file should be closed

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

  !> add a new output file
  !!
  !! \return pointer to added output file
  function add_output(oc)
    implicit none
    type(glimmer_nc_output), pointer :: add_output
    type(glimmer_nc_output), pointer :: oc !< the output file to be added

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

  !> add a new input file
  !!
  !! \return pointer to added input file
  function add_input(ic)
    implicit none
    type(glimmer_nc_input), pointer :: add_input
    type(glimmer_nc_input), pointer :: ic !< the input file to be added

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

  !> for debugging print all output files in linked list
  recursive subroutine nc_print_output(output)

    type(glimmer_nc_output),pointer :: output !< the first element in linled list

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

  !> print file stats
  subroutine nc_print_stat(stat)

    type(glimmer_nc_stat) :: stat !< the stats type

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

  !> Sets up previous points in the linked list correctly
  !!
  !! This is needed after a restart, as trying to save both
  !! next and previous pointers would cause problems
  !! Also resets some other internal components
  subroutine nc_repair_outpoint(output)

    implicit none

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

  !> Sets up previous points in the linked list correctly
  !!
  !! This is needed after a restart, as trying to save both
  !! next and previous pointers would cause problems
  subroutine nc_repair_inpoint(input)

    implicit none

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
  
  !> Adds a prefix to all the filenames stored in the linked list.
  !! Used for restarts.
  subroutine nc_prefix_outfiles(output,prefix)

    implicit none

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

!> handle netCDF error
subroutine nc_errorhandle(file,line,status)
  use netcdf
  use glimmer_log
  implicit none
  character(len=*), intent(in) :: file !< name of f90 file error occured in
  integer, intent(in) :: line          !< line number error occured at
  integer, intent(in) :: status        !< netCDF return value
  
  if (status.ne.NF90_NOERR) then
     call write_log(nf90_strerror(status),type=GM_FATAL,file=file,line=line)
  end if
end subroutine nc_errorhandle
