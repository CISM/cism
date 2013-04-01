! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glimmer_ncdf.f90 - part of the Glimmer-CISM ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2006, 2007, 2008, 2009, 2010
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

  use glimmer_global, only: fname_length, dp
  use netcdf

  implicit none

  integer, parameter :: glimmer_nc_meta_len = 100
  !*FD maximum length for meta data

  character(len=*), parameter :: glimmer_nc_mapvarname = 'mapping'
  !*FD name of the grid mapping variable

  real(dp), parameter :: glimmer_nc_max_time=1.d10
  !*FD maximum time that can be written

  !> Data structure holding netCDF file description
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
     integer :: nstaglevel = 0
     integer :: nstagwbndlevel = 0
     !*FD size of vertical and stag vertical coordinate

     integer timedim
     !*FD id of time dimension
     integer timevar
     !*FD id of time variable 
     ! TODO Create a variable for vars length so it can be made longer - apply it here for vars, vars_copy and to restart_variable_list in glimmer_ncparams.F90
     character(len=310) vars
     !*FD string containing variables to be processed
     logical :: restartfile = .false.
     !*FD Set to true if we're writing a restart file
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

     type(glimmer_nc_stat) :: nc                          !< structure containg file info
!WHLTSTEP - Changed freq, next_write, end_write and total_time to dp
!     real :: freq = 1000                                  !< frequency at which data is written to file
!     real :: next_write=0                                 !< next time step at which data is dumped
!     real :: end_write=glimmer_nc_max_time                !< stop writing after this year
     real(dp) :: freq = 1000.d0                           !< frequency at which data is written to file
     real(dp) :: next_write = 0.d0                        !< next time step at which data is dumped
     real(dp) :: end_write = glimmer_nc_max_time          !< stop writing after this year
     integer :: timecounter = 1                           !< time counter
!     real :: total_time = 0.0                             !< accumulate time steps (used for taking time averages)
     real(dp) :: total_time = 0.d0                        !< accumulate time steps (used for taking time averages)

!TODO - I'd like to understand how the default xtype is overridden
     integer :: default_xtype = NF90_REAL                 !< the default external type for storing flointing point values
     logical :: do_averages = .false.                     !< set to .true. if we need to handle averages

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
!WHLTSTEP - Changed times pointer to real(dp)
!     integer, pointer, dimension(:) :: times => NULL()     
     real(dp), pointer, dimension(:) :: times => NULL()     
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

contains

  function delete_output(oc, cf)
    !*FD remove element from linked list
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
    print*,'nstaglevel:     ',stat%nstaglevel
    print*,'nstagwbndlevel: ',stat%nstagwbndlevel
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
  
        if (status /= NF90_NOERR) then
            call write_log(nf90_strerror(status),type=GM_FATAL,file=file,line=line)
        end if
    end subroutine nc_errorhandle

!TODO - The convention now is that we have unstaggered temps for dycore = 0 (glide)
!       and staggered temps for dycore = 1 (glissade).

    subroutine check_for_tempstag(whichtemp,nc)
      !*FD check for the need to output tempstag and update the output variables if needed.
      ! When IR is used to evolve temperature, the temperature grid has an extra
      ! layer.  In that case, the netCDF output file should include a variable
      ! called tempstag(0:nz) instead of temp(1:nz). This subroutine is added for
      ! convenience to allow the variable "temp" to be specified in the config
      ! file in all cases and have it converted to "tempstag" when appropriate.
      ! MJH Nov. 2010

      use glimmer_log

      implicit none
      integer, intent(in) :: whichtemp
      type(glimmer_nc_stat) :: nc 

      ! Locals
      integer :: i

      ! Check if tempstag should be output
      ! TODO If both temp and tempstag are specfied, should one be removed?
      ! TODO Modify this to work if multiple output files are specified?

!TODO - Change from 'whichtemp' to 'whichdycore'?  At any rate, do not hardwire the whichtemp index
      !print *, "Original varstring:", varstring
      if (whichtemp == 3) then !3) then !3) then !TEMP_REMAP_ADV) then
          ! We want temp to become tempstag 
          i = index(nc%vars, " temp ")
          if (i > 0) then
            ! temp was specified - change it to tempstag
            ! If temp is listed more than once, this just changes the first instance
            nc%vars = nc%vars(1:i-1) // " tempstag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses temperature on a staggered grid.' // &
              '  The netCDF output variable "temp" has been changed to "tempstag".' )
          endif
          ! Now check if flwa needs to be changed to flwastag
          i = index(nc%vars, " flwa ") ! Look for flwa
          if (i > 0) then
            ! flwa was specified - change to flwastag
            nc%vars = nc%vars(1:i-1) // " flwastag " // nc%vars(i+6:len(nc%vars))
            call write_log('Temperature remapping option uses flwa on a staggered grid.' // &
            '  The netCDF output variable "flwa" has been changed to "flwastag".' )
          endif
      else  ! whichtemp is not IR
          ! We want tempstag to become temp
          i = index(nc%vars, " tempstag ")
          if (i > 0) then
            !Change tempstag to temp
            nc%vars = nc%vars(1:i-1) // " temp " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "tempstag" should only be used when remapping temperature.' // &
              '  The netCDF output variable "tempstag" has been changed to "temp".' )
          endif
          ! We want flwastag to become flwa
          i = index(nc%vars, " flwastag ")
          if (i > 0) then
            !Change flwastag to flwa
            nc%vars = nc%vars(1:i-1) // " flwa " // nc%vars(i+10:len(nc%vars))
            call write_log('The netCDF output variable "flwastag" should only be used when remapping temperature.' // &
              '  The netCDF output variable "flwastag" has been changed to "flwa".' )
          endif
      endif  !whichtemp == 3
      ! Copy any changes to vars_copy
      nc%vars_copy = nc%vars
    end subroutine check_for_tempstag


end module glimmer_ncdf


