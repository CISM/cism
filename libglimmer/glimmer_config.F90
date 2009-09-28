! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
! This code is taken from the Generic Mapping Tools and was 
! converted into Fortran 90 by Ian Rutt.
!
! Original code (in C) Copyright (c) 1991-2003 by P. Wessel and W. H. F. Smith
!
! Partial translation into Fortran 90 (c) 2004 Ian C. Rutt
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

module glimmer_config
  !*FD configuration file parser
  !*FD written by Magnus Hagdorn, May 2004
  !*FD everything is a singly linked list

  use glimmer_global, only : dp, msg_length
  private :: handle_section, handle_value, InsertSection, InsertValue, dp

  integer, parameter :: namelen=50
  integer, parameter :: valuelen=200
  integer, parameter :: linelen=250

  type ConfigValue
     character(len=namelen) :: name = ''
     character(len=valuelen) :: value
     type(ConfigValue), pointer :: next=>NULL()
  end type ConfigValue

  type ConfigSection
     character(len=namelen) :: name = ''
     logical :: used = .false.
     type(ConfigValue), pointer :: values=>NULL()
     type(ConfigSection), pointer :: next=>NULL()
  end type ConfigSection

  type ConfigData
     !*FD This type exists so that we can have
     !*FD arrays of config data, since f90 doesn't
     !*FD allow arrays of pointers
     type(ConfigSection), pointer :: config=>null()
  end type ConfigData

  interface GetValue
     module procedure GetValueDouble, GetValueReal, GetValueInt, GetValueChar, GetValueLogical, &
          GetValueDoubleArray, GetValueRealArray, GetValueIntArray, GetValueCharArray
  end interface

  interface ConfigSetValue
     module procedure ConfigSetValueData, ConfigSetValueSec
  end interface

  interface ConfigCombine
     module procedure ConfigCombineData, ConfigCombineSec, ConfigCombineDataSec, ConfigCombineSecData
  end interface

!MH!  !MAKE_RESTART
!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_CONFIG
!MH!#include "glimmer_rst_head.inc"
!MH!#undef RST_GLIMMER_CONFIG
!MH!#endif

contains

!MH!#ifdef RESTARTS
!MH!#define RST_GLIMMER_CONFIG
!MH!#include "glimmer_rst_body.inc"
!MH!#undef RST_GLIMMER_CONFIG
!MH!#endif

  subroutine ConfigRead(fname,config)
    !*FD read configuration file
    use glimmer_log
    implicit none
    character(len=*), intent(in) :: fname
    !*FD name of configuration file
    type(ConfigSection), pointer :: config

    ! local variables
    type(ConfigSection), pointer :: this_section
    type(ConfigValue), pointer ::  this_value
    logical there
    integer unit,ios,linenr
    character(len=linelen) :: line
    character(len=msg_length) :: message

    inquire (exist=there,file=fname)
    if (.not.there) then
       call write_log('Cannot open configuration file '//trim(fname),GM_FATAL)
    end if
    
    unit=99
    open(unit,file=trim(fname),status='old')
    ios=0
    linenr=0
    config=>NULL()
    this_section=>NULL()
    do while(ios.eq.0)
       read(unit,fmt='(a250)',iostat=ios) line
       line = adjustl(line)
       if (ios.ne.0) then
          exit
       end if
       if (.not.(line(1:1).eq.'!' .or. line(1:1).eq.'#' .or. line(1:1).eq.';' .or. line(1:1).eq.' ')) then
          ! handle comments
          if (line(1:1).eq.'[') then
             ! new section
             call handle_section(linenr,line,this_section)
             this_value=>NULL()
             if (.not.associated(config)) then
                ! this is the first section in config file
                config=>this_section
             end if
          else
             ! handle value
             if (.not.associated(this_section)) then
                call write_log('No section defined yet',GM_ERROR)
                write(message,*) trim(adjustl(fname)), linenr
                call write_log(message,GM_FATAL)
             end if
             call handle_value(linenr,line,this_value)
             if (.not.associated(this_section%values)) then
                this_section%values => this_value
             end if
          end if
       end if
       linenr = linenr + 1
    end do
    close(unit)
    return
  end subroutine ConfigRead

  subroutine PrintConfig(config)
    implicit none
    type(ConfigSection), pointer :: config

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val
    
    sec=>config
    do while(associated(sec))
       write(*,*) sec%name
       val=>sec%values
       do while(associated(val))
          write(*,*) '  ',trim(val%name),' == ', trim(val%value)
          val=>val%next
       end do
       write(*,*)
       sec=>sec%next
    end do
  end subroutine PrintConfig

  subroutine ConfigAsString(config,string)
    use glimmer_global, only: endline
    implicit none
    type(ConfigSection), pointer :: config
    character(*),intent(out) :: string

    type(ConfigSection), pointer :: sec
    type(ConfigValue), pointer ::  val
    
    string=''

    sec=>config
    do while(associated(sec))
       string=trim(string)//'['//trim(sec%name)//']'//trim(endline)
       val=>sec%values
       do while(associated(val))
          string=trim(string)//trim(val%name)//': '//trim(val%value)//trim(endline)
          val=>val%next
       end do
       sec=>sec%next
    end do
  end subroutine ConfigAsString

  subroutine ConfigSetValueData(config,secname,valname,value,tag)
    !*FD Either overwrite a given key-value pair,
    !*FD or create a new one

    type(ConfigData) :: config
    character(*) :: secname,valname,value
    character(*),optional :: tag

    call ConfigSetValueSec(config%config,secname,valname,value,tag)

  end subroutine ConfigSetValueData

  subroutine ConfigSetValueSec(config,secname,valname,value,tag)
    !*FD Either overwrite a given key-value pair,
    !*FD or create a new one
    !*FD tag is a label that you can give to a particular config section
    !*FD allowing the identification of the right section when
    !*FD several with the same name are present (e.g. [CF output])

    type(ConfigSection), pointer :: config
    character(*) :: secname,valname,value
    character(*),optional :: tag
    type(ConfigSection), pointer :: found
    type(ConfigSection), pointer :: newsec
    type(ConfigValue), pointer :: val
    type(ConfigValue), pointer :: newval
    type(ConfigValue), pointer :: newtag
    logical :: tagflag

    ! Find or create correct section

    if (.not.associated(config)) allocate(config)

    found=>config
    do
       if (associated(found)) then
          if (present(tag)) then
             tagflag=ConfigSectionHasTag(found,tag)
          else
             tagflag=.true.
          end if
          if ((trim(secname)==trim(found%name)).and.tagflag) then
             exit
          else
             if (associated(found%next)) then
                found=>found%next
             else
                allocate(newsec)
                found%next=>newsec
                found=>found%next
                found%name=trim(secname)
                if (present(tag)) then
                   allocate(newtag)
                   newtag%name='tag'
                   newtag%value=trim(tag)
                   found%values=>newtag
                end if
                exit
             end if
          end if
       else
          exit
       end if
    end do
 
    ! Add or create key-value pair

    if (.not.associated(found%values)) then
       allocate(newval)
       found%values=>newval
       found%values%name=valname
       found%values%value=value
    else
       val=>found%values
       do
          if (trim(valname)==trim(val%name)) then
             val%value=value
             exit
          else
             if (associated(val%next)) then
                val=>val%next
             else
                allocate(newval)
                val%next=>newval
                val%next%name=valname
                val%next%value=value
                exit
             end if
          end if
       end do
    end if

  end subroutine ConfigSetValueSec

  subroutine ConfigCombineDataSec(config1,config2)
    !*FD Add the contents of config2 to config1,
    !*FD overwriting if necessary

    type(ConfigData) :: config1
    type(ConfigSection),pointer :: config2

    call ConfigCombineSec(config1%config,config2)

  end subroutine ConfigCombineDataSec

  subroutine ConfigCombineSecData(config1,config2)
    !*FD Add the contents of config2 to config1,
    !*FD overwriting if necessary

    type(ConfigSection),pointer :: config1
    type(ConfigData) :: config2

    call ConfigCombineSec(config1,config2%config)

  end subroutine ConfigCombineSecData


  subroutine ConfigCombineData(config1,config2)
    !*FD Add the contents of config2 to config1,
    !*FD overwriting if necessary

    type(ConfigData) :: config1
    type(ConfigData) :: config2

    call ConfigCombineSec(config1%config,config2%config)

  end subroutine ConfigCombineData

  subroutine ConfigCombineSec(config1,config2)
    !*FD Add the contents of config2 to config1,
    !*FD overwriting if necessary

    type(ConfigSection), pointer :: config1
    type(ConfigSection), pointer :: config2

    type(ConfigSection), pointer :: thissec
    type(ConfigValue),   pointer :: thisval
    character(namelen) :: thisname

    character(150) :: tag

    thissec=>config2
    do
       if (associated(thissec)) then
          thisval=>thissec%values
          thisname=trim(thissec%name)
          do
             if (associated(thisval)) then
                if (ConfigSectionHasValue(thissec,'tag',tag)) then
                   call ConfigSetValue(config1,thisname,trim(thisval%name),trim(thisval%value),tag=tag)
                else
                   call ConfigSetValue(config1,thisname,trim(thisval%name),trim(thisval%value))
                end if
                thisval=>thisval%next
             else
                exit
             end if
          end do
          thissec=>thissec%next
       else
          exit
       end if
    end do

  end subroutine ConfigCombineSec

  logical function ConfigSectionHasTag(section,tag)
    
    type(ConfigSection), pointer :: section
    character(*) :: tag
    character(200) :: testtag

    ConfigSectionHasTag=.false.
    if (ConfigSectionHasValue(section,'tag',testtag)) then
       if (trim(tag)==trim(testtag)) then
          ConfigSectionHasTag=.true.
       end if
    end if

  end function ConfigSectionhasTag

  logical function ConfigSectionHasValue(section,valname,val)

    type(ConfigSection), pointer :: section
    type(ConfigValue), pointer :: thisval
    character(*) :: valname,val

    ConfigSectionHasValue=.false.
    val=''

    if (.not.associated(section)) return

    thisval=>section%values
    do
       if (.not.associated(thisval)) exit
       if (trim(valname)==trim(thisval%name)) then
          val=trim(thisval%value)
          ConfigSectionHasValue=.true.
          exit
       else
          thisval=>thisval%next
       end if
    end do
 
  end function ConfigSectionHasValue

  subroutine GetSection(config,found,name)
    !*FD Find and return section with name
    implicit none
    type(ConfigSection), pointer :: config
    type(ConfigSection), pointer :: found
    character(len=*),intent(in) :: name

    found=>config
    do while(associated(found))
       if (name.eq.trim(found%name)) then
          found%used = .true.
          return
       end if
       found=>found%next
    end do
  end subroutine GetSection

  subroutine CheckSections(config)
    !*FD traverse linked list and check that all sections have been used
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: config
    
    ! local variables
    type(ConfigSection), pointer :: cf

    cf=>config
    do while(associated(cf))
       if (.not.cf%used) then
          call write_log('Unused section: '//trim(cf%name),GM_WARNING)
       end if
       cf=>cf%next
    end do
  end subroutine CheckSections

  subroutine GetValueDoubleArray(section,name,val,numval)
    !*FD get real array value
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(kind=dp), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value,tmp
    real(kind=dp), dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do
    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueDoubleArray

  subroutine GetValueRealArray(section,name,val,numval)
    !*FD get real array value
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real, pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value,tmp
    real, dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueRealArray

  subroutine GetValueIntArray(section,name,val,numval)
    !*FD get integer array value
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    integer, pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value,tmp
    integer, dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tmp=value(1:ind-1)
       read(tmp,*,err=10)tempval(i)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueIntArray

  subroutine GetValueCharArray(section,name,val,numval)
    !*FD get character array value
    use glimmer_log
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    character(len=80), pointer, dimension(:) :: val
    integer,intent(in), optional :: numval

    ! local variables
    character(len=valuelen) :: value
    character(80), dimension(:),allocatable :: tempval
    integer i,numv,inds,indc,ind

    if (present(numval)) then
       numv=numval
    else
       numv=100
    end if
    allocate(tempval(numv))
    value=''
    call GetValueChar(section,name,value)
    if (value.eq.'') return

    i=1
    do
       inds=index(value,' ') ; indc=index(value,',')
       if (inds==0.and.indc==0) then
          exit
       else if (inds==1.or.indc==1) then
          value=value(2:)
          cycle
       else if (inds==0) then
          ind=indc
       else if (indc==0) then
          ind=inds
       else
          ind=min(inds,indc)
       end if
       tempval(i)=value(1:ind-1)
       value=value(ind+1:)
       if (trim(value)=='') exit
       i=i+1
    end do

    if (i.ge.1) then
       if (associated(val)) then
          deallocate(val)
       end if
       allocate(val(i))
       val = tempval(1:i)
    end if
    return
    
10  call write_log('Array error in config file - check syntax',GM_FATAL)

  end subroutine GetValueCharArray

  subroutine GetValueReal(section,name,val)
    !*FD get real value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real :: val

    ! local variables
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueReal

  subroutine GetValueDouble(section,name,val)
    !*FD get double value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    real(kind=dp) :: val

    ! local variables
    character(len=valuelen) :: value
    real temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueDouble

  subroutine GetValueInt(section,name,val)
    !*FD get integer value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    integer :: val

    ! local variables
    character(len=valuelen) :: value
    integer temp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) temp
    if (ios==0) then
       val = temp
    end if
  end subroutine GetValueInt

  subroutine GetValueChar(section,name,val)
    !*FD get character value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    character(len=*) :: val
    
    type(ConfigValue), pointer :: value

    value=>section%values
    do while(associated(value))
       if (name.eq.trim(value%name)) then
          val = value%value
          return
       end if
       value=>value%next
    end do
  end subroutine GetValueChar

  subroutine GetValueLogical(section,name,val)
    !*FD get logical value
    implicit none
    type(ConfigSection), pointer :: section
    character(len=*),intent(in) :: name
    logical :: val

    ! local variables
    character(len=valuelen) :: value
    integer itemp
    logical ltemp
    integer ios

    value=''
    call GetValueChar(section,name,value)

    read(value,*,iostat=ios) itemp
    if (ios==0) then
       val = itemp.eq.1
    end if
    read(value,*,iostat=ios) ltemp
    if (ios==0) then
       val = ltemp
    end if
  end subroutine GetValueLogical

  !==================================================================================
  ! private procedures
  !==================================================================================

  subroutine handle_section(linenr,line,section)
    use glimmer_log
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigSection), pointer :: section

    ! local variables
    integer i
    character(len=msg_length) :: message

    do i=1,linelen
       if (line(i:i).eq.']') then
          exit
       end if
    end do
    if (line(i:i).ne.']') then
       write(message,*) 'Cannot find end of section ',linenr
       call write_log(message,GM_FATAL)
    end if

    call InsertSection(trim(adjustl(line(2:i-1))),section)
  end subroutine handle_section
  
  subroutine handle_value(linenr,line,value)
    use glimmer_log
    implicit none
    integer, intent(in) :: linenr
    character(len=*) :: line
    type(ConfigValue), pointer :: value

    ! local variables
    integer i
    character(len=msg_length) :: message
    do i=1,linelen
       if (line(i:i).eq.'=' .or. line(i:i).eq.':') then
          exit
       end if
    end do
    if (.not.(line(i:i).eq.'=' .or. line(i:i).eq.':')) then
       write(message,*) 'Cannot find = or : ',linenr
       call write_log(message,GM_FATAL)
    end if

    call InsertValue(trim(adjustl(line(:i-1))), trim(adjustl(line(i+1:))),value)
  end subroutine handle_value

  subroutine InsertSection(name,section)
    !*FD add a new section
    implicit none
    character(len=*), intent(in) :: name
    type(ConfigSection), pointer :: section
    type(ConfigSection), pointer :: new_sec

    allocate(new_sec)
    new_sec%name = name

    if (associated(section)) then
       if (associated(section%next)) then
          new_sec%next => section%next
       end if
       section%next=>new_sec
    end if
    section=>new_sec
  end subroutine InsertSection

  subroutine InsertValue(name,val,value)
    !*FD add a new value
    implicit none
    character(len=*), intent(in) :: name, val
    type(ConfigValue), pointer :: value
    type(ConfigValue), pointer :: new_value

    allocate(new_value)
    
    new_value%name = name
    new_value%value = val

    if(associated(value)) then
       if (associated(value%next)) then
          new_value%next => value%next
       end if
       value%next => new_value
    end if
    value=>new_value
  end subroutine InsertValue
end module glimmer_config
