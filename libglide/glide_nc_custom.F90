! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_nc_custom.f90 - part of the Glimmer-CISM ice model + 
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

module glide_nc_custom
  !*FD module for filling in dimension variables
contains
  subroutine glide_nc_fillall(model,outfiles)
    !*FD fill dimension variables of all files
    use glide_types
    use glimmer_ncdf
    use glimmer_ncio
    implicit none
    type(glide_global_type) :: model
    type(glimmer_nc_output),pointer,optional :: outfiles

    ! local variables
    type(glimmer_nc_output), pointer :: oc

    if (present(outfiles)) then
       oc => outfiles
    else
       oc=>model%funits%out_first
    end if

    do while(associated(oc))
       if (.not.oc%append) call glide_nc_filldvars(oc,model)
       oc=>oc%next
    end do
  end subroutine glide_nc_fillall

  subroutine glide_nc_filldvars(outfile,model)
    use glide_types
    use glimmer_ncdf
    use glimmer_paramets, only : len0
    implicit none
    type(glimmer_nc_output), pointer :: outfile
    type(glide_global_type) :: model

    integer i,status,varid

    ! check if we are still in define mode and if so leave it
    if (NCO%define_mode) then
       status = nf90_enddef(NCO%id)
       call nc_errorhandle(__FILE__,__LINE__,status)
       NCO%define_mode = .FALSE.
    end if

    if (associated(model%funits%in_first)) then
    status = nf90_inq_varid(NCO%id,'x1',varid)
       status=nf90_put_var(NCO%id,varid,model%general%x1)
       call nc_errorhandle(__FILE__,__LINE__,status)
    status = nf90_inq_varid(NCO%id,'y1',varid)
       status=nf90_put_var(NCO%id,varid,model%general%y1)
       call nc_errorhandle(__FILE__,__LINE__,status)
    !create the x0 and y0 grids from x1 and y1
    status = nf90_inq_varid(NCO%id,'x0',varid)
    do i=1, model%general%ewn-1
      status=nf90_put_var(NCO%id,varid,(model%general%x1(i)+model%general%x1(i+1))/2.0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    status = nf90_inq_varid(NCO%id,'y0',varid)
    do i=1, model%general%nsn-1
       status=nf90_put_var(NCO%id,varid,(model%general%y1(i)+model%general%y1(i+1))/2.0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    
    else if(.not. associated(model%funits%in_first)) then
    ! filling coordinate variables
    status = nf90_inq_varid(NCO%id,'x0',varid)
    do i=1, model%general%ewn-1
       status=nf90_put_var(NCO%id,varid,((i-0.5)*model%numerics%dew*len0),(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    status = nf90_inq_varid(NCO%id,'y0',varid)
    do i=1, model%general%nsn-1
       status=nf90_put_var(NCO%id,varid,(i-0.5)*model%numerics%dns*len0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    status = nf90_inq_varid(NCO%id,'x1',varid)
    do i=1, model%general%ewn
       status=nf90_put_var(NCO%id,varid,(i-1.)*model%numerics%dew*len0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    status = nf90_inq_varid(NCO%id,'y1',varid)
    do i=1, model%general%nsn
       status=nf90_put_var(NCO%id,varid,(i-1.)*model%numerics%dns*len0,(/i/))
       call nc_errorhandle(__FILE__,__LINE__,status)
    end do
    end if

    status = nf90_inq_varid(NCO%id,'level',varid)
    status=nf90_put_var(NCO%id,varid,model%numerics%sigma)
    call nc_errorhandle(__FILE__,__LINE__,status)
    if (model%options%gthf.gt.0) then
       status = nf90_inq_varid(NCO%id,'lithoz',varid)
       status=nf90_put_var(NCO%id,varid,model%lithot%deltaz)
       call nc_errorhandle(__FILE__,__LINE__,status)
    end if
  end subroutine glide_nc_filldvars
end module glide_nc_custom
