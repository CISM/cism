! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_nc_custom.f90 - part of the GLIMMER ice model      + 
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
    use paramets, only : len0
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
