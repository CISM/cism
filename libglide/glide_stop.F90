! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_stop.f90 - part of the Glimmer-CISM ice model      + 
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

module glide_stop
  !*FD module containing finalisation of glide
  !*FD this subroutine had to be split out from glide.f90 to avoid
  !*FD circular dependencies

contains
  
  subroutine glide_finalise(model,crash)
    !*FD finalise GLIDE model instance
    use glimmer_ncio
    use glimmer_log
    use glide_types
    use glide_io
    use glide_lithot_io
    use profile
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    logical, optional :: crash              !*FD set to true if the model died unexpectedly
    character(len=100) :: message

    ! force last write if crashed
    if (present(crash)) then
       if (crash) then
          call glide_io_writeall(model,model,.true.)
          if (model%options%gthf.gt.0) then
             call glide_lithot_io_writeall(model,model,.true.)
          end if
       end if
    end if

    call closeall_in(model)
    call closeall_out(model)
    
    call glide_deallocarr(model)

    ! write some statistics
    call write_log('Some Stats')
    write(message,*) 'Maximum temperature iterations: ',model%temper%niter
    call write_log(message)

    ! close profile
#ifdef PROFILING
    call profile_close(model%profile)
#endif

  end subroutine glide_finalise
end module glide_stop
