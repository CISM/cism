! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_profile.f90 - part of the GLIMMER ice model        + 
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

module glide_profile
  !*FD profiling for glide

contains
  subroutine glide_prof_init(model)
    !*FD initialise glide profiling
    use profile
    use glide_types
    implicit none
    type(glide_global_type) :: model        !*FD model instance

    if (model%profile%profile_unit .eq. 0) then
       call profile_init(model%profile,'glide.profile')
       write(model%profile%profile_unit,*) '# take a profile every ',model%numerics%profile_period,' time steps'
    end if

    ! registering glide profiles
    model%glide_prof%geomderv    = profile_register(model%profile,'horizontal derivatives')
    model%glide_prof%hvelos      = profile_register(model%profile,'horizontal velocities')
    model%glide_prof%ice_mask1   = profile_register(model%profile,'ice mask 1')
    model%glide_prof%temperature = profile_register(model%profile,'temperature')
    model%glide_prof%ice_evo     = profile_register(model%profile,'ice evolution')
    model%glide_prof%ice_mask2   = profile_register(model%profile,'ice mask 2')
    model%glide_prof%isos_water  = profile_register(model%profile,'isostasy water')
    model%glide_prof%isos        = profile_register(model%profile,'isostasy')
  end subroutine glide_prof_init
  
  subroutine glide_prof_start(model,profn)
    !*FD start logging profile
    use profile
    use glide_types
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    integer, intent(in)     :: profn        !*FD profile number

    call profile_start(model%profile,profn)
  end subroutine glide_prof_start

  subroutine glide_prof_stop(model,profn)
    !*FD write message to profile
    use profile
    use glide_types
    implicit none
    type(glide_global_type) :: model        !*FD model instance
    integer, intent(in)     :: profn        !*FD profile number
    
    !local variables
    character (len=20) :: timestring

    call profile_stop(model%profile,profn)
    if (mod(model%numerics%timecounter,model%numerics%profile_period).eq.0) then
       write(timestring,*) model%numerics%time
       call profile_log(model%profile,profn,trim(timestring))
    end if
  end subroutine glide_prof_stop
end module glide_profile
