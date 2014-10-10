!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_profile.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
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

  !TODO - Are we still supporting profiling using the glide_profile module?

module glide_profile

  ! profiling for glide

  implicit none

contains

  subroutine glide_prof_init(model)

    ! initialise glide profiling
    use profile
    use glide_types
    implicit none

    type(glide_global_type) :: model        !> model instance

    if (model%profile%profile_unit == 0) then
       call profile_init(model%profile,'glide.profile')
#if (defined PROFILING && ! defined CCSMCOUPLED && ! defined CESMTIMERS)
       write(model%profile%profile_unit,*) '# take a profile every ',model%numerics%profile_period,' time steps'
#endif
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
    !> start logging profile
    use profile
    use glide_types
    implicit none
    type(glide_global_type) :: model        !> model instance
    integer, intent(in)     :: profn        !> profile number

    call profile_start(model%profile,profn)
  end subroutine glide_prof_start

  subroutine glide_prof_stop(model,profn)
    !> write message to profile
    use profile
    use glide_types
    implicit none
    type(glide_global_type) :: model        !> model instance
    integer, intent(in)     :: profn        !> profile number
    
    !local variables
    character (len=20) :: timestring

    call profile_stop(model%profile,profn)
    if (mod(model%numerics%timecounter,model%numerics%profile_period)==0) then
       write(timestring,*) real(model%numerics%time)
       call profile_log(model%profile,profn,trim(timestring))
    end if
  end subroutine glide_prof_stop
end module glide_profile
