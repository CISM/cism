! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  glide_mask.f90 - part of the Glimmer-CISM ice model      + 
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

#define MASK model%geometry%thkmask

module glide_mask
  !*FD masking ice thicknesses


  integer, parameter :: glide_mask_ocean          = -2
  integer, parameter :: glide_mask_land           = -1
  integer, parameter :: glide_mask_boundary       = 0
  integer, parameter :: glide_mask_thin_ice       = 1
  integer, parameter :: glide_mask_interior       = 2
  integer, parameter :: glide_mask_shelf          = 4
  integer, parameter :: glide_mask_stream         = 8
  integer, parameter :: glide_mask_grounding_line = 16
  integer, parameter :: glide_mask_stream_margin  = 32
  integer, parameter :: glide_mask_land_margin    = 64
  integer, parameter :: glide_mask_shelf_front    = 128
  integer, parameter :: glide_mask_marine_edge    = 256

contains
  subroutine glide_set_mask(model)
    use glide_types
    use glimmer_global, only : dp
    use glimmer_physcon, only : rhoi, rhoo
    implicit none
    type(glide_global_type) :: model        !*FD model instance

    ! local variables
    integer ew,ns
    real(dp), parameter :: con = - rhoi / rhoo

    MASK = 0
    model%geometry%iarea = 0.
    model%geometry%ivol = 0.
    do ns=1,model%general%nsn
       do ew = 1,model%general%ewn
          
          if (model%geometry%thck(ew,ns) .eq. 0.) then                               ! no ice
             if (model%geometry%topg(ew,ns) .lt. model%climate%eus) then             ! below SL
                MASK(ew,ns) = glide_mask_ocean
             else                                                                    ! above SL
                MASK(ew,ns) = glide_mask_land
             end if
          else
             model%geometry%iarea = model%geometry%iarea + 1.
             model%geometry%ivol = model%geometry%ivol + model%geometry%thck(ew,ns)
             if (model%geometry%topg(ew,ns) - model%climate%eus &                    ! ice
                  < con * model%geometry%thck(ew,ns)) then                           ! floating ice
                MASK(ew,ns) = glide_mask_shelf
             else                                                                    ! grounded ice
                MASK(ew,ns) = glide_mask_interior
             end if
             if (model%geometry%thck(ew,ns) .le. model%numerics%thklim) then         ! ice below dynamic limit
                MASK(ew,ns) = ior(MASK(ew,ns),glide_mask_thin_ice)
             end if
          end if

       end do
    end do
    model%geometry%iarea = model%geometry%iarea * model%numerics%dew * model%numerics%dns
    model%geometry%ivol = model%geometry%ivol * model%numerics%dew * model%numerics%dns

    ! finding boundaries
    do ns=2,model%general%nsn-1
       do ew = 2,model%general%ewn-1
          if (is_float(MASK(ew,ns))) then
             ! shelf front
             if (is_ocean(MASK(ew-1,ns)) .or. is_ocean(MASK(ew+1,ns)) .or. &
                  is_ocean(MASK(ew,ns-1)) .or. is_ocean(MASK(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),glide_mask_shelf_front)
             end if
          else if (is_ground(MASK(ew,ns))) then
             ! land margin
             if (is_land(MASK(ew-1,ns)) .or. is_land(MASK(ew+1,ns)) .or. &
                  is_land(MASK(ew,ns-1)) .or. is_land(MASK(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),glide_mask_land_margin)
             end if
             ! grounding line
             if (is_float(MASK(ew-1,ns)) .or. is_float(MASK(ew+1,ns)) .or. &
                  is_float(MASK(ew,ns-1)) .or. is_float(MASK(ew,ns+1))) then
                MASK(ew,ns) = ior(MASK(ew,ns),glide_mask_grounding_line)
             end if
          end if
          ! Edge of marine ice, whether floating or not
          if ((model%geometry%topg(ew,ns) .lt. model%climate%eus.and.&
               model%geometry%thck(ew,ns)>0.0).and. &
               (is_ocean(MASK(ew-1,ns)) .or. is_ocean(MASK(ew+1,ns)) .or. &
               is_ocean(MASK(ew,ns-1)) .or. is_ocean(MASK(ew,ns+1)))) then
             MASK(ew,ns) = ior(MASK(ew,ns),glide_mask_marine_edge)
          end if
       end do
    end do
  end subroutine glide_set_mask

  logical elemental function is_ocean(mask)
    !*FD returns .true. if node is ocean
    implicit none
    integer, intent(in) :: mask 

    is_ocean = mask.eq.glide_mask_ocean
  end function is_ocean

  logical elemental function is_land(mask)
    !*FD returns .true. if node is land
    implicit none
    integer, intent(in) :: mask 

    is_land = mask.eq.glide_mask_land
  end function is_land

  logical elemental function has_ice(mask)
    !*FD returns .true. if node contains ice
    implicit none
    integer, intent(in) :: mask 

    has_ice = mask .gt. 0
  end function has_ice

  logical elemental function is_thin(mask)
    !*FD returns .true. if node is below dynamic limit
    implicit none
    integer, intent(in) :: mask

    is_thin = (iand(mask,glide_mask_thin_ice) .gt. 0 .and. mask.gt.0)
  end function is_thin

  logical elemental function is_float(mask)
    !*FD returns .true. if node is floating
    implicit none
    integer, intent(in) :: mask

    is_float = (iand(mask,glide_mask_shelf) .gt. 0 .and. mask.gt.0)
  end function is_float

  logical elemental function is_ground(mask)
    !*FD returns .true. if node is grounded ice
    implicit none
    integer, intent(in) :: mask

    is_ground = (iand(mask,glide_mask_interior) .gt. 0 .and. mask.gt.0)
  end function is_ground

  logical elemental function is_calving(mask)
    !*FD return .true. if node is at the shelf front
    implicit none
    integer, intent(in) :: mask

    is_calving = (iand(mask,glide_mask_shelf_front) .gt. 0 .and. mask.gt.0)
  end function is_calving

  logical elemental function is_marine_ice_edge(mask)
    !*FD return .true. if node is at edge of ice and topgraphy is
    !*FD below sea-level
    implicit none
    integer, intent(in) :: mask

    is_marine_ice_edge = (iand(mask,glide_mask_marine_edge) .gt. 0 .and. mask.gt.0)
  end function is_marine_ice_edge
  
  logical elemental function is_grounding_line(mask)
    !*FD returns .true. if node is grounding line
    implicit none
    integer, intent(in) :: mask 

    is_grounding_line = mask.eq.glide_mask_grounding_line
  end function is_grounding_line
end module glide_mask
