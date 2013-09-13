!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   felix_dycore_interface.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


module felix_dycore_interface

   use glide_types
   use glimmer_log
   !use glissade_velo_higher  ! I assume you probably will need this one (and maybe others).
   !use glimmer_to_dycore

   implicit none
   private


   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: felix_velo_init, &
             felix_velo_driver

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************


contains


!***********************************************************************
!
!  routine felix_velo_init
!
!> \brief   Initializes the external Albany-FELIX velocity solver
!> \author  Irina Kalashnikova
!> \date    13 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine initializes the external Albany-FELIX ice velocity solver.
!
!-----------------------------------------------------------------------

   subroutine felix_velo_init(model)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(glide_global_type),intent(inout) :: model

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------


      print *, 'DEBUG: Inside felix_velo_init.'

      ! === First do any preparations needed on the CISM side (if any)


      ! === Now call the external Albany code for any init that it needs to do
      !call gtd_init_dycore(model,dycore_model_index)
      ! Doug - does this interface still make sense here?
      ! Doug - what needs to change (if anything) if the code is compiled without
      !        external Felix libraries?  Do we need a stub module?
      ! Doug - We might need to do some rearranging to make sure the call 
      !        to gtd_init_dycore_interface happens in the right place 
      !        (presumably in simple_glide/simple_felix/cism_driver).
      !        (I think I see how to do this, but will wait for now.)

   !--------------------------------------------------------------------
   end subroutine felix_velo_init




!***********************************************************************
!
!  routine felix_velo_driver
!
!> \brief   Makes preparations and calls the external Albany-FELIX velocity solver
!> \author  Irina Kalashnikova
!> \date    13 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine makes preparations and calls the external
!>  Albany-FELIX velocity solver.
!
!-----------------------------------------------------------------------

   subroutine felix_velo_driver(model)

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type(glide_global_type),intent(inout) :: model

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------


      print *, 'DEBUG: Inside felix_velo_driver.'

      ! === First do any preparations needed on the CISM side
      !call glissade_velo_higher_data(...)  or maybe this should be get_parallel_mesh_data?  ! probably makes sense to change the name


      ! === Now call the external Albany code
      !call gtd_run_dycore(dycore_model_index,cur_time,time_inc)
      ! Doug - does this interface still make sense here?
      ! Doug - what needs to change (if anything) if the code is compiled without
      !        external Felix libraries?  Do we need a stub module?


   !--------------------------------------------------------------------
   end subroutine felix_velo_driver



!***********************************************************************
! Private subroutines:
!***********************************************************************



!***********************************************************************
!
!  routine glissade_velo_higher_data  OR get_parallel_mesh_data  OR ???
!
!> \brief   ?
!> \author  Irina Kalashnikova
!> \date    13 September 2013
!> \version SVN:$Id$
!> \details
!>  This routine ...
!
!-----------------------------------------------------------------------

   subroutine glissade_velo_higher_data()
! Probably makes sense to change the name

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------



      ! your existing code here


   !--------------------------------------------------------------------
   end subroutine glissade_velo_higher_data


end module felix_dycore_interface
