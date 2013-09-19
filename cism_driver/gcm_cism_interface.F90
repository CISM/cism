!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   cism_gcm_interface.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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


module gcm_cism_interface

  use parallel
  use glimmer_commandline
  use glide
  use cism_front_end

  integer, parameter :: GCM_DATA_MODEL = 0  
  integer, parameter :: GCM_CESM = 1 

contains

subroutine gcm_init_interface(which_gcm,gcm_model,cism_model)
  use parallel
  use glimmer_commandline
  use glide
  use cism_front_end 

  integer, intent(in) :: which_gcm
  integer :: gcm_model   ! temporary decl
  type(glide_global_type) :: cism_model    ! CISM  model instance 


  call parallel_initialise
  select case (which_gcm)
    case (GCM_DATA_MODEL)
      call glimmer_GetCommandline()
    case (GCM_CESM)
      ! call gcm_glimmer_GetCommandline_proxy()
    case default
  end select

  call cism_init_dycore(cism_model)

end subroutine gcm_init_interface   


subroutine gcm_run_model(gcm_model,cism_model)
  integer :: gcm_model   ! temporary decl
  type(glide_global_type) :: cism_model    ! CISM  model instance 

  integer :: finished = 0

  do
    ! call gcm_update_model(gcm_model,cism_model)
    ! print *,"In gcm_run_model, calling cism_run_dycore"
    call cism_run_dycore(cism_model)
    ! call gcm_finished(gcm_model,finished)
    if (gcm_finished(gcm_model) .EQ. 1) exit
  end do
end subroutine gcm_run_model


function gcm_finished(gcm_model) result(finished)
  integer :: gcm_model   ! temporary decl
  integer :: finished
 
  ! do something with gcm_model
  finished = 1
end function gcm_finished
  

subroutine gcm_finalize_interface(gcm_model,cism_model)
  integer :: gcm_model   ! temporary decl
  type(glide_global_type) :: cism_model    ! CISM  model instance 

end subroutine gcm_finalize_interface


end module gcm_cism_interface
