! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! +                                                           +
! +  simple_glide.f90 - part of the Glimmer-CISM ice model    + 
! +                                                           +
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! 
! Copyright (C) 2005, 2006, 2007, 2008, 2009
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

program simple_glide
  !*FD This is a simple GLIDE test driver. It can be used to run
  !*FD the EISMINT test cases
  use parallel
  use glimmer_global, only:rk
  use glide
  use simple_forcing
  use glimmer_log
  use glimmer_config
  use glimmer_commandline
  use glimmer_writestats
  use glimmer_filenames, only : filenames_init

  use glide_diagnostics

  implicit none

#ifdef GPTL
#include <gptl.inc>
#endif
#ifdef PAPI
#include <f90papi.h>
#endif

  type(glide_global_type) :: model        ! model instance
  type(simple_climate) :: climate         ! climate
  type(ConfigSection), pointer :: config  ! configuration stuff
  real(kind=rk) time
  real(kind=dp) t1,t2
  integer clock,clock_rate,ret

  integer :: tstep_count

  call parallel_initialise

  ! start gptl
#ifdef GPTL
  ret = gptlsetoption (gptlprint_method,gptlfull_tree)
  ret = gptlsetoption (PAPI_FP_OPS, 1)
  ret = gptlsetutr (gptlnanotime)
  ret = gptlinitialize ()
  ret = gptlstart ('total')
#endif

  call glimmer_GetCommandline()


  ! start logging
  call open_log(unit=50, fname=logname(commandline_configname))
  
  ! setup paths
  call filenames_init(commandline_configname)

  ! read configuration
  call ConfigRead(commandline_configname,config)

  ! start timing
  call system_clock(clock,clock_rate)
  t1 = real(clock,kind=dp)/real(clock_rate,kind=dp)

  ! initialise GLIDE
  call glide_config(model,config)
  call simple_initialise(climate,config)
  call glide_initialise(model)
  call CheckSections(config)
  ! fill dimension variables
  call glide_nc_fillall(model)
  time = model%numerics%tstart

  call simple_massbalance(climate,model,time)
  call simple_surftemp(climate,model,time)
  call spinup_lithot(model)

  tstep_count = 0

  do while(time.le.model%numerics%tend)
     call glide_tstep_p1(model,time)
     call glide_tstep_p2(model)
     call glide_tstep_p3(model)
     ! override masking stuff for now

     tstep_count = tstep_count + 1
#ifdef JEFFORIG
!JEFF commenting out during the first pass through.  There's several direct model references and file writes in glide_write_diag.
     if (mod(tstep_count, model%numerics%ndiag) == 0) then
        call glide_write_diag(model, time, model%numerics%idiag, &
                                           model%numerics%jdiag )
     endif
#endif

     ! Redistribute calls here to spread the data back out.
     call distributed_scatter_var(model%stress%efvs, gathered_efvs)
     call distributed_scatter_var(model%velocity%uvel, gathered_uvel)
     call distributed_scatter_var(model%velocity%vvel, gathered_vvel)
     call distributed_scatter_var(model%velocity%uflx, gathered_uflx)
     call distributed_scatter_var(model%velocity%vflx, gathered_vflx)
     call distributed_scatter_var(model%velocity%velnorm, gathered_velnorm)
     call distributed_scatter_var(model%geometry%thck, gathered_thck)
     call distributed_scatter_var(model%geomderv%stagthck, gathered_stagthck)
     call distributed_scatter_var(model%climate%acab, gathered_acab)
     call distributed_scatter_var(model%temper%temp, gathered_temp)
     call distributed_scatter_var(model%geomderv%dusrfdew, gathered_dusrfdew)
     call distributed_scatter_var(model%geomderv%dusrfdns, gathered_dusrfdns)
     call distributed_scatter_var(model%geomderv%dthckdew, gathered_dthckdew)
     call distributed_scatter_var(model%geomderv%dthckdns, gathered_dthckdns)
     call distributed_scatter_var(model%stress%tau%xx, gathered_tauxx)
     call distributed_scatter_var(model%stress%tau%yy, gathered_tauyy)
     call distributed_scatter_var(model%stress%tau%xy, gathered_tauxy)
     call distributed_scatter_var(model%stress%tau%scalar, gathered_tauscalar)
     call distributed_scatter_var(model%stress%tau%xz, gathered_tauxz)
     call distributed_scatter_var(model%stress%tau%yz, gathered_tauyz)
     call distributed_scatter_var(model%geometry%topg, gathered_topg)
     call distributed_scatter_var(model%geometry%thkmask, gathered_thkmask)
     call distributed_scatter_var(model%geometry%marine_bc_normal, gathered_marine_bc_normal)
     call distributed_scatter_var(model%velocity%surfvel, gathered_surfvel)
     call distributed_scatter_var(model%ground%gline_flux, gathered_gline_flux)
     call distributed_scatter_var(model%velocity%ubas, gathered_ubas)
     call distributed_scatter_var(model%velocity%vbas, gathered_vbas)
     call distributed_scatter_var(model%isos%relx, gathered_relx)
     call distributed_scatter_var(model%temper%flwa, gathered_flwa)
     call distributed_scatter_var(model%climate%calving, gathered_calving)
     call distributed_scatter_var(model%climate%backstress, gathered_backstress)
     call distributed_scatter_var(model%geometry%usrf, gathered_usrf)
     call distributed_scatter_var(model%climate%backstressmap, gathered_backstressmap)
     call distributed_scatter_var(model%stress%tau_x, gathered_tau_x)
     call distributed_scatter_var(model%stress%tau_y, gathered_tau_y)
     call distributed_scatter_var(model%geometry%lsrf, gathered_lsrf)

     !After scattering, then update nsn and ewn back to local values
     model%general%ewn = local_ewn
     model%general%nsn = local_nsn

     ! Redistribute calls here to spread the data back out.
     call parallel_halo(model%stress%efvs)
     call parallel_halo(model%velocity%uvel)
     call parallel_halo(model%velocity%vvel)
     call parallel_halo(model%velocity%uflx)
     call parallel_halo(model%velocity%vflx)
     call parallel_halo(model%velocity%velnorm)
     call parallel_halo(model%geometry%thck)
     call parallel_halo(model%geomderv%stagthck)
     call parallel_halo(model%climate%acab)
     call parallel_halo(model%geomderv%dusrfdew)
     call parallel_halo(model%geomderv%dusrfdns)
     call parallel_halo(model%geomderv%dthckdew)
     call parallel_halo(model%geomderv%dthckdns)
     call parallel_halo(model%stress%tau%xx)
     call parallel_halo(model%stress%tau%yy)
     call parallel_halo(model%stress%tau%xy)
     call parallel_halo(model%stress%tau%scalar)
     call parallel_halo(model%stress%tau%xz)
     call parallel_halo(model%stress%tau%yz)
     call parallel_halo(model%geometry%topg)
     call parallel_halo(model%geometry%thkmask)
     call parallel_halo(model%geometry%marine_bc_normal)
     call parallel_halo(model%velocity%surfvel)
     call parallel_halo(model%ground%gline_flux)
     call parallel_halo(model%velocity%ubas)
     call parallel_halo(model%velocity%vbas)
     call parallel_halo(model%isos%relx)
     call parallel_halo(model%temper%flwa)
     call parallel_halo(model%climate%calving)
     call parallel_halo(model%climate%backstress)
     call parallel_halo(model%geometry%usrf)
     call parallel_halo(model%climate%backstressmap)
     call parallel_halo(model%stress%tau_x)
     call parallel_halo(model%stress%tau_y)
     call parallel_halo(model%geometry%lsrf)

     call parallel_halo_temperature(model%temper%temp)

     ! Perform parallel operations for restart files
     call glide_tstep_postp3(model)

     time = time + model%numerics%tinc
     call simple_massbalance(climate,model,time)
     call simple_surftemp(climate,model,time)     
  end do

  ! finalise GLIDE
  call glide_finalise(model)
  call system_clock(clock,clock_rate)
  t2 = real(clock,kind=dp)/real(clock_rate,kind=dp)
  call glimmer_write_stats(commandline_resultsname,commandline_configname,t2-t1)
  call close_log
  call parallel_finalise

  ! stop gptl
#ifdef GPTL
  ret = gptlstop ('total')
  ret = gptlpr (0)
  ret = gptlfinalize ()
#endif
  
end program simple_glide
