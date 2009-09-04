! ******************************************************************************
! verif_glide.f90
! Magnus Hagdorn
!
! testing steady state
! ******************************************************************************
!

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program verifglide

  ! load various modules
  use glimmer_global, only:rk ! precision of the model
  use glide                   ! main glide module
  use glimmer_log             ! module for logging messages
  use glimmer_config          ! module for handling configuration files
  use verif
  use verif_io
  implicit none

  ! some variables
  type(glide_global_type) :: model    
  type(ConfigSection), pointer :: config

  type(verif_type) :: veri

  character(len=50) :: fname  ! name of paramter file
  real(kind=rk) time ! current time
  
  ! Ask for configuration file
  write(*,*) 'Enter name of GLIDE configuration file to be read'
  read(*,*) fname
  
  ! start logging
  call open_log(unit=50)
  
  ! read configuration
  call ConfigRead(fname,config)
  call glide_config(model,config)
  call verif_config(config,veri)
  call verif_printconfig(veri)
  call CheckSections(config)

  ! initialise test setup
  call verif_init(model, veri)
  
  ! initialise GLIDE
  call glide_initialise(model)
  ! fill dimension variables
  ! create verif variables
  call verif_io_createall(model)
  call glide_nc_fillall(model)
  ! get current time from start time
  time = get_tstart(model)

  ! initial conditions
  call verif_update(model, veri, time)
  call verif_initthk(model, veri)

  ! loop over times
  do while(time.le.model%numerics%tend)
     call verif_update(model, veri, time)

     ! calculate temperature and velocity distribution
     call glide_tstep_p1(model,time)
     call verif_io_writeall(veri,model)
     ! write to netCDF file, move ice
     call glide_tstep_p2(model)
     ! calculate isostatic adjustment
     call glide_tstep_p3(model)
     ! increment time counter
     time = time + get_tinc(model)
  end do

  ! finalise GLIDE
  call glide_finalise(model)
end program verifglide
