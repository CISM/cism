

program simple_phaml
    use glimmer_global, only:rk ! precision of the model
    use glide                   ! main glide module
    use glimmer_log             ! module for logging messages
    use glimmer_config          ! module for handling configuration files
    use glide_mask
    
    use phaml
    use phaml_example
    use phaml_user_mod
    implicit none
    type(phaml_solution_type) :: phaml_solution
    type(glide_global_type) :: cism_model
    type(ConfigSection), pointer :: config
    character(len=50) :: fname 
    character(len=255) :: cwd
    real(kind=rk) time
    !fname = 'sm.conf'
    write(*,*) 'Enter name of GLIDE configuration file to be read'
    read(*,*) fname
    call getcwd(cwd)
    
    cwd = trim(cwd)//'/'//fname
    ! start logging
    call open_log(unit=50)
    ! read configuration
    call ConfigRead(cwd,config)
    call CheckSections(config)
    call glide_config(cism_model,config)
    ! initialise GLIDE
    call glide_initialise(cism_model)
    ! fill dimension variables
    call glide_nc_fillall(cism_model)
    ! get current time from start time
    time = get_tstart(cism_model)
    !mask set in initialise    
    call glide_set_mask(cism_model)  

    !initialize all variables needed
    call phaml_init(cism_model,phaml_solution)
    !does the evaluation and places the
    !solution in cism_model%phaml%uphaml
    call phaml_evolve(cism_model,phaml_solution)
    !close and free variables
    call phaml_close(phaml_solution)
    
    call glide_io_writeall(cism_model,cism_model)
    time = time + get_tinc(cism_model)
    cism_model%numerics%time = time
    
    call glide_finalise(cism_model)
    
end program simple_phaml
