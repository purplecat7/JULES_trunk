


PROGRAM empire

    USE empire_mod, ONLY: empi_pf_initialise, empi_pf_model_hookup_state_dims, &
    empi_pf_model_hookup_obs_size, empi_pf_assign_models, &
    empi_pf_do_work, empi_pf_cleanup
    USE pf_control, ONLY: set_pf_controls, allocate_pf, deallocate_pf   !these are DA ones - need a big tidy-up
    
    LOGICAL :: b_obssize_known = .FALSE.

    !set up the MPI runtime environment, get a channel to the models
    PRINT*, 'mpi_pf_initialise____________'
    CALL empi_pf_initialise()
    
    !distribute the model processes amongst the available filter processes
    !also sets some stuff in the pf_control structure
    PRINT*, 'mpi_pf_assign_models____________'
    CALL empi_pf_assign_models()
    
    PRINT*, 'read setup files___________'
    CALL set_pf_controls
    
    !receive information from the model (state vector sizes, time steps? etc.)
    !used to set a couple of module variables in 'sizes.f90' which get used all over the shop!
    !now moved to pf_control at least so they're more identifiable and grouped
    PRINT*, 'mpi_pf_model_hookup____________'
    CALL empi_pf_model_hookup_state_dims()
    !this next function returns a switch to decide how to do filtering
    !ie known time steps or not
    b_obssize_known = empi_pf_model_hookup_obs_size()
    
    !the first of many things which need the numbers from 'sizes.f90'
    PRINT*, 'allocate space for filtering code___________'
    CALL allocate_pf
    
    PRINT*, 'do filtering work___________'
    CALL empi_pf_do_work(b_obssize_known)
    
    !tidy up
    PRINT*, 'mpi_pf_cleanup____________'
    CALL empi_pf_cleanup()
    CALL deallocate_pf
    
    PRINT*, 'DONE!'

END PROGRAM empire



