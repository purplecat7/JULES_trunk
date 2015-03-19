


PROGRAM pf_minimal

    USE pf_mpi_mod, ONLY: empi_pf_initialise, empi_pf_model_hookup, empi_pf_assign_models, &
    empi_pf_do_filter, empi_pf_cleanup
  
    !set up the MPI runtime environment, get a channel to the models
    PRINT*, 'mpi_pf_initialise____________'
    CALL empi_pf_initialise()
    
    !distribute the model processes amongst the available filter processes
    PRINT*, 'mpi_pf_assign_models____________'
    CALL empi_pf_assign_models()
    
    !receive information from the model (state vector sizes, time steps? etc.)
    PRINT*, 'mpi_pf_model_hookup____________'
    CALL empi_pf_model_hookup()
    
    !execute the data assimilation filters and other magic stuff
    PRINT*, 'mpi_pf_do_filter____________'
    CALL empi_pf_do_filter()   
    
    !tidy up
    PRINT*, 'mpi_pf_cleanup____________'
    CALL empi_pf_cleanup()
    
    PRINT*, 'DONE!'

END PROGRAM pf_minimal



