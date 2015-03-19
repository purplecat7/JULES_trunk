#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/control/standalone/jules.F90 $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-06-21 10:27:57 +0100 (Thu, 21 Jun 2012) $
!!   $LastChangedRevision: 324 $
!!
!!****************************************************************************

PROGRAM jules
        
  USE time_varying_input_mod, ONLY :                                          &
    update_prescribed_variables => update_model_variables,                    &
    input_close_all => close_all

  USE update_mod, ONLY : update_derived_variables
  
  USE output_mod, ONLY : sample_data, output_data,                            &
                         output_close_all => close_all
  
  USE model_time_mod, ONLY : timestep, start_of_year, end_of_year, end_of_run, is_spinup

  USE switches, ONLY : l_imogen
  
#if defined(MPI_DEFINED)
  USE empi_comms_mod, ONLY : empi_initialise, empi_send_state_info, &
                            empi_send_obs_info, empi_do_filter, empi_cleanup, &
                                              empi_testing
#endif
  
  INTEGER::count
  LOGICAL ::b_first_data=.TRUE.

#if defined(MPI_DEFINED)  
!-----------------------------------------------------------------------------
! Initialise the parallel processing capability
!-----------------------------------------------------------------------------
  CALL empi_initialise()
  WRITE (*, *) "MPI initialised"
#else
  WRITE (*,*) "MPI NOT initialised"
#endif

!-----------------------------------------------------------------------------
! Initialise the model
!-----------------------------------------------------------------------------
  CALL init()
  
#if defined(MPI_DEFINED)
!-----------------------------------------------------------------------------
! Do the first send/receive calls to set up the buffer size information
! and send the number of observations.
! May wish to add sending more meta-data here eg timesteps or timeout info
!-----------------------------------------------------------------------------
  CALL empi_send_state_info()
  CALL empi_send_obs_info()
#endif
  
!-----------------------------------------------------------------------------
! Loop over timesteps.
! Note that the number of timesteps is of unknown length at the start of run,
! if the model is to determine when it has spun up.
!-----------------------------------------------------------------------------
  DO    !  timestep
#if defined(MPI_DEFINED)
        !are we finished spinup?
        !Before the model gets going SEND state vectors to the DA code, and only do it ONCE, 
        !this is effectively the send/recv before the first timestep
        IF ((b_first_data .EQV. .TRUE.) .AND. (is_spinup .EQV. .FALSE.)) THEN
            CALL empi_do_filter()
            b_first_data = .FALSE.
        END IF
#endif
!-----------------------------------------------------------------------------
! Update the IMOGEN climate variables if required
!-----------------------------------------------------------------------------
    IF ( l_imogen .AND. start_of_year ) CALL imogen_update_clim()

!-----------------------------------------------------------------------------
! The update of prescribed data is done in two phases
!  - Update variables provided by files
!  - Update variables that are derived from those updated in the first phase
!-----------------------------------------------------------------------------
    CALL update_prescribed_variables()
    CALL update_derived_variables()

!-----------------------------------------------------------------------------
! Call the main model science routine
!-----------------------------------------------------------------------------
    CALL control(timestep)

!-----------------------------------------------------------------------------
! Update IMOGEN carbon if required
!-----------------------------------------------------------------------------
    IF ( l_imogen .AND. end_of_year ) CALL imogen_update_carb()

!-----------------------------------------------------------------------------
! The gathering of data for output and the actual outputting of data are
! done in two different phases
!-----------------------------------------------------------------------------
    CALL sample_data()
    CALL output_data()

    
#if defined(MPI_DEFINED)
!-----------------------------------------------------------------------------
! Apply the particle filter to prognostic variables
!-----------------------------------------------------------------------------
    !have we finished spin-up? only call if so
    IF (is_spinup .EQV. .FALSE.) THEN
        CALL empi_do_filter()
    END IF
#endif    
    
!-----------------------------------------------------------------------------
! Move the model on to the next timestep
!-----------------------------------------------------------------------------
    CALL next_time()

    IF ( end_of_run ) THEN
        EXIT
    END IF
    
  ENDDO  !  timestep loop

!-----------------------------------------------------------------------------
! Clean up by closing all open files and closing MPI gracefully
!-----------------------------------------------------------------------------
  CALL input_close_all()
  CALL output_close_all()
 
#if defined(MPI_DEFINED)  
  PRINT *, 'ending model, cleaning MPI.'
  CALL empi_cleanup()
#endif
  
  WRITE(*,"(/,a)")'End of run.'

END PROGRAM jules
#endif
