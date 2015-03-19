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
  
  USE model_time_mod, ONLY : timestep, start_of_year, end_of_year, end_of_run

  USE switches, ONLY : l_imogen
  
  USE logging_mod, ONLY : log_init, log_info, log_debug, log_warn, log_error, log_fatal
  
  IMPLICIT NONE
  
  INCLUDE 'mpif.h'
  INTEGER mpi_err, mpi_numtasks, mpi_rank, rc, mpi_len
  INTEGER:: test_total=0
  INTEGER open_err
  CHARACTER (LEN=50):: mpi_msg, my_fmt
  CHARACTER*(MPI_MAX_PROCESSOR_NAME) mpi_name
  !variables for the mpi Send/Receive
  INTEGER, PARAMETER:: c_mpi_zero=0, c_mpi_tag=1 !constants
  INTEGER:: buffer=0, count=1, datatype=MPI_INTEGER, source=MPI_ANY_SOURCE, mpi_comm_handle=1, status(MPI_STATUS_SIZE)
  
  !OPEN (20, FILE='/home/xw904346/SatWin/jules/harvard/output/a_log_messages.txt', IOSTAT=open_err)
  !IF (open_err .eq. 0) THEN
    !CALL log_init(20, 31, 12) !args: file handle, print everything, stop for Errors and Warnings.
  !ENDIF
  CALL log_init(print_level=31, stop_level=8) !args: print everything, stop for Errors.

  !Set up an MPI environment to run the model in parallel, check for success
  !Requires a configuration file specifying servers and threads
  CALL MPI_INIT(mpi_err)
  
  IF (mpi_err .ne. MPI_SUCCESS) then
    PRINT *,'Error starting MPI parallelisation program. Terminating.'
    CALL log_error("Main",                                                    &
                   "Error starting MPI parallelisation program. Terminating MPI.")
    CALL MPI_ABORT(MPI_COMM_WORLD, rc, mpi_err)
    
  ELSE
    !Get the size of the group (mpi_numtasks) associated with the communicator (MPI_COMM_WORLD)
    CALL MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_numtasks, mpi_err)
    WRITE (mpi_msg, "(A20, I2)") "MPI process total: ", mpi_numtasks
    CALL log_info("Main", mpi_msg)
    
  END IF
  !Get the unique name of the processor node (and string length)
  CALL MPI_GET_PROCESSOR_NAME(mpi_name, mpi_len, mpi_err)
  !Get the unique id of the process (its rank) on this communicator - each process will
  !print out this onformation but the rank will be different
  CALL MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_err)
  !create the output format on the fly knowing the string length
  WRITE (my_fmt, "(A23, I2, A1)") '(A12, I3, A4, I3, A4, A' ,mpi_len, ')'
  WRITE (mpi_msg, my_fmt) "MPI process ",mpi_rank, " of ", mpi_numtasks, " on ", mpi_name
  !sadly, this notation doesn't work on the compiler...
  !WRITE (mpi_msg, "(A12, I3, A4, I3, A4, <mpi_len>A1)") "MPI process ",mpi_rank, " of ", mpi_numtasks, " on ", mpi_name

!-----------------------------------------------------------------------------
! Initialise the model
!-----------------------------------------------------------------------------
  CALL init()

!-----------------------------------------------------------------------------
! Loop over timesteps.
! Note that the number of timesteps is of unknown length at the start of run,
! if the model is to determine when it has spun up.
!-----------------------------------------------------------------------------
  DO    !  timestep

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

!-----------------------------------------------------------------------------
! Move the model on to the next timestep
!-----------------------------------------------------------------------------
    CALL next_time()

    IF ( end_of_run ) EXIT
    
!-----------------------------------------------------------------------------
! TODO Apply the particle filter to prognostic variables
!-----------------------------------------------------------------------------
    IF (mpi_err .eq. MPI_SUCCESS) THEN
      !Create info string
      WRITE (mpi_msg, "(A15, I2)") "From MPI rank: ", mpi_rank

      IF (mpi_rank .eq. 0) THEN
        !If this instance is the primary one, then get all the other's variables
        !CALL log_info("Main", mpi_msg)
        CALL MPI_RECV(buffer, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, mpi_err)
        test_total = test_total + buffer
      ELSE
        !Else, for all other instances, send the variables to the primary instance
        !CALL log_info("Main", mpi_msg)
        !WRITE (*, '(I10)'), buffer
        CALL MPI_SEND(buffer, 1, MPI_INTEGER, c_mpi_zero, c_mpi_tag, MPI_COMM_WORLD, mpi_err)
        buffer = buffer + 1

      ENDIF
    ELSE
      CALL log_warn("Main", &
                    "Model running one instance only, error in MPI process ranking.")
    ENDIF
    
    
  ENDDO  !  timestep loop
  IF (mpi_rank .eq. 0) THEN
    WRITE (mpi_msg, "(A15, I10)") "Test total: ", test_total
  ENDIF
  
  CALL log_info("Main", mpi_msg)
!-----------------------------------------------------------------------------
! Clean up by closing all open files
!-----------------------------------------------------------------------------
  CALL input_close_all()
  CALL output_close_all()
  
  !Tidy up the MPI communicator session in each process.
  CALL MPI_FINALIZE(mpi_err)

  WRITE(*,"(/,a)")'End of run.'

END PROGRAM jules
#endif
