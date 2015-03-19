!     
! File:   empi_comms_mod.F90
! Author: xw904346 (Jane Lewis, Dept of Meteorology, UREAD)
!
! Created on August 11, 2014, 3:22 PM
!
!This file contains all the wrapper functions for the model to use MPI and connected
!data assimilation processes. The model therefore only need to inset the calls in the
!correct places in the main program.
!It is also responsible for communicating the required state vectors.

MODULE empi_comms_mod

    USE logging_mod, ONLY: log_init, log_info, log_debug, log_warn, log_error, log_fatal
    USE prognostics, ONLY: smcl, t_soil, tstar_tile

    IMPLICIT NONE

    INCLUDE 'mpif.h'


    INTEGER     empi_err, &                   !holds error value
                      empi_numtasks, &          !total number of tasks running in the MPI environment
                      world_id, &                 !the id of a task instance running this code
                      empi_errcode, &            !error code used for ABORT function
                      empi_strlen, &               !processor name string size
                      empi_coupling_comms, &     !group communicator for model and particle filter processes
                      empi_total_procs, &        !total processes including model and particle filter
                      empi_model_procs, &     !total model processes
                      empi_coupling_id, &   !id of the model on the group communicator
                      filter_rank_id, &                  !process id of the particle filter this model instance will use
                      empi_request_handle        !for test purposes only
    

    !INTEGER open_err !for writing to file
    CHARACTER (LEN = 50) :: empi_msg
    INTEGER :: test_total = 0

    
    !the size of each of the 3 state vector dimensions
    INTEGER     smcl_size_0, smcl_size_1, &     !for the prognostic variable 'SMCL'
                      tsoil_size_0, tsoil_size_1, &                  !for T_SOIL
                      tstar_tile_size_0, tstar_tile_size_1, &        !for TSTAR_TILE
                      count0, count1
    REAL(KIND = KIND(1.0D0)), ALLOCATABLE, DIMENSION(:) :: state_buffer

CONTAINS

    !-----------------------------------------------------------------------------
    ! Initialise MPI and communicators.
    !
    ! Start the MPI environment and set up the communicators  for model groups and
    ! between data assimilation processes and model group head nodes.
    ! The model heads are distributed evenly amongst available DA processes.
    !
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_initialise()

        IMPLICIT NONE
        INTEGER, PARAMETER :: c_num_mdl_procs = 1 !this is how many processes run to constitute a single model
        INTEGER :: count, npfs, world_size, num_models, m_model_id_rank1, m_model_id_rank2, group_mdl_colour
        INTEGER :: coupling_colour = MPI_UNDEFINED
        INTEGER :: empi_model_comms1, empi_model_comms2
        CHARACTER*(MPI_MAX_PROCESSOR_NAME) empi_name
        CHARACTER (LEN = 50) :: my_fmt
        !This is the 'magic' number for hooking into the correct communicator group for the
        !particle filter
        INTEGER, PARAMETER :: c_coupling_colour = 9999, c_mdl_colour = 0
        
        
        !OPEN (20, FILE='/home/xw904346/SatWin/jules/harvard/output/a_log_messages.txt', IOSTAT=open_err)
        !IF (open_err .eq. 0) THEN
        !CALL log_init(20, 31, 12) !args: file handle, print everything, stop for Errors and Warnings.
        !ENDIF

        CALL log_init(print_level = 31, stop_level = 8) !args: print everything, stop for Errors.

        !Set up an MPI environment to run multiple copies of the model, check for success
        !Requires a configuration file specifying servers and threads
        CALL MPI_INIT(empi_err)

        IF (empi_err .ne. MPI_SUCCESS) then
            PRINT *, 'Error starting MPI parallelisation program. Terminating.'
            CALL log_error("Main", &
            "Error starting MPI parallelisation program. Terminating MPI.")
            CALL MPI_ABORT(MPI_COMM_WORLD, empi_errcode, empi_err)

        ELSE
            !Get the number of processes (empi_numtasks) associated with the communicator (MPI_COMM_WORLD)
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, empi_numtasks, empi_err)
            WRITE (empi_msg, "(A20, I2)") "MPI process total: ", empi_numtasks
            CALL log_info("Main", empi_msg)

        END IF

        
        PRINT*, 'mpi init, getting proc name '
        !Get the unique name of the processor node (and string length)
        CALL MPI_GET_PROCESSOR_NAME(empi_name, empi_strlen, empi_err)
        !Get the unique id of the process (its rank) on this communicator - each process will
        !print out this information but the rank will be different
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, world_id, empi_err)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, empi_err)
        !PRINT*, 'got rank: ', world_id, ' in total: ', world_size
        !create the output format on the fly knowing the string length - write the format to a string using a format! crazy shit :)
        WRITE (my_fmt, "(A23, I2, A1)") '(A12, I3, A4, I3, A4, A', empi_strlen, ')'
        WRITE (empi_msg, my_fmt) "MPI process ", world_id, " of ", empi_numtasks, " on ", empi_name
        !sadly, this notation doesn't work on the compiler...
        !WRITE (empi_msg, "(A12, I3, A4, I3, A4, <mpi_strlen>A1)") "MPI process ",mpi_rank, " of ", empi_numtasks, " on ", empi_name

        PRINT*, empi_msg
        
        !The rather pointless splits... only useful if each model is run as a bunch of processes with a head node
        !This communicator is for all model processes using 'c_mdl_colour', with them ranked using 'world_id'
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, c_mdl_colour, world_id, empi_model_comms1, empi_err)
        !We can then find the total number of model processes...
        CALL MPI_COMM_SIZE(empi_model_comms1, num_models, empi_err)
        !and each one's rank on this communicator
        CALL MPI_COMM_RANK(empi_model_comms1, m_model_id_rank1, empi_err)
        
        !This integer division will give us a number per group of models (assuming that component models are sequential)
        group_mdl_colour = m_model_id_rank1 / c_num_mdl_procs
        
        !Create communicator per group - each process in the same group will link to it - these will be the component parts of a single model
        CALL MPI_COMM_SPLIT(empi_model_comms1, group_mdl_colour, m_model_id_rank1, empi_model_comms2, empi_err)
        CALL MPI_COMM_RANK(empi_model_comms2, m_model_id_rank2, empi_err)
        
        IF (m_model_id_rank2 .EQ. 0) THEN
            !this is a head node in a model group
            coupling_colour = c_coupling_colour
        END IF
        
        !We'll join up with the filtering code's communicator, but only if we're a head node. We'll be ranked by our group number
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, coupling_colour, group_mdl_colour, empi_coupling_comms, empi_err)
        
        IF (m_model_id_rank2 .EQ. 0) THEN
            !If we're the head node with access to the coupling communicator, find the total no. of processes (head models + filters)
            CALL MPI_COMM_SIZE(empi_coupling_comms, empi_total_procs, empi_err)
            !And what our unique rank is on the coupling communicator
            CALL MPI_COMM_RANK(empi_coupling_comms, empi_coupling_id, empi_err)
            !number of filters = total number of processes in MPI - total number of models
            npfs = world_size - num_models
            !number of head model processes = no.processes on coupler - number of filters
            empi_model_procs = empi_total_procs - npfs
            !Distribute the models around the available particle filter processes
            !Knowing our own ranking, we can set the rank of the filter we'll be using
            DO count = 1, npfs
                IF (REAL(empi_coupling_id) .lt. REAL(empi_model_procs * count)/REAL(npfs)) THEN
                    filter_rank_id = count - 1 + empi_model_procs
                    EXIT
                END IF
            END DO
            PRINT*, 'particle filter rank to use: ', filter_rank_id
        ELSE
            filter_rank_id = -1
        END IF
        
        !Lastly, run each model in its own environment
        CALL empi_set_working_dir()

    END SUBROUTINE empi_initialise

    !-----------------------------------------------------------------------------
    ! Simple test for mpi without additional particle filter process
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_testing()
        IMPLICIT NONE
        !variables/constants for the simple test mpi Send/Receive
        !NOTE that these initialised values aquire the SAVE attribute and their contents would
        !therefore be carried over between calls to this function
        !http://stackoverflow.com/questions/3509208/does-fortran-preserve-the-value-of-internal-variables-through-function-and-subro

        INTEGER, PARAMETER :: c_mpi_zero = 0, c_mpi_tag = 1 !constants
        INTEGER :: buffer = 0
        LOGICAL:: interrupt, retval

        IF (empi_err .eq. MPI_SUCCESS) THEN
            !Create info string
            WRITE (empi_msg, "(A15, I2)") "From MPI rank: ", world_id

            IF (world_id .eq. 0) THEN

                !If this instance is the primary one, then get all the other's variables
                !CALL log_info("Main", empi_msg)
                CALL MPI_RECV(buffer, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE, empi_err)
                test_total = test_total + buffer

            ELSE
                !Else, for all other instances, send the variables to the primary instance
                !CALL log_info("Main", empi_msg)
                !WRITE (*, '(I10)'), buffer
                CALL MPI_SEND(buffer, 1, MPI_INTEGER, c_mpi_zero, c_mpi_tag, MPI_COMM_WORLD, empi_err)
                buffer = buffer + 1
                
            ENDIF
        ELSE
            CALL log_warn("Main", &
            "Model running one instance only, error in MPI process ranking.")
        ENDIF

    END SUBROUTINE empi_testing
    
   

    !-----------------------------------------------------------------------------
    ! Initialize the prognostic variable buffer and send its size to particle filter
    ! code so that it can initialize correctly.
    ! SMCL, TSTAR_SIZE and T_SOIL are imported from the prognostics module
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_send_state_info()
        IMPLICIT NONE
        INTEGER err_code, smcl_size, tsoil_size, tstar_tile_size, states, state_dims
        INTEGER, ALLOCATABLE, DIMENSION(:) :: state_dim_arr !an array to hold the state vector dimensions
        INTEGER, PARAMETER :: c_mpi_tag = 1 !constants

        !Set up an array to store the dimensions of the state arrays so that the particle
        !filter knows them
        !****************************************************************************
        !This is currently a 'magic' number but should ideally be read from file along with
        !the prognostic state variables to be used.
        states=3
        !****************************************************************************
        ALLOCATE (state_dim_arr(states * 2), STAT = err_code)
        !Using smcl, tsoil and tstar_tile as the prognostic variables that we're interested in.
        !Find out how big each is... FORTRAN arrays start at 1 - eugh!
        smcl_size_0 = SIZE(smcl, 1) !first dimension
        smcl_size_1 = SIZE(smcl, 2) !second dimension
        smcl_size = smcl_size_0 * smcl_size_1 !no. of elements
        tsoil_size_0 = SIZE(t_soil, 1)
        tsoil_size_1 = SIZE(t_soil, 2)
        tsoil_size = tsoil_size_0 * tsoil_size_1
        tstar_tile_size_0 = SIZE(tstar_tile, 1)
        tstar_tile_size_1 = SIZE(tstar_tile, 2)
        tstar_tile_size = tstar_tile_size_0 * tstar_tile_size_1        

        !Store the array dimensions to send to the particle filter
        state_dim_arr = (/smcl_size_0, smcl_size_1, tsoil_size_0, tsoil_size_1, &
                                  tstar_tile_size_0, tstar_tile_size_1/)
        state_dims = SIZE(state_dim_arr)
        
        PRINT *, '****************states ', states, 'state_dims ', state_dims, ' array ', state_dim_arr

        !If this was a multi-node model, we'd only send this from the group rank zero (head node)
        !Ensure the particle filter knows how big the buffer is, and how many time-steps there are
        CALL MPI_SEND(state_dim_arr, state_dims, MPI_INTEGER, &
                               filter_rank_id, c_mpi_tag, empi_coupling_comms, empi_err)
        !PRINT *, 'sent state vector dimensions '


        !Allocate the buffer needed to transfer the state variables to the particle filter
        ALLOCATE (state_buffer((smcl_size_0 * smcl_size_1) + (tsoil_size_0 * tsoil_size_1) + &
                        (tstar_tile_size_0 * tstar_tile_size_1)), STAT = err_code)
        !PRINT *, 'allocated state buffer '
        IF (err_code .ne. 0) THEN
            CALL log_error("mpi_comms", &
                                  "Unable to allocate state buffer")
        END IF

    END SUBROUTINE empi_send_state_info
    
    !-----------------------------------------------------------------------------
    ! Send the number of observations to the particle filter code
    ! Currently set at compile time with a 'magic' number, this should be read from
    ! file at runtime.
    !-----------------------------------------------------------------------------
    SUBROUTINE empi_send_obs_info()
        IMPLICIT NONE
        INTEGER :: obs_size, err_code
        INTEGER, PARAMETER :: c_mpi_tag = 1 !constant
        
        !****************************************************************************
        !This is currently a 'magic' number but should ideally be read from file.
        ! It's the 'size of the observation space'
        obs_size = 100
        !****************************************************************************  
        
        CALL MPI_SEND(obs_size, 1, MPI_INTEGER, &
                      filter_rank_id, c_mpi_tag, empi_coupling_comms, empi_err)
                      
        IF (err_code .ne. 0) THEN
            CALL log_error("mpi_comms", &
                                  "Unable to allocate send observation count.")
        END IF
        
    END SUBROUTINE empi_send_obs_info
    
    !-----------------------------------------------------------------------------
    ! To facilitate each model having independent setup parameters, change its working
    ! directory to a unique path.
    ! Note that *.nml files are picked up by file name assuming that the model is being run
    ! in the same location.
    ! TODO: collect nml files and the pf_parameters.dat file and copy in as needed - currently
    ! the working directories have to be premade and populated which isn't ideal.
    !-----------------------------------------------------------------------------  
    SUBROUTINE empi_set_working_dir()
        CHARACTER (len=255) :: path, rank_dir, new_dir
        CHARACTER (len=1) :: separator = '/'
        
        CALL GETCWD(path)
        !WRITE (*, *) TRIM(path)
        !WRITE (rank_dir, '(I2.2)') empi_coupling_id
        !PRINT *, rank_dir
        new_dir = TRIM(path)//separator//rank_dir
        !PRINT *, new_dir
        CALL CHDIR(new_dir)
        CALL GETCWD(path)
        PRINT *, 'Model working directory: ', TRIM(path)
        !WRITE (*, *) TRIM(path)

    END SUBROUTINE empi_set_working_dir
        
    !-----------------------------------------------------------------------------
    ! Call out to the particle filter(s) and get result
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_do_filter()
        
        !TODO check ordering - may need to reverse - RECV then SEND and split out to be either side of model timestep
        IMPLICIT NONE
        INTEGER, PARAMETER :: c_mpi_tag = 1 !constants
        INTEGER buffer_size

        DO count0 = 1, smcl_size_0
            DO count1 = 1, smcl_size_1
                state_buffer(count0 * count1) = smcl(count0, count1)
            END DO
        END DO
        DO count0 = 1, tsoil_size_0
            DO count1 = 1, tsoil_size_1
                state_buffer((smcl_size_0*smcl_size_1) + (count0 * count1)) = t_soil(count0, count1)
            END DO
        END DO
        DO count0 = 1, tsoil_size_0
            DO count1 = 1, tsoil_size_1
                state_buffer((smcl_size_0*smcl_size_1) + (tstar_tile_size_0 * tstar_tile_size_1) + (count0 * count1))&
                            = tstar_tile(count0, count1)
            END DO
        END DO
        
        PRINT *, 'buffer sizes:'
        PRINT *, 'smcl ', smcl_size_0, ' by ', smcl_size_1, ' total=', (smcl_size_0*smcl_size_1)
        PRINT *, 'tstar ', tstar_tile_size_0, ' by ', tstar_tile_size_1, ' total=', (tstar_tile_size_0*tstar_tile_size_1)
        PRINT *, 'tsoil ', tsoil_size_0, ' by ', tsoil_size_1, ' total=', (tsoil_size_0*tsoil_size_1)
        
        buffer_size = (smcl_size_0*smcl_size_1) + (tstar_tile_size_0 * tstar_tile_size_1) + (tsoil_size_0 * tsoil_size_1)
        PRINT *, 'grand total buffer size=', buffer_size
        
        
        PRINT *, '******************************* packed state buffer for send '
        CALL MPI_SEND (state_buffer, buffer_size, MPI_DOUBLE_PRECISION, &
                                filter_rank_id, c_mpi_tag, empi_coupling_comms, empi_err)
        PRINT *, '******************sent '
        CALL MPI_RECV (state_buffer, buffer_size, MPI_DOUBLE_PRECISION, &
                                filter_rank_id, c_mpi_tag, empi_coupling_comms, MPI_STATUS_IGNORE, empi_err)
        PRINT *, '******************received '

!        CALL MPI_SENDRECV_REPLACE (state_buffer, (smcl_size_0 * smcl_size_1), MPI_DOUBLE_PRECISION, &
!                                                      filter_rank_id, c_mpi_tag, empi_coupling_id, c_mpi_tag, &
!                                                      empi_coupling_comms, MPI_STATUS_IGNORE, empi_err)
        DO count0 = 1, smcl_size_0
            DO count1 = 1, smcl_size_1
                smcl(count0, count1) = state_buffer(count0 * count1)
            END DO
        END DO
        DO count0 = 1, tsoil_size_0
            DO count1 = 1, tsoil_size_1
                t_soil(count0, count1) = state_buffer((smcl_size_0*smcl_size_1) + (count0 * count1))
            END DO
        END DO
        DO count0 = 1, tsoil_size_0
            DO count1 = 1, tsoil_size_1
                tstar_tile(count0, count1) = &
                state_buffer((smcl_size_0*smcl_size_1) + (tstar_tile_size_0 * tstar_tile_size_1) + (count0 * count1))
            END DO
        END DO        
        PRINT *, '******************************* unpacked state buffer from receive '
        
    END SUBROUTINE empi_do_filter
    
    

    
    !-----------------------------------------------------------------------------
    ! Tidy up the mpi communicators and deallocate buffer memory
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_cleanup()

        IF (world_id .eq. 0) THEN
            WRITE (empi_msg, "(A15, I10)") "Test total: ", test_total
        ENDIF

        DEALLOCATE(state_buffer)

        CALL log_info("Main", empi_msg)
        !Tidy up the MPI communicator session in each process.
        CALL MPI_FINALIZE(empi_err)
    END SUBROUTINE empi_cleanup

END MODULE empi_comms_mod