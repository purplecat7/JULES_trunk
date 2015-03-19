!     
! File:   empi_communications.F90
! Author: xw904346 (Jane Lewis, Dept of Meteorology, UREAD)
!
! Created on August 11, 2014, 3:22 PM
!

MODULE empi_comms_mod

    USE logging_mod, ONLY: log_init, log_info, log_debug, log_warn, log_error, log_fatal
    USE prognostics, ONLY: smcl, t_soil, tstar_tile

    IMPLICIT NONE

    INCLUDE 'mpif.h'


    INTEGER     empi_err, &                   !holds error value
                      empi_numtasks, &          !total number of tasks running in the MPI environment
                      empi_rank, &                 !the id of a task instance running this code
                      empi_errcode, &            !error code used for ABORT function
                      empi_strlen, &               !processor name string size
                      empi_coupling_comms, &     !group communicator for model and particle filter processes
                      empi_total_procs, &        !total processes including model and particle filter
                      empi_model_procs, &     !total model processes
                      empi_coupling_id, &   !id of the model on the group communicator
                      filter_rank_id, &                  !process id of the particle filter this model instance will use
                      empi_request_handle        !for test purposes only
    !This is the 'magic' number for hooking into the correct communicator group for the
    !particle filter
    INTEGER, PARAMETER :: c_coupling_colour = 9999

    !INTEGER open_err !for writing to file
    CHARACTER (LEN = 50) :: empi_msg
    INTEGER :: test_total = 0

    
    !the size of each of the 2 dimensions
    INTEGER     smcl_size_0, smcl_size_1, &     !for the prognostic variable 'SMCL'
                      t_soil_0, t_soil_1, &                  !for T_SOIL
                      tstar_tile_0, tstar_tile_1, &        !for TSTAR_TILE
                      count0, count1
    REAL(KIND = KIND(1.0D0)), ALLOCATABLE, DIMENSION(:) :: state_buffer

CONTAINS

    !-----------------------------------------------------------------------------
    ! Initialise MPI and communicators.
    !
    !
    !  
    !
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_initialise()

        IMPLICIT NONE
        INTEGER :: count, npfs
        !INTEGER, PARAMETER :: c_particle_filter_flag_no = 0
        CHARACTER*(MPI_MAX_PROCESSOR_NAME) empi_name
        CHARACTER (LEN = 50) :: my_fmt

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
            !Get the size of the group (empi_numtasks) associated with the communicator (MPI_COMM_WORLD)
            CALL MPI_COMM_SIZE(MPI_COMM_WORLD, empi_numtasks, empi_err)
            WRITE (empi_msg, "(A20, I2)") "MPI process total: ", empi_numtasks
            CALL log_info("Main", empi_msg)



        END IF

        PRINT*, 'mpi init, getting proc name '
        !Get the unique name of the processor node (and string length)
        CALL MPI_GET_PROCESSOR_NAME(empi_name, empi_strlen, empi_err)
        !Get the unique id of the process (its rank) on this communicator - each process will
        !print out this information but the rank will be different
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, empi_rank, empi_err)
        PRINT*, 'got rank: ', empi_rank
        !create the output format on the fly knowing the string length - write the format to a string using a format! crazy shit :)
        WRITE (my_fmt, "(A23, I2, A1)") '(A12, I3, A4, I3, A4, A', empi_strlen, ')'
        WRITE (empi_msg, my_fmt) "MPI process ", empi_rank, " of ", empi_numtasks, " on ", empi_name
        !sadly, this notation doesn't work on the compiler...
        !WRITE (empi_msg, "(A12, I3, A4, I3, A4, <mpi_strlen>A1)") "MPI process ",mpi_rank, " of ", empi_numtasks, " on ", empi_name

        PRINT*, empi_msg
        
        !Get a communicator, empi_coupling_comms, which hooks into the particle filter
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, c_coupling_colour, empi_rank, empi_coupling_comms, empi_err)
        PRINT*, 'split done '
        !Find the number of processes on this new communicator - model AND particle filters
        CALL MPI_COMM_SIZE(empi_coupling_comms, empi_total_procs, empi_err)
        PRINT*, 'size done '
        !Find what rank we are on this new communicator
        CALL MPI_COMM_RANK(empi_coupling_comms, empi_coupling_id, empi_err)
        PRINT *, "Model: num elements ", empi_total_procs, ", coupling rank ", empi_coupling_id, &
                    ", coupling comm id ", empi_coupling_comms
        PRINT*, 'rank done '
        !Count up the total number of particle filter processes - particle filter flags will be set to 1 in them
        CALL MPI_ALLREDUCE(c_particle_filter_flag_no, npfs, 1, MPI_INTEGER, MPI_SUM, empi_coupling_comms, empi_err)
        PRINT*, 'done all reduce '
        
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
            WRITE (empi_msg, "(A15, I2)") "From MPI rank: ", empi_rank

            IF (empi_rank .eq. 0) THEN

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
    ! Initialize the prognostic variable buffer and send its size to particle filter code so that
    ! it can initialize correctly
    ! SMCL and T_SOIL are imported from the prognostics module
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_start()
        IMPLICIT NONE
        INTEGER err_code, smcl_size, tsoil_size, tstar_tile_size, states, state_dims
        INTEGER, ALLOCATABLE, DIMENSION(:) :: state_dim_arr !an array to hold the state vector dimensions
        INTEGER, PARAMETER :: c_mpi_tag = 1 !constants

        !Set up an array to store the dimensions of the state arrays so that the particle
        !filter knows them
        !This is currently a 'magic' number but should ideally be read from file along with
        !the prognostic state variables to be used.
        states=3
        ALLOCATE (state_dim_arr(states * 2), STAT = err_code)
        !Using smcl, tsoil and tstar_tile as the prognostic variables that we're interested in.
        !Find out how big each is... FORTRAN arrays start at 1 - eugh!
        smcl_size_0 = SIZE(smcl, 1) !first dimension
        smcl_size_1 = SIZE(smcl, 2) !second dimension
        smcl_size = smcl_size_0 * smcl_size_1 !no. of elements
        tsoil_size_0 = SIZE(tsoil, 1)
        tsoil_size_1 = SIZE(tsoil, 2)
        tsoil_size = tsoil_size_0 * tsoil_size_1
        tstar_tile_size_0 = SIZE(tstar_tile, 1)
        tstar_tile_size_1 = SIZE(tstar_tile, 2)
        tstar_tile_size = tstar_tile_size_0 * tstar_tile_size_1        

        !Store the array dimensions to send to the particle filter
        state_dim_arr = (/smcl_size_0, smcl_size_1, tsoil_size_0, tsoil_size_1, &
                                  tstar_tile_size_0, tstar_tile_size_1/)
        state_dims = SIZE(state_dim_arr)
        
        PRINT *, '****************states ', states, 'state_dims ', state_dims

        !If this was a multi-node model, we'd only send this from the group rank zero
        !Ensure the particle filter knows how big the buffer is, and how many time-steps there are
        CALL MPI_SEND(state_dim_arr, state_dims, MPI_INTEGER, &
                               filter_rank_id, c_mpi_tag, empi_coupling_comms, empi_err)
        PRINT *, 'sent state vector dimensions '


        !Allocate the buffer needed to transfer the state variables to the particle filter
        ALLOCATE (state_buffer((smcl_size_0 * smcl_size_1) + (tsoil_size_0 * tsoil_size_1) + &
                        (tstar_tile_size_0 * tstar_tile_size_1)), STAT = err_code)
        PRINT *, 'allocated state buffer '
        IF (err_code .ne. 0) THEN
            CALL log_error("mpi_comms", &
                                  "Unable to allocate state buffer")
        END IF

    END SUBROUTINE empi_start
    
    !-----------------------------------------------------------------------------
    ! To facilitate each model having independent setup parameters, change its working
    ! directory to a unique path.
    ! Note that *.nml files are picked up by file name assuming that the model is being run
    ! in the same location.
    !-----------------------------------------------------------------------------  
    SUBROUTINE empi_set_working_dir()
        CHARACTER (len=255) :: path, rank_dir, new_dir
        CHARACTER (len=1) :: separator = '/'
        
        CALL GETCWD(path)
        WRITE (*, *) TRIM(path)
        WRITE (rank_dir, '(I2.2)') empi_coupling_id
        PRINT *, rank_dir
        new_dir = TRIM(path)//separator//rank_dir
        PRINT *, new_dir
        CALL CHDIR(new_dir)
        CALL GETCWD(path)
        WRITE (*, *) TRIM(path)

    END SUBROUTINE empi_set_working_dir
    !-----------------------------------------------------------------------------
    ! Call out to the particle filter(s) and get result
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_do_filter()
        
        !TODO check ordering - may need to reverse - RECV then SEND and split out to be either side of model timestep
        IMPLICIT NONE
        INTEGER, PARAMETER :: c_mpi_tag = 1 !constants
        INTEGER smcl_size_0, smcl_size_1

        DO count0 = 1, smcl_size_0
            DO count1 = 1, smcl_size_1
                state_buffer(count0 * count1) = smcl(count0, count1)
            END DO
        END DO
        DO count0 = 1, tsoil_size_0
            DO count1 = 1, tsoil_size_1
                state_buffer((smcl_size_0*smcl_size_1) + (count0 * count1)) = tsoil(count0, count1)
            END DO
        END DO
        DO count0 = 1, tsoil_size_0
            DO count1 = 1, tsoil_size_1
                state_buffer((smcl_size_0*smcl_size_1) + (tstar_tile_size_0 * tstar_tile_size_1) + (count0 * count1))&
                            = tstar_tile(count0, count1)
            END DO
        END DO
        PRINT *, '******************************* packed state buffer for send '
        CALL MPI_SEND (state_buffer, (smcl_size_0 * smcl_size_1), MPI_DOUBLE_PRECISION, &
                                filter_rank_id, c_mpi_tag, empi_coupling_comms, empi_err)
        PRINT *, '******************sent '
        CALL MPI_RECV (state_buffer, (smcl_size_0 * smcl_size_1), MPI_DOUBLE_PRECISION, &
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
                tsoil(count0, count1) = state_buffer((smcl_size_0*smcl_size_1) + (count0 * count1))
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
    ! Tidy up the mpi communicators
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_cleanup()

        IF (empi_rank .eq. 0) THEN
            WRITE (empi_msg, "(A15, I10)") "Test total: ", test_total
        ENDIF

        DEALLOCATE(state_buffer)

        CALL log_info("Main", empi_msg)
        !Tidy up the MPI communicator session in each process.
        CALL MPI_FINALIZE(empi_err)
    END SUBROUTINE empi_cleanup

END MODULE empi_comms_mod