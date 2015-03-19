!     
! File:   pf_mpi_mod.F90
! Author: xw904346
!
! Created on November 6, 2014, 7:38 PM
!

MODULE pf_mpi_mod
    IMPLICIT NONE

    INCLUDE 'mpif.h'

    !Variable declarations for the module
    INTEGER :: empi_err
    INTEGER :: empi_status(MPI_STATUS_SIZE)
    
    INTEGER :: coupling_comms, pf_mpi_comm
    INTEGER :: total_procs_coupled, filter_coupling_id !found from size() and rank() calls in the setup
    INTEGER :: num_models, num_filters

    !This constant should reflect the number of state vectors which are being shared with the model
    INTEGER, PARAMETER :: states_2D = 3
    
    REAL(kind = KIND(1.0D0)), ALLOCATABLE, DIMENSION(:) :: state_vector
    INTEGER :: num_state_params !set in model_hook_up and used in the data transfer

    INTEGER :: m_filter_id_rank
    INTEGER :: filter_id        !only used in model_rank_id assignment but useful for printouts later
    INTEGER :: model_rank_id    !identifier of model instance
    INTEGER, ALLOCATABLE, DIMENSION(:) :: models !array for holding model_rank_id refs
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: end_model_runs !array for accumulating the finish flag of each model process




CONTAINS


    SUBROUTINE empi_pf_initialise()
        IMPLICIT NONE
        INTEGER, PARAMETER :: c_coupling_colour = 9999, c_pf_colour = 1
        INTEGER :: world_id, world_size

        CALL MPI_INIT(empi_err)

        PRINT *, 'PF: MPI initialising.'

        !Get the global id on the world communicator (which includes all the models too)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, world_id, empi_err)
        !Get the total number of processes running in this MPI environment
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, empi_err)
        
        ! this split() is only used to find 'm_filter_id_rank'
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, c_pf_colour, world_id, pf_mpi_comm, empi_err)
!        PRINT*, MPI_COMM_WORLD, c_pf_colour, world_id, pf_mpi_comm, empi_err
        
        !get the rank of this process on the communicator
        CALL MPI_COMM_RANK(pf_mpi_comm, m_filter_id_rank, empi_err) !m_filter_id_rank==pfrank
        PRINT*, 'rank done'
        CALL MPI_COMM_SIZE(pf_mpi_comm, num_filters, empi_err) !num_filters==npfs
        PRINT*, 'size done'
        
        ! create/get a communicator for the model-filter coupler
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, c_coupling_colour, world_id, coupling_comms, empi_err)
        PRINT*, 'split done'
        ! what is our id on the coupling communicator?
        CALL MPI_COMM_RANK(coupling_comms, filter_coupling_id, empi_err)
        PRINT*, 'rank done'
        ! and how many processes are on it?
        CALL MPI_COMM_SIZE(coupling_comms, total_procs_coupled, empi_err)
        PRINT*, 'size done'
       


        PRINT *, "PF empi_pf_initialise: Total coupled processes ", total_procs_coupled, &
                    ", coupling rank id ", filter_coupling_id, &
                    ", coupling comm id ", coupling_comms

    END SUBROUTINE empi_pf_initialise


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_assign_models()
        IMPLICIT NONE
        INTEGER :: total_models, rtmp, count

        PRINT *, 'PF: assigning models '

        total_models = total_procs_coupled - num_filters

        !let's find the models:
        !for simplicity, we can merely allocate an array that will be big enough for the largest:
        !round up from number_of_models / number_of_filters
        rtmp = CEILING(real(total_models)/real(num_filters))
        PRINT *, 'array allocation ', rtmp
        !allocate exact or one over number needed but initialize to -1 so that we don't use the extra one
        ALLOCATE(models(rtmp)) !ie the number of models attached to this pf process
        ALLOCATE(end_model_runs(rtmp))
        models = -1
        end_model_runs = .FALSE.
        PRINT*, 'models array ', models
        PRINT*, 'end_model_runs array ', end_model_runs

        rtmp = 1
        DO model_rank_id= 0, total_models-1
            IF (real(model_rank_id) .ge. real(total_models * (m_filter_id_rank))/real(num_filters) .and.&
                & real(model_rank_id) .lt. real(total_models * (m_filter_id_rank + 1))/real(num_filters)) THEN
                models(rtmp) = model_rank_id!put this before increment and then we don't need the -1 everywhere
                rtmp = rtmp + 1
            END IF
        END DO
        num_models = rtmp-1

        PRINT *, 'Total coupled processes = ', total_procs_coupled
        PRINT *, 'Number of filters = ', num_filters
        PRINT *, '-> Number of models = ', total_models
        PRINT *, 'Coupling rank id = ', filter_coupling_id
        PRINT *, "-> Particle filter local relative rank = ", m_filter_id_rank, " and I own"
        DO count=1, num_models
            PRINT*, " model particle: " , models(count)
            PRINT*, 'end_model_run: ', end_model_runs(count)
        END DO
        
    END SUBROUTINE empi_pf_assign_models
    
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_model_hookup()

        IMPLICIT NONE
        INTEGER :: count, state1_dim1, state1_dim2, &
                        state2_dim1, state2_dim2, &
                        state3_dim1, state3_dim2
        INTEGER, ALLOCATABLE, DIMENSION(:) :: state_dims
        INTEGER, PARAMETER :: c_mpi_tag = 1 !constant


        !let's get the information on the state buffer which we'll need later
        !this assumes each state will have serialised 2D data, therefore there are
        !states_2D*2 numbers in the array, state_dims.
        ALLOCATE (state_dims(states_2D * 2), stat = empi_err)
        PRINT *, "state size buffer allocated "

        !Do a loop like with the other send/receives so that we get stuff from the right number of models
        DO count = 1, num_models
            model_rank_id= models(count)
            CALL MPI_RECV(state_dims, states_2D *2, MPI_INTEGER, &
                                   model_rank_id, c_mpi_tag, coupling_comms, empi_status, empi_err)
        end DO

        PRINT *, 'PF: received state dimensions.'

        IF (states_2D .eq. 1) THEN !we've only one 2D state
            state1_dim1 = state_dims(1)
            state1_dim2 = state_dims(2)
        ENDIF

        IF (states_2D .eq. 3) THEN !we've three 2D states
            state1_dim1 = state_dims(1)
            state1_dim2 = state_dims(2)
            state2_dim1 = state_dims(3)
            state2_dim2 = state_dims(4)
            state3_dim1 = state_dims(5)
            state3_dim2 = state_dims(6)
        ENDIF
        
        num_state_params = state1_dim1 * state1_dim2
        ALLOCATE(state_vector(num_state_params), stat = empi_err)
        IF (empi_err .ne. 0) WRITE(*, *) "Particle filter could not allocate state_vector."
        PRINT *, 'State vector created '

        DEALLOCATE(state_dims)

    END SUBROUTINE empi_pf_model_hookup

    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------

   SUBROUTINE empi_pf_do_filter()
        IMPLICIT NONE

        INTEGER :: iter1, iter2, count, now
        INTEGER :: request_handle, index
        INTEGER, PARAMETER :: c_tag = 1, timeout = 3
        LOGICAL :: all_complete, msg_flag
         
        !Set up the flag for entry into the loop
        !Need a logical 'AND' but must set initial value of accumulator to first result
        all_complete = end_model_runs(1)
        PRINT *, 'START______________ all_complete: ', all_complete
        IF (num_models .GE. 2) THEN
            DO iter1 = 2, num_models
                all_complete = all_complete .AND. end_model_runs(iter1)
                PRINT *, 'checking---- all_complete: ', all_complete, 'AND end_model_runs(iter1) ', &
                        iter1, end_model_runs(iter1)
            END DO
        END IF
            
        count = 1
        
        !continue this loop as long as any of the models is still running        
        DO WHILE (all_complete .eqv. .FALSE.)
            PRINT *, '___________________LOOP__ ', count
            
            DO iter1 = 1, num_models
                
                !Test each model in turn and ensure we do the comms until it has finished
                IF (end_model_runs(iter1) .EQV. .FALSE.) THEN                        
                    model_rank_id= models(iter1)
                    now = TIME()
                    DO WHILE (TIME() .LT. now + timeout)
                        msg_flag = .FALSE.
                        CALL MPI_IPROBE(model_rank_id, c_tag, coupling_comms, &
                                    msg_flag, empi_status, empi_err)
                        IF (msg_flag .EQV. .TRUE.) THEN
                            PRINT *, 'Message waiting '
                            EXIT
                        END IF
                    END DO
                    PRINT *, 'timeout loop complete, flag= ', msg_flag

                    IF (msg_flag .EQV. .TRUE.) THEN
                        
                        CALL MPI_RECV(state_vector, num_state_params, MPI_DOUBLE_PRECISION, &
                                      model_rank_id, c_tag, coupling_comms, empi_status, empi_err)
                        PRINT*, 'Particle filter ', filter_id, 'has received state_vectors over mpi for model ',model_rank_id  

                    ELSE !no message waiting, flag should be False, i.e. the run *has* finished
                        end_model_runs(iter1) = .NOT. msg_flag
                        PRINT *, 'receive end_run flag ',end_model_runs(iter1)       
                    END IF
                END IF
            END DO
            

            !********************************************************
            !unpack the state vectors and do stuff with them HERE
            !********************************************************


            DO iter1 = 1, num_models
                !Test each model in turn and ensure we do the comms until it has finished
                IF (end_model_runs(iter1) .EQV. .FALSE.) THEN
                    model_rank_id= models(iter1)                  
                    !If it's still running, send the state vectors back
                    IF (end_model_runs(iter1) .eqv. .FALSE.) THEN
                        CALL MPI_SEND(state_vector, num_state_params, MPI_DOUBLE_PRECISION, &
                                               model_rank_id, c_tag, coupling_comms, empi_err)
                        PRINT*, 'Particle filter ', filter_id, 'has sent state_vectors over mpi for model ',model_rank_id
                    !ELSE !run is complete
    !                    PRINT *, 'getting completion message'
                        !PRINT *, 'receive end_run flag ',end_model_runs(iter1)       
 
                    END IF
                END IF
            END DO
            
                    
            !Check the flag for continuation of the loop
            !Need a logical 'AND' but must set initial value of accumulator to first result
            all_complete = end_model_runs(1)
            PRINT *, 'END______________ all_complete: ', all_complete

            IF (num_models .GE. 2) THEN
                DO iter1 = 2, num_models
                    all_complete = all_complete .AND. end_model_runs(iter1)
                    PRINT *, 'checking---- all_complete: ', all_complete, &
                        'AND end_model_runs(iter1) ', iter1, end_model_runs(iter1)
                END DO
            END IF
                       
            PRINT *, 'run complete flag ', all_complete
            count = count + 1

        END DO !main loop
        
        IF (all_complete .eqv. .TRUE.) THEN
                PRINT *, '************ MODEL RUNS COMPLETE ***************'
        END IF
        
    END SUBROUTINE empi_pf_do_filter
    
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_cleanup()
        !cleanup
        CALL MPI_FINALIZE(empi_err)
        DEALLOCATE(state_vector)
        DEALLOCATE(models)
        DEALLOCATE(end_model_runs)

    END SUBROUTINE empi_pf_cleanup


END MODULE pf_mpi_mod