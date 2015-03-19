!     
! File:   empire_mod.F90
! Author: xw904346
!
! Created on November 6, 2014, 7:38 PM
!

MODULE empire_mod
    USE pf_control
    IMPLICIT NONE
    INCLUDE 'mpif.h'
    
    !Variable declarations for the MPI specific functions
    INTEGER :: empi_err
    INTEGER :: empi_status(MPI_STATUS_SIZE)
    INTEGER :: pf_mpi_comm                                      !the communicator between DA processes
    INTEGER :: coupling_comms                                  !the comminicator between DA and models
    INTEGER :: filter_rank                                            !the rank of this process on pf_mpi_comm
    INTEGER :: filter_coupling_id                                  !the rank of this process on coupling_comms
    INTEGER :: total_models                                        !the total number of model processes (head nodes)
    INTEGER :: num_filters                                           !the total number of DA processes
    INTEGER :: total_procs_coupled                              !the number of processes on the model/filter coupler
                                                                                !i.e. all filters and all model head nodes
    INTEGER, PARAMETER :: c_mpi_tag = 1                   !constant
    INTEGER, DIMENSION(:), ALLOCATABLE :: requests    !an array of requests sent indexed by model number
    !INTEGER, DIMENSION(:, :), ALLOCATABLE :: statuses  !and one for the mpi return status for each
    INTEGER, ALLOCATABLE, DIMENSION(:) :: gblcount     !the number of ensemble members associated with each DA process
    INTEGER, ALLOCATABLE, DIMENSION(:) :: gbldisp      !the displacements of each ensemble member relative to filter_id=0.
                                                                                 !used for mpi_gatherv and mpi_scatterv on the communicator between DA processes
      
    !Variable declarations to allow compatibility with DA code, even though the names are unhelpful
    !Set them equal to the sensible ones in here once they've all been initialised
    INTEGER :: pfrank      ! the rank of this process on PF_MPI_COMM
    INTEGER :: npfs        ! the total number of DA processes
    INTEGER :: cpl_mpi_comm ! the communicator between the empire processes and the model master nodes
    
    INTEGER :: filter_models !number of models assigned to this filter process
    INTEGER :: err_code

    !This constant should reflect the number of 2D state vectors which are being shared with the model
    INTEGER, PARAMETER :: states_2D = 3
    
    REAL(kind = KIND(1.0D0)), ALLOCATABLE, DIMENSION(:) :: state_vector
    INTEGER :: num_state_params !set in model_hook_up and used in the data transfer

    INTEGER :: filter_id              !only used in model_rank_id assignment but useful for printouts later
    INTEGER, ALLOCATABLE, DIMENSION(:) :: models               !array for holding model_rank_id refs
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: end_model_runs  !array for accumulating the finish flag of each model process
    LOGICAL, DIMENSION(:), ALLOCATABLE :: received             !an array of confirmations of request success




CONTAINS
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    SUBROUTINE empi_pf_initialise()
        IMPLICIT NONE
        INTEGER :: world_id !the rank of this process on MPI_COMM_WORLD
        INTEGER :: world_size

        INTEGER, PARAMETER :: c_coupling_colour = 9999 , c_pf_colour = 10000 

        CALL MPI_INIT(empi_err)

        PRINT *, 'PF: MPI initialising.'

        !Get the global id on the world communicator (which includes all the models too)
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, world_id, empi_err)
        !Get the total number of processes running in this MPI environment
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, world_size, empi_err)

        !get a new communicator which joins up all processes with 'c_pf_colour' and key 'world_id'
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, c_pf_colour, world_id, pf_mpi_comm, empi_err)
        !get the rank of this process on the communicator
        CALL MPI_COMM_RANK(pf_mpi_comm, filter_rank, empi_err)
        !I'm not sure this is right... since the communicator is for unique colour/key combo, surely it'll only have the one filter on it?
        CALL MPI_COMM_SIZE(pf_mpi_comm, num_filters, empi_err) !and perhaps this isn't necessary since we'll be doing ALLREDUCE in assign_models???
        
        ! create/get a communicator for the model-filter coupler
        !correction to PB's code - using world_id as key not world_size
        CALL MPI_COMM_SPLIT(MPI_COMM_WORLD, c_coupling_colour, world_id, coupling_comms, empi_err)
        PRINT*, 'split2 done'
        ! what is our id on the coupling communicator?
        CALL MPI_COMM_RANK(coupling_comms, filter_coupling_id, empi_err)
        PRINT*, 'rank done'
        ! and how many processes are on it?
        CALL MPI_COMM_SIZE(coupling_comms, total_procs_coupled, empi_err)
        PRINT*, 'size done'
        
        !_____Set up the variables for compatibility_____
        pfrank = filter_rank
        npfs = num_filters
        cpl_mpi_comm = coupling_comms
        !_____end_____
        
        total_models = total_procs_coupled - num_filters

        PRINT *, "PF empi_pf_initialise: Total coupled processes ", total_procs_coupled, &
        ", coupling rank id ", filter_coupling_id, &
        ", coupling comm id ", coupling_comms

    END SUBROUTINE empi_pf_initialise


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_assign_models()
         !defines the datatype to store various pertinent information
        USE pf_control
        IMPLICIT NONE
        INTEGER :: rtmp, iter
         !identifier of model instance
        INTEGER :: model_rank_id

        PRINT *, 'PF: assigning models '

        !let's find the models:
        !for simplicity, we can merely allocate an array that will be big enough for the largest
        !round up from number_of_models / number_of_filters
        rtmp = CEILING(real(total_models)/real(num_filters))
        PRINT *, 'array allocation ', rtmp
        !allocate exact or one over number needed but initialize to -1 so that we don't use the extra one
         !ie the number of models attached to this pf process; was 'models(rtmp)'
        ALLOCATE(pf%ranks(rtmp))
        ALLOCATE(pf%particles(rtmp))
        ALLOCATE(end_model_runs(rtmp))
        pf%ranks = -1
        pf%particles = -1  !initialise array to rubbish value
        end_model_runs = .FALSE.
        PRINT*, 'models & ranks arrays ', pf%particles, pf%ranks
        PRINT*, 'end_model_runs array ', end_model_runs

        rtmp = 1
        ! use model_rank_id starting at zero, then don't need to use -1 when testing it
        DO model_rank_id= 0, total_models-1
            IF (real(model_rank_id) .ge. real(total_models * (filter_rank))/real(num_filters) .and.&
                & real(model_rank_id) .lt. real(total_models * (filter_rank + 1))/real(num_filters)) THEN
                pf%ranks(rtmp) = model_rank_id
                !and allow DA code to use horrible indexing!
                !indexing an array with the contents of another array is not a great idea!
                pf%particles(rtmp) = model_rank_id +1
                rtmp = rtmp + 1
            END IF
        END DO
        !the loop above distributes models evenly amongst filters hence 'filter_models' is the number on this filter and is less than total_models
        filter_models = rtmp-1
        pf%count = filter_models
        pf%nens = total_models
        
        !initialise array to hold the initial receive test results
        ALLOCATE (received(filter_models))
        received = .FALSE.
        !and allocate the space for the other arrays now we know numbers
        ALLOCATE (requests(filter_models))
        !ALLOCATE (statuses(mpi_status_size, filter_models))
        
         !doing the allgather...
        ALLOCATE (gblcount(num_filters))
        ALLOCATE (gbldisp(num_filters))
        !work out what the hell that ALLGATHER is doing in comms:initialise_mpi
        !there are no comments to explain :(
        !aha, it's sending the number of models each filter process has, then collecting the
        !array with each's number. This should be identical for all of course!
        CALL MPI_ALLGATHER(filter_models, 1, MPI_INTEGER, gblcount, 1, MPI_INTEGER, &
                            pf_mpi_comm, empi_err)
                 
        !PRINT *, 'ALLGATHER success.'
        !PRINT *, 'gblcount array: ', gblcount
        gbldisp = 0

        IF (filter_models .GT. 1) THEN
            DO iter = 2, num_filters
                !PRINT *, 'iterator = ', iter
                gbldisp(iter) = gbldisp(iter - 1) + gblcount(iter - 1)
            END DO
        END IF
        PRINT *, 'gbldisp array: ', gbldisp
        PRINT *, '_____________________________________________________________________'
        PRINT *, 'Total coupled processes = ', total_procs_coupled
        PRINT *, 'Number of filters = ', num_filters
        PRINT *, '-> Number of models = ', total_models, ' and num models I own = ', filter_models
        PRINT *, 'Coupling rank id = ', filter_coupling_id
       ! PRINT *, "-> Particle filter local relative rank = ", m_filter_id_rank, " and I own"
        DO rtmp=1, filter_models
            PRINT*, " model rank: " , pf%ranks(rtmp), " and particle id: ", pf%particles(rtmp) 
            PRINT*, 'end_model_run: ', end_model_runs(rtmp)
        END DO


        
    END SUBROUTINE empi_pf_assign_models
    
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_model_hookup_state_dims()
        !get state vector dimensions, 
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !TODO timeout/timestep switch, and whichever value is then needed
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        INTEGER :: count, state1_dim1, state1_dim2, &
                        state2_dim1, state2_dim2, &
                        state3_dim1, state3_dim2
        INTEGER, ALLOCATABLE, DIMENSION(:) :: state_dims
        !identifier of model instance
        INTEGER :: model_rank_id


        !let's get the information on the state buffer which we'll need later
        !this assumes each state will have serialised 2D data, therefore there are
        !states_2D*2 numbers in the array, state_dims.
        ALLOCATE (state_dims(states_2D * 2), stat = err_code)
        state_dims = 0
        PRINT *, "state size buffer allocated ", state_dims

        !Do a loop like with the other send/receives so that we get stuff from the right number of models
        !but in fact we only really need one model to tell us.
        DO count = 1, filter_models
            model_rank_id = pf%ranks(count) 
            CALL MPI_RECV(state_dims, states_2D *2, MPI_INTEGER, &
                                   model_rank_id, c_mpi_tag, coupling_comms, empi_status, empi_err)
             PRINT *, 'PF: received state dimensions individual.', state_dims
        END DO

        PRINT *, 'PF: received state dimensions.', state_dims

        IF (states_2D .eq. 1) THEN
            !we've only one 2D state
            state1_dim1 = state_dims(1)
            state1_dim2 = state_dims(2)
        ENDIF

        IF (states_2D .eq. 3) THEN
            !we've three 2D states
            state1_dim1 = state_dims(1)
            state1_dim2 = state_dims(2)
            state2_dim1 = state_dims(3)
            state2_dim2 = state_dims(4)
            state3_dim1 = state_dims(5)
            state3_dim2 = state_dims(6)
        ENDIF
        
        num_state_params = ((state1_dim1 * state1_dim2) + &
                                        (state2_dim1 * state2_dim2) +&
                                        (state3_dim1 * state3_dim2))
        ALLOCATE(state_vector(num_state_params), stat = err_code)
        model_params%state_dim = num_state_params
        PRINT *, 'model_params%state_dim: ', model_params%state_dim
        IF (err_code .ne. 0) WRITE(*, *) "Particle filter could not allocate state_vector."
        PRINT *, 'State vector created '

        DEALLOCATE(state_dims)
        
        

    END SUBROUTINE empi_pf_model_hookup_state_dims
    
    !-----------------------------------------------------------------------------
    LOGICAL FUNCTION empi_pf_model_hookup_obs_size()
        INTEGER :: obs_size, count, model_rank_id
        
        empi_pf_model_hookup_obs_size = .FALSE.
        
        DO count = 1, filter_models
            model_rank_id = pf%ranks(count) 
            CALL MPI_RECV(obs_size, 1, MPI_INTEGER, &
                                   model_rank_id, c_mpi_tag, coupling_comms, empi_status, empi_err)
        END DO
              
        
        model_params%obs_dim = obs_size
        empi_pf_model_hookup_obs_size = .TRUE.
        RETURN
        
    END FUNCTION empi_pf_model_hookup_obs_size
        
    !-----------------------------------------------------------------------------

    !-----------------------------------------------------------------------------
    SUBROUTINE empi_pf_do_work(b_obs_known)
        USE Qdata
        USE Rdata
        IMPLICIT NONE
        INTEGER :: iter
        LOGICAL :: b_obs_known

        !random seed thingy
        !this is in gen_rand.f90
        CALL random_seed_mpi(filter_id)
        !was (coupling_id)
        !setup Q and R
        IF (.NOT. pf%gen_Q) THEN
             !this just loads pf_control and sizes - seems a bit pointless
            CALL loadQ
             !loads sizes and does a bit of stuff
            CALL loadR
        END IF

        IF (.NOT. pf%gen_Q) THEN
            !receive the model's initial state at time, t, after spin-up
            PRINT*, 'empi_pf_model_first_data____________'
            CALL empi_pf_model_first_data()

            !initial perturbation
            !we have data in pf % psi
            DO iter = 1, pf%count
                CALL perturb_particle(pf%psi(:, iter))
            END DO

            !********************************************************************
            ! This check gives us the opportunity to do something else if we don't
            ! know the size of the  observation space (which shouldn't happen really), or
            ! to pass the switch on and decide within empi_pf_do_known_timestep_filter()
            !********************************************************************
            IF (b_obs_known .EQV. .TRUE.) THEN
                CALL empi_pf_do_known_timestep_filter()
            END IF
            
            PRINT*, 'PF: finished the loop - now to tidy up'
            !send the last state. This would otherwise block the looping receives in the model
            PRINT*, 'empi_pf_model_last_data____________'
            CALL empi_pf_model_last_data()

        ELSE
            CALL genQ

        END IF

    END SUBROUTINE empi_pf_do_work
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    SUBROUTINE empi_pf_model_first_data()
        IMPLICIT NONE
         !identifier of model instance
        INTEGER :: model_rank_id

        INTEGER :: iter
        LOGICAL :: b_success

        ! Launch a non_blocking receive for data from all models attached to this process
        DO iter = 1, pf%count
            model_rank_id = pf%ranks(iter)
            CALL MPI_IRECV(pf%psi(:, iter), model_params%state_dim, MPI_DOUBLE_PRECISION, &
                                    model_rank_id, c_mpi_tag, coupling_comms, requests(iter), empi_err)
            PRINT *, 'open receive sent to model number: ', model_rank_id

        END DO
        
        PRINT *, 'request handles array initally: ', requests

        !now keep testing for results
        iter = 0
        DO
            !this is a sneaky way of repeating valid array indices between 1 and pf%count
            iter = MOD(iter, pf%count) + 1
            IF (.NOT.received(iter)) THEN
                model_rank_id = pf%ranks(iter)
                PRINT *, 'testing position: ', iter, 'with request handle: ', requests(iter)
                CALL MPI_TEST(requests(iter), b_success, MPI_STATUS_IGNORE, empi_err)
                PRINT *, 'test done'
                IF (b_success) THEN
                    received(iter) = .TRUE.
                    PRINT *, 'received array: ', received
                END IF
            ENDIF
            IF (ALL(received)) EXIT
        END DO


    END SUBROUTINE empi_pf_model_first_data
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------

    SUBROUTINE empi_pf_do_known_timestep_filter()
        !Do this if we know the time step information from the model hookup
        !execute the data assimilation filters and other magic stuff
        IMPLICIT NONE

        INTEGER :: iter1, iter2
!This is the 'design':        
!    loop number obs
!        loop steps between obs
!           T+1
!           filters
!        end
!        T+1
!        filters
!    end
        CALL output_from_pf
        IF (pf%gen_data) CALL save_truth(pf%psi(:, 1))
        IF (pf%use_traj) CALL trajectories
        !start_t = mpi_wtime()

        DO iter1 = 1, pf%time_obs !number of observations to assimilate
            PRINT*, 'PF: observation counter = ', iter1
            DO iter2 = 1, pf%time_bwn_obs - 1
                pf%timestep = pf%timestep + 1

                SELECT CASE (pf%filter)
                CASE ("EW")
                    CALL proposal_filter
                CASE ("SI", "SE")
                    CALL stochastic_model 
                CASE ("ET")
                    CALL deterministic_model
                CASE DEFAULT
                    PRINT*, 'Error -555: Incorrect pf%filter'
                END SELECT

                IF (pf%use_traj) CALL trajectories
                CALL output_from_pf
            END DO !timesteps between observations

            pf%timestep = pf%timestep + 1
            PRINT*, 'starting the equal weight filter step'

                SELECT CASE (pf%filter)
                CASE ("EW")
                    CALL equal_weight_filter
                CASE ("SI")
                    CALL sir_filter 
                CASE ("SE")
                    CALL stochastic_model
                    CALL diagnostics
                CASE ("ET")
                    CALL deterministic_model
                    CALL letkf_analysis
                CASE DEFAULT
                    PRINT*, 'Error -556: Incorrect pf%filter'
                END SELECT

            PRINT*, 'PF: timestep = ', pf%timestep, 'after equal weight filter'

            IF (pf%gen_data) CALL save_truth(pf%psi(:, 1))
            IF (pf%use_traj) CALL trajectories
            CALL output_from_pf

            !CALL reconfigure_model

        END DO
        CALL diagnostics
        
        
    END SUBROUTINE empi_pf_do_known_timestep_filter
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_pf_model_last_data()
        !do last send for all models
        IMPLICIT NONE
        INTEGER :: iter
        INTEGER :: model_rank_id !identifier of model instance

        DO iter=1, pf%count
            model_rank_id = pf%ranks(iter)
            CALL MPI_ISEND(pf%psi(:,iter), model_params%state_dim, MPI_DOUBLE_PRECISION, model_rank_id, &
                    c_mpi_tag, coupling_comms, requests(model_rank_id), empi_err)
        END DO
        
        CALL MPI_WAITALL(pf%count, requests, MPI_STATUSES_IGNORE, empi_err)    
        
    END SUBROUTINE empi_pf_model_last_data
    
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_pf_cleanup()
        !cleanup
        CALL MPI_FINALIZE(empi_err)
        DEALLOCATE(state_vector)
        DEALLOCATE(pf%ranks)
        DEALLOCATE(pf%particles)
        DEALLOCATE(end_model_runs)
        DEALLOCATE(received)
        DEALLOCATE (gblcount)
        DEALLOCATE (gbldisp)

    END SUBROUTINE empi_pf_cleanup
    
    !-----------------------------------------------------------------------------
    
    
!    
!    This subroutine was written for the original specification of the DA code and may still have some
!    useful stuff.
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_pf_do_filter()
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !need to do the random seed thingy... gen_rand:random_seed_mpi(pfid)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        IMPLICIT NONE

        INTEGER :: iter1, count, now
        INTEGER :: model_rank_id
        INTEGER, PARAMETER :: c_tag = 1, timeout = 3
        LOGICAL :: all_complete, msg_flag
         
!        !Set up the flag for entry into the loop
!        !Need a logical 'AND' but must set initial value of accumulator to first result
!        all_complete = end_model_runs(1)
!        PRINT *, 'START______________ all_complete: ', all_complete
!        IF (filter_models .GE. 2) THEN
!            DO iter1 = 2, filter_models
!                all_complete = all_complete .AND. end_model_runs(iter1)
!                PRINT *, 'checking---- all_complete: ', all_complete, 'AND end_model_runs(iter1) ', &
!                        iter1, end_model_runs(iter1)
!            END DO
!        END IF
            
        !do line 118 etc in pf_couple.f90
        count = 1
        
        !continue this loop as long as any of the models is still running        
        DO WHILE (.NOT. ALL(end_model_runs))
            !(all_complete .eqv. .FALSE.)
            PRINT *, '___________________LOOP__ ', count
            
            DO iter1 = 1, filter_models
                
                !Test each model in turn and ensure we do the comms until it has finished
                IF (end_model_runs(iter1) .EQV. .FALSE.) THEN                        
                    model_rank_id= pf%ranks(iter1)
                    !models(iter1)
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


            DO iter1 = 1, filter_models
                !Test each model in turn and ensure we do the comms until it has finished
                IF (end_model_runs(iter1) .EQV. .FALSE.) THEN
                    model_rank_id= pf%ranks(iter1) !models(iter1)                  
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

            IF (filter_models .GE. 2) THEN
                DO iter1 = 2, filter_models
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



END MODULE empire_mod