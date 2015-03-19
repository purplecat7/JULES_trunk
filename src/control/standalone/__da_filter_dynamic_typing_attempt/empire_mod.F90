!     
! File:   empire_mod.F90
! Author: xw904346
!
! Created on November 6, 2014, 7:38 PM
!

MODULE empire_mod
    USE empi_wrappers
    
    IMPLICIT NONE

    !INCLUDE 'mpif.h'

    !Variable declarations for the module
    !INTEGER :: empi_err
    !INTEGER :: empi_status(MPI_STATUS_SIZE)
    
    !INTEGER :: coupling_comms
    !INTEGER :: total_procs_coupled, filter_coupling_id !found from size() and rank() calls in the initialisation
    INTEGER :: filter_models !number of models assigned to this filter process
    INTEGER :: err_code

    !This constant should reflect the number of 2D state vectors which are being shared with the model
    INTEGER, PARAMETER :: states_2D = 3
    
    REAL(kind = KIND(1.0D0)), ALLOCATABLE, DIMENSION(:) :: state_vector
    INTEGER :: num_state_params !set in model_hook_up and used in the data transfer

    INTEGER :: filter_id              !only used in model_rank_id assignment but useful for printouts later
    INTEGER, ALLOCATABLE, DIMENSION(:) :: models !array for holding model_rank_id refs
    LOGICAL, ALLOCATABLE, DIMENSION(:) :: end_model_runs !array for accumulating the finish flag of each model process
    LOGICAL, DIMENSION(:), ALLOCATABLE :: received  !an array of confirmations of request success




CONTAINS


    SUBROUTINE empi_pf_initialise()
        call empi_wrap_initialise()

    END SUBROUTINE empi_pf_initialise


    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_assign_models()
        USE pf_control !defines the datatype to store various pertinent information
        IMPLICIT NONE
        INTEGER :: rtmp !,num_particle_filters, total_models
        INTEGER :: model_rank_id    !identifier of model instance

        PRINT *, 'PF: assigning models '

        !let's find the models:
        !for simplicity, we can merely allocate an array that will be big enough for the largest
        !round up from number_of_models / number_of_filters
        rtmp = CEILING(real(total_models)/real(num_filters))
        PRINT *, 'array allocation ', rtmp
        !allocate exact or one over number needed but initialize to -1 so that we don't use the extra one
        ALLOCATE(pf%particles(rtmp)) !ie the number of models attached to this pf process; was 'models(rtmp)'
        ALLOCATE(end_model_runs(rtmp))
        pf%particles = -1 !models = -1 !initialise array to rubbish value
        end_model_runs = .FALSE.
        PRINT*, 'models array ', pf%particles !models
        PRINT*, 'end_model_runs array ', end_model_runs

        rtmp = 1
        DO model_rank_id= 0, total_models-1
            IF (real(model_rank_id) .ge. real(total_models * (filter_rank))/real(num_filters) .and.&
                & real(model_rank_id) .lt. real(total_models * (filter_rank + 1))/real(num_filters)) THEN
                pf%particles(rtmp) = model_rank_id!put this before increment and then we don't need the -1 everywhere
                rtmp = rtmp + 1
            END IF
        END DO
        !the loop above distributes models evenly amongst filters hence 'filter_models' is the number on this filter and is less than total_models
        filter_models = rtmp-1
        pf%count = filter_models
        pf%nens = total_models
        
        !initialise array to hold the initial receive test results
        ALLOCATE received(filter_models)
        received = .FALSE.

        CALL empi_wrap_set_model_count(filter_models, total_models)
        CALL empi_wrap_set_ensemble_members(filter_models) !doing the allgather
        !TODO
        !work out what the hell that ALLGATHER is doing in comms:initialise_mpi
        !there are no comments to explain as usual :(
        
        PRINT *, 'Total coupled processes = ', total_procs_coupled
        PRINT *, 'Number of filters = ', num_filters
        PRINT *, '-> Number of models = ', total_models
        PRINT *, 'Coupling rank id = ', filter_coupling_id
        PRINT *, "-> Particle filter local relative rank = ", m_filter_id_rank, " and I own"
        DO rtmp=1, filter_models
            PRINT*, " model particle: " , pf%particles(rtmp) !models(rtmp)
            PRINT*, 'end_model_run: ', end_model_runs(rtmp)
        END DO


        
    END SUBROUTINE empi_pf_assign_models
    
    
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_model_hookup()
        !get state vector dimensions, 
        !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        !TODO timeout/timestep switch, and whichever value is then needed
        !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        INTEGER :: count, state1_dim1, state1_dim2, &
                        state2_dim1, state2_dim2, &
                        state3_dim1, state3_dim2
        REAL, ALLOCATABLE, DIMENSION(:) :: state_dims
    INTEGER :: model_rank_id    !identifier of model instance
        !INTEGER, PARAMETER :: c_mpi_tag = 1 !constant


        !let's get the information on the state buffer which we'll need later
        !this assumes each state will have serialised 2D data, therefore there are
        !states_2D*2 numbers in the array, state_dims.
        ALLOCATE (state_dims(states_2D * 2), stat = err_code)
        PRINT *, "state size buffer allocated "

        !Do a loop like with the other send/receives so that we get stuff from the right number of models
        !but in fact we only really need one model to tell us.
        DO count = 1, filter_models
            model_rank_id= pf%particles(count) !models(count)
            !CALL MPI_RECV(state_dims, states_2D *2, MPI_INTEGER, &
            !                       model_rank_id, c_mpi_tag, coupling_comms, empi_status, empi_err)
                                   
            CALL empi_wrap_receive_int(state_dims, states_2D *2, model_rank_id)
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
        ALLOCATE(state_vector(num_state_params), stat = err_code)
        IF (err_code .ne. 0) WRITE(*, *) "Particle filter could not allocate state_vector."
        PRINT *, 'State vector created '

        DEALLOCATE(state_dims)

    END SUBROUTINE empi_pf_model_hookup

    !-----------------------------------------------------------------------------
    SUBROUTINE empi_pf_do_work()
        USE Qdata
        USE Rdata
        IMPLICIT NONE
        

        !random seed thingy
        CALL random_seed_mpi(coupling_id) !this is in gen_rand.f90
        !setup Q and R
        IF (.NOT. pf % gen_Q) THEN
            CALL loadQ !this just loads pf_control and sizes - seems a bit pointless
            CALL loadR !loads sizes and does a bit of stuff
        END IF

        IF (.NOT. pf % gen_Q) THEN
            !receive the model's initial state at time, t, after spin-up
            PRINT*, 'empi_pf_model_first_data____________'
            CALL empi_pf_model_first_data()

            !initial perturbation
            !we have data in pf % psi
            DO iter = 1, pf % count
                CALL perturb_particle(pf % psi(:, iter))
            END DO

            CALL empi_pf_do_known_timestep_filter()
            
            PRINT*, 'PF: finished the loop - now to tidy up'
            !send the last state. This would otherwise block the looping receives in the model
            PRINT*, 'empi_pf_model_last_data____________'
            CALL empi_pf_model_last_data()

        ELSE
            CALL genQ

        END IF

    END SUBROUTINE empi_pf_do_work
    !-----------------------------------------------------------------------------
    SUBROUTINE empi_pf_model_first_data()
        IMPLICIT NONE
        INTEGER :: model_rank_id !identifier of model instance

        INTEGER :: iter
        LOGICAL :: success

        ! Launch a non_blocking receive for all models attached to this process
        ! the 'requests' and 'statuses' arrays can stay local to the MPI module
        DO iter = 1, pf % count
            model_rank_id = pf % particles(iter)
            CALL empi_wrap_open_receive_dbl(pf % psi(:, iter), num_state_params, model_rank_id)

        END DO

        !now keep testing for results
        iter = 0
        DO
            !this is a sneaky way of repeating valid array indices between 1 and pf%count
            iter = MOD(iter, pf % count) + 1
            IF (.NOT.received(iter)) THEN
                model_rank_id = pf % particles(iter)
                CALL empi_wrap_test(model_rank_id, success)
                IF (success) THEN
                    received(iter) = .TRUE.
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
        
!    loop number obs
!        loop steps between obs
!           T+1
!           filters
!        end
!        T+1
!        filters
!    end
        CALL output_from_pf
        IF (pf % gen_data) CALL save_truth(pf % psi(:, 1))
        IF (pf % use_traj) CALL trajectories
        !start_t = mpi_wtime()

        DO iter1 = 1, pf % time_obs !number of observations to assimilate
            PRINT*, 'PF: observation counter = ', j
            DO iter2 = 1, pf % time_bwn_obs - 1
                pf % timestep = pf % timestep + 1

                SELECT CASE (pf % filter)
                CASE ("EW")
                    CALL proposal_filter
                CASE ("SI", "SE")
                    CALL stochastic_model 
                CASE ("ET")
                    CALL deterministic_model
                CASE DEFAULT
                    PRINT*, 'Error -555: Incorrect pf%filter'
                END SELECT

                IF (pf % use_traj) CALL trajectories
                CALL output_from_pf
            END DO !timesteps between observations

            pf % timestep = pf % timestep + 1
            PRINT*, 'starting the equal weight filter step'

                SELECT CASE (pf % filter)
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

            PRINT*, 'PF: timestep = ', pf % timestep, 'after equal weight filter'

            IF (pf % gen_data) CALL save_truth(pf % psi(:, 1))
            IF (pf % use_traj) CALL trajectories
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

        DO iter, pf%count
            model_rank_id = pf_particles(iter)
            CALL empi_wrap_open_send_dbl(pf%psi(:,iter), num_state_params, model_rank_id)
        END DO
        
        CALL empi_wrap_waitall(pf%count)
            
    END SUBROUTINE empi_pf_model_last_data
    !-----------------------------------------------------------------------------
    !-----------------------------------------------------------------------------
    
    SUBROUTINE empi_pf_do_filter()
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !need to do the random seed thingy... gen_rand:random_seed_mpi(pfid)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        IMPLICIT NONE

        INTEGER :: iter1, iter2, count, now
        INTEGER :: request_handle, index
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
        DO WHILE (.NOT. ALL(all_complete))!(all_complete .eqv. .FALSE.)
            PRINT *, '___________________LOOP__ ', count
            
            DO iter1 = 1, filter_models
                
                !Test each model in turn and ensure we do the comms until it has finished
                IF (end_model_runs(iter1) .EQV. .FALSE.) THEN                        
                    model_rank_id= pf%particles(iter1) !models(iter1)
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
                        
                        !CALL MPI_RECV(state_vector, num_state_params, MPI_DOUBLE_PRECISION, &
                        !              model_rank_id, c_tag, coupling_comms, empi_status, empi_err)
                        CALL empi_wrap_receive_dbl(state_vector, num_state_params, model_rank_id)
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
                    model_rank_id= pf%particles(iter1) !models(iter1)                  
                    !If it's still running, send the state vectors back
                    IF (end_model_runs(iter1) .eqv. .FALSE.) THEN
                        !CALL MPI_SEND(state_vector, num_state_params, MPI_DOUBLE_PRECISION, &
                        !                       model_rank_id, c_tag, coupling_comms, empi_err)
                        CALL empi_wrap_send_dbl(state_vector, num_state_params, model_rank_id)
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
!-----------------------------------------------------------------------------


    SUBROUTINE empi_pf_cleanup()
        !cleanup
        CALL empi_wrap_finalise !MPI_FINALIZE(empi_err)
        DEALLOCATE(state_vector)
        DEALLOCATE(pf%particles) !(models)
        DEALLOCATE(end_model_runs)
        DEALLOCATE(received)

    END SUBROUTINE empi_pf_cleanup


END MODULE empire_mod