!     
! File:   empi_wrappers.F90
! Author: Jane
!
! Created on 03 February 2015, 15:39
!
! This module's purpose is to remove the complexity of the MPI interface from the main code logic.
! It means that the MPI header file need only to be included once, the coupling details are local,
! as are constants, error codes and status codes.

MODULE empi_wrappers
    
    IMPLICIT NONE
    INCLUDE 'mpif.h'

 
! These blocks can go into the modules making use of this one.    
!    INTERFACE        
!        SUBROUTINE empi_wrap_allgatherv(sendbuf, sendcount, recvbuf, multiplier)
!            INTEGER, INTENT(IN) :: sendcount
!            INTEGER, OPTIONAL, INTENT(IN) :: multiplier
!            TYPE(C_PTR) :: sendbuf, recvbuf
!        END SUBROUTINE empi_wrap_allgatherv
!        
!        SUBROUTINE empi_wrap_gatherv(sendbuf, sendcount, recvbuf, root_proc, multiplier)
!            INTEGER, INTENT(IN) :: sendcount, root_proc
!            TYPE(C_PTR) :: sendbuf, recvbuf
!            INTEGER, OPTIONAL, INTENT(IN) :: multiplier
!        END SUBROUTINE empi_wrap_gatherv
!    
!        SUBROUTINE empi_wrap_scatterv(sendbuf, recvbuf, recvcount, root_proc, multiplier)
!            INTEGER, INTENT(IN) :: recvcount, root_proc
!            TYPE(C_PTR) :: sendbuf, recvbuf
!            INTEGER, OPTIONAL, INTENT(IN) :: multiplier
!        END SUBROUTINE empi_wrap_scatterv
!    END INTERFACE    
    
    !This is a private helper function
    !PRIVATE, INTEGER :: get_type
    
    !Constants for use by calling functions
    ENUM, BIND(C)
        ENUMERATOR :: TYPE_INT = 1, TYPE_DBL = 2
    END ENUM        
    ENUM, BIND(C)
        ENUMERATOR :: OPERATOR_SUM = 1, OPERATOR_MEAN = 2
    END ENUM
    ENUM, BIND(C)
        ENUMERATOR :: ANY_SOURCE = -1
    END ENUM
    ENUM, BIND(C)
        ENUMERATOR :: FILTER_COMMS = 1, COUPLE_COMMS = 2 !for pf_mpi_comm, coupling_comms
    END ENUM
    
    INTEGER :: MPI_MODES
    !Variable declarations for the module
    INTEGER :: empi_err
    INTEGER :: empi_status(MPI_STATUS_SIZE)
    INTEGER :: pf_mpi_comm !the communicator between DA processes
    INTEGER :: coupling_comms !the comminicator between DA and models
    INTEGER :: filter_rank !the rank of this process on pf_mpi_comm
    INTEGER :: filter_coupling_id !the rank of this process on coupling_comms
    INTEGER :: total_models !the total number of model processes (head nodes)
    INTEGER :: num_filters !the total number of DA processes
    INTEGER :: total_procs_coupled !the number of processes on the model/filter coupler i.e. all filters and all model head nodes
    INTEGER, PARAMETER :: c_mpi_tag = 1 !constant
    INTEGER, DIMENSION(:), ALLOCATABLE :: requests, requests_all  !an array of requests sent indexed by model number
    INTEGER, DIMENSION(:, :), ALLOCATABLE :: statuses, statuses_all  !and one for the mpi return status for each
    INTEGER, ALLOCATABLE, DIMENSION(:) :: gblcount !the number of ensemble members associated with each DA process
    INTEGER, ALLOCATABLE, DIMENSION(:) :: gbldisp !the displacements of each ensemble member relative to filter_id=0. used for mpi_gatherv and mpi_scatterv on the communicator between DA processes
    

    
CONTAINS

    SUBROUTINE empi_wrap_initialise()
        IMPLICIT NONE
        INTEGER :: world_id !the rank of this process on MPI_COMM_WORLD
        INTEGER :: world_size

        INTEGER, PARAMETER :: c_coupling_colour = 9999 , c_pf_colour = 10000 !unused communicator

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
        
        total_models = total_procs_coupled - num_filters

        PRINT *, "PF empi_pf_initialise: Total coupled processes ", total_procs_coupled, &
        ", coupling rank id ", filter_coupling_id, &
        ", coupling comm id ", coupling_comms

    END SUBROUTINE empi_wrap_initialise
!-----------------------------------------------------------------------------
!    SUBROUTINE empi_wrap_all_reduce_int_sum(sendbuf, recvbuf, size)
!        IMPLICIT NONE
!
!        INTEGER, INTENT(IN) :: sendbuf, size
!        INTEGER, INTENT(OUT) :: recvbuf
!
!        CALL MPI_ALLREDUCE(sendbuf, recvbuf, size, MPI_INTEGER, MPI_SUM, coupling_comms, empi_err)
!
!    END SUBROUTINE empi_wrap_all_reduce
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_set_model_count(filter_models, total_models)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: filter_models, total_models
        
        ALLOCATE (requests(filter_models))
        ALLOCATE (requests_all(total_models))
        ALLOCATE (statuses(mpi_status_size, filter_models))
        ALLOCATE (statuses_all(mpi_status_size, total_models))
    END SUBROUTINE empi_wrap_set_model_count
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_set_ensemble_members(count)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: count
        INTEGER :: iter
        ALLOCATE (gblcount(num_filters))
        ALLOCATE (gbldisp(num_filters))

        !very dubious about this - according to docs, send & recv should be same type
        !and they're not: integer vs array of int
        CALL MPI_ALLGATHER(count, 1, MPI_INTEGER, gblcount, 1, MPI_INTEGER, &
                            pf_mpi_comm, empi_err)

        gbldisp = 0
        IF (count .GT. 1) THEN
            DO iter = 2, count
                gbldisp(iter) = gbldisp(iter - 1) + gblcount(iter - 1)
            END DO
        END IF
    END SUBROUTINE empi_wrap_set_ensemble_members
!-----------------------------------------------------------------------------
    !see if this will work...
    !https://groups.google.com/forum/#!msg/comp.lang.fortran/Pid_eF9Wp6o/rODpb4fQ6YoJ
    !http://www.cs.mtu.edu/~shene/COURSES/cs201/NOTES/chap03/select.html
    !http://stackoverflow.com/questions/23995529/select-type-with-unlimited-polymorphic-pointer-to-character-variable
    !if not, then maybe define my own wrapper constants for type to pass in.
!    SELECT CASE (passedin thingy)
!        const1
!           set MPI type
!        etc
!   Given up on this tack - cannot get the pointer to be correct, tried:
!    TYPE(C_PTR) :: recvbuf  !with 'USE iso_c_binding' at the top
!    TYPE :: sendbuf(*)
!    POINTER :: sendbuf
!    POINTER, TYPE(*) :: sendbuf
!-----------------------------------------------------------------------------
!    SUBROUTINE empi_wrap_receive(recvbuf, size, source, the_type)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: size, source, the_type
!        TYPE(C_PTR) :: recvbuf
!        INTEGER :: empi_type
!        
!!        SELECT TYPE (sendbuf)
!!            TYPE is (INTEGER)
!!                empi_type = MPI_INTEGER
!!            TYPE is (DOUBLE)
!!                empi_type = MPI_DOUBLE_PRECISION
!!        END SELECT
!        SELECT CASE (the_type)
!                CASE (TYPE_INT)
!                        empi_type = MPI_INTEGER
!                CASE (TYPE_DBL)
!                        empi_type = MPI_DOUBLE_PRECISION
!        END SELECT
!        
!        CALL MPI_RECV(recvbuf, size, empi_type, &
!            source, c_mpi_tag, coupling_comms, empi_status, empi_err)
!    END SUBROUTINE empi_wrap_receive
!-----------------------------------------------------------------------------
        SUBROUTINE empi_wrap_receive_int(recvbuf, size, source)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, source
        POINTER, TYPE(INTEGER) :: recvbuf
        
        CALL MPI_RECV(recvbuf, size, MPI_INTEGER, &
            source, c_mpi_tag, coupling_comms, empi_status, empi_err)
    END SUBROUTINE empi_wrap_receive_int
!-----------------------------------------------------------------------------
        SUBROUTINE empi_wrap_receive_dbl(recvbuf, size, source)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, source
        POINTER, DOUBLE :: recvbuf
        
        CALL MPI_RECV(recvbuf, size, MPI_DOUBLE_PRECISION, &
            source, c_mpi_tag, coupling_comms, empi_status, empi_err)
    END SUBROUTINE empi_wrap_receive_dbl
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_open_receive_dbl(recvbuf, size, source)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, source
        POINTER, DOUBLE :: recvbuf               
        
        CALL MPI_IRECV(recvbuf, size, MPI_DOUBLE_PRECISION, &
            source, c_mpi_tag, coupling_comms, requests(source), empi_err)
        
    END SUBROUTINE empi_wrap_open_receive_dbl
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_open_receive_da_dbl(recvbuf, size, tag, request_index)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, request_index, tag, the_type
        POINTER, DOUBLE :: recvbuf
        
        CALL MPI_IRECV(recvbuf, size, MPI_DOUBLE_PRECISION, &
            MPI_ANY_SOURCE, tag, pf_mpi_comm, requests(request_index), empi_err)
        
    END SUBROUTINE empi_wrap_open_receive_da_dbl
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_send_dbl(sendbuf, size, tgt)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, tgt
        POINTER, DOUBLE :: sendbuf
       
        CALL MPI_SEND(sendbuf, size, MPI_DOUBLE_PRECISION, &
                tgt, c_mpi_tag, coupling_comms, empi_err)
    END SUBROUTINE empi_wrap_send_dbl
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_open_send_dbl(sendbuf, size, targetid)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, targetid
        POINTER, DOUBLE :: sendbuf

        CALL MPI_ISEND(sendbuf, size, MPI_DOUBLE_PRECISION, targetid, &
        c_mpi_tag, coupling_comms, requests(targetid), empi_err)
        
    END SUBROUTINE empi_wrap_open_send_dbl
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_open_send_da(sendbuf, size, targetid, tag, request_index, the_type)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: size, targetid, tag, request_index, the_type
        TYPE(C_PTR) :: sendbuf
        INTEGER :: empi_type
        
!        SELECT TYPE (sendbuf)
!            TYPE is (INTEGER)
!                empi_type = MPI_INTEGER
!            TYPE is (DOUBLE)
!                empi_type = MPI_DOUBLE_PRECISION
!        END SELECT
        SELECT CASE (the_type)
                CASE (TYPE_INT)
                        empi_type = MPI_INTEGER
                CASE (TYPE_DBL)
                        empi_type = MPI_DOUBLE_PRECISION
        END SELECT
        
        CALL MPI_ISEND(sendbuf, size, empi_type, targetid, &
        tag, pf_mpi_comm, requests(request_index), empi_err)
        
    END SUBROUTINE empi_wrap_open_send_da
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_test(source, success)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: source
        LOGICAL, INTENT(OUT) :: success

        CALL MPI_TEST(requests(source), success, statuses(:,source), empi_err)
        
    END SUBROUTINE empi_wrap_test
!-----------------------------------------------------------------------------
! Tests for the completion of all previously initiated communications in a list. 
    SUBROUTINE empi_wrap_testall(count, success)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: count
        LOGICAL, INTENT(OUT) :: success

        CALL MPI_TESTALL(count, requests_all(1:count), success, statuses_all, empi_err)
        
    END SUBROUTINE empi_wrap_testall
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_waitall(count)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: count
        
        CALL MPI_WAITALL(count, requests, statuses, empi_err)
    END SUBROUTINE empi_wrap_waitall
!-----------------------------------------------------------------------------
!Gathers data from all processes and delivers it to all. Each process may 
!contribute a different amount of data. 
    SUBROUTINE empi_wrap_allgatherv_dbl(sendbuf, sendcount, recvbuf, multiplier)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: sendcount
        INTEGER, OPTIONAL, INTENT(IN) :: multiplier
        POINTER, DOUBLE :: sendbuf, recvbuf
        INTEGER :: multi

        IF (present(multiplier)) THEN
            multi = multiplier
        ELSE
            multi = 1
        END IF
        
        CALL MPI_ALLGATHERV(sendbuf, sendcount, MPI_DOUBLE_PRECISION, &
            recvbuf, gblcount*multi, gbldisp*multi, empi_type, &
            pf_mpi_comm, empi_err)
            
    END SUBROUTINE empi_wrap_allgatherv_dbl
!-----------------------------------------------------------------------------
!Gathers varying amounts of data from all processes to the root process 
    SUBROUTINE empi_wrap_gatherv_dbl(sendbuf, sendcount, recvbuf, root_proc, multiplier)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: sendcount, root_proc
        POINTER:DOUBLE :: sendbuf, recvbuf
        INTEGER, OPTIONAL, INTENT(IN) :: multiplier
        INTEGER :: multi
        
        IF (present(multiplier))THEN
            multi = multiplier
        ELSE
            multi = 1
        END IF
        
        CALL MPI_GATHERV (sendbuf, sendcount, empi_type, recvbuf, &
            gblcount*multi, gbldisp*multi, MPI_DOUBLE_PRECISION, root_proc, &
            pf_mpi_comm, empi_err)
            
   END SUBROUTINE empi_wrap_gatherv_dbl
!-----------------------------------------------------------------------------            
   SUBROUTINE empi_wrap_scatterv_dbl(sendbuf, recvbuf, recvcount, root_proc, multiplier)
        INTEGER, INTENT(IN) :: recvcount, root_proc
        POINTER, DOUBLE :: sendbuf, recvbuf
        INTEGER, OPTIONAL, INTENT(IN) :: multiplier
        INTEGER :: multi
        
        IF (present(multiplier)) THEN
            multi = multiplier
        ELSE
            multi = 1
        END IF
        
        CALL MPI_SCATTERV(sendbuf, gblcount*multi, gbldisp*multi, MPI_DOUBLE_PRECISION, &
            recvbuf, recvcount, empi_type, root_proc, &
            pf_mpi_comm, empi_err)
        
    END SUBROUTINE empi_wrap_scatterv_dbl
    
!-----------------------------------------------------------------------------
!  Reduces values on all processes within a group.
    SUBROUTINE empi_wrap_reduce_dbl(sendbuf, recvbuf, count, operation, root_process)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: count, operation, root_process, the_type
        POINTER, DOUBLE :: sendbuf, recvbuf
        INTEGER :: empi_operator

        SELECT CASE (operation)
            CASE (OPERATOR_SUM)
                empi_operator = MPI_SUM
        END SELECT
                
        CALL MPI_REDUCE(sendbuf, recvbuf, count, MPI_DOUBLE_PRECISION, MPI_SUM, &
            root_process, pf_mpi_comm, empi_err)
            
    END SUBROUTINE empi_wrap_reduce_dbl
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_allreduce_in_place(buf, count, operation, the_type)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: count, operation, the_type
        TYPE(C_PTR) :: buf
        INTEGER :: empi_type, empi_operator
        
!        SELECT TYPE (sendbuf)
!            TYPE is (INTEGER)
!                empi_type = MPI_INTEGER
!            TYPE is (DOUBLE)
!                empi_type = MPI_DOUBLE_PRECISION
!        END SELECT
        SELECT CASE (the_type)
                CASE (TYPE_INT)
                        empi_type = MPI_INTEGER
                CASE (TYPE_DBL)
                        empi_type = MPI_DOUBLE_PRECISION
        END SELECT
        
        SELECT CASE (operation)
            CASE (OPERATOR_SUM)
                empi_operator = MPI_SUM
        END SELECT
        
        CALL MPI_ALLREDUCE(MPI_IN_PLACE, buf, count, empi_type, empi_operator,&
            pf_mpi_comm, empi_err)
            
    END SUBROUTINE empi_wrap_allreduce_in_place
!-----------------------------------------------------------------------------
!Broadcasts a message from the process with rank root to all other processes of the group.
    SUBROUTINE empi_wrap_broadcast_int(sendbuf, sendcount, root_proc)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: sendcount, root_proc
        POINTER, INTEGER :: sendbuf
        
        CALL MPI_BCAST(sendbuf, sendcount, MPI_INTEGER, root_proc, &
            pf_mpi_comm, empi_err)
        
    END SUBROUTINE empi_wrap_broadcast_int
!-----------------------------------------------------------------------------
    SUBROUTINE empi_wrap_barrier()
        CALL MPI_BARRIER(pf_mpi_comm, empi_err)
    END SUBROUTINE empi_wrap_barrier
!-----------------------------------------------------------------------------
!    SUBROUTINE empi_wrap_receive_int(recvbuf, size, source)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: size, source
!        POINTER, INTEGER :: recvbuf
!
!        CALL MPI_RECV(recvbuf, size, MPI_INTEGER, &
!        source, c_mpi_tag, coupling_comms, empi_status, empi_err)
!    END SUBROUTINE empi_wrap_receive_int
!!-----------------------------------------------------------------------------
!    SUBROUTINE empi_wrap_receive_dbl(recvbuf, size, source)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: size, source
!        POINTER, DOUBLE :: recvbuf
!
!        CALL MPI_RECV(recvbuf, size, MPI_DOUBLE_PRECISION, &
!                source, c_mpi_tag, coupling_comms, empi_status, empi_err)
!    END SUBROUTINE empi_wrap_receive_int
!!-----------------------------------------------------------------------------
!    SUBROUTINE empi_wrap_send_dbl(sendbuf, size, source)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: size, source
!        POINTER, DOUBLE :: sendbuf
!
!        CALL MPI_SEND(sendbuf, size, MPI_DOUBLE_PRECISION, &
!                source, c_mpi_tag, coupling_comms, empi_err)
!    END SUBROUTINE empi_wrap_receive_int
!-----------------------------------------------------------------------------    
    SUBROUTINE empi_wrap_finalise()
        CALL MPI_FINALIZE(empi_err)
        DEALLOCATE (requests)
        DEALLOCATE (statuses)
    END SUBROUTINE empi_wrap_finalise
    
    
!    PURE FUNCTION get_type(the_type)
!        IMPLICIT NONE
!        INTEGER, INTENT(IN) :: the_type
!        INTEGER :: empi_type
!! Run-time polymorphism doesn't appear to work in F90. Wretched dinosaur language, grrr.        
!!        SELECT TYPE (sendbuf)
!!            TYPE is (INTEGER)
!!                empi_type = MPI_INTEGER
!!            TYPE is (DOUBLE)
!!                empi_type = MPI_DOUBLE_PRECISION
!!        END SELECT
!        SELECT CASE (the_type)
!                CASE (TYPE_INT)
!                        empi_type = MPI_INTEGER
!                CASE (TYPE_DBL)
!                        empi_type = MPI_DOUBLE_PRECISION
!        END SELECT
!        get_type = empi_type
!        RETURN 
!    END FUNCTION get_type
    
END MODULE empi_wrappers
