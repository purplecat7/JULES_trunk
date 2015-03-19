!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-11-12 13:42:46 pbrowne>
!!!
!!!    Module and subroutine to intitalise EMPIRE coupling to models
!!!    Copyright (C) 2014  Philip A. Browne
!!!
!!!    This program is free software: you can redistribute it and/or modify
!!!    it under the terms of the GNU General Public License as published by
!!!    the Free Software Foundation, either version 3 of the License, or
!!!    (at your option) any later version.
!!!
!!!    This program is distributed in the hope that it will be useful,
!!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!!    GNU General Public License for more details.
!!!
!!!    You should have received a copy of the GNU General Public License
!!!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!!!
!!!    Email: p.browne @ reading.ac.uk
!!!    Mail:  School of Mathematical and Physical Sciences,
!!!    	      University of Reading,
!!!	      Reading, UK
!!!	      RG6 6BB
!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!> Module containing EMPIRE coupling data
module comms
  integer :: CPL_MPI_COMM !< the communicator between the empire
  !< codes and the model master nodes
  integer :: mype_id !< the rank of this process on MPI_COMM_WORLD
  integer :: myRank !< the rank of this process on CPL_MPI_COMM
  integer :: nProc !< the total number of processes
  integer :: pf_mpi_comm !< the communicator between DA processes
  integer :: pfrank      !< the rank of this process on PF_MPI_COMM
  integer :: npfs        !< the total number of DA processes
  integer, allocatable, dimension(:) :: gblcount !< the number of
  !< ensemble members associated with each DA process
  integer, allocatable, dimension(:) :: gbldisp !< the displacements
  !< of each each ensemble member relative to pfrank=0. VERY useful
  !< for mpi_gatherv and mpi_scatterv on pf_mpi_comm
  
contains
  
  subroutine allocate_data

    implicit none
    
  end subroutine allocate_data
  
  subroutine deallocate_data
    implicit none

  end subroutine deallocate_data

  !> subroutine to make EMPIRE connections and saves details into
  !! pf_control module
  subroutine initialise_mpi

    use pf_control
    implicit none
    include 'mpif.h'
    
    integer :: mpi_err!dummy_colour,mpi_err
    integer :: couple_colour !DUMMY_MPI_COMMUNICATOR,couple_colour
    !integer :: couple_mype_id,couple_root
    integer :: rtmp!,ctmp

 !   integer :: tag!,state_dim!,iter
 !   integer :: num_iters
    integer :: particle,world_id
    integer :: myrank !nproc,myrank
!    integer :: mpi_status(MPI_STATUS_SIZE)
    integer :: nens,i
    integer :: da
    integer :: count,pf_colour
    integer :: world_size
    
    pf_colour = 10000
    couple_colour=9999
    call MPI_INIT (mpi_err)
    
    da = 1
    call MPI_COMM_RANK (MPI_COMM_WORLD,world_id,     mpi_err)
    call mpi_comm_size (mpi_comm_world,world_size,   mpi_err)
    call mpi_comm_split(mpi_comm_world,da,           world_id,  pf_mpi_comm, mpi_err)
    call mpi_comm_rank (pf_mpi_comm,   pfrank,       mpi_err)
    call mpi_comm_size (pf_mpi_comm,   npfs,          mpi_err)
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,couple_colour,world_size,CPL_MPI_COMM,mpi_err)
    call MPI_COMM_RANK (CPL_MPI_COMM,  myRank,       mpi_err)
    call MPI_COMM_SIZE (CPL_MPI_COMM,  nens,         mpi_err)
    
    nens = nens-npfs
    print*,'DA'
    print*,'nens = ',nens
    print*,'npfs = ',npfs
    
    
    !lets find the particles:
    count = 0
    do particle = 1,nens
       if( real(particle-1) .ge. real(nens*(pfrank))/real(npfs) .and.&
            & real(particle-1) .lt. real(nens*(pfrank+1))/real(npfs)) then
          count = count + 1
       end if
    end do
    
    allocate(pf%particles(count))
    count = 0
    do particle = 1,nens
       if(real(particle-1) .ge. real(nens*(pfrank))/real(npfs) .and.&
            & real(particle-1) .lt. real(nens*(pfrank+1))/real(npfs))&
            & then
          count = count + 1
          pf%particles(count) = particle
       end if
    end do
    
    allocate(gblcount(npfs))
    allocate(gbldisp(npfs))
!    print*,'woohoo allgather'
!    print*,count
!    print*,gblcount
!    print*,pf_mpi_comm
    call mpi_allgather(count,1,mpi_integer,gblcount,1,mpi_integer&
         &,pf_mpi_comm,mpi_err)
!    print*,'allgather did not break'
    gbldisp = 0
    if(npfs .gt. 1) then
       do i = 2,npfs
          gbldisp(i) = gbldisp(i-1) + gblcount(i-1)
       end do
    end if
    pf%count = count

    pf%nens = nens
    PRINT*,'PF_rank = ',pfrank,' and I own particles ',pf%particles

    
  end subroutine initialise_mpi

end module comms
