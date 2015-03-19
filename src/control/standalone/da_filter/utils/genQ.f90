!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-22 15:45:36 pbrowne>
!!!
!!!    A subroutine to estimate Q from a long model run
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
!> Subroutine to estimate Q from a long model run
subroutine genQ
  !use sizes
  use pf_control
  use empire_mod !comms
  use Qdata
  use hqht_plus_r
  implicit none
  !include 'mpif.h'
  integer, parameter :: rk=(kind(1.0d0))

!  integer, dimension(a_nxn,a_nyn,a_levels) :: a_u_vec,a_v_vec,a_theta_vec,a_q_vec
!  integer, dimension(a_nxn,a_nyn) :: a_pstar_vec
!  integer, dimension(o_nxn,o_nyn,o_levels) :: o_u_vec,o_v_vec,o_theta_vec,o_sal_vec
  integer :: i,j,k,count,radius,nnz
  integer, parameter :: n = 426655238
  integer, allocatable, dimension(:) :: row,col
  real(kind=rk), allocatable, dimension(:) :: val
  integer :: ne,iter,day,tag,mpi_err,days
  integer :: mpi_status(MPI_STATUS_SIZE)
  real(kind=rk), dimension(model_params%state_dim) :: x 
  real(kind=rk) :: start_t,end_t

  STOP 'IS GENQ DESIGNED FOR THIS SPECIFIC MODEL??'

  allocate(pf%mean(model_params%state_dim),row(n),col(n),val(n))
  print*,'Going to generate Q ~ the model error covariance matrix'



  start_t = mpi_wtime()
  !let us calculate the position vectors: put everything in its place:
  count = 0
!!$  !first is PSTAR (2d)
!!$  do j = 1,a_nyn
!!$     do i = 1,a_nxn
!!$        count = count + 1
!!$        a_pstar_vec(i,j) = count
!!$     end do
!!$  end do
!!$  print*,'a_pstar finishes on the ',count,' element'
!!$  !second is  U
!!$  do k = 1,a_levels
!!$     do j = 1,a_nyn
!!$        do i = 1,a_nxn
!!$           count = count + 1
!!$           a_u_vec(i,j,k) = count
!!$        end do
!!$     end do
!!$  end do
!!$  print*,'a_u finishes on the ',count,' element'
!!$  !FINISHED COMPUTING THE VECTOR POSITIONS OF EACH COMPONENT
!!$  print*,'FINISHED COMPUTING THE VECTOR POSITIONS OF EACH COMPONENT'

  end_t = mpi_wtime()
  print*,'and it took ',end_t-start_t,' seconds'
!  stop
  !initialise the value thingy
  val = 0.0_rk

  !initialise the mean
  pf%mean = 0.0_rk
  tag = 1
  print*,'first recv...'
  call mpi_recv(pf%psi(:,1), model_params%state_dim, MPI_DOUBLE_PRECISION, &
       0, tag, CPL_MPI_COMM,mpi_status, mpi_err)
  print*,'.............recvd'
  
  !loop for 5 years:
  days = 360*5
  do day = 1,days
     start_t = mpi_wtime()
     print*,'day = ',day
     !update by a day, or 72 iterations:
     do iter = 1,72
        call mpi_send(pf%psi(:,1), model_params%state_dim , MPI_DOUBLE_PRECISION, &
             0, tag, CPL_MPI_COMM, mpi_err)
        print*,'iter = ',iter
        call mpi_recv(pf%psi(:,1), model_params%state_dim, MPI_DOUBLE_PRECISION, &
             0, tag, CPL_MPI_COMM,mpi_status, mpi_err)

     end do
     end_t = mpi_wtime()
     print*,'72 model runs on ',day,' took ',end_t-start_t,' seconds'
     !update the mean:
     pf%mean = pf%mean + pf%psi(:,1)

     x = pf%psi(:,1)

     radius = 2
     start_t = mpi_wtime()
               
     !  call genQ_at2d_at2d(ne,row,col,val,n,radius,field1,field2,x)

     print*,'finished genQ_order at the end of the day with ne = ',ne
     end_t = mpi_wtime()
     print*,'generating Q on ',day,' took ',end_t-start_t,' seconds.'

  end do !end of the daily loop



  !now we should divide the value by N-1 = 360*5-1...
  val = val/real(days-1,rk)


  !now subtract off the mean components:
  pf%mean = pf%mean/real(days,rk)


  val = -1.0_rk*val
  x =  pf%mean

!  call genQ_at2d_at2d(ne,row,col,val,n,radius,field1,field2,x)

  val = -1.0_rk*val


  print*,'MAXIMUM COVARIANCE VALUE = ',maxval(val)
  print*,'MIMINUM COVARIANCE VALUE = ',minval(val)
  print*,'MINIMUM ABSOLUTE C VALUE = ',minval(abs(val))







  !do the final mpi send to end model cleanly
  call mpi_send(pf%psi(:,1), model_params%state_dim , MPI_DOUBLE_PRECISION, &
       0, tag, CPL_MPI_COMM, mpi_err)


  !now let us scale Q
  val = val/1.0D4







  nnz = 0
  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) nnz = nnz + 1
  end do

  open(2,file='Qdata.dat',action='write',status='replace',form='unformatted')
  write(2) model_params%state_dim
  write(2) nnz

  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) write(2) val(i)
  end do
  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) write(2) row(i)
  end do
  do i = 1,ne
     if(row(i) .eq. col(i) .or. abs(val(i)) .gt. 1.0D-13) write(2) col(i)
  end do
  close(2)
  print*,'finished generating Q'

  deallocate(pf%mean,row,col,val)



  


end subroutine genQ


!!$subroutine genQ_at2d_at2d(ne,row,col,val,n,radius,field1,field2,x)
!!$  use hadcm3_config
!!$  use sizes
!!$  implicit none
!!$  integer, parameter :: rk = kind(1.0d0)
!!$  integer, intent(inout) :: ne
!!$  integer, intent(in) :: n,radius
!!$  integer, dimension(n), intent(inout) :: row,col
!!$  real(kind=rk), dimension(n), intent(inout) :: val
!!$  real(kind=rk), dimension(model_params%state_dim), intent(in) :: x
!!$  integer, dimension(a_nxn,a_nyn), intent(in) :: field1,field2
!!$  integer :: i,j,ii,jj
!!$
!!$  do j = 1,a_nyn
!!$     do i = 1,a_nxn
!!$        do jj = max(1,j-radius),min(a_nyn,j+radius)
!!$           do ii = max(1,i-radius),min(a_nxn,i+radius)
!!$              if(field1(i,j) .le. field2(ii,jj)) then
!!$                 ne = ne + 1
!!$                 row(ne) = field1(i,j)
!!$                 col(ne) = field2(ii,jj)
!!$                 val(ne) = val(ne)+x(row(ne))*x(col(ne))
!!$              end if
!!$           end do
!!$           if(i .le. radius) then
!!$              do ii = a_nxn-radius+i,a_nxn
!!$                 if(field1(i,j) .le. field2(ii,jj)) then
!!$                    ne = ne + 1
!!$                    row(ne) = field1(i,j)
!!$                    col(ne) = field2(ii,jj)
!!$                    val(ne) = val(ne)+x(row(ne))*x(col(ne))
!!$                 end if
!!$              end do
!!$           end if
!!$           if(i+radius .gt. a_nxn) then
!!$              do ii = 1,i+radius-a_nxn
!!$                 if(field1(i,j) .le. field2(ii,jj)) then
!!$                    ne = ne + 1
!!$                    row(ne) = field1(i,j)
!!$                    col(ne) = field2(ii,jj)
!!$                    val(ne) = val(ne)+x(row(ne))*x(col(ne))
!!$                 end if
!!$              end do
!!$           end if
!!$        end do
!!$     end do
!!$  end do
!!$end subroutine genQ_at2d_at2d


