!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-26 17:31:05 pbrowne>
!!!
!!!    Subroutine to perform nudging in the proposal step of EWPF
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!subroutine proposal_filter based on subroutine IntegrateModel
!given from the Particle filter code of Mel and Peter Jan
!PAB 04-02-2013

!> Subroutine to perform nudging in the proposal step of EWPF
subroutine proposal_filter
  use pf_control
 !use Sizes
  use empire_mod !comms

  IMPLICIT NONE
  !include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  real(kind=rk) :: pWeight!, qWeight
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: normaln     !vector to store uncorrelated random error
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: betan       !vector to store sqrtQ correlated random error
  real(kind=rk), dimension(model_params%obs_dim) :: y             !y, the observations
  real(kind=rk), dimension(model_params%obs_dim,pf%count) :: Hpsi          !H(psi^(n-1))
  real(kind=rk), dimension(model_params%obs_dim,pf%count) :: y_Hpsin1      !y-H(psi^(n-1))
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: fpsi        !f(psi^(n-1))
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: Qkgain
  real(kind=rk) :: t
  integer :: particle,k,tag,mpi_err
  integer :: mpi_status( MPI_STATUS_SIZE )
  real(kind=rk), dimension(7) :: ti
  logical, parameter :: time = .false.


  !get the next observations and store it in vector y
  if(.not. pf%gen_data) call get_observation_data(y)


  !compute y - H(x)
  if(.not. pf%gen_data) then
     if(time) t = mpi_wtime()
     call H(model_params%obs_dim,pf%count,pf%psi,Hpsi,pf%timestep)
     if(time) ti(1) = mpi_wtime()-t
     !$omp parallel do
     do k = 1,pf%count
        y_Hpsin1(:,k) = y - Hpsi(:,k)
     end do
     !$omp end parallel do
     if(time) ti(2) = mpi_wtime()-ti(1) -t
  else
     y_Hpsin1 = 0.0_rk
  end if

  !get the model to provide f(x)
  do k =1,pf%count
     particle = pf%particles(k)
     tag = 1
     PRINT *, 'using particle id ', particle-1, ' and communicator ', CPL_MPI_COMM
     call mpi_send(pf%psi(:,k),model_params%state_dim,MPI_DOUBLE_PRECISION&
          &,particle-1,tag,CPL_MPI_COMM,mpi_err)
  end do
  if(time) ti(3) = mpi_wtime()-ti(2)-t
  DO k = 1,pf%count
     particle = pf%particles(k)
     tag = 1
     CALL MPI_RECV(fpsi(:,k), model_params%state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
  END DO
  if(time) ti(4) = mpi_wtime()-ti(3) -t


  !draw from a Gaussian for the random noise
  call NormalRandomNumbers2D(0.0D0,1.0D0,model_params%state_dim,pf%count,normaln)
  if(time) ti(5) = mpi_wtime()-ti(4) -t
  
  !compute the relaxation term Qkgain, the intermediate
  !term kgain and apply correlation to noise
  call Bprime(y_Hpsin1,kgain,Qkgain,normaln,betan)
  if(time) ti(6) = mpi_wtime()-ti(5) -t


  !update the new state and weights based on these terms
  !$omp parallel do private(particle,pweight)
  DO k = 1,pf%count
     particle = pf%particles(k)
!     print*,'|fpsi-psi|_2 = ',dnrm2(model_params%state_dim,(fpsi(:,k)-pf%psi(:,k)),1)
     call update_state(pf%psi(:,k),fpsi(:,k),Qkgain(:,k),betan(:,k))
     pweight = sum(Qkgain(:,k)*kgain(:,k))+2.0D0*sum(betan(:,k)*kgain(:,k))
     pf%weight(particle) = pf%weight(particle) + 0.5*pWeight
  end DO
  !$omp end parallel do
  if(time) then
     ti(7) = mpi_wtime()-ti(6) -t
     print*,ti
  end if
end subroutine proposal_filter
