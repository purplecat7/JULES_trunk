!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-29 16:18:17 pbrowne>
!!!
!!!    subroutine to simply move the model forward in time one timestep
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
!> @brief
!> subroutine to simply move the model forward in time one 
!> timestep
!> PAB 21-05-2013

subroutine stochastic_model
  use pf_control
  !use Sizes
  use empire_mod !comms

  IMPLICIT NONE
  !include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  real(kind=rk), dimension(model_params%state_dim,pf%count) :: normaln     !vector to store uncorrelated random error
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: betan       !vector to store sqrtQ correlated random error
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: fpsi        !f(psi^(n-1))
!  real(kind=rk), dimension(model_params%state_dim,pf%count) :: kgain       !QH^T(HQH^T+R)^(-1)(y-H(psi^(n-1)))
  real(kind=rk), dimension(model_params%state_dim,pf%count) :: Qkgain
!  real(kind=rk) :: t,dnrm2
  integer :: particle,k,tag,mpi_err
  integer :: mpi_status( MPI_STATUS_SIZE )
  logical, parameter :: checkscaling=.true.

!  print*,'huh',pf%psi(:,1)
  do k =1,pf%count
     particle = pf%particles(k)
     tag = 1
     call mpi_send(pf%psi(:,k),model_params%state_dim,MPI_DOUBLE_PRECISION&
          &,particle-1,tag,CPL_MPI_COMM,mpi_err)
  end do
  DO k = 1,pf%count
     particle = pf%particles(k)
     tag = 1
     CALL MPI_RECV(fpsi(:,k), model_params%state_dim, MPI_DOUBLE_PRECISION, &
          particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
  END DO

  call NormalRandomNumbers2D(0.0D0,1.0D0,model_params%state_dim,pf%count,normaln)

  call Qhalf(pf%count,normaln,betan)
  Qkgain = 0.0_rk
  !$omp parallel do private(particle)
  DO k = 1,pf%count
     particle = pf%particles(k)
     call update_state(pf%psi(:,k),fpsi(:,k),Qkgain(:,k),betan(:,k))
!     if(checkscaling) call check_scaling(pf%psi(:,k),fpsi(:,k),betan(:,k))
  end DO
  !$omp end parallel do

end subroutine stochastic_model


subroutine check_scaling(x,fx,b,scales)

  use pf_control
  !use sizes
  use Qdata

  IMPLICIT NONE
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%state_dim), intent(in) :: x,fx,b
  REAL(KIND=rk), dimension(9), intent(inout) :: scales
  real(kind=rk) :: dnrm2,inc,be,rat,scal
  integer :: i
  integer, dimension(9) :: st,sp


  st = (/   1,  7009,140161,273313,406465,539617, 997128,1454639,1884535/)
  sp = (/7008,140160,273312,406464,539616,997127,1454638,1884534,2314430/)


  if(mod(pf%timestep,pf%time_obs) .gt. 0 .and. mod(pf%timestep,pf%time_obs)&
       & .le. 48) then
     do i = 1,9
        inc = dnrm2(sp(i)-st(i)+1,x(st(i):sp(i))-fx(st(i):sp(i)),1)
        be = dnrm2(sp(i)-st(i)+1,b(st(i):sp(i)),1)
        rat = be/inc
        print*,sqrt(sum( b(st(i):sp(i))**2 )),sqrt(sum( (x(st(i):sp(i))-fx(st(i):sp(i)))**2))
        scales(i) = rat
!        Qdiag(st(i):sp(i)) = Qdiag(st(i):sp(i))*scal     
     end do
  else
     do i = 1,9
        inc = dnrm2(sp(i)-st(i)+1,x(st(i):sp(i))-fx(st(i):sp(i)),1)
        be = dnrm2(sp(i)-st(i)+1,b(st(i):sp(i)),1)
        rat = be/inc
        print*,sqrt(sum( b(st(i):sp(i))**2 )),sqrt(sum( (x(st(i):sp(i))-fx(st(i):sp(i)))**2))
        scales(i) = rat
     end do
  end if
!  print*,'scales pf%timestep ',pf%timestep
!  do i = 1,9
     print*,scales
!  end do
  print*,'DIAG: ',Qdiag(st(:))

  !pstar 1:7008
  !a_u 7009:140160
  !a_v 140161:273312  
  !a_t 273313:406464
  !a_q 406465:539616
  !o_t 539617:997127
  !o_s 997128:1454638
  !o_u 1454639:1884534
  !o_v 1884535:2314430

end subroutine check_scaling
