!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-10-06 12:26:19 pbrowne>
!!!
!!!    Ensemble transform Kalman filter
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
! $RCSfile: etkf_analysis.f90,v $
! $Revision: 1.6 $
! $Date: 2013/12/12 03:58:54 $

!> subroutine to perform the ensemble transform Kalman filter as part
!! of L-ETKF
!>
!> 
subroutine letkf_analysis
use empire_mod !comms
use pf_control
!use sizes
implicit none
integer, parameter :: rk = kind(1.0d0)

!!> number of ensemble members
!integer, intent(in) :: N
!!> current size of state dimension
!integer, intent(in) :: stateDimension
!!> total number of observations
!integer, intent(in) :: obsDim
!!> Forecast ensemble on entry, analysis ensemble on exit
!real(kind=rk), dimension(stateDimension,N), intent(inout) :: x
real(kind=rk), dimension(model_params%state_dim,pf%count) :: x_loc

!> The observation
real(kind=rk), dimension(model_params%obs_dim) :: y
!real(kind=rk), dimension(obsDim) :: y
!!> Inflation parameter; forecast perturbations will be scaled by 1+rho
!real(kind=rk), intent(in) :: rho
!!> Localisation length scale
!real(kind=rk), intent(in) :: len
!!> the timestep
!integer, intent(in) :: t


! Local variables for the SVD
integer :: r
real(kind=rk), dimension(:,:), allocatable :: V
real(kind=rk), dimension(:), allocatable :: S
real(kind=rk), dimension(pf%nens,pf%nens) :: UT
integer :: LWORK,INFO
real(kind=rk), dimension(:), allocatable :: WORK

! Miscellaneous local variables
real(kind=rk), dimension(model_params%state_dim) :: mean_x
real(kind=rk), dimension(:), allocatable :: mean_xa
real(kind=rk), dimension(model_params%state_dim,pf%count) :: Xp_loc
!real(kind=rk), dimension(stateDimension,N) :: Xp
real(kind=rk), dimension(:,:), allocatable :: Xp
real(kind=rk), dimension(:,:), allocatable :: Xa
real(kind=rk), dimension(model_params%obs_dim) :: mean_yf,d,dd
real(kind=rk), dimension(model_params%obs_dim,pf%nens) :: yf,Ysf
real(kind=rk), dimension(model_params%obs_dim,pf%count) :: yf_loc
integer :: i,j,number_gridpoints

!variables for localisation
real(kind=rk) :: dist
real(kind=rk), parameter :: maxscal=1.0d12
real(kind=rk),dimension(model_params%obs_dim) :: scal
logical, dimension(model_params%obs_dim) :: yes
integer :: red_obsdim
real(kind=rk), allocatable, dimension(:,:) :: Ysf_red
real(kind=rk), allocatable, dimension(:) :: dd_red
integer :: stateDim !generally 1 for the letkf

!variables for mpi
integer :: mpi_err
integer, dimension(npfs) :: start_var,stop_var
!!$character(6) :: filename
!include 'mpif.h'

call get_observation_data(y)
!!$print*,'obs y = ',y

! Split forecast ensemble into mean and perturbation matrix, inflating
! if necessary
! mean_x will only be the sum of state vectors on this mpi process
mean_x = sum(pf%psi,dim=2)

!send mean_x to all processes and add up to get global sum
call mpi_allreduce(MPI_IN_PLACE,mean_x,model_params%state_dim,MPI_DOUBLE_PRECISION&
     &,MPI_SUM,pf_mpi_comm,mpi_err)

!now divide by the total number of ensemble members to make it a true
!mean
mean_x = mean_x/real(pf%nens,rk)

!!$open(1,file='before',action='write',status='replace')
!!$write(1,*) mean_x
!!$close(1)

! compute the ensemble perturbation matrix for those ensemble members
! stored on this local mpi process
do i = 1,pf%count
   Xp_loc(:,i) = pf%psi(:,i) - mean_x
end do

! inflate the ensemble perturbation matrix
Xp_loc = (1.0_rk + pf%rho) * Xp_loc
! store the local state vectors back in x_loc
do i = 1,pf%count
   x_loc(:,i) = mean_x + Xp_loc(:,i)
end do

! make the local ensemble perturbation matrix the correct scale
Xp_loc = Xp_loc/sqrt(real(pf%nens-1,rk))

! Calculate forecast observations, split into mean and ensemble
! perturbation matrix, scale perturbations by inverse square root of
! observation covariance

! first apply observation operator only to local state vectors
! on this mpi process
call H(model_params%obs_dim,pf%count,x_loc,yf_loc,pf%timestep)


! as yf_loc should be much smaller than x_loc, send this to mpi processes
! need to send round all yf_loc and store in yf on all processes
call mpi_allgatherv(yf_loc,pf%count*model_params%obs_dim,MPI_DOUBLE_PRECISION,yf&
     &,gblcount*model_params%obs_dim,gbldisp*model_params%obs_dim,MPI_DOUBLE_PRECISION&
     &,pf_mpi_comm,mpi_err)

!!$print*,'allgatherv 1 = ',mpi_err

! compute the mean of yf
mean_yf = sum(yf,dim=2)/real(pf%nens,rk)
do i = 1,pf%nens
   yf(:,i) = yf(:,i) - mean_yf
end do

!!$print*,'mean_yf = ',mean_yf
!!$print*,'mean_x = ',mean_x
! now make yf the forecast perturbation matrix
yf = yf/sqrt(real(pf%nens-1,rk))
!call solve_rhalf(N,obsDim,yf,Ysf)
! scale yf to become ysf
call solve_rhalf(model_params%obs_dim,pf%nens,yf,Ysf,pf%timestep)


! now let us compute which state variables will be analysed on each
! MPI process:
do i = 1,npfs
   start_var(i) = (i-1)*ceiling( real(model_params%state_dim,rk)/real(npfs,rk) )+1
   stop_var(i) = min( i*ceiling(real(model_params%state_dim,rk)/real(npfs,rk)) ,model_params%state_dim)
end do

!allocate space for Xp and Xa now that we know how many grid points we consider
allocate(Xp(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
allocate(Xa(stop_var(pfrank+1)-start_var(pfrank+1)+1,pf%nens))
allocate(mean_xa(stop_var(pfrank+1)-start_var(pfrank+1)+1))
mean_xa = mean_x(start_var(pfrank+1):stop_var(pfrank+1))

!now we have to get Xp filled with all the values from each process
!!$print*,'AHOY gbldisp = ',gbldisp

do i = 1,npfs
   number_gridpoints = stop_var(i)-start_var(i)+1
!!$   print*,'i = ',i
!!$   print*,'number_gridpoitns = ',number_gridpoints
!!$   print*,'size(xp_loc) = ',size(Xp_loc(start_var(i):stop_var(i),:))
!!$   print*,'sendcount = ',pf%count*number_gridpoints
!!$   print*,'sendtype = ',MPI_DOUBLE_PRECISION
!!$   print*,'recvbuf = ',size(Xp)
!!$   print*,'recvcounts = ',gblcount*number_gridpoints
!!$   print*,'displs = ',(gbldisp)*number_gridpoints
!!$   print*,'recvtype = ',MPI_DOUBLE_PRECISION
!!$   print*,'root = ',i-1
!!$   print*,'comm = ',pf_mpi_comm
!!$   print*,'mpi_err = ',mpi_err

   call mpi_gatherv(Xp_loc(start_var(i):stop_var(i),:),& !sendbuf
        pf%count*number_gridpoints,&                     !sendcount
        MPI_DOUBLE_PRECISION,&                           !sendtype
        Xp,&                                             !recvbuf
        gblcount*number_gridpoints,&                     !recvcounts
        gbldisp*number_gridpoints,&                      !displs
        MPI_DOUBLE_PRECISION,&                           !recvtype
        i-1,&                                            !root
        pf_mpi_comm,&                                    !comm
        mpi_err)                                         !ierror
end do

!!$print*,'gatherv = ',mpi_err


d = y - mean_yf
!!$print*,'d = ',d
!call solve_rhalf(1,obsDim,d,dd)
call solve_rhalf(model_params%obs_dim,1,d,dd,pf%timestep)

!!$print*,'dd = ',dd


!I THINK I CAN STEP IN AT THIS POINT!!!
!this is for serial processing
number_gridpoints = stop_var(pfrank+1)-start_var(pfrank+1)+1
!$OMP PARALLEL DO &
!$OMP& PRIVATE(stateDim,i,dist,scal,yes), &
!$OMP& PRIVATE(Ysf_red,red_obsDim,r,V,S,LWORK,WORK,INFO), &
!$OMP& PRIVATE(dd_red,UT,d)
do j = 1,number_gridpoints
   stateDim = 1
   yes = .false.
   !let us process the observations here:
   do i = 1,model_params%obs_dim
      ! get the distance between the current state variable
      ! j+start_var(pfrank+1)-1
      ! and the observation i
      ! store it as distance
      call dist_st_ob(j+start_var(pfrank+1)-1,i,dist,pf%timestep)
      ! compute the scaling factor based in this distance and the
      ! length scale
      scal(i) = sqrt(exp((dist**2)/(2.0_rk*pf%len**2)))
      ! determine if the scaling is such that we shall compute with
      ! this observation
      if(scal(i) .lt. maxscal) yes(i) = .true.
   end do
!!$   print*,'scal = ',scal
!!$   print*,'yes = ',yes
   
   ! count the total number of observations we shall consider for this
   ! state variable
   red_obsdim = count(yes)
   print*,'j = ',j,' red_obsdim = ',red_obsdim,scal
   PRINT*,'var = ',j+start_var(pfrank+1)-1


   ! if there are no observations in range, treat this as a special case
   if(red_obsdim .gt. 0) then

   allocate(Ysf_red(red_obsdim,pf%nens))
   !multiply by the distance matrix
   !this line only works for diagonal R matrix...
   ! reduce the forecast ensemble to only the observations in range
   ! and scale by the distance function
   do i = 1,pf%nens
      Ysf_red(:,i) = pack(Ysf(:,i),yes)/pack(scal,yes)
   end do
!!$   print*,'Ysf_red ',Ysf_red 

   ! Compute the SVD
   r = min(red_obsDim,pf%nens)
   allocate(V(red_obsDim,r))
   allocate(S(r))
   LWORK = 2*max( 3*r+max(red_obsDim,pf%nens), 5*r )
   !LWORK=400
   allocate(WORK(LWORK))
!!$   print*,'LWORK = ',LWORK
   call dgesvd('S','A',red_obsDim,pf%nens,Ysf_red,red_obsDim,S,V,red_obsDim,UT,pf%nens,WORK,LWORK,INFO)
   if(INFO .ne. 0) then
      print*,'SVD failed with INFO = ',INFO
      print*,'FYI WORK(1) = ',WORK(1)
      stop
   end if
   deallocate(WORK)

!!$   PRINT*,'size(UT)  = ',size(UT)
!!$   print*,'svd finished'
!!$
!!$   print*,'S = ',S
!!$   print*,'V = ',V
!!$   print*,'UT = ',UT
   ! Compute product of forecast ensemble perturbation matrix and U.  We
   ! store the result in xa, which we here call Xa to distinguish this
   ! secondary use of the variable.
   ! pick out the correct row of Xa and Xp here:
   call dgemm('N','T',stateDim,pf%nens,pf%nens,1.0d0,Xp(j,:),stateDim,UT,pf%nens,0.0d0,Xa(j,:),stateDim)
   !!$print*,'Xa(j,:) = ',Xa(j,:)


   ! Build up analysis ensemble mean (to be stored in mean_x); done now
   ! because we shall be reusing X shortly
   
   allocate(dd_red(red_obsDim))
   dd_red = pack(dd,yes)/pack(scal,yes)
   
   !!$print*,'dd_red = ',dd_red


   ! Only the first r elements of d are in use from here on
   call dgemv('T',red_obsDim,r,1.0d0,V,red_obsDim,dd_red,1,0.0d0,d,1)

   !!$print*,'d(1:r) = ',d(1:r)
   d(1:r) = S * d(1:r) / (1.0_rk + S**2)
   !!$print*,'d(1:r) = ',d(1:r)
   ! Only the first r columns of Xa are used in the following
   call dgemv('N',stateDim,r,1.0d0,Xa(j,:),stateDim,d,1,1.0d0,mean_xa(j),1)
   
   !!$print*,'mean_xa(j) = ',mean_xa(j)
   ! Build up analysis ensemble perturbation matrix
   do i = 1,r
      Xa(j,i) = Xa(j,i) / sqrt(1.0_rk + S(i)**2)
   end do
   call dgemm('N','N',stateDim,pf%nens,pf%nens,1.0d0,Xa(j,:),stateDim,UT,pf%nens,0.0d0,Xp(j,:),stateDim)
   
   ! Put ensemble mean and perturbation matrix back together
   xa(j,:) = sqrt(real(pf%nens-1,rk))*Xp(j,:)
   do i = 1,pf%nens
      xa(j,i) = mean_xa(j) + xa(j,i)
   end do

   !put this back into the full state vector
   !wait it is already
   deallocate(Ysf_red)
   deallocate(V)
   deallocate(S)
   deallocate(dd_red)
!   if(j ==2) stop

   else !if there are no observations near, the analysis is
        ! just the forecast
      do i = 1,pf%nens
         xa(j,i) =  mean_xa(j) + xp(j,i)
      end do
      print*,'xp ',xp(j,:)
      print*,'xa ',xa(j,:)
      print*,'mean_xa ',mean_xa(j)
   end if
end do
!$OMP END PARALLEL DO


!!$print*,'loop finished'

! now we must scatter xa back into pf%psi
do i = 1,npfs
   
   number_gridpoints = stop_var(i)-start_var(i)+1
!!$   print*,'i = ',i
!!$   print*,'number_gridpoitns = ',number_gridpoints
!!$   print*,'size(Xa) = ',size(Xa)
!!$   print*,'sendcounts = ',gblcount*number_gridpoints
!!$   print*,'displs = ',gbldisp*number_gridpoints
!!$   print*,'sendtype = ',MPI_DOUBLE_PRECISION
!!$   print*,'recvbuf = ',size(pf%psi(start_var(i):stop_var(i),:))
!!$   print*,'recvcount = ',pf%count*number_gridpoints
!!$   print*,'recvtype = ',MPI_DOUBLE_PRECISION
!!$   print*,'root = ',i-1
!!$   print*,'comm = ',pf_mpi_comm
!!$   print*,'mpi_err = ',mpi_err

   call mpi_scatterv(Xa,&                                !sendbuf  
        gblcount*number_gridpoints,&                     !sendcounts  
        gbldisp*number_gridpoints,&                      !displs
        MPI_DOUBLE_PRECISION,&                           !sendtype    
        pf%psi(start_var(i):stop_var(i),:),&             !recvbuf  
        pf%count*number_gridpoints,&                     !recvcount  
        MPI_DOUBLE_PRECISION,&                           !recvtype    
        i-1,&                                            !root     
        pf_mpi_comm,&                                    !comm     
        mpi_err)                                         !ierror 
end do

!!$write(filename,'(A,i1)') 'after',pfrank
!!$open(1,file=filename,action='write',status='replace')
!!$write(1,*) mean_xa
!!$close(1)

deallocate(mean_xa)
deallocate(Xp)
deallocate(Xa)

end subroutine letkf_analysis
