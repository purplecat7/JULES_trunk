!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-29 18:07:40 pbrowne>
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
subroutine letkf_analysis(x,N,stateDimension,obsDim,y,rho,len,t)

implicit none
integer, parameter :: rk = kind(1.0d0)

!> number of ensemble members
integer, intent(in) :: N
!> current size of state dimension
integer, intent(in) :: stateDimension
!> total number of observations
integer, intent(in) :: obsDim
!> Forecast ensemble on entry, analysis ensemble on exit
real(kind=rk), dimension(stateDimension,N), intent(inout) :: x
!> The observation
real(kind=rk), dimension(obsDim), intent(in) :: y
!real(kind=rk), dimension(obsDim) :: y
!> Inflation parameter; forecast perturbations will be scaled by 1+rho
real(kind=rk), intent(in) :: rho
!> Localisation length scale
real(kind=rk), intent(in) :: len
!> the timestep
integer, intent(in) :: t


! Local variables for the SVD
integer :: r
real(kind=rk), dimension(:,:), allocatable :: V
real(kind=rk), dimension(:), allocatable :: S
real(kind=rk), dimension(N,N) :: UT
integer :: LWORK,INFO
real(kind=rk), dimension(:), allocatable :: WORK

! Miscellaneous local variables
real(kind=rk), dimension(stateDimension) :: mean_x
real(kind=rk), dimension(stateDimension,N) :: Xp
real(kind=rk), dimension(obsDim) :: mean_yf,d,dd
real(kind=rk), dimension(obsDim,N) :: yf,Ysf
integer :: i,j,number_gridpoints

!variables for localisation
real(kind=rk) :: dist
real(kind=rk), parameter :: maxscal=1.0d12
real(kind=rk),dimension(obsDim) :: scal
logical, dimension(obsDim) :: yes
integer :: red_obsdim
real(kind=rk), allocatable, dimension(:,:) :: Ysf_red
real(kind=rk), allocatable, dimension(:) :: dd_red
integer :: stateDim !generally 1 for the letkf

print*,'obs y = ',y

! Split forecast ensemble into mean and perturbation matrix, inflating
! if necessary
mean_x = sum(x,dim=2)/real(N,rk)
open(1,file='before',action='write',status='replace')
write(1,*) mean_x
close(1)

do i = 1,N
   Xp(:,i) = x(:,i) - mean_x
end do
!if (present(rho)) then
Xp = (1.0_rk + rho) * Xp
do i = 1,N
   x(:,i) = mean_x + Xp(:,i)
end do
!end if
Xp = Xp/sqrt(real(N-1,rk))

! Calculate forecast observations, split into mean and ensemble
! perturbation matrix, scale perturbations by inverse square root of
! observation covariance
!call H(N,stateDim,obsDim,x,yf)
call H(obsDim,N,x,yf,t)
mean_yf = sum(yf,dim=2)/real(N,rk)
do i = 1,N
   yf(:,i) = yf(:,i) - mean_yf
end do

print*,'mean_yf = ',mean_yf
print*,'mean_x = ',mean_x

yf = yf/sqrt(real(N-1,rk))
!call solve_rhalf(N,obsDim,yf,Ysf)
call solve_rhalf(obsDim,N,yf,Ysf,t)

!I THINK I CAN STEP IN AT THIS POINT!!!
!this is for serial processing
number_gridpoints = stateDimension
!!!!$OMP PARALLEL DO PRIVATE(stateDim,i,dist,scal,yes,Ysf_red,red_obsDim,r,V,S,LWORK,WORK,dd_red,UT,dd,d,number_gridpoints)
do j = 1,number_gridpoints
   stateDim = 1
   yes = .false.
   !let us process the observations here:
   do i = 1,obsDim
      call dist_st_ob(j,i,dist,t)
      scal(i) = exp((dist**2)/(2.0_rk*len**2))
      if(scal(i) .lt. maxscal) yes(i) = .true.
   end do
   print*,'scal = ',scal
   print*,'yes = ',yes
   
   red_obsdim = count(yes)
   print*,'j = ',j,' red_obsdim = ',red_obsdim,scal

   allocate(Ysf_red(red_obsdim,N))
   !multiply by the distance matrix
   !this line only works for diagonal R matrix...
   !move the sqrt to the calculation of scal probably
   do i = 1,N
      Ysf_red(:,i) = pack(Ysf(:,i),yes)/sqrt(pack(scal,yes))
   end do
   print*,'Ysf_red ',Ysf_red 

   ! Compute the SVD
   r = min(red_obsDim,N)
   allocate(V(red_obsDim,r))
   allocate(S(r))
   LWORK = 2*max( 3*r+max(red_obsDim,N), 5*r )
   !LWORK=400
   allocate(WORK(LWORK))
   print*,LWORK
   call dgesvd('S','A',red_obsDim,N,Ysf_red,red_obsDim,S,V,red_obsDim,UT,N,WORK,LWORK,INFO)
   if(INFO .ne. 0) then
      print*,'SVD failed with INFO = ',INFO
      print*,'FYI WORK(1) = ',WORK(1)
      stop
   end if
   deallocate(WORK)
   

   print*,'S = ',S
   print*,'V = ',V
   print*,'UT = ',UT
   ! Compute product of forecast ensemble perturbation matrix and U.  We
   ! store the result in x, which we here call X to distinguish this
   ! secondary use of the variable.
   ! pick out the correct row of X and Xp here:
   call dgemm('N','T',stateDim,N,N,1.0d0,Xp(j,:),stateDim,UT,N,0.0d0,X(j,:),stateDim)
   print*,'X(j,:) = ',X(j,:)
   ! Build up analysis ensemble mean (to be stored in mean_x); done now
   ! because we shall be reusing X shortly
   d = y - mean_yf
   print*,'d = ',d
   !call solve_rhalf(1,obsDim,d,dd)
   call solve_rhalf(obsDim,1,d,dd,t)
   
   print*,'dd = ',dd
   allocate(dd_red(red_obsDim))
   dd_red = pack(dd,yes)/sqrt(pack(scal,yes))
   
   print*,'dd_red = ',dd_red


   ! Only the first r elements of d are in use from here on
   call dgemv('T',red_obsDim,r,1.0d0,V,red_obsDim,dd_red,1,0.0d0,d,1)

   print*,'d(1:r) = ',d(1:r)
   d(1:r) = S * d(1:r) / (1.0_rk + S**2)
   print*,'d(1:r) = ',d(1:r)
   ! Only the first r columns of X are used in the following
   call dgemv('N',stateDim,r,1.0d0,X(j,:),stateDim,d,1,1.0d0,mean_x(j),1)
   
   print*,'mean_x(j) = ',mean_x(j)
   ! Build up analysis ensemble perturbation matrix
   do i = 1,r
      X(j,i) = X(j,i) / sqrt(1.0_rk + S(i)**2)
   end do
   call dgemm('N','N',stateDim,N,N,1.0d0,X(j,:),stateDim,UT,N,0.0d0,Xp(j,:),stateDim)
   
   ! Put ensemble mean and perturbation matrix back together
   x(j,:) = sqrt(real(N-1,rk))*Xp(j,:)
   do i = 1,N
      x(j,i) = mean_x(j) + x(j,i)
   end do

   !put this back into the full state vector
   !wait it is already
   deallocate(Ysf_red)
   deallocate(V)
   deallocate(S)
   deallocate(dd_red)
!   if(j ==2) stop
end do
!!!!$OMP END PARALLEL DO


open(1,file='after',action='write',status='replace')
write(1,*) mean_x
close(1)
end subroutine letkf_analysis
