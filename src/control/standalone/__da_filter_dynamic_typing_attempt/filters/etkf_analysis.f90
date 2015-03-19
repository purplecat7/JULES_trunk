!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-25 11:46:38 pbrowne>
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

!> subroutine to perform the ensemble transform Kalman filter
!>
!> 
subroutine etkf_analysis(num_hor,num_ver,this_hor,this_ver,boundary,x,N,stateDim,obsDim,rho)

implicit none
integer, parameter :: rk = kind(1.0d0)

integer, intent(in) :: N,stateDim,obsDim
! Forecast ensemble on entry, analysis ensemble on exit
real(kind=rk), dimension(stateDim,N), intent(inout) :: x
! The observation
!real(kind=rk), dimension(obsDim), intent(in) :: y
real(kind=rk), dimension(obsDim) :: y
! Inflation parameter; forecast perturbations will be scaled by 1+rho
! if rho is present
real(kind=rk), intent(in) :: rho
! controls for localisation
integer, intent(in) :: num_hor,num_ver,this_hor,this_ver,boundary

! Local variables for the SVD
integer :: r
real(kind=rk), dimension(:,:), allocatable :: V
real(kind=rk), dimension(:), allocatable :: S
real(kind=rk), dimension(N,N) :: UT
integer :: LWORK,INFO
real(kind=rk), dimension(:), allocatable :: WORK

! Miscellaneous local variables
real(kind=rk), dimension(stateDim) :: mean_x
real(kind=rk), dimension(stateDim,N) :: Xp
real(kind=rk), dimension(obsDim) :: mean_yf,d,dd
real(kind=rk), dimension(obsDim,N) :: yf,Ysf
integer :: i

call get_local_observation_data(num_hor,num_ver,this_hor,this_ver,boundary,obsDim,y)

! Split forecast ensemble into mean and perturbation matrix, inflating
! if necessary
mean_x = sum(x,dim=2)/real(N,rk)
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
call H_local(num_hor,num_ver,this_hor,this_ver,boundary,N,stateDim,x&
     &,obsDim,yf)
mean_yf = sum(yf,dim=2)/real(N,rk)
do i = 1,N
   yf(:,i) = yf(:,i) - mean_yf
end do
yf = yf/sqrt(real(N-1,rk))
!call solve_rhalf(N,obsDim,yf,Ysf)
call solve_rhalf_local(num_hor,num_ver,this_hor,this_ver,boundary,N&
     &,obsDim,yf,Ysf)

! Compute the SVD
r = min(obsDim,N)
allocate(V(obsDim,r))
allocate(S(r))
LWORK = 2*max( 3*r+max(obsDim,N), 5*r )
allocate(WORK(LWORK))
call dgesvd('S','A',obsDim,N,Ysf,obsDim,S,V,obsDim,UT,N,WORK,LWORK,INFO)
if(INFO .ne. 0) then
   print*,'SVD failed with INFO = ',INFO
   print*,'FYI WORK(1) = ',WORK(1)
   stop
end if
deallocate(WORK)

! Compute product of forecast ensemble perturbation matrix and U.  We
! store the result in x, which we here call X to distinguish this
! secondary use of the variable.
call dgemm('N','T',stateDim,N,N,1.0d0,Xp,stateDim,UT,N,0.0d0,X,stateDim)

! Build up analysis ensemble mean (to be stored in mean_x); done now
! because we shall be reusing X shortly
d = y - mean_yf
!call solve_rhalf(1,obsDim,d,dd)
call solve_rhalf_local(num_hor,num_ver,this_hor,this_ver,boundary,1&
     &,obsDim,d,dd)
! Only the first r elements of d are in use from here on
call dgemv('T',obsDim,r,1.0d0,V,obsDim,dd,1,0.0d0,d,1)
d(1:r) = S * d(1:r) / (1.0_rk + S**2)
! Only the first r columns of X are used in the following
call dgemv('N',stateDim,r,1.0d0,X,stateDim,d,1,1.0d0,mean_x,1)

! Build up analysis ensemble perturbation matrix
do i = 1,r
   X(:,i) = X(:,i) / sqrt(1.0_rk + S(i)**2)
end do
call dgemm('N','N',stateDim,N,N,1.0d0,X,stateDim,UT,N,0.0d0,Xp,stateDim)

! Put ensemble mean and perturbation matrix back together
x = sqrt(real(N-1,rk))*Xp
do i = 1,N
   x(:,i) = mean_x + x(:,i)
end do

end subroutine etkf_analysis
