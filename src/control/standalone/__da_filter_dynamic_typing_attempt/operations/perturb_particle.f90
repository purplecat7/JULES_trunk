!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-10-01 13:45:06 pbrowne>
!!!
!!!    Collection of routines to perturb and update states
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

!> Subroutine to perturb state vector with normal random vector
!! drawn from \f$\mathcal{N}(0,Q)\f$
subroutine perturb_particle(x)
!use sizes
use comms
use pf_control
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(model_params%state_dim), intent(inout) :: x
real(kind=rk), dimension(model_params%state_dim) :: rdom,y,kgain
character(14) :: filename

if(pf%init .eq. 'P') then
   call NormalRandomNumbers1D(0.0D0,1.0D0,model_params%state_dim,rdom)
   call Qhalf(1,rdom,y)
   rdom = 0.0_rk
   kgain = 0.0_rk
   call update_state(rdom,x,kgain,y)
   x = rdom
elseif(pf%init .eq. 'N') then

!!$   x = (/4.82160531136,-1.77565265252,2.10461993106,3.46854171351&
!!$        &,7.19135143229,7.12754206422,-0.796082655036,2.45172556459&
!!$        &,-0.531641795086,-0.19366799976,-4.17678736568,6.94259660694&
!!$        &,4.67019291167,-1.63118486372,0.830259694341,1.74204790657&
!!$        &,5.63332497432,1.38456126234,-1.67858158117,2.15882451035&
!!$        &,8.24797269894,-1.18246777768,-3.58939011025,3.26473996747&
!!$        &,8.83474281013,3.34455161902,2.31141994205,-2.84763733923&
!!$        &,1.87634007581,3.3191730343,7.36150967249,3.99624410686 &
!!$        &,-7.63041979612,1.43836147006,-0.955124293143,5.38338216121&
!!$        &,2.92386476709,1.14393995662,3.17942949736,6.00540694151/)

   call NormalRandomNumbers1D(0.0D0,1.0D0,model_params%state_dim,rdom)
   kgain = 0.0_rk
   call update_state(y,x,kgain,rdom)
   x = y

elseif(pf%init .eq. 'R') then
   !get ensemble member from the restart folder                  
   write(filename,'(A,i2.2,A)') 'rstrt/',pfrank,'.state'
   call get_state(x,filename)
elseif(pf%init .eq. 'S') then
   !get ensemble member from the start folder                  
   if(pf%gen_data) then
      write(filename,'(A,i2.2,A)') 'start/',32,'.state'
   else
      write(filename,'(A,i2.2,A)') 'start/',pfrank,'.state'
      print*,'pf #',pfrank,' starting from ',filename
   end if
   call get_state(x,filename)
else
   print*,'ERROR: incorrect pf%init selected in perturb_particle'
   stop
end if


end subroutine perturb_particle

!> Subroutine to update the state
!!
!! This can be changed for the specific model
!! if it needs to be
!!
!! @param[in] fpsi deterministic model update \f$f(x^{n-1})\f$
!! @param[in] kgain nudging term
!! @param[inout] betan Stochastic term
!! @param[out] state The updated state vector 

subroutine update_state(state,fpsi,kgain,betan)
use pfcontrol!sizes
integer, parameter :: rk=kind(1.0D0)
real(kind=rk), dimension(model_params%state_dim), intent(out) :: state
real(kind=rk), dimension(model_params%state_dim), intent(in) :: fpsi,kgain
real(kind=rk), dimension(model_params%state_dim), intent(inout) :: betan
!real(kind=rk), dimension(model_params%state_dim) :: diff
!integer :: k,start
real(kind=rk) :: dnrm2
!logical, parameter :: debug=.false.
logical, parameter :: norms=.false.


!!atmosphere humidity = state_vector(406465:539616)
!do k = 406465,539616
!   ! fpsi(i)+kgain(i)+betan(i) has to lie between 0 and 1   
!   betan(k) = max(-fpsi(k)-kgain(k),min(betan(k),1.0D0-fpsi(k)-kgain(k)))
!end do


!!ocean salinity = state_vector(997128:1454638)
!do k = 997128,1454638
!   ! fpsi(i)+kgain(i)+betan(i) has to lie between 0 and 1
!   betan(k) = max(-fpsi(k)-kgain(k),min(betan(k),1.0D0-fpsi(k)-kgain(k)))
!end do


!do the addition
state = fpsi+kgain+betan
!if(debug) print*,'|fpsi|=',dnrm2(model_params%state_dim,fpsi,1),' |kgain|= ',dnrm2(model_params%state_dim,kgain&
!     &,1),' |betan| = ',dnrm2(model_params%state_dim,betan,1)
if(norms) print*,' |kgain|= ',dnrm2(model_params%state_dim,kgain,1),' |betan| = ',dnrm2(model_params%state_dim,betan,1)


!now fix the atmospheric polar values:
! PSTAR     U     V     THETA    Q

!index = 0
!do k=1,levels
!  do j=1,columns
!    do i=1,rows
!       index = index + 1
!       data_out(index)=data(i,j,k) == data(longitude,latitude,level)
!    enddo
!  enddo
!enddo

!at the poles, all the longitude values are equal
!first lets do Pstar
!state(1:a_nxn) = state(1)
!state((a_nyn-1)*a_nxn+1:a_nyn*a_nxn) = state((a_nyn-1)*a_nxn+1)

!start = a_nxn*a_nyn   !for U
!do k = 1,a_levels
   !j = 1
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test U k:',k,&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : &
!        start+(k-1)*a_nyn*a_nxn + a_nxn))-&
!        minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : &
!        start+(k-1)*a_nyn*a_nxn + a_nxn))
   !j = a_nyn
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test U k:',k,maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))
!end do
!start = a_nxn*a_nyn + 1*a_levels*a_nyn*a_nxn !for V
!do k = 1,a_levels
   !j = 1
   !        
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test V k:',k,&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & a_nxn))&
!        -minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn))
   !j = a_nyn
   !        
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test V k:',k,maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))
!end do
!start = a_nxn*a_nyn + 2*a_levels*a_nyn*a_nxn !for Theta
!do k = 1,a_levels
   !j = 1
   !        
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test Theta k:',k,&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn +&
!        & a_nxn))&
!        -minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn))
   !j = a_nyn
   !        
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test Theta k:',k,maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn))
!end do
!start = a_nxn*a_nyn + 3*a_levels*a_nyn*a_nxn !for Q
!do k = 1,a_levels
   !j = 1
   !        
!   state( start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + a_nxn )  =&
!        & state( start+(k-1)*a_nyn*a_nxn + 1 )
!if(debug)   print*,'fpsi test Q:',k,':',maxval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k-1)&
!        &*a_nyn*a_nxn + a_nxn))-minval(fpsi(start+(k-1)*a_nyn*a_nxn + 1 : start+(k&
!        &-1)*a_nyn*a_nxn + a_nxn ))
   !j = a_nyn
   !        
!   state( start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + 1 : start+(k-1)*a_nyn&
!        &*a_nxn +&
!        & (a_nyn-1)*a_nxn + a_nxn) = state(start+(k-1)*a_nyn*a_nxn +&
!        & (a_nyn-1)&
!        &*a_nxn + 1)
!if(debug)   print*,'fpsi test Q:',k,':',&
!        maxval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)&
!        &*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + a_nxn))&
!        -minval(fpsi(start+(k-1)*a_nyn*a_nxn + (a_nyn-1)&
!        &*a_nxn + 1 : start+(k-1)*a_nyn*a_nxn + (a_nyn-1)*a_nxn + a_nxn))
!end do

!if(debug)print*,'PSTAR: min fpsi= ',minval(fpsi(1:7008))
!if(debug)print*,'PSTAR: max fpsi= ',maxval(fpsi(1:7008))
!if(debug)print*,'PSTAR: min state= ',minval(state(1:7008))
!if(debug)print*,'PSTAR: max state= ',maxval(state(1:7008))




!if(debug)diff = abs(state-fpsi)
!if(debug)print*,'max perturbation is ',maxval(diff),' at ',maxloc(diff),' of '&
!     &,state(maxloc(diff))
!if(debug)diff = diff/(max(abs(fpsi),1.0d-8))
!if(debug)print*,'max prop perturb is ',maxval(diff),' at ',maxloc(diff),' of '&
!     &,fpsi(maxloc(diff))




end subroutine update_state
