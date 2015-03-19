!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-10-06 17:01:09 pbrowne>
!!!
!!!    Collection of subroutines to perform checks of user supplied routines
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
subroutine H_tests()
  !> @brief These are some tests to check that the observation
  !> operator is implemented correctly
  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%obs_dim) :: a,s
  real(kind=rk), dimension(model_params%obs_dim,10) :: aa,ss
  real(kind=rk), dimension(model_params%state_dim) :: j
  real(kind=rk), dimension(model_params%state_dim,10) :: jj
  real(kind=rk) :: dnrm2,r
  integer :: pfcount,i,k

  !let us save pf%count for later
  pfcount = pf%count


  write(6,*) 'TESTING H'
  write(6,*) 'First doing H^T then doing H'
  write(6,*) 'TESTING H WITH SINGLE RIGHT HAND SIDES'

  !FIRST TESTS: HHT
  !test one - zeros
  a = 0.0d0
  pf%count = 1
  call HT(model_params%obs_dim,pf%count,a,j,1)
  call H(model_params%obs_dim,pf%count,j,s,1)
  !s should be a
  r = dnrm2(model_params%obs_dim,s,1)
  write(6,'(A)',advance='no') 'Test 1: H H^T(0) ... '
  if(r .lt. 1.0d-15) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with r = ',r
  end if

  a = 1.0d0
  call HT(model_params%obs_dim,pf%count,a,j,1)
  call H(model_params%obs_dim,pf%count,j,s,1)
  !s should be a
  r = dnrm2(model_params%obs_dim,s-a,1)
  write(6,'(A)',advance='no') 'Test 2: H H^T(1) ... '
  if(r .lt. 1.0d-15) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with r = ',r
  end if
  
  a = -1.0d0
  call HT(model_params%obs_dim,pf%count,a,j,1)
  call H(model_params%obs_dim,pf%count,j,s,1)
  !s should be a
  r = dnrm2(model_params%obs_dim,s-a,1)
  write(6,'(A)',advance='no') 'Test 3: H H^T(-1) ... '
  if(r .lt. 1.0d-15) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with r = ',r
  end if


  call NormalRandomNumbers1D(0.0_rk,1.0_rk,model_params%obs_dim,a)
  call HT(model_params%obs_dim,pf%count,a,j,1)
  call H(model_params%obs_dim,pf%count,j,s,1)
  !s should be a                              
  r = dnrm2(model_params%obs_dim,s-a,1)
  write(6,'(A)',advance='no') 'Test 4: H H^T( N(0,1) ) ... '
  if(r .lt. 1.0d-15) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with r = ',r
  end if

  write(6,'(A)',advance='yes') 'TESTING H WITH MULTIPLE RIGHT HAND SIDES'
  do i = 1,10
     aa = 0.0_rk
     pf%count = i
     call HT(model_params%obs_dim,pf%count,aa(:,1:i),jj(:,1:i),i)
     call  H(model_params%obs_dim,pf%count,jj(:,1:i),ss(:,1:i),i)
     write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,' RHS: H H^T(0) ... '
     do k = 1,i
        !s should be a
        r = dnrm2(model_params%obs_dim,ss(:,k)-aa(:,k),1)
!        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,',',k,': H H^T(0) ... '
        if(r .lt. 1.0d-15) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        else
           !write(6,'(A)',advance='yes') ''
           write(6,'(A,es8.2)',advance='yes') 'failed with r = ',r
        end if
     end do
  end do

  
  do i = 1,10
     aa = 1.0_rk
     pf%count = i
     call HT(model_params%obs_dim,pf%count,aa(:,1:i),jj(:,1:i),i)
     call  H(model_params%obs_dim,pf%count,jj(:,1:i),ss(:,1:i),i)
     write(6,'(A,i2,A,i2,A)',advance='no') 'Test 6 ',i,' RHS: H H^T(1) ... '
     do k = 1,i
        !s should be a
        r = dnrm2(model_params%obs_dim,ss(:,k)-aa(:,k),1)
!        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,',',k,': H H^T(0) ... '
        if(r .lt. 1.0d-15) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        else
           write(6,'(A)',advance='yes') ''
           write(6,'(A,es23.17)',advance='yes') 'failed with r = ',r
        end if
     end do
  end do
  
  do i = 1,10
     aa = 0.0_rk
     pf%count = i
     call HT(model_params%obs_dim,pf%count,aa(:,1:i),jj(:,1:i),i)
     call  H(model_params%obs_dim,pf%count,jj(:,1:i),ss(:,1:i),i)
     write(6,'(A,i2,A,i2,A)',advance='no') 'Test 7 ',i,' RHS: H H^T(-1) ... '
     do k = 1,i
        !s should be a
        r = dnrm2(model_params%obs_dim,ss(:,k)-aa(:,k),1)
!        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,',',k,': H H^T(0) ... '
        if(r .lt. 1.0d-15) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        else
           write(6,'(A)',advance='yes') ''
           write(6,'(A,es23.17)',advance='yes') 'failed with r = ',r
        end if
     end do
  end do

  
  call NormalRandomNumbers2D(0.0_rk,1.0_rk,model_params%obs_dim,10,aa)
  do i = 1,10
     pf%count = i
     call HT(model_params%obs_dim,pf%count,aa(:,1:i),jj(:,1:i),i)
     call  H(model_params%obs_dim,pf%count,jj(:,1:i),ss(:,1:i),i)
     write(6,'(A,i2,A,i2,A)',advance='no') 'Test 8 ',i,' RHS: H H^T( N&
          &(0,1) ) ... '
     do k = 1,i
        !s should be a
        r = dnrm2(model_params%obs_dim,ss(:,k)-aa(:,k),1)
!        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,',',k,': H H^T(0) ... '
        if(r .lt. 1.0d-15) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        else
           write(6,'(A)',advance='yes') ''
           write(6,'(A,es8.2)',advance='yes') 'failed with r = ',r
        end if
     end do
  end do


  
  !put the count back to what it was before
  pf%count = pfcount
  

end subroutine H_tests


























subroutine R_tests()
  !> @brief These are some tests to check that the observation
  !> error covariance matrix is implemented correctly
  !use sizes
  use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%obs_dim) :: a,s,q,w,e
  real(kind=rk), dimension(model_params%obs_dim,10) :: aa,ss,qq,ww,ee
  real(kind=rk) :: dnrm2,rr
  real(kind=rk), parameter :: pass = 1.0d-15
  real(kind=rk), parameter :: warn = 1.0d-13
  integer :: pfcount,i,k,l

  !let us save pf%count for later
  pfcount = pf%count


  write(6,*) 'TESTING R'
  write(6,*) 'TESTING R WITH SINGLE RHS'


  !FIRST TESTS: HHT
  !test one - zeros
  a = 0.0d0
  pf%count = 1
  call R(model_params%obs_dim,1,a,s,1)
 
  !s should be 0
  rr = dnrm2(model_params%obs_dim,s,1)
  write(6,'(A)',advance='no') 'Test 1: R(0) ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  a = 0.0d0
  call Rhalf(model_params%obs_dim,1,a,s,1)
  !s should be a
  rr = dnrm2(model_params%obs_dim,s,1)
  write(6,'(A)',advance='no') 'Test 2: Rhalf(0) ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if
  
  a = 1.0d0
  s = 1.0d0
  call     R(model_params%obs_dim,1,a,q,1)
  call Rhalf(model_params%obs_dim,1,s,w,1)

  if(all(w .eq. 0.0d0)) then
     write(6,'(A)') 'SERIOUS ERROR: Rhalf*e = 0 i.e. Rhalf is the zero&
          & matrix'
     stop 'MEGA FAIL'
  end if

  call Rhalf(model_params%obs_dim,1,w,e,1)
  !q should be e
  rr = dnrm2(model_params%obs_dim,q-e,1)
  write(6,'(A)',advance='no') 'Test 3: [RhalfRhalf(1)]-[R(1)] ... '
  if(rr .lt. pass) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  a = -1.0d0
  s = -1.0d0
  call     R(model_params%obs_dim,1,a,q,1)
  call Rhalf(model_params%obs_dim,1,s,w,1)
  call Rhalf(model_params%obs_dim,1,w,e,1)
  !q should be e
  rr = dnrm2(model_params%obs_dim,q-e,1)
  write(6,'(A)',advance='no') 'Test 4: [RhalfRhalf(-1)]-[R(-1)] ... '
  if(rr .lt. pass) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  call NormalRandomNumbers1D(0.0d0,1.0d0,model_params%obs_dim,a)
  s = a
  call     R(model_params%obs_dim,1,a,q,1)
  call Rhalf(model_params%obs_dim,1,s,w,1)
  call Rhalf(model_params%obs_dim,1,w,e,1)
  !q should be e
  rr = dnrm2(model_params%obs_dim,q-e,1)
  write(6,'(A)',advance='no') 'Test 5: [RhalfRhalf( N(0,1) )]-[R( N(0,1) )] ... '
  if(rr .lt. pass) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if
 
  write(6,'(A)',advance='yes') 'TESTING R WITH MULTIPLE RIGHT HAND SIDES'

  
  do l = 6,9


  do i = 1,10
     if( l .eq. 6) then
        aa = 0.0_rk
        ss = 0.0_rk
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 6 ',i,' RHS: R(0)-&
             &Rhalf(Rhalf(0)) ... '
     elseif( l .eq. 7) then
        aa = 1.0_rk
        ss = 1.0_rk
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 7 ',i,' RHS: R(1)-&
             &Rhalf(Rhalf(1)) ... '        
     elseif( l .eq. 8) then
        aa = -1.0_rk
        ss = -1.0_rk
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 8 ',i,' RHS: R(-1)-Rha&
          &lf(Rhalf(-1)) ... '
     
     elseif( l .eq. 9) then
        call NormalRandomNumbers2D(0.0d0,1.0d0,model_params%obs_dim,i,aa(:,1:i))
        ss = aa
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 9 ',i,' RHS: R( N(&
             &0,1) )-&
             &Rhalf(Rhalf( N(0,1) )) ... '   
     end if
     call     R(model_params%obs_dim,i,aa(:,1:i),qq(:,1:i),1)
     call Rhalf(model_params%obs_dim,i,ss(:,1:i),ww(:,1:i),1)
     call Rhalf(model_params%obs_dim,i,ww(:,1:i),ee(:,1:i),1)
     
     do k = 1,i
        !qq should be ee
        rr = dnrm2(model_params%obs_dim,qq(:,k)-ee(:,k),1)
!        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,',',k,': H H^T(0) ... '
        if(rr .lt. pass) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        elseif(rr .lt. warn) then

           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='yes') 'passed: ',k,' warn ',rr
           elseif(k .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='no') 'passed: ',k,' warn ',rr
           elseif(k .eq. i) then
              write(6,'(A,i2,A,es8.2)',advance='yes') ' ',k,' warn ',rr
           else
              write(6,'(A,i2,A,es8.2)',advance='no') ' ',k,' warn ',rr
           end if
        else
           !write(6,'(A)',advance='yes') ''
           write(6,'(i2,A,es8.2)',advance='yes') k,' failed with rr = ',rr
        end if
     end do
  end do

  end do

  write(6,*) 'TESTING SOLVE R WITH SINGLE RHS'
  
  a = 0.0d0
  pf%count = 1
  call solve_r(model_params%obs_dim,pf%count,a,s,1)
 
  !s should be 0
  rr = dnrm2(model_params%obs_dim,s,1)
  write(6,'(A)',advance='no') 'Test 10: R^(-1)(0) ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  do l = 11,13
     if( l .eq. 11) then
        a = 1.0d0
        write(6,'(A)',advance='no') 'Test 11: R^(-1)[R(1)] ... '
     elseif (l .eq. 12) then
        a = -1.0d0
        write(6,'(A)',advance='no') 'Test 12: R^(-1)[R(-1)] ... '
     elseif (l .eq. 13) then
        call NormalRandomNumbers1D(0.0d0,1.0d0,model_params%obs_dim,a)
        write(6,'(A)',advance='no') 'Test 13: R^(-1)[R( N(0,1) )] ... '
     end if


     call     R(model_params%obs_dim,1,a,q,1)
     call solve_r(model_params%obs_dim,pf%count,q,w,1)
     !w should be a
     rr = dnrm2(model_params%obs_dim,w-a,1)

     
     if(rr .lt. pass) then
        !passed the test                                   
        write(6,'(A)',advance='yes') 'passed'
     elseif(rr .lt. warn) then
        write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
     else
        write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
     end if

  end do

  write(6,*) 'TESTING SOLVE R WITH MULTIPLE RHS'
  
  do l = 14,16

  do i = 1,10
     if( l .eq. 14) then
        aa = 1.0d0
        write(6,'(A)',advance='no') 'Test 14: R^(-1)[R(1)] ... '
     elseif (l .eq. 15) then
        aa = -1.0d0
        write(6,'(A)',advance='no') 'Test 15: R^(-1)[R(-1)] ... '
     elseif (l .eq. 16) then
        call NormalRandomNumbers2D(0.0d0,1.0d0,model_params%obs_dim,10,aa)
        write(6,'(A)',advance='no') 'Test 16: R^(-1)[R( N(0,1) )] ... '
     end if

     pf%count = i
     call     R(model_params%obs_dim,i,aa(:,1:i),qq(:,1:i),1)
     call solve_r(model_params%obs_dim,pf%count,qq(:,1:i),ww(:,1:i),1)
     !w should be a
     do k = 1,i
        rr = dnrm2(model_params%obs_dim,ww(:,k)-aa(:,k),1)

        if(rr .lt. pass) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        elseif(rr .lt. warn) then

           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='yes') 'passed: ',k,' warn ',rr
           elseif(k .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='no') 'passed: ',k,' warn ',rr
           elseif(k .eq. i) then
              write(6,'(A,i2,A,es8.2)',advance='yes') ' ',k,' warn ',rr
           else
              write(6,'(A,i2,A,es8.2)',advance='no') ' ',k,' warn ',rr
           end if
        else
           !write(6,'(A)',advance='yes') ''
           write(6,'(i2,A,es8.2)',advance='yes') k,' failed with rr = ',rr
        end if


     end do

  end do

  end do


  write(6,*) 'TESTING SOLVE Rhalf WITH SINGLE RHS'
  
  a = 0.0d0
  pf%count = 1
  call solve_rhalf(model_params%obs_dim,pf%count,a,s,1)
 
  !s should be 0
  rr = dnrm2(model_params%obs_dim,s,1)
  write(6,'(A)',advance='no') 'Test 17: R^(-1/2)(0) ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  do l = 18,20
     if( l .eq. 18) then
        a = 1.0d0
        write(6,'(A)',advance='no') 'Test 18: R^{-1/2}[R^{1/2}(1)] ... '
     elseif (l .eq. 19) then
        a = -1.0d0
        write(6,'(A)',advance='no') 'Test 19: R^{-1/2}[R^{1/2}(-1)] ... '
     elseif (l .eq. 20) then
        call NormalRandomNumbers1D(0.0d0,1.0d0,model_params%obs_dim,a)
        write(6,'(A)',advance='no') 'Test 20: R^{-1/2}[R^{1/2}( N(0,1) )] ... '
     end if


     call     Rhalf(model_params%obs_dim,1,a,q,1)
     call solve_rhalf(model_params%obs_dim,pf%count,q,w,1)
     !w should be a
     rr = dnrm2(model_params%obs_dim,w-a,1)

     
     if(rr .lt. pass) then
        !passed the test                                   
        write(6,'(A)',advance='yes') 'passed'
     elseif(rr .lt. warn) then
        write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
     else
        write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
     end if

  end do

  write(6,*) 'TESTING SOLVE Rhalf WITH MULTIPLE RHS'
  
  do l = 21,23

  do i = 1,10
     if( l .eq. 21) then
        aa = 1.0d0
        write(6,'(A)',advance='no') 'Test 21: R^{-1/2}[R^{1/2}(1)] ... '
     elseif (l .eq. 22) then
        aa = -1.0d0
        write(6,'(A)',advance='no') 'Test 22: R^{-1/2}[R^{1/2}(-1)] ... '
     elseif (l .eq. 23) then
        call NormalRandomNumbers2D(0.0d0,1.0d0,model_params%obs_dim,10,aa)
        write(6,'(A)',advance='no') 'Test 23: R^{-1/2}[R^{1/2}( N(0,1) )] ... '
     end if

     pf%count = i
     call     Rhalf(model_params%obs_dim,i,aa(:,1:i),qq(:,1:i),1)
     call solve_rhalf(model_params%obs_dim,pf%count,qq(:,1:i),ww(:,1:i),1)
     !w should be a
     do k = 1,i
        rr = dnrm2(model_params%obs_dim,ww(:,k)-aa(:,k),1)

        if(rr .lt. pass) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        elseif(rr .lt. warn) then

           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='yes') 'passed: ',k,' warn ',rr
           elseif(k .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='no') 'passed: ',k,' warn ',rr
           elseif(k .eq. i) then
              write(6,'(A,i2,A,es8.2)',advance='yes') ' ',k,' warn ',rr
           else
              write(6,'(A,i2,A,es8.2)',advance='no') ' ',k,' warn ',rr
           end if
        else
           !write(6,'(A)',advance='yes') ''
           write(6,'(i2,A,es8.2)',advance='yes') k,' failed with rr = ',rr
        end if


     end do

  end do

  end do





  pf%count=pfcount

end subroutine R_tests



















subroutine Q_tests()
  !> @brief These are some tests to check that the model error
  !> covariance matrix is implemented correctly
  !use sizes
  use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%state_dim) :: a,s,qqq,w,e
  real(kind=rk), dimension(model_params%state_dim,10) :: aa,ss,qq,ww,ee
  real(kind=rk) :: dnrm2,rr
  real(kind=rk), parameter :: pass = 1.0d-15
  real(kind=rk), parameter :: warn = 1.0d-13
  integer :: pfcount,i,k,l

  !let us save pf%count for later
  pfcount = pf%count


  write(6,*) 'TESTING Q'
  write(6,*) 'TESTING Q WITH SINGLE RHS'


  !FIRST TESTS: HHT
  !test one - zeros
  a = 0.0d0
  pf%count = 1
  call Q(1,a,s,1)
 
  !s should be 0
  rr = dnrm2(model_params%state_dim,s,1)
  write(6,'(A)',advance='no') 'Test 1: Q(0) ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  a = 0.0d0
  call Qhalf(1,a,s,1)
  !s should be a
  rr = dnrm2(model_params%state_dim,s,1)
  write(6,'(A)',advance='no') 'Test 2: Qhalf(0) ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if
  
  a = 1.0d0
  s = 1.0d0
  call     Q(1,a,qqq,1)
  call Qhalf(1,s,w,1)

  if(all(w .eq. 0.0d0)) then
     write(6,'(A)') 'SERIOUS ERROR: Qhalf*e = 0 i.e. Qhalf is the zero&
          & matrix'
     stop 'MEGA FAIL'
  end if

  call Qhalf(1,w,e,1)
  !q should be e
  rr = dnrm2(model_params%state_dim,qqq-e,1)
  write(6,'(A)',advance='no') 'Test 3: [QhalfQhalf(1)]-[Q(1)] ... '
  if(rr .lt. pass) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  a = -1.0d0
  s = -1.0d0
  call     Q(1,a,qqq,1)
  call Qhalf(1,s,w,1)
  call Qhalf(1,w,e,1)
  !q should be e
  rr = dnrm2(model_params%state_dim,qqq-e,1)
  write(6,'(A)',advance='no') 'Test 4: [QhalfQhalf(-1)]-[R(-1)] ... '
  if(rr .lt. pass) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if

  call NormalRandomNumbers1D(0.0d0,1.0d0,model_params%state_dim,a)
  s = a
  call     Q(1,a,qqq,1)
  call Qhalf(1,s,w,1)
  call Qhalf(1,w,e,1)
  !q should be e
  rr = dnrm2(model_params%state_dim,qqq-e,1)
  write(6,'(A)',advance='no') 'Test 5: [QhalfQhalf( N(0,1) )]-[Q( N(0,1) )] ... '
  if(rr .lt. pass) then
     !passed the test                                   
     write(6,'(A)',advance='yes') 'passed'
  elseif(rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if
 
  write(6,'(A)',advance='yes') 'TESTING Q WITH MULTIPLE RIGHT HAND SIDES'

  
  do l = 6,8


  do i = 1,10
     if( l .eq. 6) then
        aa = 1.0_rk
        ss = 1.0_rk
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 6 ',i,' RHS: Q(1)-&
             &Qhalf(Qhalf(1)) ... '        
     elseif( l .eq. 7) then
        aa = -1.0_rk
        ss = -1.0_rk
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 7 ',i,' RHS: Q(-1)-Qha&
          &lf(Qhalf(-1)) ... '
     
     elseif( l .eq. 8) then
        call NormalRandomNumbers2D(0.0d0,1.0d0,model_params%state_dim,i,aa(:,1:i))
        ss = aa
        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 8 ',i,' RHS: Q( N(&
             &0,1) )-Qhalf(Qhalf( N(0,1) )) ... '   
     end if
     call     Q(i,aa(:,1:i),qq(:,1:i),1)
     call Qhalf(i,ss(:,1:i),ww(:,1:i),1)
     call Qhalf(i,ww(:,1:i),ee(:,1:i),1)
     
     do k = 1,i
        !qq should be ee
        rr = dnrm2(model_params%state_dim,qq(:,k)-ee(:,k),1)
!        write(6,'(A,i2,A,i2,A)',advance='no') 'Test 5 ',i,',',k,': H H^T(0) ... '
        if(rr .lt. pass) then
           !passed the test
           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2)',advance='yes') 'passed: ',k
           elseif(k .eq. 1) then
              write(6,'(A,i2)',advance='no') 'passed: ',k
           elseif(k .eq. i) then
              write(6,'(A,i2)',advance='yes') ' ',k
           else
              write(6,'(A,i2)',advance='no') ' ',k
           end if
        elseif(rr .lt. warn) then

           if(k .eq. 1 .and. i .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='yes') 'passed: ',k,' warn ',rr
           elseif(k .eq. 1) then
              write(6,'(A,i2,A,es8.2)',advance='no') 'passed: ',k,' warn ',rr
           elseif(k .eq. i) then
              write(6,'(A,i2,A,es8.2)',advance='yes') ' ',k,' warn ',rr
           else
              write(6,'(A,i2,A,es8.2)',advance='no') ' ',k,' warn ',rr
           end if
        else
           !write(6,'(A)',advance='yes') ''
           write(6,'(i2,A,es8.2)',advance='yes') k,' failed with rr = ',rr
        end if
     end do
  end do

  end do

  pf%count = pfcount

end subroutine Q_tests


























subroutine HQHTR_tests()
  !> @brief These are some tests to check that the linear solve
  !> operator is implemented correctly
  !>
  !>
  !> This should check the operation \f$(HQH^T+R)^{-1}\f$ is working
  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%obs_dim) :: a,s,qq,hqha,hqha_r
  real(kind=rk), dimension(model_params%state_dim) :: qha,ha
  real(kind=rk) :: dnrm2,rr
  real(kind=rk), parameter :: pass = 1.0d-15
  real(kind=rk), parameter :: warn = 1.0d-13

  integer :: pfcount,l

  !let us save pf%count for later
  pfcount = pf%count


write(6,*) 'TESTING (HQH^T+R)^(-1)'

  !FIRST TESTS: HHT
  !test one - zeros
  a = 0.0d0
  pf%count = 1
  call solve_hqht_plus_r(model_params%obs_dim,a,s,1)
  !s should be a
  rr = dnrm2(model_params%obs_dim,s,1)
  write(6,'(A)',advance='no') 'Test 1: (HQH^T+R)^(-1)[0] ... '
  if(rr .lt. pass) then
     !passed the test
     write(6,'(A)',advance='yes') 'passed'
  elseif( rr .lt. warn) then
     write(6,'(A,es23.17)',advance='yes') 'passed with warning rr = ',rr   
  else
     write(6,'(A,es23.17)',advance='yes') 'failed with rr = ',rr
  end if


  do l = 2,10
     if( l .eq. 2) then
        a = 1.0d0
        write(6,'(A)',advance='no') 'Test 2: (HQH^T+R)^(-1)[(HQH^T+R)(1)] ... '
     elseif (l .eq. 3) then
        a = -1.0d0
        write(6,'(A)',advance='no') 'Test 3: (HQH^T+R)^(-1)[(HQH^T+R)(1)] ... '
     elseif (l .ge. 4) then
        call NormalRandomNumbers1D(0.0d0,1.0d0,model_params%obs_dim,a)
        write(6,'(A,i2,A)',advance='no') 'Test ',l,': (HQH^T+R)^(-1)[(&
             &HQH^T+R)( N(0,1) )] ... '
     end if


     call     R(model_params%obs_dim,1,a,qq,1)
     call   HT(model_params%obs_dim,1,a,ha,1)
     call    Q(1,ha,qha)
     call    H(model_params%obs_dim,1,qha,hqha,1)
     hqha_r = hqha+qq
     call solve_hqht_plus_r(model_params%obs_dim,hqha_r,s,1)
     !w should be a
     rr = dnrm2(model_params%obs_dim,s-a,1)

     
     if(rr .lt. pass) then
        !passed the test                                   
        write(6,'(A)',advance='yes') 'passed'
     elseif(rr .lt. warn) then
        write(6,'(A,es8.2)',advance='yes') 'warning rr = ',rr
     else
        write(6,'(A,es8.2)',advance='yes') 'failed  rr = ',rr
     end if

  end do
  




end subroutine HQHTR_tests


