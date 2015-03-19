!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-29 15:12:27 pbrowne>
!!!
!!!    Collection of subroutines to deal with i/o
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

!> Subroutine to read observation from a file
!! \n
!! Uses pf%timestep to determine which observation to read
!> @param[out] y The observation
subroutine get_observation_data(y)

  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(model_params%obs_dim), intent(out) :: y
  integer :: obs_number,ios
  character(14) :: filename

  obs_number = ((pf%timestep-1)/pf%time_bwn_obs) + 1

  write(filename,'(A,i6.6)') 'obs_num_',obs_number

  open(67,file=filename,iostat=ios,action='read',status='old',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Check it exists. I need a lie down.'
     stop
  end if
  read(67) y
  close(67)
end subroutine get_observation_data

!> Subroutine to save observation to a file              
!! \n              
!! Uses pf%timestep to determine which observation to save   
!> @param[in] y The observation
subroutine save_observation_data(y)

  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(model_params%obs_dim), intent(in) :: y
  integer :: obs_number,ios
  character(14) :: filename

  obs_number = ((pf%timestep-1)/pf%time_bwn_obs) + 1

  write(filename,'(A,i6.6)') 'obs_num_',obs_number

  open(67,file=filename,iostat=ios,action='write',status='replace',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) '1 PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
     stop
  end if
  write(67) y
  close(67)

  if(pf%human_readable) then
     open(64,file='pf_data',action='write',position='append')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_data'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
     write(64,*) pf%timestep,y
     close(64)
  end if
end subroutine save_observation_data

!> Subroutine to save truth to a file              
!! \n                 
!> @param[in] x The state vector
subroutine save_truth(x)

  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(model_params%state_dim), intent(in) :: x
  integer :: ios
  if(pf%timestep .eq. 0) then
     print*,'opening pf_truth'
     open(62,file='pf_truth',iostat=ios,action='write',status='replace')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_truth'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
  end if
  write(62,*) x
  call flush(62)
  if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
     close(62)
     print*,'closing pf_truth'
  end if
end subroutine save_truth

!>subroutine to ouput data from the filter
subroutine output_from_pf
  use pf_control
  !use sizes
  use empire_mod !comms
  implicit none
  real(kind=kind(1.0D0)), dimension(model_params%state_dim) :: mean
  integer :: ios,particle
  character(9) :: filename
  if(pf%timestep .eq. 0) then
     write(filename,'(A,i2.2)') 'pf_out_',pfrank
     open(68,file=filename,iostat=ios,action='write',status='replace')
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_out'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
  end if
!  print*,'output: ',pf%timestep,pf%weight
!  write(68,*) pf%timestep,pf%particles,pf%weight(:)
  write(68,'(i6.6,A)',advance='no') pf%timestep,' '
  do ios = 1,pf%count-1
      PRINT *, 'particle array: ', pf%particles
     write(68,'(i6.6,A,e21.15,A)',advance='no') pf%particles(ios),' ',pf&
          &%weight(pf%particles(ios)),' '
  end do
  write(68,'(i6.6,A,e21.15)',advance='yes') pf%particles(pf%count),' ',pf&
          &%weight(pf%particles(pf%count))
  call flush(68)
  if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) close(68)

  if(pf%use_mean) then
     if(pf%timestep .eq. 0) then
        open(61,file='pf_mean',iostat=ios,action='write',status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_mean'
           write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
           stop
        end if
     end if
     mean = 0.0D0
     do particle = 1,pf%nens
        mean(:) = mean(:) + pf%psi(:,particle)*exp(-pf%weight(particle))
     end do
     write(61,*) mean(:)
     call flush(61)
     if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) close(61)
  end if

  if(pf%use_weak) then
     if(pf%timestep .eq. 0) then
        open(77,file='pf_weak',iostat=ios,action='write',status='replace&
             &')
     else
        open(77,file='pf_weak',iostat=ios,action='write',status&
             &='old',position='append')
     end if
     if(ios .ne. 0)  then
        write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file pf_&
             &weak'
        write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
        stop
     end if
     
     do particle = 1,pf%nens
        write(77,*) pf%psi(:,particle),exp(-pf%weight(particle))
     end do
     close(77)
     
     if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
        open(77,file='gnuplot_weak.cfg',iostat=ios,action='write'&
             &,status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'Cannot open file gnuplot_weak'
           write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
           stop
        end if
        write(77,'(A)') 'reset'
        write(77,'(A)') 'set terminal wxt size 600,700'
        write(77,'(A)') 'set macros'
        write(77,'(A)') 'set key off'
        write(77,'(A,f0.5,A)') 'set yrange [0:',4.0D0/pf%nens,']'
        write(77,'(A)') 'Str(k)=sprintf("%d",k)'
        write(77,'(A)') 'S=0'
        write(77,'(A,i0)') 'E=S+',pf%nens-1
        write(77,'(A,i0)') 'n=',pf%nens*pf%time_bwn_obs*pf%time_obs
        write(77,'(A)') 'load "animateweak.lorenz"'
        close(77)
        open(77,file='animateweak.lorenz',iostat=ios,action='write'&
             &,status='replace')
        if(ios .ne. 0)  then
           write(*,*) 'Cannot open file animateweak.lorenz'
           write(*,*) 'Very strange that I couldnt open it. Im going to stop now.'
           stop
        end if
        write(77,'(A,i0)') 'set multiplot layout 3,1 title "Lorenz 63 usi&
             &ng equal weight particle filter. Timestep = ".S/',pf%nens
        write(77,'(A)') 'myplotoptions=''every ::''.Str(S).''::''.Str(&
             &E)'
        write(77,'(A)') 'set xrange [-20:20]'
        write(77,'(A)') 'set title "x unobserved"'
        write(77,'(A)') 'set xlabel "x"'
        write(77,'(A)') 'plot "pf_weak" using 1:4 @myplotoptions w imp&
             &ulses'
        write(77,'(A)') 'set xrange [-30:30]'
        write(77,'(A)') 'set title "y unobserved"'
        write(77,'(A)') 'set xlabel "y"'
        write(77,'(A)') 'plot "pf_weak" using 2:4 @myplotoptions w imp&
             &ulses'

        write(77,'(A)') 'set xrange [0:50]'
        write(77,'(A)') 'set title "z observed"'
        write(77,'(A)') 'set xlabel "z"'
        write(77,'(A)') 'plot "pf_weak" using 3:4 @myplotoptions w imp&
             &ulses'
        
        write(77,'(A)') 'unset multiplot'
        
        write(77,'(A,i0)') 'S=S+',pf%nens
        write(77,'(A,i0)') 'E=E+',pf%nens
        write(77,'(A)') 'if (E < n) reread'
        close(77)
     end if
  end if







end subroutine output_from_pf

!> subroutine to save the state vector to a named file
!! as an unformatted fortran file
subroutine save_state(state,filename)
  !use sizes
    use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%state_dim), intent(in) :: state !< the
  !!state vector
  character(14), intent(in) :: filename !< the name of the file to
  !!save the state vector in
  integer :: ios

  open(16,file=filename,iostat=ios,action='write',status='replace'&
       &,form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) '2 PARTICLE FILTER DATA ERROR!!!!! Cannot open file '&
          &,filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop&
          & now.'
     stop
  end if
  write(16) state
  close(16)
end subroutine save_state


!> subroutine to write the state vector to a named file
!! as an unformatted fortran file
subroutine get_state(state,filename)
  !use sizes
        use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  real(kind=rk), dimension(model_params%state_dim), intent(out) :: state !< the
  !!state vector
  character(14), intent(in) :: filename!< the name of the file to 
  !!write the state vector in
  integer :: ios

  open(16,file=filename,iostat=ios,action='read',status='old',form='un&
       &formatted')
  if(ios .ne. 0)  then
     write(*,*) '3 PARTICLE FILTER DATA ERROR!!!!! Cannot open file '&
          &,filename
     write(*,*) 'Very strange that I couldnt open it. Im going to stop&
          & now.'
     stop
  end if
  read(16) state
  close(16)
end subroutine get_state
