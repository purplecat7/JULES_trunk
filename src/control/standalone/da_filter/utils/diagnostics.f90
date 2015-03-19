!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-26 12:19:36 pbrowne>
!!!
!!!    Subroutine to give output diagnositics such as rank histograms
!!!    and trajectories
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

!> Subroutine to give output diagnositics such as rank histograms  
!!    and trajectories                                             
subroutine diagnostics
  use pf_control
  !use sizes
  use empire_mod !comms
  use histogram_data
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  real(kind=rk), dimension(rhl_n,pf%nens) :: HHpsi
  real(kind=rk), dimension(pf%nens) :: bin_marker
  real(kind=rk), dimension(rhl_n) :: y
  integer :: particle,time,i,j,mpi_err,lb,length
  logical :: placed
  character(32) :: filename
  integer, dimension(rhn_n,pf%nens+1) :: reduced_talagrand
  !include 'mpif.h'
  
  if(.not. pf%gen_data) then
     if(pf%use_talagrand) then
        
!        call H(pf%psi,Hpsi)
        inquire(iolength=length) pf%psi(1,1)

!        pf%timestep = -72
        do particle = 1,pf%count
           write(filename,'(A,i6.6,A,i5.5)') 'hist/timestep',((pf%timestep)/pf&
                &%time_bwn_obs) ,'particle',pf%particles(particle)
           open(12,file=filename,action='write',status='replace',form='unforma&
                &tted',access='direct',recl=length)
!           call H(pf%psi(:,particle),Hpsi)
           do i = 1,rhl_n
              write(12,rec=i) pf%psi(rank_hist_list(i),particle)
           end do
           close(12)
        end do

!        pf%timestep = 0
        if(pf%timestep .eq. pf%time_obs*pf%time_bwn_obs) then
           pf%talagrand = 0
           call mpi_barrier(pf_mpi_comm,mpi_err)
           
           do time = 1,pf%time_obs
              !              print*,'time = ',time
              if(mod(time,npfs) .eq. pfrank) then
                 !cock cock cock adjusted the below to make it sensible
                 pf%timestep = time*pf%time_bwn_obs
!                 print*,'pfrank = ',pfrank,'picking up truth at ',pf%timestep

                 write(filename,'(A,i6.6,A,i5.5)') 'hist/timestep',((pf%timestep)/pf&
                      &%time_bwn_obs) ,'truth'
                 !print*,'filename = ',filename
                 open(12,file=filename,action='read',status='old',form='unforma&
                      &tted',access='direct',recl=length)
                 !           call H(pf%psi(:,particle),Hpsi)
                 do i = 1,rhl_n
                    !print*,'pfrank = ',pfrank,' i = ',i
                    read(12,rec=i) y(i)
                    !print*,'pfrank = ',pfrank,' i = ',i,' y(i) = ',y(i)
                 end do
                 close(12)

                 !call get_truth(y)
                 !              print*,'got obs data at timestep ',pf%timestep
                 do particle = 1,pf%nens
                    write(filename,'(A,i6.6,A,i5.5)') 'hist/timestep',pf&
                         &%timestep/pf%time_bwn_obs,'&
                         &particle',particle
                    open(13, file=filename,action='read',status='old',form='un&
                         &for&
                         &matted',access='direct',recl=length)
                    do i = 1,rhl_n
                       read(13,rec=i) HHpsi(i,particle)
                    end do
                    close(13)
                 end do
                 !             print*,'read HHpsi'


                 do j = 1,rhn_n !for each histogram we want to make
                    if( j .eq. 1) then
                       lb = 1
                    else
                       lb = sum(rank_hist_nums(1:j-1)) + 1
                    end if
                 do i = lb,sum(rank_hist_nums(1:j))

                    do particle = 1,pf%nens
                       bin_marker(particle) = HHpsi(i,particle)
                    end do

                    call quicksort_d(bin_marker,pf%nens)
                    
                  
                    !                    print*,'quicksorted'
                    placed = .false.
                    do particle  = 1,pf%nens
                       if(y(i) .lt. bin_marker(particle)) then
                          pf%talagrand(j,particle) = pf&
                               &%talagrand(j,particle)&
                               & + 1
                          placed = .true.
                          exit
                       end if
                    end do
                    !                    print*,'did we place?'

                    if(.not. placed) then
                       if(y(i) .ge. bin_marker(pf%nens)) then
                          pf%talagrand(j,pf%nens+1) = pf&
                               &%talagrand(j,pf%nens+1) + 1
                       else
                          stop 'There was an error in the calculation of the p&
                               &lacement &
                               &in the rank histogram. Bums.'
                       end if
                    end if

                 end do
                 end do

              end if !end of mpi splitting by timestep
           end do !end of the timestep
           !           print*,'end of all the timesteps'

           !now let us reduce the information to the master processor:
           call mpi_reduce(pf%talagrand,reduced_talagrand,rhn_n*(pf%nens+1)&
                &,mpi_integer,mpi_sum,0,pf_mpi_comm,mpi_err)

           !           print*,'some reduction just happened'

           if(pfrank .eq. 0) then
              do i = 1,rhn_n
                 write(filename,'(A,i0)') 'histogram_',i
                 open(17,file=filename,action='write',status='replace')
                 do j = 1,pf%nens+1
                    write(17,'(i8.8)') reduced_talagrand(i,j)
                 end do
                 close(17)
              end do
              !now output the image
           end if
           pf%timestep = pf%time_obs*pf%time_bwn_obs
        end if !end of if we are the last step in the pf

     end if

     if(pf%use_rmse) then
        

     end if

  else
     if(pf%use_talagrand) then
        !     print*,'in diagnostics at timestep ',pf%timestep
        inquire(iolength=length) pf%psi(1,1)
        
        !        pf%timestep = -72
        do particle = 1,pf%count
           write(filename,'(A,i6.6,A)') 'hist/timestep',((pf%timestep)/pf&
                &%time_bwn_obs) ,'truth'
           open(12,file=filename,action='write',status='replace',form='unforma&
                &tted',access='direct',recl=length)
           !           call H(pf%psi(:,particle),Hpsi)
           do i = 1,rhl_n
              write(12,rec=i) pf%psi(rank_hist_list(i),particle)
           end do
           close(12)
        end do
     end if
  end if
end subroutine diagnostics

!> subroutine to output trajectories
subroutine trajectories
  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  integer :: particle,i,j
  character(28) :: filename
  integer, parameter :: n=10
  integer, dimension(n) :: trajvar

  stop 'Trajectories have not be specified for this specific model'

  !             ap   au    av    at      aq    ot       ot     os     ou      ov
  trajvar = (/2721,10099,145949,278007,410009,547011,558240,998895,1480094,1900706/)

  do i = 1,pf%count
     particle = pf%particles(i)
     
     do j = 1,n
        if(pf%gen_data) then
           write(filename,'(A,i7.7)') 'traj/truth_var',trajvar(j)
        else
           write(filename,'(A,i7.7,A,i5.5)') 'traj/var',trajvar(j),'particle',particle
        end if
        if(pf%timestep .eq. 0) then
!           print*,'i = ',i,' particle = ',particle,' trajvar(j) = '&
!                &,trajvar(j),' filename = ',filename
           open(41,file=filename,action='write',status='replace')
        else
           open(41,file=filename,action='write',status='old',position='append')
        end if
        write(41,'(es22.15)') pf%psi(trajvar(j),i)
        close(41)
        
     end do
  end do
end subroutine trajectories
