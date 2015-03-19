!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-18 10:09:54 pbrowne>
!!!
!!!    {one line to give the program's name and a brief idea of what it does.}
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
subroutine H_local(num_hor,num_ver,this_hor,this_ver,boundary,nrhs,stateDim,x&
     &,obsDim,y)
  use pf_control
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: num_hor,num_ver,this_hor,this_ver,boundary,nrhs&
       &,stateDim,obsDim
  real(kind=rk), dimension(stateDim,nrhs), intent(in) :: x
  real(kind=rk), dimension(obsDim,nrhs), intent(out) :: y
  integer :: nx,ny,nnx,nny,start_x,end_x,start_y,end_y
  integer :: i,j,k,ii,jj

  nx = 256/num_hor
  ny = 256/num_ver

  start_x = max(nx*(this_hor-1)+1-boundary,1)
  end_x = min(nx*(this_hor) + boundary,256)
  nnx = end_x-start_x+1

  start_y = max(ny*(this_ver-1)+1-boundary,1)
  end_y = min(ny*(this_ver)+boundary,256)
  nny = end_y-start_y+1


!  print*,'h local nx = ',nx,'ny = ',ny,'nnx = ',nnx,'nny = ',nny
!  print*,'h local start_x = ',start_x,'end x = ',end_x
!  print*,'h local start_y = ',start_y,'end_y = ',end_y
  k = 0
  do i = start_x,end_x
     ii = i-start_x+1
     do j = start_y,end_y
        jj = j - start_y + 1
        if(mod(i,pf%redObs) .eq. 1 .and. mod(j,pf%redObs) .eq. 1) then
           k = k + 1
           y(k,:) = x((ii-1)*nny+jj,:)
        end if
     end do
  end do

end subroutine H_local


subroutine solve_Rhalf_local(num_hor,num_ver,this_hor,this_ver,boundary,nrhs&
     &,obsDim,y,v)
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, intent(in) :: num_hor,num_ver,this_hor,this_ver,boundary,nrhs,obsDim

  real(kind=rk), dimension(obsDim,nrhs), intent(in) :: y
  real(kind=rk), dimension(obsDim,nrhs), intent(out) :: v

  v = 0.0_rk
  call daxpy(obsDim*nrhs,1.0_rk/0.05_rk,y,1,v,1)

end subroutine solve_Rhalf_local

subroutine get_local_observation_data(num_hor,num_ver,this_hor,this_ver,boundary,obsDim,y)
  use pf_control
  !use sizes
  implicit none
  integer, parameter :: rk = kind(1.0D0)
  integer, intent(in) :: num_hor,num_ver,this_hor,this_ver,boundary,obsDim
  real(kind=rk), dimension(obsDim), intent(out) :: y
  real(kind=rk), dimension(model_params%obs_dim) :: yfull
  integer :: obs_number,ios,i,j,k
  integer :: nx,ny,nnx,nny,start_x,end_x,start_y,end_y
  character(14) :: filename

  obs_number = ((pf%timestep-1)/pf%time_bwn_obs) + 1

  write(filename,'(A,i6.6)') 'obs_num_',obs_number

  open(67,file=filename,iostat=ios,action='read',status='old',form='unformatted')
  if(ios .ne. 0)  then
     write(*,*) 'PARTICLE FILTER DATA ERROR!!!!! Cannot open file ',filename
     write(*,*) 'Check it exists. I need a lie down.'
     stop
  end if
  read(67) yfull
  close(67)

  nx = (256/pf%redObs)/num_hor
  ny = (256/pf%redObs)/num_ver
  
  start_x = max(nx*(this_hor-1)+1-(boundary/pf%redObs),1)
  end_x = min(nx*(this_hor) + (boundary/pf%redObs),256/pf%redObs)
!  print*,'start_x = ',start_x,'end_x = ',end_x
!  print*,'this_nor = ',this_hor
!  print*,'nx*(this_hor-1)+1 = ',nx*(this_hor-1)+1
!  print*,'boundary/pf%redobs = ',(boundary/pf%redObs)
  nnx = end_x-start_x+1

  start_y = max(ny*(this_ver-1)+1-(boundary/pf%redObs),1)
  end_y = min(ny*(this_ver)+(boundary/pf%redObs),256/pf%redObs)
  nny = end_y-start_y+1
!  print*,'get_local_obs: nx = ',nx,'ny=',ny,'nnx = ',nnx,'nny = ',nny

  k = 0
  do i = start_x,end_x
     do j = start_y,end_y
        k = k + 1
        y(k) = yfull((i-1)*nny+j)
     end do
  end do
  if(k .ne. obsDim) then
     print*,'WOAH!! k .ne. obsdim'
     print*,k, obsdim
     print*,this_hor,this_ver,nnx,nny
     print*,start_x,end_x
     print*,start_y,end_y
     stop
  end if

end subroutine get_local_observation_data

subroutine localise_enkf(enkf_analysis)
!use sizes
use pf_control
use empire_mod !comms
implicit none
!include 'mpif.h'
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: enkf_analysis
real(kind=rk), dimension(model_params%state_dim,pf%count) :: x_analysis
real(kind=rk), dimension(:,:), allocatable :: x_local
integer :: mpi_err,particle,tag
integer :: num_hor,num_ver,boundary,stateDim,Obsdim,obsx,obsy
integer :: nx,ny,nnx,nny,start_x,end_x,start_y,end_y,k,ii,jj,i,j
integer, dimension(mpi_status_size) :: mpi_status


do k =1,pf%count
   particle = pf%particles(k)
   tag = 1
   call mpi_send(pf%psi(:,k),model_params%state_dim,MPI_DOUBLE_PRECISION&
        &,particle-1,tag,CPL_MPI_COMM,mpi_err)
end do

DO k = 1,pf%count
   particle = pf%particles(k)
   tag = 1
   CALL MPI_RECV(pf%psi(:,k), model_params%state_dim, MPI_DOUBLE_PRECISION, &
        particle-1, tag, CPL_MPI_COMM,mpi_status, mpi_err)
END DO




if(mod(pf%timestep,pf%time_bwn_obs) .eq. 0) then
num_hor = pf%enkfx
num_ver = pf%enkfy
boundary = pf%enkfb



do i = 1, num_hor
   do j = 1,num_ver

      nx = 256/num_hor
      ny = 256/num_ver
!      print*,'nx = ',nx
!      print*,'ny = ',ny
      
      start_x = max(nx*(i-1)+1-boundary,1)
      end_x = min(nx*(i) + boundary,256)
      nnx = end_x-start_x+1
!      print*,'nnx = ',nnx
      start_y = max(ny*(j-1)+1-boundary,1)
      end_y = min(ny*(j)+boundary,256)
      nny = end_y-start_y+1
!      print*,'nny = ',nny
      stateDim = nnx*nny
!      print*,'stateDim = ',stateDim
      obsx = 0
      do ii = start_x,end_x
         if(mod(ii,pf%redObs) .eq. 1) obsx = obsx + 1
      end do
      obsy = 0
      do jj = start_y,end_y
         if(mod(jj,pf%redObs) .eq. 1) obsy = obsy + 1
      end do
      obsDim = obsx*obsy
!      print*,'obsDim = ',obsdim

      allocate(x_local(stateDim,pf%count))
      k = 0
      do ii = start_x,end_x
         do jj = start_y,end_y
            k = k+1
            x_local(k,:) = pf%psi((ii-1)*256 + jj,:)
         end do
      end do

      
      if(enkf_analysis .eq. 3) then
!         print*,'before : '
!         print*,x_local(1:10,1)
         call etkf_analysis(num_hor,num_ver,i,j,boundary,x_local,pf%count&
              &,stateDim,obsDim,pf%Qscale)
!         print*,'after : '
!         print*,x_local(1:10,1)
      elseif(enkf_analysis .eq. 4) then
         call eakf_analysis(num_hor,num_ver,i,j,boundary,x_local,pf%count&
              &,stateDim,obsDim,pf%Qscale)
      else
         do k = 1,500
            print*,'I AM A FISH'
         END do
         stop 'ENKF SELECTED INCORRECTLY: 2013MONBLUE'
      END if

      k = 0
!      print*,'i = ',i,' j = ',j
!      print*,'term 1:(i-1)*nx*ny*num_ver : ',(i-1)*nx*ny*num_ver
      do ii = nx*(i-1)+1,nx*i
!         print*,'term 2: 256*mod(ii-1,nx) : ',256*mod(ii-1,nx)
         do jj = ny*(j-1)+1,ny*j
            !print*,'term 3: (j-1)*ny : ',(j-1)*ny
            k = k + 1
            !print*,k,ii,jj
            !print*,k+nnx*boundary + (ii*2-1)*boundary
            x_analysis((i-1)*nx*ny*num_ver + 256*mod(ii-1,nx) + jj,:)&
                 & = x_local(k+&
                  nny*( nx*(i-1)+1 - max(nx*(i-1)+1-boundary,1) ) +&
                  (ii-nx*(i-1)  )*(ny*(j-1)+1 - max(ny*(j-1)+1-boundary,1) ) &
                  &+&
                  (ii-1-nx*(i-1))*(min(ny*(j)+boundary,256)-ny*j) ,:)
!            x_analysis(k+nx*ny*((i-1)*num_ver+(j-1)),:) = x_local(k+nnx*boundary + (ii*2-1)*boundary,:)
         end do
      end do
      deallocate(x_local)
   end do
end do

pf%psi = x_analysis

end if


end subroutine localise_enkf
