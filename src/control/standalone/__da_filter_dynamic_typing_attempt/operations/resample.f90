!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-10-07 11:17:36 pbrowne>
!!!
!!!    Subroutine to perform Universal Importance Resampling
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
!>    Subroutine to perform Universal Importance Resampling
subroutine resample
!this takes pf%phi with corresponding weights pf%weight
!which are stored as -log(w_i)
!and returns pf%phi which has been resampled with completely equal weights



use pf_control
use random
!use sizes
use empi_wrappers!comms
implicit none
!include 'mpif.h'
integer, parameter :: rk = kind(1.0D0)
integer :: i,old,j,dupepos,dupe,particle,k
real(kind=rk), dimension(pf%nens) :: cumweights
integer, dimension(pf%nens) :: new,tobereplaced,brandnew
logical :: in
real(kind=rk) :: point,draw
integer :: mpi_err
integer, dimension(2*pf%nens) :: requests
integer :: summ,destination,source
integer, allocatable, dimension(:,:) :: statuses
logical :: flag


!print*,'in resample, pf%weight = ',pf%weight

!gather the weights onto the master processor of the particle filter
!call mpi_gatherv(pf%weight(gbldisp(pfrank+1)+1:gbldisp(pfrank+1)+pf%count&
!     &),pf%count,mpi_double_precision,pf%weight&
!     &,gblcount,gbldisp,mpi_double_precision,0,pf_mpi_comm,mpi_err)

     CALL empi_wrap_gatherv_dbl(pf%weight(gbldisp(pfrank+1)+1:gbldisp(pfrank+1)+pf%count), &
        pf%count, pf%weight, 0)
!we have to just draw one random number, so this task is best
!performed solely by the master processor...

!print*,'ahoy ',pf%weight
if(pfrank == 0 ) then
   !check for NaN in the weight
   if(.not. all(pf%weight .eq. pf%weight)) then
      do i = 1,pf%nens
         if(pf%weight(i) .ne. pf%weight(i)) then
            print*,'Particle ',i,' has weight ',pf%weight(i)
            print*,'stopping now'
         end if
      end do
      stop
   end if

!normalise the weights and store them as the actual values.
!print*,'seriously: ',pf%weight
pf%weight = exp(-pf%weight + minval(pf%weight) )
!print*,'un-normalised weights: ',pf%weight
pf%weight = pf%weight/sum(pf%weight)
!print*,'weights = ',pf%weight
!compute the cumulative weights, so cumweights(j) = sum_i=1^j weights(i)
cumweights = 0.0_rk
cumweights(1) = pf%weight(1)
do i = 2,pf%nens
   cumweights(i) = cumweights(i-1) + pf%weight(i)
end do

!print*,'cumweights'
!print*,cumweights

!now make the draw of the single random variable
!distributed normally between 0 and 1/(number of ensemble members)
call random_number(draw)
draw = draw/pf%nens



new = 0
old = 1

!for each particle, compute the POINT given by the 
!random draw.
do particle = 1,pf%nens
   point = real(particle-1,rk)/pf%nens + draw
!   print*,particle,point
   !Step through the weights of each of the particles
   !and find the 'bin' which it lies in.
   !store this in the array NEW.
   do
      if(point .lt. cumweights(old)) then
         new(particle) = old
         exit
      else
         old = old + 1
      end if
   end do
end do


!now we need to determine which particles have to be replaced:
!initialise this array to something stupid:
tobereplaced = -1

!loop over each particle
do i = 1,pf%nens
   in = .false.
   !go through the array NEW and if it does not exist 
   !in that array, then mark the array TOBEREPLACED
   !with the index of said particle
   do j= 1,pf%nens
      if(i .eq. new(j)) then
         in = .true.
         exit
      end if
   end do
   if(.not. in) tobereplaced(i) = i
end do

!dupepos is a counter storing the position of duplicated
dupepos = 1
!generate the array BRANDNEW = [1,2,3,...,N]
brandnew = (/ (i,i=1,pf%nens) /)

!loop through all the particles
do i = 1,pf%nens
   if(tobereplaced(i) .ne. -1) then
      !if the particle is to be replaced then
      !find the particle which has been duplicated
      do k = dupepos,pf%nens
         if(new(k) .eq. new(k+1)) then
            dupe = new(k+1)
            dupepos = k+1
            exit
         end if
      end do
      !place this duplicated particle into the place
      !in the array BRANDNEW that it corresponds to
      brandnew(tobereplaced(i)) = dupe
   end if
end do
!print*,'##################################'
!print*,brandnew
!print*,'##################################'


end if

!now the master processor can send the array
!BRANDNEW to all the other processors
!call mpi_bcast(brandnew,pf%nens,mpi_integer,0,pf_mpi_comm,mpi_err)
!print*,brandnew
CALL empi_wrap_broadcast_int(brandnew, pf%nens, 0)


!check to see if any of this processors particles are being
!replaced. if so, call the mpi receive
k=0
do i = 1,pf%count
   particle = pf%particles(i)
   if(particle .ne. brandnew(particle)) then
      print*,'replacing particle ',particle,'with particle ',brandnew(particle)
      !pf%psi(:,particle) = pf%psi(:,brandnew(particle))
      !pf%psi(:,i) = pf%psi(:,????)
      k=k+1
      print*,'pfrank ',pfrank,'going to receive particle ',i,' from ta&
           &g ',brandnew(particle)
      !the tag will correspond to where it is coming from.
      !who cares where the source was.
!      call mpi_irecv(pf%psi(:,i),model_params%state_dim,mpi_double_precision&
!           &,mpi_any_source,brandnew(particle),pf_mpi_comm,requests(k)&
!           &,mpi_err)
           CALL empi_wrap_open_receive_da_dbl(pf%psi(:,i),model_params%state_dim, &
            brandnew(particle), k)
   end if
end do
!print*,'all the recvs have started'

!print*,'pf%particles(1) = ',pf%particles(1),'pf%particles(pf%count) = &
!     &',pf%particles(pf%count)
!print*,'brandnew(pf%count) = ',brandnew(pf%count)

!now, go through all the particles and see if we have to send
!any of them to some process to replace that particle
do i = 1,pf%nens
   !print*,'i = ',i,'brandnew(i) = ',brandnew(i)
   !check if the new particle is in the range of the particles
   !stored on this process and if it is to be sent somewhere
   if(brandnew(i) .ge. pf%particles(1) .and. brandnew(i) .le.&
        & pf%particles(pf%count) .and. brandnew(i) .ne. i ) then
      !compute the destination!
      !this is the pfrank corresponding to ensemble member i
      destination = 0
      summ = 0
      do j = 1,npfs
!         print*,'computing destination j = ',j
         if(i .le. summ+gblcount(j)) then
!            print*,i,summ+gblcount(j),j-1
            destination = j-1
            exit
         else
            summ = summ + gblcount(j)
!            print*,'summ now is ',summ
         end if
      end do
 
      k=k+1
      print*,'destination'
      print*,destination
      print*,'brandnew(i)'
      print*,brandnew(i)

      !let us find the source to be sent...
      do j = 1,pf%count
         if(pf%particles(j) .eq. brandnew(i)) then
            source = j
            exit
         end if
      end do
      

      print*,k
      print*,'going to send ',brandnew(i),'with tag',brandnew(i)
      print*,'from pfrank ',pfrank,' to destination ',destination

      print*,'the source number is ',source

!      call mpi_isend(pf%psi(:,source),model_params%state_dim&
!           &,mpi_double_precision,destination,brandnew(i),pf_mpi_comm&
!           &,requests(k),mpi_err)
      CALL empi_wrap_open_send_dbl(pf%psi(:,source),model_params%state_dim, destination, &
        brandnew(i), k)

   end if
end do
!print*,'waiting'

!now as all the sends and receives were non-blocking, we must wait
! for them to complete.
if(k.gt.0) then
   print*,'going into mpi_waitall with pfrank = ',pfrank
   ALLOCATE( statuses(MPI_STATUS_SIZE,k) )
   print*,'going into mpi_testall'
   flag = .false.
   do
      if(flag) exit
!      call mpi_testall(k,requests(1:k),flag,statuses,mpi_err)
      CALL empi_wrap_testall(k, flag)
   end do
   print*,'mpi testall finished with flag ',flag
   !call mpi_waitall(k,requests(1:k),statuses(:,1:k),mpi_err)
   deallocate(statuses)
   print*,'did the mpi_waitall on pfrank = ',pfrank
end if
print*,'hit the barrier'
!call mpi_barrier(pf_mpi_comm,mpi_err)
CALL empi_wrap_barrier()

!we have just resampled, so by construction the particles
!have equal weights.
pf%weight = 1.0_rk/real(pf%nens,rk)
pf%weight = -log(pf%weight)
!print*,'finished resampling'
!print*,'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
end subroutine resample
