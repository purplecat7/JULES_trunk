!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-30 13:50:52 pbrowne>
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

subroutine deterministic_model
  use pf_control
  !use Sizes
  use empire_mod !comms

  IMPLICIT NONE
  !include 'mpif.h'
  integer, parameter :: rk = kind(1.0D0)

  integer :: particle,k,tag,mpi_err
  integer :: mpi_status( MPI_STATUS_SIZE )

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

end subroutine deterministic_model
