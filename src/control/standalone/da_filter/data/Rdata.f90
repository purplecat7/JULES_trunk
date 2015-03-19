!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-22 14:18:42 pbrowne>
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
!> @brief Module to hold user supplied data for \f$R\f$
!> observation error covariance matrix
module Rdata
implicit none
integer :: Rn,Rne
integer, allocatable, dimension(:) :: Rrow,Rcol
real(kind=kind(1.0D0)), allocatable, dimension(:) :: Rval,Rdiag
contains
  !> Subroutine to load data for R
  subroutine loadR
    !use sizes
      USE pf_control
    integer :: i
    Rn = model_params%obs_dim
    Rne = 1
    allocate(Rrow(Rne),Rcol(Rne),Rval(Rne))
    allocate(Rdiag(Rn))
    Rrow = (/ (i, i = 1,Rne) /)
    Rcol = (/ (i, i = 1,Rne) /)
    Rval = 0.0D0
    Rdiag = 0.3D0
  end subroutine loadR

  subroutine killR
    !> SUbroutine to deallocate R data
    if(allocated(Rrow)) deallocate(Rrow)
    if(allocated(Rcol)) deallocate(Rcol)
    if(allocated(Rval)) deallocate(Rval)
    if(allocated(Rdiag)) deallocate(Rdiag)
  end subroutine killR

end module Rdata

module hqht_plus_r

  implicit none

contains

  subroutine load_HQHTR
    call HQHTR_factor
  end subroutine load_HQHTR

  subroutine HQHTR_factor

    
  end subroutine HQHTR_factor
  
  subroutine kill_HQHTR

  end subroutine kill_HQHTR




end module hqht_plus_r
