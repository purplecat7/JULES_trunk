!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-22 15:46:43 pbrowne>
!!!
!!!    Module to control what variables are used to generate rank histograms
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

!> Module to control what variables are used to generate rank histograms
module histogram_data
  integer, allocatable, dimension(:) :: rank_hist_list
  integer, allocatable, dimension(:) :: rank_hist_nums
  integer :: rhl_n,rhn_n

contains
  !> subroutine to read from variables_hist.dat which 
  !! variables to be used to make the rank histograms
  subroutine load_histogram_data
    implicit none
    integer :: i

!    rhn_n = 9
    open(2,file='variables_hist.dat',action='read',status='old')
    read(2,'(i7.7)') rhn_n
    allocate(rank_hist_nums(rhn_n))
    do i = 1,rhn_n
       read(2,'(i7.7)') rank_hist_nums(i)
    end do
    rhl_n = sum(rank_hist_nums)
    allocate(rank_hist_list(rhl_n))
    do i = 1,rhl_n
       read(2,'(i7.7)') rank_hist_list(i)
    end do
    close(2)
  end subroutine load_histogram_data

  !>subroutine to clean up arrays used in rank histograms
  subroutine kill_histogram_data
    deallocate(rank_hist_list)
    deallocate(rank_hist_nums)
  end subroutine kill_histogram_data
end module histogram_data
