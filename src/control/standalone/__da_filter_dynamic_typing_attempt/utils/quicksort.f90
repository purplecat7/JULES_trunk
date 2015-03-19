!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Time-stamp: <2014-09-22 15:48:46 pbrowne>
!This code was taken from http://rosettacode.org/wiki/Quicksort#Fortran
!and is distributed under GNU Free Documentation License 1.2.
!see  http://www.gnu.org/licenses/fdl-1.2.html
!> subroutine to sort using the quicksort algorithm
!! @param[in,out] a array of doubles to be sorted
!! @param[in] na dimension of array a
recursive subroutine quicksort_d(a,na)
implicit none 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
real(kind=kind(1.0D0)), dimension(nA), intent(inout) :: A
 
! LOCAL VARIABLES
integer :: left, right, mid
real(kind=kind(1.0D0)) :: pivot, temp
integer :: marker
 
    if (nA > 1) then
        ! insertion sort limit of 47 seems best for sorting 10 million
       !  integers on Intel i7-980X CPU.  Derived data types that use
        ! more memory are optimized with smaller values - around 20 for a 16
       ! -byte type.
        if (nA > 47) then
            ! Do quicksort for large groups
            ! Get median of 1st, mid, & last points for pivot (helps reduce
           !  long execution time on some data sets, such as already
            ! sorted data, over simple 1st point pivot)
            mid = (nA+1)/2
            if (a(mid) >= a(1)) then
                if (a(mid) <= a(nA)) then
                    pivot = a(mid)
                else if (a(nA) > a(1)) then
                    pivot = a(nA)
                else
                    pivot = a(1)
                end if
            else if (a(1) <= a(nA)) then
                pivot = a(1)
            else if (a(nA) > a(mid)) then
                pivot = a(nA)
            else
                pivot = a(mid)
            end if
 
            left = 0
            right = nA + 1
 
            do while (left < right)
                right = right - 1
                do while (A(right) > pivot)
                    right = right - 1
                end do
                left = left + 1
                do while (A(left) < pivot)
                    left = left + 1
                end do
                if (left < right) then
                    temp = A(left)
                    A(left) = A(right)
                    A(right) = temp
                end if
            end do
 
            if (left == right) then
                marker = left + 1
            else
                marker = left
            end if
 
            call quicksort_d(A(:marker-1),marker-1)
            call quicksort_d(A(marker:),nA-marker+1)
 
        else
            call InsertionSort_d(A,nA)    ! Insertion sort for small groups is
            !  faster than Quicksort
        end if
    end if
 
  end subroutine quicksort_d

!> subroutine to sort using the insertionsort algorithm
!! @param[in,out] a array of doubles to be sorted
!! @param[in] na dimension of array a 
subroutine InsertionSort_d(A,nA)
 
! DUMMY ARGUMENTS
integer, intent(in) :: nA
real(kind=kind(1.0D0)), dimension(nA), intent(in out) :: A
 
! LOCAL VARIABLES
real(kind=kind(1.0D0)) :: temp
integer :: i, j
 
    do i = 2, nA
        j = i - 1
        temp = A(i)
        do
            if (j == 0) exit
            if (a(j) <= temp) exit
            A(j+1) = A(j)
            j = j - 1
        end do
        a(j+1) = temp
    end do
 
  end subroutine InsertionSort_d
 
