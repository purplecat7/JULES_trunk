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
subroutine eakf_analysis(num_hor,num_ver,this_hor,this_ver,boundary,x,N&
     &,stateDim,obsDim,rho)
implicit none
integer, parameter :: rk = kind(1.0d0)
integer, intent(in) :: N,stateDim,obsDim
real(kind=rk), dimension(stateDim,N), intent(inout) :: x
integer, intent(in) :: num_hor,num_ver,this_hor,this_ver,boundary
real(kind=rk), intent(in) :: rho

real(kind=rk), dimension(obsDim) :: y
real(kind=rk), dimension(obsDim) :: d,dd
real(kind=rk), dimension(stateDim) :: mean_x
real(kind=rk), dimension(obsDim) :: Hmean_xf
real(kind=rk), dimension(stateDim,N) :: Xpf
real(kind=rk), dimension(N) :: G,S,dn
real(kind=rk), dimension(stateDim,N) :: F,FG,FGU,Xpa
real(kind=rk), dimension(obsDim,N) :: HFG,Ysf
real(kind=rk), dimension(N,obsDim) :: YsfT
real(kind=rk), dimension(obsDim,obsDim) :: VT
real(kind=rk), dimension(N,N) :: WT,U
integer :: LWORK1,LWORK2,INFO
real(kind=rk), dimension(:), allocatable :: WORK1,WORK2
integer :: i

call get_local_observation_data(num_hor,num_ver,this_hor,this_ver,boundary,obsDim,y)


mean_x = 0.0_rk
do i = 1,N
   mean_x = mean_x + x(:,i)
end do
mean_x = mean_x/real(N,rk)

do i = 1,N
   Xpf(:,i) = (x(:,i) - mean_x)/sqrt(real(N-1,rk))
end do

!compute the svd
!JOBU='S' the first min(m,n) columns of U (the left singular
!                  vectors) are returned in the array U;
!JOBVT='S' the first min(m,n) rows of V**T (the right singular
!                  vectors) are returned in the array VT;
!M=stateDim The number of rows of the input matrix A.  M >= 0.
!
!N=N The number of columns of the input matrix A.  N >= 0.
!
!A=Xpf DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!               if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
!                          are destroyed.
!LDA=stateDim The leading dimension of the array A.  LDA >= max(1,M)
!
!S=G DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A, sorted so that S(i) >= S(i+1).
!
!U=F       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
!           (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!           if JOBU = 'S', U contains the first min(m,n) columns of U
!           (the left singular vectors, stored columnwise);
!
!LDU=stateDim The leading dimension of the array U.  LDU >= 1; if
!          JOBU = 'S' or 'A', LDU >= M.
!
!VT=WT DOUBLE PRECISION array, dimension (LDVT,N)
!          if JOBVT = 'S', VT contains the first min(m,n) rows of
!          V**T (the right singular vectors, stored rowwise);
!
!LDVT=N The leading dimension of the array VT.  LDVT >= 1; if
!          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N)
!
!WORK=WORK1 DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
!          superdiagonal elements of an upper bidiagonal matrix B
!          whose diagonal is in S (not necessarily sorted). B
!          satisfies A = U * B * VT, so it has the same singular values
!          as A, and singular vectors related by U and VT.
!
!LWORK=LWORK1 The dimension of the array WORK.
!          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!INFO=INFO (output) INTEGER
!          = 0:  successful exit.
!          < 0:  if INFO = -i, the i-th argument had an illegal value.
!          > 0:  if DBDSQR did not converge, INFO specifies how many
!                superdiagonals of an intermediate bidiagonal form B
!                did not converge to zero. See the description of WORK
!                above for details.

LWORK1 = 2*max( 3*N+stateDim, 5*N )
allocate(WORK1(LWORK1))

call dgesvd('S','S',stateDim,N,Xpf,stateDim,G,F,stateDim,WT,N,WORK1,LWORK1,INFO)
if(INFO .ne. 0) then
   print*,'First SVD failed with info = ',info
   print*,'FYI work(1) = ',WORK1(1)
   stop
end if
deallocate(WORK1)

!FG = F*G
do i = 1,N
   FG(:,i) = F(:,i)*G(i)
end do

call H(N,FG,HFG)
call solve_rhalf(N,HFG,Ysf)

YsfT = transpose(Ysf)

!JOBU='A'    (input) CHARACTER*1
!          Specifies options for computing all or part of the matrix U:
!          = 'A':  all M columns of U are returned in array U:
!
!JOBVT='A'   (input) CHARACTER*1
!          Specifies options for computing all or part of the matrix
!          V**T:
!          = 'A':  all N rows of V**T are returned in the array VT;
!M=N       (input) INTEGER
!          The number of rows of the input matrix A.  M >= 0.
!
!N=obsDim  (input) INTEGER
!          The number of columns of the input matrix A.  N >= 0.
!
!A=YsfT    (input/output) DOUBLE PRECISION array, dimension (LDA,N)
!          On entry, the M-by-N matrix A.
!          On exit,
!          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
!                          are destroyed.
!
!LDA=N     (input) INTEGER
!          The leading dimension of the array A.  LDA >= max(1,M).
!
!S=S       (output) DOUBLE PRECISION array, dimension (min(M,N))
!          The singular values of A, sorted so that S(i) >= S(i+1).
!
!U=U        (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
!          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
!          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
!
!LDU=N     (input) INTEGER
!          The leading dimension of the array U.  LDU >= 1; if
!          JOBU = 'S' or 'A', LDU >= M.
!
!VT=VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
!          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
!          V**T;
!          if JOBVT = 'S', VT contains the first min(m,n) rows of
!          V**T (the right singular vectors, stored rowwise);
!          if JOBVT = 'N' or 'O', VT is not referenced.
!
!LDVT=obsDim    (input) INTEGER
!          The leading dimension of the array VT.  LDVT >= 1; if
!          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
!
!WORK=WORK2    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
!          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
!          superdiagonal elements of an upper bidiagonal matrix B
!          whose diagonal is in S (not necessarily sorted). B
!          satisfies A = U * B * VT, so it has the same singular values
!          as A, and singular vectors related by U and VT.
!
!LWORK=LWORK2   (input) INTEGER
!          The dimension of the array WORK.
!          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
!          For good performance, LWORK should generally be larger.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
!  INFO    (output) INTEGER

LWORK2 = 2*max( 3*N+obsDim, 5*N )
allocate(WORK2(LWORK2))

call dgesvd('A','A',N,obsDim,YsfT,N,S,U,N,VT,ObsDim,WORK2,LWORK2,INFO)

if(INFO .ne. 0) then
   print*,'Second SVD failed with info = ',info
   print*,'FYI work(1) = ',WORK2(1)
   stop
end if
deallocate(WORK2)

!FGU = FG*U
call dgemm('N','N',stateDim,N,N,1.0_rk,FG,stateDim,U,N,0.0_rk,FGU,stateDim)

do i = 1,N
   FGU(:,i) = FGU(:,i) / sqrt(1.0_rk + S(i)**2)
end do

!Xpa = FGU * W'
call dgemm('N','N',stateDim,N,N,1.0_rk,FGU,stateDim,WT,N,0.0_rk,Xpa,stateDim)


call H(1,mean_x,Hmean_xf)
Hmean_xf = y - Hmean_xf
call solve_rhalf(1,Hmean_xf,d)

!dd = V'*d
call dgemv('N',obsDim,obsDim,1.0_rk,VT,obsDim,d,1,0.0_rk,dd,1)

do i = 1,N
   dd(i) = dd(i)/ sqrt(1.0_rk + S(i)**2)
end do

!d = S*dd
!this now makes a vector of size N, call it dn
dn = 0.0_rk
dn = S*dd(1:N)

!mean_xa = mean_x + FGU*d
!note here I am reusing mean_x to be mean_xa
call dgemv('N',stateDim,N,1.0_rk,FGU,stateDim,dn,1,1.0_rk,mean_x,1)

!xa = mean_xa*ones(1,N) + sqrt(N-1)*Xpa;
!x = mean_x*ones(1,N) + sqrt(N-1)*Xpa;

do i = 1,N
   x(:,i) = mean_x + sqrt(real(N-1,rk))*Xpa(:,i)
end do


end subroutine eakf_analysis

