#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/core/drivers/ncdf/file_ncdf_inquire_dim.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_ncdf_inquire_dim(file, dim_name, dim_id, dim_len, is_record_dim)

  USE netcdf, ONLY : nf90_inq_dimid, nf90_inquire_dimension

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Given a dimension name, returns its id and its size
!   If the returned id < 0, the dimension doesn't exist
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_ncdf), INTENT(IN) :: file  ! The file to inspect
  CHARACTER(len=*), INTENT(IN) :: dim_name  ! The name of the dimension to
                                            ! inquire about
  INTEGER, INTENT(OUT) :: dim_id  ! The id of the dimension in the file
  INTEGER, INTENT(OUT) :: dim_len  ! The length of the dimension in the
  LOGICAL, INTENT(OUT) :: is_record_dim  ! Indicates if the named dimension
                                         ! is the record dimension for the file

! Work variables
  INTEGER :: error  ! The current error code (if any)


!-----------------------------------------------------------------------------


  is_record_dim = .FALSE.

! To recover all the required information about a dimension in a NetCDF file
! needs two steps

! First, recover the dimension id from the name
  error = nf90_inq_dimid(file%id, dim_name, dim_id)
  IF ( error /= NF90_NOERR ) THEN
! If there was an error retrieving the dimension id, return a non-existent
! value
    dim_id = -1
    RETURN
  END IF

! If the dimension is the record dimension, we can return early
  IF ( dim_id == file%record_dim ) THEN
    is_record_dim = .TRUE.
    RETURN
  END IF

! If we have a non-record dim, then recover the dim size from the id
  error = nf90_inquire_dimension(file%id, dim_id, LEN=dim_len)
  IF ( error /= NF90_NOERR )                                                  &
    CALL log_fatal_ncdf("file_ncdf_inquire_dim",                              &
                        "Error recovering dimension size", error)

  RETURN

END SUBROUTINE file_ncdf_inquire_dim
#endif
