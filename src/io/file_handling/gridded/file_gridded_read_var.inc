#if !defined(UM_JULES)
! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

!!****************************************************************************
!! Version control information:
!!
!!   $HeadURL: svn://fcm2/JULES_svn/JULES/trunk/src/io/file_handling/gridded/file_gridded_read_var.inc $
!!   $Author: hadmq $
!!
!!   $LastChangedDate: 2012-08-10 12:50:04 +0100 (Fri, 10 Aug 2012) $
!!   $LastChangedRevision: 486 $
!!
!!****************************************************************************

SUBROUTINE file_gridded_read_var(file, var_id, values, extract_subgrid,       &
                                               points, x_coords, y_coords)

  USE file_mod, ONLY : file_read_var

  IMPLICIT NONE

!-----------------------------------------------------------------------------
! Description:
!   Reads data from the given variable in the given file as a 3d array
!   (x,y,levels) where dimensions can have size 1
!
! Current Code Owner: Matt Pryor
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to JULES coding standards v1.
!-----------------------------------------------------------------------------
! Argument types
  TYPE(file_gridded), INTENT(INOUT) :: file
                                 ! The file to read data from
  INTEGER, INTENT(IN) :: var_id  ! The id of the variable to read from
  REAL, INTENT(OUT), TARGET :: values(:,:,:)
                                 ! The values read from the file
                                 ! The dimensions are x,y,levels where
                                 ! the size of the y and/or levels
                                 ! dimensions can be 1, depending on the
                                 ! grid configuration, points requested
                                 ! or variable

  LOGICAL, INTENT(IN) :: extract_subgrid
                                 ! T - extract a subgrid to return from the
                                 !     full grid of the file
                                 ! F - return the full grid

  INTEGER, INTENT(IN), OPTIONAL :: points(:)
                                 ! The indices of points to select from a
                                 ! 1d grid

  INTEGER, INTENT(IN), OPTIONAL :: x_coords(:), y_coords(:)
                                 ! The x and y coordinates of points to
                                 ! select from a 2d grid


! Work variables
  INTEGER :: subgrid_npoints  ! The number of points in the subgrid

  INTEGER :: nlevs  ! The number of levels that the variable has
  INTEGER :: expected_shape(3)  ! The expected shape of values array

  REAL, POINTER :: data(:,:,:)  ! Temporary array to store the full grid
                                ! of data read from file, if required

  INTEGER :: actual_var_id  ! The id of the variable in the underlying file

  INTEGER :: p  ! Loop counter for points


!-----------------------------------------------------------------------------


  subgrid_npoints = 0

!-----------------------------------------------------------------------------
! Check if a subgrid has been requested (in the appropriate way) and that
! the points are in range of the grid
!-----------------------------------------------------------------------------
  IF ( extract_subgrid ) THEN
    IF ( file%grid_is_1d ) THEN
! If we have a 1d grid then we need points
      IF ( .NOT. PRESENT(points) )                                            &
        CALL log_fatal("file_gridded_read_var",                               &
                       "A subgrid has been requested but no points to " //    &
                       "extract specified")

! Check that all requested points are in the grid full grid
      IF ( MINVAL(points) < 1 .OR. MAXVAL(points) > file%grid_nx )            &
        CALL log_fatal("file_gridded_read_var",                               &
                       "A subgrid has been requested but points that are " // &
                       "out of range have been specified")

      subgrid_npoints = SIZE(points)
    ELSE
! If we have a 2d grid we need x_coords and y_coords
      IF ( .NOT. PRESENT(x_coords) .OR. .NOT. PRESENT(y_coords) )             &
        CALL log_fatal("file_gridded_read_var",                               &
                       "A subgrid has been requested but no points to " //    &
                       "extract specified")

! Check that all requested points are in the grid full grid
      IF ( MINVAL(x_coords) < 1 .OR. MAXVAL(x_coords) > file%grid_nx .OR.     &
           MINVAL(y_coords) < 1 .OR. MAXVAL(y_coords) > file%grid_ny )        &
        CALL log_fatal("file_gridded_read_var",                               &
                       "A subgrid has been requested but points that are " // &
                       "out of range have been specified")

      subgrid_npoints = SIZE(x_coords)

! x_coords and y_coords must have the same number of elements
      IF ( SIZE(y_coords) /= subgrid_npoints )                                &
        CALL log_fatal("file_gridded_read_var",                               &
                       "x_coords and y_coords must have the same number " //  &
                       "of elements")
    END IF  ! grid_is_1d
  END IF  ! extract_subgrid

!-----------------------------------------------------------------------------
! Calculate the shape that the array for return data should have and check
! that it actually has that shape
!-----------------------------------------------------------------------------
! Get the number of vertical levels that the variable has
  nlevs = 1
! var_id is an index in the var_ids array of the file_gridded object, as
! returned by file_gridded_def_var
  IF ( file%var_has_levs(var_id) ) nlevs = file%var_nlevs(var_id)

  expected_shape(:) = (/ file%grid_nx, file%grid_ny, nlevs /)

! If a subgrid has been asked for, overwrite the expected sizes of the grid
! dimensions
  IF ( extract_subgrid ) expected_shape(1:2) = (/ subgrid_npoints, 1 /)

  IF ( .NOT. ALL(expected_shape == SHAPE(values)) )                           &
    CALL log_fatal("file_gridded_read_var",                                   &
                   "values has the wrong shape for the requested data")

!-----------------------------------------------------------------------------
! Actually read the data
!-----------------------------------------------------------------------------
! We have to get the variable id in the underlying file from the var_ids array
  actual_var_id = file%var_ids(var_id)

! If we are extracting less than 25% of the points, we read them one at a
! time to reduce IO overhead
  IF ( extract_subgrid .AND.                                                  &
       ( 4 * subgrid_npoints < file%grid_nx * file%grid_ny ) ) THEN
    DO p = 1,subgrid_npoints
      IF ( file%grid_is_1d ) THEN
        CALL file_read_var(file%fh, actual_var_id, values(p,1,:),             &
! The file has one grid dimension, from which we want to read 1 point at the
! given offset. We want to read all the levels for that point.
! If the variable has no levels dimension, it doesn't matter that we have
! defined START and COUNT for that dimension as they will be ignored
                           (/ points(p), 1 /), (/ 1, nlevs /))
      ELSE
        CALL file_read_var(file%fh, actual_var_id, values(p,1,:),             &
! The file has two grid dimensions, from which we want to read 1 point at the
! given coordinates. We want to read all the levels for that point.
! If the variable has no levels dimension, it doesn't matter that we have
! defined START and COUNT for that dimension as they will be ignored
                           (/ x_coords(p), y_coords(p), 1 /), (/ 1, 1, nlevs /))
      END IF
    END DO

! We are done for this case, so return
    RETURN
  END IF


! We only need to allocate data if we are returning a subgrid (since in this
! case we need data to hold the full grid)
! Otherwise we just point data to values
  IF ( extract_subgrid ) THEN
    ALLOCATE(data(file%grid_nx, file%grid_ny, nlevs))
  ELSE
    data => values
  END IF

! Read a slab of data from file
! TODO : Check if these lines could be combined into one read statement
  IF ( file%grid_is_1d ) THEN
    CALL file_read_var(file%fh, actual_var_id, data(:,1,:))
  ELSE
    CALL file_read_var(file%fh, actual_var_id, data(:,:,:))
  END IF

! Populate the values array if required (note that if we are not extracting
! a subgrid, the act of filling data has already filled values, by the
! magic of pointers...!)
  IF ( extract_subgrid ) THEN
! Extract each point
    DO p = 1,subgrid_npoints
      IF ( file%grid_is_1d ) THEN
        values(p,1,:) = data(points(p),1,:)
      ELSE
        values(p,1,:) = data(x_coords(p),y_coords(p),:)
      END IF
    END DO
  END IF

  IF ( extract_subgrid ) DEALLOCATE(data)
  NULLIFY(data)

  RETURN

END SUBROUTINE file_gridded_read_var
#endif
