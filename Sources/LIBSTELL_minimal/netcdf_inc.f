!*******************************************************************************
!>  @file netcdf_inc.f
!>  @brief Contains module @ref netcdf_inc.
!
!  Note separating the Doxygen comment block here so detailed decription is
!  found in the Module not the file.
!
!>  Work around to fix some netcdf import problems with the Ezcdf. This module
!>  emulates the F77 interface using the F90 interface.
!*******************************************************************************
      MODULE netcdf_inc
      USE netcdf

      IMPLICIT NONE

!*******************************************************************************
!  model module parameters
!*******************************************************************************
!>  32bit Real type.
      INTEGER, PARAMETER :: nf_real = nf90_real
!>  32bit Real type.
      INTEGER, PARAMETER :: nf_float = nf90_float
!>  64bit Real type.
      INTEGER, PARAMETER :: nf_double = nf90_double

!>  8bit Integer type.
      INTEGER, PARAMETER :: nf_byte = nf90_byte
!>  Character type.
      INTEGER, PARAMETER :: nf_char = nf90_char
!>  32bit Integer type.
      INTEGER, PARAMETER :: nf_int = nf90_int

!>  Maximum length of strings.
      INTEGER, PARAMETER :: nf_max_name = nf90_max_name
!>  Netcdf global .
      INTEGER, PARAMETER :: nf_global = nf90_global
!>  No error responce.
      INTEGER, PARAMETER :: nf_noerr = nf90_noerr

!>  Write flag.
      INTEGER, PARAMETER :: nf_write = nf90_write
!>  Readonly flag.
      INTEGER, PARAMETER :: nf_nowrite = nf90_nowrite
!>  Clobber flag.
      INTEGER, PARAMETER :: nf_clobber = nf90_clobber
!>  64bit offset flag.
      INTEGER, PARAMETER :: nf_64bit_offset = nf90_64bit_offset

!*******************************************************************************
!  INTERFACE BLOCKS
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  Interface to open a file.
!-------------------------------------------------------------------------------
      INTERFACE nf_open
         MODULE PROCEDURE nf90_open
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to open a file.
!-------------------------------------------------------------------------------
      INTERFACE nf_create
         MODULE PROCEDURE nf90_create
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to open a file.
!-------------------------------------------------------------------------------
      INTERFACE nf_close
         MODULE PROCEDURE nf90_close
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to inquire a dimension id.
!-------------------------------------------------------------------------------
      INTERFACE nf_enddef
         MODULE PROCEDURE nf90_enddef
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to inquire a dimension id.
!-------------------------------------------------------------------------------
      INTERFACE nf_strerror
         MODULE PROCEDURE nf90_strerror
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to inquire a dimension id.
!-------------------------------------------------------------------------------
      INTERFACE nf_inq_dimid
         MODULE PROCEDURE nf90_inq_dimid
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to define a dimension.
!-------------------------------------------------------------------------------
      INTERFACE nf_def_dim
         MODULE PROCEDURE nf90_def_dim
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface to inquire a variable id.
!-------------------------------------------------------------------------------
      INTERFACE nf_inq_varid
         MODULE PROCEDURE nf90_inq_varid
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for reading integer variables.
!-------------------------------------------------------------------------------
      INTERFACE nf_get_var_int
         MODULE PROCEDURE nf_get_var_int, nf_get_var_1d_int,                   &
     &                    nf_get_var_2d_int, nf_get_var_3d_int
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for writing integer variables.
!-------------------------------------------------------------------------------
      INTERFACE nf_put_var_int
         MODULE PROCEDURE nf_put_var_int, nf_put_var_1d_int,                   &
     &                    nf_put_var_2d_int, nf_put_var_3d_int
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for writing integer variables.
!-------------------------------------------------------------------------------
      INTERFACE nf_put_vara_int
         MODULE PROCEDURE nf_put_vara_int, nf_put_vara_3d_int
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for reading real variables. Note double 2D array used by complex
!>  real array values.
!-------------------------------------------------------------------------------
      INTERFACE nf_get_var_real
         MODULE PROCEDURE nf_get_var_cpx_real, nf_get_var_real,                &
     &                    nf_get_var_1d_real, nf_get_var_1d_cpx_real,          &
     &                    nf_get_var_2d_real, nf_get_var_2d_cpx_real,          &
     &                    nf_get_var_2d_double, nf_get_var_3d_real,            &
     &                    nf_get_var_3d_cpx_real
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for reading real variables. Note double 2D array used by complex
!>  real array values.
!-------------------------------------------------------------------------------
      INTERFACE nf_put_vara_real
         MODULE PROCEDURE nf_put_vara_cpx_real, nf_put_vara_real,              &
     &                    nf_put_vara_2d_cpx_real, nf_put_vara_2d_real,        &
     &                    nf_put_vara_3d_cpx_real, nf_put_vara_3d_real
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for writing real variables. Note double 2D array used by complex
!>  real array values.
!-------------------------------------------------------------------------------
      INTERFACE nf_put_var_real
         MODULE PROCEDURE nf_put_var_cpx_real, nf_put_var_real,                &
     &                    nf_put_var_1d_real, nf_put_var_1d_cpx_real,          &
     &                    nf_put_var_2d_real, nf_put_var_2d_cpx_real,          &
     &                    nf_put_var_2d_double, nf_put_var_3d_real,            &
     &                    nf_put_var_3d_cpx_real
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for reading double variables.
!-------------------------------------------------------------------------------
      INTERFACE nf_get_var_double
         MODULE PROCEDURE nf_get_var_cpx_double, nf_get_var_double,            &
     &                    nf_get_var_1D_double,                                &
     &                    nf_get_var_1D_cpx_double,                            &
     &                    nf_get_var_2d_double,                                &
     &                    nf_get_var_2d_cpx_double,                            &
     &                    nf_get_var_3d_double,                                &
     &                    nf_get_var_3d_cpx_double
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for reading double variables. Note double 2D array used by complex
!>  double array values.
!-------------------------------------------------------------------------------
      INTERFACE nf_put_vara_double
         MODULE PROCEDURE nf_put_vara_cpx_double, nf_put_vara_double,          &
     &                    nf_put_vara_2d_cpx_double,                           &
     &                    nf_put_vara_2d_double,                               &
     &                    nf_put_vara_3d_cpx_double,                           &
     &                    nf_put_vara_3d_double
      END INTERFACE

!-------------------------------------------------------------------------------
!>  Interface for writing double variables.
!-------------------------------------------------------------------------------
      INTERFACE nf_put_var_double
         MODULE PROCEDURE nf_put_var_cpx_double, nf_put_var_double,            &
     &                    nf_put_var_1D_double,                                &
     &                    nf_put_var_1D_cpx_double,                            &
     &                    nf_put_var_2d_double,                                &
     &                    nf_put_var_2d_cpx_double,                            &
     &                    nf_put_var_3d_double,                                &
     &                    nf_put_var_3d_cpx_double
      END INTERFACE

      CONTAINS

!*******************************************************************************
!  CREATION SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Define a variable.
!>
!>  @param[in]  ncid   Netcdf file id.
!>  @param[in]  name   Variable id.
!>  @param[in]  xtype  Variable type.
!>  @param[in]  ndims  Number of dimensions.
!>  @param[in]  dimids Dimension ids.
!>  @param[out] varid  Variable id.
!-------------------------------------------------------------------------------
      FUNCTION nf_def_var(ncid, name, xtype, ndims, dimids, varid)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                      :: nf_def_var
      INTEGER, INTENT(in)                          :: ncid
      CHARACTER (len=*), INTENT(in)                :: name
      INTEGER, INTENT(in)                          :: xtype
      INTEGER, INTENT(in)                          :: ndims
      INTEGER, DIMENSION(:), INTENT(in)            :: dimids
      INTEGER, INTENT(out)                         :: varid

!  Start of executable code
      nf_def_var = nf90_def_var(ncid, name, xtype,                             &
     &                          dimids=dimids(:ndims), varid=varid)

      END FUNCTION

!*******************************************************************************
!  QUERY SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Inquire a variable.
!>
!>  @param[in]  ncid   Netcdf file id.
!>  @param[in]  varid  Variable id.
!>  @param[out] name   Variable id.
!>  @param[out] xtype  Variable type.
!>  @param[out] ndims  Number of dimensions.
!>  @param[out] dimids Dimension ids.
!>  @param[out] natts  Number of attributes.
!-------------------------------------------------------------------------------
      FUNCTION nf_inq_var(ncid, varid, name, xtype, ndims, dimids,             &
     &                    natts)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                      :: nf_inq_var
      INTEGER, INTENT(in)                          :: ncid
      INTEGER, INTENT(in)                          :: varid
      CHARACTER (len=*), OPTIONAL, INTENT(out)     :: name
      INTEGER, OPTIONAL, INTENT(out)               :: xtype
      INTEGER, OPTIONAL, INTENT(out)               :: ndims
      INTEGER, DIMENSION(:), OPTIONAL, INTENT(out) :: dimids
      INTEGER, OPTIONAL, INTENT(out)               :: natts

!  Start of executable code
      nf_inq_var = nf90_inquire_variable(ncid, varid, name, xtype,             &
     &                                   ndims, dimids, natts)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Inquire a dimension name or length.
!>
!>  @param[in]  ncid   Netcdf file id.
!>  @param[in]  dimid  Dimension id.
!>  @param[out] name   Name of the dimension.
!>  @param[out] length Length of the dimension.
!-------------------------------------------------------------------------------
      FUNCTION nf_inq_dim(ncid, dimid, name, length)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                  :: nf_inq_dim
      INTEGER, INTENT(in)                      :: ncid
      INTEGER, INTENT(in)                      :: dimid
      CHARACTER (len=*), OPTIONAL, INTENT(out) :: name
      INTEGER, OPTIONAL, INTENT(out)           :: length

!  Start of executable code
      nf_inq_dim = nf90_inquire_dimension(ncid, dimid, name, length)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Inquire a dimension length.
!>
!>  @param[in]  ncid   Netcdf file id.
!>  @param[in]  dimid  Dimension id.
!>  @param[out] length Length of the dimension.
!-------------------------------------------------------------------------------
      FUNCTION nf_inq_dimlen(ncid, dimid, length)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER              :: nf_inq_dimlen
      INTEGER, INTENT(in)  :: ncid
      INTEGER, INTENT(in)  :: dimid
      INTEGER, INTENT(out) :: length

!  Start of executable code
      nf_inq_dimlen = nf90_inquire_dimension(ncid, dimid, len=length)

      END FUNCTION

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Read an string attribute.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[in]  name  Name of the attribute.
!>  @param[out] text  String values of the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_att_text(ncid, varid, name, text)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                        :: nf_get_att_text
      INTEGER, INTENT(in)            :: ncid
      INTEGER, INTENT(in)            :: varid
      CHARACTER (len=*), INTENT(in)  :: name
      CHARACTER (len=*), INTENT(out) :: text

!  Start of executable code
      nf_get_att_text = nf90_get_att(ncid, varid, name, text)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an double attribute.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[in]  name  Name of the attribute.
!>  @param[out] dvals Double values of the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_att_double(ncid, varid, name, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                  :: nf_get_att_double
      INTEGER, INTENT(in)                      :: ncid
      INTEGER, INTENT(in)                      :: varid
      CHARACTER (len=*), INTENT(in)            :: name
      REAL (REAL64), DIMENSION(:), INTENT(out) :: dvals

!  Start of executable code
      nf_get_att_double = nf90_get_att(ncid, varid, name, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an real attribute.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[in]  name  Name of the attribute.
!>  @param[out] rvals Real values of the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_att_real(ncid, varid, name, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                  :: nf_get_att_real
      INTEGER, INTENT(in)                      :: ncid
      INTEGER, INTENT(in)                      :: varid
      CHARACTER (len=*), INTENT(in)            :: name
      REAL (REAL32), DIMENSION(:), INTENT(out) :: rvals

!  Start of executable code
      nf_get_att_real = nf90_get_att(ncid, varid, name, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an integer attribute.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[in]  name  Name of the attribute.
!>  @param[out] ivals Integer values of the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_att_int(ncid, varid, name, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                            :: nf_get_att_int
      INTEGER, INTENT(in)                :: ncid
      INTEGER, INTENT(in)                :: varid
      CHARACTER (len=*), INTENT(in)      :: name
      INTEGER, DIMENSION(:), INTENT(out) :: ivals

!  Start of executable code
      nf_get_att_int = nf90_get_att(ncid, varid, name, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a double value.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] dvals Double to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                    :: nf_get_var_double
      INTEGER, INTENT(in)        :: ncid
      INTEGER, INTENT(in)        :: varid
      REAL (REAL64), INTENT(out) :: dvals

!  Start of executable code
      nf_get_var_double = nf90_get_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a double valued array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] dvals Double array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_1d_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_1d_double
      INTEGER, INTENT(in)                      :: ncid
      INTEGER, INTENT(in)                      :: varid
      REAL (REAL64), DIMENSION(:), INTENT(out) :: dvals

!  Start of executable code
      nf_get_var_1d_double = nf90_get_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a double valued 2D array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] dvals Double array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_2d_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_2d_double
      INTEGER, INTENT(in)                        :: ncid
      INTEGER, INTENT(in)                        :: varid
      REAL (REAL64), DIMENSION(:,:), INTENT(out) :: dvals

!  Start of executable code
      nf_get_var_2d_double = nf90_get_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a double valued 3D array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] dvals Double array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_3d_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_3d_double
      INTEGER, INTENT(in)                          :: ncid
      INTEGER, INTENT(in)                          :: varid
      REAL (REAL64), DIMENSION(:,:,:), INTENT(out) :: dvals

!  Start of executable code
      nf_get_var_3d_double = nf90_get_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 64bit complex value.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: nf_get_var_cpx_double
      INTEGER, INTENT(in)           :: ncid
      INTEGER, INTENT(in)           :: varid
      COMPLEX (REAL64), INTENT(out) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(2)   :: dvals

!  Start of executable code
      nf_get_var_cpx_double = nf90_get_var(ncid, varid, dvals)
      cvals = CMPLX(dvals(1), dvals(2))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 64bit complex valued array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a double array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_1d_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_1d_cpx_double
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      COMPLEX (REAL64), DIMENSION(:), INTENT(out) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:), ALLOCATABLE    :: dvals
      INTEGER                                     :: i

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals)*2))

      nf_get_var_1d_cpx_double = nf90_get_var(ncid, varid, dvals)
      DO i=1, SIZE(dvals), 2
         cvals((i + 1)/2) = CMPLX(dvals(i), dvals(i + 1))
      END DO

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 64bit complex valued 2D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a double array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_2d_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_2d_cpx_double
      INTEGER, INTENT(in)                           :: ncid
      INTEGER, INTENT(in)                           :: varid
      COMPLEX (REAL64), DIMENSION(:,:), INTENT(out) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:,:), ALLOCATABLE    :: dvals
      INTEGER                                       :: i
      INTEGER                                       :: j

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals, 1)*2, SIZE(cvals, 2)))

      nf_get_var_2d_cpx_double = nf90_get_var(ncid, varid, dvals)
      DO j=1, SIZE(dvals, 2), 2
         DO i=1, SIZE(dvals, 1), 2
            cvals((i + 1)/2,j) = CMPLX(dvals(i,j), dvals(i + 1,j))
         END DO
      END DO

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 64bit complex valued 3D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a double array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_3d_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_3d_cpx_double
      INTEGER, INTENT(in)                             :: ncid
      INTEGER, INTENT(in)                             :: varid
      COMPLEX (REAL64), DIMENSION(:,:,:), INTENT(out) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:,:,:), ALLOCATABLE    :: dvals
      INTEGER                                         :: i
      INTEGER                                         :: j
      INTEGER                                         :: k

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals,1)*2, SIZE(cvals,2), SIZE(cvals,3)))

      nf_get_var_3d_cpx_double = nf90_get_var(ncid, varid, dvals)
      DO k=1, SIZE(dvals, 3), 2
         DO j=1, SIZE(dvals, 2), 2
            DO i=1, SIZE(dvals, 1), 2
               cvals((i + 1)/2,j,k) = CMPLX(dvals(i,j,k),                      &
     &                                      dvals(i + 1,j,k))
            END DO
         END DO
      END DO

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a real value.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] rvals Real to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                    :: nf_get_var_real
      INTEGER, INTENT(in)        :: ncid
      INTEGER, INTENT(in)        :: varid
      REAL (REAL32), INTENT(out) :: rvals

!  Start of executable code
      nf_get_var_real = nf90_get_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a real valued array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] rvals Real array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_1d_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                  :: nf_get_var_1d_real
      INTEGER, INTENT(in)                      :: ncid
      INTEGER, INTENT(in)                      :: varid
      REAL (REAL32), DIMENSION(:), INTENT(out) :: rvals

!  Start of executable code
      nf_get_var_1d_real = nf90_get_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a real valued 2D array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] rvals Real 2D array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_2d_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_2d_real
      INTEGER, INTENT(in)                        :: ncid
      INTEGER, INTENT(in)                        :: varid
      REAL (REAL32), DIMENSION(:,:), INTENT(out) :: rvals

!  Start of executable code
      nf_get_var_2d_real = nf90_get_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a real valued 3D array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] rvals Real 3D array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_3d_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_3d_real
      INTEGER, INTENT(in)                          :: ncid
      INTEGER, INTENT(in)                          :: varid
      REAL (REAL32), DIMENSION(:,:,:), INTENT(out) :: rvals

!  Start of executable code
      nf_get_var_3d_real = nf90_get_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 32bit complex value.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: nf_get_var_cpx_real
      INTEGER, INTENT(in)           :: ncid
      INTEGER, INTENT(in)           :: varid
      COMPLEX (REAL32), INTENT(out) :: cvals

      REAL (REAL32), DIMENSION(2)   :: rvals

!  Start of executable code
      nf_get_var_cpx_real = nf90_get_var(ncid, varid, rvals)
      cvals = CMPLX(rvals(1), rvals(2))

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 32bit complex valued array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_1d_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_1d_cpx_real
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      COMPLEX (REAL32), DIMENSION(:), INTENT(out) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:), ALLOCATABLE    :: rvals
      INTEGER                                     :: i

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals)*2))

      nf_get_var_1d_cpx_real = nf90_get_var(ncid, varid, rvals)
      DO i=1, SIZE(rvals), 2
         cvals((i + 1)/2) = CMPLX(rvals(i), rvals(i + 1))
      END DO

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 32bit complex valued 2D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_2d_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_2d_cpx_real
      INTEGER, INTENT(in)                           :: ncid
      INTEGER, INTENT(in)                           :: varid
      COMPLEX (REAL32), DIMENSION(:,:), INTENT(out) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:,:), ALLOCATABLE    :: rvals
      INTEGER                                       :: i
      INTEGER                                       :: j

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals, 1)*2, SIZE(cvals, 2)))

      nf_get_var_2d_cpx_real = nf90_get_var(ncid, varid, rvals)
      DO j=1, SIZE(rvals, 2)
         DO i=1, SIZE(rvals, 1), 2
            cvals((i + 1)/2, j) = CMPLX(rvals(i, j), rvals(i + 1, j))
         END DO
      END DO

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a 32bit complex valued 3D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_3d_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_get_var_3d_cpx_real
      INTEGER, INTENT(in)                             :: ncid
      INTEGER, INTENT(in)                             :: varid
      COMPLEX (REAL32), DIMENSION(:,:,:), INTENT(out) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:,:,:), ALLOCATABLE    :: rvals
      INTEGER                                         :: i
      INTEGER                                         :: j
      INTEGER                                         :: k

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals, 1)*2, SIZE(cvals, 2), SIZE(cvals, 3)))

      nf_get_var_3d_cpx_real = nf90_get_var(ncid, varid, rvals)
      DO k=1, SIZE(rvals, 3)
         DO j=1, SIZE(rvals, 2)
            DO i=1, SIZE(rvals, 1), 2
               cvals((i + 1)/2,j,k) = CMPLX(rvals(i,j,k),                      &
     &                                      rvals(i + 1,j,k))
            END DO
         END DO
      END DO

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an integer value.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] ivals Integer to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER              :: nf_get_var_int
      INTEGER, INTENT(in)  :: ncid
      INTEGER, INTENT(in)  :: varid
      INTEGER, INTENT(out) :: ivals

!  Start of executable code
      nf_get_var_int = nf90_get_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an integer valued array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] ivals Integer array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_1d_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                            :: nf_get_var_1d_int
      INTEGER, INTENT(in)                :: ncid
      INTEGER, INTENT(in)                :: varid
      INTEGER, DIMENSION(:), INTENT(out) :: ivals

!  Start of executable code
      nf_get_var_1d_int = nf90_get_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an integer valued 2D array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] ivals Integer array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_2d_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                              :: nf_get_var_2d_int
      INTEGER, INTENT(in)                  :: ncid
      INTEGER, INTENT(in)                  :: varid
      INTEGER, DIMENSION(:,:), INTENT(out) :: ivals

!  Start of executable code
      nf_get_var_2d_int = nf90_get_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read an integer valued 3D array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] ivals Integer array to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_3d_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                :: nf_get_var_3d_int
      INTEGER, INTENT(in)                    :: ncid
      INTEGER, INTENT(in)                    :: varid
      INTEGER, DIMENSION(:,:,:), INTENT(out) :: ivals

!  Start of executable code
      nf_get_var_3d_int = nf90_get_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a text value.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] text  String to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_var_text(ncid, varid, text)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                        :: nf_get_var_text
      INTEGER, INTENT(in)            :: ncid
      INTEGER, INTENT(in)            :: varid
      CHARACTER (len=*), INTENT(out) :: text

!  Start of executable code
      nf_get_var_text = nf90_get_var(ncid, varid, text)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a text valued array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[in]  start Starting index of the string array.
!>  @param[in]  count Number of elements to count.
!>  @param[out] text  String to store the value.
!-------------------------------------------------------------------------------
      FUNCTION nf_get_vara_text(ncid, varid, start, count, text)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                           :: nf_get_vara_text
      INTEGER, INTENT(in)               :: ncid
      INTEGER, INTENT(in)               :: varid
      INTEGER, DIMENSION(:), INTENT(in) :: start
      INTEGER, DIMENSION(:), INTENT(in) :: count
      CHARACTER (len=*), INTENT(out)    :: text

!  Start of executable code
      nf_get_vara_text = nf90_get_var(ncid, varid, text, start=start,          &
     &                                count=count)

      END FUNCTION

!*******************************************************************************
!  SETTER SUBROUTINES
!*******************************************************************************
!-------------------------------------------------------------------------------
!>  @brief Put an text attribute.
!>
!>  @param[in] ncid   Netcdf file id.
!>  @param[in] varid  Variable id.
!>  @param[in] name   Name of the atribute.
!>  @param[in] length Length of the real array.
!>  @param[in] text   String value to asign to the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_att_text(ncid, varid, name, length, text)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: nf_put_att_text
      INTEGER, INTENT(in)           :: ncid
      INTEGER, INTENT(in)           :: varid
      CHARACTER (len=*), INTENT(in) :: name
      INTEGER, INTENT(in)           :: length
      CHARACTER (len=*), INTENT(in) :: text

!  Start of executable code
      nf_put_att_text = nf90_put_att(ncid, varid, name, text)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Put an double attribute.
!>
!>  @param[in] ncid   Netcdf file id.
!>  @param[in] varid  Variable id.
!>  @param[in] name   Name of the atribute.
!>  @param[in] xtype  Type of the value. This is not used.
!>  @param[in] length Length of the real array.
!>  @param[in] dvals  Real values to asign to the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_att_double(ncid, varid, name, xtype, length,             &
     &                           dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                 :: nf_put_att_double
      INTEGER, INTENT(in)                     :: ncid
      INTEGER, INTENT(in)                     :: varid
      CHARACTER (len=*), INTENT(in)           :: name
      INTEGER, INTENT(in)                     :: xtype
      INTEGER, INTENT(in)                     :: length
      REAL (REAL64), DIMENSION(:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_att_double = nf90_put_att(ncid, varid, name, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Put an real attribute.
!>
!>  @param[in] ncid   Netcdf file id.
!>  @param[in] varid  Variable id.
!>  @param[in] name   Name of the atribute.
!>  @param[in] xtype  Type of the value. This is not used.
!>  @param[in] length Length of the real array.
!>  @param[in] rvals  Real values to asign to the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_att_real(ncid, varid, name, xtype, length, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                 :: nf_put_att_real
      INTEGER, INTENT(in)                     :: ncid
      INTEGER, INTENT(in)                     :: varid
      CHARACTER (len=*), INTENT(in)           :: name
      INTEGER, INTENT(in)                     :: xtype
      INTEGER, INTENT(in)                     :: length
      REAL (REAL32), DIMENSION(:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_att_real = nf90_put_att(ncid, varid, name, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Put an integer attribute.
!>
!>  @param[in] ncid   Netcdf file id.
!>  @param[in] varid  Variable id.
!>  @param[in] name   Name of the atribute.
!>  @param[in] xtype  Type of the value. This is not used.
!>  @param[in] length Length of the integer array.
!>  @param[in] ivals  Integer values to asign to the attribute.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_att_int(ncid, varid, name, xtype, length, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                           :: nf_put_att_int
      INTEGER, INTENT(in)               :: ncid
      INTEGER, INTENT(in)               :: varid
      CHARACTER (len=*), INTENT(in)     :: name
      INTEGER, INTENT(in)               :: xtype
      INTEGER, INTENT(in)               :: length
      INTEGER, DIMENSION(:), INTENT(in) :: ivals

!  Start of executable code
      nf_put_att_int = nf90_put_att(ncid, varid, name, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write an integer value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] ivals Integer value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER             :: nf_put_var_int
      INTEGER, INTENT(in) :: ncid
      INTEGER, INTENT(in) :: varid
      INTEGER, INTENT(in) :: ivals

!  Start of executable code
      nf_put_var_int = nf90_put_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write an integer valued array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] ivals Integer value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_1d_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                           :: nf_put_var_1d_int
      INTEGER, INTENT(in)               :: ncid
      INTEGER, INTENT(in)               :: varid
      INTEGER, DIMENSION(:), INTENT(in) :: ivals

!  Start of executable code
      nf_put_var_1d_int = nf90_put_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write an integer valued 2D array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] ivals Integer value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_2d_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                             :: nf_put_var_2d_int
      INTEGER, INTENT(in)                 :: ncid
      INTEGER, INTENT(in)                 :: varid
      INTEGER, DIMENSION(:,:), INTENT(in) :: ivals

!  Start of executable code
      nf_put_var_2d_int = nf90_put_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write an integer valued 3D array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] ivals Integer value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_3d_int(ncid, varid, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                               :: nf_put_var_3d_int
      INTEGER, INTENT(in)                   :: ncid
      INTEGER, INTENT(in)                   :: varid
      INTEGER, DIMENSION(:,:,:), INTENT(in) :: ivals

!  Start of executable code
      nf_put_var_3d_int = nf90_put_var(ncid, varid, ivals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a integer value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] ivals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_int(ncid, varid, start, count, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                             :: nf_put_vara_int
      INTEGER, INTENT(in)                 :: ncid
      INTEGER, INTENT(in)                 :: varid
      INTEGER, DIMENSION(:), INTENT(in)   :: start
      INTEGER, DIMENSION(:), INTENT(in)   :: count
      INTEGER, DIMENSION(:,:), INTENT(in) :: ivals

!  Start of executable code
      nf_put_vara_int = nf90_put_var(ncid, varid, ivals, start=start,          &
     &                               count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a integer value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] ivals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_3d_int(ncid, varid, start, count, ivals)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                               :: nf_put_vara_3d_int
      INTEGER, INTENT(in)                   :: ncid
      INTEGER, INTENT(in)                   :: varid
      INTEGER, DIMENSION(:), INTENT(in)     :: start
      INTEGER, DIMENSION(:), INTENT(in)     :: count
      INTEGER, DIMENSION(:,:,:), INTENT(in) :: ivals

!  Start of executable code
      nf_put_vara_3d_int = nf90_put_var(ncid, varid, ivals, start=start,       &
     &                                  count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] dvals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                   :: nf_put_var_double
      INTEGER, INTENT(in)       :: ncid
      INTEGER, INTENT(in)       :: varid
      REAL (REAL64), INTENT(in) :: dvals

!  Start of executable code
      nf_put_var_double = nf90_put_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double valued array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] dvals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_1d_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_1d_double
      INTEGER, INTENT(in)                     :: ncid
      INTEGER, INTENT(in)                     :: varid
      REAL (REAL64), DIMENSION(:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_var_1d_double = nf90_put_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a double valued 2D array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] dvals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_2d_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_2d_double
      INTEGER, INTENT(in)                       :: ncid
      INTEGER, INTENT(in)                       :: varid
      REAL (REAL64), DIMENSION(:,:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_var_2d_double = nf90_put_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double valued 3D array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] dvals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_3d_double(ncid, varid, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_3d_double
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      REAL (REAL64), DIMENSION(:,:,:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_var_3d_double = nf90_put_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 64bit complex value.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                      :: nf_put_var_cpx_double
      INTEGER, INTENT(in)          :: ncid
      INTEGER, INTENT(in)          :: varid
      COMPLEX (REAL64), INTENT(in) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(2)  :: dvals

!  Start of executable code
      dvals(1) = REAL(cvals)
      dvals(2) = AIMAG(cvals)
      nf_put_var_cpx_double = nf90_put_var(ncid, varid, dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 64bit complex valued array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a double array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_1d_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_1d_cpx_double
      INTEGER, INTENT(in)                        :: ncid
      INTEGER, INTENT(in)                        :: varid
      COMPLEX (REAL64), DIMENSION(:), INTENT(in) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:), ALLOCATABLE   :: dvals
      INTEGER                                    :: i

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals)*2))

      DO i=1, SIZE(dvals, 1), 2
         dvals(i) = REAL(cvals((i + 1)/2))
         dvals(i + 1) = AIMAG(cvals((i + 1)/2))
      END DO

      nf_put_var_1d_cpx_double = nf90_put_var(ncid, varid, dvals)

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 64bit complex valued 2D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a double array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] cvals Complex value to write.
!-------------------------------------------------------------------------------
       FUNCTION nf_put_var_2d_cpx_double(ncid, varid, cvals)
       USE iso_fortran_env, only: REAL64

       IMPLICIT NONE

!  Declare Arguments
       INTEGER :: nf_put_var_2d_cpx_double
       INTEGER, INTENT(in)                          :: ncid
       INTEGER, INTENT(in)                          :: varid
       COMPLEX (REAL64), DIMENSION(:,:), INTENT(in) :: cvals

!  local variables
       REAL (REAL64), DIMENSION(:,:), ALLOCATABLE   :: dvals
       INTEGER                                      :: i
       INTEGER                                      :: j

!  Start of executable code
       ALLOCATE(dvals(SIZE(cvals, 1)*2, SIZE(cvals, 2)))

       DO j=1, SIZE(dvals, 2)
          DO i=1, SIZE(dvals, 1), 2
              dvals(i,j) = REAL(cvals((i + 1)/2,j))
              dvals(i + 1,j) = AIMAG(cvals((i + 1)/2,j))
          END DO
       END DO

       nf_put_var_2d_cpx_double = nf90_put_var(ncid, varid, dvals)

       DEALLOCATE(dvals)

       END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 64bit complex valued 3D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a double array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_3d_cpx_double(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_3d_cpx_double
      INTEGER, INTENT(in)                            :: ncid
      INTEGER, INTENT(in)                            :: varid
      COMPLEX (REAL64), DIMENSION(:,:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:,:,:), ALLOCATABLE   :: dvals
      INTEGER                                        :: i
      INTEGER                                        :: j
      INTEGER                                        :: k

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals,1)*2, SIZE(cvals,2), SIZE(cvals,3)))

      DO k=1, SIZE(dvals, 3)
         DO j=1, SIZE(dvals, 2)
            DO i=1, SIZE(dvals, 1), 2
               dvals(i,j,k) = REAL(cvals((i + 1)/2,j,k))
               dvals(i + 1,j,k) = AIMAG(cvals((i + 1)/2,j,k))
            END DO
         END DO
      END DO

      nf_put_var_3d_cpx_double = nf90_put_var(ncid, varid, dvals)

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] dvals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_double(ncid, varid, start, count, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                 :: nf_put_vara_double
      INTEGER, INTENT(in)                     :: ncid
      INTEGER, INTENT(in)                     :: varid
      INTEGER, DIMENSION(:), INTENT(in)       :: start
      INTEGER, DIMENSION(:), INTENT(in)       :: count
      REAL (REAL64), DIMENSION(:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_vara_double = nf90_put_var(ncid, varid, dvals, start=start,       &
     &                                  count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] dvals Double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_2d_double(ncid, varid, start, count, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                   :: nf_put_vara_2d_double
      INTEGER, INTENT(in)                       :: ncid
      INTEGER, INTENT(in)                       :: varid
      INTEGER, DIMENSION(:), INTENT(in)         :: start
      INTEGER, DIMENSION(:), INTENT(in)         :: count
      REAL (REAL64), DIMENSION(:,:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_vara_2d_double = nf90_put_var(ncid, varid, dvals,                 &
     &                                     start=start, count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] dvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_3d_double(ncid, varid, start, count, dvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_vara_3d_double
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      INTEGER, DIMENSION(:), INTENT(in)           :: start
      INTEGER, DIMENSION(:), INTENT(in)           :: count
      REAL (REAL64), DIMENSION(:,:,:), INTENT(in) :: dvals

!  Start of executable code
      nf_put_vara_3d_double = nf90_put_var(ncid, varid, dvals,                 &
     &                                     start=start, count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a complex double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] cvals Complex double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_cpx_double(ncid, varid, start, count, cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_vara_cpx_double
      INTEGER, INTENT(in)                        :: ncid
      INTEGER, INTENT(in)                        :: varid
      INTEGER, DIMENSION(:), INTENT(in)          :: start
      INTEGER, DIMENSION(:), INTENT(in)          :: count
      COMPLEX (REAL64), DIMENSION(:), INTENT(in) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:), ALLOCATABLE   :: dvals
      INTEGER                                    :: i

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals)*2))

      DO i = 1, SIZE(dvals), 2
         dvals(i) = REAL(cvals((i + 1)/2))
         dvals(i + 1) = AIMAG(cvals((i + 1)/2))
      END DO

      nf_put_vara_cpx_double = nf90_put_var(ncid, varid, dvals,                &
     &                                      start=start, count=count)

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a complex double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] cvals Complex double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_2d_cpx_double(ncid, varid, start, count,            &
     &                                   cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_vara_2d_cpx_double
      INTEGER, INTENT(in)                          :: ncid
      INTEGER, INTENT(in)                          :: varid
      INTEGER, DIMENSION(:), INTENT(in)            :: start
      INTEGER, DIMENSION(:), INTENT(in)            :: count
      COMPLEX (REAL64), DIMENSION(:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:,:), ALLOCATABLE   :: dvals
      INTEGER                                      :: i
      INTEGER                                      :: j

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals, 1)*2, SIZE(cvals, 2)))

      DO j = 1, SIZE(dvals, 2)
         DO i = 1, SIZE(dvals, 1), 2
            dvals(i,j) = REAL(cvals((i + 1)/2, j))
            dvals(i + 1,j) = AIMAG(cvals((i + 1)/2, j))
         END DO
      END DO

      nf_put_vara_2d_cpx_double = nf90_put_var(ncid, varid, dvals,             &
     &                                         start=start, count=count)

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a double value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] cvals Complex double value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_3d_cpx_double(ncid, varid, start, count,            &
     &                                   cvals)
      USE iso_fortran_env, only: REAL64

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_vara_3d_cpx_double
      INTEGER, INTENT(in)                            :: ncid
      INTEGER, INTENT(in)                            :: varid
      INTEGER, DIMENSION(:), INTENT(in)              :: start
      INTEGER, DIMENSION(:), INTENT(in)              :: count
      COMPLEX (REAL64), DIMENSION(:,:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL64), DIMENSION(:,:,:), ALLOCATABLE   :: dvals
      INTEGER                                        :: i
      INTEGER                                        :: j
      INTEGER                                        :: k

!  Start of executable code
      ALLOCATE(dvals(SIZE(cvals, 1)*2, SIZE(cvals, 2), SIZE(cvals, 3)))

      DO k = 1, SIZE(dvals, 3)
         DO j = 1, SIZE(dvals, 2)
            DO i = 1, SIZE(dvals, 1), 2
               dvals(i,j,k) = REAL(cvals((i + 1)/2,j,k))
               dvals(i + 1,j,k) = AIMAG(cvals((i + 1)/2,j,k))
            END DO
         END DO
      END DO

      nf_put_vara_3d_cpx_double = nf90_put_var(ncid, varid, dvals,             &
     &                                         start=start, count=count)

      DEALLOCATE(dvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                   :: nf_put_var_real
      INTEGER, INTENT(in)       :: ncid
      INTEGER, INTENT(in)       :: varid
      REAL (REAL32), INTENT(in) :: rvals

!  Start of executable code
      nf_put_var_real = nf90_put_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a real valued array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_1d_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_1d_real
      INTEGER, INTENT(in)                     :: ncid
      INTEGER, INTENT(in)                     :: varid
      REAL (REAL32), DIMENSION(:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_var_1d_real = nf90_put_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Read a real valued 2D array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_2d_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_2d_real
      INTEGER, INTENT(in)                       :: ncid
      INTEGER, INTENT(in)                       :: varid
      REAL (REAL32), DIMENSION(:,:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_var_2d_real = nf90_put_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a real valued 3D array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_3d_real(ncid, varid, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_3d_real
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      REAL (REAL32), DIMENSION(:,:,:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_var_3d_real = nf90_put_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 32bit complex value.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                      :: nf_put_var_cpx_real
      INTEGER, INTENT(in)          :: ncid
      INTEGER, INTENT(in)          :: varid
      COMPLEX (REAL32), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(2)  :: rvals

!  Start of executable code
      rvals(1) = REAL(cvals)
      rvals(2) = AIMAG(cvals)
      nf_put_var_cpx_real = nf90_put_var(ncid, varid, rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 32bit complex valued array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_1d_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_1d_cpx_real
      INTEGER, INTENT(in)                        :: ncid
      INTEGER, INTENT(in)                        :: varid
      COMPLEX (REAL32), DIMENSION(:), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:), ALLOCATABLE   :: rvals
      INTEGER                                    :: i

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals)*2))

      DO i=1, SIZE(rvals, 1), 2
         rvals(i) = REAL(cvals((i + 1)/2))
         rvals(i + 1) = AIMAG(cvals((i + 1)/2))
      END DO

      nf_put_var_1d_cpx_real = nf90_put_var(ncid, varid, rvals)

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 32bit complex valued 2D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_2d_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_2d_cpx_real
      INTEGER, INTENT(in)                          :: ncid
      INTEGER, INTENT(in)                          :: varid
      COMPLEX (REAL32), DIMENSION(:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:,:), ALLOCATABLE   :: rvals
      INTEGER                                      :: i
      INTEGER                                      :: j

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals, 1)*2, SIZE(cvals, 2)))

      DO j=1, SIZE(rvals, 2)
         DO i=1, SIZE(rvals, 1), 2
            rvals(i,j) = REAL(cvals((i + 1)/2,j))
            rvals(i + 1,j) = AIMAG(cvals((i + 1)/2,j))
         END DO
      END DO

      nf_put_var_2d_cpx_real = nf90_put_var(ncid, varid, rvals)

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a 32bit complex valued 3D array.
!>
!>  Netcdf does not support complex types directly. Instead these values are
!>  interweaved into a real array.
!>
!>  @param[in]  ncid  Netcdf file id.
!>  @param[in]  varid Variable id.
!>  @param[out] cvals Complex value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_3d_cpx_real(ncid, varid, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_var_3d_cpx_real
      INTEGER, INTENT(in)                            :: ncid
      INTEGER, INTENT(in)                            :: varid
      COMPLEX (REAL32), DIMENSION(:,:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:,:,:), ALLOCATABLE   :: rvals
      INTEGER                                        :: i
      INTEGER                                        :: j
      INTEGER                                        :: k

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals,1)*2, SIZE(cvals,2), SIZE(cvals,3)))

      DO k=1, SIZE(rvals, 3)
         DO j=1, SIZE(rvals, 2)
            DO i=1, SIZE(rvals, 1), 2
               rvals(i,j,k) = REAL(cvals((i + 1)/2,j,k))
               rvals(i + 1,j,k) = AIMAG(cvals((i + 1)/2,j,k))
            END DO
         END DO
      END DO

      nf_put_var_3d_cpx_real = nf90_put_var(ncid, varid, rvals)

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_real(ncid, varid, start, count, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                 :: nf_put_vara_real
      INTEGER, INTENT(in)                     :: ncid
      INTEGER, INTENT(in)                     :: varid
      INTEGER, DIMENSION(:), INTENT(in)       :: start
      INTEGER, DIMENSION(:), INTENT(in)       :: count
      REAL (REAL32), DIMENSION(:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_vara_real = nf90_put_var(ncid, varid, rvals, start=start,         &
     &                                count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_2d_real(ncid, varid, start, count, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                   :: nf_put_vara_2d_real
      INTEGER, INTENT(in)                       :: ncid
      INTEGER, INTENT(in)                       :: varid
      INTEGER, DIMENSION(:), INTENT(in)         :: start
      INTEGER, DIMENSION(:), INTENT(in)         :: count
      REAL (REAL32), DIMENSION(:,:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_vara_2d_real = nf90_put_var(ncid, varid, rvals,                   &
     &                                   start=start, count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] rvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_3d_real(ncid, varid, start, count, rvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                     :: nf_put_vara_3d_real
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      INTEGER, DIMENSION(:), INTENT(in)           :: start
      INTEGER, DIMENSION(:), INTENT(in)           :: count
      REAL (REAL32), DIMENSION(:,:,:), INTENT(in) :: rvals

!  Start of executable code
      nf_put_vara_3d_real = nf90_put_var(ncid, varid, rvals,                   &
     &                                   start=start, count=count)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a complex real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] cvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_cpx_real(ncid, varid, start, count, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                    :: nf_put_vara_cpx_real
      INTEGER, INTENT(in)                        :: ncid
      INTEGER, INTENT(in)                        :: varid
      INTEGER, DIMENSION(:), INTENT(in)          :: start
      INTEGER, DIMENSION(:), INTENT(in)          :: count
      COMPLEX (REAL32), DIMENSION(:), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:), ALLOCATABLE   :: rvals
      INTEGER                                    :: i

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals)*2))

      DO i = 1, SIZE(rvals), 2
         rvals(i) = REAL(cvals((i + 1)/2))
         rvals(i + 1) = AIMAG(cvals((i + 1)/2))
      END DO

      nf_put_vara_cpx_real = nf90_put_var(ncid, varid, rvals,                  &
     &                                    start=start, count=count)

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a complex real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] cvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_2d_cpx_real(ncid, varid, start, count, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_vara_2d_cpx_real
      INTEGER, INTENT(in)                          :: ncid
      INTEGER, INTENT(in)                          :: varid
      INTEGER, DIMENSION(:), INTENT(in)            :: start
      INTEGER, DIMENSION(:), INTENT(in)            :: count
      COMPLEX (REAL32), DIMENSION(:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:,:), ALLOCATABLE   :: rvals
      INTEGER                                      :: i
      INTEGER                                      :: j

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals, 1)*2, SIZE(cvals, 2)))

      DO j = 1, SIZE(rvals, 2)
         DO i = 1, SIZE(rvals, 1), 2
            rvals(i,j) = REAL(cvals((i + 1)/2, j))
            rvals(i + 1,j) = AIMAG(cvals((i + 1)/2, j))
         END DO
      END DO

      nf_put_vara_2d_cpx_real = nf90_put_var(ncid, varid, rvals,               &
     &                                       start=start, count=count)

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a complex real value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting indices.
!>  @param[in] count Number of elements to count.
!>  @param[in] cvals Real value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_3d_cpx_real(ncid, varid, start, count, cvals)
      USE iso_fortran_env, only: REAL32

      IMPLICIT NONE

!  Declare Arguments
      INTEGER :: nf_put_vara_3d_cpx_real
      INTEGER, INTENT(in)                            :: ncid
      INTEGER, INTENT(in)                            :: varid
      INTEGER, DIMENSION(:), INTENT(in)              :: start
      INTEGER, DIMENSION(:), INTENT(in)              :: count
      COMPLEX (REAL32), DIMENSION(:,:,:), INTENT(in) :: cvals

!  local variables
      REAL (REAL32), DIMENSION(:,:,:), ALLOCATABLE   :: rvals
      INTEGER                                        :: i
      INTEGER                                        :: j
      INTEGER                                        :: k

!  Start of executable code
      ALLOCATE(rvals(SIZE(cvals, 1)*2, SIZE(cvals, 2), SIZE(cvals, 3)))

      DO k = 1, SIZE(rvals, 3)
         DO j = 1, SIZE(rvals, 2)
            DO i = 1, SIZE(rvals, 1), 2
               rvals(i,j,k) = REAL(cvals((i + 1)/2,j,k))
               rvals(i + 1,j,k) = AIMAG(cvals((i + 1)/2,j,k))
            END DO
         END DO
      END DO

      nf_put_vara_3d_cpx_real = nf90_put_var(ncid, varid, rvals,               &
     &                                       start=start, count=count)

      DEALLOCATE(rvals)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a text value.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] text  String value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_var_text(ncid, varid, text)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                       :: nf_put_var_text
      INTEGER, INTENT(in)           :: ncid
      INTEGER, INTENT(in)           :: varid
      CHARACTER (len=*), INTENT(in) :: text

!  Start of executable code
      nf_put_var_text = nf90_put_var(ncid, varid, text)

      END FUNCTION

!-------------------------------------------------------------------------------
!>  @brief Write a text valued array.
!>
!>  @param[in] ncid  Netcdf file id.
!>  @param[in] varid Variable id.
!>  @param[in] start Starting index of the string array.
!>  @param[in] count Number of elements to count.
!>  @param[in] text  String value to write.
!-------------------------------------------------------------------------------
      FUNCTION nf_put_vara_text(ncid, varid, start, count, text)

      IMPLICIT NONE

!  Declare Arguments
      INTEGER                                     :: nf_put_vara_text
      INTEGER, INTENT(in)                         :: ncid
      INTEGER, INTENT(in)                         :: varid
      INTEGER, DIMENSION(:), INTENT(in)           :: start
      INTEGER, DIMENSION(:), INTENT(in)           :: count
      CHARACTER (len=*), DIMENSION(:), INTENT(in) :: text

!  Start of executable code
      nf_put_vara_text = nf90_put_var(ncid, varid, text, start=start,          &
     &                                count=count)

      END FUNCTION

      END MODULE
