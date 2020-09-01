      SUBROUTINE getcarg(narg, arg, numargs)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)  :: narg
      INTEGER, INTENT(out) :: numargs
      CHARACTER(LEN=*), INTENT(out) :: arg
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: numchars
!-----------------------------------------------
#if defined (LINUX)
      INTEGER iargc
      numargs = iargc()
      CALL getarg(narg, arg)
#elif defined (CRAY)
      INTEGER :: ier
      INTEGER ipxfargc
      numargs = ipxfargc()
      CALL pxfgetarg(narg, arg, numchars, ier)
#elif defined (MACOSX)
      numargs = command_argument_count()
      CALL get_command_argument(narg,arg)
#else
      INTEGER iargc, getarg
      numargs = iargc()
      numchars = getarg(narg, arg)
#endif
      END SUBROUTINE getcarg
