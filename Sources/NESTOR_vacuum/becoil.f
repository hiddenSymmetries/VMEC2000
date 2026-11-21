      SUBROUTINE becoil (rad, zee, br, bp, bz, brvac, bpvac, bzvac,            &
     &                   lscreen)
      USE vparams, ONLY: nthreed
      USE vacmod
      USE parallel_include_module
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      REAL(dp), DIMENSION(nuv3), INTENT(in) :: rad, zee
      REAL(dp), DIMENSION(nuv3), INTENT(out) :: br, bp, bz
      REAL(dp), DIMENSION(nr0b,nz0b,np0b), INTENT(in) ::
     1   brvac, bpvac, bzvac
      LOGICAL :: lscreen
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      CHARACTER(LEN=50), PARAMETER :: warning =
     1   'Plasma Boundary exceeded Vacuum Grid Size'
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER, SAVE :: icount = 0
      INTEGER :: i, kv, ir, jz, ii, jj, nrw, nzw
      REAL(dp) :: rad0, zee0, ri, zj, pr, qz, w, tbecon, tbecoff
      REAL(dp), DIMENSION(4) :: wr, wz
C-----------------------------------------------
!
!     DETERMINE THE CYLINDRICAL COMPONENTS OF THE EXTERNAL
!     MAGNETIC FIELD (BR, BP, BZ) AT A FIXED PHI PLANE BY 
!     USING 2-D TENSOR-PRODUCT 4-POINT (CUBIC) LAGRANGE INTERPOLATION
!     ON A CENTRED STENCIL, SHIFTED INWARD AT THE GRID EDGES.
!     THE FORMER BILINEAR FORMULA (ABRAMOWITZ AND STEGUN 25.2.66) LEFT
!     AN O(H**2) FIELD ERROR THAT MOVED THE FREE BOUNDARY.
!
!     BRVAC, BPVAC, BZVAC:       CYLINDRICAL COMPONENTS OF VACUUM B-FIELD 
!                                STORED ON R, Z, PHI GRID
!
!     RAD, ZEE:                  ARRAY OF R, Z VALUES (ON FLUX SURFACE)
!                                AT WHICH B-FIELD IS TO BE DETERMINED
!
C-----------------------------------------------
      CALL second0(tbecon)

      icount = icount + 1

      DO i = nuv3min, nuv3max
!
!       CHECK THAT BOUNDARY POINTS ARE INSIDE VACUUM GRID.  IF NOT,
!       SET THEM EQUAL TO LARGEST (OR SMALLEST) VACUUM GRID POINTS
!
         rad0 = MIN(rad(i), rmaxb)
         rad0 = MAX(rad0,   rminb)
         zee0 = MIN(zee(i), zmaxb)
         zee0 = MAX(zee0,   zminb)
!
!       DETERMINE PHI-PLANE, KV (MUST LIE IN FIRST FIELD PERIOD)
!
         kv = 1 + MOD(i - 1,nv)
         kv = MIN (kv, np0b)                       !!Axi-symmetric special case
!
!
!       DETERMINE THE FIRST STENCIL NODE (IR,JZ) AND THE OFFSETS PR, QZ
!       OF THE POINT FROM IT IN UNITS OF THE GRID SPACING
!
         nrw = MIN(4, nr0b)
         nzw = MIN(4, nz0b)
         ri = (rad0 - rminb)/delrb
         zj = (zee0 - zminb)/delzb
         ir = MIN(MAX(INT(ri) - 1, 0), nr0b - nrw) + 1
         jz = MIN(MAX(INT(zj) - 1, 0), nz0b - nzw) + 1
         pr = ri - (ir - 1)
         qz = zj - (jz - 1)
         CALL lagrange_weights (pr, nrw, wr)
         CALL lagrange_weights (qz, nzw, wz)
!
!       COMPUTE B FIELD AT R, PHI, Z BY INTERPOLATION
!
         br(i) = 0; bp(i) = 0; bz(i) = 0
         DO jj = 1, nzw
            DO ii = 1, nrw
               w = wr(ii)*wz(jj)
               br(i) = br(i) + w*brvac(ir+ii-1,jz+jj-1,kv)
               bp(i) = bp(i) + w*bpvac(ir+ii-1,jz+jj-1,kv)
               bz(i) = bz(i) + w*bzvac(ir+ii-1,jz+jj-1,kv)
            END DO
         END DO

      END DO

!
!     PRINT INFO IF R, Z OUT OF BOX
!
      IF (MOD(icount,25).EQ.0 .AND. rank.EQ.0) THEN
         i = 0
         rad0 = MAXVAL(rad)
         zee0 = MAXVAL(zee)
         IF (rad0 .gt. rmaxb) i = 1
         IF (zee0 .gt. zmaxb) i = i + 2
         ri   = MINVAL(rad)
         zj   = MINVAL(zee)
         IF (ri .lt. rminb) i = i + 4
         IF (zj .lt. zminb) i = i + 8
         IF (i .ne. 0 .and. lscreen) THEN
            PRINT *, warning
            WRITE (nthreed, *) warning
            IF (i/8 .ne. 0) PRINT *,' zmin = ', zj
            i = MOD(i,8)
            IF (i/4 .ne. 0) PRINT *,' rmin = ', ri
            i = MOD(i,4)
            IF (i/2 .ne. 0) PRINT *,' zmax = ', zee0
            i = MOD(i,2)
            IF (i .ne. 0) PRINT *,' rmax = ', rad0
         END IF
      ENDIF

      CALL second0(tbecoff)
      becoil_time = becoil_time + (tbecoff - tbecon)

      CONTAINS

      SUBROUTINE lagrange_weights (u, n, wt)
!     WEIGHTS OF THE N-NODE LAGRANGE POLYNOMIAL ON NODES 0, ..., N-1 AT U
      REAL(dp), INTENT(in) :: u
      INTEGER, INTENT(in) :: n
      REAL(dp), DIMENSION(4), INTENT(out) :: wt
      INTEGER :: k, m

      wt = 0
      DO k = 1, n
         wt(k) = 1
         DO m = 1, n
            IF (m .ne. k) wt(k) = wt(k)*(u - (m - 1))/(k - m)
         END DO
      END DO

      END SUBROUTINE lagrange_weights

      END SUBROUTINE becoil
