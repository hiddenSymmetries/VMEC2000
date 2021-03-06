      SUBROUTINE analysum2 (grpmn, bvec, slp, tlp, slm, tlm,
     1    m, n, l, ivacskip, ndim)
      USE vacmod
      USE parallel_include_module
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER, INTENT(in) :: m, n, l, ivacskip, ndim
      REAL(dp), INTENT(inout) :: grpmn(0:mf,-nf:nf,ndim,nuv3)
      REAL(dp), INTENT(inout) :: bvec(0:mf,-nf:nf,ndim)
      REAL(dp), DIMENSION(nuv3), INTENT(in) ::
     1   slp, tlp, slm, tlm
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i
      REAL(dp) :: sinp, sinm, cosp, cosm, temp, ton, toff
C-----------------------------------------------
      CALL second0(ton)
    
      IF (n .LT. 0) STOP 'error calling analysum2!'

      DO i = nuv3min, nuv3max
         sinp =  sinu1(i,m)*cosv1(i,n)*cmns(l,m,n)
         temp = -cosu1(i,m)*sinv1(i,n)*cmns(l,m,n)
         sinm = sinp - temp                 !SIN(mu + |n|v) * cmns (l,m,|n|)
         sinp = sinp + temp                 !SIN(mu - |n|v) * cmns (l,m,|n|)
         bvec(m,n,1)  = bvec(m,n,1)  + tlp(i)*bexni(i)*sinp
         bvec(m,-n,1) = bvec(m,-n,1) + tlm(i)*bexni(i)*sinm

         IF (ivacskip .EQ. 0) THEN
            grpmn(m,n,1,i)  = grpmn(m,n,1,i)  + slp(i)*sinp
            grpmn(m,-n,1,i) = grpmn(m,-n,1,i) + slm(i)*sinm
         END IF


         IF (lasym) THEN
         cosp = cosu1(i,m)*cosv1(i,n)*cmns(l,m,n)
         temp = sinu1(i,m)*sinv1(i,n)*cmns(l,m,n)
         cosm = cosp - temp                !COS(mu + |n|v) * cmns (l,m,|n|)
         cosp = cosp + temp                !COS(mu - |n|v) * cmns (l,m,|n|)
         bvec(m,n,2)  = bvec(m,n,2)  + tlp(i)*bexni(i)*cosp
         bvec(m,-n,2) = bvec(m,-n,2) + tlm(i)*bexni(i)*cosm

         IF (ivacskip .EQ. 0) THEN
            grpmn(m,n,2,i)  = grpmn(m,n,2,i)  + slp(i)*cosp
            grpmn(m,-n,2,i) = grpmn(m,-n,2,i) + slm(i)*cosm
         END IF
         END IF
      END DO

      CALL second0(toff)
      timer_vac(tasum2) = timer_vac(tasum2) + (toff-ton)

      END SUBROUTINE analysum2
