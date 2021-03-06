      SUBROUTINE profil3d_par(rmn, zmn, lreset, linterp)
      USE vmec_main
      USE vmec_params
      USE vspline, ONLY: sknots, pknots, hstark, hthom
      USE realspace
      USE xstuff
#ifdef _HBANGLE
      USE angle_constraints, ONLY: store_init_array
#endif
      USE parallel_include_module
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(0:ntor,0:mpol1,ns,ntmax), INTENT(inout) ::
     &   rmn, zmn
      LOGICAL, INTENT(in) :: lreset, linterp

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
      REAL(dp), DIMENSION(0:ntor,ntmax) :: rold, zold
      REAL(dp) :: sm0, t1, facj, si, rax1, zax1
      INTEGER :: jcount, jk, k
      INTEGER :: i, j, nsmin, nsmax, lpar
      REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: bcast_buf

!-----------------------------------------------
!                INDEX OF LOCAL VARIABLES
!
!     phip     radial derivative of phi/(2*pi) on half-grid
!     chip     radial derivative of chi/(2*pi) on half-grid
!     shalf    sqrt(s) ,two-dimensional array on half-grid
!     sqrts    sqrt(s), two-dimensional array on full-grid
!     wint     two-dimensional array for normalizing angle integrations
!     ireflect two-dimensional array for computing 2pi-v angle
!-----------------------------------------------
!      CALL second0(tprofon)

      nsmin = t1lglob; nsmax = t1rglob
      DO js = nsmin, nsmax
         pphip(:,js) = phips(js)
         pchip(:,js) = chips(js)
      END DO

      faclam = 0
      sigma_an = 1         !INITIALIZE sigma FOR ANISOTROPIC PLASMA
      pwint(:,1) = 0
      pfaclam(:,:,nsmin:nsmax,:)=0
!
!     COMPUTE ARRAY FOR REFLECTING v = -v (ONLY needed for lasym)
!
      DO k = 1, nzeta
         jk = nzeta + 2 - k
         IF (k .eq. 1) THEN
            jk = 1
         END IF
         ireflect_par(k) = jk
      END DO

      lk = 0
      DO lt = 1, ntheta3
         DO lz = 1, nzeta
            lk = lk + 1
            pwint_ns(lk) = cosmui3(lt,0)/mscale(0)
            DO js = MAX(2, t1lglob), t1rglob
               pwint(lk,js) = pwint_ns(lk)
            END DO
         END DO
      END DO

!     INDEX FOR u = -u (need for lasym integration in wrout)
      lk = 0
      IF (.NOT.ALLOCATED(uminus)) THEN
         ALLOCATE(uminus(nznt))
      END IF
      DO lt = 1, ntheta2
         k = ntheta1 - lt + 2
         IF (lt .eq. 1) THEN
            k = 1             !u=-0 => u=0
         END IF
         DO lz = 1, nzeta
            lk = lk + 1
            uminus(lk) = k                !(-u), for u = 0,pi
         END DO
      END DO

!
!     COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
!     FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
!     (1/SQRTS FACTOR FOR ODD M VALUES)
!

      nsmin = t1lglob
      nsmax = t1rglob

      rold(0:ntor,1:ntmax) = rmn(0:ntor,0,1,1:ntmax)
      zold(0:ntor,1:ntmax) = zmn(0:ntor,0,1,1:ntmax)

      IF (nranks.GT.1) THEN
         ALLOCATE(bcast_buf(0:2*ntor+1,1:ntmax))
         bcast_buf(0:ntor,1:ntmax)          = rold(0:ntor,1:ntmax)
         bcast_buf(ntor+1:2*ntor+1,1:ntmax) = zold(0:ntor,1:ntmax)
         CALL MPI_Bcast(bcast_buf, 2*(ntor + 1)*ntmax, MPI_REAL8, 0,
     &                  NS_COMM, MPI_ERR)
         rold(0:ntor,1:ntmax) = bcast_buf(0:ntor,1:ntmax)
         zold(0:ntor,1:ntmax) = bcast_buf(ntor+1:2*ntor+1,1:ntmax)
         DEALLOCATE(bcast_buf)
      END IF

      nsmin = t1lglob
      nsmax = t1rglob
      DO js = nsmin, nsmax
         si = psqrts(1,js)*psqrts(1,js)
         sm0 = one - si
         DO ntype = 1, ntmax
            DO m = 0, mpol1
               DO n = 0, ntor
                  t1 = one/(mscale(m)*nscale(n))
                  mn = n + ntor1*m
                  lpar = mn+mnsize*(js - 1) + (ntype - 1)*mns + 1
                  IF (MOD(m,2) .eq. 0) THEN
                     pscalxc(lpar) = one
                  ELSE
                     pscalxc(lpar) = one/psqrts(1,MAX(2,js))
                  END IF

                  pscalxc(lpar+irzloff)=pscalxc(lpar)
                  pscalxc(lpar+2*irzloff)=pscalxc(lpar)

!             Do not overwrite r,z if read in from wout file AND in free bdy mode
!             For fixed boundary, edge values MAY have been perturbed, so must execute this loop
                  IF (.not.lreset .and. lfreeb) CYCLE
                  IF (m .eq. 0) THEN
                     IF (.not.lreset) CYCLE        !Freeze axis if read in from wout file

                     rmn(n,m,js,ntype) = rmn(n,m,js,ntype)
     &                                 + si*(rmn_bdy(n,m,ntype)*t1 -
     &                                       rmn(n,m,ns,ntype))
                     zmn(n,m,js,ntype) = zmn(n,m,js,ntype)
     &                                 + si*(zmn_bdy(n,m,ntype)*t1 -
     &                                       zmn(n,m,ns,ntype))

                     IF (ntype .eq. rcc) rax1 = raxis_cc(n)
                     IF (ntype .eq. zcs) zax1 =-zaxis_cs(n)
                     IF (ntype .eq. rcs) rax1 =-raxis_cs(n)
                     IF (ntype .eq. zcc) zax1 = zaxis_cc(n)

                     IF (ntype.eq.rcc .or. ntype.eq.rcs) THEN
                        rmn(n,m,js,ntype) = rmn(n,m,js,ntype)
     &                                    + sm0*(rax1*t1 -
     &                                           rold(n,ntype))
                     END IF
                     IF (ntype.eq.zcs .or. ntype.eq.zcc) THEN
                        zmn(n,m,js,ntype) = zmn(n,m,js,ntype)
     &                                    + sm0*(zax1*t1 -
     &                                           zold(n,ntype))
                     END IF
                  ELSE
                     facj = psqrts(1,js)**m        !!TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
                     rmn(n,m,js,ntype) = rmn(n,m,js,ntype)
     &                                 + (rmn_bdy(n,m,ntype)*t1 -
     &                                    rmn(n,m,ns,ntype))*facj
                     zmn(n,m,js,ntype) = zmn(n,m,js,ntype)
     &                                 + (zmn_bdy(n,m,ntype)*t1 -
     &                                    zmn(n,m,ns,ntype))*facj
                  END IF
               END DO
            END DO
         END DO
      END DO

#ifdef _HBANGLE
      IF (.NOT.linterp) THEN
         CALL store_init_array(xc)
      END IF
#endif

      END SUBROUTINE profil3d_par

      SUBROUTINE profil3d(rmn, zmn, lreset, linterp)
      USE vmec_main
      USE vmec_params
      USE vspline, ONLY: sknots, pknots, hstark, hthom
      USE realspace
      USE xstuff
#ifdef _HBANGLE
      USE angle_constraints, ONLY: store_init_array
#endif

      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) ::
     &   rmn, zmn
      LOGICAL, INTENT(in) :: lreset, linterp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: js, l, lk, lt, lz, ntype, m, n, mn
      REAL(dp), DIMENSION(0:ntor,ntmax) :: rold, zold
      REAL(dp) :: sm0, t1, facj, si, rax1, zax1
      INTEGER :: jcount, jk, k

!-----------------------------------------------

!
!                INDEX OF LOCAL VARIABLES
!
!        phip     radial derivative of phi/(2*pi) on half-grid
!        chip     radial derivative of chi/(2*pi) on half-grid
!        shalf    sqrt(s) ,two-dimensional array on half-grid
!        sqrts    sqrt(s), two-dimensional array on full-grid
!        wint     two-dimensional array for normalizing angle integrations
!        ireflect two-dimensional array for computing 2pi-v angle

      DO js = 1, ns
         phip(js:nrzt:ns) = phips(js)
         chip(js:nrzt:ns) = chips(js)
      END DO

      phip(nrzt + 1) = 0
      faclam = 0
      sigma_an = 1         !INITIALIZE sigma FOR ANISOTROPIC PLASMA
      wint(1:nrzt:ns) = 0

      lk = 0
      DO lt = 1, ntheta3
         DO lz = 1, nzeta
            lk = lk + 1
            DO js = 2,ns
               wint(js+ns*(lk-1)) = cosmui3(lt,0)/mscale(0)
            END DO
         END DO
      END DO

!
!     COMPUTE ARRAY FOR REFLECTING v = -v (ONLY needed for lasym)
!
      jcount = 0
      DO k = 1, nzeta
         jk = nzeta + 2 - k
         IF (k .eq. 1) jk = 1
         DO js = 1, ns
            jcount = jcount+1
            ireflect(jcount) = js + ns*(jk - 1)           !Index for -zeta[k]
         END DO
      END DO

!     INDEX FOR u = -u (need for lasym integration in wrout)
      lk = 0
      IF (.NOT.ALLOCATED(uminus)) THEN
         ALLOCATE (uminus(nznt))
      END IF
      DO lt = 1, ntheta2
         k = ntheta1 - lt + 2
         IF (lt .eq. 1) THEN
            k = 1             !u=-0 => u=0
         END IF
         DO lz = 1, nzeta
            lk = lk + 1
            uminus(lk) = k                !(-u), for u = 0,pi
         END DO
      END DO

!
!     COMPUTE INITIAL R AND Z FOURIER COEFFICIENTS,
!     FROM SCALED BOUNDARY VALUES, AND SCALXC ARRAY
!     (1/SQRTS FACTOR FOR ODD M VALUES)
!

      DO js = 1, ns
         si = sqrts(js)*sqrts(js)
         sm0 = one - si
         DO ntype = 1, ntmax
            DO m = 0, mpol1
               DO n = 0, ntor
                  t1 = one/(mscale(m)*nscale(n))
                  mn = n + ntor1*m
                  l = js + ns*mn + (ntype - 1)*mns
                  IF (MOD(m,2) .eq. 0) THEN
                     scalxc(l) = one
                  ELSE
                     scalxc(l) = one/MAX(sqrts(js),sqrts(2))
                  END IF
!                 Do not overwrite r,z if read in from wout file AND in free bdy mode
!                 For fixed boundary, edge values MAY have been perturbed, so must execute this loop
                  IF (.not.lreset .and. lfreeb) CYCLE
                  IF (m .eq. 0) THEN
                     IF (.not.lreset) CYCLE        !Freeze axis if read in from wout file
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype)
     &                                    + si*(rmn_bdy(n,m,ntype)*t1 -
     &                                          rmn(ns,n,m,ntype))
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype)
     &                                    + si*(zmn_bdy(n,m,ntype)*t1 -
     &                                          zmn(ns,n,m,ntype))
                     IF (js .eq. 1) THEN
                        rold(n,ntype) = rmn(1,n,0,ntype)
                        zold(n,ntype) = zmn(1,n,0,ntype)
                     END IF
                     IF (ntype .eq. rcc) rax1 = raxis_cc(n)
                     IF (ntype .eq. zcs) zax1 =-zaxis_cs(n)
                     IF (ntype .eq. rcs) rax1 =-raxis_cs(n)
                     IF (ntype .eq. zcc) zax1 = zaxis_cc(n)
                     IF (ntype.eq.rcc .or. ntype.eq.rcs) THEN
                        rmn(js,n,m,ntype) = rmn(js,n,m,ntype)
     &                                    + sm0*(rax1*t1 -
     &                                           rold(n,ntype))
                     END IF
                     IF (ntype.eq.zcs .or. ntype.eq.zcc) THEN
                        zmn(js,n,m,ntype) = zmn(js,n,m,ntype)
     &                                    + sm0*(zax1*t1 -
     &                                           zold(n,ntype))
                     END IF
                  ELSE
                     facj = sqrts(js)**m        !!TURN OFF NEXT 3 LINES IF THIS ONE ACTIVATED
                     rmn(js,n,m,ntype) = rmn(js,n,m,ntype)
     &                                 + (rmn_bdy(n,m,ntype)*t1 -
     &                                    rmn(ns,n,m,ntype))*facj
                     zmn(js,n,m,ntype) = zmn(js,n,m,ntype)
     &                                 + (zmn_bdy(n,m,ntype)*t1 -
     &                                    zmn(ns,n,m,ntype))*facj
                  END IF
               END DO
            END DO
         END DO
      END DO

      scalxc(1+irzloff:2*irzloff)   = scalxc(:irzloff)              !Z-components
      scalxc(1+2*irzloff:3*irzloff) = scalxc(:irzloff)              !Lamda-components

#ifdef _HBANGLE
      IF (.NOT.linterp) CALL store_init_array(xc)
#endif
      
      END SUBROUTINE profil3d
