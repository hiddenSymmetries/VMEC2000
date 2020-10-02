      MODULE vmec_ext_interface
        USE, INTRINSIC :: ISO_C_BINDING
        USE stel_kinds
        USE stel_constants
        IMPLICIT NONE
        
        PUBLIC :: vmec_output_data

        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rmnc, rmns, zmns, 
     &   zmnc, lmns, lmnc
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: gmnc, bmnc,
     &   gmns, bmns,
     &   bsubumnc, bsubvmnc, bsubsmns, bsubumns, bsubvmns, bsubsmnc,
     &   currumnc, currvmnc, currumns, currvmns
        REAL(dp), DIMENSION(:), ALLOCATABLE :: gmn, bmn,
     &   bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn

        REAL(dp), DIMENSION(:), POINTER ::   xm_nyq0, xn_nyq0
        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: bsupumnc, bsupumns,
     &   bsupvmnc, bsupvmns

        LOGICAL :: lcurr

        INTERFACE set_vmec_data
          MODULE PROCEDURE set_vmec_data_real
          MODULE PROCEDURE set_vmec_data_int
          MODULE PROCEDURE set_vmec_data_bool
          MODULE PROCEDURE set_vmec_data_char
        END INTERFACE

        INTERFACE retrieve_vmec_data
          MODULE PROCEDURE set_vmec_data_real
          MODULE PROCEDURE set_vmec_data_int
          MODULE PROCEDURE set_vmec_data_bool
        END INTERFACE

      CONTAINS

      SUBROUTINE runvmec_ext(input_file0,EXT_MPI_COMM) 
!     & BIND(C,name='runvmec_c')
        USE vmec_input
        USE vmec_seq
        USE safe_open_mod
        USE vparams, ONLY: nlog, nlog0, nthreed
        USE vmec_params, ONLY: more_iter_flag,
     &                       bad_jacobian_flag,
     &    restart_flag, readin_flag, timestep_flag,
     &    output_flag, cleanup_flag,
     &    norm_term_flag, successful_term_flag ! J Geiger: for more iterations and full 3D1-output
        USE parallel_include_module, ONLY: grank, 
     &                                   MPI_ERR
        USE parallel_vmec_module, ONLY: MyEnvVariables, 
     &                                InitializeParallel,
     &                                FinalizeParallel
        IMPLICIT NONE
        TYPE(C_PTR), INTENT(IN) :: input_file0
        INTEGER, INTENT(IN) :: EXT_MPI_COMM
        CHARACTER(LEN=:), ALLOCATABLE :: file_string

C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
        INTEGER, PARAMETER :: nseq0 = 12
        CHARACTER(LEN=*), PARAMETER ::
     &    increase_niter = "Try increasing NITER",
     &    bad_jacobian = "The jacobian was non-definite!",
     &    full_3d1output_request = "Full threed1-output request!"
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
        INTEGER :: numargs, ierr_vmec, index_end,
     &   iopen, isnml, iread, iseq, index_seq,
     &   index_dat, iunit, ncount, nsteps, i
        INTEGER :: ictrl(5)
        CHARACTER(LEN=120) :: input_file, seq_ext, reset_file_name, arg
        CHARACTER(LEN=120) :: log_file
        CHARACTER(LEN=120), DIMENSION(10) :: command_arg
        LOGICAL :: lscreen

C-----------------------------------------------
!***
!                              D   I   S   C   L   A   I   M   E   R
!
!       You are using a beta version of the PROGRAM VMEC, which is currently
!       under development by S. P. Hirshman at the Fusion Energy Division,
!       Oak Ridge National Laboratory.  Please report any problems or comments
!       to him.  As a beta version, this program is subject to change
!       and improvement without notice.
!
!       1. CODE SYNOPSIS
!
!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see SUBROUTINE BECOIL)
!
!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
!
!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
!
!                  iota(s) * grad(v) X grad(phiT)
!
!       WHERE phiT is the toroidal flux (called phi in code) and
!       u,v are the poloidal, toroidal angles, respectively.
!
!       2. ADDITIONAL CODES REQUIRED
!       For the fixed boundary calculation, the user must provide the Fourier
!       coefficients for the plasma boundary (the last surface outside of which
!       the pressure gradient vanishes). For ALL but the simplest geometry, the
!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
!       code, can be used to produce the optimized VMEC Fourier representation for
!       an arbritrary closed boundary (it need not be a 'star-like' DOmain, nor
!       need it possess vertical, or 'stellarator', symmetry).
!
!       For the free boundary calculation, the MAKEGRID code (available upon
!       request) is needed to create a binary Green''s FUNCTION table for the
!       vacuum magnetic field(s) and, IF data analysis is to be done, flux and
!       field loops as well. The user provides a SUBROUTINE (BFIELD) which can be
!       called at an arbitrary spatial location and which should RETURN the three
!       cylindrical components of the vacuum field at that point. (Similary,
!       locations of diagnostic flux loops, Rogowski coils, etc. are required IF
!       equilibrium reconstruction is to be done.)
!
!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
!       file, WOUT.EXT, WHERE 'EXT' is the command-line extension of the INPUT file.
!
!
!       3. UNIX SCRIPT SETUP PARAMETERS
!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
!       the C-precompiler to produce both the machine-specific Fortran source and a
!       make-file specific to ANY one of the following platforms:
!
!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
!       WINDOWS-NT, DEC-VMS
!
!       Additional platforms are easy to add to the existing script as required.
!
!
!       4. FORTRAN PARAMETER STATEMENTS set by user
!       In the Fortran-90 version of VMEC these PARAMETER statements have
!       been replaced by dynamic memory allocation. So the user should set the
!       run-time parameters ns (through ns_array), mpol, ntor in the NAMELIST INDATA.
!
!
!       Added features since last edition (see vmec_params for revision history list)
!       1. Implemented preconditioning algorithm for R,Z
!       2. The physical (unpreconditioned) residuals are used
!          to determine the level of convergence
!       3. The original (MOMCON) scaling of lambda is used, i.e.,
!          Bsupu = phip*(iota - lamda[sub]v)/SQRT(g). This is needed to
!          maintain consistency with the time-stepper for arbitrary PHIP.
!
!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
!

!     Local variables
!
!     ictrl:   array(5) of control variables for running "runvmec" routine
!              see "runvmec" for a description
!

!
!     Read in command-line arguments to get input file or sequence file,
!     screen display information, and restart information
        INTERFACE
           SUBROUTINE runvmec(ictrl_array, input_file0, 
     &                      lscreen, EXT_MPI_COMM, reset_file_name)
             IMPLICIT NONE
             INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
             LOGICAL, INTENT(in) :: lscreen
             CHARACTER(LEN=*), INTENT(in) :: input_file0
             INTEGER, INTENT(in), OPTIONAL :: EXT_MPI_COMM
             CHARACTER(LEN=*), OPTIONAL :: reset_file_name
           END SUBROUTINE runvmec
        END INTERFACE

        CALL MyEnvVariables
        CALL InitializeParallel(EXT_MPI_COMM)

        lscreen = .false.
        IF(grank.EQ.0) lscreen = .true.
        reset_file_name = " "

!
!     Determine type of file opened (sequential or input-data)
!     ARG1 (char var)
!          By DEFAULT, ARG1 obtained from the command
!          line is parsed as follows to determine the input data file(s):
!               a. Attempt to OPEN file ARG1 (full path + file name).
!                  Look for the VSEQ NAMELIST to obtain nseq, nseq_select, and
!                  extension array. If they exist and nseq>0, VMEC will run
!                  sequentially using input determined from the array EXTENSION[i]
!                  or input.EXTENSION[i]
!               b. If the command argument is not a sequence NAMELIST, THEN the data file
!                  ARG1 or input.ARG1 is READ directly, with NSEQ=1.
!
        CALL c2f_string_array_1d(input_file0,file_string)
        arg = TRIM(file_string)
        index_dat = INDEX(arg,'.')
        index_end = LEN_TRIM(arg)
        IF (index_dat .gt. 0) THEN
           seq_ext  = arg(index_dat + 1:index_end)
           input_file = TRIM(arg)
        END IF
!
        nseq = 1
        nseq_select(1) = 1
        extension(1) = input_file
!
!        iunit = nseq0
!        IF (iopen .eq. 0) THEN
!          CALL read_namelist (iunit, isnml, 'indata')
!
!        END IF
!        CLOSE (iunit)

!
!     CALL EQUILIBRIUM SOLVER
!
!     nseq_select:      If sequence file (VSEQ NAMELIST given with nseq >0)
!                       array giving indices into EXTENSION array prescribing
!                       the order in which the input files are run by VMEC
!     nseq:             number of sequential VMEC runs to make
!
!
!     CALL VMEC WITH POSSIBLE SEQUENCE EXTENSION (SEQ_EXT)
!     AND ARRAY OF INPUT FILE EXTENSIONS (EXTENSION)
!
        ictrl = 0


        iseq = 1
        index_seq = nseq_select(iseq)
        ictrl(1) = restart_flag + readin_flag + timestep_flag
     &            + output_flag + cleanup_flag                !Sets all flags
        ictrl(2) = 0
!         ictrl(3) = 100
!         ictrl(4) = 2
        ictrl(5) = iseq - 1
        ncount = 0

        print *, ictrl
        CALL runvmec(ictrl, extension(index_seq), lscreen, EXT_MPI_COMM,
     &                reset_file_name)

                print *, ictrl
        ierr_vmec = ictrl(2)
        print *, ierr_vmec

        SELECT CASE (ierr_vmec)
          CASE (more_iter_flag)                                !Need a few more iterations to converge
            IF (grank .EQ. 0) THEN
              IF(lscreen) WRITE (6, '(1x,a)') increase_niter
              WRITE (nthreed, '(1x,a)') increase_niter
              WRITE (nthreed, '(1x,a)') "PARVMEC aborting..."
              CALL FLUSH(nthreed)
            END IF
! J Geiger: if lmoreiter and lfull3d1out are false
!           the o-lines (original) are the only
!           ones to be executed.
            IF (lmoreiter) THEN                                 ! J Geiger: --start--
              DO i = 2, max_main_iterations                    ! Changes to run
                ictrl(1) = timestep_flag                      ! some more iterations if requested
                ictrl(3) = niter                              ! - this is the number of iterations
                CALL runvmec(ictrl, extension(1), lscreen,
     &                            EXT_MPI_COMM, reset_file_name)       ! - the second iteration run with ictrl(3) iterations
                IF (ictrl(2) .EQ. more_iter_flag .and.
     &                   grank    .EQ. 0) THEN
                  WRITE (nthreed, '(1x,a)') increase_niter
                  IF(lscreen) WRITE (6, '(1x,a)') increase_niter
                END IF
              END DO
              ictrl(1) = output_flag + cleanup_flag            ! - Output, cleanup
              IF (ictrl(2) .ne. successful_term_flag) THEN
                ictrl(2)=successful_term_flag                 ! - force success flag to get full threed1-output!
              END IF
              ictrl(3) = 0                                     ! - this is the number of iterations
              CALL runvmec(ictrl, extension(1), lscreen, EXT_MPI_COMM,
     &                         reset_file_name)
            ELSE                                                ! else-branch contains original code.
#if defined(MPI_OPT)
              CALL MPI_Barrier(EXT_MPI_COMM, MPI_ERR)
#endif

              ictrl(1) = output_flag + cleanup_flag            !Output, cleanup  ! o-lines
              ictrl(2) = 0                                     ! o-lines
!              IF (lfull3d1out) THEN
!                ictrl(2) = successful_term_flag
!                IF (grank .EQ. 0) THEN
!                  WRITE(6,'(1x,a)') full_3d1output_request
!                  WRITE(nthreed,'(1x,a)') full_3d1output_request
!                END IF
!              END IF

              CALL runvmec(ictrl, extension(1), lscreen, EXT_MPI_COMM, ! o-lines
     &                         reset_file_name)
            END IF                                              ! J Geiger: -- end --

          CASE (bad_jacobian_flag)                               !Bad jacobian even after axis reset and ns->3
            IF (grank .EQ. 0) THEN
              IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
                WRITE (nthreed, '(/,1x,a)') bad_jacobian
              END IF
          CASE DEFAULT
        END SELECT


        CALL FinalizeParallel

      END SUBROUTINE runvmec_ext

      SUBROUTINE c2f_string_array_1d(c_pointer,f_string)
        USE, INTRINSIC :: ISO_C_BINDING
        TYPE(C_PTR), INTENT(IN) :: c_pointer
        CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: f_string
        CHARACTER, DIMENSION(:), POINTER :: f_pointer
        CHARACTER(LEN=:), POINTER :: f_char_pointer

        INTEGER :: i, j
        LOGICAL :: null_char_found
        i = 0
        null_char_found = .FALSE.
        DO WHILE (null_char_found .EQV. .FALSE.)
          i = i + 1
          CALL C_F_POINTER(c_pointer,f_pointer,[i])
          IF (f_pointer(i) == C_NULL_CHAR) THEN
            null_char_found = .TRUE.
          END IF
        END DO
        i = i - 1
        
        ALLOCATE(CHARACTER(LEN=i)::f_string)
        CALL C_F_POINTER(C_LOC(f_pointer),f_char_pointer)

        f_string = f_char_pointer(1:i)

      END SUBROUTINE c2f_string_array_1d
    

      SUBROUTINE vmec_output_data(bsq, gsqrt,
     &  bsubu, bsubv, bsubs, bsupv, bsupu,rzl_array, gc_array)
        USE vmec_input, ONLY: ns_array, ftol_array, lwouttxt
        USE vmec_params
        USE vmec_main
        USE vmercier
        USE vmec_persistent
        USE vparams, p5 => cp5, two => c2p0
        USE vac_persistent
        USE vspline
        USE xstuff
        USE vmec_io
        USE realspace, ONLY: phip, chip, gsqrta=>z1, z1=>z1
        USE totzsp_mod
        USE vforces, ONLY: bsupua=>brmn_e, bsupva=>czmn_o, 
     &                   bsqa=>bzmn_e, bsubsa=>armn_e,
     &                   bsubua=>azmn_e, bsubva=>armn_o
        USE vacmod, ONLY: potvac, mnpd, xmpot, xnpot       !added for diagno, J.Geiger
        USE read_wout_mod, ONLY: Compute_Currents
        USE mgrid_mod

        IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
        REAL(dp), DIMENSION(mnmax,ns,3*MAX(ntmax/2,1)),           !reverse ns, mnmax for backwards compatibility
     &   INTENT(inout), TARGET :: rzl_array, gc_array
        REAL(dp), DIMENSION(ns,nznt), INTENT(inout) ::
     &   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu
        REAL(dp) :: qfact(ns)
        REAL(dp), PARAMETER :: c1p5 = 1.5_dp
        LOGICAL :: lnyquist = .TRUE.                               !=false, suppress nyquist stuff
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
        INTEGER :: j, js, jlk, mn, lk, iasym,
     &           m, n, k, n1, istat, i, indx1(1),
     &           mnmax_nyq0, mnyq0, nnyq0
     &          ,isgn, js2, nfort      !for diagno 1.5
        REAL(dp) :: dmult, tcosi, tsini, vversion, sgn, tmult, 
     &            presfactor, ftolx1, d_bsupumn, d_bsupvmn   ! diagno 1.5
        REAL(dp), DIMENSION(mnmax) :: rmnc1, zmns1, lmns1,
     &   rmns1, zmnc1, lmnc1, bmodmn, bmodmn1
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubumnc_sur  !MRC 10-15-15
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubvmnc_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupumnc_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupvmnc_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubumns_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubvmns_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupumns_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupvmns_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubua_sur, bsubva_sur
        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupua_sur, bsupva_sur

        REAL(dp), DIMENSION(:), ALLOCATABLE :: xfinal

!
!     THIS SUBROUTINE CREATES THE FILE WOUT.IT CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
!     COEFFICIENTS RMN,ZMN (full), LMN (half_mesh - CONVERTED FROM
!     INTERNAL full REPRESENTATION), AS WELL AS COEFFICIENTS (ON NYQ MESH) FOR COMPUTED
!     QUANTITIES:
!
!     BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF); BSUBS (FULL-CONVERTED IN JXBFORCE)
!
        IF (lnyquist) THEN
           mnmax_nyq0 = mnmax_nyq
           mnyq0 = mnyq
           nnyq0 = nnyq
           xm_nyq0 => xm_nyq; xn_nyq0 => xn_nyq
        ELSE
           mnmax_nyq0 = mnmax
           mnyq0 = mpol1
           nnyq0 = ntor
           xm_nyq0 => xm; xn_nyq0 => xn
        END IF

        ALLOCATE (gmn(mnmax_nyq0), bmn(mnmax_nyq0),
     &   bsubumn(mnmax_nyq0), bsubvmn(mnmax_nyq0), bsubsmn(mnmax_nyq0),
     &   bsupumn(mnmax_nyq0), bsupvmn(mnmax_nyq0), stat=istat) 
        ALLOCATE(rmnc(mnmax,ns),zmns(mnmax,ns),lmns(mnmax,ns),
     &   rmns(mnmax,ns),zmnc(mnmax,ns),lmnc(mnmax,ns))

        IF (lfreeb) THEN        !MRC 10-15-15
           ALLOCATE 
     & (bsubua_sur(nzeta*ntheta2), bsubva_sur(nzeta*ntheta2))
           ALLOCATE 
     & (bsupua_sur(nzeta*ntheta2), bsupva_sur(nzeta*ntheta2))

           ALLOCATE 
     & (bsubumnc_sur(mnmax_nyq0), bsubvmnc_sur(mnmax_nyq0))
           ALLOCATE 
     & (bsupumnc_sur(mnmax_nyq0), bsupvmnc_sur(mnmax_nyq0))
           IF (lasym) THEN
              ALLOCATE (bsubumns_sur(mnmax_nyq0),
     &                bsubvmns_sur(mnmax_nyq0))
              ALLOCATE (bsupumns_sur(mnmax_nyq0),
     &                bsupvmns_sur(mnmax_nyq0))
           END IF
        END IF

        ALLOCATE (gmnc(mnmax_nyq0,ns), bmnc(mnmax_nyq0,ns),
     &          bsubumnc(mnmax_nyq0,ns), bsubvmnc(mnmax_nyq0,ns),
     &          bsubsmns(mnmax_nyq0,ns), bsupumnc(mnmax_nyq0,ns),
     &          bsupvmnc(mnmax_nyq0,ns), 
     &          currumnc(mnmax_nyq0,ns), currvmnc(mnmax_nyq0,ns))

        IF (lasym) THEN
          ALLOCATE (gmns(mnmax_nyq0,ns), bmns(mnmax_nyq0,ns),
     &          bsubumns(mnmax_nyq0,ns), bsubvmns(mnmax_nyq0,ns),
     &          bsubsmnc(mnmax_nyq0,ns), bsupumns(mnmax_nyq0,ns),
     &          bsupvmns(mnmax_nyq0,ns),
     &          currumns(mnmax_nyq0,ns), currvmns(mnmax_nyq0,ns))
        END IF
        IF (istat .ne. 0) STOP 
     & 'Error allocating arrays in compute_output_representation'
        n1 = MAX(1,ntmax/2)
        rmnc = rzl_array(:,:,1)            !!store COS(mu-nv) components
        zmns = rzl_array(:,:,1+n1)         !!store SIN(mu-nv)
        lmns = rzl_array(:,:,1+2*n1)       !!store SIN(mu-nv)

        IF (lasym) THEN
           rmns = gc_array(:,:,1)            !!store SIN(mu-nv)
           zmnc = gc_array(:,:,1+n1)         !!store COS(mu-nv)
           lmnc = gc_array(:,:,1+2*n1)       !!store COS(mu-nv)
        END IF


        indx1=MAXLOC(ns_array)
        ftolx1=ftol_array(indx1(1))

!     NYQUIST FREQUENCY REQUIRES FACTOR OF 1/2
        IF (lnyquist) THEN
           IF (mnyq .ne. 0) cosmui(:,mnyq) = p5*cosmui(:,mnyq)
           IF (nnyq .ne. 0) cosnv (:,nnyq) = p5*cosnv (:,nnyq)
        END IF

        ALLOCATE (xfinal(neqs), stat=js)
        IF (js .NE. 0) STOP 
     &  'Allocation error for xfinal in compute_output_representation'
        xfinal = xc
!
!     MUST CONVERT m=1 MODES... FROM INTERNAL TO PHYSICAL FORM
!     Extrapolation of m=0 Lambda (cs) modes, which are not evolved at j=1, done in CONVERT
!
        lk = ns*ntor1
        IF (lthreed) CALL convert_sym  (xfinal(1+mns*(rss-1)+lk), 
     &                                xfinal(1+irzloff+mns*(zcs-1)+lk))
        IF (lasym)   CALL convert_asym (xfinal(1+mns*(rsc-1)+lk), 
     &                                xfinal(1+irzloff+mns*(zcc-1)+lk))

!
!     CONVERT TO rmnc, zmns, lmns, etc EXTERNAL representation (without internal mscale, nscale)
!     IF B^v ~ phip + lamu, MUST DIVIDE BY phipf(js) below to maintain old-style format
!     THIS COULD BE A PROBLEM FOR RFP WHERE PHIPF->0 INSIDE THE PLASMA!
!
        RADIUS1: DO js = 1, ns

           CALL convert (rmnc1, zmns1, lmns1, rmns1, zmnc1, lmnc1, 
     &                         xfinal, js)

           rmnc(:,js) = rmnc1(:)
           zmns(:,js) = zmns1(:)
           lmns(:,js) = (lmns1(:)/phipf(js)) * lamscale
           IF (lasym) THEN
              rmns(:,js) = rmns1(:)
              zmnc(:,js) = zmnc1(:)
              lmnc(:,js) = (lmnc1(:)/phipf(js)) * lamscale
           END IF

        END DO RADIUS1

        DEALLOCATE (xfinal)

!
!     INTERPOLATE LAMBDA ONTO HALF-MESH FOR BACKWARDS CONSISTENCY WITH EARLIER VERSIONS OF VMEC
!     AND SMOOTHS POSSIBLE UNPHYSICAL "WIGGLE" ON RADIAL MESH
!

        WHERE (NINT(xm) .le. 1) lmns(:,1) = lmns(:,2)
        DO js = ns,2,-1
           WHERE (MOD(NINT(xm),2) .eq. 0) 
              lmns(:,js) = p5*(lmns(:,js) + lmns(:,js-1))
           ELSEWHERE
              lmns(:,js) = 
     &  p5*(sm(js)*lmns(:,js) + sp(js-1)*lmns(:,js-1))
           END WHERE
        END DO

        lmns(:,1) = 0  
        raxis_cc(0:ntor) = rmnc(1:ntor+1,1)
        zaxis_cs(0:ntor) = zmns(1:ntor+1,1)
        
        IF (.NOT.lasym) GOTO 900

        WHERE (NINT(xm) .le. 1) lmnc(:,1) = lmnc(:,2)
        DO js = ns,2,-1
           WHERE (MOD(NINT(xm),2) .eq. 0) 
              lmnc(:,js) = p5*(lmnc(:,js) + lmnc(:,js-1))
           ELSEWHERE
              lmnc(:,js) = 
     &  p5*(sm(js)*lmnc(:,js) + sp(js-1)*lmnc(:,js-1))
           END WHERE
        END DO

        lmnc(:,1) = 0;   
        raxis_cs(0:ntor) = rmns(1:ntor+1,1)
        zaxis_cc(0:ntor) = zmnc(1:ntor+1,1)

 900  CONTINUE
        DO js = 2, ns
           bsq(js,:nznt) = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
        END DO

        tmult = p5/r0scale**2
!SPH: FIXED THIS 03-05-07 TO CALL symmetrization routine
        IF (lasym) THEN
!Changed integration norm in fixaray, SPH012314
           tmult = 2*tmult
           bsubs(1,:) = 0
           CALL symoutput (bsq,   gsqrt,  bsubu,  bsubv,  bsupu,
     &                   bsupv,  bsubs, 
     &                   bsqa,  gsqrta, bsubua, bsubva, bsupua,
     &                   bsupva, bsubsa)

           IF (lfreeb) THEN     !MRC  10-15-15
              CALL symoutput_sur(bsubu_sur, bsubv_sur,
     &                         bsupu_sur, bsupv_sur,
     &                         bsubua_sur, bsubva_sur,
     &                         bsupua_sur, bsupva_sur)
           END IF
        END IF

!         DO js = 2, ns
!            WRITE (200, *) 'JS: ', js, 'BSUBU, BSUBV'
!            WRITE (200, '(1p,6e12.4)') bsubu(js,:), bsubv(js,:)
!         END DO

        RADIUS2: DO js = 2, ns
           gmn = 0
           bmn = 0
           bsubumn = 0
           bsubvmn = 0
           bsubsmn = 0
           bsupumn = 0
           bsupvmn = 0

           MN2: DO mn = 1, mnmax_nyq0
              n = NINT(xn_nyq0(mn))/nfp
              m = NINT(xm_nyq0(mn))
              n1 = ABS(n)
              dmult = mscale(m)*nscale(n1)*tmult
              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
              sgn = SIGN(1, n)
              lk = 0
              DO j = 1, ntheta2
                 DO k = 1, nzeta
                    lk = lk + 1 
                    tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     &                       sgn*sinmui(j,m)*sinnv(k,n1))
                    tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     &                       sgn*cosmui(j,m)*sinnv(k,n1))
                    bmn(mn) = bmn(mn) + tcosi*bsq(js,lk)
                    gmn(mn) = gmn(mn) + tcosi*gsqrt(js,lk)
                    bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lk)
                    bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lk)
                    bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lk)
                    bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lk) 
                    bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lk) 
                 END DO
              END DO
           END DO MN2

           IF (js .eq. ns/2) bmodmn = bmn(1:mnmax)
           IF (js .eq. ns) bmodmn1 = bmn(1:mnmax)
           gmnc(:,js) = gmn(:)
           bmnc(:,js) = bmn(:)
           bsubumnc(:,js) = bsubumn(:)
           bsubvmnc(:,js) = bsubvmn(:)
           bsubsmns(:,js) = bsubsmn(:)
           bsupumnc(:,js) = bsupumn(:)
           bsupvmnc(:,js) = bsupvmn(:)
        END DO RADIUS2

        IF (lfreeb) THEN    !MRC    10-15-15
           bsubumnc_sur = 0
           bsubvmnc_sur = 0
           bsupumnc_sur = 0
           bsupvmnc_sur = 0
           DO mn = 1, mnmax_nyq0
              n = NINT(xn_nyq0(mn))/nfp
              m = NINT(xm_nyq0(mn))
              n1 = ABS(n)
              dmult = mscale(m)*nscale(n1)*tmult
              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
              sgn = SIGN(1, n)
              lk = 0
              DO j = 1, ntheta2
                 DO k = 1, nzeta
                    lk = lk + 1
                    tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     &                         sgn*sinmui(j,m)*sinnv(k,n1))
                    bsubumnc_sur(mn) = bsubumnc_sur(mn)
     &                               + tcosi*bsubu_sur(lk)
                    bsubvmnc_sur(mn) = bsubvmnc_sur(mn)
     &                               + tcosi*bsubv_sur(lk)
                    bsupumnc_sur(mn) = bsupumnc_sur(mn)
     &                               + tcosi*bsupu_sur(lk)
                    bsupvmnc_sur(mn) = bsupvmnc_sur(mn)
     &                               + tcosi*bsupv_sur(lk)
                 END DO
              END DO
           END DO
        END IF

        gmnc(:,1) = 0; bmnc(:,1) = 0;
        bsubumnc(:,1) = 0
        bsubvmnc(:,1) = 0
        bsubsmns(:,1) = 2*bsubsmns(:,2) - bsubsmns(:,3)
        bsupumnc(:,1) = 0;  bsupvmnc(:,1) = 0

        IF (.not.lasym) GO TO 200

        RADIUS3: DO js = 2, ns
           gmn = 0
           bmn = 0
           bsubumn = 0
           bsubvmn = 0
           bsubsmn = 0
           bsupumn = 0
           bsupvmn = 0

           MN3: DO mn = 1, mnmax_nyq0
              n = NINT(xn_nyq0(mn))/nfp
              m = NINT(xm_nyq0(mn))
              n1 = ABS(n)
              dmult = mscale(m)*nscale(n1)*tmult
              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
              sgn = SIGN(1, n)
              lk = 0
              jlk = js
              DO j = 1, ntheta2
                 DO k = 1, nzeta
                    lk = lk + 1
                    tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
     &                         sgn*sinmui(j,m)*sinnv(k,n1))
                    tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     &                         sgn*cosmui(j,m)*sinnv(k,n1))
                    bmn(mn) = bmn(mn) + tsini*bsqa(jlk)
                    gmn(mn) = gmn(mn) + tsini*gsqrta(jlk,0)
                    bsubumn(mn) = bsubumn(mn) + tsini*bsubua(jlk)
                    bsubvmn(mn) = bsubvmn(mn) + tsini*bsubva(jlk)
                    bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubsa(jlk)
                    bsupumn(mn) = bsupumn(mn) + tsini*bsupua(jlk)
                    bsupvmn(mn) = bsupvmn(mn) + tsini*bsupva(jlk)

                    jlk = jlk+ns
                 END DO
              END DO
           END DO MN3
     
           gmns(:,js) = gmn(:)
           bmns(:,js) = bmn(:)
           bsubumns(:,js) = bsubumn(:)
           bsubvmns(:,js) = bsubvmn(:)
           bsubsmnc(:,js) = bsubsmn(:)
           bsupumns(:,js) = bsupumn(:)
           bsupvmns(:,js) = bsupvmn(:)
        END DO RADIUS3

        gmns(:,1) = 0; bmns(:,1) = 0
        bsubumns(:,1) = 0
        bsubvmns(:,1) = 0
        bsubsmnc(:,1) = 2*bsubsmnc(:,2) - bsubsmnc(:,3)
        bsupumns(:,1) = 0;  bsupvmns(:,1) = 0

        IF (lfreeb) THEN        !MRC  10-15-15
           bsubumns_sur = 0
           bsubvmns_sur = 0
           bsupumns_sur = 0
           bsupvmns_sur = 0

           DO mn = 1, mnmax_nyq0
              n = NINT(xn_nyq0(mn))/nfp
              m = NINT(xm_nyq0(mn))
              n1 = ABS(n)
              dmult = mscale(m)*nscale(n1)*tmult
              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
              sgn = SIGN(1, n)
              lk = 0
              DO j = 1, ntheta2
                 DO k = 1, nzeta
                    lk = lk + 1
                    tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
     &                         sgn*cosmui(j,m)*sinnv(k,n1))
                    bsubumns_sur(mn) = bsubumns_sur(mn)
     &                               + tsini*bsubua_sur(lk)
                    bsubvmns_sur(mn) = bsubvmns_sur(mn)
     &                               + tsini*bsubva_sur(lk)
                    bsupumns_sur(mn) = bsupumns_sur(mn)
     &                               + tsini*bsupua_sur(lk)
                    bsupvmns_sur(mn) = bsupvmns_sur(mn)
     &                               + tsini*bsupva_sur(lk)
                 END DO
              END DO
           END DO
        END IF

 200  CONTINUE

        CALL Compute_Currents(bsubsmnc, bsubsmns, bsubumnc, bsubumns,
     &                        bsubvmnc, bsubvmns,
     &                        xm_nyq0, xn_nyq0, mnmax_nyq0, lasym, ns,
     &                        currumnc, currvmnc, currumns, currvmns)


!     RESTORE nyq ENDPOINT VALUES

        IF (lnyquist) THEN
           IF (mnyq .ne. 0) cosmui(:,mnyq) = 2*cosmui(:,mnyq)
           IF (nnyq .ne. 0) cosnv (:,nnyq) = 2*cosnv (:,nnyq)
        END IF

        IF (ALLOCATED(bsubumnc_sur)) THEN
           DEALLOCATE(bsubumnc_sur, bsubvmnc_sur)
           DEALLOCATE(bsupumnc_sur, bsupvmnc_sur)
        END IF
        IF (ALLOCATED(bsubumns_sur)) THEN
           DEALLOCATE(bsubumns_sur, bsubvmns_sur)
           DEALLOCATE(bsupumns_sur, bsupvmns_sur)
        END IF
        IF (ALLOCATED(bsubua_sur)) THEN
           DEALLOCATE(bsubua_sur, bsubva_sur)
           DEALLOCATE(bsupua_sur, bsupva_sur)
        END IF

        CALL freeb_data(rmnc1, zmns1, rmns1, zmnc1, bmodmn, bmodmn1)

        rzl_array = 0
      END SUBROUTINE vmec_output_data

      SUBROUTINE initialize(m_pol,n_tor)
!     &  BIND(C,name='initialize_vmec')
        USE vmec_input
        USE vparams
        INTEGER(C_INT), INTENT(IN) :: m_pol, n_tor

        mpol = m_pol
        ntor = n_tor
        omp_num_threads = 8
        gamma = 0
        spres_ped = 1
        ntheta = 0;  nzeta = 0
        ns_array = 0;  ns_array(1) = ns_default
        niter_array = -1;
        bloat = 1
        rbc = 0;  rbs = 0; zbs = 0; zbc = 0
        time_slice = 0
        nfp = 1
        ncurr = 0
        nsin = ns_default
        niter = 100
        nstep = 10
        nvacskip = 1
        delt = 1
        ftol = 1.E-10_dp
        ftol_array = 0;  ftol_array(1) = ftol
        am = 0; ai = 0; ac = 0; aphi = 0; aphi(1) = 1
        pres_scale = 1
        raxis_cc = 0; zaxis_cs = 0; raxis_cs = 0; zaxis_cc = 0;
        mfilter_fbdy = -1; nfilter_fbdy = -1
        tcon0 = 1
        precon_type = 'NONE'; prec2d_threshold = 1.E-30_dp
        curtor = 0; 
        extcur = 0;  phiedge = 1;
        mgrid_file = 'NONE'
        trip3d_file = 'NONE' ! SAL - TRIP3D
        lfreeb = .true.
        lmove_axis = .true.
        lmac = .false.
        lforbal = .false.
        lasym = .false.
        lrfp = .false.
        loldout = .false.        ! J Geiger 2010-05-04 start
        ldiagno = .false.
        lgiveup = .false.        ! inserted M.Drevlak
        fgiveup = 3.E+01_dp      ! inserted M.Drevlak
        lbsubs = .false.         ! J Hanson. See jxbforce coding
        lfull3d1out = .false.
        lmovie = .false.         ! S Lazerson for making movie files
        lmoreiter = .false.      ! default value if no max_main_iterations given.
        max_main_iterations = 1  ! to keep a presumably expected standard behavior.
!DEC$ IF DEFINED (NETCDF)
        lwouttxt = .false.       ! to keep functionality as expected with netcdf
!DEC$ ELSE
        lwouttxt = .true.        ! and without netcdf
!DEC$ ENDIF

        pcurr_type = 'power_series'
        piota_type = 'power_series'
        pmass_type = 'power_series'

!     ANISTROPY PARAMETERS
        bcrit = 1
        at(0) = 1;  at(1:) = 0
        ah = 0
        ph_type = 'power_series'
        pt_type = 'power_series'
        ah_aux_s(:) = -1
        at_aux_s(:) = -1
        am_aux_s(:) = -1
        ac_aux_s(:) = -1
        ai_aux_s(:) = -1
      

!
!     BACKWARDS COMPATIBILITY
!
        raxis = 0;  zaxis = 0

      END SUBROUTINE initialize

      SUBROUTINE finalize() !BIND(C,name='finalize_vmec')
        USE vmec_persistent, ONLY: xm, xn, xm_nyq, xn_nyq
        USE vmec_main
        IF (ALLOCATED(xm)) DEALLOCATE(xm,xn,xm_nyq,xn_nyq)
        IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc,zmns,lmns,rmns,zmnc,lmnc)
        IF (ALLOCATED(gmnc)) DEALLOCATE(gmnc, bmnc, bsubumnc, bsubvmnc,
     &                                bsubsmns, bsupumnc, bsupvmnc,
     &                                currumnc, currvmnc)
        IF (ALLOCATED(gmns)) DEALLOCATE(gmns, bmns, bsubumns, bsubvmns,
     &                                  bsubsmnc, bsupumns, bsupvmns,
     &                                  currumns, currvmns)
! J Geiger: check also for allocation.
        IF (ALLOCATED(gmn)) DEALLOCATE (gmn, bmn, bsubumn, bsubvmn,
     &               bsubsmn, bsupumn, bsupvmn)
        IF (ALLOCATED(iotaf)) THEN
         DEALLOCATE(iotaf, mass, phi, presf, jcuru, jcurv, jdotb, buco,
     &              bvco, bucof, bvcof, chi,
     &              bdotgradv, equif, specw, tcon, psi, yellip, yinden,
     &              ytrian, yshift, ygeo, overr, faclam, iotas, phips,
     &              chips, pres, vp, beta_vol, jperp2, jpar2, bdotb,
     &              clam, blam, dlam, phipf, chipf)
        END IF
!        IF (ALLOCATED(cosmu))
!     &   DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
!     &             sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn,
!     &             cosmui3, cosmumi3, cos01, sin01)

!        IF (ALLOCATED(tanu))
!     &   DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, sinu,
!     &              cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
!     &              cosu1, sinv1, cosv1, imirr, xmpot, xnpot)

      END SUBROUTINE finalize

      SUBROUTINE set_vmec_data_real(data_size,data_id,data_ptr)
!     &  BIND(C,name='set_vmec_data_real')
        USE vmec_input
        USE vparams, ONLY: ntord, ndatafmax
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        REAL(C_DOUBLE), DIMENSION(data_size), INTENT(IN) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_string
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: data_2d
        REAL(dp), DIMENSION(:), ALLOCATABLE :: data_1d
        INTEGER :: ix, ix1, ix2

        CALL c2f_string_array_1d(data_id,data_string)

        IF (TRIM(data_string) .eq. "rbc") THEN
          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
            DO ix = 1,data_size
              ix1 = MOD(ix-1,2*ntor+1) - ntor
              ix2 = (ix-1)/(2*ntor+1)
              data_2d(ix1, ix2) = data_ptr(ix)
            END DO
            rbc(-ntor:ntor,0:mpol-1) = data_2d
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "rbs") THEN
          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
            DO ix = 1,data_size
              ix1 = MOD(ix-1,2*ntor+1) - ntor
              ix2 = (ix-1)/(2*ntor+1)
              data_2d(ix1, ix2) = data_ptr(ix)
            END DO
            rbs(-ntor:ntor,0:mpol-1) = data_2d
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "zbc") THEN
          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
            DO ix = 1,data_size
              ix1 = MOD(ix-1,2*ntor+1) - ntor
              ix2 = (ix-1)/(2*ntor+1)
              data_2d(ix1, ix2) = data_ptr(ix)
            END DO
            zbc(-ntor:ntor,0:mpol-1) = data_2d
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "zbs") THEN
          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
            DO ix = 1,data_size
              ix1 = MOD(ix-1,2*ntor+1) - ntor
              ix2 = (ix-1)/(2*ntor+1)
              data_2d(ix1, ix2) = data_ptr(ix)
            END DO
            zbs(-ntor:ntor,0:mpol-1) = data_2d
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "am_aux_s") THEN
          IF (data_size .gt. ndatafmax) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          am_aux_s(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "am_aux_f") THEN
          IF (data_size .gt. ndatafmax) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          am_aux_f(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ai_aux_s") THEN
          IF (data_size .gt. ndatafmax) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ai_aux_s(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ai_aux_f") THEN
          IF (data_size .gt. ndatafmax) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ai_aux_f(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ac_aux_s") THEN
          IF (data_size .gt. ndatafmax) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ac_aux_s(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ac_aux_f") THEN
          IF (data_size .gt. ndatafmax) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ac_aux_f(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ai") THEN
          IF (data_size .gt. 21) THEN
            WRITE(6,"(A)") "Fatal VMEC input error!"
            WRITE(6,"(A)") "Must specify <= 21 expansion coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ai(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ac") THEN
          IF (data_size .gt. 21) THEN
            WRITE(6,"(A)") "Fatal VMEC input error!"
            WRITE(6,"(A)") "Must specify <= 21 expansion coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ac(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "am") THEN
          IF (data_size .gt. 21) THEN
            WRITE(6,"(A)") "Fatal VMEC input error!"
            WRITE(6,"(A)") "Must specify <= 21 expansion coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          am(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "aphi") THEN
          IF (data_size .gt. 20) THEN
            WRITE(6,"(A)") "Fatal VMEC input error!"
            WRITE(6,"(A)") "Must specify <= 20 expansion coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          aphi(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "raxis") THEN
          IF (data_size .gt. ntord+1) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          raxis(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "zaxis") THEN
          IF (data_size .gt. ntord+1) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          zaxis(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "raxis_cc") THEN
          IF (data_size .gt. ntord+1) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          raxis_cc(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "raxis_cs") THEN
          IF (data_size .gt. ntord+1) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          raxis_cs(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "zaxis_cc") THEN
          IF (data_size .gt. ntord+1) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          zaxis_cc(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "zaxis_cs") THEN
          IF (data_size .gt. ntord+1) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
     &        "coefficients!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          zaxis_cs(0:data_size-1) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "ftol_array") THEN
          IF (data_size .gt. 100) THEN
            WRITE(6,"(A)") "Fatal VMEC Input Error" 
            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= 100 coefficients"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ftol_array(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "delt") THEN
          delt = data_ptr(1) 
        ELSE IF (TRIM(data_string) .eq. "prec2d_threshold") THEN
          prec2d_threshold = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "phiedge") THEN
          phiedge = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "gamma") THEN
          gamma = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "bloat") THEN
          bloat = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "spres_ped") THEN
          spres_ped = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "pres_scale") THEN
          pres_scale = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "curtor") THEN
          curtor = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "tcon0") THEN
          tcon0 = data_ptr(1)
        ELSE
          WRITE(6,"(A)") "Fatal VMEC input error!"
          WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_string)//" found!"
          STOP 
        END IF
      END SUBROUTINE

      SUBROUTINE set_vmec_data_int(data_size,data_id,data_ptr)
!     &  BIND(C,name='set_vmec_data_int')
        USE vmec_input
        USE vparams
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        INTEGER(C_INT), DIMENSION(data_size), INTENT(IN) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_string
        INTEGER, DIMENSION(:,:), ALLOCATABLE :: data_2d
        INTEGER, DIMENSION(:), ALLOCATABLE :: data_1d
        INTEGER :: ix1, ix2

        CALL c2f_string_array_1d(data_id,data_string)

        IF (TRIM(data_string) .eq. "ns_array") THEN
          IF (data_size .gt. 100) THEN
            WRITE(6,"(A)") "Fatal VMEC input error!"
            WRITE(6,"(A)") "ns_array cannot be larger than 100!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          ns_array(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "niter_array") THEN
          IF (data_size .gt. 100) THEN
            WRITE(6,"(A)") "Fatal VMEC input error!"
            WRITE(6,"(A)") "niter_array cannot be larger than 100!"
            STOP
          END IF
          ALLOCATE(data_1d(data_size))
          data_1d = data_ptr
          niter_array(1:data_size) = data_1d
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "mpol") THEN
          mpol = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "ntor") THEN
          ntor = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "niter") THEN
          niter = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "nfp") THEN
          nfp = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "nsin") THEN
          nsin = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "nstep") THEN
          nstep = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "ncurr") THEN
          ncurr = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "nvacskip") THEN
          nvacskip = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "ntheta") THEN
          ntheta = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "nzeta") THEN
          nzeta = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "mfilter_fbdy") THEN
          mfilter_fbdy = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "nfilter_fbdy") THEN
          nfilter_fbdy = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "max_main_iterations") THEN
          max_main_iterations = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "omp_num_threads") THEN
          omp_num_threads = data_ptr(1)
        ELSE
          WRITE(6,"(A)") "Fatal VMEC input error!"
          WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_string)//" found!"
          STOP 
        END IF

      END SUBROUTINE set_vmec_data_int

      SUBROUTINE set_vmec_data_bool(data_size,data_id,data_ptr)
!     &  BIND(C,name='set_vmec_data_bool')
        USE vmec_input
        USE vparams
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        LOGICAL(C_BOOL), DIMENSION(data_size), INTENT(IN) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_string
        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: data_2d
        LOGICAL, DIMENSION(:), ALLOCATABLE :: data_1d
        INTEGER :: ix1, ix2

        CALL c2f_string_array_1d(data_id,data_string)

        IF (TRIM(data_string) .eq. "lasym") THEN
          lasym = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "lfreeb") THEN
          lfreeb = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "lwouttxt") THEN
          lwouttxt = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "lpofr") THEN
          lpofr = data_ptr(1)
        ELSE IF (TRIM(data_string) .eq. "full3d1out") THEN
          lfull3d1out = data_ptr(1)
        ELSE
          WRITE(6,"(A)") "Fatal VMEC input error!"
          WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_string)//" found!"
          STOP 
        END IF

      END SUBROUTINE set_vmec_data_bool

      SUBROUTINE set_vmec_data_char(data_size,data_id,data_ptr)
!     &  BIND(C,name='set_vmec_data_char')
        USE vmec_input
        USE vparams
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        TYPE(C_PTR), INTENT(IN) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_name, data_string

        CALL c2f_string_array_1d(data_id,data_name)
        CALL c2f_string_array_1d(data_ptr,data_string)

        IF (TRIM(data_name) .eq. "precon_type") THEN
          precon_type = data_string 
        ELSE IF (TRIM(data_name) .eq. "pmass_type") THEN
          pmass_type = data_string
        ELSE IF (TRIM(data_name) .eq. "pcurr_type") THEN
          pcurr_type = data_string
        ELSE IF (TRIM(data_name) .eq. "piota_type") THEN
          piota_type = data_string
        ELSE IF (TRIM(data_name) .eq. "mgrid_file") THEN
          mgrid_file = data_string
        ELSE
          WRITE(6,"(A)") "Fatal VMEC input error!"
          WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_name)//" found!"
          STOP 
        END IF

      END SUBROUTINE set_vmec_data_char

      SUBROUTINE retrieve_vmec_data_real(data_size, data_id, data_ptr)
!     & BIND(C,name="retrieve_vmec_data_real")
        USE vmec_main
        USE vmec_persistent, ONLY: xm, xn, xm_nyq, xn_nyq
        USE vmec_io
        USE vmercier
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        REAL(C_DOUBLE), DIMENSION(data_size), INTENT(OUT) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_string
        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: data_2d
        REAL(dp), DIMENSION(:), ALLOCATABLE :: data_1d
        INTEGER :: ix1, ix2, cix

        CALL c2f_string_array_1d(data_id,data_string)

        IF (TRIM(data_string) .eq. "rmnc") THEN 
          ALLOCATE(data_2d(mnmax,ns))
          data_2d = rmnc
          DO ix1 = 1,mnmax
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "rmns") THEN
          ALLOCATE(data_2d(mnmax,ns))
          data_2d = rmns
          DO ix1 = 1,mnmax
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "zmns") THEN
          ALLOCATE(data_2d(mnmax,ns))
          data_2d = zmns
          DO ix1 = 1,mnmax
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (trim(data_string) .eq. "zmnc") THEN
          ALLOCATE(data_2d(mnmax,ns))
          data_2d = zmnc
          DO ix1 = 1,mnmax
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "lmns") THEN
          ALLOCATE(data_2d(mnmax,ns))
          data_2d = lmns
          DO ix1 = 1,mnmax
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (trim(data_string) .eq. "lmnc") THEN
          ALLOCATE(data_2d(mnmax,0:ns))
          data_2d = lmnc
          DO ix1 = 1,mnmax
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "gmnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = gmnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "gmns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = gmns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = ix2*mnmax_nyq + ix1-1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bmnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bmnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bmns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bmns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsubsmnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsubsmnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsubsmns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsubsmns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsubumnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsubumnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsubumns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsubumns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsubvmnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsubvmnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsubvmns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsubvmns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsupumnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsupumnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsupumns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsupumns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsupvmnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsupvmnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "bsupvmns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = bsupvmns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "currumnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = currumnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "currumns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = currumns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "currvmnc") THEN 
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = currvmnc
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "currvmns") THEN
          ALLOCATE(data_2d(mnmax_nyq,ns))
          data_2d = currvmns
          DO ix1 = 1,mnmax_nyq
            DO ix2 = 1,ns 
              cix = (ix2-1)*mnmax_nyq + ix1;
              data_ptr(cix) = data_2d(ix1,ix2); 
            END DO 
          END DO
          DEALLOCATE(data_2d)
        ELSE IF (TRIM(data_string) .eq. "xm") THEN 
          ALLOCATE(data_1d(mnmax))
          data_1d = xm
          DO ix1 = 1,mnmax
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "xn") THEN
          ALLOCATE(data_1d(mnmax))
          data_1d = xn
          DO ix1 = 1,mnmax
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "xm_nyq") THEN
          ALLOCATE(data_1d(mnmax_nyq))
          data_1d = xm_nyq
          DO ix1 = 1,mnmax_nyq
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "xn_nyq") THEN
          ALLOCATE(data_1d(mnmax_nyq))
          data_1d = xn_nyq
          DO ix1 = 1,mnmax_nyq
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "phi") THEN
          ALLOCATE(data_1d(ns))
          data_1d = phi
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "phipf") THEN
          ALLOCATE(data_1d(ns))
          data_1d = phipf
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "phips") THEN
          ALLOCATE(data_1d(ns+1))
          data_1d = phips
          DO ix1 = 1,ns+1
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "iotaf") THEN
          ALLOCATE(data_1d(ns+1))
          data_1d = iotaf
          DO ix1 = 1,ns+1
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "iotas") THEN
          ALLOCATE(data_1d(ns+1))
          data_1d = iotas
          DO ix1 = 1,ns+1
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "pres") THEN
          ALLOCATE(data_1d(ns+1))
          data_1d = pres
          DO ix1 = 1,ns+1
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "presf") THEN
          ALLOCATE(data_1d(ns))
          data_1d = presf
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "DCurr") THEN
          ALLOCATE(data_1d(ns))
          data_1d = DCurr
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "DGeod") THEN
          ALLOCATE(data_1d(ns))
          data_1d = DGeod
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "DMerc") THEN
          ALLOCATE(data_1d(ns))
          data_1d = DMerc
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "DShear") THEN
          ALLOCATE(data_1d(ns))
          data_1d = DShear
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "bdotb") THEN
          ALLOCATE(data_1d(ns))
          data_1d = bdotb
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "beta_vol") THEN
          ALLOCATE(data_1d(ns))
          data_1d = beta_vol
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "buco") THEN
          ALLOCATE(data_1d(ns))
          data_1d = buco
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "bvco") THEN
          ALLOCATE(data_1d(ns))
          data_1d = bvco
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "chi") THEN
          ALLOCATE(data_1d(ns))
          data_1d = chi
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "chipf") THEN
          ALLOCATE(data_1d(ns))
          data_1d = chipf
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "equif") THEN
          ALLOCATE(data_1d(ns))
          data_1d = equif
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "jcuru") THEN
          ALLOCATE(data_1d(ns))
          data_1d = jcuru
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "jcurv") THEN
          ALLOCATE(data_1d(ns))
          data_1d = jcurv
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "jdotb") THEN
          ALLOCATE(data_1d(ns))
          data_1d = jdotb
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "mass") THEN
          ALLOCATE(data_1d(ns))
          data_1d = mass
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "specw") THEN
          ALLOCATE(data_1d(ns))
          data_1d = specw
          DO ix1 = 1,ns
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "vp") THEN
          ALLOCATE(data_1d(ns+1))
          data_1d = vp
          DO ix1 = 1,ns+1
            data_ptr(ix1) = data_1d(ix1)
          END DO
          DEALLOCATE(data_1d)
        ELSE IF (TRIM(data_string) .eq. "Aminor_p") THEN
          data_ptr(1) = Aminor_p
        ELSE IF (TRIM(data_string) .eq. "Rmajor_p") THEN
          data_ptr(1) = Rmajor_p
        ELSE IF (TRIM(data_string) .eq. "IonLarmor") THEN
          data_ptr(1) = IonLarmor
        ELSE IF (TRIM(data_string) .eq. "aspect") THEN
          data_ptr(1) = aspect
        ELSE IF (TRIM(data_string) .eq. "b0") THEN
          data_ptr(1) = b0
        ELSE IF (TRIM(data_string) .eq. "betapol") THEN
          data_ptr(1) = betapol
        ELSE IF (TRIM(data_string) .eq. "betator") THEN
          data_ptr(1) = betator
        ELSE IF (TRIM(data_string) .eq. "betatotal") THEN
          data_ptr(1) = betatot
        ELSE IF (TRIM(data_string) .eq. "betaxis") THEN
          data_ptr(1) = betaxis
        ELSE IF (TRIM(data_string) .eq. "gamma") THEN
          data_ptr(1) = gamma
        ELSE IF (TRIM(data_string) .eq. "volavgB") THEN
          data_ptr(1) = volavgB
        ELSE IF (TRIM(data_string) .eq. "volume_p") THEN
          data_ptr(1) = volume_p
        ELSE IF (TRIM(data_string) .eq. "wp") THEN
          data_ptr(1) = wp
        ELSE
         WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_string)//" found!"
        END IF

      END SUBROUTINE retrieve_vmec_data_real

      SUBROUTINE retrieve_vmec_data_int(data_size, data_id, data_ptr)
!     &  BIND(C,name='retrieve_vmec_data_int')
        USE vmec_dim, ONLY: ns, mnmax
        USE vmec_main
        USE vmec_input
        USE vmec_params, ONLY: signgs
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        INTEGER(C_INT), DIMENSION(data_size), INTENT(OUT) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_string
        INTEGER :: ix1, ix2, cix

        CALL c2f_string_array_1d(data_id,data_string)

        IF (TRIM(data_string) .eq. "nfp") THEN
          data_ptr(1) = nfp
        ELSE IF (TRIM(data_string) .eq. "mpol") THEN
          data_ptr(1) = mpol
        ELSE IF (TRIM(data_string) .eq. "ntor") THEN
          data_ptr(1) = ntor
        ELSE IF (TRIM(data_string) .eq. "ns") THEN
          data_ptr(1) = ns
        ELSE IF (TRIM(data_string) .eq. "mnmax") THEN
          data_ptr(1) = mnmax
        ELSE IF (TRIM(data_string) .eq. "mnmax_nyq") THEN
          data_ptr(1) = mnmax_nyq
        ELSE IF (TRIM(data_string) .eq. "signgs") THEN
          data_ptr(1) = NINT(signgs)
          print *, data_ptr(1)
        ELSE
          WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_string)//" found!"
        END IF

      END SUBROUTINE retrieve_vmec_data_int

      SUBROUTINE retrieve_vmec_data_bool(data_size, data_id, data_ptr)
!     &  BIND(C,name='retrieve_vmec_data_bool')
        USE vmec_main
        INTEGER(C_INT), INTENT(IN) :: data_size
        TYPE(C_PTR), INTENT(IN) :: data_id
        LOGICAL(C_BOOL), DIMENSION(data_size), INTENT(OUT) :: data_ptr
        CHARACTER(LEN=:), ALLOCATABLE :: data_string
        INTEGER :: ix1, ix2, cix

        CALL c2f_string_array_1d(data_id,data_string)

        IF (TRIM(data_string) .eq. "lasym") THEN 
          data_ptr(1) = lasym
        ELSE IF (TRIM(data_string) .eq. "lfreeb") THEN
          data_ptr(1) = lfreeb
        ELSE IF (TRIM(data_string) .eq. "lrecon") THEN
          data_ptr(1) = lrecon
        ELSE IF (TRIM(data_string) .eq. "lrfp") THEN
          data_ptr(1) = lrfp
        ELSE 
          WRITE(6,"(A)") "No data field with label "
     &      //TRIM(data_string)//" found!"
        END IF

      END SUBROUTINE retrieve_vmec_data_bool


      END MODULE vmec_ext_interface
