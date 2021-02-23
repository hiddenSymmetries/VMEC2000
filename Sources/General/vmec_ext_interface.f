c$$$      MODULE vmec_ext_interface
c$$$        USE, INTRINSIC :: ISO_C_BINDING
c$$$        USE stel_kinds
c$$$        USE stel_constants
c$$$        IMPLICIT NONE
c$$$        
c$$$        PUBLIC :: vmec_output_data
c$$$
c$$$        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: rmnc, rmns, zmns, 
c$$$     &   zmnc, lmns, lmnc
c$$$        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: gmnc, bmnc,
c$$$     &   gmns, bmns,
c$$$     &   bsubumnc, bsubvmnc, bsubsmns, bsubumns, bsubvmns, bsubsmnc,
c$$$     &   currumnc, currvmnc, currumns, currvmns
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: gmn, bmn,
c$$$     &   bsubumn, bsubvmn, bsubsmn, bsupumn, bsupvmn
c$$$
c$$$        REAL(dp), DIMENSION(:), POINTER ::   xm_nyq0, xn_nyq0
c$$$        REAL(dp), ALLOCATABLE, DIMENSION(:,:) :: bsupumnc, bsupumns,
c$$$     &   bsupvmnc, bsupvmns
c$$$
c$$$        LOGICAL :: lcurr
c$$$
c$$$        INTERFACE set_vmec_data
c$$$          MODULE PROCEDURE set_vmec_data_real
c$$$          MODULE PROCEDURE set_vmec_data_int
c$$$          MODULE PROCEDURE set_vmec_data_bool
c$$$          MODULE PROCEDURE set_vmec_data_char
c$$$        END INTERFACE
c$$$
c$$$        INTERFACE retrieve_vmec_data
c$$$          MODULE PROCEDURE set_vmec_data_real
c$$$          MODULE PROCEDURE set_vmec_data_int
c$$$          MODULE PROCEDURE set_vmec_data_bool
c$$$        END INTERFACE
c$$$
c$$$      CONTAINS
c$$$
c$$$      SUBROUTINE runvmec_ext(input_file0,EXT_MPI_COMM) 
c$$$!     & BIND(C,name='runvmec_c')
c$$$        USE vmec_input
c$$$        USE vmec_seq
c$$$        USE safe_open_mod
c$$$        USE vparams, ONLY: nlog, nlog0, nthreed
c$$$        USE vmec_params, ONLY: more_iter_flag,
c$$$     &                       bad_jacobian_flag,
c$$$     &    restart_flag, readin_flag, timestep_flag,
c$$$     &    output_flag, cleanup_flag,
c$$$     &    norm_term_flag, successful_term_flag ! J Geiger: for more iterations and full 3D1-output
c$$$        USE parallel_include_module, ONLY: grank, 
c$$$     &                                   MPI_ERR
c$$$        USE parallel_vmec_module, ONLY: MyEnvVariables, 
c$$$     &                                InitializeParallel,
c$$$     &                                FinalizeParallel
c$$$        IMPLICIT NONE
c$$$        TYPE(C_PTR), INTENT(IN) :: input_file0
c$$$        INTEGER, INTENT(IN) :: EXT_MPI_COMM
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: file_string
c$$$
c$$$C-----------------------------------------------
c$$$C   L o c a l   P a r a m e t e r s
c$$$C-----------------------------------------------
c$$$        INTEGER, PARAMETER :: nseq0 = 12
c$$$        CHARACTER(LEN=*), PARAMETER ::
c$$$     &    increase_niter = "Try increasing NITER",
c$$$     &    bad_jacobian = "The jacobian was non-definite!",
c$$$     &    full_3d1output_request = "Full threed1-output request!"
c$$$C-----------------------------------------------
c$$$C   L o c a l   V a r i a b l e s
c$$$C-----------------------------------------------
c$$$        INTEGER :: numargs, ierr_vmec, index_end,
c$$$     &   iopen, isnml, iread, iseq, index_seq,
c$$$     &   index_dat, iunit, ncount, nsteps, i
c$$$        INTEGER :: ictrl(5)
c$$$        CHARACTER(LEN=120) :: input_file, seq_ext, reset_file_name, arg
c$$$        CHARACTER(LEN=120) :: log_file
c$$$        CHARACTER(LEN=120), DIMENSION(10) :: command_arg
c$$$        LOGICAL :: lscreen
c$$$
c$$$C-----------------------------------------------
c$$$!***
c$$$!                              D   I   S   C   L   A   I   M   E   R
c$$$!
c$$$!       You are using a beta version of the PROGRAM VMEC, which is currently
c$$$!       under development by S. P. Hirshman at the Fusion Energy Division,
c$$$!       Oak Ridge National Laboratory.  Please report any problems or comments
c$$$!       to him.  As a beta version, this program is subject to change
c$$$!       and improvement without notice.
c$$$!
c$$$!       1. CODE SYNOPSIS
c$$$!
c$$$!       THIS PROGRAM - VMEC (Variational Moments Equilibrium Code)  -
c$$$!       SOLVES THREE-DIMENSIONAL MHD EQUILIBRIUM EQUATIONS USING
c$$$!       FOURIER SPECTRAL (MOMENTS) METHODS. A CYLINDRICAL COORDINATE
c$$$!       REPRESENTATION IS USED (R-Z COORDINATES). THE POLOIDAL
c$$$!       ANGLE VARIABLE IS RENORMALIZED THROUGH THE STREAM FUNCTION
c$$$!       LAMBDA, WHICH IS SELF-CONSISTENTLY DETERMINED AND DIFFERENCED
c$$$!       VARIATIONALLY ON THE HALF-RADIAL MESH. THE POLOIDAL ANGLE IS
c$$$!       DETERMINED BY MINIMIZING <M> = m**2 S(m) , WHERE S(m) =
c$$$!       Rm**2 + Zm**2 . AN EVEN-ODD DECOMPOSITION IN THE POLOIDAL MODE
c$$$!       NO. OF R,Z, AND LAMDA IS USED TO IMPROVE RADIAL RESOLUTION.
c$$$!       A FREE-BOUNDARY OPTION IS AVAILABLE (FOR lfreeb=T), WITH A
c$$$!       USER-SUPPLIED DATA-FILE "MGRID" NEEDED TO COMPUTE THE PLASMA
c$$$!       VACUUM FIELD COMPONENTS BR, BPHI, BZ (see SUBROUTINE BECOIL)
c$$$!
c$$$!       THE MAGNETIC FIELD IS REPRESENTED INTERNALLY AS FOLLOWS:
c$$$!
c$$$!       B(s,u,v) = grad(phiT) X ( grad(u) + grad(lambda) ) +
c$$$!
c$$$!                  iota(s) * grad(v) X grad(phiT)
c$$$!
c$$$!       WHERE phiT is the toroidal flux (called phi in code) and
c$$$!       u,v are the poloidal, toroidal angles, respectively.
c$$$!
c$$$!       2. ADDITIONAL CODES REQUIRED
c$$$!       For the fixed boundary calculation, the user must provide the Fourier
c$$$!       coefficients for the plasma boundary (the last surface outside of which
c$$$!       the pressure gradient vanishes). For ALL but the simplest geometry, the
c$$$!       SCRUNCH code (available from R. Wieland), based on the DESCUR curve-fitting
c$$$!       code, can be used to produce the optimized VMEC Fourier representation for
c$$$!       an arbritrary closed boundary (it need not be a 'star-like' DOmain, nor
c$$$!       need it possess vertical, or 'stellarator', symmetry).
c$$$!
c$$$!       For the free boundary calculation, the MAKEGRID code (available upon
c$$$!       request) is needed to create a binary Green''s FUNCTION table for the
c$$$!       vacuum magnetic field(s) and, IF data analysis is to be done, flux and
c$$$!       field loops as well. The user provides a SUBROUTINE (BFIELD) which can be
c$$$!       called at an arbitrary spatial location and which should RETURN the three
c$$$!       cylindrical components of the vacuum field at that point. (Similary,
c$$$!       locations of diagnostic flux loops, Rogowski coils, etc. are required IF
c$$$!       equilibrium reconstruction is to be done.)
c$$$!
c$$$!       Plotting is handled by a stand-alone package, PROUT.NCARG (written by
c$$$!       R. M. Wieland). It uses NCAR-graphics calls and reads the primary VMEC output
c$$$!       file, WOUT.EXT, WHERE 'EXT' is the command-line extension of the INPUT file.
c$$$!
c$$$!
c$$$!       3. UNIX SCRIPT SETUP PARAMETERS
c$$$!       The VMEC source code (vmec.lsqh) is actually a UNIX script file which uses
c$$$!       the C-precompiler to produce both the machine-specific Fortran source and a
c$$$!       make-file specific to ANY one of the following platforms:
c$$$!
c$$$!       IBM-RISC6000, CRAY, ALPHA (DEC-STATION), HP-UX WORKSTATION,
c$$$!       WINDOWS-NT, DEC-VMS
c$$$!
c$$$!       Additional platforms are easy to add to the existing script as required.
c$$$!
c$$$!
c$$$!       4. FORTRAN PARAMETER STATEMENTS set by user
c$$$!       In the Fortran-90 version of VMEC these PARAMETER statements have
c$$$!       been replaced by dynamic memory allocation. So the user should set the
c$$$!       run-time parameters ns (through ns_array), mpol, ntor in the NAMELIST INDATA.
c$$$!
c$$$!
c$$$!       Added features since last edition (see vmec_params for revision history list)
c$$$!       1. Implemented preconditioning algorithm for R,Z
c$$$!       2. The physical (unpreconditioned) residuals are used
c$$$!          to determine the level of convergence
c$$$!       3. The original (MOMCON) scaling of lambda is used, i.e.,
c$$$!          Bsupu = phip*(iota - lamda[sub]v)/SQRT(g). This is needed to
c$$$!          maintain consistency with the time-stepper for arbitrary PHIP.
c$$$!
c$$$!       WRITTEN BY S. P. HIRSHMAN (8/28/85 - REVISED 3/1/86) BASED ON
c$$$!       1. S. P. Hirshman and J. C. Whitson, Phys. Fluids 26, 3553 (1983).
c$$$!       2. S. P. Hirshman and H. K. Meier, Phys. Fluids 28, 1387 (1985).
c$$$!       3. S. P. Hirshman and D. K. Lee, Comp. Phys. Comm. 39, 161 (1986).
c$$$!
c$$$
c$$$!     Local variables
c$$$!
c$$$!     ictrl:   array(5) of control variables for running "runvmec" routine
c$$$!              see "runvmec" for a description
c$$$!
c$$$
c$$$!
c$$$!     Read in command-line arguments to get input file or sequence file,
c$$$!     screen display information, and restart information
c$$$        INTERFACE
c$$$           SUBROUTINE runvmec(ictrl_array, input_file0, 
c$$$     &                      lscreen, EXT_MPI_COMM, reset_file_name)
c$$$             IMPLICIT NONE
c$$$             INTEGER, INTENT(inout), TARGET :: ictrl_array(5)
c$$$             LOGICAL, INTENT(in) :: lscreen
c$$$             CHARACTER(LEN=*), INTENT(in) :: input_file0
c$$$             INTEGER, INTENT(in), OPTIONAL :: EXT_MPI_COMM
c$$$             CHARACTER(LEN=*), OPTIONAL :: reset_file_name
c$$$           END SUBROUTINE runvmec
c$$$        END INTERFACE
c$$$
c$$$        CALL MyEnvVariables
c$$$        CALL InitializeParallel(EXT_MPI_COMM)
c$$$
c$$$        lscreen = .false.
c$$$        IF(grank.EQ.0) lscreen = .true.
c$$$        reset_file_name = " "
c$$$
c$$$!
c$$$!     Determine type of file opened (sequential or input-data)
c$$$!     ARG1 (char var)
c$$$!          By DEFAULT, ARG1 obtained from the command
c$$$!          line is parsed as follows to determine the input data file(s):
c$$$!               a. Attempt to OPEN file ARG1 (full path + file name).
c$$$!                  Look for the VSEQ NAMELIST to obtain nseq, nseq_select, and
c$$$!                  extension array. If they exist and nseq>0, VMEC will run
c$$$!                  sequentially using input determined from the array EXTENSION[i]
c$$$!                  or input.EXTENSION[i]
c$$$!               b. If the command argument is not a sequence NAMELIST, THEN the data file
c$$$!                  ARG1 or input.ARG1 is READ directly, with NSEQ=1.
c$$$!
c$$$        CALL c2f_string_array_1d(input_file0,file_string)
c$$$        arg = TRIM(file_string)
c$$$        index_dat = INDEX(arg,'.')
c$$$        index_end = LEN_TRIM(arg)
c$$$        IF (index_dat .gt. 0) THEN
c$$$           seq_ext  = arg(index_dat + 1:index_end)
c$$$           input_file = TRIM(arg)
c$$$        END IF
c$$$!
c$$$        nseq = 1
c$$$        nseq_select(1) = 1
c$$$        extension(1) = input_file
c$$$!
c$$$!        iunit = nseq0
c$$$!        IF (iopen .eq. 0) THEN
c$$$!          CALL read_namelist (iunit, isnml, 'indata')
c$$$!
c$$$!        END IF
c$$$!        CLOSE (iunit)
c$$$
c$$$!
c$$$!     CALL EQUILIBRIUM SOLVER
c$$$!
c$$$!     nseq_select:      If sequence file (VSEQ NAMELIST given with nseq >0)
c$$$!                       array giving indices into EXTENSION array prescribing
c$$$!                       the order in which the input files are run by VMEC
c$$$!     nseq:             number of sequential VMEC runs to make
c$$$!
c$$$!
c$$$!     CALL VMEC WITH POSSIBLE SEQUENCE EXTENSION (SEQ_EXT)
c$$$!     AND ARRAY OF INPUT FILE EXTENSIONS (EXTENSION)
c$$$!
c$$$        ictrl = 0
c$$$
c$$$
c$$$        iseq = 1
c$$$        index_seq = nseq_select(iseq)
c$$$        ictrl(1) = restart_flag + readin_flag + timestep_flag
c$$$     &            + output_flag + cleanup_flag                !Sets all flags
c$$$        ictrl(2) = 0
c$$$!         ictrl(3) = 100
c$$$!         ictrl(4) = 2
c$$$        ictrl(5) = iseq - 1
c$$$        ncount = 0
c$$$
c$$$        print *, ictrl
c$$$        CALL runvmec(ictrl, extension(index_seq), lscreen, EXT_MPI_COMM,
c$$$     &                reset_file_name)
c$$$
c$$$                print *, ictrl
c$$$        ierr_vmec = ictrl(2)
c$$$        print *, ierr_vmec
c$$$
c$$$        SELECT CASE (ierr_vmec)
c$$$          CASE (more_iter_flag)                                !Need a few more iterations to converge
c$$$            IF (grank .EQ. 0) THEN
c$$$              IF(lscreen) WRITE (6, '(1x,a)') increase_niter
c$$$              WRITE (nthreed, '(1x,a)') increase_niter
c$$$              WRITE (nthreed, '(1x,a)') "PARVMEC aborting..."
c$$$              CALL FLUSH(nthreed)
c$$$            END IF
c$$$! J Geiger: if lmoreiter and lfull3d1out are false
c$$$!           the o-lines (original) are the only
c$$$!           ones to be executed.
c$$$            IF (lmoreiter) THEN                                 ! J Geiger: --start--
c$$$              DO i = 2, max_main_iterations                    ! Changes to run
c$$$                ictrl(1) = timestep_flag                      ! some more iterations if requested
c$$$                ictrl(3) = niter                              ! - this is the number of iterations
c$$$                CALL runvmec(ictrl, extension(1), lscreen,
c$$$     &                            EXT_MPI_COMM, reset_file_name)       ! - the second iteration run with ictrl(3) iterations
c$$$                IF (ictrl(2) .EQ. more_iter_flag .and.
c$$$     &                   grank    .EQ. 0) THEN
c$$$                  WRITE (nthreed, '(1x,a)') increase_niter
c$$$                  IF(lscreen) WRITE (6, '(1x,a)') increase_niter
c$$$                END IF
c$$$              END DO
c$$$              ictrl(1) = output_flag + cleanup_flag            ! - Output, cleanup
c$$$              IF (ictrl(2) .ne. successful_term_flag) THEN
c$$$                ictrl(2)=successful_term_flag                 ! - force success flag to get full threed1-output!
c$$$              END IF
c$$$              ictrl(3) = 0                                     ! - this is the number of iterations
c$$$              CALL runvmec(ictrl, extension(1), lscreen, EXT_MPI_COMM,
c$$$     &                         reset_file_name)
c$$$            ELSE                                                ! else-branch contains original code.
c$$$#if defined(MPI_OPT)
c$$$              CALL MPI_Barrier(EXT_MPI_COMM, MPI_ERR)
c$$$#endif
c$$$
c$$$              ictrl(1) = output_flag + cleanup_flag            !Output, cleanup  ! o-lines
c$$$              ictrl(2) = 0                                     ! o-lines
c$$$!              IF (lfull3d1out) THEN
c$$$!                ictrl(2) = successful_term_flag
c$$$!                IF (grank .EQ. 0) THEN
c$$$!                  WRITE(6,'(1x,a)') full_3d1output_request
c$$$!                  WRITE(nthreed,'(1x,a)') full_3d1output_request
c$$$!                END IF
c$$$!              END IF
c$$$
c$$$              CALL runvmec(ictrl, extension(1), lscreen, EXT_MPI_COMM, ! o-lines
c$$$     &                         reset_file_name)
c$$$            END IF                                              ! J Geiger: -- end --
c$$$
c$$$          CASE (bad_jacobian_flag)                               !Bad jacobian even after axis reset and ns->3
c$$$            IF (grank .EQ. 0) THEN
c$$$              IF (lscreen) WRITE (6, '(/,1x,a)') bad_jacobian
c$$$                WRITE (nthreed, '(/,1x,a)') bad_jacobian
c$$$              END IF
c$$$          CASE DEFAULT
c$$$        END SELECT
c$$$
c$$$
c$$$        CALL FinalizeParallel
c$$$
c$$$      END SUBROUTINE runvmec_ext
c$$$
c$$$      SUBROUTINE c2f_string_array_1d(c_pointer,f_string)
c$$$        USE, INTRINSIC :: ISO_C_BINDING
c$$$        TYPE(C_PTR), INTENT(IN) :: c_pointer
c$$$        CHARACTER(LEN=:), ALLOCATABLE, INTENT(OUT) :: f_string
c$$$        CHARACTER, DIMENSION(:), POINTER :: f_pointer
c$$$        CHARACTER(LEN=:), POINTER :: f_char_pointer
c$$$
c$$$        INTEGER :: i, j
c$$$        LOGICAL :: null_char_found
c$$$        i = 0
c$$$        null_char_found = .FALSE.
c$$$        DO WHILE (null_char_found .EQV. .FALSE.)
c$$$          i = i + 1
c$$$          CALL C_F_POINTER(c_pointer,f_pointer,[i])
c$$$          IF (f_pointer(i) == C_NULL_CHAR) THEN
c$$$            null_char_found = .TRUE.
c$$$          END IF
c$$$        END DO
c$$$        i = i - 1
c$$$        
c$$$        ALLOCATE(CHARACTER(LEN=i)::f_string)
c$$$        CALL C_F_POINTER(C_LOC(f_pointer),f_char_pointer)
c$$$
c$$$        f_string = f_char_pointer(1:i)
c$$$
c$$$      END SUBROUTINE c2f_string_array_1d
c$$$    
c$$$
c$$$      SUBROUTINE vmec_output_data(bsq, gsqrt,
c$$$     &  bsubu, bsubv, bsubs, bsupv, bsupu,rzl_array, gc_array)
c$$$        USE vmec_input, ONLY: ns_array, ftol_array, lwouttxt
c$$$        USE vmec_params
c$$$        USE vmec_main
c$$$        USE vmercier
c$$$        USE vmec_persistent
c$$$        USE vparams, p5 => cp5, two => c2p0
c$$$        USE vac_persistent
c$$$        USE vspline
c$$$        USE xstuff
c$$$        USE vmec_io
c$$$        USE realspace, ONLY: phip, chip, gsqrta=>z1, z1=>z1
c$$$        USE totzsp_mod
c$$$        USE vforces, ONLY: bsupua=>brmn_e, bsupva=>czmn_o, 
c$$$     &                   bsqa=>bzmn_e, bsubsa=>armn_e,
c$$$     &                   bsubua=>azmn_e, bsubva=>armn_o
c$$$        USE vacmod, ONLY: potvac, mnpd, xmpot, xnpot       !added for diagno, J.Geiger
c$$$        USE read_wout_mod, ONLY: Compute_Currents
c$$$        USE mgrid_mod
c$$$
c$$$        IMPLICIT NONE
c$$$!-----------------------------------------------
c$$$!   D u m m y   A r g u m e n t s
c$$$!-----------------------------------------------
c$$$        REAL(dp), DIMENSION(mnmax,ns,3*MAX(ntmax/2,1)),           !reverse ns, mnmax for backwards compatibility
c$$$     &   INTENT(inout), TARGET :: rzl_array, gc_array
c$$$        REAL(dp), DIMENSION(ns,nznt), INTENT(inout) ::
c$$$     &   bsq, gsqrt, bsubu, bsubv, bsubs, bsupv, bsupu
c$$$        REAL(dp) :: qfact(ns)
c$$$        REAL(dp), PARAMETER :: c1p5 = 1.5_dp
c$$$        LOGICAL :: lnyquist = .TRUE.                               !=false, suppress nyquist stuff
c$$$!-----------------------------------------------
c$$$!   L o c a l   V a r i a b l e s
c$$$!-----------------------------------------------
c$$$        INTEGER :: j, js, jlk, mn, lk, iasym,
c$$$     &           m, n, k, n1, istat, i, indx1(1),
c$$$     &           mnmax_nyq0, mnyq0, nnyq0
c$$$     &          ,isgn, js2, nfort      !for diagno 1.5
c$$$        REAL(dp) :: dmult, tcosi, tsini, vversion, sgn, tmult, 
c$$$     &            presfactor, ftolx1, d_bsupumn, d_bsupvmn   ! diagno 1.5
c$$$        REAL(dp), DIMENSION(mnmax) :: rmnc1, zmns1, lmns1,
c$$$     &   rmns1, zmnc1, lmnc1, bmodmn, bmodmn1
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubumnc_sur  !MRC 10-15-15
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubvmnc_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupumnc_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupvmnc_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubumns_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubvmns_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupumns_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupvmns_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsubua_sur, bsubva_sur
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: bsupua_sur, bsupva_sur
c$$$
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: xfinal
c$$$
c$$$!
c$$$!     THIS SUBROUTINE CREATES THE FILE WOUT.IT CONTAINS THE CYLINDRICAL COORDINATE SPECTRAL
c$$$!     COEFFICIENTS RMN,ZMN (full), LMN (half_mesh - CONVERTED FROM
c$$$!     INTERNAL full REPRESENTATION), AS WELL AS COEFFICIENTS (ON NYQ MESH) FOR COMPUTED
c$$$!     QUANTITIES:
c$$$!
c$$$!     BSQ, BSUPU,V, BSUBU,V, GSQRT (HALF); BSUBS (FULL-CONVERTED IN JXBFORCE)
c$$$!
c$$$        IF (lnyquist) THEN
c$$$           mnmax_nyq0 = mnmax_nyq
c$$$           mnyq0 = mnyq
c$$$           nnyq0 = nnyq
c$$$           xm_nyq0 => xm_nyq; xn_nyq0 => xn_nyq
c$$$        ELSE
c$$$           mnmax_nyq0 = mnmax
c$$$           mnyq0 = mpol1
c$$$           nnyq0 = ntor
c$$$           xm_nyq0 => xm; xn_nyq0 => xn
c$$$        END IF
c$$$
c$$$        ALLOCATE (gmn(mnmax_nyq0), bmn(mnmax_nyq0),
c$$$     &   bsubumn(mnmax_nyq0), bsubvmn(mnmax_nyq0), bsubsmn(mnmax_nyq0),
c$$$     &   bsupumn(mnmax_nyq0), bsupvmn(mnmax_nyq0), stat=istat) 
c$$$        ALLOCATE(rmnc(mnmax,ns),zmns(mnmax,ns),lmns(mnmax,ns),
c$$$     &   rmns(mnmax,ns),zmnc(mnmax,ns),lmnc(mnmax,ns))
c$$$
c$$$        IF (lfreeb) THEN        !MRC 10-15-15
c$$$           ALLOCATE 
c$$$     & (bsubua_sur(nzeta*ntheta2), bsubva_sur(nzeta*ntheta2))
c$$$           ALLOCATE 
c$$$     & (bsupua_sur(nzeta*ntheta2), bsupva_sur(nzeta*ntheta2))
c$$$
c$$$           ALLOCATE 
c$$$     & (bsubumnc_sur(mnmax_nyq0), bsubvmnc_sur(mnmax_nyq0))
c$$$           ALLOCATE 
c$$$     & (bsupumnc_sur(mnmax_nyq0), bsupvmnc_sur(mnmax_nyq0))
c$$$           IF (lasym) THEN
c$$$              ALLOCATE (bsubumns_sur(mnmax_nyq0),
c$$$     &                bsubvmns_sur(mnmax_nyq0))
c$$$              ALLOCATE (bsupumns_sur(mnmax_nyq0),
c$$$     &                bsupvmns_sur(mnmax_nyq0))
c$$$           END IF
c$$$        END IF
c$$$
c$$$        ALLOCATE (gmnc(mnmax_nyq0,ns), bmnc(mnmax_nyq0,ns),
c$$$     &          bsubumnc(mnmax_nyq0,ns), bsubvmnc(mnmax_nyq0,ns),
c$$$     &          bsubsmns(mnmax_nyq0,ns), bsupumnc(mnmax_nyq0,ns),
c$$$     &          bsupvmnc(mnmax_nyq0,ns), 
c$$$     &          currumnc(mnmax_nyq0,ns), currvmnc(mnmax_nyq0,ns))
c$$$
c$$$        IF (lasym) THEN
c$$$          ALLOCATE (gmns(mnmax_nyq0,ns), bmns(mnmax_nyq0,ns),
c$$$     &          bsubumns(mnmax_nyq0,ns), bsubvmns(mnmax_nyq0,ns),
c$$$     &          bsubsmnc(mnmax_nyq0,ns), bsupumns(mnmax_nyq0,ns),
c$$$     &          bsupvmns(mnmax_nyq0,ns),
c$$$     &          currumns(mnmax_nyq0,ns), currvmns(mnmax_nyq0,ns))
c$$$        END IF
c$$$        IF (istat .ne. 0) STOP 
c$$$     & 'Error allocating arrays in compute_output_representation'
c$$$        n1 = MAX(1,ntmax/2)
c$$$        rmnc = rzl_array(:,:,1)            !!store COS(mu-nv) components
c$$$        zmns = rzl_array(:,:,1+n1)         !!store SIN(mu-nv)
c$$$        lmns = rzl_array(:,:,1+2*n1)       !!store SIN(mu-nv)
c$$$
c$$$        IF (lasym) THEN
c$$$           rmns = gc_array(:,:,1)            !!store SIN(mu-nv)
c$$$           zmnc = gc_array(:,:,1+n1)         !!store COS(mu-nv)
c$$$           lmnc = gc_array(:,:,1+2*n1)       !!store COS(mu-nv)
c$$$        END IF
c$$$
c$$$
c$$$        indx1=MAXLOC(ns_array)
c$$$        ftolx1=ftol_array(indx1(1))
c$$$
c$$$!     NYQUIST FREQUENCY REQUIRES FACTOR OF 1/2
c$$$        IF (lnyquist) THEN
c$$$           IF (mnyq .ne. 0) cosmui(:,mnyq) = p5*cosmui(:,mnyq)
c$$$           IF (nnyq .ne. 0) cosnv (:,nnyq) = p5*cosnv (:,nnyq)
c$$$        END IF
c$$$
c$$$        ALLOCATE (xfinal(neqs), stat=js)
c$$$        IF (js .NE. 0) STOP 
c$$$     &  'Allocation error for xfinal in compute_output_representation'
c$$$        xfinal = xc
c$$$!
c$$$!     MUST CONVERT m=1 MODES... FROM INTERNAL TO PHYSICAL FORM
c$$$!     Extrapolation of m=0 Lambda (cs) modes, which are not evolved at j=1, done in CONVERT
c$$$!
c$$$        lk = ns*ntor1
c$$$        IF (lthreed) CALL convert_sym  (xfinal(1+mns*(rss-1)+lk), 
c$$$     &                                xfinal(1+irzloff+mns*(zcs-1)+lk))
c$$$        IF (lasym)   CALL convert_asym (xfinal(1+mns*(rsc-1)+lk), 
c$$$     &                                xfinal(1+irzloff+mns*(zcc-1)+lk))
c$$$
c$$$!
c$$$!     CONVERT TO rmnc, zmns, lmns, etc EXTERNAL representation (without internal mscale, nscale)
c$$$!     IF B^v ~ phip + lamu, MUST DIVIDE BY phipf(js) below to maintain old-style format
c$$$!     THIS COULD BE A PROBLEM FOR RFP WHERE PHIPF->0 INSIDE THE PLASMA!
c$$$!
c$$$        RADIUS1: DO js = 1, ns
c$$$
c$$$           CALL convert (rmnc1, zmns1, lmns1, rmns1, zmnc1, lmnc1, 
c$$$     &                         xfinal, js)
c$$$
c$$$           rmnc(:,js) = rmnc1(:)
c$$$           zmns(:,js) = zmns1(:)
c$$$           lmns(:,js) = (lmns1(:)/phipf(js)) * lamscale
c$$$           IF (lasym) THEN
c$$$              rmns(:,js) = rmns1(:)
c$$$              zmnc(:,js) = zmnc1(:)
c$$$              lmnc(:,js) = (lmnc1(:)/phipf(js)) * lamscale
c$$$           END IF
c$$$
c$$$        END DO RADIUS1
c$$$
c$$$        DEALLOCATE (xfinal)
c$$$
c$$$!
c$$$!     INTERPOLATE LAMBDA ONTO HALF-MESH FOR BACKWARDS CONSISTENCY WITH EARLIER VERSIONS OF VMEC
c$$$!     AND SMOOTHS POSSIBLE UNPHYSICAL "WIGGLE" ON RADIAL MESH
c$$$!
c$$$
c$$$        WHERE (NINT(xm) .le. 1) lmns(:,1) = lmns(:,2)
c$$$        DO js = ns,2,-1
c$$$           WHERE (MOD(NINT(xm),2) .eq. 0) 
c$$$              lmns(:,js) = p5*(lmns(:,js) + lmns(:,js-1))
c$$$           ELSEWHERE
c$$$              lmns(:,js) = 
c$$$     &  p5*(sm(js)*lmns(:,js) + sp(js-1)*lmns(:,js-1))
c$$$           END WHERE
c$$$        END DO
c$$$
c$$$        lmns(:,1) = 0  
c$$$        raxis_cc(0:ntor) = rmnc(1:ntor+1,1)
c$$$        zaxis_cs(0:ntor) = zmns(1:ntor+1,1)
c$$$        
c$$$        IF (.NOT.lasym) GOTO 900
c$$$
c$$$        WHERE (NINT(xm) .le. 1) lmnc(:,1) = lmnc(:,2)
c$$$        DO js = ns,2,-1
c$$$           WHERE (MOD(NINT(xm),2) .eq. 0) 
c$$$              lmnc(:,js) = p5*(lmnc(:,js) + lmnc(:,js-1))
c$$$           ELSEWHERE
c$$$              lmnc(:,js) = 
c$$$     &  p5*(sm(js)*lmnc(:,js) + sp(js-1)*lmnc(:,js-1))
c$$$           END WHERE
c$$$        END DO
c$$$
c$$$        lmnc(:,1) = 0;   
c$$$        raxis_cs(0:ntor) = rmns(1:ntor+1,1)
c$$$        zaxis_cc(0:ntor) = zmnc(1:ntor+1,1)
c$$$
c$$$ 900  CONTINUE
c$$$        DO js = 2, ns
c$$$           bsq(js,:nznt) = SQRT(2*ABS(bsq(js,:nznt)-pres(js)))
c$$$        END DO
c$$$
c$$$        tmult = p5/r0scale**2
c$$$!SPH: FIXED THIS 03-05-07 TO CALL symmetrization routine
c$$$        IF (lasym) THEN
c$$$!Changed integration norm in fixaray, SPH012314
c$$$           tmult = 2*tmult
c$$$           bsubs(1,:) = 0
c$$$           CALL symoutput (bsq,   gsqrt,  bsubu,  bsubv,  bsupu,
c$$$     &                   bsupv,  bsubs, 
c$$$     &                   bsqa,  gsqrta, bsubua, bsubva, bsupua,
c$$$     &                   bsupva, bsubsa)
c$$$
c$$$           IF (lfreeb) THEN     !MRC  10-15-15
c$$$              CALL symoutput_sur(bsubu_sur, bsubv_sur,
c$$$     &                         bsupu_sur, bsupv_sur,
c$$$     &                         bsubua_sur, bsubva_sur,
c$$$     &                         bsupua_sur, bsupva_sur)
c$$$           END IF
c$$$        END IF
c$$$
c$$$!         DO js = 2, ns
c$$$!            WRITE (200, *) 'JS: ', js, 'BSUBU, BSUBV'
c$$$!            WRITE (200, '(1p,6e12.4)') bsubu(js,:), bsubv(js,:)
c$$$!         END DO
c$$$
c$$$        RADIUS2: DO js = 2, ns
c$$$           gmn = 0
c$$$           bmn = 0
c$$$           bsubumn = 0
c$$$           bsubvmn = 0
c$$$           bsubsmn = 0
c$$$           bsupumn = 0
c$$$           bsupvmn = 0
c$$$
c$$$           MN2: DO mn = 1, mnmax_nyq0
c$$$              n = NINT(xn_nyq0(mn))/nfp
c$$$              m = NINT(xm_nyq0(mn))
c$$$              n1 = ABS(n)
c$$$              dmult = mscale(m)*nscale(n1)*tmult
c$$$              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
c$$$              sgn = SIGN(1, n)
c$$$              lk = 0
c$$$              DO j = 1, ntheta2
c$$$                 DO k = 1, nzeta
c$$$                    lk = lk + 1 
c$$$                    tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
c$$$     &                       sgn*sinmui(j,m)*sinnv(k,n1))
c$$$                    tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
c$$$     &                       sgn*cosmui(j,m)*sinnv(k,n1))
c$$$                    bmn(mn) = bmn(mn) + tcosi*bsq(js,lk)
c$$$                    gmn(mn) = gmn(mn) + tcosi*gsqrt(js,lk)
c$$$                    bsubumn(mn) = bsubumn(mn) + tcosi*bsubu(js,lk)
c$$$                    bsubvmn(mn) = bsubvmn(mn) + tcosi*bsubv(js,lk)
c$$$                    bsubsmn(mn) = bsubsmn(mn) + tsini*bsubs(js,lk)
c$$$                    bsupumn(mn) = bsupumn(mn) + tcosi*bsupu(js,lk) 
c$$$                    bsupvmn(mn) = bsupvmn(mn) + tcosi*bsupv(js,lk) 
c$$$                 END DO
c$$$              END DO
c$$$           END DO MN2
c$$$
c$$$           IF (js .eq. ns/2) bmodmn = bmn(1:mnmax)
c$$$           IF (js .eq. ns) bmodmn1 = bmn(1:mnmax)
c$$$           gmnc(:,js) = gmn(:)
c$$$           bmnc(:,js) = bmn(:)
c$$$           bsubumnc(:,js) = bsubumn(:)
c$$$           bsubvmnc(:,js) = bsubvmn(:)
c$$$           bsubsmns(:,js) = bsubsmn(:)
c$$$           bsupumnc(:,js) = bsupumn(:)
c$$$           bsupvmnc(:,js) = bsupvmn(:)
c$$$        END DO RADIUS2
c$$$
c$$$        IF (lfreeb) THEN    !MRC    10-15-15
c$$$           bsubumnc_sur = 0
c$$$           bsubvmnc_sur = 0
c$$$           bsupumnc_sur = 0
c$$$           bsupvmnc_sur = 0
c$$$           DO mn = 1, mnmax_nyq0
c$$$              n = NINT(xn_nyq0(mn))/nfp
c$$$              m = NINT(xm_nyq0(mn))
c$$$              n1 = ABS(n)
c$$$              dmult = mscale(m)*nscale(n1)*tmult
c$$$              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
c$$$              sgn = SIGN(1, n)
c$$$              lk = 0
c$$$              DO j = 1, ntheta2
c$$$                 DO k = 1, nzeta
c$$$                    lk = lk + 1
c$$$                    tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
c$$$     &                         sgn*sinmui(j,m)*sinnv(k,n1))
c$$$                    bsubumnc_sur(mn) = bsubumnc_sur(mn)
c$$$     &                               + tcosi*bsubu_sur(lk)
c$$$                    bsubvmnc_sur(mn) = bsubvmnc_sur(mn)
c$$$     &                               + tcosi*bsubv_sur(lk)
c$$$                    bsupumnc_sur(mn) = bsupumnc_sur(mn)
c$$$     &                               + tcosi*bsupu_sur(lk)
c$$$                    bsupvmnc_sur(mn) = bsupvmnc_sur(mn)
c$$$     &                               + tcosi*bsupv_sur(lk)
c$$$                 END DO
c$$$              END DO
c$$$           END DO
c$$$        END IF
c$$$
c$$$        gmnc(:,1) = 0; bmnc(:,1) = 0;
c$$$        bsubumnc(:,1) = 0
c$$$        bsubvmnc(:,1) = 0
c$$$        bsubsmns(:,1) = 2*bsubsmns(:,2) - bsubsmns(:,3)
c$$$        bsupumnc(:,1) = 0;  bsupvmnc(:,1) = 0
c$$$
c$$$        IF (.not.lasym) GO TO 200
c$$$
c$$$        RADIUS3: DO js = 2, ns
c$$$           gmn = 0
c$$$           bmn = 0
c$$$           bsubumn = 0
c$$$           bsubvmn = 0
c$$$           bsubsmn = 0
c$$$           bsupumn = 0
c$$$           bsupvmn = 0
c$$$
c$$$           MN3: DO mn = 1, mnmax_nyq0
c$$$              n = NINT(xn_nyq0(mn))/nfp
c$$$              m = NINT(xm_nyq0(mn))
c$$$              n1 = ABS(n)
c$$$              dmult = mscale(m)*nscale(n1)*tmult
c$$$              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
c$$$              sgn = SIGN(1, n)
c$$$              lk = 0
c$$$              jlk = js
c$$$              DO j = 1, ntheta2
c$$$                 DO k = 1, nzeta
c$$$                    lk = lk + 1
c$$$                    tcosi = dmult*(cosmui(j,m)*cosnv(k,n1) +
c$$$     &                         sgn*sinmui(j,m)*sinnv(k,n1))
c$$$                    tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
c$$$     &                         sgn*cosmui(j,m)*sinnv(k,n1))
c$$$                    bmn(mn) = bmn(mn) + tsini*bsqa(jlk)
c$$$                    gmn(mn) = gmn(mn) + tsini*gsqrta(jlk,0)
c$$$                    bsubumn(mn) = bsubumn(mn) + tsini*bsubua(jlk)
c$$$                    bsubvmn(mn) = bsubvmn(mn) + tsini*bsubva(jlk)
c$$$                    bsubsmn(mn) = bsubsmn(mn) + tcosi*bsubsa(jlk)
c$$$                    bsupumn(mn) = bsupumn(mn) + tsini*bsupua(jlk)
c$$$                    bsupvmn(mn) = bsupvmn(mn) + tsini*bsupva(jlk)
c$$$
c$$$                    jlk = jlk+ns
c$$$                 END DO
c$$$              END DO
c$$$           END DO MN3
c$$$     
c$$$           gmns(:,js) = gmn(:)
c$$$           bmns(:,js) = bmn(:)
c$$$           bsubumns(:,js) = bsubumn(:)
c$$$           bsubvmns(:,js) = bsubvmn(:)
c$$$           bsubsmnc(:,js) = bsubsmn(:)
c$$$           bsupumns(:,js) = bsupumn(:)
c$$$           bsupvmns(:,js) = bsupvmn(:)
c$$$        END DO RADIUS3
c$$$
c$$$        gmns(:,1) = 0; bmns(:,1) = 0
c$$$        bsubumns(:,1) = 0
c$$$        bsubvmns(:,1) = 0
c$$$        bsubsmnc(:,1) = 2*bsubsmnc(:,2) - bsubsmnc(:,3)
c$$$        bsupumns(:,1) = 0;  bsupvmns(:,1) = 0
c$$$
c$$$        IF (lfreeb) THEN        !MRC  10-15-15
c$$$           bsubumns_sur = 0
c$$$           bsubvmns_sur = 0
c$$$           bsupumns_sur = 0
c$$$           bsupvmns_sur = 0
c$$$
c$$$           DO mn = 1, mnmax_nyq0
c$$$              n = NINT(xn_nyq0(mn))/nfp
c$$$              m = NINT(xm_nyq0(mn))
c$$$              n1 = ABS(n)
c$$$              dmult = mscale(m)*nscale(n1)*tmult
c$$$              IF (m.eq.0 .or. n.eq.0) dmult = 2*dmult
c$$$              sgn = SIGN(1, n)
c$$$              lk = 0
c$$$              DO j = 1, ntheta2
c$$$                 DO k = 1, nzeta
c$$$                    lk = lk + 1
c$$$                    tsini = dmult*(sinmui(j,m)*cosnv(k,n1) -
c$$$     &                         sgn*cosmui(j,m)*sinnv(k,n1))
c$$$                    bsubumns_sur(mn) = bsubumns_sur(mn)
c$$$     &                               + tsini*bsubua_sur(lk)
c$$$                    bsubvmns_sur(mn) = bsubvmns_sur(mn)
c$$$     &                               + tsini*bsubva_sur(lk)
c$$$                    bsupumns_sur(mn) = bsupumns_sur(mn)
c$$$     &                               + tsini*bsupua_sur(lk)
c$$$                    bsupvmns_sur(mn) = bsupvmns_sur(mn)
c$$$     &                               + tsini*bsupva_sur(lk)
c$$$                 END DO
c$$$              END DO
c$$$           END DO
c$$$        END IF
c$$$
c$$$ 200  CONTINUE
c$$$
c$$$        CALL Compute_Currents(bsubsmnc, bsubsmns, bsubumnc, bsubumns,
c$$$     &                        bsubvmnc, bsubvmns,
c$$$     &                        xm_nyq0, xn_nyq0, mnmax_nyq0, lasym, ns,
c$$$     &                        currumnc, currvmnc, currumns, currvmns)
c$$$
c$$$
c$$$!     RESTORE nyq ENDPOINT VALUES
c$$$
c$$$        IF (lnyquist) THEN
c$$$           IF (mnyq .ne. 0) cosmui(:,mnyq) = 2*cosmui(:,mnyq)
c$$$           IF (nnyq .ne. 0) cosnv (:,nnyq) = 2*cosnv (:,nnyq)
c$$$        END IF
c$$$
c$$$        IF (ALLOCATED(bsubumnc_sur)) THEN
c$$$           DEALLOCATE(bsubumnc_sur, bsubvmnc_sur)
c$$$           DEALLOCATE(bsupumnc_sur, bsupvmnc_sur)
c$$$        END IF
c$$$        IF (ALLOCATED(bsubumns_sur)) THEN
c$$$           DEALLOCATE(bsubumns_sur, bsubvmns_sur)
c$$$           DEALLOCATE(bsupumns_sur, bsupvmns_sur)
c$$$        END IF
c$$$        IF (ALLOCATED(bsubua_sur)) THEN
c$$$           DEALLOCATE(bsubua_sur, bsubva_sur)
c$$$           DEALLOCATE(bsupua_sur, bsupva_sur)
c$$$        END IF
c$$$
c$$$        CALL freeb_data(rmnc1, zmns1, rmns1, zmnc1, bmodmn, bmodmn1)
c$$$
c$$$        rzl_array = 0
c$$$      END SUBROUTINE vmec_output_data
c$$$
c$$$      SUBROUTINE initialize(m_pol,n_tor)
c$$$!     &  BIND(C,name='initialize_vmec')
c$$$        USE vmec_input
c$$$        USE vparams
c$$$        INTEGER(C_INT), INTENT(IN) :: m_pol, n_tor
c$$$
c$$$        mpol = m_pol
c$$$        ntor = n_tor
c$$$        omp_num_threads = 8
c$$$        gamma = 0
c$$$        spres_ped = 1
c$$$        ntheta = 0;  nzeta = 0
c$$$        ns_array = 0;  ns_array(1) = ns_default
c$$$        niter_array = -1;
c$$$        bloat = 1
c$$$        rbc = 0;  rbs = 0; zbs = 0; zbc = 0
c$$$        time_slice = 0
c$$$        nfp = 1
c$$$        ncurr = 0
c$$$        nsin = ns_default
c$$$        niter = 100
c$$$        nstep = 10
c$$$        nvacskip = 1
c$$$        delt = 1
c$$$        ftol = 1.E-10_dp
c$$$        ftol_array = 0;  ftol_array(1) = ftol
c$$$        am = 0; ai = 0; ac = 0; aphi = 0; aphi(1) = 1
c$$$        pres_scale = 1
c$$$        raxis_cc = 0; zaxis_cs = 0; raxis_cs = 0; zaxis_cc = 0;
c$$$        mfilter_fbdy = -1; nfilter_fbdy = -1
c$$$        tcon0 = 1
c$$$        precon_type = 'NONE'; prec2d_threshold = 1.E-30_dp
c$$$        curtor = 0; 
c$$$        extcur = 0;  phiedge = 1;
c$$$        mgrid_file = 'NONE'
c$$$        trip3d_file = 'NONE' ! SAL - TRIP3D
c$$$        lfreeb = .true.
c$$$        lmove_axis = .true.
c$$$        lmac = .false.
c$$$        lforbal = .false.
c$$$        lasym = .false.
c$$$        lrfp = .false.
c$$$        loldout = .false.        ! J Geiger 2010-05-04 start
c$$$        ldiagno = .false.
c$$$        lgiveup = .false.        ! inserted M.Drevlak
c$$$        fgiveup = 3.E+01_dp      ! inserted M.Drevlak
c$$$        lbsubs = .false.         ! J Hanson. See jxbforce coding
c$$$        lfull3d1out = .false.
c$$$        lmovie = .false.         ! S Lazerson for making movie files
c$$$        lmoreiter = .false.      ! default value if no max_main_iterations given.
c$$$        max_main_iterations = 1  ! to keep a presumably expected standard behavior.
c$$$!DEC$ IF DEFINED (NETCDF)
c$$$        lwouttxt = .false.       ! to keep functionality as expected with netcdf
c$$$!DEC$ ELSE
c$$$        lwouttxt = .true.        ! and without netcdf
c$$$!DEC$ ENDIF
c$$$
c$$$        pcurr_type = 'power_series'
c$$$        piota_type = 'power_series'
c$$$        pmass_type = 'power_series'
c$$$
c$$$!     ANISTROPY PARAMETERS
c$$$        bcrit = 1
c$$$        at(0) = 1;  at(1:) = 0
c$$$        ah = 0
c$$$        ph_type = 'power_series'
c$$$        pt_type = 'power_series'
c$$$        ah_aux_s(:) = -1
c$$$        at_aux_s(:) = -1
c$$$        am_aux_s(:) = -1
c$$$        ac_aux_s(:) = -1
c$$$        ai_aux_s(:) = -1
c$$$      
c$$$
c$$$!
c$$$!     BACKWARDS COMPATIBILITY
c$$$!
c$$$        raxis = 0;  zaxis = 0
c$$$
c$$$      END SUBROUTINE initialize
c$$$
c$$$      SUBROUTINE finalize() !BIND(C,name='finalize_vmec')
c$$$        USE vmec_persistent, ONLY: xm, xn, xm_nyq, xn_nyq
c$$$        USE vmec_main
c$$$        IF (ALLOCATED(xm)) DEALLOCATE(xm,xn,xm_nyq,xn_nyq)
c$$$        IF (ALLOCATED(rmnc)) DEALLOCATE(rmnc,zmns,lmns,rmns,zmnc,lmnc)
c$$$        IF (ALLOCATED(gmnc)) DEALLOCATE(gmnc, bmnc, bsubumnc, bsubvmnc,
c$$$     &                                bsubsmns, bsupumnc, bsupvmnc,
c$$$     &                                currumnc, currvmnc)
c$$$        IF (ALLOCATED(gmns)) DEALLOCATE(gmns, bmns, bsubumns, bsubvmns,
c$$$     &                                  bsubsmnc, bsupumns, bsupvmns,
c$$$     &                                  currumns, currvmns)
c$$$! J Geiger: check also for allocation.
c$$$        IF (ALLOCATED(gmn)) DEALLOCATE (gmn, bmn, bsubumn, bsubvmn,
c$$$     &               bsubsmn, bsupumn, bsupvmn)
c$$$        IF (ALLOCATED(iotaf)) THEN
c$$$         DEALLOCATE(iotaf, mass, phi, presf, jcuru, jcurv, jdotb, buco,
c$$$     &              bvco, bucof, bvcof, chi,
c$$$     &              bdotgradv, equif, specw, tcon, psi, yellip, yinden,
c$$$     &              ytrian, yshift, ygeo, overr, faclam, iotas, phips,
c$$$     &              chips, pres, vp, beta_vol, jperp2, jpar2, bdotb,
c$$$     &              clam, blam, dlam, phipf, chipf)
c$$$        END IF
c$$$!        IF (ALLOCATED(cosmu))
c$$$!     &   DEALLOCATE(cosmu, sinmu, cosmum, sinmum, cosmui, cosmumi,
c$$$!     &             sinmui, sinmumi, cosnv, sinnv, cosnvn, sinnvn,
c$$$!     &             cosmui3, cosmumi3, cos01, sin01)
c$$$
c$$$!        IF (ALLOCATED(tanu))
c$$$!     &   DEALLOCATE(tanu, tanv, sinper, cosper, sinuv, cosuv, sinu,
c$$$!     &              cosu, sinv, cosv, sinui, cosui, cmns, csign, sinu1,
c$$$!     &              cosu1, sinv1, cosv1, imirr, xmpot, xnpot)
c$$$
c$$$      END SUBROUTINE finalize
c$$$
c$$$      SUBROUTINE set_vmec_data_real(data_size,data_id,data_ptr)
c$$$!     &  BIND(C,name='set_vmec_data_real')
c$$$        USE vmec_input
c$$$        USE vparams, ONLY: ntord, ndatafmax
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        REAL(C_DOUBLE), DIMENSION(data_size), INTENT(IN) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_string
c$$$        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: data_2d
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: data_1d
c$$$        INTEGER :: ix, ix1, ix2
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_string)
c$$$
c$$$        IF (TRIM(data_string) .eq. "rbc") THEN
c$$$          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
c$$$            DO ix = 1,data_size
c$$$              ix1 = MOD(ix-1,2*ntor+1) - ntor
c$$$              ix2 = (ix-1)/(2*ntor+1)
c$$$              data_2d(ix1, ix2) = data_ptr(ix)
c$$$            END DO
c$$$            rbc(-ntor:ntor,0:mpol-1) = data_2d
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "rbs") THEN
c$$$          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
c$$$            DO ix = 1,data_size
c$$$              ix1 = MOD(ix-1,2*ntor+1) - ntor
c$$$              ix2 = (ix-1)/(2*ntor+1)
c$$$              data_2d(ix1, ix2) = data_ptr(ix)
c$$$            END DO
c$$$            rbs(-ntor:ntor,0:mpol-1) = data_2d
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "zbc") THEN
c$$$          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
c$$$            DO ix = 1,data_size
c$$$              ix1 = MOD(ix-1,2*ntor+1) - ntor
c$$$              ix2 = (ix-1)/(2*ntor+1)
c$$$              data_2d(ix1, ix2) = data_ptr(ix)
c$$$            END DO
c$$$            zbc(-ntor:ntor,0:mpol-1) = data_2d
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "zbs") THEN
c$$$          ALLOCATE(data_2d(-ntor:ntor,0:mpol-1))
c$$$            DO ix = 1,data_size
c$$$              ix1 = MOD(ix-1,2*ntor+1) - ntor
c$$$              ix2 = (ix-1)/(2*ntor+1)
c$$$              data_2d(ix1, ix2) = data_ptr(ix)
c$$$            END DO
c$$$            zbs(-ntor:ntor,0:mpol-1) = data_2d
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "am_aux_s") THEN
c$$$          IF (data_size .gt. ndatafmax) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          am_aux_s(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "am_aux_f") THEN
c$$$          IF (data_size .gt. ndatafmax) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          am_aux_f(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ai_aux_s") THEN
c$$$          IF (data_size .gt. ndatafmax) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ai_aux_s(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ai_aux_f") THEN
c$$$          IF (data_size .gt. ndatafmax) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ai_aux_f(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ac_aux_s") THEN
c$$$          IF (data_size .gt. ndatafmax) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ac_aux_s(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ac_aux_f") THEN
c$$$          IF (data_size .gt. ndatafmax) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ndatafmax,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ac_aux_f(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ai") THEN
c$$$          IF (data_size .gt. 21) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$            WRITE(6,"(A)") "Must specify <= 21 expansion coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ai(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ac") THEN
c$$$          IF (data_size .gt. 21) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$            WRITE(6,"(A)") "Must specify <= 21 expansion coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ac(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "am") THEN
c$$$          IF (data_size .gt. 21) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$            WRITE(6,"(A)") "Must specify <= 21 expansion coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          am(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "aphi") THEN
c$$$          IF (data_size .gt. 20) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$            WRITE(6,"(A)") "Must specify <= 20 expansion coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          aphi(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "raxis") THEN
c$$$          IF (data_size .gt. ntord+1) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          raxis(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "zaxis") THEN
c$$$          IF (data_size .gt. ntord+1) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          zaxis(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "raxis_cc") THEN
c$$$          IF (data_size .gt. ntord+1) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          raxis_cc(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "raxis_cs") THEN
c$$$          IF (data_size .gt. ntord+1) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          raxis_cs(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "zaxis_cc") THEN
c$$$          IF (data_size .gt. ntord+1) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          zaxis_cc(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "zaxis_cs") THEN
c$$$          IF (data_size .gt. ntord+1) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= ",ntord,
c$$$     &        "coefficients!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          zaxis_cs(0:data_size-1) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "ftol_array") THEN
c$$$          IF (data_size .gt. 100) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC Input Error" 
c$$$            WRITE(6,"(A,2X,I6,2X,A)") "Must specify <= 100 coefficients"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ftol_array(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "delt") THEN
c$$$          delt = data_ptr(1) 
c$$$        ELSE IF (TRIM(data_string) .eq. "prec2d_threshold") THEN
c$$$          prec2d_threshold = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "phiedge") THEN
c$$$          phiedge = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "gamma") THEN
c$$$          gamma = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "bloat") THEN
c$$$          bloat = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "spres_ped") THEN
c$$$          spres_ped = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "pres_scale") THEN
c$$$          pres_scale = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "curtor") THEN
c$$$          curtor = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "tcon0") THEN
c$$$          tcon0 = data_ptr(1)
c$$$        ELSE
c$$$          WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$          WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_string)//" found!"
c$$$          STOP 
c$$$        END IF
c$$$      END SUBROUTINE
c$$$
c$$$      SUBROUTINE set_vmec_data_int(data_size,data_id,data_ptr)
c$$$!     &  BIND(C,name='set_vmec_data_int')
c$$$        USE vmec_input
c$$$        USE vparams
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        INTEGER(C_INT), DIMENSION(data_size), INTENT(IN) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_string
c$$$        INTEGER, DIMENSION(:,:), ALLOCATABLE :: data_2d
c$$$        INTEGER, DIMENSION(:), ALLOCATABLE :: data_1d
c$$$        INTEGER :: ix1, ix2
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_string)
c$$$
c$$$        IF (TRIM(data_string) .eq. "ns_array") THEN
c$$$          IF (data_size .gt. 100) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$            WRITE(6,"(A)") "ns_array cannot be larger than 100!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          ns_array(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "niter_array") THEN
c$$$          IF (data_size .gt. 100) THEN
c$$$            WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$            WRITE(6,"(A)") "niter_array cannot be larger than 100!"
c$$$            STOP
c$$$          END IF
c$$$          ALLOCATE(data_1d(data_size))
c$$$          data_1d = data_ptr
c$$$          niter_array(1:data_size) = data_1d
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "mpol") THEN
c$$$          mpol = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "ntor") THEN
c$$$          ntor = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "niter") THEN
c$$$          niter = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "nfp") THEN
c$$$          nfp = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "nsin") THEN
c$$$          nsin = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "nstep") THEN
c$$$          nstep = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "ncurr") THEN
c$$$          ncurr = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "nvacskip") THEN
c$$$          nvacskip = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "ntheta") THEN
c$$$          ntheta = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "nzeta") THEN
c$$$          nzeta = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "mfilter_fbdy") THEN
c$$$          mfilter_fbdy = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "nfilter_fbdy") THEN
c$$$          nfilter_fbdy = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "max_main_iterations") THEN
c$$$          max_main_iterations = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "omp_num_threads") THEN
c$$$          omp_num_threads = data_ptr(1)
c$$$        ELSE
c$$$          WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$          WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_string)//" found!"
c$$$          STOP 
c$$$        END IF
c$$$
c$$$      END SUBROUTINE set_vmec_data_int
c$$$
c$$$      SUBROUTINE set_vmec_data_bool(data_size,data_id,data_ptr)
c$$$!     &  BIND(C,name='set_vmec_data_bool')
c$$$        USE vmec_input
c$$$        USE vparams
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        LOGICAL(C_BOOL), DIMENSION(data_size), INTENT(IN) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_string
c$$$        LOGICAL, DIMENSION(:,:), ALLOCATABLE :: data_2d
c$$$        LOGICAL, DIMENSION(:), ALLOCATABLE :: data_1d
c$$$        INTEGER :: ix1, ix2
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_string)
c$$$
c$$$        IF (TRIM(data_string) .eq. "lasym") THEN
c$$$          lasym = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "lfreeb") THEN
c$$$          lfreeb = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "lwouttxt") THEN
c$$$          lwouttxt = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "lpofr") THEN
c$$$          lpofr = data_ptr(1)
c$$$        ELSE IF (TRIM(data_string) .eq. "full3d1out") THEN
c$$$          lfull3d1out = data_ptr(1)
c$$$        ELSE
c$$$          WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$          WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_string)//" found!"
c$$$          STOP 
c$$$        END IF
c$$$
c$$$      END SUBROUTINE set_vmec_data_bool
c$$$
c$$$      SUBROUTINE set_vmec_data_char(data_size,data_id,data_ptr)
c$$$!     &  BIND(C,name='set_vmec_data_char')
c$$$        USE vmec_input
c$$$        USE vparams
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        TYPE(C_PTR), INTENT(IN) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_name, data_string
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_name)
c$$$        CALL c2f_string_array_1d(data_ptr,data_string)
c$$$
c$$$        IF (TRIM(data_name) .eq. "precon_type") THEN
c$$$          precon_type = data_string 
c$$$        ELSE IF (TRIM(data_name) .eq. "pmass_type") THEN
c$$$          pmass_type = data_string
c$$$        ELSE IF (TRIM(data_name) .eq. "pcurr_type") THEN
c$$$          pcurr_type = data_string
c$$$        ELSE IF (TRIM(data_name) .eq. "piota_type") THEN
c$$$          piota_type = data_string
c$$$        ELSE IF (TRIM(data_name) .eq. "mgrid_file") THEN
c$$$          mgrid_file = data_string
c$$$        ELSE
c$$$          WRITE(6,"(A)") "Fatal VMEC input error!"
c$$$          WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_name)//" found!"
c$$$          STOP 
c$$$        END IF
c$$$
c$$$      END SUBROUTINE set_vmec_data_char
c$$$
c$$$      SUBROUTINE retrieve_vmec_data_real(data_size, data_id, data_ptr)
c$$$!     & BIND(C,name="retrieve_vmec_data_real")
c$$$        USE vmec_main
c$$$        USE vmec_persistent, ONLY: xm, xn, xm_nyq, xn_nyq
c$$$        USE vmec_io
c$$$        USE vmercier
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        REAL(C_DOUBLE), DIMENSION(data_size), INTENT(OUT) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_string
c$$$        REAL(dp), DIMENSION(:,:), ALLOCATABLE :: data_2d
c$$$        REAL(dp), DIMENSION(:), ALLOCATABLE :: data_1d
c$$$        INTEGER :: ix1, ix2, cix
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_string)
c$$$
c$$$        IF (TRIM(data_string) .eq. "rmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax,ns))
c$$$          data_2d = rmnc
c$$$          DO ix1 = 1,mnmax
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "rmns") THEN
c$$$          ALLOCATE(data_2d(mnmax,ns))
c$$$          data_2d = rmns
c$$$          DO ix1 = 1,mnmax
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "zmns") THEN
c$$$          ALLOCATE(data_2d(mnmax,ns))
c$$$          data_2d = zmns
c$$$          DO ix1 = 1,mnmax
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (trim(data_string) .eq. "zmnc") THEN
c$$$          ALLOCATE(data_2d(mnmax,ns))
c$$$          data_2d = zmnc
c$$$          DO ix1 = 1,mnmax
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "lmns") THEN
c$$$          ALLOCATE(data_2d(mnmax,ns))
c$$$          data_2d = lmns
c$$$          DO ix1 = 1,mnmax
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (trim(data_string) .eq. "lmnc") THEN
c$$$          ALLOCATE(data_2d(mnmax,0:ns))
c$$$          data_2d = lmnc
c$$$          DO ix1 = 1,mnmax
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "gmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = gmnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "gmns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = gmns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = ix2*mnmax_nyq + ix1-1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bmnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bmns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bmns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsubsmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsubsmnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsubsmns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsubsmns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsubumnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsubumnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsubumns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsubumns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsubvmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsubvmnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsubvmns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsubvmns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsupumnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsupumnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsupumns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsupumns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsupvmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsupvmnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bsupvmns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = bsupvmns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "currumnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = currumnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "currumns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = currumns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "currvmnc") THEN 
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = currvmnc
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "currvmns") THEN
c$$$          ALLOCATE(data_2d(mnmax_nyq,ns))
c$$$          data_2d = currvmns
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            DO ix2 = 1,ns 
c$$$              cix = (ix2-1)*mnmax_nyq + ix1;
c$$$              data_ptr(cix) = data_2d(ix1,ix2); 
c$$$            END DO 
c$$$          END DO
c$$$          DEALLOCATE(data_2d)
c$$$        ELSE IF (TRIM(data_string) .eq. "xm") THEN 
c$$$          ALLOCATE(data_1d(mnmax))
c$$$          data_1d = xm
c$$$          DO ix1 = 1,mnmax
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "xn") THEN
c$$$          ALLOCATE(data_1d(mnmax))
c$$$          data_1d = xn
c$$$          DO ix1 = 1,mnmax
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "xm_nyq") THEN
c$$$          ALLOCATE(data_1d(mnmax_nyq))
c$$$          data_1d = xm_nyq
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "xn_nyq") THEN
c$$$          ALLOCATE(data_1d(mnmax_nyq))
c$$$          data_1d = xn_nyq
c$$$          DO ix1 = 1,mnmax_nyq
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "phi") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = phi
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "phipf") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = phipf
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "phips") THEN
c$$$          ALLOCATE(data_1d(ns+1))
c$$$          data_1d = phips
c$$$          DO ix1 = 1,ns+1
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "iotaf") THEN
c$$$          ALLOCATE(data_1d(ns+1))
c$$$          data_1d = iotaf
c$$$          DO ix1 = 1,ns+1
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "iotas") THEN
c$$$          ALLOCATE(data_1d(ns+1))
c$$$          data_1d = iotas
c$$$          DO ix1 = 1,ns+1
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "pres") THEN
c$$$          ALLOCATE(data_1d(ns+1))
c$$$          data_1d = pres
c$$$          DO ix1 = 1,ns+1
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "presf") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = presf
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "DCurr") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = DCurr
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "DGeod") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = DGeod
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "DMerc") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = DMerc
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "DShear") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = DShear
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bdotb") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = bdotb
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "beta_vol") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = beta_vol
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "buco") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = buco
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "bvco") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = bvco
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "chi") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = chi
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "chipf") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = chipf
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "equif") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = equif
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "jcuru") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = jcuru
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "jcurv") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = jcurv
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "jdotb") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = jdotb
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "mass") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = mass
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "specw") THEN
c$$$          ALLOCATE(data_1d(ns))
c$$$          data_1d = specw
c$$$          DO ix1 = 1,ns
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "vp") THEN
c$$$          ALLOCATE(data_1d(ns+1))
c$$$          data_1d = vp
c$$$          DO ix1 = 1,ns+1
c$$$            data_ptr(ix1) = data_1d(ix1)
c$$$          END DO
c$$$          DEALLOCATE(data_1d)
c$$$        ELSE IF (TRIM(data_string) .eq. "Aminor_p") THEN
c$$$          data_ptr(1) = Aminor_p
c$$$        ELSE IF (TRIM(data_string) .eq. "Rmajor_p") THEN
c$$$          data_ptr(1) = Rmajor_p
c$$$        ELSE IF (TRIM(data_string) .eq. "IonLarmor") THEN
c$$$          data_ptr(1) = IonLarmor
c$$$        ELSE IF (TRIM(data_string) .eq. "aspect") THEN
c$$$          data_ptr(1) = aspect
c$$$        ELSE IF (TRIM(data_string) .eq. "b0") THEN
c$$$          data_ptr(1) = b0
c$$$        ELSE IF (TRIM(data_string) .eq. "betapol") THEN
c$$$          data_ptr(1) = betapol
c$$$        ELSE IF (TRIM(data_string) .eq. "betator") THEN
c$$$          data_ptr(1) = betator
c$$$        ELSE IF (TRIM(data_string) .eq. "betatotal") THEN
c$$$          data_ptr(1) = betatot
c$$$        ELSE IF (TRIM(data_string) .eq. "betaxis") THEN
c$$$          data_ptr(1) = betaxis
c$$$        ELSE IF (TRIM(data_string) .eq. "gamma") THEN
c$$$          data_ptr(1) = gamma
c$$$        ELSE IF (TRIM(data_string) .eq. "volavgB") THEN
c$$$          data_ptr(1) = volavgB
c$$$        ELSE IF (TRIM(data_string) .eq. "volume_p") THEN
c$$$          data_ptr(1) = volume_p
c$$$        ELSE IF (TRIM(data_string) .eq. "wp") THEN
c$$$          data_ptr(1) = wp
c$$$        ELSE
c$$$         WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_string)//" found!"
c$$$        END IF
c$$$
c$$$      END SUBROUTINE retrieve_vmec_data_real
c$$$
c$$$      SUBROUTINE retrieve_vmec_data_int(data_size, data_id, data_ptr)
c$$$!     &  BIND(C,name='retrieve_vmec_data_int')
c$$$        USE vmec_dim, ONLY: ns, mnmax
c$$$        USE vmec_main
c$$$        USE vmec_input
c$$$        USE vmec_params, ONLY: signgs
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        INTEGER(C_INT), DIMENSION(data_size), INTENT(OUT) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_string
c$$$        INTEGER :: ix1, ix2, cix
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_string)
c$$$
c$$$        IF (TRIM(data_string) .eq. "nfp") THEN
c$$$          data_ptr(1) = nfp
c$$$        ELSE IF (TRIM(data_string) .eq. "mpol") THEN
c$$$          data_ptr(1) = mpol
c$$$        ELSE IF (TRIM(data_string) .eq. "ntor") THEN
c$$$          data_ptr(1) = ntor
c$$$        ELSE IF (TRIM(data_string) .eq. "ns") THEN
c$$$          data_ptr(1) = ns
c$$$        ELSE IF (TRIM(data_string) .eq. "mnmax") THEN
c$$$          data_ptr(1) = mnmax
c$$$        ELSE IF (TRIM(data_string) .eq. "mnmax_nyq") THEN
c$$$          data_ptr(1) = mnmax_nyq
c$$$        ELSE IF (TRIM(data_string) .eq. "signgs") THEN
c$$$          data_ptr(1) = NINT(signgs)
c$$$          print *, data_ptr(1)
c$$$        ELSE
c$$$          WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_string)//" found!"
c$$$        END IF
c$$$
c$$$      END SUBROUTINE retrieve_vmec_data_int
c$$$
c$$$      SUBROUTINE retrieve_vmec_data_bool(data_size, data_id, data_ptr)
c$$$!     &  BIND(C,name='retrieve_vmec_data_bool')
c$$$        USE vmec_main
c$$$        INTEGER(C_INT), INTENT(IN) :: data_size
c$$$        TYPE(C_PTR), INTENT(IN) :: data_id
c$$$        LOGICAL(C_BOOL), DIMENSION(data_size), INTENT(OUT) :: data_ptr
c$$$        CHARACTER(LEN=:), ALLOCATABLE :: data_string
c$$$        INTEGER :: ix1, ix2, cix
c$$$
c$$$        CALL c2f_string_array_1d(data_id,data_string)
c$$$
c$$$        IF (TRIM(data_string) .eq. "lasym") THEN 
c$$$          data_ptr(1) = lasym
c$$$        ELSE IF (TRIM(data_string) .eq. "lfreeb") THEN
c$$$          data_ptr(1) = lfreeb
c$$$        ELSE IF (TRIM(data_string) .eq. "lrecon") THEN
c$$$          data_ptr(1) = lrecon
c$$$        ELSE IF (TRIM(data_string) .eq. "lrfp") THEN
c$$$          data_ptr(1) = lrfp
c$$$        ELSE 
c$$$          WRITE(6,"(A)") "No data field with label "
c$$$     &      //TRIM(data_string)//" found!"
c$$$        END IF
c$$$
c$$$      END SUBROUTINE retrieve_vmec_data_bool
c$$$
c$$$
c$$$      END MODULE vmec_ext_interface
