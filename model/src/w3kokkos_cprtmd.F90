!/ ------------------------------------------------------------------- /
!> @file w3kokkos_cprtmd.F90
!> @brief ww3-gpu wiring unit cprt: W3CPRT's spectral partitioning on the
!>        device copy of VA.
!>
!> Binds kokkos/src/fortran_iface/cprt_c.hpp of h0ffmann/ww3-gpu; PATCH.md
!> section 27 there is the caller side. Nothing is computed here but the packing
!> of what W3CPRT (w3iosfmd.F90) reads and the allocation of ICPRT and DTPRT it
!> makes. W3KOKKOS_CPRT replaces the host read of VA and CALL W3CPRT ( IMOD )
!> of W3WAVE: the Kokkos translation (iosf_w3cprt.hpp, part_w3part.hpp) reads
!> the resident VA, and only ICPRT and the FINSIZ columns of DTPRT come back.
!> It calls the Fortran W3CPRT whole, after a host read of VA, when a sea point
!> has a NaN depth (W3CPRT stops there, w3iosfmd.F90:270-274) or a spectrum
!> the translation cannot take (PTMEAN's out-of-bounds read, cprt_c.hpp).
!> Opt-in: WW_KOKKOS_CPRT=1; WW_KOKKOS_ALL does not turn it on.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_CPRTMD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  USE W3KOKKOSMD, ONLY: KOKKOS_PART, W3KOKKOS_REGISTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: KOKKOS_CPRT, W3KOKKOS_CPRT, W3KOKKOS_CPRT_SETUP
  !
  !> WW_KOKKOS_CPRT=1 on a build the translation covers (W3KOKKOS_CPRT_SETUP).
  LOGICAL :: KOKKOS_CPRT = .FALSE.
  !
  !> cprt_c.hpp: WW_CPRT_NIPAR, WW_CPRT_NRPAR.
  INTEGER, PARAMETER :: NIPAR = 13, NRPAR = 10
  !
  INTERFACE
    INTEGER(C_INT) FUNCTION WW_CPRT_VA ( IGRID, IPAR, NIPAR, RPAR,      &
         NRPAR, VA, NSPEC, NSEALM, MAPSF, MAPSTA, U10, U10D, ASF, DW,   &
         CG, WN, SIG, DSIP, DSII, TH, ECOS, ESIN, EWN1, ECG1, FINSIZ,   &
         NXPFULL, FALLBACK ) BIND(C, NAME='ww_cprt_va')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE        :: IGRID, NIPAR, NRPAR, NSPEC, NSEALM
      INTEGER(C_INT), INTENT(IN)   :: IPAR(*), MAPSF(*), MAPSTA(*)
      REAL(C_FLOAT), INTENT(IN)    :: RPAR(*), U10(*), U10D(*), ASF(*), &
           DW(*), CG(*), WN(*), SIG(*), DSIP(*), DSII(*), TH(*),        &
           ECOS(*), ESIN(*), EWN1(*), ECG1(*)
      REAL(C_FLOAT), INTENT(INOUT) :: VA(*)
      INTEGER(C_INT), INTENT(OUT)  :: FINSIZ, NXPFULL, FALLBACK
    END FUNCTION WW_CPRT_VA
    INTEGER(C_INT) FUNCTION WW_CPRT_FETCH ( IGRID, NSEALM, ICPRT,       &
         FINSIZ, DTPRT ) BIND(C, NAME='ww_cprt_fetch')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE         :: IGRID, NSEALM, FINSIZ
      INTEGER(C_INT), INTENT(INOUT) :: ICPRT(*)
      REAL(C_FLOAT), INTENT(INOUT)  :: DTPRT(*)
    END FUNCTION WW_CPRT_FETCH
  END INTERFACE
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> The flag. The translation covers W3CPRT as W3_SHRD compiles it, with or
  !> without W3_OMPG (W3CPRT has no OpenMP), and not W3_MPI (NSEAL /= NSEA,
  !> INIT_GET_ISEA) or W3_T, whose statements it does not hold.
  SUBROUTINE W3KOKKOS_CPRT_SETUP ( IMOD )
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    LOGICAL                 :: OK
    IF ( IMOD < 1 ) RETURN
    OK = .FALSE.
#ifdef W3_SHRD
    OK = .TRUE.
#endif
#if defined(W3_MPI) || defined(W3_T)
    OK = .FALSE.
#endif
    KOKKOS_CPRT = KOKKOS_PART ( 'CPRT', OK )
    CALL W3KOKKOS_REGISTER ( 'CPRT', KOKKOS_CPRT )
  END SUBROUTINE W3KOKKOS_CPRT_SETUP
  !/ ------------------------------------------------------------------- /
  !> CALL W3CPRT ( IMOD ) of W3WAVE, on Kokkos, or the Fortran W3CPRT after a
  !> host read of VA for the two cases above. Refused inside an OpenMP
  !> parallel region, as every VA sync.
  SUBROUTINE W3KOKKOS_CPRT ( IMOD )
    USE CONSTANTS, ONLY: UNDEF
    USE W3GDATMD, ONLY: NK, NTH, NSEA, NSEAL, MAPSF, MAPSTA, SIG, DSIP,  &
         DSII, TH, ECOS, ESIN, DTH, XFR, FACHFE, FTE, IGRID
    USE W3ADATMD, ONLY: WN, CG, U10, U10D, DW, NSEALM
    USE W3WDATMD, ONLY: VA, ASF
    USE W3ODATMD, ONLY: OUTPTS, O6INIT, ICPRT, DTPRT, DIMP, PTMETH,      &
         IHMAX, HSPMIN, WSMULT, WSCUT, PTFCUT, FLCOMB, IAPROC, NAPERR,   &
         NDSE
    USE W3DISPMD, ONLY: NAR1D, N1MAX, EWN1, ECG1, DSIE
    USE W3IOSFMD, ONLY: W3CPRT
    USE W3KOKKOS_RESIDENTMD, ONLY: W3KOKKOS_VA_HOST_READ
    USE W3SERVMD, ONLY: EXTCDE
!$  USE OMP_LIB, ONLY: OMP_IN_PARALLEL
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    INTEGER                 :: JSEA, ISEA, IX, IY, I, DIMXP
    INTEGER(C_INT)          :: IPAR(NIPAR), IERR, FINSIZ, NXPFULL,       &
         FALLBACK
    REAL(C_FLOAT)           :: RPAR(NRPAR)
    LOGICAL                 :: FORT
    !
!$  IF ( OMP_IN_PARALLEL() ) THEN
!$    WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_CPRTMD : '// &
!$         'VA SYNC CALLED INSIDE AN OPENMP PARALLEL REGION'
!$    CALL EXTCDE ( 1 )
!$  END IF
    !
    ! A NaN depth at a point W3CPRT partitions (w3iosfmd.F90:253, :260,
    ! :270-274): the Fortran W3CPRT writes its lines and stops.
    FORT = .FALSE.
    DO JSEA=1, NSEAL
      ISEA   = JSEA
      IX     = MAPSF(ISEA,1)
      IY     = MAPSF(ISEA,2)
      IF ( MAPSTA(IY,IX) .LT. 0 ) CYCLE
      IF ( DW(ISEA) .NE. DW(ISEA) ) FORT = .TRUE.
    END DO
    !
    ! The shapes the C side reads (cprt_c.hpp).
    IF ( SIZE(CG,1).NE.NK+2 .OR. SIZE(WN,1).NE.NK+2 .OR.                 &
         SIZE(CG,2).NE.NSEA+1 .OR. SIZE(WN,2).NE.NSEA+1 .OR.             &
         SIZE(DW).NE.NSEA+1 .OR. SIZE(U10).NE.NSEA .OR.                  &
         SIZE(U10D).NE.NSEA .OR. SIZE(MAPSF,1).NE.NSEA .OR.              &
         SIZE(SIG).NE.NK+2 .OR. SIZE(DSIP).NE.NK+2 ) THEN
      WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_CPRTMD : '// &
           'AN ARRAY W3CPRT READS HAS AN UNEXPECTED SHAPE'
      CALL EXTCDE ( 1 )
    END IF
    !
    IF ( .NOT. FORT ) THEN
      IPAR = 0
      IPAR(1)  = NK
      IPAR(2)  = NTH
      IPAR(3)  = NSEAL
      IPAR(4)  = NSEA
      IPAR(5)  = SIZE(MAPSTA,2)
      IPAR(6)  = SIZE(MAPSTA,1)
      IPAR(7)  = SIZE(ASF)
      IPAR(8)  = PTMETH
      IPAR(9)  = IHMAX
      IPAR(10) = MERGE ( 1, 0, FLCOMB )
      IPAR(11) = DIMP
      IPAR(12) = NAR1D
      IPAR(13) = N1MAX
      RPAR = (/ DTH, XFR, FACHFE, FTE, HSPMIN, WSMULT, WSCUT, PTFCUT,    &
           UNDEF, DSIE /)
      IERR = WW_CPRT_VA ( IGRID, IPAR, INT(NIPAR,C_INT), RPAR,          &
           INT(NRPAR,C_INT), VA, SIZE(VA,1), UBOUND(VA,2), MAPSF,        &
           MAPSTA, U10, U10D, ASF, DW, CG, WN, SIG, DSIP, DSII, TH,      &
           ECOS, ESIN, EWN1, ECG1, FINSIZ, NXPFULL, FALLBACK )
      IF ( IERR /= 0_C_INT ) THEN
        WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_CPRTMD : '// &
             'KOKKOS W3CPRT FAILED, CODE ', IERR
        CALL EXTCDE ( 1 )
      END IF
      FORT = FALLBACK /= 0_C_INT
    END IF
    !
    IF ( FORT ) THEN
      CALL W3KOKKOS_VA_HOST_READ
      CALL W3CPRT ( IMOD )
      RETURN
    END IF
    !
    ! W3CPRT's storage (w3iosfmd.F90:221-229, :320-326).
    IF ( O6INIT ) THEN
      DEALLOCATE ( OUTPTS(IMOD)%OUT6%DTPRT )
    ELSE
      ALLOCATE ( OUTPTS(IMOD)%OUT6%ICPRT(NSEALM+1,2) )
      ICPRT => OUTPTS(IMOD)%OUT6%ICPRT
      O6INIT = .TRUE.
    END IF
    !
    ! PTMEAN's format 1000 (w3partmd.F90:1415-1418, :1533-1534), once per
    ! PTMEAN call that met a full XP array, NPO+1 = DIMXP+1 each time.
    IF ( PTMETH .EQ. 4 .OR. PTMETH .EQ. 5 ) THEN
      DIMXP = 2
    ELSE
      DIMXP = ((NK+1)/2) * ((NTH-1)/2)
    END IF
    DO I=1, NXPFULL
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) DIMXP+1
    END DO
    !
    ALLOCATE ( OUTPTS(IMOD)%OUT6%DTPRT(DIMP,MAX(1,FINSIZ)) )
    DTPRT => OUTPTS(IMOD)%OUT6%DTPRT
    IERR = WW_CPRT_FETCH ( IGRID, INT(NSEALM,C_INT), ICPRT, FINSIZ, DTPRT )
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_CPRTMD : '// &
           'KOKKOS W3CPRT FETCH FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
    IF ( FINSIZ .LE. 0 ) DTPRT = 0.
    !
1000 FORMAT (/' *** WAVEWATCH III ERROR IN PTMEAN :'/                &
         '     XP ARRAY TOO SMALL AT PARTITION',I6/)
  END SUBROUTINE W3KOKKOS_CPRT
  !/ ------------------------------------------------------------------- /
END MODULE W3KOKKOS_CPRTMD
