!/ ------------------------------------------------------------------- /
!> @file w3kokkos_prop3md.F90
!> @brief ww3-gpu wiring unit prop3: the PR3 propagation (W3XYP3, W3KTP3, with
!>        W3QCK1, W3QCK2, W3QCK3 of UQ) on Kokkos, one call per loop of W3WAVE.
!>
!> Binds kokkos/src/fortran_iface/prop3_c.hpp of h0ffmann/ww3-gpu; PATCH.md
!> section 11 there is the caller side. Nothing is computed here but the packing
!> of what the Fortran routines read: W3KOKKOS_XYP3 replaces the DO ISPEC loop of
!> W3WAVE section 3.6.3 (W3GATH, W3XYP3, W3SCAT per component, W3_SHRD), and
!> W3KOKKOS_KTP3 one intra-spectral part, the DO ITLOC loop of section 3.6.2 or
!> 3.6.4 with its DO JSEA loop. Each part is opt-in: WW_KOKKOS_XYP3=1 and
!> WW_KOKKOS_KTP3=1; WW_KOKKOS_ALL does not turn them on.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_PROP3MD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  USE W3KOKKOSMD, ONLY: KOKKOS_PART, W3KOKKOS_REGISTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: KOKKOS_XYP3, KOKKOS_KTP3, W3KOKKOS_XYP3, W3KOKKOS_KTP3,      &
       W3KOKKOS_PROP3_SETUP
  !
  !> WW_KOKKOS_XYP3=1, WW_KOKKOS_KTP3=1 on a grid the translation covers
  !> (W3KOKKOS_PROP3_SETUP).
  LOGICAL :: KOKKOS_XYP3 = .FALSE.
  LOGICAL :: KOKKOS_KTP3 = .FALSE.
  !
  INTERFACE
    INTEGER(C_INT) FUNCTION WW_XYP3_BATCH ( RPAR, NRPAR, IPAR, NIPAR,   &
         IMAPS, NIMAPS, TSCAL, NTSCAL, DTG, MAPSF, CLATS, SIG, ECOS,     &
         ESIN, DPDX, DPDY, DQDX, DQDY, GSQRT, ATRNX, ATRNY, MAPX2,       &
         MAPY2, MAPAXY, MAPCXY, MAPTRN, MAPSTA, CG, CX, CY, FLBPI, NBI,  &
         ISBPI, BBPI0, BBPIN, ON, NCOL, VA )                             &
         BIND(C, NAME='ww_xyp3_batch')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE        :: NRPAR, NIPAR, NIMAPS, NTSCAL,     &
           FLBPI, NBI, NCOL
      REAL(C_FLOAT), VALUE         :: DTG
      REAL(C_FLOAT), INTENT(IN)    :: RPAR(*), TSCAL(*), CLATS(*),      &
           SIG(*), ECOS(*), ESIN(*), DPDX(*), DPDY(*), DQDX(*), DQDY(*), &
           GSQRT(*), ATRNX(*), ATRNY(*), CG(*), CX(*), CY(*), BBPI0(*),  &
           BBPIN(*)
      INTEGER(C_INT), INTENT(IN)   :: IPAR(*), IMAPS(*), MAPSF(*),      &
           MAPX2(*), MAPY2(*), MAPAXY(*), MAPCXY(*), MAPTRN(*),          &
           MAPSTA(*), ISBPI(*), ON(*)
      REAL(C_FLOAT), INTENT(INOUT) :: VA(*)
    END FUNCTION WW_XYP3_BATCH
    INTEGER(C_INT) FUNCTION WW_KTP3_BATCH ( RPAR, NRPAR, IPAR, NIPAR,   &
         SIG, DSIP, ECOS, ESIN, EC2, ESC, ES2, MAPWN, MAPTH2, MAPWN2,    &
         ITSTEP, NTIMES, FACTH, FACK, NSEA, CG, WN, NPTS, JSEA, ISEA,    &
         PSCAL, NCOL, VA, CFLTHMAX, CFLKMAX )                            &
         BIND(C, NAME='ww_ktp3_batch')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE        :: NRPAR, NIPAR, ITSTEP, NTIMES,     &
           NSEA, NPTS, NCOL
      REAL(C_FLOAT), VALUE         :: FACTH, FACK
      REAL(C_FLOAT), INTENT(IN)    :: RPAR(*), SIG(*), DSIP(*),         &
           ECOS(*), ESIN(*), EC2(*), ESC(*), ES2(*), CG(*), WN(*),       &
           PSCAL(*)
      INTEGER(C_INT), INTENT(IN)   :: IPAR(*), MAPWN(*), MAPTH2(*),     &
           MAPWN2(*), JSEA(*), ISEA(*)
      REAL(C_FLOAT), INTENT(INOUT) :: VA(*), CFLTHMAX(*), CFLKMAX(*)
    END FUNCTION WW_KTP3_BATCH
  END INTERFACE
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> The flags. XYP3 needs W3_SHRD (the W3GATH/W3SCAT it replaces are the
  !> shared-memory ones) and a grid that is not tripole; both need PR3 with UQ
  !> and a regular or curvilinear grid, and none of PDLIB, MGP, MGG (moving grid)
  !> or REFRX, which the translation does not cover.
  SUBROUTINE W3KOKKOS_PROP3_SETUP ( IMOD )
#if defined(W3_PR3) && defined(W3_UQ)
    USE W3GDATMD, ONLY: GRIDS, RLGTYPE, CLGTYPE, ICLOSE_TRPL
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
#if defined(W3_PR3) && defined(W3_UQ)
    LOGICAL                 :: OKG, OKX, OKK
#endif
    IF ( IMOD < 1 ) RETURN
#if defined(W3_PR3) && defined(W3_UQ)
    OKG = GRIDS(IMOD)%GTYPE .EQ. RLGTYPE .OR. GRIDS(IMOD)%GTYPE .EQ. CLGTYPE
#if defined(W3_PDLIB) || defined(W3_MGP) || defined(W3_MGG)
    OKG = .FALSE.
#endif
    OKX = OKG .AND. GRIDS(IMOD)%ICLOSE .NE. ICLOSE_TRPL
#ifndef W3_SHRD
    OKX = .FALSE.
#endif
    OKK = OKG
#ifdef W3_REFRX
    OKK = .FALSE.
#endif
    KOKKOS_XYP3 = KOKKOS_PART ( 'XYP3', OKX )
    KOKKOS_KTP3 = KOKKOS_PART ( 'KTP3', OKK )
    CALL W3KOKKOS_REGISTER ( 'XYP3', KOKKOS_XYP3 )
    CALL W3KOKKOS_REGISTER ( 'KTP3', KOKKOS_KTP3 )
#endif
  END SUBROUTINE W3KOKKOS_PROP3_SETUP
  !/ ------------------------------------------------------------------- /
  !> W3WAVE section 3.6.3 on Kokkos: for every ISPEC with IAPPRO(ISPEC) ==
  !> IAPROC, W3GATH, W3XYP3 (ISPEC, DTG, MAPSTA, MAPFS, FIELD, VGX, VGY) and
  !> W3SCAT, W3_SHRD. MAPFS is not read by W3XYP3, VGX and VGY only under MGP.
  SUBROUTINE W3KOKKOS_XYP3 ( DTG, MAPSTA )
#if defined(W3_PR3) && defined(W3_UQ)
    USE W3GDATMD, ONLY: NX, NY, NSEA, NK, NTH, NSPEC, MAPSF, DTCFL,     &
         CLATS, ICLOSE, FLCX, FLCY, DTH, XFR, ECOS, ESIN, SIG, WDCG,     &
         WDTH, FLAGLL, DPDX, DPDY, DQDX, DQDY, GSQRT, GTYPE, RLGTYPE,    &
         CLGTYPE, ICLOSE_TRPL
    USE W3WDATMD, ONLY: TIME, VA
    USE W3ADATMD, ONLY: NMX0, NMX1, NMX2, NMY0, NMY1, NMY2, NACT,       &
         NCENT, MAPX2, MAPY2, MAPAXY, MAPCXY, MAPTRN, CG, CX, CY,        &
         ATRNX, ATRNY, IAPPRO
    USE W3IDATMD, ONLY: FLCUR
    USE W3ODATMD, ONLY: NDSE, IAPROC, FLBPI, NBI, TBPI0, TBPIN, ISBPI,  &
         BBPI0, BBPIN
    USE W3TIMEMD, ONLY: DSEC21
    USE W3SERVMD, ONLY: EXTCDE
#endif
    IMPLICIT NONE
    REAL, INTENT(IN)        :: DTG
    INTEGER, INTENT(IN)     :: MAPSTA(*)
#if defined(W3_PR3) && defined(W3_UQ)
    REAL(C_FLOAT)           :: RPAR(5), TSCAL(3), RDUM(1)
    INTEGER(C_INT)          :: IPAR(11), IMAPS(8), IDUM(1), IERR
    INTEGER(C_INT), ALLOCATABLE, SAVE :: MTRN(:), ON(:)
    INTEGER                 :: TTEST(2), ISP
    !
    ! The set-up checked the grid; a later grid of a multi-grid run may differ.
    IF ( SIZE(MAPSF,1) .NE. NSEA .OR. ICLOSE .EQ. ICLOSE_TRPL .OR.       &
         ( GTYPE .NE. RLGTYPE .AND. GTYPE .NE. CLGTYPE ) ) THEN
      WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PROP3MD : '// &
           'W3XYP3 ON KOKKOS NEEDS A REGULAR OR CURVILINEAR, NOT '//     &
           'TRIPOLE, GRID AND MAPSF WITH NSEA ROWS'
      CALL EXTCDE ( 1 )
    END IF
    IF ( ALLOCATED(MTRN) ) THEN
      IF ( SIZE(MTRN) .NE. NX*NY .OR. SIZE(ON) .NE. NSPEC )             &
           DEALLOCATE ( MTRN, ON )
    END IF
    IF ( .NOT. ALLOCATED(MTRN) ) ALLOCATE ( MTRN(NX*NY), ON(NSPEC) )
    MTRN = MERGE ( 1_C_INT, 0_C_INT, MAPTRN(1:NX*NY) )
    DO ISP = 1, NSPEC
      ON(ISP) = MERGE ( 1_C_INT, 0_C_INT, IAPPRO(ISP) .EQ. IAPROC )
    END DO
    !
    RPAR = (/ DTCFL, DTH, XFR, WDCG, WDTH /)
    IPAR = (/ NX, NY, NSEA, NK, NTH, NSPEC, ICLOSE,                      &
         MERGE(1,0,FLCX), MERGE(1,0,FLCY), MERGE(1,0,FLAGLL),            &
         MERGE(1,0,FLCUR) /)
    IMAPS = (/ NMX0, NMX1, NMX2, NMY0, NMY1, NMY2, NACT, NCENT /)
    ! W3XYP3's DSEC21 calls (w3pro3md.F90:928-930, :1396-1398), REAL as there.
    TTEST(1) = TIME(1)
    TTEST(2) = 0
    TSCAL = 0.
    TSCAL(1) = DSEC21 ( TTEST, TIME )
    IF ( FLBPI ) THEN
      TSCAL(2) = DSEC21 ( TBPI0, TIME )
      TSCAL(3) = DSEC21 ( TBPI0, TBPIN )
      IERR = WW_XYP3_BATCH ( RPAR, 5_C_INT, IPAR, 11_C_INT, IMAPS,      &
           8_C_INT, TSCAL, 3_C_INT, DTG, MAPSF, CLATS, SIG, ECOS, ESIN,  &
           DPDX, DPDY, DQDX, DQDY, GSQRT, ATRNX, ATRNY, MAPX2, MAPY2,    &
           MAPAXY, MAPCXY, MTRN, MAPSTA, CG, CX, CY, 1_C_INT, NBI,       &
           ISBPI, BBPI0, BBPIN, ON, SIZE(VA,2), VA )
    ELSE
      IDUM = 0
      RDUM = 0.
      IERR = WW_XYP3_BATCH ( RPAR, 5_C_INT, IPAR, 11_C_INT, IMAPS,      &
           8_C_INT, TSCAL, 3_C_INT, DTG, MAPSF, CLATS, SIG, ECOS, ESIN,  &
           DPDX, DPDY, DQDX, DQDY, GSQRT, ATRNX, ATRNY, MAPX2, MAPY2,    &
           MAPAXY, MAPCXY, MTRN, MAPSTA, CG, CX, CY, 0_C_INT, 0_C_INT,   &
           IDUM, RDUM, RDUM, ON, SIZE(VA,2), VA )
    END IF
    CALL PROP3_CHECK ( 'W3XYP3', IERR )
#else
    CALL PROP3_CHECK ( 'W3XYP3', 5_C_INT )
#endif
  END SUBROUTINE W3KOKKOS_XYP3
  !/ ------------------------------------------------------------------- /
  !> W3WAVE section 3.6.2 or 3.6.4 on Kokkos: NTIMES sub-steps of W3KTP3
  !> (ISEA, FACTH, FACK, CTHG0S(ISEA), CG(:,ISEA), WN(:,ISEA), DEPTH,
  !> DDDX(IY,IX), ..., VA(:,JSEA), CFLTHMAX(JSEA), CFLKMAX(JSEA)) on every JSEA
  !> with MAPSTA(IY,IX) == 1, the arguments packed as W3WAVE passes them.
  SUBROUTINE W3KOKKOS_KTP3 ( FACTH, FACK, NTIMES )
#if defined(W3_PR3) && defined(W3_UQ)
    USE W3GDATMD, ONLY: NK, NK2, NTH, NSPEC, NSEA, NSEAL, SIG, DSIP,    &
         ECOS, ESIN, EC2, ESC, ES2, FACHFA, MAPWN, FLCTH, FLCK, CTMAX,   &
         DMIN, MAPSF, MAPSTA, CTHG0S, GTYPE, RLGTYPE, CLGTYPE
    USE W3WDATMD, ONLY: VA
    USE W3ADATMD, ONLY: MAPTH2, MAPWN2, ITSTEP, CG, WN, DW, CX, CY,     &
         DDDX, DDDY, DCXDX, DCXDY, DCYDX, DCYDY, CFLTHMAX, CFLKMAX
    USE W3IDATMD, ONLY: FLCUR
    USE W3PARALL, ONLY: INIT_GET_ISEA
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#endif
    IMPLICIT NONE
    REAL, INTENT(IN)        :: FACTH, FACK
    INTEGER, INTENT(IN)     :: NTIMES
#if defined(W3_PR3) && defined(W3_UQ)
    REAL(C_FLOAT)           :: RPAR(3)
    INTEGER(C_INT)          :: IPAR(7), IERR
    INTEGER(C_INT), ALLOCATABLE, SAVE :: PJSEA(:), PISEA(:)
    REAL(C_FLOAT), ALLOCATABLE, SAVE  :: PS(:,:)
    INTEGER                 :: JSEA, ISEA, IX, IY, NP
    !
    IF ( NTIMES .LT. 1 ) RETURN
    IF ( GTYPE .NE. RLGTYPE .AND. GTYPE .NE. CLGTYPE ) THEN
      WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PROP3MD : '// &
           'W3KTP3 ON KOKKOS NEEDS A REGULAR OR CURVILINEAR GRID'
      CALL EXTCDE ( 1 )
    END IF
    IF ( ALLOCATED(PJSEA) ) THEN
      IF ( SIZE(PJSEA) .NE. NSEAL ) DEALLOCATE ( PJSEA, PISEA, PS )
    END IF
    IF ( .NOT. ALLOCATED(PJSEA) )                                       &
         ALLOCATE ( PJSEA(NSEAL), PISEA(NSEAL), PS(10,NSEAL) )
    NP = 0
    DO JSEA = 1, NSEAL
      CALL INIT_GET_ISEA(ISEA, JSEA)
      IX     = MAPSF(ISEA,1)
      IY     = MAPSF(ISEA,2)
      IF ( MAPSTA(IY,IX) .EQ. 1 ) THEN
        NP = NP + 1
        PJSEA(NP) = JSEA
        PISEA(NP) = ISEA
        PS(1,NP) = CTHG0S(ISEA)
        PS(2,NP) = MAX ( DMIN , DW(ISEA) )
        PS(3,NP) = DDDX(IY,IX)
        PS(4,NP) = DDDY(IY,IX)
        PS(5,NP) = CX(ISEA)
        PS(6,NP) = CY(ISEA)
        PS(7,NP) = DCXDX(IY,IX)
        PS(8,NP) = DCXDY(IY,IX)
        PS(9,NP) = DCYDX(IY,IX)
        PS(10,NP) = DCYDY(IY,IX)
      END IF
    END DO
    !
    RPAR = (/ DMIN, CTMAX, FACHFA /)
    IPAR = (/ NK, NTH, NSPEC, NK2, MERGE(1,0,FLCTH), MERGE(1,0,FLCK),    &
         MERGE(1,0,FLCUR) /)
    IERR = WW_KTP3_BATCH ( RPAR, 3_C_INT, IPAR, 7_C_INT, SIG, DSIP,     &
         ECOS, ESIN, EC2, ESC, ES2, MAPWN, MAPTH2, MAPWN2, ITSTEP,       &
         NTIMES, FACTH, FACK, NSEA, CG, WN, NP, PJSEA, PISEA, PS,       &
         SIZE(VA,2), VA, CFLTHMAX, CFLKMAX )
    CALL PROP3_CHECK ( 'W3KTP3', IERR )
#else
    CALL PROP3_CHECK ( 'W3KTP3', 5_C_INT )
#endif
  END SUBROUTINE W3KOKKOS_KTP3
  !/ ------------------------------------------------------------------- /
  SUBROUTINE PROP3_CHECK ( WHAT, IERR )
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WHAT
    INTEGER(C_INT), INTENT(IN)   :: IERR
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PROP3MD : '// &
           'KOKKOS '//WHAT//' FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
  END SUBROUTINE PROP3_CHECK
  !/ ------------------------------------------------------------------- /
END MODULE W3KOKKOS_PROP3MD
