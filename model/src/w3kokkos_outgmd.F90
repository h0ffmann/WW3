!/ ------------------------------------------------------------------- /
!> @file w3kokkos_outgmd.F90
!> @brief ww3-gpu wiring unit outg (WFIP-0005 unit 5): W3OUTG's integrals on
!>        the device copy of VA.
!>
!> Binds kokkos/src/fortran_iface/outg_c.hpp of h0ffmann/ww3-gpu; PATCH.md
!> section 15 there is the caller side. Nothing is computed here but the packing
!> of what W3OUTG (w3iogomd.F90) reads. W3KOKKOS_OUTG replaces CALL W3OUTG
!> (VA, FLPFLD, FLOUTG, FLOUTG2) of W3WAVE section 4.b: the Kokkos translation
!> (iogo_w3outg.hpp) reads the resident VA, so a field-output time copies no
!> VA to the host for it. When a field the translation leaves to the Fortran
!> is requested, the partitions of section 6 (FLPART), CALC_U3STOKES
!> (FLOLOC(6,8), (6,12)), SKEWNESS (FLOLOC(8,7:9)) or CALC_WBT (FLOLOC(2,17)),
!> or a requested 3-D field is not allocated, it calls the Fortran W3OUTG
!> whole, after a host read of VA. Opt-in: WW_KOKKOS_OUTG=1; WW_KOKKOS_ALL
!> does not turn it on.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_OUTGMD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT, C_PTR, C_LOC, C_NULL_PTR
  USE W3KOKKOSMD, ONLY: KOKKOS_PART, W3KOKKOS_REGISTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: KOKKOS_OUTG, W3KOKKOS_OUTG, W3KOKKOS_OUTG_SETUP
  !
  !> WW_KOKKOS_OUTG=1 on a build the translation covers (W3KOKKOS_OUTG_SETUP).
  LOGICAL :: KOKKOS_OUTG = .FALSE.
  !
  !> outg_c.hpp: WW_OUTG_NIPAR, WW_OUTG_NRPAR, WW_OUTG_NFIELD, WW_OUTG_N1D.
  INTEGER, PARAMETER :: NIPAR = 25, NRPAR = 10, NFIELD = 46, N1D = 40
  !
  INTERFACE
    INTEGER(C_INT) FUNCTION WW_OUTG_VA ( IGRID, IPAR, NIPAR, RPAR,      &
         NRPAR, FLOGRD, FLOGR2, FLOGRR, VA, NSPEC, NSEALM, CG, WN, DW,  &
         SIG, DDEN, DSII, ECOS, ESIN, EC2, ES2, ESC, MAPSF, MAPSTA,     &
         FIELDS, NFIELDS ) BIND(C, NAME='ww_outg_va')
      IMPORT :: C_INT, C_FLOAT, C_PTR
      INTEGER(C_INT), VALUE        :: IGRID, NIPAR, NRPAR, NSPEC,       &
           NSEALM, NFIELDS
      INTEGER(C_INT), INTENT(IN)   :: IPAR(*), FLOGRD(*), FLOGR2(*),    &
           FLOGRR(*), MAPSF(*), MAPSTA(*)
      REAL(C_FLOAT), INTENT(IN)    :: RPAR(*), CG(*), WN(*), DW(*),     &
           SIG(*), DDEN(*), DSII(*), ECOS(*), ESIN(*), EC2(*), ES2(*),  &
           ESC(*)
      REAL(C_FLOAT), INTENT(INOUT) :: VA(*)
      TYPE(C_PTR), INTENT(IN)      :: FIELDS(*)
    END FUNCTION WW_OUTG_VA
  END INTERFACE
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> The flag. The translation covers W3OUTG as W3_SHRD compiles it, with or
  !> without W3_OMPG, and none of MPI (NSEAL /= NSEA), IG1, O8, O9 or T, whose
  !> statements it does not hold.
  SUBROUTINE W3KOKKOS_OUTG_SETUP ( IMOD )
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    LOGICAL                 :: OK
    IF ( IMOD < 1 ) RETURN
    OK = .FALSE.
#ifdef W3_SHRD
    OK = .TRUE.
#endif
#if defined(W3_MPI) || defined(W3_IG1) || defined(W3_O8) || defined(W3_O9) || defined(W3_T)
    OK = .FALSE.
#endif
    KOKKOS_OUTG = KOKKOS_PART ( 'OUTG', OK )
    CALL W3KOKKOS_REGISTER ( 'OUTG', KOKKOS_OUTG )
  END SUBROUTINE W3KOKKOS_OUTG_SETUP
  !/ ------------------------------------------------------------------- /
  !> CALL W3OUTG ( VA, FLPART, FLOUTG, FLOUTG2 ) of W3WAVE, on Kokkos when the
  !> requested fields allow it, else the Fortran W3OUTG after a host read of VA.
  SUBROUTINE W3KOKKOS_OUTG ( FLPART, FLOUTG, FLOUTG2 )
    USE W3GDATMD, ONLY: NK, NTH, NSEAL, DTH, XFR, FTE, FTWL, FTTR,      &
         STEXU, STEYU, STEDU, E3DF, P2MSF, SIG, DDEN, DSII, ECOS, ESIN, &
         EC2, ES2, ESC, MAPSF, MAPSTA, IGRID
    USE W3WDATMD, ONLY: VA
    USE W3ADATMD, ONLY: CG, WN, DW, HS, WLM, T02, T0M1, T01, FP0, THM,  &
         THS, THP0, HSIG, QKK, QP, WBT, STMAXE, STMAXD, HCMAXE, HMAXE,  &
         HCMAXD, HMAXD, ABA, ABD, UBA, UBD, SXX, SYY, SXY, USSX, USSY,  &
         TUSX, TUSY, MSSX, MSSY, MSSD, MSCX, MSCY, MSCD, PRMS, TPMS,    &
         BHD, CGE, EF, TH1M, STH1M, TH2M, STH2M, P2SMS
    USE W3ODATMD, ONLY: UNDEF, WSCUT, FLOGRD, FLOGR2, FLOGRR, NOGRP,    &
         NGRPP, NDSE
    USE W3IORSMD, ONLY: OARST
    USE W3IOGOMD, ONLY: W3OUTG
    USE W3KOKKOS_RESIDENTMD, ONLY: W3KOKKOS_VA_HOST_READ
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    LOGICAL, INTENT(IN)     :: FLPART, FLOUTG, FLOUTG2
    LOGICAL                 :: FLOLOC(NOGRP,NGRPP), FORT
    INTEGER                 :: I, J, K
    INTEGER(C_INT)          :: IPAR(NIPAR), IERR
    INTEGER(C_INT)          :: LGRD(NOGRP*NGRPP), LGR2(NOGRP*NGRPP),    &
         LGRR(NOGRP*NGRPP)
    REAL(C_FLOAT)           :: RPAR(NRPAR)
    TYPE(C_PTR)             :: FPTR(NFIELD)
    !
    ! W3OUTG's FLOLOC (w3iogomd.F90:1381-1386).
    DO I=1, NOGRP
      DO J=1, NGRPP
        FLOLOC(I,J) = (FLOUTG.AND.FLOGRD(I,J)).OR.(FLOUTG2.AND.FLOGR2(I,J))
      END DO
    END DO
    FORT = FLPART .OR. FLOLOC(6,8) .OR. FLOLOC(6,12) .OR. FLOLOC(8,7)   &
         .OR. FLOLOC(8,8) .OR. FLOLOC(8,9) .OR. FLOLOC(2,17)
    DO K=1, 5
      FORT = FORT .OR. ( FLOLOC(3,K) .AND. E3DF(1,K).LE.0 )
    END DO
    FORT = FORT .OR. ( FLOLOC(6,9) .AND. P2MSF(1).LE.0 )
    IF ( FORT ) THEN
      CALL W3KOKKOS_VA_HOST_READ
      CALL W3OUTG ( VA, FLPART, FLOUTG, FLOUTG2 )
      RETURN
    END IF
    !
    DO J=1, NGRPP
      DO I=1, NOGRP
        K = I + (J-1)*NOGRP
        LGRD(K) = MERGE ( 1_C_INT, 0_C_INT, FLOGRD(I,J) )
        LGR2(K) = MERGE ( 1_C_INT, 0_C_INT, FLOGR2(I,J) )
        LGRR(K) = MERGE ( 1_C_INT, 0_C_INT, FLOGRR(I,J) )
      END DO
    END DO
    IPAR = 0
    IPAR(1)  = NK
    IPAR(2)  = NTH
    IPAR(3)  = NSEAL
    IPAR(4)  = SIZE(MAPSF,1)
    IPAR(5)  = SIZE(MAPSTA,2)
    IPAR(6)  = SIZE(MAPSTA,1)
    IPAR(7)  = SIZE(HS)
    IPAR(8)  = NOGRP
    IPAR(9)  = NGRPP
    IPAR(10) = MERGE ( 1, 0, FLOUTG )
    IPAR(11) = MERGE ( 1, 0, FLOUTG2 )
    IPAR(12) = MERGE ( 1, 0, OARST )
#ifdef W3_OMPG
    IPAR(13) = 1
#endif
    IPAR(14) = P2MSF(2)
    IPAR(15) = P2MSF(3)
    IPAR(16:20) = E3DF(2,1:5)
    IPAR(21:25) = E3DF(3,1:5)
    RPAR = (/ DTH, XFR, FTE, FTWL, FTTR, UNDEF, STEXU, STEYU, STEDU,    &
         WSCUT /)
    !
    ! CG and WN are (0:NK+1,0:NSEA), DW (0:NSEA): the C side reads them so.
    IF ( SIZE(CG,1).NE.NK+2 .OR. SIZE(WN,1).NE.NK+2 .OR.                 &
         SIZE(CG,2).NE.SIZE(MAPSF,1)+1 .OR. SIZE(WN,2).NE.SIZE(MAPSF,1)+1 &
         .OR. SIZE(DW).NE.SIZE(MAPSF,1)+1 ) THEN
      WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_OUTGMD : '// &
           'CG, WN OR DW IS NOT (0:NK+1,0:NSEA)'
      CALL EXTCDE ( 1 )
    END IF
    !
    ! The fields, in the order of WwOutgField (outg_c.hpp).
    FPTR = C_NULL_PTR
    FPTR(1)  = C_LOC(HS(1))
    FPTR(2)  = C_LOC(WLM(1))
    FPTR(3)  = C_LOC(T02(1))
    FPTR(4)  = C_LOC(T0M1(1))
    FPTR(5)  = C_LOC(T01(1))
    FPTR(6)  = C_LOC(FP0(1))
    FPTR(7)  = C_LOC(THM(1))
    FPTR(8)  = C_LOC(THS(1))
    FPTR(9)  = C_LOC(THP0(1))
    FPTR(10) = C_LOC(HSIG(1))
    FPTR(11) = C_LOC(QKK(1))
    FPTR(12) = C_LOC(QP(1))
    FPTR(13) = C_LOC(WBT(1))
    FPTR(14) = C_LOC(STMAXE(1))
    FPTR(15) = C_LOC(STMAXD(1))
    FPTR(16) = C_LOC(HCMAXE(1))
    FPTR(17) = C_LOC(HMAXE(1))
    FPTR(18) = C_LOC(HCMAXD(1))
    FPTR(19) = C_LOC(HMAXD(1))
    FPTR(20) = C_LOC(ABA(1))
    FPTR(21) = C_LOC(ABD(1))
    FPTR(22) = C_LOC(UBA(1))
    FPTR(23) = C_LOC(UBD(1))
    FPTR(24) = C_LOC(SXX(1))
    FPTR(25) = C_LOC(SYY(1))
    FPTR(26) = C_LOC(SXY(1))
    FPTR(27) = C_LOC(USSX(1))
    FPTR(28) = C_LOC(USSY(1))
    FPTR(29) = C_LOC(TUSX(1))
    FPTR(30) = C_LOC(TUSY(1))
    FPTR(31) = C_LOC(MSSX(1))
    FPTR(32) = C_LOC(MSSY(1))
    FPTR(33) = C_LOC(MSSD(1))
    FPTR(34) = C_LOC(MSCX(1))
    FPTR(35) = C_LOC(MSCY(1))
    FPTR(36) = C_LOC(MSCD(1))
    FPTR(37) = C_LOC(PRMS(1))
    FPTR(38) = C_LOC(TPMS(1))
    FPTR(39) = C_LOC(BHD(1))
    FPTR(40) = C_LOC(CGE(1))
    ! W3DIMA allocates a 3-D field only when E3DF(1,K) > 0 (w3adatmd.F90:1072-1101).
    IF ( E3DF(1,1).GT.0 ) FPTR(N1D+1) = C_LOC(EF(1,E3DF(2,1)))
    IF ( E3DF(1,2).GT.0 ) FPTR(N1D+2) = C_LOC(TH1M(1,E3DF(2,2)))
    IF ( E3DF(1,3).GT.0 ) FPTR(N1D+3) = C_LOC(STH1M(1,E3DF(2,3)))
    IF ( E3DF(1,4).GT.0 ) FPTR(N1D+4) = C_LOC(TH2M(1,E3DF(2,4)))
    IF ( E3DF(1,5).GT.0 ) FPTR(N1D+5) = C_LOC(STH2M(1,E3DF(2,5)))
    IF ( P2MSF(1).GT.0 )  FPTR(N1D+6) = C_LOC(P2SMS(1,P2MSF(2)))
    !
    IERR = WW_OUTG_VA ( IGRID, IPAR, INT(NIPAR,C_INT), RPAR,             &
         INT(NRPAR,C_INT), LGRD, LGR2, LGRR, VA, SIZE(VA,1),             &
         UBOUND(VA,2), CG, WN, DW, SIG, DDEN, DSII, ECOS, ESIN, EC2, ES2,&
         ESC, MAPSF, MAPSTA, FPTR, INT(NFIELD,C_INT) )
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_OUTGMD : '// &
           'KOKKOS W3OUTG FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
  END SUBROUTINE W3KOKKOS_OUTG
  !/ ------------------------------------------------------------------- /
END MODULE W3KOKKOS_OUTGMD
