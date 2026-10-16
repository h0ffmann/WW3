!/ ------------------------------------------------------------------- /
!> @file w3kokkosmd.F90
!> @brief WAVEWATCH III's interface to the C++/Kokkos source-term port.
!>
!> Binds kokkos/src/fortran_iface/ww_kokkos_c.hpp, one INTERFACE block per
!> declaration in that header. Nothing is computed here: the module is a contract,
!> plus the one switch the model reads.
!>
!> Written to be dropped into WW3/model/src unchanged -- see
!> kokkos/src/fortran_iface/PATCH.md for the caller side (w3srcemd.F90, w3initmd.F90,
!> src_list.cmake, switches.json). Hence the WW3 house style, the LGPL header and
!> the absence of any dependency on the rest of this lab.
!>
!> Calling order:
!>
!>     IERR = WW_KOKKOS_INIT ( -1 )          ! once per process, W3INIT
!>     CALL W3KOKKOS_SETUP                   ! sets KOKKOS_SNL1
!>     IERR = WW_SNL1_INIT ( NK, NTH, ..., SIG(1) )
!>                                           ! once per process -- ONE grid; note
!>                                           ! SIG(1), not SIG (see remark 4)
!>     CALL WW_SNL1 ( NPTS, A, CG, KDMEAN, S, D )   ! per source-term call
!>     ...
!>     CALL WW_KOKKOS_FINALIZE               ! once, at W3WAVE teardown
!>
!> Error handling: WW_SNL1 is a SUBROUTINE because a Fortran CALL cannot inspect a
!> return value; its status is WW_SNL1_LAST_ERROR(), 0 for success. The C++ side
!> guarantees that no exception ever crosses the boundary, so a non-zero code is
!> the only way a failure is reported and it must not be ignored.
!>
!> Precision: every REAL here is REAL(C_FLOAT), matching WW3's default REAL. A WW3
!> built with an 8-byte default REAL needs a different ABI, not a KIND change.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!>
MODULE W3KOKKOSMD
  !/
  !/                  +-----------------------------------+
  !/                  | WAVEWATCH III           NOAA/NCEP |
  !/                  |          Kokkos port              |
  !/                  |                        FORTRAN 90 |
  !/                  | Last update :         15-Sep-2026 |
  !/                  +-----------------------------------+
  !/
  !/    15-Sep-2026 : Origination.                        ( version 7.14 )
  !/
  !  1. Purpose :
  !
  !     ISO_C_BINDING interface to the Kokkos port of the source terms, and the
  !     run-time switch that selects it.
  !
  !  2. Variables and types :
  !
  !      Name        Type  Scope    Description
  !     ----------------------------------------------------------------
  !      KOKKOS_SNL1 Log.  Public   Use the Kokkos DIA instead of W3SNL1.
  !     ----------------------------------------------------------------
  !
  !  3. Subroutines and functions :
  !
  !      Name              Type  Scope    Description
  !     ----------------------------------------------------------------
  !      W3KOKKOS_SETUP    Subr. Public   Read the switches from the shim.
  !      WW_KOKKOS_INIT    Func. Public   bind(C) ww_kokkos_init
  !      WW_KOKKOS_FINALIZE Subr. Public  bind(C) ww_kokkos_finalize
  !      WW_SNL1_INIT      Func. Public   bind(C) ww_snl1_init
  !      WW_SNL1           Subr. Public   bind(C) ww_snl1
  !      WW_SNL1_ENABLED   Func. Public   bind(C) ww_snl1_enabled
  !      WW_SNL1_LAST_ERROR Func. Public  bind(C) ww_snl1_last_error
  !     ----------------------------------------------------------------
  !
  !  4. Remarks :
  !
  !     The arrays passed to WW_SNL1 are assumed-size, so what reaches C is the
  !     address of the first element and the Fortran column-major order is kept.
  !     A(NSPEC,NPTS), CG(NK,NPTS), KDMEAN(NPTS), S(NSPEC,NPTS), D(NSPEC,NPTS).
  !
  !     Assumed-size dummies also mean the *lower bound of the actual argument is
  !     lost*. W3GDATMD allocates SIG(0:MK+1) (w3gdatmd.F90:2066), so passing the
  !     bare name SIG to WW_SNL1_INIT hands over SIG(0) and builds every
  !     quadruplet one frequency bin low. Pass SIG(1), which is WW3's own idiom.
  !     The same applies to any array whose lower bound is not 1; CG1 inside
  !     W3SRCE is declared CG1(NK) and is safe.
  !
  !     Three limits of phase 1, all from the shim holding one context:
  !
  !       a) One spectral grid per process. WW_SNL1_INIT replaces the previous
  !          grid's tables rather than adding a second grid, so a ww3_multi run
  !          must not call it once per grid -- with different NK*NTH the later
  !          grids would read out of bounds. The caller must refuse that case.
  !       b) One caller at a time. Nothing in the shim is locked, so concurrent
  !          WW_SNL1 calls are undefined behaviour. W3SRCE is called from inside
  !          !$OMP PARALLEL regions under W3_OMPG and W3_OMP0, so a phase-1
  !          caller either serialises the call (!$OMP CRITICAL) or runs without
  !          those switches. See PATCH.md.
  !       c) Every call copies host->device and back, so calling it once per sea
  !          point pays that cost once per sea point. Phase 2 removes it.
  !
  !/ ------------------------------------------------------------------- /
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  !
  IMPLICIT NONE
  !
  PUBLIC
  !
  !/ ------------------------------------------------------------------- /
  !/ Parameter list
  !/
  !> Use the Kokkos DIA in place of W3SNL1. Set by W3KOKKOS_SETUP from the
  !> environment variable WW_KOKKOS_SNL1 as read by WW_KOKKOS_INIT; .FALSE. until
  !> then, so a model that never calls the set-up keeps its Fortran source terms.
  LOGICAL :: KOKKOS_SNL1 = .FALSE.
  !> Kokkos W3SDS4, W3SIN4, W3SPR4: WW_KOKKOS_SDS4/SIN4/SPR4=1, read likewise.
  LOGICAL :: KOKKOS_SDS4 = .FALSE., KOKKOS_SIN4 = .FALSE.,             &
             KOKKOS_SPR4 = .FALSE.
  !/
  !/ ST4 packet indices, 0-based as in ww_kokkos/st4_params.hpp and
  !/ st4_arrays.hpp: Fortran RPAR(R_X+1). Change them only with those headers.
  !/
  INTEGER, PARAMETER :: R_DTH = 0, R_XFR = 1, R_FTE = 2, R_FTF = 3,    &
       R_AAIRCMIN = 4, R_AAIRGB = 5, R_IGPARS5 = 6,                    &
       R_WWNMEANPTAIL = 7, R_SSTXFTFTAIL = 8, R_AALPHA = 9,            &
       R_BBETA = 10, R_ZZ0MAX = 11, R_ZZ0RAT = 12, R_ZZALP = 13,       &
       R_SSINTHP = 14, R_TTAUWSHELTER = 15, R_SSWELLF = 16,            &
       R_SSDSC = R_SSWELLF + 7, R_SSDSBR = R_SSDSC + 21,               &
       R_SINTAILPAR = R_SSDSBR + 1, R_SSDSP = R_SINTAILPAR + 5,        &
       R_WWNMEANP = R_SSDSP + 1, R_SSTXFTF = R_SSDSP + 2,              &
       R_SSTXFTWN = R_SSDSP + 3, R_FFXPM = R_SSDSP + 4,                &
       R_FFXFM = R_SSDSP + 5, R_FFXFA = R_SSDSP + 6,                   &
       R_SSDSBRF1 = R_SSDSP + 7, R_SSDSBRF2 = R_SSDSP + 8,             &
       R_SSDSBINT = R_SSDSP + 9, R_SSDSBCK = R_SSDSP + 10,             &
       R_SSDSHCK = R_SSDSP + 11, R_SSDSABK = R_SSDSP + 12,             &
       R_SSDSPBK = R_SSDSP + 13, R_SSINBR = R_SSDSP + 14,              &
       R_CAPCHNK = R_SSDSP + 15, R_ZZWND = R_CAPCHNK + 10,             &
       R_SSDSCOS = R_ZZWND + 1, R_SSDSDTH = R_ZZWND + 2,               &
       R_SSDSBT = R_ZZWND + 3, R_SSDSBM = R_ZZWND + 4,                 &
       R_DELUST = R_SSDSBM + 5, R_DELALP = R_DELUST + 1,               &
       R_DELTAUW = R_DELUST + 2, R_DELU = R_DELUST + 3,                &
       R_DELTAIL = R_DELUST + 4, R_DELAB = R_DELUST + 5,               &
       NRPAR = R_DELUST + 6
  INTEGER, PARAMETER :: I_NK = 0, I_NTH = 1, I_NSPEC = 2, I_MK = 3,    &
       I_MTH = 4, I_MSPEC = 5, I_SSWELLFPAR = 6, I_SSDSISO = 7,        &
       I_SSDSBRFDF = 8, I_STAB3 = 9, I_IG1 = 10, I_FLX5 = 11,          &
       NIPAR = 12
  INTEGER, PARAMETER :: A_SIG = 0, A_SIG2 = 1, A_DSIP = 2, A_DDEN = 3, &
       A_DDEN2 = 4, A_ESIN = 5, A_ECOS = 6, A_TH = 7,                  &
       A_SATWEIGHTS = 8, A_CUMULW = 9, A_DCKI = 10, A_QBI = 11,        &
       A_SATINDICES = 12, A_IKTAB = 13, A_TAUT = 14, A_TAUHFT = 15,    &
       A_TAUHFT2 = 16, A_FWTABLE = 17, NARRAYS = 18
  !/
  !/ ------------------------------------------------------------------- /
  !/ The C ABI. One block per declaration of ww_kokkos_c.hpp.
  !/
  INTERFACE
     !
     !> Start the Kokkos runtime. COMM_F is an MPI_Comm_c2f handle or -1;
     !> phase 1 ignores it and takes the device from WW_KOKKOS_DEVICE_ID.
     !> Idempotent. Returns 0 on success.
     INTEGER(C_INT) FUNCTION WW_KOKKOS_INIT ( COMM_F )                  &
          BIND(C, NAME='ww_kokkos_init')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE :: COMM_F
     END FUNCTION WW_KOKKOS_INIT
     !
     !> Release the shim's buffers and shut Kokkos down if, and only if,
     !> WW_KOKKOS_INIT started it.
     SUBROUTINE WW_KOKKOS_FINALIZE ( ) BIND(C, NAME='ww_kokkos_finalize')
     END SUBROUTINE WW_KOKKOS_FINALIZE
     !
     !> Build the INSNL1 quadruplet tables for one spectral grid. SIG is
     !> SIG(1:NK) -- pass SIG(1), not SIG, when the actual argument's lower
     !> bound is 0 -- and is copied. One grid per process in phase 1: a second
     !> call replaces the first grid's tables. Returns 0 on success.
     INTEGER(C_INT) FUNCTION WW_SNL1_INIT                                &
          ( IGRID, NK, NTH, XFR, DTH, LAM, SNLC1, KDCON, KDMN,           &
            SNLS1, SNLS2, SNLS3, FACHFE, SIG )                           &
          BIND(C, NAME='ww_snl1_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NK, NTH
       REAL(C_FLOAT),  VALUE           :: XFR, DTH, LAM, SNLC1, KDCON,   &
                                          KDMN, SNLS1, SNLS2, SNLS3,     &
                                          FACHFE
       REAL(C_FLOAT),  INTENT(IN)      :: SIG(*)
     END FUNCTION WW_SNL1_INIT
     !
     !> W3SNL1 for NPTS sea points in one launch. S and D are overwritten, not
     !> accumulated. Status in WW_SNL1_LAST_ERROR().
     SUBROUTINE WW_SNL1 ( IGRID, NPTS, A, CG, KDMEAN, S, D )             &
          BIND(C, NAME='ww_snl1')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), KDMEAN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
     END SUBROUTINE WW_SNL1
     !
     !> 1 if WW_KOKKOS_SNL1=1 was set when WW_KOKKOS_INIT ran, else 0.
     INTEGER(C_INT) FUNCTION WW_SNL1_ENABLED ( )                         &
          BIND(C, NAME='ww_snl1_enabled')
       IMPORT :: C_INT
     END FUNCTION WW_SNL1_ENABLED
     !
     !> Status of the last WW_SNL1_INIT or WW_SNL1 call; 0 means success.
     INTEGER(C_INT) FUNCTION WW_SNL1_LAST_ERROR ( )                      &
          BIND(C, NAME='ww_snl1_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_SNL1_LAST_ERROR
     !
     !> ST4. Every per-point scalar is an array of NPTS; every OUT array is
     !> copied in as well (INTENT(INOUT)); LLWS crosses as INTEGER 0/1.
     INTEGER(C_INT) FUNCTION WW_ST4_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_st4_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_ST4_INIT
     !
     INTEGER(C_INT) FUNCTION WW_ST4_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_st4_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_ST4_SET_REAL
     !
     INTEGER(C_INT) FUNCTION WW_ST4_SET_INT ( IGRID, ID, X, RANK, LB,    &
          EXT ) BIND(C, NAME='ww_st4_set_int')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       INTEGER(C_INT), INTENT(IN)      :: X(*), LB(*), EXT(*)
     END FUNCTION WW_ST4_SET_INT
     !
     INTEGER(C_INT) FUNCTION WW_ST4_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_st4_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_ST4_COMMIT
     !
     SUBROUTINE WW_SDS4 ( IGRID, NPTS, A, K, CG, USTAR, USDIR, DEPTH,    &
          DAIR, DLWMEAN, FLOGRD57, FLOGRD58, SRHS, DDIAG, BRLAMBDA,      &
          WHITECAP ) BIND(C, NAME='ww_sds4')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS, FLOGRD57, FLOGRD58
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), K(*), CG(*), USTAR(*),   &
                                          USDIR(*), DEPTH(*), DAIR(*),   &
                                          DLWMEAN(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: SRHS(*), DDIAG(*),             &
                                          BRLAMBDA(*), WHITECAP(*)
     END SUBROUTINE WW_SDS4
     !
     SUBROUTINE WW_SIN4 ( IGRID, NPTS, A, CG, K, U, USTAR, DRAT, AS,     &
          USDIR, Z0, CD, TAUWX, TAUWY, TAUWNX, TAUWNY, S, D, LLWS )      &
          BIND(C, NAME='ww_sin4')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), K(*), U(*),       &
                                          USTAR(*), DRAT(*), AS(*),      &
                                          USDIR(*), Z0(*), CD(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: TAUWX(*), TAUWY(*), TAUWNX(*), &
                                          TAUWNY(*), S(*), D(*)
       INTEGER(C_INT), INTENT(INOUT)   :: LLWS(*)
     END SUBROUTINE WW_SIN4
     !
     SUBROUTINE WW_SPR4 ( IGRID, NPTS, A, CG, WN, U, UDIR, TAUWX, TAUWY, &
          LLWS, USTAR, USDIR, EMEAN, FMEAN, FMEAN1, WNMEAN, AMAX, CD,    &
          Z0, CHARN, FMEANWS, DLWMEAN ) BIND(C, NAME='ww_spr4')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*), U(*),      &
                                          UDIR(*), TAUWX(*), TAUWY(*)
       INTEGER(C_INT), INTENT(IN)      :: LLWS(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: USTAR(*), USDIR(*), EMEAN(*),  &
                                          FMEAN(*), FMEAN1(*), WNMEAN(*),&
                                          AMAX(*), CD(*), Z0(*),         &
                                          CHARN(*), FMEANWS(*),          &
                                          DLWMEAN(*)
     END SUBROUTINE WW_SPR4
     !
     !> WHICH = 1 SDS4, 2 SIN4, 3 SPR4, 4 ST4_TABLES.
     INTEGER(C_INT) FUNCTION WW_ST4_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_st4_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_ST4_ENABLED
     !
     INTEGER(C_INT) FUNCTION WW_ST4_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_st4_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_ST4_LAST_ERROR
     !
  END INTERFACE
  !/
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> @brief Read the run-time switches from the shim into module variables.
  !>
  !> Call once, after WW_KOKKOS_INIT and before the first source-term call. Kept
  !> separate from WW_KOKKOS_INIT so the model reads a LOGICAL in its inner loop
  !> rather than crossing the C boundary for every sea point.
  !>
  SUBROUTINE W3KOKKOS_SETUP
    IMPLICIT NONE
    !
    KOKKOS_SNL1 = ( WW_SNL1_ENABLED() == 1_C_INT )
    KOKKOS_SDS4 = ( WW_ST4_ENABLED(1_C_INT) == 1_C_INT )
    KOKKOS_SIN4 = ( WW_ST4_ENABLED(2_C_INT) == 1_C_INT )
    KOKKOS_SPR4 = ( WW_ST4_ENABLED(3_C_INT) == 1_C_INT )
    !
  END SUBROUTINE W3KOKKOS_SETUP
#ifdef W3_ST4
  !/ ------------------------------------------------------------------- /
  !> @brief Give the shim grid IMOD's ST4 scalars and arrays.
  !>
  !> Call from W3INIT after W3IOGR, with W3SETG/W3SETO done for IMOD, so the
  !> W3GDATMD aliases, the W3SRC4MD tables and FWTABLE hold this grid's values.
  !>
  SUBROUTINE W3KOKKOS_ST4_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, XFR, FTE, FTF, AAIRCMIN,  &
         AAIRGB, WWNMEANPTAIL, SSTXFTFTAIL, AALPHA, BBETA, ZZ0MAX,     &
         ZZ0RAT, ZZALP, SSINTHP, TTAUWSHELTER, SSWELLF, SSDSC, SSDSBR, &
         SINTAILPAR, SSDSP, WWNMEANP, SSTXFTF, SSTXFTWN, FFXPM, FFXFM, &
         FFXFA, SSDSBRF1, SSDSBRF2, SSDSBINT, SSDSBCK, SSDSHCK,        &
         SSDSABK, SSDSPBK, SSINBR, CAPCHNK, ZZWND, SSDSCOS, SSDSDTH,   &
         SSDSBT, SSDSBM, SSWELLFPAR, SSDSISO, SSDSBRFDF,               &
         SIG, SIG2, DSIP, DDEN, DDEN2, ESIN, ECOS, TH, SATWEIGHTS,     &
         CUMULW, DCKI, QBI, SATINDICES, IKTAB
#ifdef W3_IG1
    USE W3GDATMD, ONLY: IGPARS
#endif
    USE W3SRC4MD, ONLY: TAUT, TAUHFT, TAUHFT2, DELUST, DELALP,        &
         DELTAUW, DELU, DELTAIL
    USE CONSTANTS, ONLY: DELAB, FWTABLE
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(NRPAR)
    INTEGER(C_INT)          :: IPAR(NIPAR)
    INTEGER                 :: J
    !
    RPAR = 0.
    IPAR = 0
    RPAR(R_DTH+1)          = DTH
    RPAR(R_XFR+1)          = XFR
    RPAR(R_FTE+1)          = FTE
    RPAR(R_FTF+1)          = FTF
    RPAR(R_AAIRCMIN+1)     = AAIRCMIN
    RPAR(R_AAIRGB+1)       = AAIRGB
#ifdef W3_IG1
    RPAR(R_IGPARS5+1)      = IGPARS(5)
#endif
    RPAR(R_WWNMEANPTAIL+1) = WWNMEANPTAIL
    RPAR(R_SSTXFTFTAIL+1)  = SSTXFTFTAIL
    RPAR(R_AALPHA+1)       = AALPHA
    RPAR(R_BBETA+1)        = BBETA
    RPAR(R_ZZ0MAX+1)       = ZZ0MAX
    RPAR(R_ZZ0RAT+1)       = ZZ0RAT
    RPAR(R_ZZALP+1)        = ZZALP
    RPAR(R_SSINTHP+1)      = SSINTHP
    RPAR(R_TTAUWSHELTER+1) = TTAUWSHELTER
    DO J = 1, 7
      RPAR(R_SSWELLF+J)    = SSWELLF(J)
    END DO
    DO J = 1, 21
      RPAR(R_SSDSC+J)      = SSDSC(J)
    END DO
    RPAR(R_SSDSBR+1)       = SSDSBR
    DO J = 1, 5
      RPAR(R_SINTAILPAR+J) = SINTAILPAR(J)
    END DO
    RPAR(R_SSDSP+1)        = SSDSP
    RPAR(R_WWNMEANP+1)     = WWNMEANP
    RPAR(R_SSTXFTF+1)      = SSTXFTF
    RPAR(R_SSTXFTWN+1)     = SSTXFTWN
    RPAR(R_FFXPM+1)        = FFXPM
    RPAR(R_FFXFM+1)        = FFXFM
    RPAR(R_FFXFA+1)        = FFXFA
    RPAR(R_SSDSBRF1+1)     = SSDSBRF1
    RPAR(R_SSDSBRF2+1)     = SSDSBRF2
    RPAR(R_SSDSBINT+1)     = SSDSBINT
    RPAR(R_SSDSBCK+1)      = SSDSBCK
    RPAR(R_SSDSHCK+1)      = SSDSHCK
    RPAR(R_SSDSABK+1)      = SSDSABK
    RPAR(R_SSDSPBK+1)      = SSDSPBK
    RPAR(R_SSINBR+1)       = SSINBR
    DO J = 1, 10
      RPAR(R_CAPCHNK+J)    = CAPCHNK(J)
    END DO
    RPAR(R_ZZWND+1)        = ZZWND
    RPAR(R_SSDSCOS+1)      = SSDSCOS
    RPAR(R_SSDSDTH+1)      = SSDSDTH
    RPAR(R_SSDSBT+1)       = SSDSBT
    ! SSDSBM is SSDSBM(0:4) (w3gdatmd.F90:915); the pointer keeps that base.
    DO J = 0, 4
      RPAR(R_SSDSBM+1+J)   = SSDSBM(J)
    END DO
    RPAR(R_DELUST+1)       = DELUST
    RPAR(R_DELALP+1)       = DELALP
    RPAR(R_DELTAUW+1)      = DELTAUW
    RPAR(R_DELU+1)         = DELU
    RPAR(R_DELTAIL+1)      = DELTAIL
    RPAR(R_DELAB+1)        = DELAB
    !
    ! MK, MTH, MSPEC are W3DIMS arguments, not module variables: read them off
    ! the arrays W3DIMS allocated, SIG(0:MK+1), TH(MTH), SIG2(MSPEC).
    IPAR(I_NK+1)           = NK
    IPAR(I_NTH+1)          = NTH
    IPAR(I_NSPEC+1)        = NSPEC
    IPAR(I_MK+1)           = UBOUND(SIG,1) - 1
    IPAR(I_MTH+1)          = SIZE(TH)
    IPAR(I_MSPEC+1)        = SIZE(SIG2)
    IPAR(I_SSWELLFPAR+1)   = SSWELLFPAR
    IPAR(I_SSDSISO+1)      = SSDSISO
    IPAR(I_SSDSBRFDF+1)    = SSDSBRFDF
#ifdef W3_STAB3
    IPAR(I_STAB3+1)        = 1
#endif
#ifdef W3_IG1
    IPAR(I_IG1+1)          = 1
#endif
#ifdef W3_FLX5
    IPAR(I_FLX5+1)         = 1
#endif
    !
    CALL ST4_CHECK ( 'WW_ST4_INIT',                                     &
         WW_ST4_INIT ( IMOD, RPAR, NRPAR, IPAR, NIPAR ) )
    CALL ST4_SET_R1 ( IMOD, A_SIG, SIG, LBOUND(SIG) )
    CALL ST4_SET_R1 ( IMOD, A_SIG2, SIG2, LBOUND(SIG2) )
    CALL ST4_SET_R1 ( IMOD, A_DSIP, DSIP, LBOUND(DSIP) )
    CALL ST4_SET_R1 ( IMOD, A_DDEN, DDEN, LBOUND(DDEN) )
    CALL ST4_SET_R1 ( IMOD, A_DDEN2, DDEN2, LBOUND(DDEN2) )
    CALL ST4_SET_R1 ( IMOD, A_ESIN, ESIN, LBOUND(ESIN) )
    CALL ST4_SET_R1 ( IMOD, A_ECOS, ECOS, LBOUND(ECOS) )
    CALL ST4_SET_R1 ( IMOD, A_TH, TH, LBOUND(TH) )
    CALL ST4_SET_R2 ( IMOD, A_SATWEIGHTS, SATWEIGHTS, LBOUND(SATWEIGHTS) )
    CALL ST4_SET_R2 ( IMOD, A_CUMULW, CUMULW, LBOUND(CUMULW) )
    CALL ST4_SET_R2 ( IMOD, A_DCKI, DCKI, LBOUND(DCKI) )
    CALL ST4_SET_R2 ( IMOD, A_QBI, QBI, LBOUND(QBI) )
    CALL ST4_SET_I2 ( IMOD, A_SATINDICES, SATINDICES, LBOUND(SATINDICES) )
    CALL ST4_SET_I2 ( IMOD, A_IKTAB, IKTAB, LBOUND(IKTAB) )
    ! Module tables: INSIN4 allocates them only when SINTAILPAR(1) > 0.5
    ! (TAUHFT2 also needs TTAUWSHELTER > 0), w3src4md.F90:1118-1127.
    IF ( ALLOCATED(TAUT) )                                              &
         CALL ST4_SET_R2 ( IMOD, A_TAUT, TAUT, LBOUND(TAUT) )
    IF ( ALLOCATED(TAUHFT) )                                            &
         CALL ST4_SET_R2 ( IMOD, A_TAUHFT, TAUHFT, LBOUND(TAUHFT) )
    IF ( ALLOCATED(TAUHFT2) )                                           &
         CALL ST4_SET_R3 ( IMOD, A_TAUHFT2, TAUHFT2, LBOUND(TAUHFT2) )
    CALL ST4_SET_R1 ( IMOD, A_FWTABLE, FWTABLE, LBOUND(FWTABLE) )
    CALL ST4_CHECK ( 'WW_ST4_COMMIT', WW_ST4_COMMIT ( IMOD ) )
    !
  END SUBROUTINE W3KOKKOS_ST4_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the model on a non-zero status of the shim.
  SUBROUTINE ST4_CHECK ( WHAT, IERR )
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: WHAT
    INTEGER(C_INT), INTENT(IN)   :: IERR
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOSMD : '// &
           'KOKKOS '//WHAT//' FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
  END SUBROUTINE ST4_CHECK
  !/ ------------------------------------------------------------------- /
  !> Register one array. X is assumed-shape, so its own bounds are 1-based:
  !> the actual's LBOUND travels in LB, its SHAPE in SHAPE(X).
  SUBROUTINE ST4_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL ST4_CHECK ( 'WW_ST4_SET_REAL', WW_ST4_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST4_SET_R1
  SUBROUTINE ST4_SET_R2 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(2)
    REAL, INTENT(IN)        :: X(:,:)
    CALL ST4_CHECK ( 'WW_ST4_SET_REAL', WW_ST4_SET_REAL ( IMOD, ID, X, &
         2, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST4_SET_R2
  SUBROUTINE ST4_SET_R3 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(3)
    REAL, INTENT(IN)        :: X(:,:,:)
    CALL ST4_CHECK ( 'WW_ST4_SET_REAL', WW_ST4_SET_REAL ( IMOD, ID, X, &
         3, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST4_SET_R3
  SUBROUTINE ST4_SET_I2 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(2)
    INTEGER, INTENT(IN)     :: X(:,:)
    CALL ST4_CHECK ( 'WW_ST4_SET_INT', WW_ST4_SET_INT ( IMOD, ID, X,   &
         2, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST4_SET_I2
  !/ ------------------------------------------------------------------- /
  !> W3SDS4 on Kokkos for one point, the argument list of W3SDS4 less IX, IY.
  !> Serialised: the shim's contexts are unlocked and W3SRCE runs inside
  !> !$OMP PARALLEL under W3_OMPG / W3_OMP0.
  SUBROUTINE W3KOKKOS_SDS4 ( A, K, CG, USTAR, USDIR, DEPTH, DAIR, SRHS, &
       DDIAG, BRLAMBDA, WHITECAP, DLWMEAN )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    USE W3ODATMD, ONLY: FLOGRD
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), K(NK), CG(NK), USTAR, USDIR,  &
                               DEPTH, DAIR, DLWMEAN
    REAL, INTENT(INOUT)     :: SRHS(NSPEC), DDIAG(NSPEC),               &
                               BRLAMBDA(NSPEC), WHITECAP(4)
    REAL(C_FLOAT)           :: US(1), UD(1), DP(1), DA(1), DL(1)
    INTEGER(C_INT)          :: IERR
    US(1) = USTAR
    UD(1) = USDIR
    DP(1) = DEPTH
    DA(1) = DAIR
    DL(1) = DLWMEAN
    !$OMP CRITICAL (WW_KOKKOS_ST4)
    CALL WW_SDS4 ( IGRID, 1, A, K, CG, US, UD, DP, DA, DL,              &
         MERGE(1_C_INT, 0_C_INT, FLOGRD(5,7)),                          &
         MERGE(1_C_INT, 0_C_INT, FLOGRD(5,8)),                          &
         SRHS, DDIAG, BRLAMBDA, WHITECAP )
    IERR = WW_ST4_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST4)
    CALL ST4_CHECK ( 'W3SDS4', IERR )
  END SUBROUTINE W3KOKKOS_SDS4
  !/ ------------------------------------------------------------------- /
  !> W3SIN4 on Kokkos for one point, the argument list of W3SIN4 less IX, IY
  !> and BRLAMBDA, which W3SIN4 does not use.
  SUBROUTINE W3KOKKOS_SIN4 ( A, CG, K, U, USTAR, DRAT, AS, USDIR, Z0,  &
       CD, TAUWX, TAUWY, TAUWNX, TAUWNY, S, D, LLWS )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), K(NSPEC), U, USTAR,   &
                               DRAT, AS, USDIR, Z0, CD
    REAL, INTENT(INOUT)     :: TAUWX, TAUWY, TAUWNX, TAUWNY,            &
                               S(NSPEC), D(NSPEC)
    LOGICAL, INTENT(INOUT)  :: LLWS(NSPEC)
    REAL(C_FLOAT)           :: U1(1), US(1), DR(1), AS1(1), UD(1),      &
                               Z01(1), CD1(1), TX(1), TY(1), TNX(1),    &
                               TNY(1)
    INTEGER(C_INT)          :: LW(NSPEC), IERR
    U1(1) = U
    US(1) = USTAR
    DR(1) = DRAT
    AS1(1) = AS
    UD(1) = USDIR
    Z01(1) = Z0
    CD1(1) = CD
    TX(1) = TAUWX
    TY(1) = TAUWY
    TNX(1) = TAUWNX
    TNY(1) = TAUWNY
    LW = MERGE(1_C_INT, 0_C_INT, LLWS)
    !$OMP CRITICAL (WW_KOKKOS_ST4)
    CALL WW_SIN4 ( IGRID, 1, A, CG, K, U1, US, DR, AS1, UD, Z01, CD1,   &
         TX, TY, TNX, TNY, S, D, LW )
    IERR = WW_ST4_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST4)
    CALL ST4_CHECK ( 'W3SIN4', IERR )
    TAUWX = TX(1)
    TAUWY = TY(1)
    TAUWNX = TNX(1)
    TAUWNY = TNY(1)
    LLWS = ( LW /= 0_C_INT )
  END SUBROUTINE W3KOKKOS_SIN4
  !/ ------------------------------------------------------------------- /
  !> W3SPR4 on Kokkos for one point, the argument list of W3SPR4 without
  !> W3_FLX5 (refused by WW_ST4_INIT and at configure time).
  SUBROUTINE W3KOKKOS_SPR4 ( A, CG, WN, EMEAN, FMEAN, FMEAN1, WNMEAN,  &
       AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0, CHARN, LLWS,  &
       FMEANWS, DLWMEAN )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), WN(NK), U, UDIR,      &
                               TAUWX, TAUWY
    LOGICAL, INTENT(IN)     :: LLWS(NSPEC)
    REAL, INTENT(INOUT)     :: USTAR, USDIR, EMEAN, FMEAN, FMEAN1,      &
                               WNMEAN, AMAX, CD, Z0, CHARN, FMEANWS,    &
                               DLWMEAN
    REAL(C_FLOAT)           :: U1(1), UD1(1), TX(1), TY(1), US(1),      &
                               UD(1), EM(1), FM(1), FM1(1), WM(1),      &
                               AM(1), CD1(1), Z01(1), CH(1), FWS(1),    &
                               DL(1)
    INTEGER(C_INT)          :: LW(NSPEC), IERR
    U1(1) = U
    UD1(1) = UDIR
    TX(1) = TAUWX
    TY(1) = TAUWY
    US(1) = USTAR
    UD(1) = USDIR
    EM(1) = EMEAN
    FM(1) = FMEAN
    FM1(1) = FMEAN1
    WM(1) = WNMEAN
    AM(1) = AMAX
    CD1(1) = CD
    Z01(1) = Z0
    CH(1) = CHARN
    FWS(1) = FMEANWS
    DL(1) = DLWMEAN
    LW = MERGE(1_C_INT, 0_C_INT, LLWS)
    !$OMP CRITICAL (WW_KOKKOS_ST4)
    CALL WW_SPR4 ( IGRID, 1, A, CG, WN, U1, UD1, TX, TY, LW, US, UD,    &
         EM, FM, FM1, WM, AM, CD1, Z01, CH, FWS, DL )
    IERR = WW_ST4_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST4)
    CALL ST4_CHECK ( 'W3SPR4', IERR )
    USTAR = US(1)
    USDIR = UD(1)
    EMEAN = EM(1)
    FMEAN = FM(1)
    FMEAN1 = FM1(1)
    WNMEAN = WM(1)
    AMAX = AM(1)
    CD = CD1(1)
    Z0 = Z01(1)
    CHARN = CH(1)
    FMEANWS = FWS(1)
    DLWMEAN = DL(1)
  END SUBROUTINE W3KOKKOS_SPR4
#endif
  !/ ------------------------------------------------------------------- /
  !/
END MODULE W3KOKKOSMD
