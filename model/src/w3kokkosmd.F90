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
  !> Kokkos W3SRCE (with every source term it calls): WW_KOKKOS_SRCE=1.
  LOGICAL :: KOKKOS_SRCE = .FALSE.
#ifdef W3_ST0
  !> Kokkos W3SPR0: WW_KOKKOS_SPR0=1, read by W3KOKKOS_SETUP.
  LOGICAL :: KOKKOS_SPR0 = .FALSE.
  !/
  !/ ST0 packet indices and array ids, 0-based as WwSt0Rpar, WwSt0Ipar and
  !/ WwSt0Array of ww_kokkos_c.hpp.
  !/
  INTEGER, PARAMETER :: R0_FTE = 0, R0_FTF = 1, R0_FTWN = 2, NR0PAR = 3
  INTEGER, PARAMETER :: I0_NK = 0, I0_NTH = 1, I0_NSPEC = 2, NI0PAR = 3
  INTEGER, PARAMETER :: A0_SIG = 0, A0_DDEN = 1
#endif
#if defined(W3_FLX1) || defined(W3_FLX2) || defined(W3_FLX3) || defined(W3_FLX4) || defined(W3_FLX5)
  !> Kokkos W3FLX1..W3FLX5: WW_KOKKOS_FLX1..FLX5=1. Only the flag of the
  !> compiled W3_FLXn has an effect.
  LOGICAL :: KOKKOS_FLX1 = .FALSE., KOKKOS_FLX2 = .FALSE.,             &
             KOKKOS_FLX3 = .FALSE., KOKKOS_FLX4 = .FALSE.,             &
             KOKKOS_FLX5 = .FALSE.
#endif
#ifdef W3_ST6
  !> Kokkos W3SWL6: WW_KOKKOS_SWL6=1, read by W3KOKKOS_SETUP.
  LOGICAL :: KOKKOS_SWL6 = .FALSE.
  !/
  !/ SWL packet indices and array ids, 0-based as WwSwlRpar, WwSwlIpar and
  !/ WwSwlArray of ww_kokkos_c.hpp.
  !/
  INTEGER, PARAMETER :: SWL_R_DTH = 0, SWL_R_SWL6B1 = 1, SWL_NRPAR = 2
  INTEGER, PARAMETER :: SWL_I_NK = 0, SWL_I_NTH = 1, SWL_I_NSPEC = 2,  &
       SWL_I_SWL6CSTB1 = 3, SWL_NIPAR = 4
  INTEGER, PARAMETER :: SWL_A_SIG = 0, SWL_A_DDEN = 1
#endif
  !/
  !/ W3SRCE packet indices, 0-based as in ww_kokkos/srce_params.hpp.
  !/
  INTEGER, PARAMETER :: RS_DMIN = 0, RS_DTMAX = 1, RS_DTMIN = 2,       &
       RS_FACTI1 = 3, RS_FACTI2 = 4, RS_FACHFA = 5, RS_FACP = 6,        &
       RS_XFLT = 7, RS_XREL = 8, RS_FHMAX = 9, RS_FFXFM = 10,           &
       RS_FFXPM = 11, RS_FFXFA = 12, RS_ICESCALES = 13,                 &
       RS_SLNC1 = RS_ICESCALES + 4, RS_FSPM = RS_SLNC1 + 1,             &
       RS_FSHF = RS_SLNC1 + 2, RS_SBTC1 = RS_SLNC1 + 3,                 &
       RS_SDBC1 = RS_SLNC1 + 4, RS_SDBC2 = RS_SLNC1 + 5,                &
       NRSPAR = RS_SLNC1 + 6
  INTEGER, PARAMETER :: IS_NK = 0, IS_NTH = 1, IS_NSPEC = 2,           &
       IS_IQTPE = 3, IS_IICEDISP = 4, IS_IC_NUMERICS = 5, IS_FDONLY = 6, &
       IS_LAB_SWITCH = 7, NISPAR = 8
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
       R_WWNMEANP = R_SSDSP + 1, R_SSTXFTWN = R_SSDSP + 2,             &
       R_SSDSBRF1 = R_SSDSP + 3, R_SSDSBINT = R_SSDSP + 4,             &
       R_SSDSBCK = R_SSDSP + 5, R_SSDSHCK = R_SSDSP + 6,               &
       R_SSDSABK = R_SSDSP + 7, R_SSDSPBK = R_SSDSP + 8,               &
       R_CAPCHNK = R_SSDSP + 9, R_ZZWND = R_CAPCHNK + 10,              &
       R_SSDSCOS = R_ZZWND + 1, R_SSDSDTH = R_ZZWND + 2,               &
       R_SSDSBT = R_ZZWND + 3, R_SSDSBM = R_ZZWND + 4,                 &
       R_DELUST = R_SSDSBM + 5, R_DELALP = R_DELUST + 1,               &
       R_DELTAUW = R_DELUST + 2, R_DELU = R_DELUST + 3,                &
       R_DELTAIL = R_DELUST + 4, R_DELAB = R_DELUST + 5,               &
       NRPAR = R_DELUST + 6
  INTEGER, PARAMETER :: I_NK = 0, I_NTH = 1, I_NSPEC = 2, I_MK = 3,    &
       I_MTH = 4, I_MSPEC = 5, I_STAB3 = 6, I_IG1 = 7, I_FLX5 = 8,     &
       NIPAR = 9
  INTEGER, PARAMETER :: A_SIG = 0, A_SIG2 = 1, A_DSIP = 2, A_DDEN = 3, &
       A_DDEN2 = 4, A_ESIN = 5, A_ECOS = 6, A_TH = 7,                  &
       A_SATWEIGHTS = 8, A_CUMULW = 9, A_DCKI = 10, A_QBI = 11,        &
       A_SATINDICES = 12, A_IKTAB = 13, A_TAUT = 14, A_TAUHFT = 15,    &
       A_TAUHFT2 = 16, A_FWTABLE = 17, NARRAYS = 18
  !/
#if defined(W3_FLX1) || defined(W3_FLX2) || defined(W3_FLX3) || defined(W3_FLX4) || defined(W3_FLX5)
  !/ FLX packet indices, 0-based as WwFlxRpar, WwFlxIpar, WwFlxArray in
  !/ ww_kokkos_c.hpp ("FLX"): Fortran RPAR(FR_X+1).
  !/
  INTEGER, PARAMETER :: FR_CINXSI = 0, FR_CD_MAX = 1, FR_FLX4A0 = 2,   &
       FR_DSIE = 3, NFRPAR = 4
  INTEGER, PARAMETER :: FI_FLX = 0, FI_NITTIN = 1, FI_CAP_ID = 2,      &
       FI_N1MAX = 3, NFIPAR = 4
  INTEGER, PARAMETER :: FA_EWN1 = 0
  !/
#endif
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
     !> W3SRCE. Needs WW_SNL1_INIT and WW_ST4_COMMIT done for IGRID. SIG,
     !> DDEN are (1:NK), ECOS, ESIN, MAPWN (1:NSPEC); copied.
     INTEGER(C_INT) FUNCTION WW_SRCE_INIT ( IGRID, RPAR, NR, IPAR, NI,  &
          SIG, DDEN, ECOS, ESIN, MAPWN ) BIND(C, NAME='ww_srce_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*), SIG(*), DDEN(*),      &
                                          ECOS(*), ESIN(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*), MAPWN(*)
     END FUNCTION WW_SRCE_INIT
     !
     !> W3SRCE(srce_direct, ...) for one point, the arguments the lab switch
     !> reads. SHAVE crosses as INTEGER 0/1. Status in WW_SRCE_LAST_ERROR().
     SUBROUTINE WW_SRCE ( IGRID, IT, INFLAGS24, FLOGRD57, FLOGRD58,      &
          SPEC, WN1, CG1, D_INP, U10ABS, U10DIR, AS, ICE, DTG, DAIR,     &
          COEF, USTAR, USTDIR, TAUWX, TAUWY, DTDYN, FCUT, TAUOX, TAUOY,  &
          TAUWIX, TAUWIY, TAUWNX, TAUWNY, PHIAW, CHARN, TWS, PHIOC,      &
          WHITECAP, PHIBBL, TAUBBL, TAUICE, PHICE, TAUOCX, TAUOCY,       &
          WNMEAN, SHAVE ) BIND(C, NAME='ww_srce')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, IT, INFLAGS24,          &
                                          FLOGRD57, FLOGRD58
       REAL(C_FLOAT),  VALUE           :: D_INP, U10ABS, U10DIR, AS,     &
                                          ICE, DTG, DAIR, COEF
       REAL(C_FLOAT),  INTENT(IN)      :: WN1(*), CG1(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: SPEC(*), USTAR, USTDIR, TAUWX, &
                                          TAUWY, DTDYN, FCUT, TAUOX,     &
                                          TAUOY, TAUWIX, TAUWIY, TAUWNX, &
                                          TAUWNY, PHIAW, CHARN, TWS,     &
                                          PHIOC, WHITECAP(*), PHIBBL,    &
                                          TAUBBL(*), TAUICE(*), PHICE,   &
                                          TAUOCX, TAUOCY, WNMEAN
       INTEGER(C_INT), INTENT(INOUT)   :: SHAVE
     END SUBROUTINE WW_SRCE
     !
     INTEGER(C_INT) FUNCTION WW_SRCE_ENABLED ( )                         &
          BIND(C, NAME='ww_srce_enabled')
       IMPORT :: C_INT
     END FUNCTION WW_SRCE_ENABLED
     !
     INTEGER(C_INT) FUNCTION WW_SRCE_LAST_ERROR ( )                      &
          BIND(C, NAME='ww_srce_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_SRCE_LAST_ERROR
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
#ifdef W3_ST0
     !> ST0 configuration of grid IGRID (W3INIT's IMOD): scalars, then each
     !> array with its LBOUND and SHAPE, then the commit.
     INTEGER(C_INT) FUNCTION WW_ST0_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_st0_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_ST0_INIT
     !
     INTEGER(C_INT) FUNCTION WW_ST0_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_st0_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_ST0_SET_REAL
     !
     INTEGER(C_INT) FUNCTION WW_ST0_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_st0_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_ST0_COMMIT
     !
     !> W3SPR0. A (NSPEC,NPTS); CG, WN (NK,NPTS); the four results (NPTS).
     SUBROUTINE WW_SPR0 ( IGRID, NPTS, A, CG, WN, EMEAN, FMEAN, WNMEAN,  &
          AMAX ) BIND(C, NAME='ww_spr0')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: EMEAN(*), FMEAN(*), WNMEAN(*), &
                                          AMAX(*)
     END SUBROUTINE WW_SPR0
     !
     !> WHICH = 1 SPR0.
     INTEGER(C_INT) FUNCTION WW_ST0_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_st0_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_ST0_ENABLED
     !
     INTEGER(C_INT) FUNCTION WW_ST0_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_st0_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_ST0_LAST_ERROR
#endif
     !
#if defined(W3_FLX1) || defined(W3_FLX2) || defined(W3_FLX3) || defined(W3_FLX4) || defined(W3_FLX5)
     !> FLX. Every per-point scalar is an array of NPTS; every OUT argument is
     !> copied in as well (INTENT(INOUT)). IEXT is the code of the EXTCDE the
     !> Fortran would call (ZWND not 10 m), 0 otherwise.
     INTEGER(C_INT) FUNCTION WW_FLX_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_flx_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_FLX_INIT
     !
     INTEGER(C_INT) FUNCTION WW_FLX_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_flx_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_FLX_SET_REAL
     !
     INTEGER(C_INT) FUNCTION WW_FLX_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_flx_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_FLX_COMMIT
     !
     SUBROUTINE WW_FLX1 ( IGRID, NPTS, ZWND, U10, U10D, UST, USTD, Z0,   &
          CD, IEXT ) BIND(C, NAME='ww_flx1')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: ZWND(*), U10(*), U10D(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: UST(*), USTD(*), Z0(*), CD(*)
       INTEGER(C_INT), INTENT(INOUT)   :: IEXT(*)
     END SUBROUTINE WW_FLX1
     !
     SUBROUTINE WW_FLX2 ( IGRID, NPTS, ZWIND, DEPTH, FP, U, UDIR, UST,   &
          USTD, Z0, CD ) BIND(C, NAME='ww_flx2')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: ZWIND(*), DEPTH(*), FP(*),     &
                                          U(*), UDIR(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: UST(*), USTD(*), Z0(*), CD(*)
     END SUBROUTINE WW_FLX2
     !
     SUBROUTINE WW_FLX3 ( IGRID, NPTS, ZWIND, DEPTH, FP, U, UDIR, UST,   &
          USTD, Z0, CD ) BIND(C, NAME='ww_flx3')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: ZWIND(*), DEPTH(*), FP(*),     &
                                          U(*), UDIR(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: UST(*), USTD(*), Z0(*), CD(*)
     END SUBROUTINE WW_FLX3
     !
     SUBROUTINE WW_FLX4 ( IGRID, NPTS, ZWND, U10, U10D, UST, USTD, Z0,   &
          CD, IEXT ) BIND(C, NAME='ww_flx4')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: ZWND(*), U10(*), U10D(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: UST(*), USTD(*), Z0(*), CD(*)
       INTEGER(C_INT), INTENT(INOUT)   :: IEXT(*)
     END SUBROUTINE WW_FLX4
     !
     SUBROUTINE WW_FLX5 ( IGRID, NPTS, ZWND, U10, U10D, TAUA, TAUADIR,   &
          RHOAIR, UST, USTD, Z0, CD, CHARN, IEXT )                       &
          BIND(C, NAME='ww_flx5')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: ZWND(*), U10(*), U10D(*),      &
                                          TAUA(*), TAUADIR(*), RHOAIR(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: UST(*), USTD(*), Z0(*), CD(*), &
                                          CHARN(*)
       INTEGER(C_INT), INTENT(INOUT)   :: IEXT(*)
     END SUBROUTINE WW_FLX5
     !
     !> WHICH = 1..5 FLX1..FLX5, 6 FLX_TABLES.
     INTEGER(C_INT) FUNCTION WW_FLX_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_flx_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_FLX_ENABLED
     !
     INTEGER(C_INT) FUNCTION WW_FLX_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_flx_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_FLX_LAST_ERROR
#endif
#ifdef W3_ST6
     !> SWL configuration of grid IGRID (W3INIT's IMOD): scalars, then each
     !> array with its LBOUND and SHAPE, then the commit.
     INTEGER(C_INT) FUNCTION WW_SWL_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_swl_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_SWL_INIT
     !
     INTEGER(C_INT) FUNCTION WW_SWL_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_swl_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_SWL_SET_REAL
     !
     INTEGER(C_INT) FUNCTION WW_SWL_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_swl_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_SWL_COMMIT
     !
     !> W3SWL6. A, S, D (NSPEC,NPTS); CG, WN (NK,NPTS). Returns 0, or 2
     !> where W3SWL6 stops in EXTCDE(2), with NIKN = SIZE(IKN).
     INTEGER(C_INT) FUNCTION WW_SWL6 ( IGRID, NPTS, A, CG, WN, S, D,     &
          NIKN ) BIND(C, NAME='ww_swl6')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
       INTEGER(C_INT), INTENT(OUT)     :: NIKN
     END FUNCTION WW_SWL6
     !
     !> WHICH = 1 SWL6.
     INTEGER(C_INT) FUNCTION WW_SWL_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_swl_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_SWL_ENABLED
     !
     INTEGER(C_INT) FUNCTION WW_SWL_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_swl_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_SWL_LAST_ERROR
#endif
  END INTERFACE
#ifdef W3_ST1
  !/
  !> Kokkos W3SPR1, W3SIN1, W3SDS1: WW_KOKKOS_SPR1/SIN1/SDS1=1. Set by
  !> W3KOKKOS_ST1_SETUP, .FALSE. until then.
  LOGICAL :: KOKKOS_SPR1 = .FALSE., KOKKOS_SIN1 = .FALSE.,             &
             KOKKOS_SDS1 = .FALSE.
  !/
  !/ ST1 packet indices and array ids, 0-based as in ww_kokkos_c.hpp
  !/ ("ST1"): Fortran RPAR(ST1_R_X+1).
  !/
  INTEGER, PARAMETER :: ST1_R_FTE = 0, ST1_R_FTF = 1, ST1_R_FTWN = 2,  &
       ST1_R_SINC1 = 3, ST1_R_SDSC1 = 4, ST1_NRPAR = 5
  INTEGER, PARAMETER :: ST1_I_NK = 0, ST1_I_NTH = 1, ST1_I_NSPEC = 2,  &
       ST1_NIPAR = 3
  INTEGER, PARAMETER :: ST1_A_SIG = 0, ST1_A_DDEN = 1, ST1_A_SIG2 = 2, &
       ST1_A_ECOS = 3, ST1_A_ESIN = 4
  !/
  INTERFACE
     !
     INTEGER(C_INT) FUNCTION WW_ST1_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_st1_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_ST1_INIT
     !
     INTEGER(C_INT) FUNCTION WW_ST1_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_st1_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_ST1_SET_REAL
     !
     INTEGER(C_INT) FUNCTION WW_ST1_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_st1_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_ST1_COMMIT
     !
     SUBROUTINE WW_SPR1 ( IGRID, NPTS, A, CG, WN, EMEAN, FMEAN, WNMEAN,  &
          AMAX ) BIND(C, NAME='ww_spr1')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: EMEAN(*), FMEAN(*), WNMEAN(*), &
                                          AMAX(*)
     END SUBROUTINE WW_SPR1
     !
     SUBROUTINE WW_SIN1 ( IGRID, NPTS, A, K, USTAR, USDIR, S, D )        &
          BIND(C, NAME='ww_sin1')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), K(*), USTAR(*), USDIR(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
     END SUBROUTINE WW_SIN1
     !
     SUBROUTINE WW_SDS1 ( IGRID, NPTS, A, K, EMEAN, FMEAN, WNMEAN, S, D ) &
          BIND(C, NAME='ww_sds1')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), K(*), EMEAN(*), FMEAN(*),&
                                          WNMEAN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
     END SUBROUTINE WW_SDS1
     !
     !> WHICH = 1 SPR1, 2 SIN1, 3 SDS1.
     INTEGER(C_INT) FUNCTION WW_ST1_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_st1_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_ST1_ENABLED
     !
     INTEGER(C_INT) FUNCTION WW_ST1_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_st1_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_ST1_LAST_ERROR
     !
  END INTERFACE
#endif
#ifdef W3_ST2
  !/
  !> Kokkos W3SPR2, W3SIN2, W3SDS2: WW_KOKKOS_SPR2/SIN2/SDS2=1, read by
  !> W3KOKKOS_ST2_SETUP. .FALSE. until then: the Fortran runs.
  LOGICAL :: KOKKOS_SPR2 = .FALSE., KOKKOS_SIN2 = .FALSE.,             &
             KOKKOS_SDS2 = .FALSE.
  !/
  !/ ST2 packet indices, 0-based as WwSt2Rpar, WwSt2Ipar and WwSt2Array in
  !/ ww_kokkos_c.hpp: Fortran RPAR(ST2_R_X+1). Change them only with it.
  !/
  INTEGER, PARAMETER :: ST2_R_DTH = 0, ST2_R_XFR = 1, ST2_R_FTE = 2,   &
       ST2_R_FTF = 3, ST2_R_FTWN = 4, ST2_R_FTTR = 5,                   &
       ST2_R_FPIMIN = 6, ST2_R_ZWIND = 7, ST2_R_FSWELL = 8,             &
       ST2_R_FACTI1 = 9, ST2_R_FACTI2 = 10, ST2_R_XF1 = 11,             &
       ST2_R_XF2 = 12, ST2_R_XFH = 13, ST2_R_SDSALN = 14,               &
       ST2_R_CDSA0 = 15, ST2_R_CDSA1 = 16, ST2_R_CDSA2 = 17,            &
       ST2_R_CDSB0 = 18, ST2_R_CDSB1 = 19, ST2_R_CDSB2 = 20,            &
       ST2_R_CDSB3 = 21, ST2_NRPAR = 22
  INTEGER, PARAMETER :: ST2_I_NK = 0, ST2_I_NTH = 1, ST2_I_NSPEC = 2,  &
       ST2_NIPAR = 3
  INTEGER, PARAMETER :: ST2_A_SIG = 0, ST2_A_SIG2 = 1, ST2_A_DDEN = 2, &
       ST2_A_ECOS = 3, ST2_A_ESIN = 4, ST2_NARRAYS = 5
  !/
  INTERFACE
     INTEGER(C_INT) FUNCTION WW_ST2_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_st2_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_ST2_INIT
     INTEGER(C_INT) FUNCTION WW_ST2_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_st2_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_ST2_SET_REAL
     INTEGER(C_INT) FUNCTION WW_ST2_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_st2_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_ST2_COMMIT
     SUBROUTINE WW_SPR2 ( IGRID, NPTS, A, CG, WN, DEPTH, FPI, U, USTAR,  &
          EMEAN, FMEAN, WNMEAN, AMAX, ALFA, FP ) BIND(C, NAME='ww_spr2')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*), DEPTH(*),  &
                                          FPI(*), U(*), USTAR(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: EMEAN(*), FMEAN(*), WNMEAN(*), &
                                          AMAX(*), ALFA(*), FP(*)
     END SUBROUTINE WW_SPR2
     SUBROUTINE WW_SIN2 ( IGRID, NPTS, A, CG, K, U, UDIR, CD, Z0, FPI,   &
          S, D ) BIND(C, NAME='ww_sin2')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), K(*), U(*),       &
                                          UDIR(*), CD(*), Z0(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: FPI(*), S(*), D(*)
     END SUBROUTINE WW_SIN2
     SUBROUTINE WW_SDS2 ( IGRID, NPTS, A, CG, K, FPI, USTAR, ALFA, S, D ) &
          BIND(C, NAME='ww_sds2')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), K(*), FPI(*),     &
                                          USTAR(*), ALFA(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
     END SUBROUTINE WW_SDS2
     INTEGER(C_INT) FUNCTION WW_ST2_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_st2_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_ST2_ENABLED
     INTEGER(C_INT) FUNCTION WW_ST2_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_st2_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_ST2_LAST_ERROR
  END INTERFACE
#endif
#ifdef W3_ST3
  !/
  !> Kokkos W3SPR3, W3SIN3, W3SDS3: WW_KOKKOS_SPR3/SIN3/SDS3=1, read by
  !> W3KOKKOS_ST3_SETUP. .FALSE. until then: the Fortran runs.
  LOGICAL :: KOKKOS_SPR3 = .FALSE., KOKKOS_SIN3 = .FALSE.,             &
             KOKKOS_SDS3 = .FALSE.
  !/
  !/ ST3 packet indices, 0-based as WwSt3Rpar, WwSt3Ipar and WwSt3Array in
  !/ ww_kokkos_c.hpp: Fortran RPAR(ST3_R_X+1). Change them only with it.
  !/
  INTEGER, PARAMETER :: ST3_R_DTH = 0, ST3_R_ZZWND = 1,                 &
       ST3_R_AALPHA = 2, ST3_R_BBETA = 3, ST3_R_ZZ0MAX = 4,             &
       ST3_R_ZZALP = 5, ST3_R_SSINTHP = 6, ST3_R_SSWELLF1 = 7,          &
       ST3_R_SSDSC1 = 8, ST3_R_DDELTA1 = 9, ST3_R_DDELTA2 = 10,         &
       ST3_R_WWNMEANP = 11, ST3_R_WWNMEANPTAIL = 12, ST3_R_FTE = 13,    &
       ST3_R_SSTXFTF = 14, ST3_R_SSTXFTWN = 15,                         &
       ST3_R_SSTXFTFTAIL = 16, ST3_NRPAR = 17
  INTEGER, PARAMETER :: ST3_I_NK = 0, ST3_I_NTH = 1, ST3_I_NSPEC = 2,  &
       ST3_I_STAB3 = 3, ST3_NIPAR = 4
  INTEGER, PARAMETER :: ST3_A_SIG = 0, ST3_A_SIG2 = 1, ST3_A_DDEN = 2, &
       ST3_A_DDEN2 = 3, ST3_A_ECOS = 4, ST3_A_ESIN = 5, ST3_NARRAYS = 6
  !/
  INTERFACE
     INTEGER(C_INT) FUNCTION WW_ST3_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_st3_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_ST3_INIT
     INTEGER(C_INT) FUNCTION WW_ST3_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_st3_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_ST3_SET_REAL
     INTEGER(C_INT) FUNCTION WW_ST3_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_st3_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_ST3_COMMIT
     SUBROUTINE WW_SPR3 ( IGRID, NPTS, A, CG, WN, EMEAN, FMEAN, FMEANS,  &
          WNMEAN, AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0,     &
          CHARN, LLWS, FMEANWS ) BIND(C, NAME='ww_spr3')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*), U(*),      &
                                          UDIR(*), TAUWX(*), TAUWY(*)
       INTEGER(C_INT), INTENT(IN)      :: LLWS(*)
       REAL(C_FLOAT),  INTENT(INOUT)   :: USTAR(*), USDIR(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: EMEAN(*), FMEAN(*), FMEANS(*), &
                                          WNMEAN(*), AMAX(*), CD(*),     &
                                          Z0(*), CHARN(*), FMEANWS(*)
     END SUBROUTINE WW_SPR3
     SUBROUTINE WW_SIN3 ( IGRID, NPTS, A, CG, K, U, USTAR, DRAT, AS,     &
          USDIR, Z0, TAUWX, TAUWY, TAUWNX, TAUWNY, S, D, LLWS )          &
          BIND(C, NAME='ww_sin3')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), K(*), U(*),       &
                                          USTAR(*), DRAT(*), AS(*),      &
                                          USDIR(*), Z0(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: TAUWX(*), TAUWY(*), TAUWNX(*), &
                                          TAUWNY(*), S(*), D(*)
       INTEGER(C_INT), INTENT(OUT)     :: LLWS(*)
     END SUBROUTINE WW_SIN3
     SUBROUTINE WW_SDS3 ( IGRID, NPTS, A, K, EMEAN, FMEAN, WNMEAN, S, D ) &
          BIND(C, NAME='ww_sds3')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), K(*), EMEAN(*), FMEAN(*), &
                                          WNMEAN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
     END SUBROUTINE WW_SDS3
     INTEGER(C_INT) FUNCTION WW_ST3_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_st3_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_ST3_ENABLED
     INTEGER(C_INT) FUNCTION WW_ST3_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_st3_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_ST3_LAST_ERROR
  END INTERFACE
#endif
#ifdef W3_ST6
  !/
  !> Kokkos W3SPR6, W3SIN6, W3SDS6: WW_KOKKOS_SPR6/SIN6/SDS6=1, read by
  !> W3KOKKOS_ST6_SETUP. .FALSE. until then: the Fortran runs.
  LOGICAL :: KOKKOS_SPR6 = .FALSE., KOKKOS_SIN6 = .FALSE.,             &
             KOKKOS_SDS6 = .FALSE.
  !/
  !/ ST6 packet indices, 0-based as WwSt6Rpar, WwSt6Ipar, WwSt6Array and
  !/ WwSt6Code in ww_kokkos_c.hpp: Fortran RPAR(ST6_R_X+1). Change them only
  !/ with it.
  !/
  INTEGER, PARAMETER :: ST6_R_DTH = 0, ST6_R_XFR = 1, ST6_R_FTE = 2,   &
       ST6_R_FTF = 3, ST6_R_FTWN = 4, ST6_R_SIN6A0 = 5,                 &
       ST6_R_SIN6WS = 6, ST6_R_SIN6CHKMIN = 7, ST6_R_SIN6CHKCAP = 8,    &
       ST6_R_SIN6CHKINF = 9, ST6_R_SIN6CHKSIG = 10, ST6_R_SDS6A1 = 11,  &
       ST6_R_SDS6A2 = 12, ST6_NRPAR = 13
  INTEGER, PARAMETER :: ST6_I_NK = 0, ST6_I_NTH = 1, ST6_I_NSPEC = 2,  &
       ST6_I_SIN6FLCAP = 3, ST6_I_SDS6P1 = 4, ST6_I_SDS6P2 = 5,         &
       ST6_I_SDS6ET = 6, ST6_NIPAR = 7
  INTEGER, PARAMETER :: ST6_A_SIG = 0, ST6_A_DDEN = 1, ST6_A_DSII = 2, &
       ST6_A_SIG2 = 3, ST6_A_DDEN2 = 4, ST6_A_ECOS = 5, ST6_A_ESIN = 6
  INTEGER, PARAMETER :: ST6_C_IKN = 1, ST6_C_LF_IK10HZ = 2,            &
       ST6_C_LF_MK10HZ = 3, ST6_C_TWA_MK10HZ = 4, ST6_C_LF_ITERMAX = 5
  !/
  INTERFACE
     INTEGER(C_INT) FUNCTION WW_ST6_INIT ( IGRID, RPAR, NR, IPAR, NI )   &
          BIND(C, NAME='ww_st6_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NR, NI
       REAL(C_FLOAT),  INTENT(IN)      :: RPAR(*)
       INTEGER(C_INT), INTENT(IN)      :: IPAR(*)
     END FUNCTION WW_ST6_INIT
     INTEGER(C_INT) FUNCTION WW_ST6_SET_REAL ( IGRID, ID, X, RANK, LB,   &
          EXT ) BIND(C, NAME='ww_st6_set_real')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, ID, RANK
       REAL(C_FLOAT),  INTENT(IN)      :: X(*)
       INTEGER(C_INT), INTENT(IN)      :: LB(*), EXT(*)
     END FUNCTION WW_ST6_SET_REAL
     INTEGER(C_INT) FUNCTION WW_ST6_COMMIT ( IGRID )                     &
          BIND(C, NAME='ww_st6_commit')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: IGRID
     END FUNCTION WW_ST6_COMMIT
     SUBROUTINE WW_SPR6 ( IGRID, NPTS, A, CG, WN, EMEAN, FMEAN, WNMEAN,  &
          AMAX, FP ) BIND(C, NAME='ww_spr6')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: EMEAN(*), FMEAN(*), WNMEAN(*), &
                                          AMAX(*), FP(*)
     END SUBROUTINE WW_SPR6
     INTEGER(C_INT) FUNCTION WW_SIN6 ( IGRID, NPTS, A, CG, WN2, UABS,    &
          USTAR, USDIR, CD, DAIR, TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN,   &
          S, D, INFO ) BIND(C, NAME='ww_sin6')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN2(*), UABS(*),  &
                                          USTAR(*), USDIR(*), CD(*),     &
                                          DAIR(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: TAUWX(*), TAUWY(*), TAUNWX(*), &
                                          TAUNWY(*), CHARN(*), S(*), D(*)
       INTEGER(C_INT), INTENT(OUT)     :: INFO(2)
     END FUNCTION WW_SIN6
     INTEGER(C_INT) FUNCTION WW_SDS6 ( IGRID, NPTS, A, CG, WN, S, D,     &
          INFO ) BIND(C, NAME='ww_sds6')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: IGRID, NPTS
       REAL(C_FLOAT),  INTENT(IN)      :: A(*), CG(*), WN(*)
       REAL(C_FLOAT),  INTENT(OUT)     :: S(*), D(*)
       INTEGER(C_INT), INTENT(OUT)     :: INFO(2)
     END FUNCTION WW_SDS6
     !> WHICH = 1 SPR6, 2 SIN6, 3 SDS6.
     INTEGER(C_INT) FUNCTION WW_ST6_ENABLED ( WHICH )                    &
          BIND(C, NAME='ww_st6_enabled')
       IMPORT :: C_INT
       INTEGER(C_INT), VALUE           :: WHICH
     END FUNCTION WW_ST6_ENABLED
     INTEGER(C_INT) FUNCTION WW_ST6_LAST_ERROR ( )                       &
          BIND(C, NAME='ww_st6_last_error')
       IMPORT :: C_INT
     END FUNCTION WW_ST6_LAST_ERROR
  END INTERFACE
#endif
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
    KOKKOS_SRCE = ( WW_SRCE_ENABLED() == 1_C_INT )
#ifdef W3_ST0
    KOKKOS_SPR0 = ( WW_ST0_ENABLED(1_C_INT) == 1_C_INT )
#endif
#if defined(W3_FLX1) || defined(W3_FLX2) || defined(W3_FLX3) || defined(W3_FLX4) || defined(W3_FLX5)
    KOKKOS_FLX1 = ( WW_FLX_ENABLED(1_C_INT) == 1_C_INT )
    KOKKOS_FLX2 = ( WW_FLX_ENABLED(2_C_INT) == 1_C_INT )
    KOKKOS_FLX3 = ( WW_FLX_ENABLED(3_C_INT) == 1_C_INT )
    KOKKOS_FLX4 = ( WW_FLX_ENABLED(4_C_INT) == 1_C_INT )
    KOKKOS_FLX5 = ( WW_FLX_ENABLED(5_C_INT) == 1_C_INT )
#endif
#ifdef W3_ST6
    KOKKOS_SWL6 = ( WW_SWL_ENABLED(1_C_INT) == 1_C_INT )
#endif
    !
  END SUBROUTINE W3KOKKOS_SETUP
#ifdef W3_ST0
  !/ ------------------------------------------------------------------- /
  !> @brief Give the shim grid IMOD's ST0 scalars and arrays.
  !>
  !> Call from W3INIT after W3IOGR, with W3SETG done for IMOD, so the
  !> W3GDATMD aliases hold this grid's values.
  !>
  SUBROUTINE W3KOKKOS_ST0_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, FTE, FTF, FTWN, SIG, DDEN
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(NR0PAR)
    INTEGER(C_INT)          :: IPAR(NI0PAR)
    !
    RPAR(R0_FTE+1)   = FTE
    RPAR(R0_FTF+1)   = FTF
    RPAR(R0_FTWN+1)  = FTWN
    IPAR(I0_NK+1)    = NK
    IPAR(I0_NTH+1)   = NTH
    IPAR(I0_NSPEC+1) = NSPEC
    CALL ST0_CHECK ( 'WW_ST0_INIT',                                     &
         WW_ST0_INIT ( IMOD, RPAR, NR0PAR, IPAR, NI0PAR ) )
    ! SIG is SIG(0:MK+1) (w3gdatmd.F90:2066), DDEN is DDEN(MK): the shim
    ! takes (1:NK) of each from the LBOUND it is given.
    CALL ST0_SET_R1 ( IMOD, A0_SIG, SIG, LBOUND(SIG) )
    CALL ST0_SET_R1 ( IMOD, A0_DDEN, DDEN, LBOUND(DDEN) )
    CALL ST0_CHECK ( 'WW_ST0_COMMIT', WW_ST0_COMMIT ( IMOD ) )
    !
  END SUBROUTINE W3KOKKOS_ST0_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the model on a non-zero status of the shim.
  SUBROUTINE ST0_CHECK ( WHAT, IERR )
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
  END SUBROUTINE ST0_CHECK
  !/ ------------------------------------------------------------------- /
  !> Register one array. X is assumed-shape, so its own bounds are 1-based:
  !> the actual's LBOUND travels in LB, its SHAPE in SHAPE(X).
  SUBROUTINE ST0_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL ST0_CHECK ( 'WW_ST0_SET_REAL', WW_ST0_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST0_SET_R1
  !/ ------------------------------------------------------------------- /
  !> W3SPR0 on Kokkos for one point, the argument list of W3SPR0.
  !> Serialised: the shim's contexts are unlocked and W3SRCE runs inside
  !> !$OMP PARALLEL under W3_OMPG / W3_OMP0.
  SUBROUTINE W3KOKKOS_SPR0 ( A, CG, WN, EMEAN, FMEAN, WNMEAN, AMAX )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), WN(NK)
    REAL, INTENT(OUT)       :: EMEAN, FMEAN, WNMEAN, AMAX
    REAL(C_FLOAT)           :: EM(1), FM(1), WM(1), AM(1)
    INTEGER(C_INT)          :: IERR
    !$OMP CRITICAL (WW_KOKKOS_ST0)
    CALL WW_SPR0 ( IGRID, 1, A, CG, WN, EM, FM, WM, AM )
    IERR = WW_ST0_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST0)
    CALL ST0_CHECK ( 'W3SPR0', IERR )
    EMEAN = EM(1)
    FMEAN = FM(1)
    WNMEAN = WM(1)
    AMAX = AM(1)
  END SUBROUTINE W3KOKKOS_SPR0
#endif
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
         SINTAILPAR, SSDSP, WWNMEANP, SSTXFTWN, SSDSBRF1, SSDSBINT,    &
         SSDSBCK, SSDSHCK, SSDSABK, SSDSPBK, CAPCHNK, ZZWND, SSDSCOS,  &
         SSDSDTH, SSDSBT, SSDSBM,                                      &
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
    RPAR(R_SSTXFTWN+1)     = SSTXFTWN
    RPAR(R_SSDSBRF1+1)     = SSDSBRF1
    RPAR(R_SSDSBINT+1)     = SSDSBINT
    RPAR(R_SSDSBCK+1)      = SSDSBCK
    RPAR(R_SSDSHCK+1)      = SSDSHCK
    RPAR(R_SSDSABK+1)      = SSDSABK
    RPAR(R_SSDSPBK+1)      = SSDSPBK
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
#if defined(W3_ST4) && defined(W3_NL1)
  !/ ------------------------------------------------------------------- /
  !> @brief Give the shim grid IMOD's W3SRCE scalars and arrays.
  !>
  !> Call from W3INIT after WW_SNL1_INIT and W3KOKKOS_ST4_SETUP for IMOD: the
  !> Kokkos W3SRCE calls the Kokkos DIA and ST4 with those contexts.
  !>
  SUBROUTINE W3KOKKOS_SRCE_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, SIG, DDEN, ECOS, ESIN, MAPWN,  &
         DMIN, DTMAX, DTMIN, FACTI1, FACTI2, FACHFA, FACP, XFLT, XREL,  &
         FHMAX, FFXFM, FFXPM, FFXFA, ICESCALES, IICEDISP, IC_NUMERICS,  &
         SLNC1, FSPM, FSHF, SBTC1, SDBC1, SDBC2, FDONLY, IQTPE
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(NRSPAR)
    INTEGER(C_INT)          :: IPAR(NISPAR)
    INTEGER                 :: J
    !
    RPAR = 0.
    IPAR = 0
    RPAR(RS_DMIN+1)   = DMIN
    RPAR(RS_DTMAX+1)  = DTMAX
    RPAR(RS_DTMIN+1)  = DTMIN
    RPAR(RS_FACTI1+1) = FACTI1
    RPAR(RS_FACTI2+1) = FACTI2
    RPAR(RS_FACHFA+1) = FACHFA
    RPAR(RS_FACP+1)   = FACP
    RPAR(RS_XFLT+1)   = XFLT
    RPAR(RS_XREL+1)   = XREL
    RPAR(RS_FHMAX+1)  = FHMAX
    RPAR(RS_FFXFM+1)  = FFXFM
    RPAR(RS_FFXPM+1)  = FFXPM
    RPAR(RS_FFXFA+1)  = FFXFA
    DO J = 1, 4
      RPAR(RS_ICESCALES+J) = ICESCALES(J)
    END DO
    RPAR(RS_SLNC1+1)  = SLNC1
    RPAR(RS_FSPM+1)   = FSPM
    RPAR(RS_FSHF+1)   = FSHF
    RPAR(RS_SBTC1+1)  = SBTC1
    RPAR(RS_SDBC1+1)  = SDBC1
    RPAR(RS_SDBC2+1)  = SDBC2
    IPAR(IS_NK+1)     = NK
    IPAR(IS_NTH+1)    = NTH
    IPAR(IS_NSPEC+1)  = NSPEC
    IPAR(IS_IQTPE+1)  = IQTPE
    IPAR(IS_IICEDISP+1)    = MERGE(1, 0, IICEDISP)
    IPAR(IS_IC_NUMERICS+1) = MERGE(1, 0, IC_NUMERICS)
    IPAR(IS_FDONLY+1)      = MERGE(1, 0, FDONLY)
    ! The kernel translates W3SRCE as switch_lab_shrd compiles it: these six
    ! switches on, none of the other switches W3SRCE, W3SLN1, W3SBT1, W3SDB1
    ! branch on (W3_OMPG changes no arithmetic and is allowed).
#if defined(W3_BT1) && defined(W3_DB1) && defined(W3_LN1) && defined(W3_MLIM)
#if !defined(W3_BS1) && !defined(W3_BT4) && !defined(W3_BT8) && !defined(W3_BT9) && !defined(W3_DEBUGSRC) && !defined(W3_DEBUGRUN)
#if !defined(W3_FLD1) && !defined(W3_FLD2) && !defined(W3_FLX1) && !defined(W3_FLX2) && !defined(W3_FLX3) && !defined(W3_FLX4) && !defined(W3_FLX5)
#if !defined(W3_IC1) && !defined(W3_IC2) && !defined(W3_IC3) && !defined(W3_IC4) && !defined(W3_IC5) && !defined(W3_IG1) && !defined(W3_IS1) && !defined(W3_IS2)
#if !defined(W3_LN0) && !defined(W3_NL0) && !defined(W3_NL2) && !defined(W3_NL3) && !defined(W3_NL4) && !defined(W3_NL5) && !defined(W3_NLS) && !defined(W3_NNT)
#if !defined(W3_PDLIB) && !defined(W3_REF1) && !defined(W3_S) && !defined(W3_SEED) && !defined(W3_T) && !defined(W3_T0) && !defined(W3_T1) && !defined(W3_TR1) && !defined(W3_UOST)
#if !defined(W3_ST0) && !defined(W3_ST1) && !defined(W3_ST2) && !defined(W3_ST3) && !defined(W3_ST6) && !defined(W3_STAB3)
    IPAR(IS_LAB_SWITCH+1)  = 1
#endif
#endif
#endif
#endif
#endif
#endif
#endif
    !
    CALL SRCE_CHECK ( 'WW_SRCE_INIT', WW_SRCE_INIT ( IMOD, RPAR, NRSPAR, &
         IPAR, NISPAR, SIG(1:NK), DDEN(1:NK), ECOS(1:NSPEC),             &
         ESIN(1:NSPEC), MAPWN(1:NSPEC) ) )
    !
  END SUBROUTINE W3KOKKOS_SRCE_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the model on a non-zero status of the shim.
  SUBROUTINE SRCE_CHECK ( WHAT, IERR )
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
  END SUBROUTINE SRCE_CHECK
  !/ ------------------------------------------------------------------- /
  !> .TRUE. when W3SRCE(srce_direct, ...) runs on Kokkos for this call: the
  !> flag is on and the call takes no branch the kernel leaves to the Fortran:
  !> GQM (IQTPE <= 0, w3srcemd.F90:1259); Liu's ice dispersion (INFLAGS2(4)
  !> and IICEDISP, w3srcemd.F90:2124-2126; ICE > 0 is not needed for the test,
  !> without it the kernel skips the ice block as well); TAUWBUG = 0
  !> (SINTAILPAR(4) <= 0.5) with IT > 0, where W3SPR4 at :1119 reads the local
  !> LLWS before anything sets it (w3src4md.F90:302), which no translation can
  !> reproduce.
  LOGICAL FUNCTION W3KOKKOS_SRCE_ON ( IT )
    USE W3GDATMD, ONLY: IQTPE, IICEDISP, SINTAILPAR
    USE W3IDATMD, ONLY: INFLAGS2
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IT
    W3KOKKOS_SRCE_ON = KOKKOS_SRCE .AND. IQTPE > 0 .AND.                 &
         .NOT. ( INFLAGS2(4) .AND. IICEDISP ) .AND.                      &
         ( IT == 0 .OR. SINTAILPAR(4) > 0.5 )
  END FUNCTION W3KOKKOS_SRCE_ON
  !/ ------------------------------------------------------------------- /
  !> W3SRCE(srce_direct, ...) on Kokkos for one point: the argument list of
  !> W3SRCE as the lab switch reads it (w3srcemd.F90:190-2412). VSIO, VDIO
  !> are zeroed as W3SRCE zeroes them (:879-880); the arguments W3SRCE does not
  !> read under that switch are not passed.
  SUBROUTINE W3KOKKOS_SRCE ( IT, SPEC, VSIO, VDIO, SHAVEIO, WN1, CG1,  &
       D_INP, U10ABS, U10DIR, AS, USTAR, USTDIR, ICE, DTDYN, FCUT, DTG, &
       TAUWX, TAUWY, TAUOX, TAUOY, TAUWIX, TAUWIY, TAUWNX, TAUWNY,      &
       PHIAW, CHARN, TWS, PHIOC, WHITECAP, PHIBBL, TAUBBL, TAUICE,      &
       PHICE, TAUOCX, TAUOCY, WNMEAN, DAIR, COEF )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    USE W3IDATMD, ONLY: INFLAGS2
    USE W3ODATMD, ONLY: FLOGRD
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IT
    REAL, INTENT(INOUT)     :: SPEC(NSPEC)
    REAL, INTENT(OUT)       :: VSIO(NSPEC), VDIO(NSPEC)
    LOGICAL, INTENT(OUT)    :: SHAVEIO
    REAL, INTENT(IN)        :: WN1(NK), CG1(NK), D_INP, U10ABS, U10DIR, &
                               AS, ICE, DTG, DAIR, COEF
    REAL, INTENT(INOUT)     :: USTAR, USTDIR, TAUWX, TAUWY, TAUOX,      &
                               TAUOY, TAUWIX, TAUWIY, TAUWNX, TAUWNY,   &
                               PHIAW, CHARN, TWS, PHIOC, WHITECAP(4),   &
                               PHIBBL, TAUBBL(2), TAUICE(2), PHICE,     &
                               TAUOCX, TAUOCY, WNMEAN
    REAL, INTENT(OUT)       :: DTDYN, FCUT
    INTEGER(C_INT)          :: SH, IERR
    VSIO = 0.
    VDIO = 0.
    SH = 0
    !$OMP CRITICAL (WW_KOKKOS_SRCE)
    CALL WW_SRCE ( IGRID, IT, MERGE(1_C_INT, 0_C_INT, INFLAGS2(4)),     &
         MERGE(1_C_INT, 0_C_INT, FLOGRD(5,7)),                          &
         MERGE(1_C_INT, 0_C_INT, FLOGRD(5,8)),                          &
         SPEC, WN1, CG1, D_INP, U10ABS, U10DIR, AS, ICE, DTG, DAIR,     &
         COEF, USTAR, USTDIR, TAUWX, TAUWY, DTDYN, FCUT, TAUOX, TAUOY,  &
         TAUWIX, TAUWIY, TAUWNX, TAUWNY, PHIAW, CHARN, TWS, PHIOC,      &
         WHITECAP, PHIBBL, TAUBBL, TAUICE, PHICE, TAUOCX, TAUOCY,       &
         WNMEAN, SH )
    IERR = WW_SRCE_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_SRCE)
    CALL SRCE_CHECK ( 'W3SRCE', IERR )
    SHAVEIO = ( SH /= 0_C_INT )
  END SUBROUTINE W3KOKKOS_SRCE
#endif
#ifdef W3_ST1
  !/ ------------------------------------------------------------------- /
  !> @brief Read the ST1 run-time switches and, when one is on, give the
  !> shim grid IMOD's ST1 scalars and arrays.
  !>
  !> Call from W3INIT after W3IOGR and WW_KOKKOS_INIT, with W3SETG done for
  !> IMOD. The flags are read here, not in W3KOKKOS_SETUP, so the ST1 block
  !> stays under one switch.
  !>
  SUBROUTINE W3KOKKOS_ST1_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, FTE, FTF, FTWN, SINC1, SDSC1,  &
         SIG, DDEN, SIG2, ECOS, ESIN
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(ST1_NRPAR)
    INTEGER(C_INT)          :: IPAR(ST1_NIPAR)
    !
    KOKKOS_SPR1 = ( WW_ST1_ENABLED(1_C_INT) == 1_C_INT )
    KOKKOS_SIN1 = ( WW_ST1_ENABLED(2_C_INT) == 1_C_INT )
    KOKKOS_SDS1 = ( WW_ST1_ENABLED(3_C_INT) == 1_C_INT )
    IF ( .NOT. ( KOKKOS_SPR1 .OR. KOKKOS_SIN1 .OR. KOKKOS_SDS1 ) ) RETURN
    !
    RPAR(ST1_R_FTE+1)   = FTE
    RPAR(ST1_R_FTF+1)   = FTF
    RPAR(ST1_R_FTWN+1)  = FTWN
    RPAR(ST1_R_SINC1+1) = SINC1
    RPAR(ST1_R_SDSC1+1) = SDSC1
    IPAR(ST1_I_NK+1)    = NK
    IPAR(ST1_I_NTH+1)   = NTH
    IPAR(ST1_I_NSPEC+1) = NSPEC
    !
    CALL ST1_CHECK ( 'WW_ST1_INIT',                                     &
         WW_ST1_INIT ( IMOD, RPAR, ST1_NRPAR, IPAR, ST1_NIPAR ) )
    ! Whole arrays with their LBOUND: SIG is SIG(0:MK+1) (w3gdatmd.F90:2066);
    ! the shim keeps elements 1..NK or 1..NSPEC.
    CALL ST1_SET_R1 ( IMOD, ST1_A_SIG, SIG, LBOUND(SIG) )
    CALL ST1_SET_R1 ( IMOD, ST1_A_DDEN, DDEN, LBOUND(DDEN) )
    CALL ST1_SET_R1 ( IMOD, ST1_A_SIG2, SIG2, LBOUND(SIG2) )
    CALL ST1_SET_R1 ( IMOD, ST1_A_ECOS, ECOS, LBOUND(ECOS) )
    CALL ST1_SET_R1 ( IMOD, ST1_A_ESIN, ESIN, LBOUND(ESIN) )
    CALL ST1_CHECK ( 'WW_ST1_COMMIT', WW_ST1_COMMIT ( IMOD ) )
    !
  END SUBROUTINE W3KOKKOS_ST1_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the model on a non-zero status of the shim.
  SUBROUTINE ST1_CHECK ( WHAT, IERR )
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
  END SUBROUTINE ST1_CHECK
  !/ ------------------------------------------------------------------- /
  !> Register one array. X is assumed-shape, so its own bounds are 1-based:
  !> the actual's LBOUND travels in LB, its SHAPE in SHAPE(X).
  SUBROUTINE ST1_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL ST1_CHECK ( 'WW_ST1_SET_REAL', WW_ST1_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST1_SET_R1
  !/ ------------------------------------------------------------------- /
  !> W3SPR1 on Kokkos for one point, the argument list of W3SPR1
  !> (w3src1md.F90:87). Serialised: the shim's contexts are unlocked and
  !> W3SRCE runs inside !$OMP PARALLEL under W3_OMPG / W3_OMP0.
  SUBROUTINE W3KOKKOS_SPR1 ( A, CG, WN, EMEAN, FMEAN, WNMEAN, AMAX )
    USE W3GDATMD, ONLY: NK, NTH, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK)
    REAL, INTENT(OUT)       :: EMEAN, FMEAN, WNMEAN, AMAX
    REAL(C_FLOAT)           :: EM(1), FM(1), WM(1), AM(1)
    INTEGER(C_INT)          :: IERR
    !$OMP CRITICAL (WW_KOKKOS_ST1)
    CALL WW_SPR1 ( IGRID, 1, A, CG, WN, EM, FM, WM, AM )
    IERR = WW_ST1_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST1)
    CALL ST1_CHECK ( 'W3SPR1', IERR )
    EMEAN = EM(1)
    FMEAN = FM(1)
    WNMEAN = WM(1)
    AMAX = AM(1)
  END SUBROUTINE W3KOKKOS_SPR1
  !/ ------------------------------------------------------------------- /
  !> W3SIN1 on Kokkos for one point, the argument list of W3SIN1
  !> (w3src1md.F90:255).
  SUBROUTINE W3KOKKOS_SIN1 ( A, K, USTAR, USDIR, S, D )
    USE W3GDATMD, ONLY: NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), K(NSPEC), USTAR, USDIR
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    REAL(C_FLOAT)           :: US(1), UD(1)
    INTEGER(C_INT)          :: IERR
    US(1) = USTAR
    UD(1) = USDIR
    !$OMP CRITICAL (WW_KOKKOS_ST1)
    CALL WW_SIN1 ( IGRID, 1, A, K, US, UD, S, D )
    IERR = WW_ST1_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST1)
    CALL ST1_CHECK ( 'W3SIN1', IERR )
  END SUBROUTINE W3KOKKOS_SIN1
  !/ ------------------------------------------------------------------- /
  !> W3SDS1 on Kokkos for one point, the argument list of W3SDS1
  !> (w3src1md.F90:432).
  SUBROUTINE W3KOKKOS_SDS1 ( A, K, EMEAN, FMEAN, WNMEAN, S, D )
    USE W3GDATMD, ONLY: NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), K(NSPEC), EMEAN, FMEAN, WNMEAN
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    REAL(C_FLOAT)           :: EM(1), FM(1), WM(1)
    INTEGER(C_INT)          :: IERR
    EM(1) = EMEAN
    FM(1) = FMEAN
    WM(1) = WNMEAN
    !$OMP CRITICAL (WW_KOKKOS_ST1)
    CALL WW_SDS1 ( IGRID, 1, A, K, EM, FM, WM, S, D )
    IERR = WW_ST1_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST1)
    CALL ST1_CHECK ( 'W3SDS1', IERR )
  END SUBROUTINE W3KOKKOS_SDS1
#endif
#ifdef W3_ST2
  !/ ------------------------------------------------------------------- /
  !> Read WW_KOKKOS_SPR2/SIN2/SDS2 and, when any is set, configure grid
  !> IMOD of the ST2 shim from W3GDATMD (after W3IOGR). BETATB, DSIGA and
  !> DDRAG are PRIVATE to W3SRC2MD: the shim builds them itself with its
  !> INPTAB at the first commit of the process.
  SUBROUTINE W3KOKKOS_ST2_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, XFR, FTE, FTF, FTWN, FTTR, &
         FPIMIN, ZWIND, FSWELL, FACTI1, FACTI2, XF1, XF2, XFH, SDSALN,  &
         CDSA0, CDSA1, CDSA2, CDSB0, CDSB1, CDSB2, CDSB3,               &
         SIG, SIG2, DDEN, ECOS, ESIN
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(ST2_NRPAR)
    INTEGER(C_INT)          :: IPAR(ST2_NIPAR)
    KOKKOS_SPR2 = ( WW_ST2_ENABLED(1_C_INT) == 1_C_INT )
    KOKKOS_SIN2 = ( WW_ST2_ENABLED(2_C_INT) == 1_C_INT )
    KOKKOS_SDS2 = ( WW_ST2_ENABLED(3_C_INT) == 1_C_INT )
    IF ( .NOT. ( KOKKOS_SPR2 .OR. KOKKOS_SIN2 .OR. KOKKOS_SDS2 ) ) RETURN
    RPAR = 0.
    IPAR = 0
    RPAR(ST2_R_DTH+1)    = DTH
    RPAR(ST2_R_XFR+1)    = XFR
    RPAR(ST2_R_FTE+1)    = FTE
    RPAR(ST2_R_FTF+1)    = FTF
    RPAR(ST2_R_FTWN+1)   = FTWN
    RPAR(ST2_R_FTTR+1)   = FTTR
    RPAR(ST2_R_FPIMIN+1) = FPIMIN
    RPAR(ST2_R_ZWIND+1)  = ZWIND
    RPAR(ST2_R_FSWELL+1) = FSWELL
    RPAR(ST2_R_FACTI1+1) = FACTI1
    RPAR(ST2_R_FACTI2+1) = FACTI2
    RPAR(ST2_R_XF1+1)    = XF1
    RPAR(ST2_R_XF2+1)    = XF2
    RPAR(ST2_R_XFH+1)    = XFH
    RPAR(ST2_R_SDSALN+1) = SDSALN
    RPAR(ST2_R_CDSA0+1)  = CDSA0
    RPAR(ST2_R_CDSA1+1)  = CDSA1
    RPAR(ST2_R_CDSA2+1)  = CDSA2
    RPAR(ST2_R_CDSB0+1)  = CDSB0
    RPAR(ST2_R_CDSB1+1)  = CDSB1
    RPAR(ST2_R_CDSB2+1)  = CDSB2
    RPAR(ST2_R_CDSB3+1)  = CDSB3
    IPAR(ST2_I_NK+1)     = NK
    IPAR(ST2_I_NTH+1)    = NTH
    IPAR(ST2_I_NSPEC+1)  = NSPEC
    CALL ST2_CHECK ( 'WW_ST2_INIT',                                     &
         WW_ST2_INIT ( IMOD, RPAR, ST2_NRPAR, IPAR, ST2_NIPAR ) )
    ! The whole arrays with their LBOUND: SIG is (0:MK+1) and the shim
    ! takes SIG(1:NK+1) from it (W3SDS2 reads SIG(NK+1)).
    CALL ST2_SET_R1 ( IMOD, ST2_A_SIG, SIG, LBOUND(SIG) )
    CALL ST2_SET_R1 ( IMOD, ST2_A_SIG2, SIG2, LBOUND(SIG2) )
    CALL ST2_SET_R1 ( IMOD, ST2_A_DDEN, DDEN, LBOUND(DDEN) )
    CALL ST2_SET_R1 ( IMOD, ST2_A_ECOS, ECOS, LBOUND(ECOS) )
    CALL ST2_SET_R1 ( IMOD, ST2_A_ESIN, ESIN, LBOUND(ESIN) )
    CALL ST2_CHECK ( 'WW_ST2_COMMIT', WW_ST2_COMMIT ( IMOD ) )
  END SUBROUTINE W3KOKKOS_ST2_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the run with the shim's status code, as ST4_CHECK does.
  SUBROUTINE ST2_CHECK ( WHAT, IERR )
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
  END SUBROUTINE ST2_CHECK
  !/ ------------------------------------------------------------------- /
  SUBROUTINE ST2_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL ST2_CHECK ( 'WW_ST2_SET_REAL', WW_ST2_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST2_SET_R1
  !/ ------------------------------------------------------------------- /
  !> W3SPR2 (w3src2md.F90:125) on Kokkos, same dummy list. The shim's
  !> contexts are unlocked and W3SRCE may run inside an OpenMP parallel
  !> region (W3_OMPG): the call and its status read are serialised.
  SUBROUTINE W3KOKKOS_SPR2 ( A, CG, WN, DEPTH, FPI, U, USTAR, EMEAN,   &
       FMEAN, WNMEAN, AMAX, ALFA, FP )
    USE W3GDATMD, ONLY: NK, NTH, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK), DEPTH, FPI,  &
                               U, USTAR
    REAL, INTENT(OUT)       :: EMEAN, FMEAN, WNMEAN, AMAX, ALFA(NK), FP
    REAL(C_FLOAT)           :: DP1(1), FPI1(1), U1(1), US1(1), EM1(1),  &
                               FM1(1), WM1(1), AM1(1), FP1(1)
    INTEGER(C_INT)          :: IERR
    DP1(1) = DEPTH
    FPI1(1) = FPI
    U1(1) = U
    US1(1) = USTAR
    !$OMP CRITICAL (WW_KOKKOS_ST2)
    CALL WW_SPR2 ( IGRID, 1, A, CG, WN, DP1, FPI1, U1, US1, EM1, FM1,  &
         WM1, AM1, ALFA, FP1 )
    IERR = WW_ST2_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST2)
    CALL ST2_CHECK ( 'W3SPR2', IERR )
    EMEAN = EM1(1)
    FMEAN = FM1(1)
    WNMEAN = WM1(1)
    AMAX = AM1(1)
    FP = FP1(1)
  END SUBROUTINE W3KOKKOS_SPR2
  !/ ------------------------------------------------------------------- /
  !> W3SIN2 (w3src2md.F90:308) on Kokkos, same dummy list.
  SUBROUTINE W3KOKKOS_SIN2 ( A, CG, K, U, UDIR, CD, Z0, FPI, S, D )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), K(NSPEC), U, UDIR,    &
                               CD, Z0
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), FPI
    REAL(C_FLOAT)           :: U1(1), UD1(1), CD1(1), Z01(1), FPI1(1)
    INTEGER(C_INT)          :: IERR
    U1(1) = U
    UD1(1) = UDIR
    CD1(1) = CD
    Z01(1) = Z0
    !$OMP CRITICAL (WW_KOKKOS_ST2)
    CALL WW_SIN2 ( IGRID, 1, A, CG, K, U1, UD1, CD1, Z01, FPI1, S, D )
    IERR = WW_ST2_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST2)
    CALL ST2_CHECK ( 'W3SIN2', IERR )
    FPI = FPI1(1)
  END SUBROUTINE W3KOKKOS_SIN2
  !/ ------------------------------------------------------------------- /
  !> W3SDS2 (w3src2md.F90:582) on Kokkos, same dummy list.
  SUBROUTINE W3KOKKOS_SDS2 ( A, CG, K, FPI, USTAR, ALFA, S, D )
    USE W3GDATMD, ONLY: NK, NTH, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), K(NK), FPI, USTAR,   &
                               ALFA(NK)
    REAL, INTENT(OUT)       :: S(NTH,NK), D(NTH,NK)
    REAL(C_FLOAT)           :: FPI1(1), US1(1)
    INTEGER(C_INT)          :: IERR
    FPI1(1) = FPI
    US1(1) = USTAR
    !$OMP CRITICAL (WW_KOKKOS_ST2)
    CALL WW_SDS2 ( IGRID, 1, A, CG, K, FPI1, US1, ALFA, S, D )
    IERR = WW_ST2_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST2)
    CALL ST2_CHECK ( 'W3SDS2', IERR )
  END SUBROUTINE W3KOKKOS_SDS2
#endif
#ifdef W3_ST3
  !/ ------------------------------------------------------------------- /
  !> Read WW_KOKKOS_SPR3/SIN3/SDS3 and, when any is set, configure grid
  !> IMOD of the ST3 shim from W3GDATMD (after W3IOGR). The first commit
  !> of the process builds TAUT, TAUHFT and their steps with the shim's
  !> INSIN3 from this grid, as W3IOGR's first read does for W3SRC3MD.
  SUBROUTINE W3KOKKOS_ST3_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, ZZWND, AALPHA, BBETA,      &
         ZZ0MAX, ZZALP, SSINTHP, SSWELLF, SSDSC1, DDELTA1, DDELTA2,     &
         WWNMEANP, WWNMEANPTAIL, FTE, SSTXFTF, SSTXFTWN, SSTXFTFTAIL,   &
         SIG, SIG2, DDEN, DDEN2, ECOS, ESIN
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(ST3_NRPAR)
    INTEGER(C_INT)          :: IPAR(ST3_NIPAR)
    KOKKOS_SPR3 = ( WW_ST3_ENABLED(1_C_INT) == 1_C_INT )
    KOKKOS_SIN3 = ( WW_ST3_ENABLED(2_C_INT) == 1_C_INT )
    KOKKOS_SDS3 = ( WW_ST3_ENABLED(3_C_INT) == 1_C_INT )
    IF ( .NOT. ( KOKKOS_SPR3 .OR. KOKKOS_SIN3 .OR. KOKKOS_SDS3 ) ) RETURN
    RPAR = 0.
    IPAR = 0
    RPAR(ST3_R_DTH+1)          = DTH
    RPAR(ST3_R_ZZWND+1)        = ZZWND
    RPAR(ST3_R_AALPHA+1)       = AALPHA
    RPAR(ST3_R_BBETA+1)        = BBETA
    RPAR(ST3_R_ZZ0MAX+1)       = ZZ0MAX
    RPAR(ST3_R_ZZALP+1)        = ZZALP
    RPAR(ST3_R_SSINTHP+1)      = SSINTHP
    RPAR(ST3_R_SSWELLF1+1)     = SSWELLF(1)
    RPAR(ST3_R_SSDSC1+1)       = SSDSC1
    RPAR(ST3_R_DDELTA1+1)      = DDELTA1
    RPAR(ST3_R_DDELTA2+1)      = DDELTA2
    RPAR(ST3_R_WWNMEANP+1)     = WWNMEANP
    RPAR(ST3_R_WWNMEANPTAIL+1) = WWNMEANPTAIL
    RPAR(ST3_R_FTE+1)          = FTE
    RPAR(ST3_R_SSTXFTF+1)      = SSTXFTF
    RPAR(ST3_R_SSTXFTWN+1)     = SSTXFTWN
    RPAR(ST3_R_SSTXFTFTAIL+1)  = SSTXFTFTAIL
    IPAR(ST3_I_NK+1)           = NK
    IPAR(ST3_I_NTH+1)          = NTH
    IPAR(ST3_I_NSPEC+1)        = NSPEC
    ! The shim is one object for every switch file: W3SIN3's STAB3 blocks
    ! are chosen by this flag at run time.
#ifdef W3_STAB3
    IPAR(ST3_I_STAB3+1)        = 1
#endif
    CALL ST3_CHECK ( 'WW_ST3_INIT',                                     &
         WW_ST3_INIT ( IMOD, RPAR, ST3_NRPAR, IPAR, ST3_NIPAR ) )
    ! The whole arrays with their LBOUND: SIG is (0:MK+1) and the shim
    ! takes SIG(1:NK) from it.
    CALL ST3_SET_R1 ( IMOD, ST3_A_SIG, SIG, LBOUND(SIG) )
    CALL ST3_SET_R1 ( IMOD, ST3_A_SIG2, SIG2, LBOUND(SIG2) )
    CALL ST3_SET_R1 ( IMOD, ST3_A_DDEN, DDEN, LBOUND(DDEN) )
    CALL ST3_SET_R1 ( IMOD, ST3_A_DDEN2, DDEN2, LBOUND(DDEN2) )
    CALL ST3_SET_R1 ( IMOD, ST3_A_ECOS, ECOS, LBOUND(ECOS) )
    CALL ST3_SET_R1 ( IMOD, ST3_A_ESIN, ESIN, LBOUND(ESIN) )
    CALL ST3_CHECK ( 'WW_ST3_COMMIT', WW_ST3_COMMIT ( IMOD ) )
  END SUBROUTINE W3KOKKOS_ST3_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the run with the shim's status code, as ST4_CHECK does.
  SUBROUTINE ST3_CHECK ( WHAT, IERR )
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
  END SUBROUTINE ST3_CHECK
  !/ ------------------------------------------------------------------- /
  SUBROUTINE ST3_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL ST3_CHECK ( 'WW_ST3_SET_REAL', WW_ST3_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST3_SET_R1
  !/ ------------------------------------------------------------------- /
  !> W3SPR3 (w3src3md.F90:134) on Kokkos, same dummy list. The shim's
  !> contexts are unlocked and W3SRCE may run inside an OpenMP parallel
  !> region (W3_OMPG): the call and its status read are serialised.
  SUBROUTINE W3KOKKOS_SPR3 ( A, CG, WN, EMEAN, FMEAN, FMEANS, WNMEAN,  &
       AMAX, U, UDIR, USTAR, USDIR, TAUWX, TAUWY, CD, Z0, CHARN, LLWS,  &
       FMEANWS )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK), U, UDIR
    REAL, INTENT(IN)        :: TAUWX, TAUWY
    LOGICAL, INTENT(IN)     :: LLWS(NSPEC)
    REAL, INTENT(INOUT)     :: USTAR, USDIR
    REAL, INTENT(OUT)       :: EMEAN, FMEAN, FMEANS, WNMEAN, AMAX,      &
                               CD, Z0, CHARN, FMEANWS
    REAL(C_FLOAT)           :: U1(1), UD1(1), US1(1), USD1(1), TX1(1),  &
                               TY1(1), EM1(1), FM1(1), FMS1(1), WM1(1), &
                               AM1(1), CD1(1), Z01(1), CH1(1), FWS1(1)
    INTEGER(C_INT)          :: LW(NSPEC), IERR
    U1(1) = U
    UD1(1) = UDIR
    US1(1) = USTAR
    USD1(1) = USDIR
    TX1(1) = TAUWX
    TY1(1) = TAUWY
    LW = MERGE(1_C_INT, 0_C_INT, LLWS)
    !$OMP CRITICAL (WW_KOKKOS_ST3)
    CALL WW_SPR3 ( IGRID, 1, A, CG, WN, EM1, FM1, FMS1, WM1, AM1, U1,  &
         UD1, US1, USD1, TX1, TY1, CD1, Z01, CH1, LW, FWS1 )
    IERR = WW_ST3_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST3)
    CALL ST3_CHECK ( 'W3SPR3', IERR )
    USTAR = US1(1)
    USDIR = USD1(1)
    EMEAN = EM1(1)
    FMEAN = FM1(1)
    FMEANS = FMS1(1)
    WNMEAN = WM1(1)
    AMAX = AM1(1)
    CD = CD1(1)
    Z0 = Z01(1)
    CHARN = CH1(1)
    FMEANWS = FWS1(1)
  END SUBROUTINE W3KOKKOS_SPR3
  !/ ------------------------------------------------------------------- /
  !> W3SIN3 (w3src3md.F90:384) on Kokkos, same dummy list. CD, ICE, IX
  !> and IY are read by no statement of W3SIN3 and do not cross.
  SUBROUTINE W3KOKKOS_SIN3 ( A, CG, K, U, USTAR, DRAT, AS, USDIR, Z0,  &
       CD, TAUWX, TAUWY, TAUWNX, TAUWNY, ICE, S, D, LLWS, IX, IY )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC)
    REAL, INTENT(IN)        :: CG(NK), K(NSPEC), Z0, U, CD
    REAL, INTENT(IN)        :: USTAR, USDIR, AS, DRAT, ICE
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC), TAUWX, TAUWY,        &
                               TAUWNX, TAUWNY
    LOGICAL, INTENT(OUT)    :: LLWS(NSPEC)
    INTEGER, INTENT(IN)     :: IX, IY
    REAL(C_FLOAT)           :: U1(1), US1(1), DR1(1), AS1(1), UD1(1),   &
                               Z01(1), TX1(1), TY1(1), TNX1(1), TNY1(1)
    INTEGER(C_INT)          :: LW(NSPEC), IERR
    U1(1) = U
    US1(1) = USTAR
    DR1(1) = DRAT
    AS1(1) = AS
    UD1(1) = USDIR
    Z01(1) = Z0
    !$OMP CRITICAL (WW_KOKKOS_ST3)
    CALL WW_SIN3 ( IGRID, 1, A, CG, K, U1, US1, DR1, AS1, UD1, Z01,     &
         TX1, TY1, TNX1, TNY1, S, D, LW )
    IERR = WW_ST3_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST3)
    CALL ST3_CHECK ( 'W3SIN3', IERR )
    TAUWX = TX1(1)
    TAUWY = TY1(1)
    TAUWNX = TNX1(1)
    TAUWNY = TNY1(1)
    LLWS = ( LW /= 0_C_INT )
  END SUBROUTINE W3KOKKOS_SIN3
  !/ ------------------------------------------------------------------- /
  !> W3SDS3 (w3src3md.F90:1253) on Kokkos, same dummy list. CG, USTAR,
  !> USDIR, DEPTH, IX and IY are read by no statement of W3SDS3 and do
  !> not cross.
  SUBROUTINE W3KOKKOS_SDS3 ( A, K, CG, EMEAN, FMEAN, WNMEAN, USTAR,    &
       USDIR, DEPTH, S, D, IX, IY )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), K(NK), CG(NK),                 &
                               DEPTH, USTAR, USDIR, EMEAN, FMEAN, WNMEAN
    INTEGER, INTENT(IN)     :: IX, IY
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    REAL(C_FLOAT)           :: EM1(1), FM1(1), WM1(1)
    INTEGER(C_INT)          :: IERR
    EM1(1) = EMEAN
    FM1(1) = FMEAN
    WM1(1) = WNMEAN
    !$OMP CRITICAL (WW_KOKKOS_ST3)
    CALL WW_SDS3 ( IGRID, 1, A, K, EM1, FM1, WM1, S, D )
    IERR = WW_ST3_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST3)
    CALL ST3_CHECK ( 'W3SDS3', IERR )
  END SUBROUTINE W3KOKKOS_SDS3
#endif
#if defined(W3_FLX1) || defined(W3_FLX2) || defined(W3_FLX3) || defined(W3_FLX4) || defined(W3_FLX5)
  !/ ------------------------------------------------------------------- /
  !> @brief Give the shim grid IMOD's W3FLXn scalars (and W3DISPMD's table).
  !>
  !> Call from W3INIT after W3IOGR, with W3SETG done for IMOD, so the W3GDATMD
  !> pointers NITTIN, CINXSI, CD_MAX, CAP_ID, FLX4A0 are this grid's and
  !> W3IOGR has run DISTAB (w3iogrmd.F90:2015-2021). Does nothing unless the
  !> flag of the compiled W3_FLXn is on.
  !>
  SUBROUTINE W3KOKKOS_FLX_SETUP ( IMOD )
#if defined(W3_FLX2)
    USE W3GDATMD, ONLY: NITTIN, CINXSI
#endif
#if defined(W3_FLX3)
    USE W3GDATMD, ONLY: NITTIN, CINXSI, CD_MAX, CAP_ID
#endif
#if defined(W3_FLX4)
    USE W3GDATMD, ONLY: FLX4A0
#endif
#if defined(W3_FLX2) || defined(W3_FLX3)
    USE W3DISPMD, ONLY: EWN1, DSIE, N1MAX
#endif
    USE W3ODATMD, ONLY: NDSO
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(NFRPAR)
    INTEGER(C_INT)          :: IPAR(NFIPAR)
    LOGICAL                 :: ON
    INTEGER                 :: N
    !
    RPAR = 0.
    IPAR = 0
#ifdef W3_FLX1
    ON = KOKKOS_FLX1
    N  = 1
#endif
#ifdef W3_FLX2
    ON = KOKKOS_FLX2
    N  = 2
    IPAR(FI_NITTIN+1) = NITTIN
    RPAR(FR_CINXSI+1) = CINXSI
#endif
#ifdef W3_FLX3
    ON = KOKKOS_FLX3
    N  = 3
    IPAR(FI_NITTIN+1) = NITTIN
    IPAR(FI_CAP_ID+1) = CAP_ID
    RPAR(FR_CINXSI+1) = CINXSI
    RPAR(FR_CD_MAX+1) = CD_MAX
#endif
#ifdef W3_FLX4
    ON = KOKKOS_FLX4
    N  = 4
    RPAR(FR_FLX4A0+1) = FLX4A0
#endif
#ifdef W3_FLX5
    ON = KOKKOS_FLX5
    N  = 5
#endif
    IF ( .NOT. ON ) RETURN
    IPAR(FI_FLX+1) = N
#if defined(W3_FLX2) || defined(W3_FLX3)
    RPAR(FR_DSIE+1)  = DSIE
    IPAR(FI_N1MAX+1) = N1MAX
#endif
    CALL FLX_CHECK ( 'WW_FLX_INIT',                                     &
         WW_FLX_INIT ( IMOD, RPAR, NFRPAR, IPAR, NFIPAR ) )
#if defined(W3_FLX2) || defined(W3_FLX3)
    ! EWN1 is EWN1(0:NAR1D), a fixed-size module array: the bare name passes
    ! EWN1(0), which is what the shim expects with LBOUND 0.
    CALL FLX_CHECK ( 'WW_FLX_SET_REAL', WW_FLX_SET_REAL ( IMOD, FA_EWN1, &
         EWN1, 1_C_INT, INT(LBOUND(EWN1),C_INT),                        &
         INT(SHAPE(EWN1),C_INT) ) )
#endif
    CALL FLX_CHECK ( 'WW_FLX_COMMIT', WW_FLX_COMMIT ( IMOD ) )
    WRITE (NDSO,'(A,I1,A)') '  Kokkos W3FLX', N, ' enabled.'
    !
  END SUBROUTINE W3KOKKOS_FLX_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the model on a non-zero status of the shim.
  SUBROUTINE FLX_CHECK ( WHAT, IERR )
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
  END SUBROUTINE FLX_CHECK
#endif
#ifdef W3_FLX1
  !/ ------------------------------------------------------------------- /
  !> W3FLX1 on Kokkos for one point, the argument list of W3FLX1. Where
  !> W3FLX1 stops (w3flx1md.F90:178-181) this writes its message and stops
  !> with the same code. Serialised: W3SRCE runs inside !$OMP PARALLEL
  !> under W3_OMPG / W3_OMP0 and the shim's contexts are unlocked.
  SUBROUTINE W3KOKKOS_FLX1 ( ZWND, U10, U10D, UST, USTD, Z0, CD )
    USE W3GDATMD, ONLY: IGRID
    USE W3ODATMD, ONLY: NDSE, IAPROC, NAPERR
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    REAL, INTENT(IN)        :: ZWND, U10, U10D
    REAL, INTENT(INOUT)     :: UST, USTD, Z0, CD
    REAL(C_FLOAT)           :: ZW(1), U1(1), UD(1), US(1), USD(1),      &
                               Z01(1), CD1(1)
    INTEGER(C_INT)          :: IEXT(1), IERR
    ZW(1) = ZWND
    U1(1) = U10
    UD(1) = U10D
    US(1) = UST
    USD(1) = USTD
    Z01(1) = Z0
    CD1(1) = CD
    IEXT(1) = 0
    !$OMP CRITICAL (WW_KOKKOS_FLX)
    CALL WW_FLX1 ( IGRID, 1, ZW, U1, UD, US, USD, Z01, CD1, IEXT )
    IERR = WW_FLX_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_FLX)
    CALL FLX_CHECK ( 'W3FLX1', IERR )
    IF ( IEXT(1) /= 0_C_INT ) THEN
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) ZWND
      CALL EXTCDE ( INT(IEXT(1)) )
    END IF
    UST = US(1)
    USTD = USD(1)
    Z0 = Z01(1)
    CD = CD1(1)
    ! The text of w3flx1md.F90:194-196, as W3FLX1 writes it.
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3STR1 : '/               &
         '     HIGHT OF WIND SHOULD BE 10m IN THIS APPRACH '/   &
         '     ZWND =',F8.2,'m'/)
  END SUBROUTINE W3KOKKOS_FLX1
#endif
#ifdef W3_FLX2
  !/ ------------------------------------------------------------------- /
  !> W3FLX2 on Kokkos for one point, the argument list of W3FLX2.
  SUBROUTINE W3KOKKOS_FLX2 ( ZWIND, DEPTH, FP, U, UDIR, UST, USTD, Z0, &
       CD )
    USE W3GDATMD, ONLY: IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: ZWIND, DEPTH, FP, U, UDIR
    REAL, INTENT(INOUT)     :: UST, USTD, Z0, CD
    REAL(C_FLOAT)           :: ZW(1), DP(1), FP1(1), U1(1), UD(1),      &
                               US(1), USD(1), Z01(1), CD1(1)
    INTEGER(C_INT)          :: IERR
    ZW(1) = ZWIND
    DP(1) = DEPTH
    FP1(1) = FP
    U1(1) = U
    UD(1) = UDIR
    US(1) = UST
    USD(1) = USTD
    Z01(1) = Z0
    CD1(1) = CD
    !$OMP CRITICAL (WW_KOKKOS_FLX)
    CALL WW_FLX2 ( IGRID, 1, ZW, DP, FP1, U1, UD, US, USD, Z01, CD1 )
    IERR = WW_FLX_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_FLX)
    CALL FLX_CHECK ( 'W3FLX2', IERR )
    UST = US(1)
    USTD = USD(1)
    Z0 = Z01(1)
    CD = CD1(1)
  END SUBROUTINE W3KOKKOS_FLX2
#endif
#ifdef W3_FLX3
  !/ ------------------------------------------------------------------- /
  !> W3FLX3 on Kokkos for one point, the argument list of W3FLX3.
  SUBROUTINE W3KOKKOS_FLX3 ( ZWIND, DEPTH, FP, U, UDIR, UST, USTD, Z0, &
       CD )
    USE W3GDATMD, ONLY: IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: ZWIND, DEPTH, FP, U, UDIR
    REAL, INTENT(INOUT)     :: UST, USTD, Z0, CD
    REAL(C_FLOAT)           :: ZW(1), DP(1), FP1(1), U1(1), UD(1),      &
                               US(1), USD(1), Z01(1), CD1(1)
    INTEGER(C_INT)          :: IERR
    ZW(1) = ZWIND
    DP(1) = DEPTH
    FP1(1) = FP
    U1(1) = U
    UD(1) = UDIR
    US(1) = UST
    USD(1) = USTD
    Z01(1) = Z0
    CD1(1) = CD
    !$OMP CRITICAL (WW_KOKKOS_FLX)
    CALL WW_FLX3 ( IGRID, 1, ZW, DP, FP1, U1, UD, US, USD, Z01, CD1 )
    IERR = WW_FLX_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_FLX)
    CALL FLX_CHECK ( 'W3FLX3', IERR )
    UST = US(1)
    USTD = USD(1)
    Z0 = Z01(1)
    CD = CD1(1)
  END SUBROUTINE W3KOKKOS_FLX3
#endif
#ifdef W3_FLX4
  !/ ------------------------------------------------------------------- /
  !> W3FLX4 on Kokkos for one point, the argument list of W3FLX4. Where
  !> W3FLX4 stops (w3flx4md.F90:199-202) this writes its message and stops
  !> with the same code.
  SUBROUTINE W3KOKKOS_FLX4 ( ZWND, U10, U10D, UST, USTD, Z0, CD )
    USE W3GDATMD, ONLY: IGRID
    USE W3ODATMD, ONLY: NDSE, IAPROC, NAPERR
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    REAL, INTENT(IN)        :: ZWND, U10, U10D
    REAL, INTENT(INOUT)     :: UST, USTD, Z0, CD
    REAL(C_FLOAT)           :: ZW(1), U1(1), UD(1), US(1), USD(1),      &
                               Z01(1), CD1(1)
    INTEGER(C_INT)          :: IEXT(1), IERR
    ZW(1) = ZWND
    U1(1) = U10
    UD(1) = U10D
    US(1) = UST
    USD(1) = USTD
    Z01(1) = Z0
    CD1(1) = CD
    IEXT(1) = 0
    !$OMP CRITICAL (WW_KOKKOS_FLX)
    CALL WW_FLX4 ( IGRID, 1, ZW, U1, UD, US, USD, Z01, CD1, IEXT )
    IERR = WW_FLX_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_FLX)
    CALL FLX_CHECK ( 'W3FLX4', IERR )
    IF ( IEXT(1) /= 0_C_INT ) THEN
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) ZWND
      CALL EXTCDE ( INT(IEXT(1)) )
    END IF
    UST = US(1)
    USTD = USD(1)
    Z0 = Z01(1)
    CD = CD1(1)
    ! The text of w3flx4md.F90:225-227, as W3FLX4 writes it.
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLX4 : '/               &
         '     HIGHT OF WIND SHOULD BE 10m IN THIS APPROACH '/   &
         '     ZWND =',F8.2,'m'/)
  END SUBROUTINE W3KOKKOS_FLX4
#endif
#ifdef W3_FLX5
  !/ ------------------------------------------------------------------- /
  !> W3FLX5 on Kokkos for one point, the argument list of W3FLX5. Where
  !> W3FLX5 stops (w3flx5md.F90:228-231) this writes its message and stops
  !> with the same code.
  SUBROUTINE W3KOKKOS_FLX5 ( ZWND, U10, U10D, TAUA, TAUADIR, RHOAIR,   &
       UST, USTD, Z0, CD, CHARN )
    USE W3GDATMD, ONLY: IGRID
    USE W3ODATMD, ONLY: NDSE, IAPROC, NAPERR
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    REAL, INTENT(IN)        :: ZWND, U10, U10D, TAUA, TAUADIR, RHOAIR
    REAL, INTENT(INOUT)     :: UST, USTD, Z0, CD, CHARN
    REAL(C_FLOAT)           :: ZW(1), U1(1), UD(1), TA(1), TD(1),       &
                               RA(1), US(1), USD(1), Z01(1), CD1(1),    &
                               CH(1)
    INTEGER(C_INT)          :: IEXT(1), IERR
    ZW(1) = ZWND
    U1(1) = U10
    UD(1) = U10D
    TA(1) = TAUA
    TD(1) = TAUADIR
    RA(1) = RHOAIR
    US(1) = UST
    USD(1) = USTD
    Z01(1) = Z0
    CD1(1) = CD
    CH(1) = CHARN
    IEXT(1) = 0
    !$OMP CRITICAL (WW_KOKKOS_FLX)
    CALL WW_FLX5 ( IGRID, 1, ZW, U1, UD, TA, TD, RA, US, USD, Z01, CD1, &
         CH, IEXT )
    IERR = WW_FLX_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_FLX)
    CALL FLX_CHECK ( 'W3FLX5', IERR )
    IF ( IEXT(1) /= 0_C_INT ) THEN
      IF ( IAPROC .EQ. NAPERR ) WRITE (NDSE,1000) ZWND
      CALL EXTCDE ( INT(IEXT(1)) )
    END IF
    UST = US(1)
    USTD = USD(1)
    Z0 = Z01(1)
    CD = CD1(1)
    CHARN = CH(1)
    ! The text of w3flx5md.F90:254-256, as W3FLX5 writes it.
1000 FORMAT (/' *** WAVEWATCH III ERROR IN W3FLX5 : '/                 &
         '     HEIGHT OF WIND SHOULD BE 10m IN THIS APPROACH '/   &
         '     ZWND =',F8.2,'m'/)
  END SUBROUTINE W3KOKKOS_FLX5
#endif
#ifdef W3_ST6
  !/ ------------------------------------------------------------------- /
  !> @brief Give the shim grid IMOD's W3SWL6 scalars and arrays.
  !>
  !> Call from W3INIT after W3IOGR, with W3SETG done for IMOD, so the
  !> W3GDATMD aliases hold this grid's values.
  !>
  SUBROUTINE W3KOKKOS_SWL_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, SIG, DDEN, SWL6B1,       &
         SWL6CSTB1
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(SWL_NRPAR)
    INTEGER(C_INT)          :: IPAR(SWL_NIPAR)
    !
    RPAR(SWL_R_DTH+1)       = DTH
    RPAR(SWL_R_SWL6B1+1)    = SWL6B1
    IPAR(SWL_I_NK+1)        = NK
    IPAR(SWL_I_NTH+1)       = NTH
    IPAR(SWL_I_NSPEC+1)     = NSPEC
    IPAR(SWL_I_SWL6CSTB1+1) = MERGE(1_C_INT, 0_C_INT, SWL6CSTB1)
    CALL SWL_CHECK ( 'WW_SWL_INIT',                                     &
         WW_SWL_INIT ( IMOD, RPAR, SWL_NRPAR, IPAR, SWL_NIPAR ) )
    ! SIG is SIG(0:MK+1) (w3gdatmd.F90:2066), DDEN is DDEN(MK): the shim
    ! takes (1:NK) of each from the LBOUND it is given.
    CALL SWL_SET_R1 ( IMOD, SWL_A_SIG, SIG, LBOUND(SIG) )
    CALL SWL_SET_R1 ( IMOD, SWL_A_DDEN, DDEN, LBOUND(DDEN) )
    CALL SWL_CHECK ( 'WW_SWL_COMMIT', WW_SWL_COMMIT ( IMOD ) )
    !
  END SUBROUTINE W3KOKKOS_SWL_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the model on a non-zero status of the shim.
  SUBROUTINE SWL_CHECK ( WHAT, IERR )
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
  END SUBROUTINE SWL_CHECK
  !/ ------------------------------------------------------------------- /
  !> Register one array. X is assumed-shape, so its own bounds are 1-based:
  !> the actual's LBOUND travels in LB, its SHAPE in SHAPE(X).
  SUBROUTINE SWL_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL SWL_CHECK ( 'WW_SWL_SET_REAL', WW_SWL_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE SWL_SET_R1
  !/ ------------------------------------------------------------------- /
  !> W3SWL6 on Kokkos for one point, the argument list of W3SWL6.
  !> Serialised: the shim's contexts and its one IKN per process are
  !> unlocked, and W3SRCE runs inside !$OMP PARALLEL under W3_OMPG / W3_OMP0.
  !> Where W3SWL6 stops (SIZE(IKN) /= NK, w3swldmd.F90:374-377) this writes
  !> its FORMAT 6001 and stops with the same EXTCDE(2).
  SUBROUTINE W3KOKKOS_SWL6 ( A, CG, WN, S, D )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), WN(NK)
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    INTEGER(C_INT)          :: RC, NIKN, IERR
    !$OMP CRITICAL (WW_KOKKOS_SWL)
    RC = WW_SWL6 ( IGRID, 1, A, CG, WN, S, D, NIKN )
    IERR = WW_SWL_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_SWL)
    CALL SWL_CHECK ( 'W3SWL6', IERR )
    IF ( RC /= 0_C_INT ) THEN
      WRITE (NDSE,6001) NIKN, NK
      CALL EXTCDE ( 2 )
    END IF
6001 FORMAT (' *** WAVEWATCH III ERROR IN W3SWL6 : '/ &
       '     INDEX ARRAY DIM MISMATCH SIZE(IKN)=',I4,' .NE. NK=',I4)
  END SUBROUTINE W3KOKKOS_SWL6
#endif
#ifdef W3_ST6
  !/ ------------------------------------------------------------------- /
  !> Read WW_KOKKOS_SPR6/SIN6/SDS6 and, when any is set, configure grid
  !> IMOD of the ST6 shim from W3GDATMD (after W3IOGR, W3SETG done for
  !> IMOD). The SAVE tables of W3SRC6MD (IKN, MK10Hz, SIG10Hz, DSII10Hz)
  !> are not built here: the shim builds them at the first W3SIN6/W3SDS6
  !> call of the process, from that call's grid, as W3SRC6MD does.
  SUBROUTINE W3KOKKOS_ST6_SETUP ( IMOD )
    USE W3GDATMD, ONLY: NK, NTH, NSPEC, DTH, XFR, FTE, FTF, FTWN,       &
         SIN6A0, SIN6WS, SIN6CHKMIN, SIN6CHKCAP, SIN6CHKINF,            &
         SIN6CHKSIG, SIN6FLCAP, SDS6A1, SDS6A2, SDS6P1, SDS6P2, SDS6ET, &
         SIG, DDEN, DSII, SIG2, DDEN2, ECOS, ESIN
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    REAL(C_FLOAT)           :: RPAR(ST6_NRPAR)
    INTEGER(C_INT)          :: IPAR(ST6_NIPAR)
    KOKKOS_SPR6 = ( WW_ST6_ENABLED(1_C_INT) == 1_C_INT )
    KOKKOS_SIN6 = ( WW_ST6_ENABLED(2_C_INT) == 1_C_INT )
    KOKKOS_SDS6 = ( WW_ST6_ENABLED(3_C_INT) == 1_C_INT )
    IF ( .NOT. ( KOKKOS_SPR6 .OR. KOKKOS_SIN6 .OR. KOKKOS_SDS6 ) ) RETURN
    RPAR(ST6_R_DTH+1)         = DTH
    RPAR(ST6_R_XFR+1)         = XFR
    RPAR(ST6_R_FTE+1)         = FTE
    RPAR(ST6_R_FTF+1)         = FTF
    RPAR(ST6_R_FTWN+1)        = FTWN
    RPAR(ST6_R_SIN6A0+1)      = SIN6A0
    RPAR(ST6_R_SIN6WS+1)      = SIN6WS
    RPAR(ST6_R_SIN6CHKMIN+1)  = SIN6CHKMIN
    RPAR(ST6_R_SIN6CHKCAP+1)  = SIN6CHKCAP
    RPAR(ST6_R_SIN6CHKINF+1)  = SIN6CHKINF
    RPAR(ST6_R_SIN6CHKSIG+1)  = SIN6CHKSIG
    RPAR(ST6_R_SDS6A1+1)      = SDS6A1
    RPAR(ST6_R_SDS6A2+1)      = SDS6A2
    IPAR(ST6_I_NK+1)          = NK
    IPAR(ST6_I_NTH+1)         = NTH
    IPAR(ST6_I_NSPEC+1)       = NSPEC
    IPAR(ST6_I_SIN6FLCAP+1)   = MERGE(1_C_INT, 0_C_INT, SIN6FLCAP)
    IPAR(ST6_I_SDS6P1+1)      = SDS6P1
    IPAR(ST6_I_SDS6P2+1)      = SDS6P2
    IPAR(ST6_I_SDS6ET+1)      = MERGE(1_C_INT, 0_C_INT, SDS6ET)
    CALL ST6_CHECK ( 'WW_ST6_INIT',                                     &
         WW_ST6_INIT ( IMOD, RPAR, ST6_NRPAR, IPAR, ST6_NIPAR ) )
    ! The whole arrays with their LBOUND: SIG is (0:MK+1), the others are
    ! allocated for MK, MSPEC (w3gdatmd.F90:2062-2069); the shim takes
    ! (1:NK) or (1:NSPEC) of each.
    CALL ST6_SET_R1 ( IMOD, ST6_A_SIG, SIG, LBOUND(SIG) )
    CALL ST6_SET_R1 ( IMOD, ST6_A_DDEN, DDEN, LBOUND(DDEN) )
    CALL ST6_SET_R1 ( IMOD, ST6_A_DSII, DSII, LBOUND(DSII) )
    CALL ST6_SET_R1 ( IMOD, ST6_A_SIG2, SIG2, LBOUND(SIG2) )
    CALL ST6_SET_R1 ( IMOD, ST6_A_DDEN2, DDEN2, LBOUND(DDEN2) )
    CALL ST6_SET_R1 ( IMOD, ST6_A_ECOS, ECOS, LBOUND(ECOS) )
    CALL ST6_SET_R1 ( IMOD, ST6_A_ESIN, ESIN, LBOUND(ESIN) )
    CALL ST6_CHECK ( 'WW_ST6_COMMIT', WW_ST6_COMMIT ( IMOD ) )
  END SUBROUTINE W3KOKKOS_ST6_SETUP
  !/ ------------------------------------------------------------------- /
  !> Stop the run with the shim's status code, as ST4_CHECK does.
  SUBROUTINE ST6_CHECK ( WHAT, IERR )
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
  END SUBROUTINE ST6_CHECK
  !/ ------------------------------------------------------------------- /
  SUBROUTINE ST6_SET_R1 ( IMOD, ID, X, LB )
    INTEGER, INTENT(IN)     :: IMOD, ID, LB(1)
    REAL, INTENT(IN)        :: X(:)
    CALL ST6_CHECK ( 'WW_ST6_SET_REAL', WW_ST6_SET_REAL ( IMOD, ID, X, &
         1, INT(LB,C_INT), INT(SHAPE(X),C_INT) ) )
  END SUBROUTINE ST6_SET_R1
  !/ ------------------------------------------------------------------- /
  !> W3SPR6 (w3src6md.F90:129) on Kokkos, same dummy list. The shim's
  !> contexts are unlocked and W3SRCE may run inside an OpenMP parallel
  !> region (W3_OMPG): the call and its status read are serialised.
  SUBROUTINE W3KOKKOS_SPR6 ( A, CG, WN, EMEAN, FMEAN, WNMEAN, AMAX, FP )
    USE W3GDATMD, ONLY: NK, NTH, IGRID
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NTH,NK), CG(NK), WN(NK)
    REAL, INTENT(OUT)       :: EMEAN, FMEAN, WNMEAN, AMAX, FP
    REAL(C_FLOAT)           :: EM1(1), FM1(1), WM1(1), AM1(1), FP1(1)
    INTEGER(C_INT)          :: IERR
    !$OMP CRITICAL (WW_KOKKOS_ST6)
    CALL WW_SPR6 ( IGRID, 1, A, CG, WN, EM1, FM1, WM1, AM1, FP1 )
    IERR = WW_ST6_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST6)
    CALL ST6_CHECK ( 'W3SPR6', IERR )
    EMEAN = EM1(1)
    FMEAN = FM1(1)
    WNMEAN = WM1(1)
    AMAX = AM1(1)
    FP = FP1(1)
  END SUBROUTINE W3KOKKOS_SPR6
  !/ ------------------------------------------------------------------- /
  !> W3SIN6 (w3src6md.F90:298) on Kokkos, same dummy list. Where W3SIN6,
  !> LFACTOR or TAU_WAVE_ATMOS stop, this writes their FORMAT to NDSE and
  !> stops with the same EXTCDE code. Where LFACTOR reaches ITERMAX
  !> (:1093-1094) the Fortran writes a warning to NDST from values the
  !> kernel does not return (TAU, TAU_TOT, ERR, TAU_NND): that point is
  !> computed again by W3SIN6 itself, which writes the warning and gives
  !> the outputs.
  SUBROUTINE W3KOKKOS_SIN6 ( A, CG, WN2, UABS, USTAR, USDIR, CD, DAIR, &
       TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN, S, D )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    USE W3SRC6MD, ONLY: W3SIN6
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), WN2(NSPEC)
    REAL, INTENT(IN)        :: UABS, USTAR, USDIR, CD, DAIR
    REAL, INTENT(OUT)       :: TAUWX, TAUWY, TAUNWX, TAUNWY, CHARN
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    REAL(C_FLOAT)           :: UA1(1), US1(1), UD1(1), CD1(1), DA1(1),  &
                               TX1(1), TY1(1), TNX1(1), TNY1(1), CH1(1)
    INTEGER(C_INT)          :: RC, INFO(2), IERR
    UA1(1) = UABS
    US1(1) = USTAR
    UD1(1) = USDIR
    CD1(1) = CD
    DA1(1) = DAIR
    !$OMP CRITICAL (WW_KOKKOS_ST6)
    RC = WW_SIN6 ( IGRID, 1, A, CG, WN2, UA1, US1, UD1, CD1, DA1, TX1,  &
         TY1, TNX1, TNY1, CH1, S, D, INFO )
    IERR = WW_ST6_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST6)
    CALL ST6_CHECK ( 'W3SIN6', IERR )
    SELECT CASE ( RC )
    CASE ( ST6_C_IKN )
      WRITE (NDSE,6061) INFO(1), INFO(2)
      CALL EXTCDE ( 2 )
    CASE ( ST6_C_LF_IK10HZ )
      WRITE (NDSE,6062) INFO(1), INFO(2)
      CALL EXTCDE ( 2 )
    CASE ( ST6_C_LF_MK10HZ )
      WRITE (NDSE,6063) INFO(1), INFO(2)
      CALL EXTCDE ( 3 )
    CASE ( ST6_C_TWA_MK10HZ )
      WRITE (NDSE,6064) INFO(1), INFO(2)
      CALL EXTCDE ( 3 )
    CASE ( ST6_C_LF_ITERMAX )
      CALL W3SIN6 ( A, CG, WN2, UABS, USTAR, USDIR, CD, DAIR, TAUWX,    &
           TAUWY, TAUNWX, TAUNWY, CHARN, S, D )
      RETURN
    END SELECT
    TAUWX = TX1(1)
    TAUWY = TY1(1)
    TAUNWX = TNX1(1)
    TAUNWY = TNY1(1)
    CHARN = CH1(1)
    ! The FORMATs of w3src6md.F90:534-535, :1127-1130, :1322-1323, verbatim.
6061 FORMAT (' *** WAVEWATCH III ERROR IN W3SIN6 : '/ &
       '     INDEX ARRAY DIM MISMATCH SIZE(IKN)=',I4,' .NE. NK=',I4)
6062 FORMAT (' *** WAVEWATCH III ERROR IN W3SIN6 LFACTOR: '/ &
       '     INDEX ARRAY DIM MISMATCH SIZE(IKN)=',I4,' .NE. NK=',I4)
6063 FORMAT (' *** WAVEWATCH III ERROR IN W3SIN6 LFACTOR: '/ &
       '     ARRAY SIZE MISMATCH MK=',I4,' .NE. NK=',I4)
6064 FORMAT (' *** WAVEWATCH III ERROR IN W3SIN6 TAU_WAV_ATMOS : '/ &
       '     ARRAY SIZE MISMATCH MK=',I4,' .NE. NK=',I4)
  END SUBROUTINE W3KOKKOS_SIN6
  !/ ------------------------------------------------------------------- /
  !> W3SDS6 (w3src6md.F90:566) on Kokkos, same dummy list. Where W3SDS6
  !> stops (SIZE(IKN) /= NK, :699-702) this writes its FORMAT 601 and
  !> stops with the same EXTCDE(2). FORMAT of w3src6md.F90:765-766, verbatim.
  SUBROUTINE W3KOKKOS_SDS6 ( A, CG, WN, S, D )
    USE W3GDATMD, ONLY: NK, NSPEC, IGRID
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    REAL, INTENT(IN)        :: A(NSPEC), CG(NK), WN(NK)
    REAL, INTENT(OUT)       :: S(NSPEC), D(NSPEC)
    INTEGER(C_INT)          :: RC, INFO(2), IERR
    !$OMP CRITICAL (WW_KOKKOS_ST6)
    RC = WW_SDS6 ( IGRID, 1, A, CG, WN, S, D, INFO )
    IERR = WW_ST6_LAST_ERROR()
    !$OMP END CRITICAL (WW_KOKKOS_ST6)
    CALL ST6_CHECK ( 'W3SDS6', IERR )
    IF ( RC == ST6_C_IKN ) THEN
      WRITE (NDSE,6065) INFO(1), INFO(2)
      CALL EXTCDE ( 2 )
    END IF
6065 FORMAT (' *** WAVEWATCH III ERROR IN W3SDS6 : '/ &
       '     INDEX ARRAY DIM MISMATCH SIZE(IKN)=',I4,' .NE. NK=',I4)
  END SUBROUTINE W3KOKKOS_SDS6
#endif
  !/ ------------------------------------------------------------------- /
  !/
END MODULE W3KOKKOSMD
