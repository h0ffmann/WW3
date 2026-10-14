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
          ( NK, NTH, XFR, DTH, LAM, SNLC1, KDCON, KDMN,                  &
            SNLS1, SNLS2, SNLS3, FACHFE, SIG )                           &
          BIND(C, NAME='ww_snl1_init')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: NK, NTH
       REAL(C_FLOAT),  VALUE           :: XFR, DTH, LAM, SNLC1, KDCON,   &
                                          KDMN, SNLS1, SNLS2, SNLS3,     &
                                          FACHFE
       REAL(C_FLOAT),  INTENT(IN)      :: SIG(*)
     END FUNCTION WW_SNL1_INIT
     !
     !> W3SNL1 for NPTS sea points in one launch. S and D are overwritten, not
     !> accumulated. Status in WW_SNL1_LAST_ERROR().
     SUBROUTINE WW_SNL1 ( NPTS, A, CG, KDMEAN, S, D )                    &
          BIND(C, NAME='ww_snl1')
       IMPORT :: C_INT, C_FLOAT
       INTEGER(C_INT), VALUE           :: NPTS
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
    !
  END SUBROUTINE W3KOKKOS_SETUP
  !/ ------------------------------------------------------------------- /
  !/
END MODULE W3KOKKOSMD
