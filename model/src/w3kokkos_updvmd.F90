!/ ------------------------------------------------------------------- /
!> @file w3kokkos_updvmd.F90
!> @brief ww3-gpu wiring unit updv (WFIP-0005): W3UICE's write of VA repeated
!>        on the device copy of VA, so an ice update copies no VA.
!>
!> Binds ww_va_host_update and ww_va_device_zero_cols of
!> kokkos/src/fortran_iface/resident_c.hpp of h0ffmann/ww3-gpu; PATCH.md
!> section 17 there is the caller side. W3UICE (w3updtmd.F90:1806-2052)
!> writes VA in one statement form only, VA(:,JSEA) = 0. at :1965 (a sea
!> point that becomes ice) and :2001 (an ice point that becomes sea and is
!> not dry), and reads no VA. W3KOKKOS_UICE runs the Fortran W3UICE unchanged
!> on the host VA, finds the columns it zeroed from MAPST2 before and after,
!> and has the same columns zeroed in the device copy. No VA crosses, and the
!> resident state is left as it was. Opt-in: WW_KOKKOS_UPDV=1; WW_KOKKOS_ALL
!> does not turn it on. W3ULEV is not covered (PATCH.md section 17).
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_UPDVMD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  USE W3KOKKOSMD, ONLY: KOKKOS_PART, W3KOKKOS_REGISTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: KOKKOS_UPDV, W3KOKKOS_UICE, W3KOKKOS_UPDV_SETUP
  !
  !> WW_KOKKOS_UPDV=1 on a build where JSEA = ISEA (W3KOKKOS_UPDV_SETUP).
  LOGICAL :: KOKKOS_UPDV = .FALSE.
  !
  INTERFACE
    !> Before the host statement: no copy (a host read under
    !> WW_KOKKOS_VA_RESIDENT=1).
    INTEGER(C_INT) FUNCTION WW_VA_HOST_UPDATE ( IGRID, VA, NSPEC,      &
         NSEALM ) BIND(C, NAME='ww_va_host_update')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE         :: IGRID, NSPEC, NSEALM
      REAL(C_FLOAT), INTENT(INOUT)  :: VA(*)
    END FUNCTION WW_VA_HOST_UPDATE
    !> After it: the device columns COLS (0-based JSEA) set to 0.
    INTEGER(C_INT) FUNCTION WW_VA_DEVICE_ZERO_COLS ( IGRID, VA, NSPEC, &
         NSEALM, NCOL, COLS ) BIND(C, NAME='ww_va_device_zero_cols')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE         :: IGRID, NSPEC, NSEALM, NCOL
      REAL(C_FLOAT), INTENT(INOUT)  :: VA(*)
      INTEGER(C_INT), INTENT(IN)    :: COLS(*)
    END FUNCTION WW_VA_DEVICE_ZERO_COLS
  END INTERFACE
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> The flag. W3_SHRD (NAPROC = 1, so JSEA = ISEA, w3parall.F90:1236) and
  !> neither W3_MPI nor W3_PDLIB, whose JSEA this wrapper does not compute.
  SUBROUTINE W3KOKKOS_UPDV_SETUP ( IMOD )
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    LOGICAL                 :: OK
    IF ( IMOD < 1 ) RETURN
    OK = .FALSE.
#ifdef W3_SHRD
    OK = .TRUE.
#endif
#if defined(W3_MPI) || defined(W3_PDLIB)
    OK = .FALSE.
#endif
    KOKKOS_UPDV = KOKKOS_PART ( 'UPDV', OK )
    CALL W3KOKKOS_REGISTER ( 'UPDV', KOKKOS_UPDV )
  END SUBROUTINE W3KOKKOS_UPDV_SETUP
  !/ ------------------------------------------------------------------- /
  !> CALL W3UICE ( VA ) of W3WAVE, with its VA write repeated on the device.
  !> Which columns W3UICE zeroes follows from its ice bit of MAPST2 (W3_IC0
  !> only; without it W3UICE writes no VA): MAPICE = MOD(MAPST2,2) and MAPST2
  !> without it on entry (:1938-1939), MAPST2 + MAPICE on exit (:2035). The
  !> bit goes 0 -> 1 exactly where :1956 holds, and :1965 zeroes the column;
  !> it goes 1 -> 0 exactly where :1987 holds, and :2001 zeroes the column
  !> when MAPST2 without the bit is 0 (:1992), that is when MAPST2 on exit is
  !> 0. The two cases exclude each other (ICEI >= FICEN against < FICEN), and
  !> MAPSF visits each sea cell once. So the columns are the sea points whose
  !> bit changed and whose MAPST2 on exit is odd or 0.
  SUBROUTINE W3KOKKOS_UICE ( VA )
    USE W3GDATMD, ONLY: NX, NY, NSEA, NSPEC, MAPSF, MAPST2, IGRID
    USE W3ADATMD, ONLY: NSEALM
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
    USE W3UPDTMD, ONLY: W3UICE
!$  USE OMP_LIB, ONLY: OMP_IN_PARALLEL
    IMPLICIT NONE
    REAL, INTENT(INOUT)     :: VA(NSPEC,0:NSEALM)
    INTEGER, ALLOCATABLE    :: MAPOLD(:,:)
    INTEGER(C_INT), ALLOCATABLE :: COLS(:)
    INTEGER(C_INT)          :: IERR
    INTEGER                 :: ISEA, IX, IY, N
!$  IF ( OMP_IN_PARALLEL() ) THEN
!$    WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_UPDVMD : '// &
!$         'VA SYNC CALLED INSIDE AN OPENMP PARALLEL REGION'
!$    CALL EXTCDE ( 1 )
!$  END IF
    ALLOCATE ( MAPOLD(NY,NX), COLS(MAX(NSEA,1)) )
    MAPOLD = MAPST2
    IERR = WW_VA_HOST_UPDATE ( IGRID, VA, INT(NSPEC,C_INT),             &
         INT(NSEALM,C_INT) )
    IF ( IERR /= 0_C_INT ) CALL FAILED ( IERR )
    CALL W3UICE ( VA )
    N = 0
#ifdef W3_IC0
    DO ISEA=1, NSEA
      IX     = MAPSF(ISEA,1)
      IY     = MAPSF(ISEA,2)
      IF ( MOD(MAPOLD(IY,IX),2) .NE. MOD(MAPST2(IY,IX),2) .AND.          &
           ( MOD(MAPST2(IY,IX),2) .EQ. 1 .OR. MAPST2(IY,IX) .EQ. 0 ) ) THEN
        N       = N + 1
        COLS(N) = INT ( ISEA, C_INT )
      END IF
    END DO
#endif
    IERR = WW_VA_DEVICE_ZERO_COLS ( IGRID, VA, INT(NSPEC,C_INT),        &
         INT(NSEALM,C_INT), INT(N,C_INT), COLS )
    IF ( IERR /= 0_C_INT ) CALL FAILED ( IERR )
    DEALLOCATE ( MAPOLD, COLS )
  CONTAINS
    SUBROUTINE FAILED ( CODE )
      INTEGER(C_INT), INTENT(IN) :: CODE
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_UPDVMD : '// &
           'KOKKOS VA UPDATE FAILED, CODE ', CODE
      CALL EXTCDE ( 1 )
    END SUBROUTINE FAILED
  END SUBROUTINE W3KOKKOS_UICE
END MODULE W3KOKKOS_UPDVMD
