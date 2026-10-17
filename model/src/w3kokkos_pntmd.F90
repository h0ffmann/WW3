!/ ------------------------------------------------------------------- /
!> @file w3kokkos_pntmd.F90
!> @brief ww3-gpu wiring unit pnt (WFIP-0005): the point output reads only
!>        the spectra it needs from the device copy of VA.
!>
!> Binds ww_va_host_read_cols of kokkos/src/fortran_iface/resident_c.hpp of
!> h0ffmann/ww3-gpu; PATCH.md section 16 there is the caller side. Under
!> W3_SHRD, W3IOPE (w3iopomd.F90) reads A(:,:,IS(J)) only, IS(J) =
!> MAPFS(IPTINT(2,J,I),IPTINT(1,J,I)) for the four interpolation corners J of
!> each output point I (w3iopomd.F90:1051-1052, :1069, :1214). W3KOKKOS_IOPE_READ
!> passes that column list, so a point-output time copies those columns of VA
!> back from the device and not the whole array; the other columns of the host
!> VA stay stale, and the resident state is left as it was. Opt-in:
!> WW_KOKKOS_IOPE=1; WW_KOKKOS_ALL does not turn it on.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_PNTMD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  USE W3KOKKOSMD, ONLY: KOKKOS_PART, W3KOKKOS_REGISTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: KOKKOS_IOPE, W3KOKKOS_IOPE_READ, W3KOKKOS_PNT_SETUP
  !
  !> WW_KOKKOS_IOPE=1 on a build where W3IOPE reads VA as above.
  LOGICAL :: KOKKOS_IOPE = .FALSE.
  !
  INTERFACE
    !> resident_c.hpp: COLS are 0-based JSEA, 0 .. NSEALM.
    INTEGER(C_INT) FUNCTION WW_VA_HOST_READ_COLS ( IGRID, VA, NSPEC,   &
         NSEALM, NCOL, COLS ) BIND(C, NAME='ww_va_host_read_cols')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE         :: IGRID, NSPEC, NSEALM, NCOL
      REAL(C_FLOAT), INTENT(INOUT)  :: VA(*)
      INTEGER(C_INT), INTENT(IN)    :: COLS(*)
    END FUNCTION WW_VA_HOST_READ_COLS
  END INTERFACE
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> The flag. Only W3_SHRD's copy loop of W3IOPE (w3iopomd.F90:1210-1217)
  !> reads VA by MAPFS; under W3_MPI the spectra arrive through the
  !> requests of W3MPIO, and W3WAVE keeps its single host read there.
  SUBROUTINE W3KOKKOS_PNT_SETUP ( IMOD )
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    LOGICAL                 :: OK
    IF ( IMOD < 1 ) RETURN
    OK = .FALSE.
#ifdef W3_SHRD
    OK = .TRUE.
#endif
#ifdef W3_MPI
    OK = .FALSE.
#endif
    KOKKOS_IOPE = KOKKOS_PART ( 'IOPE', OK )
    CALL W3KOKKOS_REGISTER ( 'IOPE', KOKKOS_IOPE )
  END SUBROUTINE W3KOKKOS_PNT_SETUP
  !/ ------------------------------------------------------------------- /
  !> Call in place of W3KOKKOS_VA_HOST_READ right before CALL W3IOPE ( VA ):
  !> the columns of VA W3IOPE reads, current on the host. Refused inside an
  !> OpenMP parallel region, as every VA sync.
  SUBROUTINE W3KOKKOS_IOPE_READ ( )
    USE W3GDATMD, ONLY: IGRID, MAPFS
    USE W3WDATMD, ONLY: VA
    USE W3ODATMD, ONLY: NOPTS, IPTINT, NDSE
    USE W3SERVMD, ONLY: EXTCDE
!$  USE OMP_LIB, ONLY: OMP_IN_PARALLEL
    IMPLICIT NONE
    INTEGER(C_INT)          :: COLS(4*MAX(NOPTS,1)), IERR
    INTEGER                 :: I, J, N
!$  IF ( OMP_IN_PARALLEL() ) THEN
!$    WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PNTMD : '// &
!$         'VA SYNC CALLED INSIDE AN OPENMP PARALLEL REGION'
!$    CALL EXTCDE ( 1 )
!$  END IF
    ! IS(J) of W3IOPE (w3iopomd.F90:1051-1052, :1069); JSEA = ISEA under
    ! W3_SHRD, and a land corner gives column 0, which W3IOPE reads too.
    N = 0
    DO I=1, NOPTS
      DO J=1, 4
        N = N + 1
        COLS(N) = INT ( MAPFS(IPTINT(2,J,I),IPTINT(1,J,I)), C_INT )
      END DO
    END DO
    IERR = WW_VA_HOST_READ_COLS ( IGRID, VA, SIZE(VA,1), UBOUND(VA,2), &
         INT(N,C_INT), COLS )
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PNTMD : '// &
           'KOKKOS VA PARTIAL READ FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
  END SUBROUTINE W3KOKKOS_IOPE_READ
END MODULE W3KOKKOS_PNTMD
