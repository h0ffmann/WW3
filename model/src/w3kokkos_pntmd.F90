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
!> VA stay stale, and the resident state is left as it was. The same flag
!> covers the boundary output, W3IOBC ('WRITE') reading VA(:,ISBPO(ISOUT))
!> (w3iobcmd.F90:660-668), and the track output, W3IOTR reading
!> A(:,:,MAPFS(IY,IX)) where MASK1 and MAPSTA > 0 (w3iotrmd.F90:728-819).
!> Opt-in: WW_KOKKOS_IOPE=1; WW_KOKKOS_ALL does not turn it on.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_PNTMD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  USE W3KOKKOSMD, ONLY: KOKKOS_PART, W3KOKKOS_REGISTER
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: KOKKOS_IOPE, W3KOKKOS_IOPE_READ, W3KOKKOS_PNT_SETUP
  PUBLIC :: W3KOKKOS_IOBC_READ, W3KOKKOS_IOTR_READ
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
    USE W3GDATMD, ONLY: MAPFS
    USE W3ODATMD, ONLY: NOPTS, IPTINT
    IMPLICIT NONE
    INTEGER(C_INT)          :: COLS(4*MAX(NOPTS,1))
    INTEGER                 :: I, J, N
    ! IS(J) of W3IOPE (w3iopomd.F90:1051-1052, :1069); JSEA = ISEA under
    ! W3_SHRD, and a land corner gives column 0, which W3IOPE reads too.
    N = 0
    DO I=1, NOPTS
      DO J=1, 4
        N = N + 1
        COLS(N) = INT ( MAPFS(IPTINT(2,J,I),IPTINT(1,J,I)), C_INT )
      END DO
    END DO
    CALL READ_COLS ( N, COLS )
  END SUBROUTINE W3KOKKOS_IOPE_READ
  !/ ------------------------------------------------------------------- /
  !> Call in place of W3KOKKOS_VA_HOST_READ right before CALL W3IOBC
  !> ( 'WRITE', ... ): the columns ISBPO(1:NBO2(NFBPO)) of VA, the only ones
  !> its W3_SHRD gather reads (w3iobcmd.F90:660-668, ISEA = ISBPO(ISOUT) for
  !> every ISOUT of every file). ISBPO is set by ww3_grid (w3gridmd.F90:6023,
  !> a sea point, never 0) and read with the model definition.
  SUBROUTINE W3KOKKOS_IOBC_READ ( )
    USE W3ODATMD, ONLY: NFBPO, NBO2, ISBPO
    IMPLICIT NONE
    INTEGER(C_INT)          :: COLS(MAX(NBO2(NFBPO),1))
    INTEGER                 :: I, N
    N = NBO2(NFBPO)
    DO I=1, N
      COLS(I) = INT ( ISBPO(I), C_INT )
    END DO
    CALL READ_COLS ( N, COLS )
  END SUBROUTINE W3KOKKOS_IOBC_READ
  !/ ------------------------------------------------------------------- /
  !> Called inside W3IOTR, after its section 3 has read the track file and
  !> set MASK1, before its section 4 reads A (= VA): the columns
  !> MAPFS(IY,IX) of every grid point with MASK1(IY,IX) and MAPSTA(IY,IX) > 0,
  !> the points W3IOTR writes as 'SEA' and the only ones whose spectrum it
  !> reads (w3iotrmd.F90:728-747, :774, :808-819; JSEA = ISEA under W3_SHRD).
  SUBROUTINE W3KOKKOS_IOTR_READ ( )
    USE W3GDATMD, ONLY: NX, NY, MAPFS, MAPSTA
    USE W3ODATMD, ONLY: MASK1
    IMPLICIT NONE
    INTEGER(C_INT), ALLOCATABLE :: COLS(:)
    INTEGER                 :: IX, IY, N
    N = 0
    DO IY=1, NY
      DO IX=1, NX
        IF ( MASK1(IY,IX) .AND. MAPSTA(IY,IX) .GT. 0 ) N = N + 1
      END DO
    END DO
    ALLOCATE ( COLS(MAX(N,1)) )
    N = 0
    DO IY=1, NY
      DO IX=1, NX
        IF ( MASK1(IY,IX) .AND. MAPSTA(IY,IX) .GT. 0 ) THEN
          N = N + 1
          COLS(N) = INT ( MAPFS(IY,IX), C_INT )
        END IF
      END DO
    END DO
    CALL READ_COLS ( N, COLS )
    DEALLOCATE ( COLS )
  END SUBROUTINE W3KOKKOS_IOTR_READ
  !/ ------------------------------------------------------------------- /
  !> The N columns COLS of VA current on the host (ww_va_host_read_cols).
  !> Refused inside an OpenMP parallel region, as every VA sync.
  SUBROUTINE READ_COLS ( N, COLS )
    USE W3GDATMD, ONLY: IGRID
    USE W3WDATMD, ONLY: VA
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
!$  USE OMP_LIB, ONLY: OMP_IN_PARALLEL
    IMPLICIT NONE
    INTEGER, INTENT(IN)        :: N
    INTEGER(C_INT), INTENT(IN) :: COLS(*)
    INTEGER(C_INT)             :: IERR
!$  IF ( OMP_IN_PARALLEL() ) THEN
!$    WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PNTMD : '// &
!$         'VA SYNC CALLED INSIDE AN OPENMP PARALLEL REGION'
!$    CALL EXTCDE ( 1 )
!$  END IF
    IERR = WW_VA_HOST_READ_COLS ( IGRID, VA, SIZE(VA,1), UBOUND(VA,2), &
         INT(N,C_INT), COLS )
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_PNTMD : '// &
           'KOKKOS VA PARTIAL READ FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
  END SUBROUTINE READ_COLS
END MODULE W3KOKKOS_PNTMD
