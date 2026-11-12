!/ ------------------------------------------------------------------- /
!> @file w3kokkos_exchmd.F90
!> @brief ww3-gpu wiring unit exch: the MPI messages of the aggregated
!>        exchange of the resident VA (ADR-0020, plan unit M2).
!>
!> The C++ side packs and unpacks (kokkos/src/fortran_iface/prop3_c.hpp of
!> h0ffmann/ww3-gpu, PATCH.md section 33); this module only moves the packed
!> REALs between ranks, so ww_kokkos stays free of MPI (ADR-0020 decision 3).
!> One message per pair of ranks and direction replaces the persistent
!> requests IRQSG1 and IRQSG2 of W3MPIP for a propagation step. The messages
!> run on a duplicate of MPI_COMM_WAVE, so no tag of WW3's own requests can
!> match them. A wiring unit of PATCH.md section 9 with no part of its own:
!> W3KOKKOS_PROP3_SETUP calls W3KOKKOS_EXCH_INIT when XYP3 runs under MPI.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!/ ------------------------------------------------------------------- /
MODULE W3KOKKOS_EXCHMD
#ifdef W3_MPI
  USE MPI_F08
#endif
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: W3KOKKOS_EXCH_SETUP, W3KOKKOS_EXCH_INIT
#ifdef W3_MPI
  PUBLIC :: W3KOKKOS_EXCH_PAIRS
  !> The duplicate of each grid's MPI_COMM_WAVE, by IMOD.
  TYPE(MPI_COMM), ALLOCATABLE, SAVE :: XCOMM(:)
  LOGICAL, ALLOCATABLE, SAVE        :: XCOMM_SET(:)
  INTEGER, PARAMETER                :: XTAG = 1
#endif
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> Nothing to set up for every grid: the communicator is duplicated only
  !> when a part asks for it (W3KOKKOS_EXCH_INIT).
  SUBROUTINE W3KOKKOS_EXCH_SETUP ( IMOD )
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
    IF ( IMOD < 1 ) RETURN
  END SUBROUTINE W3KOKKOS_EXCH_SETUP
  !/ ------------------------------------------------------------------- /
  !> Duplicate MPI_COMM_WAVE of grid IMOD. Collective over that
  !> communicator: W3INIT calls it, through the set-up of the units, on every
  !> rank of the grid, IAPROC > NAPROC included. A second call for the same
  !> grid keeps the first duplicate.
  SUBROUTINE W3KOKKOS_EXCH_INIT ( IMOD )
#ifdef W3_MPI
    USE W3GDATMD, ONLY: NGRIDS
    USE W3ADATMD, ONLY: MPI_COMM_WAVE
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
#endif
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IMOD
#ifdef W3_MPI
    TYPE(MPI_COMM), ALLOCATABLE :: TC(:)
    LOGICAL, ALLOCATABLE        :: TS(:)
    INTEGER                     :: N, IERR_MPI
    IF ( IMOD < 1 ) RETURN
    N = MAX ( IMOD, NGRIDS, 1 )
    IF ( .NOT. ALLOCATED(XCOMM) ) THEN
      ALLOCATE ( XCOMM(N), XCOMM_SET(N) )
      XCOMM_SET = .FALSE.
    ELSE IF ( SIZE(XCOMM) < N ) THEN
      ALLOCATE ( TC(N), TS(N) )
      TS = .FALSE.
      TC(1:SIZE(XCOMM)) = XCOMM
      TS(1:SIZE(XCOMM)) = XCOMM_SET
      CALL MOVE_ALLOC ( TC, XCOMM )
      CALL MOVE_ALLOC ( TS, XCOMM_SET )
    END IF
    IF ( XCOMM_SET(IMOD) ) RETURN
    CALL MPI_COMM_DUP ( MPI_COMM_WAVE, XCOMM(IMOD), IERR_MPI )
    IF ( IERR_MPI /= MPI_SUCCESS ) THEN
      WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_EXCHMD : '// &
           'MPI_COMM_DUP OF MPI_COMM_WAVE FAILED'
      CALL EXTCDE ( 1 )
    END IF
    XCOMM_SET(IMOD) = .TRUE.
#else
    IF ( IMOD < 1 ) RETURN
#endif
  END SUBROUTINE W3KOKKOS_EXCH_INIT
#ifdef W3_MPI
  !/ ------------------------------------------------------------------- /
  !> One message each way with every other rank of grid IGRID's
  !> MPI_COMM_WAVE: rank IP (1 .. NAPROC) gets SCNT(IP) REALs from
  !> SBUF(SOFF(IP)+1), and its RCNT(IP) REALs land at RBUF(ROFF(IP)+1). A zero
  !> count sends or posts nothing; both ranks of a pair compute the same
  !> counts, so a message is posted on both sides or on neither. The entries
  !> of this rank, IAPROC, are not read. Every receive is posted before the
  !> sends, then all are waited for: the buffers are free when this returns.
  SUBROUTINE W3KOKKOS_EXCH_PAIRS ( IGRID, SBUF, SOFF, SCNT, RBUF, ROFF, &
       RCNT )
    USE W3ODATMD, ONLY: NDSE, NAPROC, IAPROC
    USE W3SERVMD, ONLY: EXTCDE
    IMPLICIT NONE
    INTEGER, INTENT(IN)     :: IGRID
    REAL, INTENT(IN)        :: SBUF(*)
    REAL, INTENT(INOUT)     :: RBUF(*)
    INTEGER, INTENT(IN)     :: SOFF(*), SCNT(*), ROFF(*), RCNT(*)
    TYPE(MPI_REQUEST), ALLOCATABLE :: RQ(:)
    INTEGER                 :: IP, NRQ, IERR_MPI
    IF ( .NOT. ALLOCATED(XCOMM_SET) ) CALL NOT_SET
    IF ( IGRID < 1 .OR. IGRID > SIZE(XCOMM_SET) ) CALL NOT_SET
    IF ( .NOT. XCOMM_SET(IGRID) ) CALL NOT_SET
    ALLOCATE ( RQ(2*MAX(1,NAPROC)) )
    NRQ = 0
    DO IP = 1, NAPROC
      IF ( IP .EQ. IAPROC .OR. RCNT(IP) .LE. 0 ) CYCLE
      NRQ = NRQ + 1
      CALL MPI_IRECV ( RBUF(ROFF(IP)+1), RCNT(IP), MPI_REAL, IP-1, XTAG, &
           XCOMM(IGRID), RQ(NRQ), IERR_MPI )
    END DO
    DO IP = 1, NAPROC
      IF ( IP .EQ. IAPROC .OR. SCNT(IP) .LE. 0 ) CYCLE
      NRQ = NRQ + 1
      CALL MPI_ISEND ( SBUF(SOFF(IP)+1), SCNT(IP), MPI_REAL, IP-1, XTAG, &
           XCOMM(IGRID), RQ(NRQ), IERR_MPI )
    END DO
    IF ( NRQ .GT. 0 ) CALL MPI_WAITALL ( NRQ, RQ(1:NRQ),               &
         MPI_STATUSES_IGNORE, IERR_MPI )
    DEALLOCATE ( RQ )
  CONTAINS
    SUBROUTINE NOT_SET
      WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_EXCHMD : '// &
           'W3KOKKOS_EXCH_PAIRS BEFORE W3KOKKOS_EXCH_INIT OF THE GRID'
      CALL EXTCDE ( 1 )
    END SUBROUTINE NOT_SET
  END SUBROUTINE W3KOKKOS_EXCH_PAIRS
#endif
  !/ ------------------------------------------------------------------- /
END MODULE W3KOKKOS_EXCHMD
