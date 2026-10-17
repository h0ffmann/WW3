!/ ------------------------------------------------------------------- /
!> @file w3kokkos_residentmd.F90
!> @brief The resident spectrum (ww3-gpu WFIP-0005, units 1 and 2): the host
!>        side of the sync protocol for VA.
!>
!> Binds kokkos/src/fortran_iface/resident_c.hpp of h0ffmann/ww3-gpu. The C++
!> side keeps a device copy of VA(NSPEC,0:NSEALM) per grid and a state that
!> says which copy holds the newer bits; W3WAVE calls W3KOKKOS_VA_HOST_READ
!> before a host statement that reads VA and W3KOKKOS_VA_HOST_WRITE before one
!> that writes it (PATCH.md section 11 lists the sites). With no Kokkos kernel
!> on the resident VA the state never leaves HOST and both calls copy
!> nothing. A wiring unit of PATCH.md section 9: no part of its own.
!>
!> SPDX-License-Identifier: LGPL-3.0-or-later
!>
MODULE W3KOKKOS_RESIDENTMD
  USE ISO_C_BINDING, ONLY: C_INT, C_FLOAT
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: W3KOKKOS_VA_HOST_READ, W3KOKKOS_VA_HOST_WRITE,             &
       W3KOKKOS_RESIDENT_SETUP
  !
  INTERFACE
    !> Before a host read of VA: a copy back when the device has the newer
    !> bits. VA is VA(1,0) of grid IGRID, NSPEC = SIZE(VA,1),
    !> NSEALM = UBOUND(VA,2). Returns 0, or a status of ww_kokkos_c.hpp.
    INTEGER(C_INT) FUNCTION WW_VA_HOST_READ ( IGRID, VA, NSPEC,        &
         NSEALM ) BIND(C, NAME='ww_va_host_read')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE         :: IGRID, NSPEC, NSEALM
      REAL(C_FLOAT), INTENT(INOUT)  :: VA(*)
    END FUNCTION WW_VA_HOST_READ
    !> Before a host write of VA: WW_VA_HOST_READ, then the host is newer.
    INTEGER(C_INT) FUNCTION WW_VA_HOST_WRITE ( IGRID, VA, NSPEC,       &
         NSEALM ) BIND(C, NAME='ww_va_host_write')
      IMPORT :: C_INT, C_FLOAT
      INTEGER(C_INT), VALUE         :: IGRID, NSPEC, NSEALM
      REAL(C_FLOAT), INTENT(INOUT)  :: VA(*)
    END FUNCTION WW_VA_HOST_WRITE
  END INTERFACE
  !
CONTAINS
  !/ ------------------------------------------------------------------- /
  !> Nothing to set up: VA does not exist yet when W3INIT calls this
  !> (W3DIMW comes later), so the C++ side binds it on the first sync call.
  SUBROUTINE W3KOKKOS_RESIDENT_SETUP ( IMOD )
    INTEGER, INTENT(IN)     :: IMOD
    IF ( IMOD < 1 ) RETURN
  END SUBROUTINE W3KOKKOS_RESIDENT_SETUP
  !/ ------------------------------------------------------------------- /
  !> Call before a statement of W3WAVE that reads VA on the host.
  SUBROUTINE W3KOKKOS_VA_HOST_READ ( )
    CALL VA_SYNC ( .FALSE. )
  END SUBROUTINE W3KOKKOS_VA_HOST_READ
  !/ ------------------------------------------------------------------- /
  !> Call before a statement of W3WAVE that writes VA on the host (or reads
  !> and writes it).
  SUBROUTINE W3KOKKOS_VA_HOST_WRITE ( )
    CALL VA_SYNC ( .TRUE. )
  END SUBROUTINE W3KOKKOS_VA_HOST_WRITE
  !/ ------------------------------------------------------------------- /
  !> The two calls above. Refused inside an OpenMP parallel region: every
  !> site is outside one, and the C++ contexts are not locked.
  SUBROUTINE VA_SYNC ( LWRITE )
    USE W3GDATMD, ONLY: IGRID
    USE W3WDATMD, ONLY: VA
    USE W3ODATMD, ONLY: NDSE
    USE W3SERVMD, ONLY: EXTCDE
!$  USE OMP_LIB, ONLY: OMP_IN_PARALLEL
    LOGICAL, INTENT(IN)     :: LWRITE
    INTEGER(C_INT)          :: IERR
!$  IF ( OMP_IN_PARALLEL() ) THEN
!$    WRITE (NDSE,'(A)') '*** WAVEWATCH III ERROR IN W3KOKKOS_RESIDENTMD : '// &
!$         'VA SYNC CALLED INSIDE AN OPENMP PARALLEL REGION'
!$    CALL EXTCDE ( 1 )
!$  END IF
    IF ( LWRITE ) THEN
      IERR = WW_VA_HOST_WRITE ( IGRID, VA, SIZE(VA,1), UBOUND(VA,2) )
    ELSE
      IERR = WW_VA_HOST_READ ( IGRID, VA, SIZE(VA,1), UBOUND(VA,2) )
    END IF
    IF ( IERR /= 0_C_INT ) THEN
      WRITE (NDSE,'(A,I0)') '*** WAVEWATCH III ERROR IN W3KOKKOS_RESIDENTMD : '// &
           'KOKKOS VA SYNC FAILED, CODE ', IERR
      CALL EXTCDE ( 1 )
    END IF
  END SUBROUTINE VA_SYNC
END MODULE W3KOKKOS_RESIDENTMD
