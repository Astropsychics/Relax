
MODULE flux_mapping
	! 1 LY = 6.3e4 AU
	! 1 PC = 3.26 PC
!	REAL(KIND=8),PARAMETER :: R_map = 1.0D0*6.3D4 ! [AU]
	REAL(KIND=8),PARAMETER :: R_map = 1.0D5 ! [AU]

CONTAINS

SUBROUTINE flux_map( x_now, u_now, E )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Determine if particle crosses mapping surface
! and has sunward velocity component
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 	IMPLICIT NONE

  !! Inputs
  REAL(KIND=8)  :: x_now(3), u_now(3)
	REAL(KIND=8)	:: E

	!! Internal
	REAL(KIND=8)	:: R_now

	R_now = SQRT( x_now(1)**2 + x_now(2)**2 + x_now(3)**2 )

!	WRITE(*,*) 'rnow: R: ', R_now, R_map
 
	IF (R_now .LE. R_map) THEN
		WRITE(420,*) x_now(1), x_now(2), x_now(3), E*u_now(1), E*u_now(2), E*u_now(3)
	END IF

END SUBROUTINE flux_map


END MODULE flux_mapping

