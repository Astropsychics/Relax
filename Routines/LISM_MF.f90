
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!! Gives magnetic field vector at current 
!! location 
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE LISM_MF( R, B )
	USE ENA, ONLY : B0
	USE physics_constants, ONLY : PI

	IMPLICIT NONE

	!! Input
	REAL(KIND=8)	:: R(3)

	!! Output
	REAL(KIND=8)	:: B(3)

	!! Internal
	REAL(KIND=8)	:: Bt, Bp

  Bt   = 38.0D0*PI/180.0D0
  Bp   = 23.0D0*PI/180.0D0
  B(1) = B0*COS(Bp)*COS(Bt)
  B(2) = B0*COS(Bp)*SIN(Bt)
  B(3) = B0*SIN(Bp)		


END SUBROUTINE LISM_MF


