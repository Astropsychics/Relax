
PROGRAM distance_finder


	IMPLICIT NONE

	INTEGER :: N, i, j, cbin, rbin, bc

	REAL(KIND=8) :: rL, rU, cL, cU, x, y, z, dr, r1, r2, dc, c1, c2, rr, cc
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)	:: r, c

	rbin = 100 
	cbin = 100
	N    = 1000000
	ALLOCATE(r(N),c(N))

	OPEN(555, FILE="./Distance_Histogram.dat", ACCESS="APPEND")
	OPEN(666, FILE="./Collision_Histogram.dat", ACCESS="APPEND")
!	OPEN(UNIT=55,FILE="./Dist_LISM_ENA.dat",ACCESS="APPEND")
	OPEN(UNIT=66,FILE="../LISM_ENA_Relax.dat",STATUS="OLD",ACTION="READ")

	rL = 10000.0D0
	rU = 0.0D0
	cL = 1000.0D0
	cU = 0.0D0

	WRITE(*,*) "Starting calculation of distances"

	DO i=1,N

		READ(66,*) x, y, z, c(i)

		r(i) = SQRT( x*x + y*y + z*z )

		IF (c(i) .LT. cL) THEN
			cL = c(i)
		END IF
		IF (c(i) .GT. cU) THEN
			cU = c(i)
		END IF

		IF (r(i) .LT. rL) THEN
			rL = r(i)
		END IF
		IF (r(i) .GT. rU) THEN
			rU = r(i)
		END IF
	
	END DO

	WRITE(*,*)
	WRITE(*,*) "r bins, lower upper: ", rL, rU 
	WRITE(*,*) "c bins, lower upper: ", cL, cU 
	
	CLOSE(55)
	CLOSE(66)

	dr = (rU - rL)/REAL(rbin)
	dc = (cU - cL)/REAL(cbin)

	WRITE(*,*) "Starting binning of distances"

	DO i=1,rbin
		r1 = rL + (i-1)*dr 
		r2 = rL + i*dr
		rr = r1 + REAL(i)*dr/2.0D0		
		bc = 0
		DO j=1,N
			IF ( r(j) .GT. r1 .AND. r(j) .LT. r2 ) THEN
				bc = bc + 1
!				WRITE(*,*) "i, j, bc ", i, j, bc
			END IF
		END DO	
		WRITE(*,*) "Bin: ", i, "Left: ", r1, "Right: ", r2, "Count: ", bc
		WRITE(555,*) rr*1.5e-5, bc
	END DO

	CLOSE(555)

	WRITE(*,*) "Starting binning of collision numbers"
				
	DO i=1,cbin
		c1 = cL + (i-1)*dc 
		c2 = cL + i*dc
		cc = c1 + REAL(i)*dc/2.0D0		
		bc = 0
		DO j=1,N
			IF ( c(j) .GT. c1 .AND. c(j) .LT. c2 ) THEN
				bc = bc + 1
			END IF
		END DO
		WRITE(666,*) cc, bc
		WRITE(*,*) "Bin: ", i, "Left: ", c1, "Right: ", c2, "Count: ", bc
	END DO
			
	CLOSE(666)
	

END PROGRAM

