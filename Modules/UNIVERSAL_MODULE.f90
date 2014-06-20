
MODULE universal

	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE	:: table_X_H, table_X_He, table_X_O, table_X_Ar	
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE	:: table_X_H2, table_X_N2, table_X_CO, table_X_CO2	
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE	:: table_X_N, table_X_S, table_X_O2	
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE	:: table_X_SO, table_X_SO2
	REAL(KIND=8),DIMENSION(:), ALLOCATABLE	:: F_ENERGY,  F_ANGLE
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_H,  F_PD_X_He
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_O,  F_PD_X_Ar
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_H2, F_PD_X_N2
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_CO, F_PD_X_CO2
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_N,  F_PD_X_S
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_SO, F_PD_X_O2
	REAL(KIND=8),DIMENSION(:,:),ALLOCATABLE :: F_PD_X_SO2

	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_ENERGY
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_H,  TCS_X_He
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_O,  TCS_X_Ar
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_H2, TCS_X_N2
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_CO, TCS_X_CO2
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_S,  TCS_X_N
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_O2, TCS_X_SO
	REAL(KIND=8),DIMENSION(:),ALLOCATABLE		:: TCS_X_SO2

	INTEGER																	:: NL_TCS
	INTEGER																	:: NL_X_H, NL_X_He, NL_X_O, NL_X_Ar
	INTEGER																	:: NL_X_H2, NL_X_N2, NL_X_CO, NL_X_CO2
	INTEGER																	:: NL_X_N, NL_X_S, NL_X_O2, NL_X_SO, NL_X_SO2
	INTEGER																	:: NUM_ENG, NUM_ANG

  REAL,PARAMETER                          :: MB = 8.0D-6
  REAL,PARAMETER                          :: GB = 8.0D-9

	REAL(KIND=8)		:: tau_0, tau_b

CONTAINS


!################################################
!################################################

SUBROUTINE read_uni_tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in Angular Probability data table from file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE planet, ONLY : PROJ, ATMOSPHERE, VENUS_ATMOSPHERE
  USE mpi_info

  IMPLICIT NONE

	INCLUDE 'mpif.h'

  REAL(KIND=8)                            :: E0
  REAL(KIND=8)                            :: Enow

  INTEGER                                 :: N
  INTEGER                                 :: i
  INTEGER                                 :: j
  INTEGER                                 :: k, s1, s2, s3, s4

	IF (myid .EQ. 0) THEN
		IF (ATMOSPHERE .NE. 99) THEN
			IF (PROJ .EQ. 'H') THEN
				!! open tables files to read from
 				OPEN(UNIT=10, FILE="../Tables/H_O_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=11, FILE="../Tables/H_Ar_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=12, FILE="../Tables/H_H2_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=13, FILE="../Tables/H_N2_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=14, FILE="../Tables/H_CO_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=15, FILE="../Tables/H_CO2_ANG.dat", STATUS="old", ACTION="read")
 		 		READ(10,*) NL_X_O
 		 		READ(11,*) NL_X_Ar
 		 		READ(12,*) NL_X_H2
 		 		READ(13,*) NL_X_N2
 		 		READ(14,*) NL_X_CO
 		 		READ(15,*) NL_X_CO2
 		 		ALLOCATE(table_X_O(NL_X_O ,3))
  			ALLOCATE(table_X_Ar(NL_X_Ar ,3))
  			ALLOCATE(table_X_H2(NL_X_H2,3))
  			ALLOCATE(table_X_N2(NL_X_N2,3))
  			ALLOCATE(table_X_CO(NL_X_CO ,3))
  			ALLOCATE(table_X_CO2(NL_X_CO2 ,3))
  			DO i=1,NL_X_O
					READ(10,*) table_X_O(i,1), table_X_O(i,2), table_X_O(i,3)
					READ(11,*) table_X_Ar(i,1), table_X_Ar(i,2), table_X_Ar(i,3)
					READ(12,*) table_X_H2(i,1), table_X_H2(i,2), table_X_H2(i,3)
					READ(13,*) table_X_N2(i,1), table_X_N2(i,2), table_X_N2(i,3)
					READ(14,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
					READ(15,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
  			END DO
			ELSE IF (PROJ .EQ. 'He') THEN
 				OPEN(UNIT=10, FILE="../Tables/He_Ar_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=11, FILE="../Tables/He_H2_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=12, FILE="../Tables/He_N2_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=13, FILE="../Tables/He_CO_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=14, FILE="../Tables/He_CO2_ANG.dat", STATUS="old", ACTION="read")
  			READ(10,*) NL_X_Ar
  			READ(11,*) NL_X_H2
  			READ(12,*) NL_X_N2
  			READ(13,*) NL_X_CO
  			READ(14,*) NL_X_CO2
  			ALLOCATE(table_X_Ar(NL_X_Ar ,3))
  			ALLOCATE(table_X_H2(NL_X_H2,3))
  			ALLOCATE(table_X_N2(NL_X_N2,3))
  			ALLOCATE(table_X_CO(NL_X_CO ,3))
  			ALLOCATE(table_X_CO2(NL_X_CO2 ,3))
  			DO i=1,NL_X_Ar
					READ(10,*) table_X_Ar(i,1), table_X_Ar(i,2), table_X_Ar(i,3)
					READ(11,*) table_X_H2(i,1), table_X_H2(i,2), table_X_H2(i,3)
					READ(12,*) table_X_N2(i,1), table_X_N2(i,2), table_X_N2(i,3)
					READ(13,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
					READ(14,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
  			END DO
			ELSE IF (PROJ .EQ. 'O') THEN
 				OPEN(UNIT=10, FILE="../Tables/O_Ar_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=11, FILE="../Tables/O_H2_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=12, FILE="../Tables/O_N2_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=13, FILE="../Tables/O_CO_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=14, FILE="../Tables/O_CO2_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=15, FILE="../Tables/O_H_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=16, FILE="../Tables/O_He_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=17, FILE="../Tables/O_O_ANG.dat", STATUS="old", ACTION="read")
  			READ(10,*) NL_X_Ar
  			READ(11,*) NL_X_H2
  			READ(12,*) NL_X_N2
  			READ(13,*) NL_X_CO
  			READ(14,*) NL_X_CO2
  			READ(15,*) NL_X_H
  			READ(16,*) NL_X_He
  			READ(17,*) NL_X_O
  			ALLOCATE(table_X_H(NL_X_H ,3))
  			ALLOCATE(table_X_He(NL_X_He ,3))
  			ALLOCATE(table_X_O(NL_X_O ,3))
  			ALLOCATE(table_X_Ar(NL_X_Ar ,3))
  			ALLOCATE(table_X_H2(NL_X_H2,3))
  			ALLOCATE(table_X_N2(NL_X_N2,3))
  			ALLOCATE(table_X_CO(NL_X_CO ,3))
  			ALLOCATE(table_X_CO2(NL_X_CO2 ,3))
  			DO i=1,NL_X_Ar
					READ(10,*) table_X_Ar(i,1), table_X_Ar(i,2), table_X_Ar(i,3)
					READ(11,*) table_X_H2(i,1), table_X_H2(i,2), table_X_H2(i,3)
					READ(12,*) table_X_N2(i,1), table_X_N2(i,2), table_X_N2(i,3)
					READ(13,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
					READ(14,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
					READ(15,*) table_X_H(i,1), table_X_H(i,2), table_X_H(i,3)
					READ(16,*) table_X_He(i,1), table_X_He(i,2), table_X_He(i,3)
					READ(17,*) table_X_O(i,1), table_X_O(i,2), table_X_O(i,3)
  			END DO
			ELSE IF (PROJ .EQ. 'H2O') THEN
 				OPEN(UNIT=10, FILE="../Tables/H2O_H_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=11, FILE="../Tables/H2O_He_ANG.dat", STATUS="old", ACTION="read")
 			 	OPEN(UNIT=12, FILE="../Tables/H2O_O_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=13, FILE="../Tables/H2O_Ar_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=14, FILE="../Tables/H2O_H2_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=15, FILE="../Tables/H2O_N2_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=16, FILE="../Tables/H2O_CO_ANG.dat", STATUS="old", ACTION="read")
 				OPEN(UNIT=17, FILE="../Tables/H2O_CO2_ANG.dat", STATUS="old", ACTION="read")
  			READ(10,*) NL_X_H
  			READ(11,*) NL_X_He
  			READ(12,*) NL_X_O
  			READ(13,*) NL_X_Ar
  			READ(14,*) NL_X_H2
  			READ(15,*) NL_X_N2
  			READ(16,*) NL_X_CO
  			READ(17,*) NL_X_CO2
  			ALLOCATE(table_X_H(NL_X_H ,3))
  			ALLOCATE(table_X_He(NL_X_He ,3))
  			ALLOCATE(table_X_O(NL_X_O ,3))
  			ALLOCATE(table_X_Ar(NL_X_Ar ,3))
  			ALLOCATE(table_X_H2(NL_X_H2,3))
  			ALLOCATE(table_X_N2(NL_X_N2,3))
  			ALLOCATE(table_X_CO(NL_X_CO ,3))
  			ALLOCATE(table_X_CO2(NL_X_CO2 ,3))
  			DO i=1,NL_X_Ar
					READ(10,*) table_X_H(i,1), table_X_H(i,2), table_X_H(i,3)
					READ(11,*) table_X_He(i,1), table_X_He(i,2), table_X_He(i,3)
					READ(12,*) table_X_O(i,1), table_X_O(i,2), table_X_O(i,3)
					READ(13,*) table_X_Ar(i,1), table_X_Ar(i,2), table_X_Ar(i,3)
					READ(14,*) table_X_H2(i,1), table_X_H2(i,2), table_X_H2(i,3)
					READ(15,*) table_X_N2(i,1), table_X_N2(i,2), table_X_N2(i,3)
					READ(16,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
					READ(17,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
  			END DO
			ELSE 
				WRITE(*,*) 'Current Projectile: ', PROJ, ' not currently supported for Mars atmosphere'
			END IF
		END IF ! if Mars atmosphere
		IF (VENUS_ATMOSPHERE .NE. 99) THEN
			IF (PROJ .EQ. 'H') THEN
        OPEN(UNIT=10, FILE="../Tables/H_N_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=11, FILE="../Tables/H_O_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=12, FILE="../Tables/H_S_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=13, FILE="../Tables/H_CO_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=14, FILE="../Tables/H_CO2_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=15, FILE="../Tables/H_O2_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=16, FILE="../Tables/H_SO_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=17, FILE="../Tables/H_SO2_ANG.dat", STATUS="old",ACTION="read")
        READ(10,*) NL_X_N
        READ(11,*) NL_X_O
        READ(12,*) NL_X_S
        READ(13,*) NL_X_CO
        READ(14,*) NL_X_CO2
        READ(15,*) NL_X_O2
        READ(16,*) NL_X_SO
        READ(17,*) NL_X_SO2
        ALLOCATE(table_X_N(NL_X_N ,3))
        ALLOCATE(table_X_O(NL_X_O ,3))
        ALLOCATE(table_X_S(NL_X_S ,3))
        ALLOCATE(table_X_CO(NL_X_CO ,3))
        ALLOCATE(table_X_CO2(NL_X_CO2,3))
        ALLOCATE(table_X_O2(NL_X_O2,3))
        ALLOCATE(table_X_SO(NL_X_SO ,3))
        ALLOCATE(table_X_SO2(NL_X_SO2 ,3))
        DO i=1,NL_X_N
          READ(10,*) table_X_N(i,1), table_X_N(i,2), table_X_N(i,3)
          READ(11,*) table_X_O(i,1), table_X_O(i,2), table_X_O(i,3)
          READ(12,*) table_X_S(i,1), table_X_S(i,2), table_X_S(i,3)
          READ(13,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
          READ(14,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
          READ(15,*) table_X_O2(i,1), table_X_O2(i,2), table_X_O2(i,3)
          READ(16,*) table_X_SO(i,1), table_X_SO(i,2), table_X_SO(i,3)
          READ(17,*) table_X_SO2(i,1), table_X_SO2(i,2), table_X_SO2(i,3)
        END DO
			ELSE IF (PROJ .EQ. 'He') THEN
        OPEN(UNIT=10, FILE="../Tables/He_N_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=12, FILE="../Tables/He_S_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=13, FILE="../Tables/He_CO_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=14, FILE="../Tables/He_CO2_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=15, FILE="../Tables/He_O2_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=16, FILE="../Tables/He_SO_ANG.dat", STATUS="old",ACTION="read")
        OPEN(UNIT=17, FILE="../Tables/He_SO2_ANG.dat", STATUS="old",ACTION="read")
        READ(10,*) NL_X_N
        READ(12,*) NL_X_S
        READ(13,*) NL_X_CO
        READ(14,*) NL_X_CO2
        READ(15,*) NL_X_O2
        READ(16,*) NL_X_SO
        READ(17,*) NL_X_SO2
        ALLOCATE(table_X_N(NL_X_N ,3))
        ALLOCATE(table_X_S(NL_X_S ,3))
        ALLOCATE(table_X_CO(NL_X_CO ,3))
        ALLOCATE(table_X_CO2(NL_X_CO2,3))
        ALLOCATE(table_X_O2(NL_X_O2,3))
        ALLOCATE(table_X_SO(NL_X_SO ,3))
        ALLOCATE(table_X_SO2(NL_X_SO2 ,3))
        DO i=1,NL_X_N
          READ(10,*) table_X_N(i,1), table_X_N(i,2), table_X_N(i,3)
          READ(12,*) table_X_S(i,1), table_X_S(i,2), table_X_S(i,3)
          READ(13,*) table_X_CO(i,1), table_X_CO(i,2), table_X_CO(i,3)
          READ(14,*) table_X_CO2(i,1), table_X_CO2(i,2), table_X_CO2(i,3)
          READ(15,*) table_X_O2(i,1), table_X_O2(i,2), table_X_O2(i,3)
          READ(16,*) table_X_SO(i,1), table_X_SO(i,2), table_X_SO(i,3)
          READ(17,*) table_X_SO2(i,1), table_X_SO2(i,2), table_X_SO2(i,3)
        END DO
			END IF
		END IF ! Venus atmosphere

		IF (ATMOSPHERE .NE. 99) THEN
		 	E0    = table_X_Ar(1,1)
 		 	Enow  = E0
 		 	N     = 0
 			DO WHILE (Enow == E0)
 	  		N    = N + 1
 		   	Enow = table_X_Ar(N,1)
 		 	END DO
  		NUM_ANG = N-1
  		NUM_ENG = NL_X_Ar/NUM_ANG
			ALLOCATE(F_ENERGY(NUM_ENG), F_ANGLE(NUM_ANG))
  		DO i=1,NUM_ANG
   		 	F_ANGLE(i)  = table_X_Ar(i,2)
  		END DO
		ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
		 	E0    = table_X_N(1,1)
 		 	Enow  = E0
 		 	N     = 0
 			DO WHILE (Enow == E0)
 	  		N    = N + 1
 		   	Enow = table_X_N(N,1)
 		 	END DO
  		NUM_ANG = N-1
  		NUM_ENG = NL_X_N/NUM_ANG
			ALLOCATE(F_ENERGY(NUM_ENG), F_ANGLE(NUM_ANG))
  		DO i=1,NUM_ANG
   		 	F_ANGLE(i)  = table_X_N(i,2)
  		END DO
		END IF

		IF (ATMOSPHERE .NE. 99) THEN
			IF (PROJ .EQ. 'H') THEN
				ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
				s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 6*SIZE(F_PD_X_Ar)
  			DO i=1,NUM_ENG
    			DO j=1,NUM_ANG
    	  		k = (i-1)*NUM_ANG + j
    	  		F_PD_X_O(i,j) = table_X_O(k,3)
    	  		F_PD_X_Ar(i,j) = table_X_Ar(k,3)
    	  		F_PD_X_H2(i,j) = table_X_H2(k,3)
    	  		F_PD_X_N2(i,j) = table_X_N2(k,3)
    	  		F_PD_X_CO(i,j) = table_X_CO(k,3)
    	  		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    			END DO
    			F_ENERGY(i) = table_X_Ar(k,1)
  			END DO
  			DEALLOCATE(table_X_O, table_X_Ar, table_X_H2, table_X_N2, table_X_CO, table_X_CO2)
  			CLOSE(10)
  			CLOSE(11)
  			CLOSE(12)
  			CLOSE(13)
  			CLOSE(14)
  			CLOSE(15)
			ELSE IF (PROJ .EQ. 'He') THEN
				ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
				s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 5*SIZE(F_PD_X_Ar)
  			DO i=1,NUM_ENG
   		 		DO j=1,NUM_ANG
      			k = (i-1)*NUM_ANG + j
     		 		F_PD_X_Ar(i,j) = table_X_Ar(k,3)
     		 		F_PD_X_H2(i,j) = table_X_H2(k,3)
     		 		F_PD_X_N2(i,j) = table_X_N2(k,3)
      			F_PD_X_CO(i,j) = table_X_CO(k,3)
      			F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    			END DO
    			F_ENERGY(i) = table_X_Ar(k,1)
  			END DO
  			DEALLOCATE(table_X_Ar, table_X_H2, table_X_N2, table_X_CO, table_X_CO2)
  			CLOSE(10)
  			CLOSE(11)
  			CLOSE(12)
  			CLOSE(13)
  			CLOSE(14)
			ELSE IF (PROJ .EQ. 'O') THEN
				ALLOCATE(F_PD_X_H(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_He(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
				s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 5*SIZE(F_PD_X_Ar)
  			DO i=1,NUM_ENG
    			DO j=1,NUM_ANG
      			k = (i-1)*NUM_ANG + j
      			F_PD_X_H(i,j)   = table_X_H(k,3)
      			F_PD_X_He(i,j)  = table_X_He(k,3)
      			F_PD_X_O(i,j)   = table_X_O(k,3)
      			F_PD_X_Ar(i,j)  = table_X_Ar(k,3)
      			F_PD_X_H2(i,j)  = table_X_H2(k,3)
      			F_PD_X_N2(i,j)  = table_X_N2(k,3)
      			F_PD_X_CO(i,j)  = table_X_CO(k,3)
     		 		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    			END DO
    			F_ENERGY(i) = table_X_Ar(k,1)
  			END DO
  			DEALLOCATE(table_X_H, table_X_He, table_X_O, table_X_Ar)
				DEALLOCATE(table_X_H2, table_X_N2, table_X_CO, table_X_CO2)
  			CLOSE(10)
  			CLOSE(11)
  			CLOSE(12)
  			CLOSE(13)
  			CLOSE(14)
  			CLOSE(15)
  			CLOSE(16)
  			CLOSE(17)
			ELSE IF (PROJ .EQ. 'H2O') THEN
				ALLOCATE(F_PD_X_H(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_He(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
				s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 5*SIZE(F_PD_X_Ar)
  			DO i=1,NUM_ENG
    			DO j=1,NUM_ANG
      			k = (i-1)*NUM_ANG + j
      			F_PD_X_H(i,j) = table_X_H(k,3)
      			F_PD_X_He(i,j) = table_X_He(k,3)
      			F_PD_X_O(i,j) = table_X_O(k,3)
      			F_PD_X_Ar(i,j) = table_X_Ar(k,3)
      			F_PD_X_H2(i,j) = table_X_H2(k,3)
    	  		F_PD_X_N2(i,j) = table_X_N2(k,3)
    	  		F_PD_X_CO(i,j) = table_X_CO(k,3)
    	  		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    			END DO
    			F_ENERGY(i) = table_X_Ar(k,1)
  			END DO
  			DEALLOCATE(table_X_H, table_X_He, table_X_O, table_X_Ar) 
				DEALLOCATE(table_X_H2, table_X_N2, table_X_CO, table_X_CO2)
  			CLOSE(10)
  			CLOSE(11)
  			CLOSE(12)
  			CLOSE(13)
  			CLOSE(14)
  			CLOSE(15)
  			CLOSE(16)
 		 		CLOSE(17)
			END IF
		ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
			IF (PROJ .EQ. 'H') THEN
				ALLOCATE(F_PD_X_N(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_S(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_O2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_SO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_SO2(NUM_ENG,NUM_ANG))	
				s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 8*SIZE(F_PD_X_N)
  			DO i=1,NUM_ENG
    			DO j=1,NUM_ANG
    	  		k = (i-1)*NUM_ANG + j
    	  		F_PD_X_N(i,j)   = table_X_N(k,3)
    	  		F_PD_X_O(i,j)   = table_X_O(k,3)
    	  		F_PD_X_S(i,j)   = table_X_S(k,3)
    	  		F_PD_X_CO(i,j)  = table_X_CO(k,3)
    	  		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    	  		F_PD_X_O2(i,j)  = table_X_O2(k,3)
    	  		F_PD_X_SO(i,j)  = table_X_SO(k,3)
    	  		F_PD_X_SO2(i,j) = table_X_SO2(k,3)
    			END DO
    			F_ENERGY(i) = table_X_N(k,1)
  			END DO
  			DEALLOCATE(table_X_N, table_X_O, table_X_S, table_X_CO, table_X_CO2)
  			DEALLOCATE(table_X_O2, table_X_SO, table_X_SO2)
  			CLOSE(10)
  			CLOSE(11)
  			CLOSE(12)
  			CLOSE(13)
  			CLOSE(14)
  			CLOSE(15)
  			CLOSE(16)
  			CLOSE(17)
			ELSE IF (PROJ .EQ. 'He') THEN
				ALLOCATE(F_PD_X_N(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_S(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_O2(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_SO(NUM_ENG,NUM_ANG))	
				ALLOCATE(F_PD_X_SO2(NUM_ENG,NUM_ANG))	
				s1 = SIZE(F_ENERGY) + SIZE(F_ANGLE) + 7*SIZE(F_PD_X_N)
  			DO i=1,NUM_ENG
    			DO j=1,NUM_ANG
    	  		k = (i-1)*NUM_ANG + j
    	  		F_PD_X_N(i,j)   = table_X_N(k,3)
    	  		F_PD_X_S(i,j)   = table_X_S(k,3)
    	  		F_PD_X_CO(i,j)  = table_X_CO(k,3)
    	  		F_PD_X_CO2(i,j) = table_X_CO2(k,3)
    	  		F_PD_X_O2(i,j)  = table_X_O2(k,3)
    	  		F_PD_X_SO(i,j)  = table_X_SO(k,3)
    	  		F_PD_X_SO2(i,j) = table_X_SO2(k,3)
    			END DO
    			F_ENERGY(i) = table_X_N(k,1)
  			END DO
  			DEALLOCATE(table_X_N, table_X_S, table_X_CO, table_X_CO2)
  			DEALLOCATE(table_X_O2, table_X_SO, table_X_SO2)
  			CLOSE(10)
  			CLOSE(12)
  			CLOSE(13)
  			CLOSE(14)
  			CLOSE(15)
  			CLOSE(16)
  			CLOSE(17)
			END IF ! Proj = He
		END IF ! Venus

  	34 FORMAT(A, ES10.2)
  	35 FORMAT(A, I6, A, I5)

    WRITE(*,'(A)')      'UNIVERSAL TABLES READ'
		WRITE(*,'(a)')      '#################################################'
    WRITE(*,35)         "N_ENERGIES: ", NUM_ENG, " N_ANGLES: ", NUM_ANG
    WRITE(*,34)         "Minimum Energy [eV]: ", F_ENERGY(1)
    WRITE(*,34)         "Maximum Energy [eV]: ", F_ENERGY(NUM_ENG)
		WRITE(*,'(a)')      '#################################################'
    WRITE(*,'(A,F6.2)') 'Universal Memory Footprint [GB]: ', np*s1*MB*1.0D-3
		WRITE(*,'(a)')      '#################################################'
  END IF ! rank = 0

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(NUM_ENG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(NUM_ANG,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	IF (myid .NE. 0) THEN
		ALLOCATE(F_ENERGY(NUM_ENG), F_ANGLE(NUM_ANG))
		IF (ATMOSPHERE .NE. 99) THEN
			IF (PROJ .EQ. 'H') THEN
 	    	ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
			ELSE IF (PROJ .EQ. 'He') THEN
      	ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
			ELSE IF (PROJ .EQ. 'O') THEN
      	ALLOCATE(F_PD_X_H(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_He(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))
     		ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
			ELSE IF (PROJ .EQ. 'H2O') THEN
      	ALLOCATE(F_PD_X_H(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_He(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_Ar(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_H2(NUM_ENG,NUM_ANG))
     		ALLOCATE(F_PD_X_N2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
			END IF
		ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
			IF (PROJ .EQ. 'H') THEN
      	ALLOCATE(F_PD_X_N(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_O(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_S(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_O2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_SO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_SO2(NUM_ENG,NUM_ANG))
			ELSE IF (PROJ .EQ. 'He') THEN
      	ALLOCATE(F_PD_X_N(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_S(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_CO2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_O2(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_SO(NUM_ENG,NUM_ANG))
      	ALLOCATE(F_PD_X_SO2(NUM_ENG,NUM_ANG))
			END IF
		END IF
	END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(F_ENERGY,NUM_ENG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)	
	CALL MPI_BCAST(F_ANGLE,NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)	

	IF (ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H') THEN
			CALL MPI_BCAST(F_PD_X_O,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_Ar,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_H2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_N2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		ELSE IF (PROJ .EQ. 'He') THEN
			CALL MPI_BCAST(F_PD_X_Ar,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_H2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_N2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		ELSE IF (PROJ .EQ. 'O') THEN
			CALL MPI_BCAST(F_PD_X_H,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_He,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_O,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_Ar,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_H2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_N2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		ELSE IF (PROJ .EQ. 'H2O') THEN
			CALL MPI_BCAST(F_PD_X_H,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_He,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_O,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_Ar,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_H2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_N2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		END IF
	ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H') THEN
			CALL MPI_BCAST(F_PD_X_N,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_O,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_S,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_O2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_SO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_SO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		ELSE IF (PROJ .EQ. 'He') THEN
			CALL MPI_BCAST(F_PD_X_N,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_S,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_CO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_O2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_SO,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
			CALL MPI_BCAST(F_PD_X_SO2,NUM_ENG*NUM_ANG,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		END IF
	END IF

END SUBROUTINE read_uni_tables

!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE clean_uni_tables
	USE planet, ONLY : PROJ, ATMOSPHERE, VENUS_ATMOSPHERE

	IF (ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H') THEN
  		DEALLOCATE(F_PD_X_O)
  		DEALLOCATE(F_PD_X_Ar)
  		DEALLOCATE(F_PD_X_H2)
  		DEALLOCATE(F_PD_X_N2)
  		DEALLOCATE(F_PD_X_CO)
  		DEALLOCATE(F_PD_X_CO2)
		ELSE IF (PROJ .EQ. 'He') THEN
  		DEALLOCATE(F_PD_X_Ar)
  		DEALLOCATE(F_PD_X_H2)
  		DEALLOCATE(F_PD_X_N2)
  		DEALLOCATE(F_PD_X_CO)
  		DEALLOCATE(F_PD_X_CO2)
		ELSE IF (PROJ .EQ. 'O') THEN
  		DEALLOCATE(F_PD_X_H)
  		DEALLOCATE(F_PD_X_He)
  		DEALLOCATE(F_PD_X_O)
  		DEALLOCATE(F_PD_X_Ar)
  		DEALLOCATE(F_PD_X_H2)
  		DEALLOCATE(F_PD_X_N2)
  		DEALLOCATE(F_PD_X_CO)
  		DEALLOCATE(F_PD_X_CO2)
		ELSE IF (PROJ .EQ. 'H2O') THEN
  		DEALLOCATE(F_PD_X_H)
  		DEALLOCATE(F_PD_X_He)
  		DEALLOCATE(F_PD_X_O)
  		DEALLOCATE(F_PD_X_Ar)
  		DEALLOCATE(F_PD_X_H2)
  		DEALLOCATE(F_PD_X_N2)
  		DEALLOCATE(F_PD_X_CO)
  		DEALLOCATE(F_PD_X_CO2)
		END IF ! proj type
	ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H') THEN
  		DEALLOCATE(F_PD_X_N)
  		DEALLOCATE(F_PD_X_O2)
  		DEALLOCATE(F_PD_X_O)
  		DEALLOCATE(F_PD_X_S)
  		DEALLOCATE(F_PD_X_SO)
  		DEALLOCATE(F_PD_X_SO2)
  		DEALLOCATE(F_PD_X_CO)
  		DEALLOCATE(F_PD_X_CO2)
		ELSE IF (PROJ .EQ. 'He') THEN
  		DEALLOCATE(F_PD_X_N)
  		DEALLOCATE(F_PD_X_O2)
  		DEALLOCATE(F_PD_X_S)
  		DEALLOCATE(F_PD_X_SO)
  		DEALLOCATE(F_PD_X_SO2)
  		DEALLOCATE(F_PD_X_CO)
  		DEALLOCATE(F_PD_X_CO2)
		END IF
	END IF ! atmosphere models

	DEALLOCATE(F_ENERGY,F_ANGLE)

END SUBROUTINE clean_uni_tables

!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE read_uni_tcs_tables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read in Angular Probability data table from file
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE planet, ONLY : PROJ, ATMOSPHERE, VENUS_ATMOSPHERE
  USE mpi_info

  IMPLICIT NONE

  INCLUDE 'mpif.h'

	INTEGER						:: i
	REAL(KIND=8)			:: dummy

	IF (myid .EQ. 0) THEN
	IF (ATMOSPHERE .NE. 99) THEN ! Mars atmosphere model
		IF (PROJ .EQ. 'H ') THEN
      OPEN(UNIT=10, FILE="../Tables/H_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=11, FILE="../Tables/H_Ar_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/H_H2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/H_N2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/H_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/H_CO2_TCS.dat", STATUS="old", ACTION="read")
      READ(10,*) NL_TCS
      READ(11,*) 
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			DO i=1,NL_TCS
				READ(10,*) TCS_ENERGY(i), TCS_X_O(i)
				READ(11,*) dummy,         TCS_X_Ar(i)
				READ(12,*) dummy,         TCS_X_H2(i)
				READ(13,*) dummy,         TCS_X_N2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
			END DO ! i
			CLOSE(10)
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
		ELSE IF (PROJ .EQ. 'He') THEN
      OPEN(UNIT=11, FILE="../Tables/He_Ar_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/He_H2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/He_N2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/He_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/He_CO2_TCS.dat", STATUS="old", ACTION="read")
      READ(11,*) NL_TCS
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			DO i=1,NL_TCS
				READ(11,*) TCS_ENERGY(i), TCS_X_Ar(i)
				READ(12,*) dummy,         TCS_X_H2(i)
				READ(13,*) dummy,         TCS_X_N2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
			END DO ! i
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
		ELSE IF (PROJ .EQ. 'O') THEN
      OPEN(UNIT=8, FILE="../Tables/H_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=9, FILE="../Tables/He_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=10, FILE="../Tables/O_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=11, FILE="../Tables/O_Ar_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/O_H2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/O_N2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/O_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/O_CO2_TCS.dat", STATUS="old", ACTION="read")
      READ(8,*) 
      READ(9,*) 
      READ(10,*) 
      READ(11,*) NL_TCS
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_H(NL_TCS))
			ALLOCATE(TCS_X_He(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			DO i=1,NL_TCS
				READ(8,*) dummy,         TCS_X_H(i)
				READ(9,*) dummy,         TCS_X_He(i)
				READ(10,*) dummy,         TCS_X_O(i)
				READ(11,*) TCS_ENERGY(i), TCS_X_Ar(i)
				READ(12,*) dummy,         TCS_X_H2(i)
				READ(13,*) dummy,         TCS_X_N2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
			END DO ! i
			CLOSE(8)
			CLOSE(9)
			CLOSE(10)
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
		ELSE IF (PROJ .EQ. 'H2O') THEN
      OPEN(UNIT=8, FILE="../Tables/H2O_H_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=9, FILE="../Tables/H2O_He_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=10, FILE="../Tables/H2O_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=11, FILE="../Tables/H2O_Ar_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/H2O_H2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/H2O_N2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/H2O_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/H2O_CO2_TCS.dat", STATUS="old", ACTION="read")
      READ(8,*) 
      READ(9,*) 
      READ(10,*) 
      READ(11,*) NL_TCS
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_H(NL_TCS))
			ALLOCATE(TCS_X_He(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			DO i=1,NL_TCS
				READ(8,*) dummy,         TCS_X_H(i)
				READ(9,*) dummy,         TCS_X_He(i)
				READ(10,*) dummy,         TCS_X_O(i)
				READ(11,*) TCS_ENERGY(i), TCS_X_Ar(i)
				READ(12,*) dummy,         TCS_X_H2(i)
				READ(13,*) dummy,         TCS_X_N2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
			END DO ! i
			CLOSE(8)
			CLOSE(9)
			CLOSE(10)
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
		ELSE
			WRITE(*,*) 'Projectile ', PROJ, ' not currently supported in tables'
		END IF ! PROJ
	ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN 
		IF (PROJ .EQ. 'H ') THEN
      OPEN(UNIT=10, FILE="../Tables/H_N_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=11, FILE="../Tables/H_O_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/H_S_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/H_O2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/H_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/H_CO2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=16, FILE="../Tables/H_SO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=17, FILE="../Tables/H_SO2_TCS.dat", STATUS="old", ACTION="read")
      READ(10,*) NL_TCS
      READ(11,*) 
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
      READ(16,*) 
      READ(17,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_N(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_S(NL_TCS))
			ALLOCATE(TCS_X_O2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			ALLOCATE(TCS_X_SO(NL_TCS))
			ALLOCATE(TCS_X_SO2(NL_TCS))
			DO i=1,NL_TCS
				READ(10,*) TCS_ENERGY(i), TCS_X_N(i)
				READ(11,*) dummy,         TCS_X_O(i)
				READ(12,*) dummy,         TCS_X_S(i)
				READ(13,*) dummy,         TCS_X_O2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
				READ(16,*) dummy,         TCS_X_SO(i)
				READ(17,*) dummy,         TCS_X_SO2(i)
			END DO ! i
			CLOSE(10)
			CLOSE(11)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
			CLOSE(16)
			CLOSE(17)
		ELSE IF (PROJ .EQ. 'He') THEN
      OPEN(UNIT=10, FILE="../Tables/He_N_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=12, FILE="../Tables/He_S_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=13, FILE="../Tables/He_O2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=14, FILE="../Tables/He_CO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=15, FILE="../Tables/He_CO2_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=16, FILE="../Tables/He_SO_TCS.dat", STATUS="old", ACTION="read")
      OPEN(UNIT=17, FILE="../Tables/He_SO2_TCS.dat", STATUS="old", ACTION="read")
      READ(10,*) NL_TCS
      READ(12,*) 
      READ(13,*) 
      READ(14,*) 
      READ(15,*) 
      READ(16,*) 
      READ(17,*) 
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_N(NL_TCS))
			ALLOCATE(TCS_X_S(NL_TCS))
			ALLOCATE(TCS_X_O2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			ALLOCATE(TCS_X_SO(NL_TCS))
			ALLOCATE(TCS_X_SO2(NL_TCS))
			DO i=1,NL_TCS
				READ(10,*) TCS_ENERGY(i), TCS_X_N(i)
				READ(12,*) dummy,         TCS_X_S(i)
				READ(13,*) dummy,         TCS_X_O2(i)
				READ(14,*) dummy,         TCS_X_CO(i)
				READ(15,*) dummy,         TCS_X_CO2(i)
				READ(16,*) dummy,         TCS_X_SO(i)
				READ(17,*) dummy,         TCS_X_SO2(i)
			END DO ! i
			CLOSE(10)
			CLOSE(12)
			CLOSE(13)
			CLOSE(14)
			CLOSE(15)
			CLOSE(16)
			CLOSE(17)
		END IF ! projectile
	END IF ! atmosphere 
	END IF ! myid = 0

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( NL_TCS, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF ( myid .NE. 0) THEN
	IF (ATMOSPHERE .NE. 99) THEN 
		IF (PROJ .EQ. 'H ') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
		ELSE IF (PROJ .EQ. 'He') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
		ELSE IF (PROJ .EQ. 'O') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_H(NL_TCS))
			ALLOCATE(TCS_X_He(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
		ELSE IF (PROJ .EQ. 'H2O') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_H(NL_TCS))
			ALLOCATE(TCS_X_He(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_Ar(NL_TCS))
			ALLOCATE(TCS_X_H2(NL_TCS))
			ALLOCATE(TCS_X_N2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
		END IF
	ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H ') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_O(NL_TCS))
			ALLOCATE(TCS_X_N(NL_TCS))
			ALLOCATE(TCS_X_S(NL_TCS))
			ALLOCATE(TCS_X_O2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			ALLOCATE(TCS_X_SO(NL_TCS))
			ALLOCATE(TCS_X_SO2(NL_TCS))
		ELSE IF (PROJ .EQ. 'He') THEN
			ALLOCATE(TCS_ENERGY(NL_TCS))
			ALLOCATE(TCS_X_N(NL_TCS))
			ALLOCATE(TCS_X_S(NL_TCS))
			ALLOCATE(TCS_X_O2(NL_TCS))
			ALLOCATE(TCS_X_CO(NL_TCS))
			ALLOCATE(TCS_X_CO2(NL_TCS))
			ALLOCATE(TCS_X_SO(NL_TCS))
			ALLOCATE(TCS_X_SO2(NL_TCS))
		END IF ! proj
	END IF ! atmosphere
	END IF ! myid
		
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	CALL MPI_BCAST(TCS_ENERGY, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

	IF (ATMOSPHERE .NE. 99) THEN
	IF (PROJ .EQ. 'H ') THEN
		CALL MPI_BCAST(TCS_X_O,  NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_Ar, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_H2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_N2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	ELSE IF (PROJ .EQ. 'He') THEN	
		CALL MPI_BCAST(TCS_X_Ar, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_H2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_N2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	ELSE IF (PROJ .EQ. 'O') THEN	
		CALL MPI_BCAST(TCS_X_H, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_He, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_O, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_Ar, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_H2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_N2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	ELSE IF (PROJ .EQ. 'H2O') THEN	
		CALL MPI_BCAST(TCS_X_H, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_He, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_O, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_Ar, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_H2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_N2, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO, NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	END IF ! proj
	ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H') THEN
		CALL MPI_BCAST(TCS_X_N,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_O,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_S,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_SO,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_SO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_O2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		ELSE IF (PROJ .EQ. 'He') THEN
		CALL MPI_BCAST(TCS_X_N,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_S,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_SO,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_SO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_O2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		CALL MPI_BCAST(TCS_X_CO2,NL_TCS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
		END IF ! proj
	END IF ! atmosphere

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE read_uni_tcs_tables

!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!

SUBROUTINE clean_uni_tcs_tables
	USE planet, ONLY : PROJ, ATMOSPHERE, VENUS_ATMOSPHERE

	IMPLICIT NONE

	IF (ATMOSPHERE .NE. 99) THEN
	IF (PROJ .EQ. 'H ') THEN
		DEALLOCATE( TCS_X_O, TCS_X_Ar, TCS_X_H2, TCS_X_N2, TCS_X_CO, TCS_X_CO2 )
	ELSE IF (PROJ .EQ. 'He') THEN
		DEALLOCATE( TCS_X_Ar, TCS_X_H2, TCS_X_N2, TCS_X_CO, TCS_X_CO2 )
	ELSE IF (PROJ .EQ. 'O') THEN
		DEALLOCATE( TCS_X_H,TCS_X_He,TCS_X_O,TCS_X_Ar,TCS_X_H2,TCS_X_N2,TCS_X_CO,TCS_X_CO2 )
	ELSE IF (PROJ .EQ. 'H2O') THEN
		DEALLOCATE( TCS_X_H,TCS_X_He,TCS_X_O,TCS_X_Ar,TCS_X_H2,TCS_X_N2,TCS_X_CO,TCS_X_CO2 )
	END IF
	ELSE IF (VENUS_ATMOSPHERE .NE. 99) THEN
		IF (PROJ .EQ. 'H') THEN
			DEALLOCATE( TCS_X_N, TCS_X_O, TCS_X_S, TCS_X_O2, TCS_X_CO, TCS_X_CO2, TCS_X_SO, TCS_X_SO2 )
		ELSE IF (PROJ .EQ. 'He') THEN
			DEALLOCATE( TCS_X_N, TCS_X_S, TCS_X_O2, TCS_X_CO, TCS_X_CO2, TCS_X_SO, TCS_X_SO2 )
		END IF
	END IF 

END SUBROUTINE clean_uni_tcs_tables

END MODULE universal


