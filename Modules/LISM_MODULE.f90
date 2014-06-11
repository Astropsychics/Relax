!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Module to be used for LISM calculations 
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE lism

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: H_ENA_PROB
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: H_ENA_DIST
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: He_ENA_PROB
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: He_ENA_DIST

	INTEGER																	:: LB_MB_N
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: LB_MB_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:) 	:: LB_MB_C

	REAL(KIND=8),PARAMETER									:: H_COMP  = 0.92D0					! percentage of SW that is H+
	REAL(KIND=8),PARAMETER									:: He_COMP = 1.0D0-H_COMP		! percentage of SW that is He++

	INTEGER				:: Max_Click   ! total number of clicks to transport
	INTEGER				:: Tot_Click   ! total number of clicks to write to files
	INTEGER				:: Click_Skip	 ! number of clicks to skip for output
  INTEGER		    :: N_Part      ! Total number of test particles
  INTEGER		    :: Num_Hist    ! Total number of histogram bins
	INTEGER				:: ENA_TYPE	   ! Type of ENA 0 => H | 1 => He
	INTEGER				:: ION_METH    ! Type of ion energy loss
														   ! 0 => none | 1 => Butler | 2 => Bethe
	INTEGER				:: STAR_METH   ! Type of star initial conditions to use
	INTEGER				:: ENG_METH    ! 0 => mono | 1 => Prob Den
	REAL(KIND=8)	:: MONO_ENG    ! Mono Energy [eV]
	INTEGER				:: N_r_zone    ! Number of radial zones
	INTEGER				:: N_Engy_Dist ! Number of energy zones
	REAL(KIND=8)	:: Max_T_zone  ! Maximum time zone [days]
	REAL(KIND=8)	:: Max_R_zone  ! Maximum radial zone [PC]
	REAL(KIND=8)	:: Max_SR_zone ! Maximum small grid radial zone [PC]
	REAL(KIND=8)	:: Engy_Dist_Fn! Maximum energy [eV]
	REAL(KIND=8)	:: dr_zone     ! differential for zone distance [AU]
	REAL(KIND=8)	:: d_Engy_Dist ! differential for energy [eV]
	INTEGER				:: DO_RND_In_EN! Do random initial energy test
	INTEGER				:: DO_RND_AN   ! Do random scattering angle test
	INTEGER				:: QM_ON			 ! Use QM cross sections if = 1
															 ! Use HS cross sections if = 0

	CHARACTER(Len=3)	:: Proj
	REAL(KIND=8)			:: MP			 ! projectils mass

CONTAINS

SUBROUTINE read_LB_MB
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read file containing Maxwell-Boltzmann
! energy distributions for the local 
! bubble with 10^6 K plasma
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info	

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER					:: i

	IF ( myid .EQ. 0 ) THEN
		OPEN(31, FILE="../Tables/Local_Bubble_MaxBoltz_Energy_Dist.dat", STATUS="old", ACTION="read")
		READ(31,*) LB_MB_N  ! number of lines in file
		ALLOCATE( LB_MB_E(LB_MB_N), LB_MB_C(LB_MB_N) )
		DO i=1,LB_MB_N
			READ(31,*) LB_MB_E(i), LB_MB_C(i)
		END DO ! i
		CLOSE(31)
	END IF ! myid = 0

	CALL MPI_BCAST( LB_MB_N, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	IF ( myid .NE. 0 ) ALLOCATE( LB_MB_E(LB_MB_N), LB_MB_C(LB_MB_N) )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  CALL MPI_BCAST( LB_MB_E, LB_MB_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )	
  CALL MPI_BCAST( LB_MB_E, LB_MB_N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )	
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE read_LB_MB

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_LB_MB_Energy(LB_Energy)

	IMPLICIT NONE

	! Output
	REAL(KIND=8)		:: LB_Energy ! [eV]

	! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: randy, C_now, x0, y0, x1, y1
	INTEGER					:: i

	randy = lfg()
	i     = 1
	C_now = LB_MB_C(i)

	DO WHILE (C_now .LT. randy)
		i 		= i + 1
		C_now = LB_MB_C(i)	
	END DO ! C_now < randy

	x0 = LB_MB_C(i-1)
	y0 = LB_MB_E(i-1)
	x1 = LB_MB_C(i)
	y1 = LB_MB_E(i)

	LB_Energy = y0 + (y1-y0)*((randy-x0)/(x1-x0))

END SUBROUTINE get_LB_MB_Energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE lism_read_inputs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read input file ../LISM_keys.in and 
! save parameters to module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info
	USE physics_constants

	IMPLICIT NONE

  INCLUDE 'mpif.h'

	INTEGER 					:: N_head, i

	IF (myid .EQ. 0) THEN

		OPEN(UNIT=33,FILE='../Inputs/LISM_keys.in',STATUS='OLD')
		N_head = 5
		DO i=1,N_head
			READ(33,*)
		END DO

		READ(33,*) ENA_TYPE 					! Int
		READ(33,*) 
		READ(33,*) N_Part
		READ(33,*) 
		READ(33,*) ION_METH						! Int
		READ(33,*) 
		READ(33,*) ENG_METH						! Int
		READ(33,*) MONO_ENG						! Real
		READ(33,*) 
		READ(33,*) STAR_METH
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) Max_Click
		READ(33,*) Click_Skip
		READ(33,*) 
		READ(33,*) Num_Hist				! Int
		READ(33,*) 
		READ(33,*) Max_T_zone					! Max time [days]
		READ(33,*) Max_R_zone					! Max distance [PC]
		READ(33,*) Max_SR_zone					! Max distance [PC]
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) 
		READ(33,*) DO_RND_In_EN				! Int
		READ(33,*) 
		READ(33,*) DO_RND_AN					! Int
		READ(33,*) 

		CLOSE(33)	
	
		CALL lism_inputs_write

	END IF

	!! Broadcast the input arrays to rest of ranks
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( MONO_ENG,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Max_T_zone,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Max_R_zone,   1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( ENA_TYPE,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( ION_METH,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( ENG_METH,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( STAR_METH,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Num_Hist,     1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( DO_RND_In_EN, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( DO_RND_AN,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( N_Part,       1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Max_Click,    1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( Click_Skip,   1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	Tot_Click     = Max_Click/Click_Skip

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (ENA_TYPE .EQ. 0) THEN
		Proj = 'H  '
		MP   = 1.0D0
	ELSE IF (ENA_TYPE .EQ. 1) THEN
		Proj = 'He '
		MP   = 4.0D0
	ELSE
		WRITE(*,*) 'IN lism_read_inputs ENA TYPE: ', ENA_TYPE, ' UNKNOWN!!'
	END IF

END SUBROUTINE lism_read_inputs

SUBROUTINE lism_inputs_write
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Write parameters from LISM_keys.in
! to screen
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE mpi_info	
	
	IMPLICIT NONE

	IF ( myid .EQ. 0 ) THEN

	50 format( A )
	51 format( A,ES10.2,A )
	52 format( A,I7,A )

	WRITE(*,50) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'	
	WRITE(*,50) '$$$$$        LISM INPUT PARAMETERS         $$$$$'	
	WRITE(*,50) '$$$$$--------------------------------------$$$$$'	
	IF (ENA_TYPE .EQ. 0) THEN
		WRITE(*,50) '$$$$$        ENA TYPE:      H              $$$$$'
	ELSE 
		WRITE(*,50) '$$$$$        ENA TYPE:      He             $$$$$'
	END IF
	IF (ION_METH .EQ. 0) THEN
		WRITE(*,50) '$$$$$        ION LOSS:      NONE           $$$$$'
	ELSE IF (ION_METH .EQ. 1) THEN
		WRITE(*,50) '$$$$$        ION LOSS:      BUTLER         $$$$$'
	ELSE IF (ION_METH .EQ. 2) THEN
		WRITE(*,50) '$$$$$        ION LOSS:      BETHE          $$$$$'
	END IF	
	IF (ENG_METH .EQ. 0) THEN
		WRITE(*,51) ' $$$$$        MONO ENERGY: ', MONO_ENG, '       $$$$$'
	ELSE 
		WRITE(*,50) '$$$$$        INIT ENERGY: PROB DIST        $$$$$'
	END IF
	IF (STAR_METH .EQ. 0) THEN
		WRITE(*,50) '$$$$$        STAR METHOD: SUN              $$$$$'
	ELSE IF (STAR_METH .EQ. 1) THEN
		WRITE(*,50) '$$$$$        STAR METHOD: SUN+60           $$$$$'
	ELSE IF (STAR_METH .EQ. 2) THEN
		WRITE(*,50) '$$$$$        STAR METHOD: Sun+60+Boundary  $$$$$'
	ELSE IF (STAR_METH .EQ. 3) THEN
		WRITE(*,50) '$$$$$        STAR METHOD: 60 Nearest       $$$$$'
	ELSE IF (STAR_METH .EQ. 4) THEN
		WRITE(*,50) '$$$$$        STAR METHOD: Boundary         $$$$$'
	END IF
	WRITE(*,50) '$$$$$--------------------------------------$$$$$'	
	WRITE(*,50) '$$$$$        GRID PARAMETERS               $$$$$'
	WRITE(*,50) '$$$$$--------------------------------------$$$$$'	
	WRITE(*,52) '$$$$$  MAX CLICK        :   ', Max_Click, '        $$$$$'
	WRITE(*,52) '$$$$$  CLICK SKIP       :   ', Click_Skip, '        $$$$$'
	WRITE(*,52) '$$$$$  N GRID           :   ', Num_Hist, '        $$$$$'
	WRITE(*,51) '$$$$$  R GRID END [PC]  :   ', Max_R_zone, '     $$$$$'
	WRITE(*,51) '$$$$$  T GRID END [yrs] :   ', Max_T_zone, '     $$$$$'
	WRITE(*,50) '$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$'	

	END IF

END SUBROUTINE lism_inputs_write

SUBROUTINE lism_grid_write
	USE mpi_info
	USE physics_constants

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	IF ( myid .EQ. 0 ) THEN
		!! print grid info to screen
    WRITE(*,*)
    WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    WRITE(*,*) '%%%%%              GRID INFO              %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       R_end [LY]: ', Max_r_zone*AUtoLY, '        %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       dR    [LY]: ', dr_zone*AUtoLY, '        %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       E_end [eV]: ', Engy_Dist_Fn, '        %%%%%'
    WRITE(*,'(A,ES10.2,A)') ' %%%%%       dE    [eV]: ', d_Engy_Dist, '        %%%%%'
    WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    WRITE(*,*)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE lism_grid_write

END MODULE lism

