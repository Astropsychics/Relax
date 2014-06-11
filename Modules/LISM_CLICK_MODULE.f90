
MODULE lism_click
	USE lism, ONLY : N_Part, Click_Skip, Max_Click, Proj, Num_Hist, get_LB_MB_Energy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CLICK_0/1 holds all current and next data 
! CLICK( MC, z )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z = 1   time [sec]
! z = 2 	energy [eV]
! z = 3   x position [AU]
! z = 4   y position [AU]
! z = 5   z position [AU]
! z = 6   r distance from star [AU]
! z = 7   ux
! z = 8   uy
! z = 9   uz
! z = 10  theta [rad]
! z = 11  phi [rad]
! z = 12  total number of collisions for MC particle
! z = 13  total number of elastic atom-atom collisions for MC particle
! z = 14  total number of CX atom-ion collisions for MC particle
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: CLICK_0
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: CLICK_1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! X_CLICK hols the reference vectors for histgram data
! X_CLICK( histogram_bin, z )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! z = 1  time [sec]
! z = 2  energy [eV]
! z = 3  x position [AU]
! z = 4  y position [AU]
! z = 5  z position [AU]
! z = 6  r position [AU]
! z = 7  ux unit velocity
! z = 8  uy unit velocity
! z = 9  uz unit velocity
! z = 10  theta range
! z = 11 phi range
! z = 12 Num Collisions
! z = 13 energy loss [eV]
! z = 14  small x position [AU]
! z = 15  small y position [AU]
! z = 16  small z position [AU]
! z = 17  small r position [AU]
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: X_CLICK
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER,ALLOCATABLE,DIMENSION(:)								:: TRANS
	INTEGER,ALLOCATABLE,DIMENSION(:)								:: R_TRANS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data save arrays for all clicks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_X_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_Y_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_Z_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_R_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_N_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_NE_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_NC_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_R_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_X_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_Y_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_E_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_E_X_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_E_Y_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_C_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_C_X_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: A_C_Y_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_ALL
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)					:: ROOT_ALL2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Data save arrays for individual clicks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_X_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_Y_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_Z_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_R_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_N_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_NE_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_NC_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_R_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_X_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: C_Y_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: Cs_R_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: Cs_E_T
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: Cs_X_Y
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: Cs_X_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: Cs_Y_Z
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:,:)				:: ROOT_CLICK


CONTAINS

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE allocate_lism_click
	USE lism 
	USE mpi_info

	IMPLICIT NONE

  REAL(KIND=8),PARAMETER 	:: GB_per_real = 8.0D-9
	INTEGER									:: N

	ALLOCATE( CLICK_0( N_Part,   14 ) )
	ALLOCATE( CLICK_1( N_Part,   14 ) )
	ALLOCATE( X_CLICK( Num_Hist, 17 ) )
	ALLOCATE( TRANS( N_Part ) )
	ALLOCATE( R_TRANS( N_Part ) )
	ALLOCATE( A_E_T( Num_Hist, Num_Hist ) )
	ALLOCATE( A_X_T( Num_Hist, Num_Hist ) )
	ALLOCATE( A_Y_T( Num_Hist, Num_Hist ) )
	ALLOCATE( A_Z_T( Num_Hist, Num_Hist ) )
	ALLOCATE( A_R_T( Num_Hist, Num_Hist ) )
	ALLOCATE( A_E_X_Y( Num_Hist, Num_Hist ) )
	ALLOCATE( A_E_X_Z( Num_Hist, Num_Hist ) )
	ALLOCATE( A_E_Y_Z( Num_Hist, Num_Hist ) )
	ALLOCATE( A_C_X_Y( Num_Hist, Num_Hist ) )
	ALLOCATE( A_C_X_Z( Num_Hist, Num_Hist ) )
	ALLOCATE( A_C_Y_Z( Num_Hist, Num_Hist ) )
	ALLOCATE( C_E_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_X_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_Y_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_Z_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_R_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_X_Y( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_X_Z( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( C_Y_Z( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( Cs_E_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( Cs_R_T( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( Cs_X_Y( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( Cs_X_Z( Tot_Click, Num_Hist, Num_Hist ) )
	ALLOCATE( Cs_Y_Z( Tot_Click, Num_Hist, Num_Hist ) )

	N = 2*N_Part*14 + Num_Hist*17 + N_Part + 5*Num_Hist*Num_Hist + &
	&   8*Tot_Click*Num_Hist*Num_Hist

	TRANS(:)   = 0
	R_TRANS(:) = 0

	A_E_X_Y = 0.0D0
	A_E_X_Z = 0.0D0
	A_E_Y_Z = 0.0D0
	A_C_X_Y = 0.0D0
	A_C_X_Z = 0.0D0
	A_C_Y_Z = 0.0D0
	A_E_T = 0.0D0
	C_E_T = 0.0D0
	C_X_T = 0.0D0
	C_Y_T = 0.0D0
	C_Z_T = 0.0D0
	C_R_T = 0.0D0
	C_X_Y = 0.0D0
	C_X_Z = 0.0D0
	C_Y_Z = 0.0D0
	Cs_E_T = 0.0D0
	Cs_R_T = 0.0D0
	Cs_X_Y = 0.0D0
	Cs_X_Z = 0.0D0
	Cs_Y_Z = 0.0D0

	IF (myid .EQ. 0) THEN
		WRITE(*,'(A,I7,A)') '####  ', Tot_Click, ' Write Clicks Written           ####'
		WRITE(*,'(A,ES10.2,A)') '####  ', N*GB_per_real, ' GB Used For Clicks          ####'
	END IF

END SUBROUTINE allocate_lism_click

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

SUBROUTINE clean_lism_click

	DEALLOCATE( CLICK_0 )
	DEALLOCATE( CLICK_1 )
	DEALLOCATE( X_CLICK )
	DEALLOCATE( TRANS   )
	DEALLOCATE( R_TRANS   )
	DEALLOCATE( A_E_X_Y )
	DEALLOCATE( A_E_X_Z )
	DEALLOCATE( A_E_Y_Z )
	DEALLOCATE( A_C_X_Y )
	DEALLOCATE( A_C_X_Z )
	DEALLOCATE( A_C_Y_Z )
	DEALLOCATE( A_E_T ) 
	DEALLOCATE( A_X_T ) 
	DEALLOCATE( A_Y_T ) 
	DEALLOCATE( A_Z_T ) 
	DEALLOCATE( A_R_T ) 
	DEALLOCATE( C_E_T ) 
	DEALLOCATE( C_X_T ) 
	DEALLOCATE( C_Y_T ) 
	DEALLOCATE( C_Z_T ) 
	DEALLOCATE( C_R_T ) 
	DEALLOCATE( C_X_Y ) 
	DEALLOCATE( C_X_Z ) 
	DEALLOCATE( C_Y_Z ) 
	DEALLOCATE( Cs_E_T ) 
	DEALLOCATE( Cs_R_T ) 
	DEALLOCATE( Cs_X_Y ) 
	DEALLOCATE( Cs_X_Z ) 
	DEALLOCATE( Cs_Y_Z ) 

END SUBROUTINE clean_lism_click

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set initial value for click data

SUBROUTINE set_init_lism_click_data_Vsw
	USE physics_constants, ONLY : PI, PCTOAU
	USE lism, ONLY : Proj, N_Part, STAR_METH
	USE ena
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Input
	REAL(KIND=8)		:: h		!! Initial height [m]
	REAL(KIND=8)		:: SZA	!! Solar Zenith Angle [rad]

	!! Internal
	INTEGER					:: i
	REAL(KIND=8)		:: v		!! velocity [m/s]
	REAL(KIND=8)		:: lfg	!! random number generator
	REAL(KIND=8)		:: theta, phi, r_phi, r_theta, E0, MP, rh
	REAL(KIND=8)		:: x_now(3), u_now, C1, C2, C3, Norm, dr1, dr2, dr3

	IF (myid .EQ. 0) THEN
		IF (TRIM(PROJ) .EQ. 'H') THEN
			MP = 1.0D0
		ELSE IF (TRIM(PROJ) .EQ. 'He') THEN
			MP = 4.0D0
		ELSE
			WRITE(*,*) 'Proj: ', PROJ, ' not known!!'
		END IF
		OPEN(UNIT=65,FILE='../Data/Initial_Energy.dat',ACCESS='APPEND')
		OPEN(UNIT=66,FILE='../Data/Initial_Position.dat',ACCESS='APPEND')
		DO i=1,N_Part

			!! get initial staring location
			CALL lism_ena_start( x_now )

			!! get intial energy
			IF ( STAR_METH .EQ. 4) THEN
				CALL get_LB_MB_Energy( E0 )
			ELSE 
				CALL rand_init_energy( MP, E0 )
			END IF

			WRITE(65,*) E0
			WRITE(66,*) x_now(1), x_now(2), x_now(3)

			r_phi 					 = lfg()
			r_theta 				 = lfg()
			phi 	= r_phi*2.0D0*PI
			theta = ACOS(2.0D0*r_theta-1.0D0)
			CLICK_0( i, 1 )  = 0.0D0								! Time [sec]
			CLICK_0( i, 2 )  = E0									  ! Energy [eV] 		
			CLICK_0( i, 3 )  = x_now(1)						  ! x0		
			CLICK_0( i, 4 )  = x_now(2)						  ! y0		
			CLICK_0( i, 5 )  = x_now(3)						  ! z0		
			CLICK_0( i, 6 )  = SQRT(x_now(1)**2+x_now(2)**2+x_now(3)**2) ! r0
			CLICK_0( i, 7 )  = SIN(theta)*COS(phi)	! ux0
			CLICK_0( i, 8 )  = SIN(theta)*SIN(phi)	! uy0
			CLICK_0( i, 9 )  = COS(theta)	 					! uz0
			CLICK_0( i, 10 ) = 0.0D0 							  ! theta scattering angle [rad]
			CLICK_0( i, 11)  = 0.0D0								! phi scattering angle [rad]
			CLICK_0( i, 12)  = 0.0D0								! number of total collisions
			CLICK_0( i, 13)  = 0.0D0								! number of elastic collisions
			CLICK_0( i, 14)  = 0.0D0								! number of cx collisions
		END DO
		CLOSE(65)
		CLOSE(66)
	END IF
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(CLICK_0,N_Part*14,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE set_init_lism_click_data_Vsw

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set initial value for click data

SUBROUTINE set_init_lism_click_data
	USE physics_constants, ONLY : PI
	USE lism, ONLY : Proj, N_Part
	USE ena
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Input
	REAL(KIND=8)		:: h		!! Initial height [m]
	REAL(KIND=8)		:: SZA	!! Solar Zenith Angle [rad]

	!! Internal
	INTEGER					:: i
	REAL(KIND=8)		:: v		!! velocity [m/s]
	REAL(KIND=8)		:: lfg	!! random number generator
	REAL(KIND=8)		:: theta, phi, r_phi, r_theta, E0, MP, rh
	REAL(KIND=8)		:: x_now(3), u_now

	IF (myid .EQ. 0) THEN
		IF (TRIM(PROJ) .EQ. 'H') THEN
			MP = 1.0D0
		ELSE IF (TRIM(PROJ) .EQ. 'He') THEN
			MP = 4.0D0
		ELSE
			WRITE(*,*) 'Proj: ', PROJ, ' not known!!'
		END IF
		OPEN(UNIT=65,FILE='../Data/Initial_Energy.dat',ACCESS='APPEND')
		OPEN(UNIT=66,FILE='../Data/Initial_Position.dat',ACCESS='APPEND')
		DO i=1,N_Part

			CALL rand_init_energy( MP, E0 )
			CALL lism_ena_start( x_now )
			WRITE(65,*) E0
			WRITE(66,*) x_now(1), x_now(2), x_now(3)

			r_phi 					 = lfg()
			r_theta 				 = lfg()
			phi 						 = r_phi*2.0D0*PI
			theta            = r_theta*PI 
			CLICK_0( i, 1 )  = 0.0D0								! Time [sec]
			CLICK_0( i, 2 )  = E0									  ! Energy [eV] 		
			CLICK_0( i, 3 )  = x_now(1)						  ! x0		
			CLICK_0( i, 4 )  = x_now(2)						  ! y0		
			CLICK_0( i, 5 )  = x_now(3)						  ! z0		
			CLICK_0( i, 6 )  = SQRT(x_now(1)**2+x_now(2)**2+x_now(3)**2) ! r0
			CLICK_0( i, 7 )  = SIN(theta)*COS(phi)  ! ux0
			CLICK_0( i, 8 )  = SIN(theta)*SIN(phi)  ! uy0
			CLICK_0( i, 9 )  = COS(theta) 					! uz0
			CLICK_0( i, 10 ) = 0.0D0 							  ! theta scattering angle [rad]
			CLICK_0( i, 11)  = 0.0D0								! phi scattering angle [rad]
			CLICK_0( i, 12)  = 0.0D0								! number of total collisions
			CLICK_0( i, 13)  = 0.0D0								! number of elastic collisions
			CLICK_0( i, 14)  = 0.0D0								! number of cx collisions
		END DO
		CLOSE(65)
		CLOSE(66)
	END IF
	
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(CLICK_0,N_Part*14,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

END SUBROUTINE set_init_lism_click_data

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! set reference vectors for click data

SUBROUTINE lism_click_ref_vectors
	USE physics_constants, ONLY : PI, MTOAU, PCTOAU, YRStoSEC
	USE lism , ONLY : Max_T_zone, Max_R_zone, Max_SR_zone
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER				:: i, MONO_ENGY
	REAL(KIND=8)	:: TEnd, EEnd, REnd, NEnd, dEEnd, INIT_ENGY, SREnd

	MONO_ENGY = 2
	INIT_ENGY = 1.0D3
	TEnd 	= Max_T_zone*YRStoSEC	! [sec]
	REnd 	= Max_R_zone*PCTOAU 	! [AU] 
	SREnd	= Max_SR_zone				 	! [AU] 
	NEnd 	= Max_Click/4
	dEEnd = 5.0D0
	IF (MONO_ENGY .EQ. 1) THEN
		EEnd = INIT_ENGY
	ELSE
		IF (TRIM(Proj) .EQ. 'H') THEN
			EEnd = 3.0D3
		ELSE IF (TRIM(Proj) .EQ. 'He') THEN
			EEnd = 4.0D0*3.0D3
		ELSE
			WRITE(*,*) 'PROJ: ', PROJ, ' not known! '
		END IF
	END IF

	IF (myid .EQ. 0) THEN
		DO i=1,Num_Hist
			X_CLICK(i,1)  = 0.0D0 + (i-1)*(TEnd/REAL(Num_Hist-1))					! time [sec]
			X_CLICK(i,2)  = 1.0D0 + (i-1)*(EEnd/REAL(Num_Hist-1))					! energy [eV]
			X_CLICK(i,3)  = -REnd + (i-1)*(2.0D0*REnd/REAL(Num_Hist-1))		! x position [AU]
			X_CLICK(i,4)  = -REnd + (i-1)*(2.0D0*REnd/REAL(Num_Hist-1))		! y position [AU]
			X_CLICK(i,5)  = -REnd + (i-1)*(2.0D0*REnd/REAL(Num_Hist-1))		! z position [AU]
			X_CLICK(i,6)  = 0.0D0 + (i-1)*(REnd/REAL(Num_Hist-1))	        ! r position [m]
			X_CLICK(i,7)  = -1.0D0+ (i-1)*(2.0D0/REAL(Num_Hist-1))				! ux unit velocity
			X_CLICK(i,8)  = -1.0D0+ (i-1)*(2.0D0/REAL(Num_Hist-1))				! uy unit velocity
			X_CLICK(i,9)  = -1.0D0+ (i-1)*(2.0D0/REAL(Num_Hist-1))				! uz unit velocity
			X_CLICK(i,10) = 0.0D0 + (i-1)*(PI/REAL(Num_Hist-1))						! theta range
			X_CLICK(i,11) = 0.0D0 + (i-1)*(2.0D0*PI/REAL(Num_Hist-1))			! phi range
			X_CLICK(i,12) = 0.0D0 + (i-1)*(NEnd/REAL(Num_Hist-1))					! number of collisions
			X_CLICK(i,13) = 0.0D0 + (i-1)*(dEEnd/REAL(Num_Hist-1))				! energy loss [eV]
			X_CLICK(i,14) = -SREnd+ (i-1)*(2.0D0*SREnd/REAL(Num_Hist-1))	! sm x position [m]
			X_CLICK(i,15) = -SREnd+ (i-1)*(2.0D0*SREnd/REAL(Num_Hist-1))	! sm y position [m]
			X_CLICK(i,16) = -SREnd+ (i-1)*(2.0D0*SREnd/REAL(Num_Hist-1))	! sm z position [m]
			X_CLICK(i,17) = 0.0D0 + (i-1)*(SREnd/REAL(Num_Hist-1))	      ! sm r position [m]
		END DO
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( X_CLICK, Num_Hist*13, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )

END SUBROUTINE lism_click_ref_vectors

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! write the reference vectors to files
! in data directory
SUBROUTINE write_lism_click_ref_vectors

	IMPLICIT NONE

	INTEGER		:: i

	OPEN(UNIT=90, FILE="../Data/lism_X_time.dat", 				ACCESS="APPEND")
	OPEN(UNIT=91, FILE="../Data/lism_X_energy.dat", 			ACCESS="APPEND")
	OPEN(UNIT=92, FILE="../Data/lism_X_xyz.dat", 					ACCESS="APPEND")
	OPEN(UNIT=93, FILE="../Data/lism_X_r.dat", 						ACCESS="APPEND")
	OPEN(UNIT=94, FILE="../Data/lism_X_unit_vel.dat", 		ACCESS="APPEND")
	OPEN(UNIT=95, FILE="../Data/lism_X_theta.dat", 				ACCESS="APPEND")
	OPEN(UNIT=96, FILE="../Data/lism_X_phi.dat", 					ACCESS="APPEND")
	OPEN(UNIT=97, FILE="../Data/lism_X_Ncoll.dat",       	ACCESS="APPEND")
	OPEN(UNIT=98, FILE="../Data/lism_X_energy_loss.dat", 	ACCESS="APPEND")
	OPEN(UNIT=99, FILE="../Data/lism_small_X_xyz.dat",		ACCESS="APPEND")
	OPEN(UNIT=100,FILE="../Data/lism_small_X_r.dat", 			ACCESS="APPEND")

	DO i=1,Num_Hist
		WRITE(90,*) X_CLICK(i,1)
		WRITE(91,*) X_CLICK(i,2)
		WRITE(92,*) X_CLICK(i,3)
		WRITE(93,*) X_CLICK(i,6)
		WRITE(94,*) X_CLICK(i,7)
		WRITE(95,*) X_CLICK(i,10)
		WRITE(96,*) X_CLICK(i,11)
		WRITE(97,*) X_CLICK(i,12)
		WRITE(98,*) X_CLICK(i,13)
		WRITE(99,*) X_CLICK(i,16)
		WRITE(100,*) X_CLICK(i,17)
	END DO

	CLOSE(90)
	CLOSE(91)
	CLOSE(92)
	CLOSE(93)
	CLOSE(94)
	CLOSE(95)
	CLOSE(96)
	CLOSE(97)
	CLOSE(98)
	CLOSE(99)
	CLOSE(100)

END SUBROUTINE write_lism_click_ref_vectors

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



END MODULE lism_click

