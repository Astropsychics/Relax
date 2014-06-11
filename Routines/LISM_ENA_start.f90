
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Read parameters for ENA start probabilies and 
! save them in the lism module
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lism_ena_table_read
	USE lism
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	INTEGER					:: i, N_H, N_He

	IF (myid .EQ. 0) THEN
		OPEN(UNIT=50,FILE='../Tables/LISM_H_ENA_Start.dat')	
		OPEN(UNIT=60,FILE='../Tables/LISM_He_ENA_Start.dat')	

		READ(50,*) N_H	
		READ(60,*) N_He	

		ALLOCATE( H_ENA_PROB(N_H), H_ENA_DIST(N_H) ) 
		ALLOCATE( He_ENA_PROB(N_He), He_ENA_DIST(N_He) ) 

		DO i=1,N_H
			READ(50,*) H_ENA_PROB(i), H_ENA_DIST(i)	
		END DO

		DO i=1,N_He
			READ(60,*) He_ENA_PROB(i), He_ENA_DIST(i)	
		END DO

		CLOSE(50)
		CLOSE(60)
	END IF

	CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

	CALL MPI_BCAST(N_H,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST(N_He,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

	IF (myid .NE. 0) THEN
		ALLOCATE( H_ENA_PROB(N_H), H_ENA_DIST(N_H) ) 
		ALLOCATE( He_ENA_PROB(N_He), He_ENA_DIST(N_He) ) 
	END IF

	CALL MPI_BCAST( H_ENA_PROB, N_H, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST( He_ENA_PROB, N_He, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST( H_ENA_DIST, N_H, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)
	CALL MPI_BCAST( He_ENA_DIST, N_He, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD,ierr)

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END SUBROUTINE lism_ena_table_read

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Determine starting location and velocity
! for ENA in LISM
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE lism_ena_start( R )
	USE physics_constants
	USE lism

	IMPLICIT NONE

	!! Outputs
	REAL(KIND=8)		:: R(3)		! starting location [x,y,z] in [m]

	!! Internal
	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: phi, theta, dist, randy, randy1, randy2, Star(3), P_Star, randy3
	REAL(KIND=8)		:: N_Near_Stars, N_Boundary_Stars, R_Boundary, Norm
	REAL(KIND=8)		:: dr1, dr2, dr3, Cr1, Cr2, Cr3, randy_phi, randy_theta, u, v
	INTEGER					:: i, N

	N_Near_Stars 			= 61.0D0
	N_Boundary_Stars 	= 100.0D0
	P_Star						= N_Near_Stars/(N_Near_Stars+N_Boundary_Stars)
	R_Boundary				= 4.5D0*PCTOAU

	u 					= lfg()
	v 					= lfg()
	randy_phi 	= 2.0D0*PI*u
	randy_theta = ACOS(2.0D0*v-1.0D0)

	IF (ENA_TYPE .EQ. 0) THEN
		Proj   = 'H  '
		randy2 = lfg()
		N      = 0
		i      = 1
		DO WHILE (N .EQ. 0)
			IF ( (H_ENA_PROB(i) .GE. randy2) .OR. (i .GE. SIZE(H_ENA_PROB(:))) ) THEN
				dist = H_ENA_DIST(i)
				N    = 1
			END IF
			i = i+1	
		END DO
	ELSE
		Proj   = 'He '
		randy2 = lfg()
		N      = 0
		i      = 1
		DO WHILE (N .EQ. 0)
			IF ( (He_ENA_PROB(i) .GE. randy2) .OR. (i .GE. SIZE(He_ENA_PROB(:))) ) THEN
				dist = He_ENA_DIST(i)
				N    = 1
			END IF
			i = i+1	
		END DO
	END IF

	IF (STAR_METH .EQ. 0) THEN
		R(1) = dist*SIN(randy_theta)*COS(randy_phi)
		R(2) = dist*SIN(randy_theta)*SIN(randy_phi)	
		R(3) = dist*COS(randy_theta)
	ELSE IF ( (STAR_METH .EQ. 1) .OR. (STAR_METH .EQ. 3) ) THEN
		!! star start
		CALL lism_star_start(Star)
		R(1) = Star(1) + dist*SIN(randy_theta)*COS(randy_phi)
		R(2) = Star(2) + dist*SIN(randy_theta)*SIN(randy_phi)	
		R(3) = Star(3) + dist*COS(randy_theta)
	ELSE IF ( STAR_METH .EQ. 2 ) THEN
		randy3  = lfg()
		IF ( randy3 .LE. P_Star ) THEN
			!! star start
			CALL lism_star_start(Star)
			R(1) = Star(1) + dist*SIN(randy_theta)*COS(randy_phi)
			R(2) = Star(2) + dist*SIN(randy_theta)*SIN(randy_phi)	
			R(3) = Star(3) + dist*COS(randy_theta)
		ELSE 
			!! boundary start
			R(1) = R_Boundary*SIN(randy_theta)*COS(randy_phi)
			R(2) = R_Boundary*SIN(randy_theta)*SIN(randy_phi)
			R(3) = R_Boundary*COS(randy_theta)
		END IF
	ELSE IF (STAR_METH .EQ. 4) THEN
		!! boundary start
		R(1) = R_Boundary*SIN(randy_theta)*COS(randy_phi)
		R(2) = R_Boundary*SIN(randy_theta)*SIN(randy_phi)
		R(3) = R_Boundary*COS(randy_theta)
	END IF	

END SUBROUTINE lism_ena_start

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE lism_star_start( R )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Start the ENA from a rand star
! within 4 PC of Earth
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE lism, ONLY : STAR_METH

	IMPLICIT NONE

	REAL(KIND=8) 		:: R(3) 	! [AU]
	REAL(KIND=8)		:: S(60,3) ! S( [star] , [x,y,z] )

	REAL(KIND=8)		:: lfg
	REAL(KIND=8)		:: randy	
	INTEGER					:: star_i

	S(1,:)  = (/ -2.095962D05, -1.624248D05, 4.072926D04 /)
	S(2,:)  = (/ -1.143598D05, 1.075290D03, 2.512421D05 /)
	S(3,:)  = (/ -1.135135D05, 1.729598D03, 2.516220D05 /)
	S(4,:)  = (/ -5.378940D03, 4.764384D03, -3.770376D05 /)
	S(5,:)  = (/ 2.657614D05, 2.520116D05, 3.287786D05 /)
	S(6,:)  = (/ 6.477696D04, -5.134554D04, -5.177428D05 /)
	S(7,:)  = (/ -2.106394D05, -1.986520D05, 4.591120D05 /)
	S(8,:)  = (/ 3.192950D05, -1.265410D05, 4.321677D05 /)
	S(9,:)  = (/ 3.839154D05, -5.982054D04, 4.415209D05 /)
	S(10,:) = (/ 1.589607D05, -4.700438D04, 5.893881D05 /)
	S(11,:) = (/ -5.688958D05, -2.946903D05, 1.250142D05 /)
	S(12,:) = (/ 6.549388D05, -1.155119D05, 2.232429D04 /)
	S(13,:) = (/ -1.160037D05, -1.436507D05, 6.537610D05 /)
	S(14,:) = (/ 2.563989D05, 4.044685D05, 4.974881D05 /)
	S(15,:) = (/ -6.155008D05, -2.226965D05, -2.813948D05 /)
	S(16,:) = (/ -5.539531D04, 3.492884D05, -6.283920D05 /)
	S(17,:) = (/ -3.006930D05, 1.954074D05, 6.256504D05 /)
	S(18,:) = (/ -3.057686D05, 1.961785D05, 6.229433D05 /)
	S(19,:) = (/ 3.373347D05, 6.446102D05, 4.369229D04 /)
	S(20,:) = (/ 3.340387D05, 6.461558D05, 4.611718D04 /)
	S(21,:) = (/ -8.572071D04, -7.294837D05, 2.995375D04 /)
	S(22,:) = (/ -4.089358D05, -6.008182D05, -1.758873D05 /)
	S(23,:) = (/ -2.794086D05, -6.715568D05, -1.734638D05 /)
	S(24,:) = (/ 1.206360D04, -5.350780D04, 7.458739D05 /)
	S(25,:) = (/ -4.639588D05, -5.662230D05, 1.710420D05 /)
	S(26,:) = (/ -5.429981D05, -3.648861D05, -3.834902D05 /)
	S(27,:) = (/ -1.590468D05, 1.400744D05, 7.373874D05 /)
	S(28,:) = (/ 1.228515D05, -3.639465D05, -6.812005D05 /)
	S(29,:) = (/ -2.285190D05, 2.039085D05, -7.297335D05 /)
	S(30,:) = (/ 3.360228D04, -3.364572D05, -7.195114D05 /)
	S(31,:) = (/ 2.974458D04, -3.360174D05, -7.198866D05 /)
	S(32,:) = (/ -3.349998D05, 2.394881D05, -6.952219D05 /)
	S(33,:) = (/ 1.363751D05, -2.881103D05, -7.488968D05 /)
	S(34,:) = (/ -5.338416D05, -4.126150D05, 4.697496D05 /)
	S(35,:) = (/ -2.262428D05, -2.562461D05, 7.580484D05 /)
	S(36,:) = (/ -1.366862D05, 4.953699D05, -6.552139D05 /)
	S(37,:) = (/ 7.985783D05, -3.345784D04, -2.717238D05 /)
	S(38,:) = (/ -7.128564D05, 4.986565D05, -8.388684D04 /)
	S(39,:) = (/ 5.357874D05, -1.512414D05, -6.937964D05 /)
	S(40,:) = (/ 1.843336D05, 8.277003D05, 3.014998D05 /)
	S(41,:) = (/ -8.197087D05, 1.439956D05, 3.557503D05 /)
	S(42,:) = (/ 1.690371D05, -7.935777D05, 4.284670D05 /)
	S(43,:) = (/ 6.525355D05, 1.419055D05, -6.551836D05 /)
	S(44,:) = (/ 8.078924D04, -3.031107D05, 8.821694D05 /)
	S(45,:) = (/ 1.423856D05, 9.009877D05, -2.122433D05 /)
	S(46,:) = (/ -8.500280D05, -2.085643D04, 3.930297D05 /)
	S(47,:) = (/ -8.498246D05, -2.439489D04, 3.932659D05 /)
	S(48,:) = (/ -8.500280D05, -2.085643D04, 3.930297D05 /)
	S(49,:) = (/ 7.834772D04, 6.858229D05, 6.506758D05 /)
	S(50,:) = (/ -3.551394D05, -1.938113D05, -8.622047D05 /)
	S(51,:) = (/ -3.339720D04, 2.977358D05, -9.208989D05 /)
	S(52,:) = (/ 7.739290D04, 9.476856D04, -9.624976D05 /)
	S(53,:) = (/ 3.259818D04, 3.841348D04, 9.859056D05 /)
	S(54,:) = (/ -8.908796D05, 1.212261D05, -4.405705D05 /)
	S(55,:) = (/ -8.895677D05, 1.142591D05, -4.450604D05 /)
	S(56,:) = (/ -3.155060D05, 6.137422D05, -7.268207D05 /)
	S(57,:) = (/ -3.011036D05, -4.334208D05, 8.590263D05 /)
	S(58,:) = (/ -3.085843D05, 6.010776D04, 9.674310D05 /)
	S(59,:) = (/ 5.439944D05, -8.669371D05, 3.809411D04 /)
	S(60,:) = (/ -9.943333D05, 2.190975D05, -1.118828D05 /)
	
	randy  = lfg()

	IF (STAR_METH .EQ. 3) THEN
		star_i = FLOOR(60.0D0*randy) + 1
		R      = S(star_i, :)
	ELSE IF ( (STAR_METH .EQ. 1) .OR. (STAR_METH .EQ. 2) ) THEN
		star_i = FLOOR(61.0D0*randy)			
		IF (star_i .EQ. 0) THEN
			R = (/ 0.0D0, 0.0D0, 0.0D0 /)
		ELSE
			R      = S(star_i, :)
		END IF
	END IF 

END SUBROUTINE lism_star_start




