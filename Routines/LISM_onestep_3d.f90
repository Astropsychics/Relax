
SUBROUTINE LISM_onestep_3d
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Experiment to calculate relaxation of ENAs in the LISM which
! orginate throught solar wind charge exchange collisions. 
! Made to be 'main' routine for experiment, in that most 
! computation takes place in other subroutines for ease of 
! addition of higher complexity to the problem in the future. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	USE physics_constants
	USE lism
	USE lism_click
	USE ENA, ONLY : L_ion, L_targ
	USE flux_mapping
	USE planet, ONLY : M_H, M_He
	USE elastic_atom_ion
	USE mpi_info	
	
	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Internal
	INTEGER																	:: N_MC, My_N_MC, MC, My_MC_Start, My_MC_End
	INTEGER																	:: q, zz, z, i, j, tot, ck
	INTEGER																	:: N_hist
	INTEGER																	:: max_collisions 
	INTEGER																	:: N_Coll_zone, Trans_tot
	INTEGER																	:: My_10_per_done
	INTEGER																	:: Thermalized, Root_Thermalized
	INTEGER,DIMENSION(np) 									:: buff_tot

	REAL,PARAMETER													:: MB_per_real = 8.0D-6
	REAL,PARAMETER													:: GB_per_real = 8.0D-9

	REAL(KIND=8)														:: lfg
	REAL(KIND=8)														:: N_Coll, N_CX_Coll_H, N_Atom_Coll_H, N_Coll_H
	REAL(KIND=8)														:: N_CX_Coll_He, N_Atom_Coll_He, N_Coll_He
	REAL(KIND=8)														:: root_N_CX, root_N_AT, root_N
	REAL(KIND=8)														:: ScattAng, Lab_ScattAng, randy
	REAL(KIND=8)														:: CX_CS, TCS, MFP, Den_now, ion_Den_now
	REAL(KIND=8)														:: E_now, E_nxt, E_temp, ds, dt
	REAL(KIND=8)														:: r1, r2, temp_t, temp_p, r0
	REAL(KIND=8)														:: theta_now, phi_now, theta_prev, phi_prev, phi
	REAL(KIND=8)														:: max_diss, tot_per, Click_Elapsed
	REAL(KIND=8)														:: Click_End, Click_hr, Click_min, Click_sec
	REAL(KIND=8)														:: Start_t, End_t, Full_Start, Full_End
	REAL(KIND=8)														:: Elapsed, Full_min, Full_sec
	REAL(KIND=8)														:: E0, x0, y0, z0, ux0, uy0, uz0, P0, T0
	REAL(KIND=8)														:: E1, x1, y1, z1, ux1, uy1, uz1, P1, T1
	REAL(KIND=8)														:: Den_tot, MR_H, MR_Hp, force, ion_path
	REAL(KIND=8)														:: Engy_Dist_In
	REAL(KIND=8)														:: Prob_atom, Prob_ion, Prob_tot, tt, tau

	REAL(KIND=8),DIMENSION(3)								:: u_now, u_nxt
	REAL(KIND=8),DIMENSION(3)								:: x_now, x_nxt
	REAL(KIND=8),DIMENSION(4)								:: XXX

	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: My_Data
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:) :: Root_Data
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: My_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:,:)	:: Root_Engy_Dist
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_CX_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: My_N_EL_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Diss
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Hist_x
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Hist_y 
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone_C
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: r_zone
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: E_zone
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_r_zone_E
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_r_zone_C
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_C_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_E_H
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_CX_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_C_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: Root_N_EL_E_He
	REAL(KIND=8),ALLOCATABLE,DIMENSION(:)		:: R_Final_Energy

	QM_ON = 1

	!! get inptus from file LISM_keys.in
	CALL lism_read_inputs	! parallel check

	!! read LB Maxwell Boltzmann tables
	CALL read_LB_MB

	!! set up initial energy tables
	CALL read_Vsw_table ! parallel check

	!! elastic atom-ion collisions
	CALL read_elastic_atom_ion_tables
!	CALL test_elastic_atom_ion

	!! Run tests if input file dictates such
	IF ( DO_RND_IN_EN .EQ. 1) CALL test_rand_energy
	IF ( DO_RND_AN    .EQ. 1) CALL test_lin_rand_angle

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Use number of MC particles as specified in keys.in input file
	N_MC = N_Part

	!! Output formats
	11 FORMAT(I7, A)
	12 FORMAT(A, I4, A, I9, A)
	13 FORMAT(A, I4, A, F5.1, A)
	
	!! Set up number of MC particles per rank
	My_N_MC 			 = N_MC / np
	My_10_per_done = My_N_MC / 10

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Collect My_N_MC to root to make sure full number of particles are being computed
	CALL MPI_REDUCE( My_N_MC, tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BCAST( tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
	
	!! Add 1 MC particle to all ranks N_MC if not enough
	IF (tot .LT. N_MC) THEN
		My_N_MC = My_N_MC + 1 
		CALL MPI_REDUCE( My_N_MC, tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
		CALL MPI_BCAST( tot, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
		N_MC = tot
	END IF	

	N_Part      = N_MC
	My_MC_Start = myid*My_N_MC + 1
	My_MC_End   = (myid + 1)*My_N_MC

	Thermalized 			= 0	
	Root_Thermalized 	= 0

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Write welcome messages to screen
	IF ( myid == 0 ) THEN
		WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		WRITE(*,*) 'STARTING LISM ENA SIMULATION'
		WRITE(*,'(I3,A)') np, ' Processors Being Used'
		WRITE(*,'(ES9.2,A)') REAL(N_MC), ' Total Particles Being Simulated'
		WRITE(*,'(A,F5.2)') ' Total Storage Array Memory Footprint [GB]: ', &
		&  (12.0D0*N_MC + np*(N_r_zone*N_Engy_Dist+N_Engy_Dist) + &
		&  N_r_zone*N_Engy_Dist + np*3.0D0*N_r_zone+2.0D0*N_r_zone)*GB_per_real
		WRITE(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
		WRITE(*,*)
	END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Write number of MC particles per processor to screen
	IF (myid .EQ. 0) THEN
		WRITE(*,'(A,ES9.2,A)') ' Computing ', REAL(My_N_MC), ' MC Particles per MPI rank' 
		WRITE(*,*)
	END IF

  CALL allocate_lism_click
  CALL lism_click_ref_vectors
  IF (myid .EQ. 0) CALL write_lism_click_ref_vectors

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! read starting location parameter tables
	CALL lism_ena_table_read ! parallel check

	!! initialized all test particles
	CALL set_init_lism_click_data_Vsw

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF (myid .EQ. 0) Full_Start = MPI_WTIME()

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!! Main DO loop
DO ck=1,Max_Click
!	WRITE(*,*) 'click: ', ck
	CLICK_1 = 0.0D0

	! set Trans to 1 if first click for particle
	IF (ck .EQ. 1) THEN
!		WRITE(*,*) 'rank: ', myid, ' st: ', My_MC_Start, ' end: ', My_MC_End
		DO MC=My_MC_Start,My_MC_End
			Trans(MC) = 1
		END DO
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
		R_Trans = 0
		CALL MPI_REDUCE (Trans,R_Trans,N_MC,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		IF (myid .EQ. 0) THEN
			Trans_tot = 0
			DO i=1,N_MC
				Trans_tot = Trans_tot + R_Trans(i)
			END DO
			WRITE(*,*) 'INITIAL ', Trans_tot, ' particles of ', N_MC, ' being transported'
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	END IF ! ck = 1

	DO MC=My_MC_Start,My_MC_End
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!! Initialize all parameters
		!! for current MC test particle	
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		E0  = CLICK_0( MC, 2 )
		x0  = CLICK_0( MC, 3 )
		y0  = CLICK_0( MC, 4 )
		z0  = CLICK_0( MC, 5 )
		r0  = CLICK_0( MC, 6 )
		ux0 = CLICK_0( MC, 7 )
		uy0 = CLICK_0( MC, 8 )
		uz0 = CLICK_0( MC, 9 )

		!! Start energy dependent thermalization loop
		IF ( TRANS(MC) .EQ. 1 ) THEN

			CALL lism_onestep_transport( E0, x0, y0, z0, ux0, uy0, uz0, &
			&														 E1, x1, y1, z1, ux1, uy1, uz1, dt )

			CLICK_1( MC, 1 ) = CLICK_0( MC, 1) + dt
			CLICK_1( MC, 2 ) = E1
			CLICK_1( MC, 3 ) = x1
			CLICK_1( MC, 4 ) = y1
			CLICK_1( MC, 5 ) = z1
			CLICK_1( MC, 6 ) = SQRT(x1**2+y1**2+z1**2)
			CLICK_1( MC, 7 ) = ux1
			CLICK_1( MC, 8 ) = uy1
			CLICK_1( MC, 9 ) = uz1
			CLICK_1( MC, 10:14) = 0.0D0

			!! If particle lost all energy in 
			!! transport, stop transporting particle
			!! Thermal in LIC is 6000 K = 0.5 eV
			!! Cutoff at 1 eV to account for upper MB tail
			IF (E1 .LT. 1.0D0) TRANS(MC) = 0

      IF ( MOD(ck,Click_Skip) .EQ. 0 ) THEN
	      CALL fill_lism_phase_spaces(ck/Click_Skip,MC)
      END IF

      CALL fill_all_lism_phase_spaces(MC)

		ELSE !! particle is not being transported
			CLICK_1( MC, : ) = CLICK_0( MC, : )
		END IF

		x_now(:) = CLICK_0(MC,3:5)
		u_now(:) = CLICK_0(MC,7:9)
		x_nxt(:) = CLICK_1(MC,3:5)
		u_nxt(:) = CLICK_1(MC,7:9)

		!! Check for particle crossing mapping surface
		CALL flux_map( x_nxt(:), u_nxt(:), CLICK_0(MC,2) )

	END DO ! MC

	CLICK_0 = CLICK_1

	IF ( myid .EQ. 0 ) THEN

	  IF ( ck .EQ. 10 ) THEN
  		Click_End     = MPI_WTIME()
  		Click_Elapsed = Max_Click*(Click_End - Full_Start)/10.0D0
  		Click_hr      = FLOOR(Click_Elapsed/3600.0D0)
  		Click_min     = FLOOR(MOD((Click_Elapsed/60.0D0),60.0))
  		Click_sec     = MOD(Click_Elapsed,60.0)
  		WRITE(*,*) 'Projected Simulation Time: '
 			WRITE(*,'(I4,A,I4,A,I4,A)') INT(Click_hr), ' hours  ', INT(Click_min), ' minutes  ', INT(Click_sec), ' seconds '
	  	WRITE(*,*)
  	END IF


	  IF ( MOD(ck,(Max_click/10)) .EQ. 0) THEN
!			WRITE(*,*) 'num_therm: ck: ', N_MC-ROOT_Thermalized, ck
 	 		IF ( (ck/REAL(Max_click)) .EQ. 0.1 ) THEN
    		WRITE(*,*) 'Percent Complete'
      	WRITE(*,*) '----------------'
      	WRITE(*,*) '[*         ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.2) THEN
      	WRITE(*,*) '[**        ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.3) THEN
      	WRITE(*,*) '[***       ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.4) THEN
      	WRITE(*,*) '[****      ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.5) THEN
      	WRITE(*,*) '[*****     ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.6) THEN
      	WRITE(*,*) '[******    ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.7) THEN
      	WRITE(*,*) '[*******   ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.8) THEN
      	WRITE(*,*) '[********  ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 0.9) THEN
      	WRITE(*,*) '[********* ]'
    	ELSE IF ( (ck/REAL(Max_click)) .EQ. 1.0) THEN
      	WRITE(*,*) '[**********]'
    	END IF
  	END IF
	END IF

  IF ( MOD(ck,(Max_click/10)) .EQ. 0) THEN
		CALL MPI_REDUCE (Trans,R_Trans,N_MC,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,ierr)
		IF (myid .EQ. 0) THEN
			Trans_tot = 0
			DO i=1,N_MC
				Trans_tot = Trans_tot + R_Trans(i)
			END DO
			WRITE(*,*) Trans_tot, ' particles of ', N_MC, ' being transported'
		END IF
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	END IF

!	IF (ROOT_THERMALIZED .GE. N_MC) THEN
!		WRITE(*,*) 'ALL PARTICLES THERMALIZED WITH CLICKS: ', ck, ROOT_THERMALIZED	
!		EXIT
!	END IF
	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

END DO ! CLICK

!WRITE(*,*) 'rank: ', myid, 'st: ', My_MC_Start, 'end: ', My_MC_End,' My_N_MC: ', My_N_MC, ' DFF: ', My_MC_End-My_MC_Start
!WRITE(*,*) 'randk: ', myid, 'L: ', SIZE(CLICK_0(My_MC_Start:My_MC_End,2))
!WRITE(*,*) 'rank: ', myid, ' CE: ', Click_0(My_MC_Start:My_MC_End,2)
IF (myid .EQ. 0) THEN
	ALLOCATE( R_Final_Energy( My_N_MC*np ) ) 
	R_Final_Energy(:) = 0.0D0
!	WRITE(*,*) 'Size R_Fin: ', SIZE(R_Final_Energy)
END IF
CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
CALL MPI_GATHER(CLICK_0(My_MC_Start:My_MC_End,2),My_N_MC,MPI_DOUBLE_PRECISION, &
&               R_Final_Energy,My_N_MC,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
IF (myid .EQ. 0) THEN
	OPEN(UNIT=55,FILE="../Data/Final_Energy.dat",ACCESS="APPEND")
	DO i=1,N_MC
		WRITE(55,*) R_Final_Energy(i)
	END DO
	CLOSE(55)
END IF

IF (myid .EQ. 0) WRITE(*,*) 'SIMULATION ENDED: ', ck

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	!! Write timing data from run to screen
	IF (myid .EQ. 0) THEN
		Full_End = MPI_WTIME()
		Elapsed  = Full_End - Full_Start
		Full_Min = Elapsed/60
		Full_Sec = MOD(Elapsed,60.0)
		WRITE(*,*) 
    WRITE(*,FMT="(A,I5,A,I2,A)") "LISM ENA RELAX complete in ",INT(Full_Min), " min ", INT(Full_Sec)," sec"
    WRITE(*,FMT="(A,I3,A,I5,A,I2,A)") "Time saved by using", np, " processors: ", &
    & INT((Elapsed*(np-1))/60)," min ", INT(MOD((Elapsed*(np-1)),60.0)), " sec"
    WRITE(*,*)
  END IF

	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	CALL write_lism_phase_spaces

	CALL clean_lism_click

END SUBROUTINE LISM_onestep_3d

