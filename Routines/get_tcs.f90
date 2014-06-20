
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! gather required tcs for 
! transport calculations
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_cx_tcs(E,CX_X_CO2,CX_X_CO,CX_X_N2,CX_X_H2,CX_X_O,CX_X_Ar,CX_X_He,CX_X_H)
	USE physics_constants, ONLY : BOHRTOM
	USE planet

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E		 			! Lab Frame [eV]
	
	!! Outputs
	REAL(KIND=8)	:: CX_X_CO2	! [m^2]
	REAL(KIND=8)	:: CX_X_CO		! [m^2]
	REAL(KIND=8)	:: CX_X_N2		! [m^2]
	REAL(KIND=8)	:: CX_X_H2		! [m^2]
	REAL(KIND=8)	:: CX_X_O		! [m^2]
	REAL(KIND=8)	:: CX_X_Ar		! [m^2]
	REAL(KIND=8)	:: CX_X_He		! [m^2]
	REAL(KIND=8)	:: CX_X_H		! [m^2]

	!! Internal
	REAL(KIND=8)			:: dummy
	CHARACTER(LEN=5)	:: ion

	IF ( PROJ .EQ. 'H ') ion = 'H_p  '
	IF ( PROJ .EQ. 'He') ion = 'He_pp'

	CALL cx_cross_sections( ion, 'H  ', 1, E, CX_X_H   )
	CALL cx_cross_sections( ion, 'He ', 1, E, CX_X_He  )
	CALL cx_cross_sections( ion, 'O  ', 1, E, CX_X_O   )
	CALL cx_cross_sections( ion, 'O  ', 1, E, CX_X_Ar  )
	CALL cx_cross_sections( ion, 'H  ', 1, E, CX_X_H2  )
	CALL cx_cross_sections( ion, 'CO ', 1, E, CX_X_N2  )
	CALL cx_cross_sections( ion, 'CO ', 1, E, CX_X_CO  )
	CALL cx_cross_sections( ion, 'CO2', 1, E, CX_X_CO2 )

END SUBROUTINE get_cx_tcs

!######################
!######################

SUBROUTINE get_ion_atom_tcs(E,TCS_X_CO2,TCS_X_CO,TCS_X_N2,TCS_X_H2,TCS_X_O,TCS_X_Ar,TCS_X_He,TCS_X_H)
	USE physics_constants, ONLY : BOHRTOM
	USE planet

	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E					! Lab Frame [eV]
	
	!! Outputs
	REAL(KIND=8)	:: TCS_X_CO2	! [m^2]
	REAL(KIND=8)	:: TCS_X_CO		! [m^2]
	REAL(KIND=8)	:: TCS_X_N2		! [m^2]
	REAL(KIND=8)	:: TCS_X_H2		! [m^2]
	REAL(KIND=8)	:: TCS_X_O		! [m^2]
	REAL(KIND=8)	:: TCS_X_Ar		! [m^2]
	REAL(KIND=8)	:: TCS_X_He		! [m^2]
	REAL(KIND=8)	:: TCS_X_H		! [m^2]

	!! Internal
	REAL(KIND=8)	:: EC_CO2, EC_CO, EC_N2, EC_H2, EC_O
	REAL(KIND=8)	:: EC_Ar, EC_He, EC_H, MassP

	IF ( Proj .EQ. 'H ') MassP = 1.0D0
	IF ( Proj .EQ. 'He') MassP = 4.0D0

  MT     = M_CO2
  EC_CO2 = E*MT/(MP+MT)
  MT     = M_CO
  EC_CO  = E*MT/(MP+MT)
  MT     = M_N2
  EC_N2  = E*MT/(MP+MT)
  MT     = M_H2
  EC_H2  = E*MT/(MP+MT)
  MT     = M_O
  EC_O   = E*MT/(MP+MT)
  MT     = M_Ar
  EC_Ar  = E*MT/(MP+MT)
  MT     = M_He
  EC_He  = E*MT/(MP+MT)
  MT     = M_H
  EC_H   = E*MT/(MP+MT)
	
	CALL Hp_H_tcs(EC_H,   TCS_X_H)
	CALL Hp_H_tcs(EC_He,  TCS_X_He)
	CALL Hp_H_tcs(EC_O,   TCS_X_O)
	CALL Hp_H_tcs(EC_Ar,  TCS_X_Ar)
	CALL Hp_H_tcs(EC_H2,  TCS_X_H2)
	CALL Hp_H_tcs(EC_N2,  TCS_X_N2)
	CALL Hp_H_tcs(EC_CO,  TCS_X_CO)
	CALL Hp_H_tcs(EC_CO2, TCS_X_CO2)

END SUBROUTINE get_ion_atom_tcs

!######################
!######################

SUBROUTINE write_qm_tcs
  USE physics_constants, ONLY : BOHRTOM
  USE planet, ONLY : PROJ

	IMPLICIT NONE

	REAL(KIND=8)	:: Ei, Ef, dE, E, T_H, T_He, T_O, T_Ar, T_H2, T_N2, T_CO, T_CO2
	INTEGER				:: NE, i

	Ei = 1.0D0
	Ef = 1.0D4
	NE = 1000
	dE = (Ef-Ei)/REAL(NE-1)

	OPEN(UNIT=99,FILE='../Data/QM_Total_Cross_Sections.dat',ACCESS='APPEND')	

	IF ( PROJ .EQ. 'H ' ) THEN
		DO i=1,NE
			E = Ei + REAL(i-1)*dE
  		CALL TCS_HH( E, T_H )
 	   	CALL TCS_HeH( E, T_He )
 	   	CALL universal_table_tcs( E, 'O  ', T_O   )
 	   	CALL universal_table_tcs( E, 'Ar ', T_Ar  )
 	   	CALL universal_table_tcs( E, 'H2 ', T_H2  )
 	   	CALL universal_table_tcs( E, 'N2 ', T_N2  )
 	   	CALL universal_table_tcs( E, 'CO ', T_CO  )
 	   	CALL universal_table_tcs( E, 'CO2', T_CO2 )
  		T_H   = T_H*BOHRTOM*BOHRTOM
  		T_He  = T_He*BOHRTOM*BOHRTOM
  		T_O   = T_O*BOHRTOM*BOHRTOM
  		T_Ar  = T_Ar*BOHRTOM*BOHRTOM
  		T_CO  = T_CO*BOHRTOM*BOHRTOM
  		T_CO2 = T_CO2*BOHRTOM*BOHRTOM
  		T_H2  = T_H2*BOHRTOM*BOHRTOM
  		T_N2  = T_N2*BOHRTOM*BOHRTOM
			WRITE(99,'(9ES10.2)') E, T_H, T_He, T_O, T_Ar, T_H2, T_N2, T_CO, T_CO2
		END DO
	ELSE IF ( PROJ .EQ. 'He' ) THEN
		DO i=1,NE
			E = Ei + REAL(i-1)*dE
    	CALL TCS_HeH( E, T_H )
    	CALL TCS_HeHe( E, T_He )
    	CALL TCS_HeO( E, T_O )
    	CALL universal_table_tcs( E, 'Ar ', T_Ar  )
    	CALL universal_table_tcs( E, 'H2 ', T_H2  )
    	CALL universal_table_tcs( E, 'N2 ', T_N2  )
    	CALL universal_table_tcs( E, 'CO ', T_CO  )
    	CALL universal_table_tcs( E, 'CO2', T_CO2 )
  		T_H   = T_H*BOHRTOM*BOHRTOM
  		T_He  = T_He*BOHRTOM*BOHRTOM
  		T_O   = T_O*BOHRTOM*BOHRTOM
  		T_Ar  = T_Ar*BOHRTOM*BOHRTOM
  		T_CO  = T_CO*BOHRTOM*BOHRTOM
  		T_CO2 = T_CO2*BOHRTOM*BOHRTOM
  		T_H2  = T_H2*BOHRTOM*BOHRTOM
  		T_N2  = T_N2*BOHRTOM*BOHRTOM
			WRITE(99,'(9ES10.2)') E, T_H, T_He, T_O, T_Ar, T_H2, T_N2, T_CO, T_CO2
		END DO
	END IF

	CLOSE(99)

END SUBROUTINE write_qm_tcs

!######################
!######################

SUBROUTINE get_qm_tcs(E_now,TCS_X_CO2,TCS_X_CO,TCS_X_N2,TCS_X_H2,TCS_X_O,TCS_X_Ar,TCS_X_He,TCS_X_H)
	use physics_constants, only : bohrtom
	use planet

	implicit none

	!! inputs
	real(kind=8)	:: e_now 			! lab frame [ev]
	
	!! outputs
	real(kind=8)	:: tcs_x_co2	! [m^2]
	real(kind=8)	:: tcs_x_co		! [m^2]
	real(kind=8)	:: tcs_x_n2		! [m^2]
	real(kind=8)	:: tcs_x_h2		! [m^2]
	real(kind=8)	:: tcs_x_o		! [m^2]
	real(kind=8)	:: tcs_x_ar		! [m^2]
	real(kind=8)	:: tcs_x_he		! [m^2]
	real(kind=8)	:: tcs_x_h		! [m^2]

	!! internal
	real(kind=8)	:: ec_co2, ec_co, ec_n2, ec_h2, ec_o
	real(kind=8)	:: ec_ar, ec_he, ec_h

	!! dummy targ getting the proj
	IF (Targ .EQ. '') Targ = 'H  '
	call mass_finder

	!! convert energy to CM frame
	MT     = M_CO2	
	EC_CO2 = E_now*MT/(MP+MT)
	MT     = M_CO	
	EC_CO  = E_now*MT/(MP+MT)
	MT     = M_N2	
	EC_N2  = E_now*MT/(MP+MT)
	MT     = M_H2	
	EC_H2  = E_now*MT/(MP+MT)
	MT     = M_O	
	EC_O   = E_now*MT/(MP+MT)
	MT     = M_Ar	
	EC_Ar  = E_now*MT/(MP+MT)
	MT     = M_He	
	EC_He  = E_now*MT/(MP+MT)
	MT     = M_H	
	EC_H   = E_now*MT/(MP+MT)

	IF ( TRIM(Proj) .EQ. 'H') THEN
		CALL TCS_HH( EC_H, TCS_X_H )
		CALL TCS_HeH( EC_He, TCS_X_He )
		CALL universal_table_tcs( EC_O,  'O  ', TCS_X_O   )
		CALL universal_table_tcs( EC_Ar, 'Ar ', TCS_X_Ar  )
		CALL universal_table_tcs( EC_H2, 'H2 ', TCS_X_H2  )
		CALL universal_table_tcs( EC_N2, 'N2 ', TCS_X_N2  )
		CALL universal_table_tcs( EC_CO, 'CO ', TCS_X_CO  )
		CALL universal_table_tcs( EC_CO2,'CO2', TCS_X_CO2 )
	ELSE IF ( TRIM(Proj) .EQ. 'He' ) THEN
		CALL TCS_HeH( EC_H, TCS_X_H )
		CALL TCS_HeHe( EC_He, TCS_X_He )
		CALL TCS_HeO( EC_O, TCS_X_O )
		CALL universal_table_tcs( EC_Ar, 'Ar ', TCS_X_Ar  )
		CALL universal_table_tcs( EC_H2, 'H2 ', TCS_X_H2  )
		CALL universal_table_tcs( EC_N2, 'N2 ', TCS_X_N2  )
		CALL universal_table_tcs( EC_CO, 'CO ', TCS_X_CO  )
		CALL universal_table_tcs( EC_CO2,'CO2', TCS_X_CO2 )
	ELSE IF ( TRIM(Proj) .EQ. 'O' ) THEN
		CALL universal_table_tcs( EC_H, 'H  ', TCS_X_H  )
		CALL TCS_HeO( EC_He, TCS_X_He )
		CALL universal_table_tcs( EC_Ar, 'Ar ', TCS_X_Ar  )
		CALL universal_table_tcs( EC_O,  'O  ', TCS_X_O   )
		CALL universal_table_tcs( EC_H2, 'H2 ', TCS_X_H2  )
		CALL universal_table_tcs( EC_N2, 'N2 ', TCS_X_N2  )
		CALL universal_table_tcs( EC_CO, 'CO ', TCS_X_CO  )
		CALL universal_table_tcs( EC_CO2,'CO2', TCS_X_CO2 )
	ELSE IF ( TRIM(Proj) .EQ. 'H2O' ) THEN
		CALL TCS_HeH( EC_H, TCS_X_H )
		CALL TCS_HeHe( EC_He, TCS_X_He )
		CALL TCS_HeO( EC_O, TCS_X_O )
		CALL universal_table_tcs( EC_Ar, 'Ar ', TCS_X_Ar  )
		CALL universal_table_tcs( EC_H2, 'H2 ', TCS_X_H2  )
		CALL universal_table_tcs( EC_N2, 'N2 ', TCS_X_N2  )
		CALL universal_table_tcs( EC_CO, 'CO ', TCS_X_CO  )
		CALL universal_table_tcs( EC_CO2,'CO2', TCS_X_CO2 )
	ELSE
		CALL universal_tcs_int( EC_H, MP, M_H, 	1, TCS_X_H )
		CALL universal_tcs_int( EC_He, MP, M_He, 	1, TCS_X_He )
		CALL universal_tcs_int( EC_Ar, MP, M_Ar, 	1, TCS_X_Ar )
		CALL universal_tcs_int( EC_O, MP, M_O, 	1, TCS_X_O )
		CALL universal_tcs_int( EC_CO2, MP, M_CO2, 2, TCS_X_CO2 )
		CALL universal_tcs_int( EC_CO, MP, M_CO, 	2, TCS_X_CO )
		CALL universal_tcs_int( EC_N2, MP, M_N2, 	2, TCS_X_N2 )
		CALL universal_tcs_int( EC_H2, MP, M_H2, 	2, TCS_X_H2 )
	END IF	
	TCS_X_H   = TCS_X_H*BOHRTOM*BOHRTOM
	TCS_X_He  = TCS_X_He*BOHRTOM*BOHRTOM
	TCS_X_O   = TCS_X_O*BOHRTOM*BOHRTOM
	TCS_X_Ar  = TCS_X_Ar*BOHRTOM*BOHRTOM
	TCS_X_CO  = TCS_X_CO*BOHRTOM*BOHRTOM
	TCS_X_CO2 = TCS_X_CO2*BOHRTOM*BOHRTOM
	TCS_X_H2  = TCS_X_H2*BOHRTOM*BOHRTOM
	TCS_X_N2  = TCS_X_N2*BOHRTOM*BOHRTOM

END SUBROUTINE get_qm_tcs

!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!

SUBROUTINE get_qm_tcs_venus(E_now,TCS_X_H,TCS_X_N,TCS_X_O,TCS_X_S,TCS_X_O2,TCS_X_CO,TCS_X_CO2,TCS_X_SO,TCS_X_SO2)
	use physics_constants, only : bohrtom
	use planet

	implicit none

	!! inputs
	real(kind=8)	:: e_now 			! lab frame [ev]
	
	!! outputs
	real(kind=8)	:: tcs_x_co2	! [m^2]
	real(kind=8)	:: tcs_x_co		! [m^2]
	real(kind=8)	:: tcs_x_n		! [m^2]
	real(kind=8)	:: tcs_x_h		! [m^2]
	real(kind=8)	:: tcs_x_o		! [m^2]
	real(kind=8)	:: tcs_x_s		! [m^2]
	real(kind=8)	:: tcs_x_O2		! [m^2]
	real(kind=8)	:: tcs_x_so		! [m^2]
	real(kind=8)	:: tcs_x_so2	! [m^2]

	!! internal
	real(kind=8)	:: ec_co2, ec_co, ec_n, ec_h, ec_o
	real(kind=8)	:: ec_S, ec_o2, ec_so, ec_so2

	!! dummy targ getting the proj
!	targ = 'h  '
	call mass_finder

	!! convert energy to CM frame
	MT     = M_H
	EC_H   = E_now*MT/(MP+MT)
	MT     = M_CO2	
	EC_CO2 = E_now*MT/(MP+MT)
	MT     = M_CO	
	EC_CO  = E_now*MT/(MP+MT)
	MT     = M_N	
	EC_N   = E_now*MT/(MP+MT)
	MT     = M_O	
	EC_O   = E_now*MT/(MP+MT)
	MT     = M_O2	
	EC_O2  = E_now*MT/(MP+MT)
	MT     = M_S
	EC_S   = E_now*MT/(MP+MT)
	MT     = M_SO	
	EC_SO  = E_now*MT/(MP+MT)
	MT     = M_SO2
	EC_SO2 = E_now*MT/(MP+MT)

	IF ( TRIM(Proj) .EQ. 'H') THEN
		CALL TCS_HH( EC_H, TCS_X_H )
		CALL universal_table_tcs( EC_O,  'O  ', TCS_X_O   )
		CALL universal_table_tcs( EC_H, 'H  ', TCS_X_H  )
		CALL universal_table_tcs( EC_N, 'N  ', TCS_X_N  )
		CALL universal_table_tcs( EC_S, 'S  ', TCS_X_S  )
		CALL universal_table_tcs( EC_CO, 'CO ', TCS_X_CO  )
		CALL universal_table_tcs( EC_CO2,'CO2', TCS_X_CO2 )
		CALL universal_table_tcs( EC_SO, 'SO ', TCS_X_SO  )
		CALL universal_table_tcs( EC_SO2,'SO2', TCS_X_SO2 )
		CALL universal_table_tcs( EC_O2,'O2 ', TCS_X_O2 )
	ELSE IF ( TRIM(Proj) .EQ. 'He' ) THEN
		CALL TCS_HH( EC_H, TCS_X_H )
		CALL universal_table_tcs( EC_H, 'H  ', TCS_X_H  )
		CALL universal_table_tcs( EC_N, 'N  ', TCS_X_N  )
		CALL universal_table_tcs( EC_S, 'S  ', TCS_X_S  )
		CALL universal_table_tcs( EC_CO, 'CO ', TCS_X_CO  )
		CALL universal_table_tcs( EC_CO2,'CO2', TCS_X_CO2 )
		CALL universal_table_tcs( EC_SO, 'SO ', TCS_X_SO  )
		CALL universal_table_tcs( EC_SO2,'SO2', TCS_X_SO2 )
		CALL universal_table_tcs( EC_O2,'O2 ', TCS_X_O2 )
		CALL TCS_HeO( EC_O, TCS_X_O )
	ELSE
		CALL universal_tcs_int( EC_H, MP, M_H, 	1, TCS_X_H )
		CALL universal_tcs_int( EC_N, MP, M_N, 	1, TCS_X_N )
		CALL universal_tcs_int( EC_O, MP, M_O, 	1, TCS_X_O )
		CALL universal_tcs_int( EC_S, MP, M_S, 	1, TCS_X_S )
		CALL universal_tcs_int( EC_CO, MP, M_CO, 	1, TCS_X_CO )
		CALL universal_tcs_int( EC_CO2, MP, M_CO2, 2, TCS_X_CO2 )
		CALL universal_tcs_int( EC_SO, MP, M_SO, 	2, TCS_X_SO )
		CALL universal_tcs_int( EC_SO2, MP, M_SO2, 	2, TCS_X_SO2 )
		CALL universal_tcs_int( EC_O2, MP, M_O2, 	2, TCS_X_O2 )
	END IF	
	TCS_X_H   = TCS_X_H*BOHRTOM*BOHRTOM
	TCS_X_N   = TCS_X_N*BOHRTOM*BOHRTOM
	TCS_X_O   = TCS_X_O*BOHRTOM*BOHRTOM
	TCS_X_S   = TCS_X_S*BOHRTOM*BOHRTOM
	TCS_X_O2  = TCS_X_O2*BOHRTOM*BOHRTOM
	TCS_X_CO  = TCS_X_CO*BOHRTOM*BOHRTOM
	TCS_X_CO2 = TCS_X_CO2*BOHRTOM*BOHRTOM
	TCS_X_SO  = TCS_X_SO*BOHRTOM*BOHRTOM
	TCS_X_SO2 = TCS_X_SO2*BOHRTOM*BOHRTOM

END SUBROUTINE get_qm_tcs_venus

!!
!!

SUBROUTINE get_hs_tcs(TCS_X_CO2,TCS_X_CO,TCS_X_N2,TCS_X_H2,TCS_X_O,TCS_X_Ar,TCS_X_He,TCS_X_H)
  USE physics_constants, ONLY : BOHRTOM
	USE planet

  IMPLICIT NONE

  !! Outputs
	REAL(KIND=8)	:: TCS_X_CO2	! [m^2]
	REAL(KIND=8)	:: TCS_X_CO		! [m^2]
	REAL(KIND=8)	:: TCS_X_N2		! [m^2]
	REAL(KIND=8)	:: TCS_X_H2		! [m^2]
	REAL(KIND=8)	:: TCS_X_O		! [m^2]
	REAL(KIND=8)	:: TCS_X_Ar		! [m^2]
	REAL(KIND=8)	:: TCS_X_He		! [m^2]
	REAL(KIND=8)	:: TCS_X_H		! [m^2]

	!! Internal
	REAL(KIND=8)	:: RP

	IF ( Proj .EQ. 'H  ' ) THEN
		RP = HS_R_H
	ELSE IF ( Proj .EQ. 'He ' ) THEN
		RP = HS_R_He
	ELSE IF ( Proj .EQ. 'O  ' ) THEN
		RP = HS_R_O
	ELSE IF ( Proj .EQ. 'Ar ' ) THEN
		RP = HS_R_Ar
	ELSE IF ( Proj .EQ. 'H2 ' ) THEN
		RP = HS_R_H2
	ELSE IF ( Proj .EQ. 'N2 ' ) THEN
		RP = HS_R_N2
	ELSE IF ( Proj .EQ. 'CO ' ) THEN
		RP = HS_R_CO
	ELSE IF ( Proj .EQ. 'CO2' ) THEN
		RP = HS_R_CO2
	ELSE 
		WRITE(*,*) 'Projectile: ', Proj, ' Not found in HS TCS '
	END IF	

	CALL HS_tcs( RP, HS_R_CO2, TCS_X_CO2 )
	CALL HS_tcs( RP, HS_R_CO , TCS_X_CO  )
	CALL HS_tcs( RP, HS_R_N2 , TCS_X_N2  )
	CALL HS_tcs( RP, HS_R_H2 , TCS_X_H2  )
	CALL HS_tcs( RP, HS_R_O  , TCS_X_O   )
	CALL HS_tcs( RP, HS_R_He , TCS_X_He  )
	CALL HS_tcs( RP, HS_R_Ar , TCS_X_Ar  )
	CALL HS_tcs( RP, HS_R_H  , TCS_X_H   )

END SUBROUTINE get_hs_tcs

!!!!!!!!!!!!!!
!!!!!!!!!!!!!!

SUBROUTINE get_hs_tcs_venus(TCS_X_H,TCS_X_N,TCS_X_O,TCS_X_S,TCS_X_O2,TCS_X_CO,TCS_X_CO2,TCS_X_SO,TCS_X_SO2)
  USE physics_constants, ONLY : BOHRTOM
	USE planet

  IMPLICIT NONE

  !! Outputs
	REAL(KIND=8)	:: TCS_X_CO2	! [m^2]
	REAL(KIND=8)	:: TCS_X_CO		! [m^2]
	REAL(KIND=8)	:: TCS_X_O2		! [m^2]
	REAL(KIND=8)	:: TCS_X_SO		! [m^2]
	REAL(KIND=8)	:: TCS_X_O		! [m^2]
	REAL(KIND=8)	:: TCS_X_SO2  ! [m^2]
	REAL(KIND=8)	:: TCS_X_S		! [m^2]
	REAL(KIND=8)	:: TCS_X_H		! [m^2]
	REAL(KIND=8)	:: TCS_X_N		! [m^2]

	!! Internal
	REAL(KIND=8)	:: RP

	IF ( Proj .EQ. 'H  ' ) THEN
		RP = HS_R_H
	ELSE IF ( Proj .EQ. 'He ' ) THEN
		RP = HS_R_He
	ELSE IF ( Proj .EQ. 'O  ' ) THEN
		RP = HS_R_O
	ELSE IF ( Proj .EQ. 'Ar ' ) THEN
		RP = HS_R_Ar
	ELSE IF ( Proj .EQ. 'H2 ' ) THEN
		RP = HS_R_H2
	ELSE IF ( Proj .EQ. 'N2 ' ) THEN
		RP = HS_R_N2
	ELSE IF ( Proj .EQ. 'N  ' ) THEN
		RP = HS_R_N
	ELSE IF ( Proj .EQ. 'CO ' ) THEN
		RP = HS_R_CO
	ELSE IF ( Proj .EQ. 'CO2' ) THEN
		RP = HS_R_CO2
	ELSE IF ( Proj .EQ. 'SO2' ) THEN
		RP = HS_R_SO2
	ELSE IF ( Proj .EQ. 'SO ' ) THEN
		RP = HS_R_SO
	ELSE IF ( Proj .EQ. 'S  ' ) THEN
		RP = HS_R_S
	ELSE IF ( Proj .EQ. 'O2 ' ) THEN
		RP = HS_R_O2
	ELSE 
		WRITE(*,*) 'Projectile: ', Proj, ' Not found in HS TCS '
	END IF	

	CALL HS_tcs( RP, HS_R_CO2, TCS_X_CO2 )
	CALL HS_tcs( RP, HS_R_CO , TCS_X_CO  )
	CALL HS_tcs( RP, HS_R_N , TCS_X_N  )
	CALL HS_tcs( RP, HS_R_H , TCS_X_H  )
	CALL HS_tcs( RP, HS_R_O  , TCS_X_O   )
	CALL HS_tcs( RP, HS_R_O2  , TCS_X_O2   )
	CALL HS_tcs( RP, HS_R_S , TCS_X_S  )
	CALL HS_tcs( RP, HS_R_SO , TCS_X_SO  )
	CALL HS_tcs( RP, HS_R_SO2  , TCS_X_SO2   )

END SUBROUTINE get_hs_tcs_venus


