
SUBROUTINE  lism_onestep_transport( E0, x0, y0, z0, ux0, uy0, uz0, &
      &                            E1, x1, y1, z1, ux1, uy1, uz1, dt )
	USE lism
	USE planet, ONLY : M_H, M_He
	USE physics_constants
	USE elastic_atom_ion
	USE mpi_info
	
	IMPLICIT NONE

	!! Inputs
	REAL(KIND=8)	:: E0, x0, y0, z0, ux0, uy0, uz0

	!! Outputs
	REAL(KIND=8)	:: E1, x1, y1, z1, ux1, uy1, uz1, dt

	!! Internal
	REAL(KIND=8)	:: lfg
	REAL(KIND=8)	:: Den_now, E_temp, r0, TCS, ion_Den_now, CX_CS, theta_now
	REAL(KIND=8)	:: Den_tot, Prob_tot, ds, MFP, phi, theta, phi_now
	REAL(KIND=8)	:: Prob_atom, Prob_ion, r1, ScattAng, ion_path, Lab_ScattAng, tt
	REAL(KIND=8)	:: atom_den, C1, ddt, Length, velocity, r2, rat
	REAL(KIND=8)	:: x_now(3), u_now(3), x_nxt(3), u_nxt(3)
	REAL(KIND=8)	:: He_neu_Den_now, He_ion_Den_now, H_neu_Den_now, H_ion_Den_now
	REAL(KIND=8)	:: TCS_EL_HeH, TCS_EL_HeHe, TCS_EL_HeHp, TCS_EL_HeHep
	REAL(KIND=8)	:: TCS_IE_HeHp, TCS_IE_HeHep, Prob_ion_el, Prob_ion_ie
	REAL(KIND=8)	:: P_EL_HeH, P_EL_HeHe, P_EL_HeHp, P_EL_HeHep
	REAL(KIND=8)	:: P_IE_HeHp, P_IE_HeHep
	REAL(KIND=8)	:: MR_H, MR_He, MR_Hp, MR_Hep
	REAL(KIND=8)	:: pp1, pp2, pp3, pp4, pp5, pp6, pp_tot
	REAL(KIND=8)	:: pp_prob(6)
	INTEGER				:: i	

	x_now(1) = x0
	x_now(2) = y0
	x_now(3) = z0
	u_now(1) = ux0
	u_now(2) = uy0
	u_now(3) = uz0
	r0       = SQRT(x0**2+y0**2+z0**2)

	!! Get density [1/m^3] and temperature [eV] in current position
  CALL LISM_density( r0, H_neu_Den_now )
  CALL LISM_density_He( r0, He_neu_Den_now )
  CALL LISM_ion_density( r0, H_ion_Den_now )
  CALL LISM_ion_density_He( r0, He_ion_Den_now )
  CALL LISM_temp( r0, E_temp )

!	WRITE(*,'(A,3ES10.2)') 'TRANSPORT R: den: temp: ', r0, Den_now, E_temp

  !! Find atom-atom TCS for current energy
  IF (TRIM(Proj) .EQ. 'He') THEN
  	IF (QM_ON .EQ. 1) THEN
    	IF (E0 .GE. 1.0D4) CALL TCS_HeH( 9.99D3, TCS_EL_HeH )
    	IF (E0 .LT. 1.0D4) CALL TCS_HeH( E0, TCS_EL_HeH )
			CALL TCS_HeHe( E0, TCS_EL_HeHe )
			CALL elastic_He_Hp_tcs( E0, TCS_EL_HeHp )
			CALL elastic_He_Hep_tcs( E0, TCS_EL_HeHep )
			TCS_EL_HeH  = TCS_EL_HeH*BOHRTOM*BOHRTOM
			TCS_EL_HeHe = TCS_EL_HeHe*BOHRTOM*BOHRTOM
			CALL cx_cross_sections( 'H_p  ', 'He ', 1, E0, TCS_IE_HeHp )
			CALL cx_cross_sections( 'He_p ', 'He ', 1, E0, TCS_IE_HeHep )
    ELSE IF (QM_ON .EQ. 0) THEN
    	CALL HS_TCS( M_He, M_H, TCS )
    END IF
  ELSE IF (TRIM(Proj) .EQ. 'H') THEN
  	IF (QM_ON .EQ. 1) THEN
    	IF (E0 .GE. 1.0D4) CALL TCS_HH( 1.0D4, TCS )
      IF (E0 .LT. 1.0D4) CALL TCS_HH( E0, TCS )
    ELSE IF (QM_ON .EQ. 0) THEN
    	CALL HS_TCS( M_H, M_He, TCS )
    END IF
  END IF

	P_EL_HeH   = H_neu_Den_now*TCS_EL_HeH
	P_EL_HeHe  = He_neu_Den_now*TCS_EL_HeHe
	P_EL_HeHp  = H_ion_Den_now*TCS_EL_HeHp
	P_EL_HeHep = He_ion_Den_now*TCS_EL_HeHep
	P_IE_HeHp  = H_ion_Den_now*TCS_IE_HeHp
	P_IE_HeHep = He_ion_Den_now*TCS_IE_HeHep

  !! Get CX cross sections for proj (H or He) plus H+      
	CALL cx_cross_sections( 'H_p  ', Proj, 1, E0, CX_CS )

  !! Get mixing ratio for neutrals and ions
	Den_tot  = He_neu_Den_now + He_ion_Den_now + H_neu_Den_now + H_ion_Den_now 
	MR_H     = H_neu_Den_now/Den_tot
	MR_He    = He_neu_Den_now/Den_tot
	MR_Hp    = H_ion_Den_now/Den_tot
	MR_Hep   = He_ion_Den_now/Den_tot

	!! total collisional probability including elastic and inelastic channels
	Prob_tot = P_EL_HeH + P_EL_HeHe + P_EL_HeHp + P_EL_HeHep + P_IE_HeHp + P_IE_HeHep

  MFP      = MTOAU/Prob_tot ! [AU]
  ds       = 0.01D0*MFP			! [AU]
	C1       = EXP(-ds/MFP)

  !! Determine if collision is atom-atom or atom-ion
  Prob_atom    = (P_EL_HeH + P_EL_HeHe)/Prob_tot
  Prob_ion_el  = (P_EL_HeHp + P_EL_HeHep)/Prob_tot
	Prob_ion_ie  = (P_IE_HeHp + P_IE_HeHep)/Prob_tot

  !! get random number to deterimine if collision is with atom or ion
  r1 = lfg()
	

	IF (r1 .GT. C1) THEN
 		!#######################################
		!! Collision occurs
 		!#######################################

		CALL energy_to_velocity( E0, MP, velocity )
		Length 	= -MFP*LOG(r1) 						! [m]
		ddt 		= Length/(MTOAU*Velocity) ! [sec]
		r2  		= lfg()

		pp1     = MR_H*TCS_EL_HeH
		pp2			= MR_He*TCS_EL_HeHe
		pp3			= MR_Hp*TCS_EL_HeHp
		pp4			= MR_Hep*TCS_EL_HeHep
		pp5			= MR_Hp*TCS_IE_HeHp
		pp6			= MR_Hep*TCS_IE_HeHep
		pp_tot  = pp1 + pp2 + pp3 + pp4 + pp5 + pp6

		pp_prob(1) = pp1/pp_tot
		pp_prob(2) = pp_prob(1) + pp2/pp_tot		
		pp_prob(3) = pp_prob(2) + pp3/pp_tot		
		pp_prob(4) = pp_prob(3) + pp4/pp_tot		
		pp_prob(5) = pp_prob(4) + pp5/pp_tot		
		pp_prob(6) = pp_prob(5) + pp6/pp_tot		

		IF (myid .eq. 0) WRITE(44,*) E0, pp_prob(1), pp_prob(2), pp_prob(3), pp_prob(4), pp_prob(5), pp_prob(6)

  	IF ( r2 .LT. pp_prob(4) ) THEN
	  	!#######################################
 	   	!## collision is elastic atom-atom #####
 	   	!## or elastic atom-ion            #####
    	!#######################################

    	!! Find random theta scattering angle
    	IF (TRIM(Proj) .EQ. 'H') THEN
    		IF (E0 .GT. 2.5D3) THEN
      		IF (QM_ON .EQ. 0) CALL HS_rand_angle( ScattAng )
        	IF (QM_ON .EQ. 1) CALL lin_rand_angle( 2.5D3, 'HH    ', ScattAng )
      	ELSE
        	IF (QM_ON .EQ. 0) CALL HS_rand_angle( ScattAng )
        	IF (QM_ON .EQ. 1) CALL lin_rand_angle( E0, 'HH    ', ScattAng )
      	END IF
				rat = 1.0D0	! mass ratio MP/MT for energy-angle conversions
			!! He projectile
    	ELSE
				!! He + H elastic
				IF ( r2 .LE. pp_prob(1) ) THEN
					CALL lin_rand_angle( E0, 'HeH   ', ScattAng )
					if (myid .eq. 0) WRITE(55,*) E0, ScattAng, 1
				!! He + He elastic
				ELSE IF ( (r2 .GT. pp_prob(1)) .and. (r2 .LE. pp_prob(2)) ) THEN
					CALL lin_rand_angle( E0, 'HeHe  ', ScattAng )
					if (myid .eq. 0) WRITE(55,*) E0, ScattAng, 2
				!! He + Hp elastic
				ELSE IF ( (r2 .GT. pp_prob(2)) .and. (r2 .LE. pp_prob(3)) ) THEN
					CALL elastic_He_Hp_rand_angle( E0, ScattAng )
					if (myid .eq. 0) WRITE(55,*) E0, ScattAng, 3
				!! He + Hep elastic
				ELSE IF ( (r2 .GT. pp_prob(3)) .and. (r2 .LE. pp_prob(4)) ) THEN
					CALL elastic_He_Hep_rand_angle( E0, ScattAng )
					if (myid .eq. 0) WRITE(55,*) E0, ScattAng, 4
				END IF

				rat = 4.0D0	! mass ratio MP/MT for energy-angle conversions
    	END IF

    	!! Find new energy after collision at angle ScattAng
    	CALL find_new_energy( E0, ScattAng, rat, E1 )

    	!! Convert ScattAng to lab frame
    	CALL angle_to_lab( ScattAng, rat, Lab_ScattAng )

    	!! Find random phi scattering angle
    	r1     = lfg()
    	phi    = r1*2.0D0*PI

    	!! Update new scattering angles
    	theta_now = Lab_ScattAng
    	phi_now   = phi

    	!! Transport particle and get new unit velocity and position vectors
			CALL particle_dir_cos_transport(u_now,x_now,theta_now,phi_now,Length,u_nxt,x_nxt)

			x1 	= x_nxt(1)	
			y1 	= x_nxt(2)	
			z1 	= x_nxt(3)	
			ux1	= u_nxt(1)
			uy1	= u_nxt(2)
			uz1	= u_nxt(3)
			dt  = ddt

  	ELSE
	  	!#######################################
    	!## collision is atom-ion CX ###########
    	!#######################################

    	!! Update position location to CX position before transporting ion
    	x_now(:) = x_now(:) + u_now(:)*Length
			r0       = SQRT(x_now(1)**2+x_now(2)**2+x_now(3)**2)
			CALL LISM_density( r0, atom_den )

			if (myid .eq. 0) WRITE(66,*) E0

			IF (atom_den .NE. 0.0D0) THEN
	    	!! transport newly created ion 
 		   	CALL ion_transport(MP,1,E0,atom_den,CX_CS,x_now/MTOAU,u_now,x_nxt,u_nxt,tt,ion_path)
				dt = ddt + tt
	    	!! convert x_nxt from [m] back to [AU]
 		   	x_nxt = x_nxt * MTOAU
				x1 		= x_nxt(1)	
				y1 		= x_nxt(2)	
				z1 		= x_nxt(3)	
				ux1		= u_nxt(1)
				uy1		= u_nxt(2)
				uz1		= u_nxt(3)

	    	!! find energy loss due to ion transport
 		   	IF ( TRIM(Proj) .EQ. 'H') THEN
 		   		IF ( ION_METH .EQ. 1) THEN
						CALL butler_drag( 1, ion_path, ion_Den_now, E_temp*11604.505D0, E0, E1 )
      		ELSE IF (ION_METH .EQ. 2) THEN
						CALL bethe_drag( 1, ion_path, ion_Den_now, E0, E1 )
					END IF
    		ELSE IF ( TRIM(Proj) .EQ. 'He') THEN
	    		IF ( ION_METH .EQ. 1) THEN
						CALL butler_drag( 2, ion_path, ion_Den_now, E_temp*11604.505D0, E0, E1 )
      		ELSE IF ( ION_METH .EQ. 2) THEN
						CALL bethe_drag( 2, ion_path, ion_Den_now, E0, E1 )
					END IF ! ION_METH
    		END IF ! Proj type
     		IF (E1 .LT. 0.01D0) E1 = 0.0D0
			ELSE IF (atom_den .EQ. 0.0D0) THEN 
				dt 		= ddt
				E1 		= 0.0D0
				u_nxt = u_now
				x_nxt = x_now
				x1 		= x_nxt(1)	
				y1 		= x_nxt(2)	
				z1 		= x_nxt(3)	
				ux1		= u_nxt(1)
				uy1		= u_nxt(2)
				uz1		= u_nxt(3)
			END IF ! atom_den is zero or not

		END IF  ! collision is atom-ion CX
	ELSE ! no collision occurs
	 	!#######################################
   	!######## NO COLLISION #################
   	!#######################################
!		WRITE(*,*) 'rank: NO collision', myid

		x_nxt = x_now + u_now*ds
		u_nxt = u_now	
		E1    = E0
    CALL energy_to_velocity( E1, MP, velocity )
    dt    = ds/(MTOAU*Velocity) ! [sec]
		x1  	= x_nxt(1)	
		y1  	= x_nxt(2)	
		z1  	= x_nxt(3)	
		ux1 	= u_nxt(1)
		uy1 	= u_nxt(2)
		uz1 	= u_nxt(3)


!		WRITE(*,'(A,6ES10.1)') 'NC: x: u: ', x1, y1, z1, ux1, uy1, uz1

	END IF ! collision occurs or not


END SUBROUTINE lism_onestep_transport

