
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! ion transport routine
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE ion_onestep_transport( E0, R0, U0, E1, R1, U1, tt )
	USE ENA
	USE physics_constants, ONLY : PI
	USE mpi_info 

	IMPLICIT NONE

	INCLUDE 'mpif.h'

	!! Inputs
	REAL(KIND=8)			:: E0				! energy of ion [eV]
	REAL(KIND=8)			:: R0(3)		! initial position vector [x y z] [AU]
	REAL(KIND=8)			:: U0(3)		! initial directional cosines [ux uy uz]

	!! Outputs
	REAL(KIND=8)			:: E1				! energy of ion [eV]
	REAL(KIND=8)			:: R1(3) 		! final position vector [x y z] [m]
	REAL(KIND=8)			:: U1(3) 		! final unit velocity vector [vx vy vz]
	REAL(KIND=8)			:: tt				! time during transport

	!! Internal
	INTEGER						:: GOOD, NumL, i, j
	REAL(KIND=8)			:: M, M_p, q, JeV, Bt, Bp, B(3), eps, dt, L, lam, Mvout
	REAL(KIND=8)			:: v_perp, v_parr, VV0(3), v, pathL, d, RR, r, theta
	REAL(KIND=8)			:: AU, Rad, w0, v0_perp(3), v0_parr(3), C(3), n1(3), n2(3), MC
	REAL(KIND=8)			:: RC(3), rnow(3), rnxt(3), dL, t, P, nowL, rt(3), rp(3)
	REAL(KIND=8)			:: lfg
	REAL(KIND=8)			:: t_start, t_end, t_elap

	AU  = 1.5D11					! astronomical unit [m]
	M_p = 1.672D-27				! mass of proton [kg]
	q   = -1.602D-19			! charge of electron [C]
	JeV = 6.24150974D18		! Joules to eV  1 [J] = JeV [eV]

	!!! Get absolute distance from sun
	RR = SQRT(R0(1)**2 + R0(2)**2 + R0(3)**2)

	!!! Get magnetic field vector
	CALL LISM_MF( R0, B )

	!!! Get mass of ion [kg]
	M    = amu*M_p 

	!!! Percent of rad freq per time step
	eps  = 0.01D0

	!!! Time step
	dt   = eps*M/(Z*ABS(q)*B0)

	!!! Mean free path [m]	
	lam  = 1.0D0/(den*CX) 
	
	!!! Step length for transport [m]
	L    = 0.2D0*lam
	
	!!! Get random path length
	P    = EXP(-L/lam)
	GOOD = 0
	NumL = 0
	
!	WRITE(*,'(A,4ES10.2)') 'IT - Den: CX: MFP: Prob: ', den, CX, lam/AU, P

	DO WHILE (GOOD .EQ. 0)
		r = lfg()
		IF ( r .GE. P ) THEN 	!! CX collision back to neutral
			d    = -lam*LOG(r)
			GOOD = 1	
		ELSE 									!! No CX collision, keep transporting
			NumL = NumL + 1
		END IF
	END DO 	

	pathL = L*NumL + d	
	ion_path = pathL

!	WRITE(*,*) 'IT - MFP: Path: ', lam/AU, ion_path/AU

	!!! get initial velocity
	v   = SQRT(2.0D0*E/(M*JeV))		! magnitude of velocity [m/s]
	VV0 = V0*v										! initial velocity vector [m/s]

	!!! get angle between B and VV0
	theta = ACOS( (B(1)*VV0(1)+B(2)*VV0(2)+B(3)*VV0(3))/(v*B0) )
	
	!!! get perp and parr velocity vectors
	v_perp  = v*SIN(theta)
	v_parr  = v*COS(theta)	
	v0_perp = VV0 - ( (VV0(1)*B(1)+VV0(2)*B(2)+VV0(3)*B(3))/(B0*B0) )*B
	v0_parr =  ( (VV0(1)*B(1)+VV0(2)*B(2)+VV0(3)*B(3))/(B0*B0) )*B

	!!! get natural frequency and radius of oscillation
	w0      = Z*q*B0/M
	Rad     = M*v_perp/ABS(Z*q*B0)	
	
	!!! get unit vector for perp plane
	C(1)    = v0_perp(2)*B(3) - v0_perp(3)*B(2)
	C(2)    = v0_perp(3)*B(1) - v0_perp(1)*B(3)
	C(3)    = v0_perp(1)*B(2) - v0_perp(2)*B(1)
	MC      = SQRT( C(1)**2 + C(2)**2 + C(3)**2 )
	n1      = C/MC
	C(1)    = n1(2)*B(3) - n1(3)*B(2)
	C(2)    = n1(3)*B(1) - n1(1)*B(3)
	C(3)    = n1(1)*B(2) - n1(2)*B(1)
	MC      = SQRT( C(1)**2 + C(2)**2 + C(3)**2 )
	n2      = C/MC

	!!! get constant for position vector
	RC      = R0 - Rad*n1

	!!! get time as ion	
	t       = pathL/v

	!!! get position and velocity at time t	
	rp   	= v0_parr*t
	rt   	= Rad*COS(w0*t)*n1 + Rad*SIN(w0*t)*n2
	rout 	= RC + rp + rt
	vout 	= v0_parr + Rad*w0*( -SIN(w0*t)*n1 + COS(w0*t)*n2 )
	Mvout = SQRT( vout(1)**2 + vout(2)**2 + vout(3)**2 )	

	vout(1) = vout(1)/Mvout	
	vout(2) = vout(2)/Mvout	
	vout(3) = vout(3)/Mvout	
	tt      = t

END SUBROUTINE ion_onestep_transport

