
!! Fill all phase spaces which will
!! be written to files for data
!! analysis
!! ALL CLICKS AVERAGE

SUBROUTINE fill_all_lism_phase_spaces( MC )
	USE lism_click

	IMPLICIT NONE

	!! Input
	INTEGER					:: MC		! current MC particle

	!! Internal
	INTEGER					:: i, j, iX, iY
	REAL(KIND=8)		:: ddT, ddE, ddN, dEL, ddXYZ, ddU, ddR, X, Y, dX, dY, X1, X2, Y1, Y2
	INTEGER					:: W_A_E_T, W_A_X_T, W_A_Y_T, W_A_Z_T, W_A_R_T

	W_A_E_T = 1
	W_A_X_T = 1
	W_A_Y_T = 1
	W_A_Z_T = 1
	W_A_R_T = 1

	ddT  = (X_CLICK(2,1) - X_CLICK(1,1))/2.0D0
	ddE  = (X_CLICK(2,2) - X_CLICK(1,2))/2.0D0
	ddXYZ= (X_CLICK(2,3) - X_CLICK(1,3))/2.0D0
	ddR  = (X_CLICK(2,6) - X_CLICK(1,6))/2.0D0
	ddU  = (X_CLICK(2,7) - X_CLICK(1,7))/2.0D0
	ddN  = (X_CLICK(2,12) - X_CLICK(1,12))/2.0D0
	dEL  = (X_CLICK(2,13) - X_CLICK(1,13))/2.0D0

	IF (W_A_E_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <E> vs T
		! Energy vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 2 ! energy
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddE
		dX = ddT
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						A_E_T(i,j) = A_E_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF

	IF (W_A_X_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <X> vs T
		! X-pos vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 3 ! x-pos
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddXYZ
		dX = ddT
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						A_X_T(i,j) = A_X_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF

	IF (W_A_Y_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <Y> vs T
		! Y-pos vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 4 ! y-pos
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddXYZ
		dX = ddT
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						A_Y_T(i,j) = A_Y_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF

	IF (W_A_Z_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <Z> vs T
		! z-pos vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 5 ! z-pos
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddXYZ
		dX = ddT
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						A_Z_T(i,j) = A_Z_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF

	IF (W_A_R_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		!  y  vs x
		! <R> vs T
		! R-pos vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 6 ! R-pos
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dY = ddR
		dX = ddT
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						A_R_T(i,j) = A_R_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF

END SUBROUTINE fill_all_lism_phase_spaces

!#######################################
!#######################################
!#######################################

!! Fill all phase spaces which will
!! be written to files for data
!! analysis

SUBROUTINE fill_lism_phase_spaces( ck, MC )

	USE lism_click
	USE lism, ONLY : N_Part

	IMPLICIT NONE

	!! Input
	INTEGER					:: ck		! click
	INTEGER					:: MC		! current MC particle

	!! Internal
	INTEGER					:: i, j, iX, iY, iz
	REAL(KIND=8)		:: ddT, ddE, ddN, dEL, ddXYZ, ddR, ddU, X, Y, dX, dY, X1, X2, Y1, Y2
	REAL(KIND=8)		:: ddSX, ddSR, Z
  INTEGER         :: W_C_E_T, W_C_X_T, W_C_Y_T, W_C_Z_T, W_C_R_T, W_SC_R_T
	INTEGER					:: W_C_X_Y, W_C_X_Z, W_C_Y_Z, W_SC_X_Y, W_SC_X_Z, W_SC_Y_Z
	INTEGER					:: W_A_E_T, AVE_ENG

	AVE_ENG  = 1

	W_A_E_T  = 1
  W_C_E_T  = 1
  W_C_X_T  = 1
  W_C_Y_T  = 1
  W_C_Z_T  = 1
  W_C_R_T  = 1
  W_C_X_Y  = 1
  W_C_X_Z  = 1
  W_C_Y_Z  = 1
  W_SC_R_T = 1
  W_SC_X_Y = 1
  W_SC_X_Z = 1
  W_SC_Y_Z = 1

  ddT  = (X_CLICK(2,1) - X_CLICK(1,1))/2.0D0
  ddE  = (X_CLICK(2,2) - X_CLICK(1,2))/2.0D0
  ddXYZ= (X_CLICK(2,3) - X_CLICK(1,3))/2.0D0
  ddR  = (X_CLICK(2,6) - X_CLICK(1,6))/2.0D0
  ddU  = (X_CLICK(2,7) - X_CLICK(1,7))/2.0D0
  ddN  = (X_CLICK(2,12) - X_CLICK(1,12))/2.0D0
  dEL  = (X_CLICK(2,13) - X_CLICK(1,13))/2.0D0
	ddSX = (X_CLICK(2,14) - X_CLICK(1,14))/2.0D0
	ddSR = (X_CLICK(2,17) - X_CLICK(1,17))/2.0D0

	IF ( AVE_ENG .EQ. 1 ) THEN
		X = CLICK_1(MC,3)
		Y = CLICK_1(MC,4)		
		Z = CLICK_1(MC,5)
		ix = 0
		iy = 0
		iz = 0
		DO i=1,NUM_Hist
			X1 = X_CLICK(i,3) - ddXYZ
			X2 = X_CLICK(i,3) + ddXYZ
			IF ((X1 .LE. X) .AND. (X .LE. X2)) THEN
				ix = i
			END IF
		END DO		
		DO i=1,NUM_Hist
			X1 = X_CLICK(i,4) - ddXYZ
			X2 = X_CLICK(i,4) + ddXYZ
			IF ((X1 .LE. Y) .AND. (Y .LE. X2)) THEN
				iy = i
			END IF
		END DO		
		DO i=1,NUM_Hist
			X1 = X_CLICK(i,5) - ddXYZ
			X2 = X_CLICK(i,5) + ddXYZ
			IF ((X1 .LE. Z) .AND. (Z .LE. X2)) THEN
				iz = i
			END IF
		END DO		
		IF ( (ix .GT. 0) .AND. (iy .GT. 0) .AND. (ABS(Z) .LE. 2.0D0*ddXYZ) ) THEN
			A_E_X_Y(ix,iy) = A_E_X_Y(ix,iy) + CLICK_1(MC,2)
			A_C_X_Y(ix,iy) = A_C_X_Y(ix,iy) + 1.0D0	
		END IF
		IF ( (ix .GT. 0) .AND. (iz .GT. 0) .AND. (ABS(Y) .LE. 2.0D0*ddXYZ) ) THEN
			A_E_X_Z(ix,iz) = A_E_X_Z(ix,iz) + CLICK_1(MC,2)
			A_C_X_Z(ix,iz) = A_C_X_Z(ix,iz) + 1.0D0	
		END IF
		IF ( (iy .GT. 0) .AND. (iz .GT. 0) .AND. (ABS(X) .LE. 2.0D0*ddXYZ) ) THEN
			A_E_Y_Z(iy,iz) = A_E_Y_Z(iy,iz) + CLICK_1(MC,2)
			A_C_Y_Z(iy,iz) = A_C_Y_Z(iy,iz) + 1.0D0	
		END IF
	END IF

	IF ( W_A_E_T .EQ. 1 ) THEN
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !  y  vs x
    ! <E> vs T
    ! Energy vs Time Prob Den
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iY = 2 ! energy
    iX = 1 ! time
    X  = CLICK_1(MC,iX)
    Y  = CLICK_1(MC,iY)
    dY = ddE
    dX = ddT
    DO i=1,Num_Hist
      X1 = X_CLICK(i,iX)-dX
      X2 = X_CLICK(i,iX)+dX
      DO j=1,Num_Hist
        Y1 = X_CLICK(j,iY)-dY
        Y2 = X_CLICK(j,iY)+dY
        IF ((X.GE.X1).AND.(X.LE.X2)) THEN
          IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
            A_E_T(i,j) = A_E_T(i,j) + 1.0D0/REAL(N_Part*Max_Click)
          END IF ! Y
        END IF ! X 
      END DO ! j 
    END DO ! i
	END IF

  IF ( W_SC_R_T .EQ. 1) THEN
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! r-pos vs Time Prob Den
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iY = 6 ! r-pos
    iX = 1 ! time
    X  = CLICK_1(MC,iX)
    Y  = CLICK_1(MC,iY)
    dY = ddSR
    dX = ddT
    DO i=1,Num_Hist
      X1 = X_CLICK(i,iX)-dX
      X2 = X_CLICK(i,iX)+dX
      DO j=1,Num_Hist
        Y1 = X_CLICK(j,17)-dY
        Y2 = X_CLICK(j,17)+dY
        IF ((X.GE.X1).AND.(X.LE.X2)) THEN
          IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
            Cs_R_T(ck,i,j) = Cs_R_T(ck,i,j) + 1.0D0/REAL(N_Part)
          END IF ! Y
        END IF ! X 
      END DO ! j 
    END DO ! i
  END IF ! test data write on

	IF ( W_SC_X_Y .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! X vs Y Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iZ = 5 ! z-pos
		iY = 4 ! y-pos
		iX = 3 ! x-pos
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		Z  = CLICK_1(MC,iZ)
		dX = ddSX
		dY = ddSX
		IF ( ABS(Z) .LE. 2.0D0*dX ) THEN
			DO i=1,Num_Hist
				X1 = X_CLICK(i,14)-dX
				X2 = X_CLICK(i,14)+dX
				DO j=1,Num_Hist
					Y1 = X_CLICK(j,15)-dY
					Y2 = X_CLICK(j,15)+dY
					IF ((X.GE.X1).AND.(X.LE.X2)) THEN
						IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
							Cs_X_Y(ck,i,j) = Cs_X_Y(ck,i,j) + 1.0D0/REAL(N_Part)
						END IF ! Y
					END IF ! X 
				END DO ! j 
			END DO ! i
		END IF ! within z = 0 layer
	END IF ! test data write on

	IF ( W_SC_X_Z .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! X vs Z Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 5 ! z-pos
		iX = 3 ! x-pos
		iZ = 4 ! y-pos
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		Z  = CLICK_1(MC,iZ)
		dX = ddSX
		dY = ddSX
		IF ( ABS(Z) .LE. 2.0D0*dX ) THEN
			DO i=1,Num_Hist
				X1 = X_CLICK(i,14)-dX
				X2 = X_CLICK(i,14)+dX
				DO j=1,Num_Hist
					Y1 = X_CLICK(j,16)-dY
					Y2 = X_CLICK(j,16)+dY
					IF ((X.GE.X1).AND.(X.LE.X2)) THEN
						IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
							Cs_X_Z(ck,i,j) = Cs_X_Z(ck,i,j) + 1.0D0/REAL(N_Part)
						END IF ! Y
					END IF ! X 
				END DO ! j 
			END DO ! i
		END IF ! within z = 0 layer
	END IF ! test data write on

	IF ( W_SC_Y_Z .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! Y vs Z Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 5 ! z-pos
		iX = 4 ! y-pos
		iZ = 3 ! x-pos
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		Z  = CLICK_1(MC,iZ)
		dX = ddSX
		dY = ddSX
		IF ( ABS(Z) .LE. 2.0D0*dX ) THEN
			DO i=1,Num_Hist
				X1 = X_CLICK(i,15)-dX
				X2 = X_CLICK(i,15)+dX
				DO j=1,Num_Hist
					Y1 = X_CLICK(j,16)-dY
					Y2 = X_CLICK(j,16)+dY
					IF ((X.GE.X1).AND.(X.LE.X2)) THEN
						IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
							Cs_Y_Z(ck,i,j) = Cs_Y_Z(ck,i,j) + 1.0D0/REAL(N_Part)
						END IF ! Y
					END IF ! X 
				END DO ! j 
			END DO ! i
		END IF ! within z = 0 layer
	END IF ! test data write on

	IF ( W_C_X_Y .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! X vs Y Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 4 ! y-pos
		iX = 3 ! x-pos
		iZ = 5 ! z-pos
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		Z  = CLICK_1(MC,iZ)
		dX = ddXYZ
		dY = ddXYZ
		IF ( ABS(Z) .LE. 2.0D0*dX ) THEN
			DO i=1,Num_Hist
				X1 = X_CLICK(i,iX)-dX
				X2 = X_CLICK(i,iX)+dX
				DO j=1,Num_Hist
					Y1 = X_CLICK(j,iY)-dY
					Y2 = X_CLICK(j,iY)+dY
					IF ((X.GE.X1).AND.(X.LE.X2)) THEN
						IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
							C_X_Y(ck,i,j) = C_X_Y(ck,i,j) + 1.0D0/REAL(N_Part)
						END IF ! Y
					END IF ! X 
				END DO ! j 
			END DO ! i
		END IF ! within z = 0 layer
	END IF ! test data write on

	IF ( W_C_X_Z .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! X vs Z Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 5 ! z-pos
		iX = 3 ! x-pos
		iZ = 4 ! y-pos
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		Z  = CLICK_1(MC,iZ)
		dX = ddXYZ
		dY = ddXYZ
		IF ( ABS(Z) .LE. 2.0D0*dX ) THEN
			DO i=1,Num_Hist
				X1 = X_CLICK(i,iX)-dX
				X2 = X_CLICK(i,iX)+dX
				DO j=1,Num_Hist
					Y1 = X_CLICK(j,iY)-dY
					Y2 = X_CLICK(j,iY)+dY
					IF ((X.GE.X1).AND.(X.LE.X2)) THEN
						IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
							C_X_Z(ck,i,j) = C_X_Z(ck,i,j) + 1.0D0/REAL(N_Part)
						END IF ! Y
					END IF ! X 
				END DO ! j 
			END DO ! i
		END IF ! within z = 0 layer
	END IF ! test data write on

	IF ( W_C_Y_Z .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! Y vs Z Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 5 ! z-pos
		iX = 4 ! y-pos
		iZ = 3 ! x-pos
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		Z  = CLICK_1(MC,iZ)
		dX = ddXYZ
		dY = ddXYZ
		IF ( ABS(Z) .LE. 2.0D0*dX ) THEN
			DO i=1,Num_Hist
				X1 = X_CLICK(i,iX)-dX
				X2 = X_CLICK(i,iX)+dX
				DO j=1,Num_Hist
					Y1 = X_CLICK(j,iY)-dY
					Y2 = X_CLICK(j,iY)+dY
					IF ((X.GE.X1).AND.(X.LE.X2)) THEN
						IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
							C_Y_Z(ck,i,j) = C_Y_Z(ck,i,j) + 1.0D0/REAL(N_Part)
						END IF ! Y
					END IF ! X 
				END DO ! j 
			END DO ! i
		END IF ! within z = 0 layer
	END IF ! test data write on

	IF ( W_C_E_T .EQ. 1) THEN
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		! Energy vs Time Prob Den
		!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		iY = 2 ! energy
		iX = 1 ! time
		X  = CLICK_1(MC,iX)
		Y  = CLICK_1(MC,iY)
		dX = ddT
		dY = ddE
		DO i=1,Num_Hist
			X1 = X_CLICK(i,iX)-dX
			X2 = X_CLICK(i,iX)+dX
			DO j=1,Num_Hist
				Y1 = X_CLICK(j,iY)-dY
				Y2 = X_CLICK(j,iY)+dY
				IF ((X.GE.X1).AND.(X.LE.X2)) THEN
					IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
						C_E_T(ck,i,j) = C_E_T(ck,i,j) + 1.0D0/REAL(N_Part)
					END IF ! Y
				END IF ! X 
			END DO ! j 
		END DO ! i
	END IF ! test data write on

  IF ( W_C_X_T .EQ. 1) THEN
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! X-pos vs Time Prob Den
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iY = 3 ! x-pos
    iX = 1 ! time
    X  = CLICK_1(MC,iX)
    Y  = CLICK_1(MC,iY)
    dY = ddXYZ
    dX = ddT
    DO i=1,Num_Hist
      X1 = X_CLICK(i,iX)-dX
      X2 = X_CLICK(i,iX)+dX
      DO j=1,Num_Hist
        Y1 = X_CLICK(j,iY)-dY
        Y2 = X_CLICK(j,iY)+dY
        IF ((X.GE.X1).AND.(X.LE.X2)) THEN
          IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
            C_X_T(ck,i,j) = C_X_T(ck,i,j) + 1.0D0/REAL(N_Part)
          END IF ! Y
        END IF ! X 
      END DO ! j 
    END DO ! i
  END IF ! test data write on

  IF ( W_C_Y_T .EQ. 1) THEN
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! Y-pos vs Time Prob Den
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iY = 4 ! y-pos
    iX = 1 ! time
    X  = CLICK_1(MC,iX)
    Y  = CLICK_1(MC,iY)
    dY = ddXYZ
    dX = ddT
    DO i=1,Num_Hist
      X1 = X_CLICK(i,iX)-dX
      X2 = X_CLICK(i,iX)+dX
      DO j=1,Num_Hist
        Y1 = X_CLICK(j,iY)-dY
        Y2 = X_CLICK(j,iY)+dY
        IF ((X.GE.X1).AND.(X.LE.X2)) THEN
          IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
            C_Y_T(ck,i,j) = C_Y_T(ck,i,j) + 1.0D0/REAL(N_Part)
          END IF ! Y
        END IF ! X 
      END DO ! j 
    END DO ! i
  END IF ! test data write on

  IF ( W_C_Z_T .EQ. 1) THEN
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! z-pos vs Time Prob Den
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iY = 5 ! z-pos
    iX = 1 ! time
    X  = CLICK_1(MC,iX)
    Y  = CLICK_1(MC,iY)
    dY = ddXYZ
    dX = ddT
    DO i=1,Num_Hist
      X1 = X_CLICK(i,iX)-dX
      X2 = X_CLICK(i,iX)+dX
      DO j=1,Num_Hist
        Y1 = X_CLICK(j,iY)-dY
        Y2 = X_CLICK(j,iY)+dY
        IF ((X.GE.X1).AND.(X.LE.X2)) THEN
          IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
            C_Z_T(ck,i,j) = C_Z_T(ck,i,j) + 1.0D0/REAL(N_Part)
          END IF ! Y
        END IF ! X 
      END DO ! j 
    END DO ! i
  END IF ! test data write on

  IF ( W_C_R_T .EQ. 1) THEN
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    ! r-pos vs Time Prob Den
    !^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    iY = 6 ! r-pos
    iX = 1 ! time
    X  = CLICK_1(MC,iX)
    Y  = CLICK_1(MC,iY)
    dY = ddR
    dX = ddT
    DO i=1,Num_Hist
      X1 = X_CLICK(i,iX)-dX
      X2 = X_CLICK(i,iX)+dX
      DO j=1,Num_Hist
        Y1 = X_CLICK(j,iY)-dY
        Y2 = X_CLICK(j,iY)+dY
        IF ((X.GE.X1).AND.(X.LE.X2)) THEN
          IF ((Y.GE.Y1).AND.(Y.LE.Y2)) THEN
            C_R_T(ck,i,j) = C_R_T(ck,i,j) + 1.0D0/REAL(N_Part)
          END IF ! Y
        END IF ! X 
      END DO ! j 
    END DO ! i
  END IF ! test data write on


END SUBROUTINE fill_lism_phase_spaces

!#######################################
!#######################################

SUBROUTINE write_lism_phase_spaces
  USE lism_click
  USE lism, ONLY : N_Part, Tot_Click, Click_Skip
	USE mpi_info

	IMPLICIT NONE

	INCLUDE 'mpif.h'

  INTEGER         :: W_C_X_Y, W_C_X_Z, W_C_Y_Z
  INTEGER         :: W_C_E_T, W_C_X_T, W_C_Y_T, W_C_Z_T, W_C_R_T
  INTEGER         :: W_A_E_T, W_A_X_T, W_A_Y_T, W_A_Z_T, W_A_R_T
	INTEGER					:: W_SC_X_Y, W_SC_X_Z, W_SC_Y_Z, W_SC_R_T
	INTEGER					:: AVE_ENG
	INTEGER					:: ck, i, j

	AVE_ENG = 1

  W_A_E_T = 1
  W_A_X_T = 1
  W_A_Y_T = 1
  W_A_Z_T = 1
  W_A_R_T = 1

  W_C_X_Y = 1
  W_C_X_Z = 1
  W_C_Y_Z = 1
  W_C_E_T = 1
  W_C_X_T = 1
  W_C_Y_T = 1
  W_C_Z_T = 1
  W_C_R_T = 1
  W_SC_X_Y = 1
  W_SC_X_Z = 1
  W_SC_Y_Z = 1
  W_SC_R_T = 1

	IF ( AVE_ENG .EQ. 1 ) THEN
		IF ( myid .EQ. 0) THEN
			ALLOCATE( ROOT_ALL( Num_Hist, Num_Hist) )
			ALLOCATE( ROOT_ALL2( Num_Hist, Num_Hist) )
		END IF
    CALL MPI_REDUCE( A_E_X_Y(:,:), ROOT_ALL(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_REDUCE( A_C_X_Y(:,:), ROOT_ALL2(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_all_avg_energy_XY.dat",ACCESS="APPEND")
      OPEN(UNIT=42,FILE="../Data/lism_all_count_XY.dat",ACCESS="APPEND")
      DO i=1,Num_Hist
        DO j=1,Num_Hist
        	IF (ROOT_ALL2(i,j) .GT. 0.0D0) THEN 
						WRITE(41,*) i, j, ROOT_ALL(i,j)/ROOT_ALL2(i,j)
					ELSE
						WRITE(41,*) i, j, 0.0D0
					END IF
					WRITE(42,*) i, j, ROOT_ALL2(i,j)
        END DO ! j
      END DO ! i      
			CLOSE(41)
			CLOSE(42)
			ROOT_ALL 	= 0.0D0
			ROOT_ALL2 = 0.0D0
		END IF ! root
  	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_REDUCE( A_E_X_Z(:,:), ROOT_ALL(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_REDUCE( A_C_X_Z(:,:), ROOT_ALL2(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_all_avg_energy_XZ.dat",ACCESS="APPEND")
      OPEN(UNIT=42,FILE="../Data/lism_all_count_XZ.dat",ACCESS="APPEND")
      DO i=1,Num_Hist
        DO j=1,Num_Hist
        	IF (ROOT_ALL2(i,j) .GT. 0.0D0) THEN 
						WRITE(41,*) i, j, ROOT_ALL(i,j)/ROOT_ALL2(i,j)
					ELSE
						WRITE(41,*) i, j, 0.0D0
					END IF
					WRITE(42,*) i, j, ROOT_ALL2(i,j)
        END DO ! j
      END DO ! i      
			CLOSE(41)
			CLOSE(42)
			ROOT_ALL 	= 0.0D0
			ROOT_ALL2 = 0.0D0
		END IF ! root
  	CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_REDUCE( A_E_Y_Z(:,:), ROOT_ALL(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    CALL MPI_REDUCE( A_C_Y_Z(:,:), ROOT_ALL2(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_all_avg_energy_YZ.dat",ACCESS="APPEND")
      OPEN(UNIT=42,FILE="../Data/lism_all_count_YZ.dat",ACCESS="APPEND")
      DO i=1,Num_Hist
        DO j=1,Num_Hist
        	IF (ROOT_ALL2(i,j) .GT. 0.0D0) THEN 
						WRITE(41,*) i, j, ROOT_ALL(i,j)/ROOT_ALL2(i,j)
					ELSE
						WRITE(41,*) i, j, 0.0D0
					END IF
					WRITE(42,*) i, j, ROOT_ALL2(i,j)
        END DO ! j
      END DO ! i      
			CLOSE(41)
			CLOSE(42)
			DEALLOCATE(ROOT_ALL)
			DEALLOCATE(ROOT_ALL2)
		END IF ! root
		CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
	END IF

	IF ( W_A_E_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_ALL( Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( A_E_T(:,:), ROOT_ALL(:,:), Num_Hist*Num_Hist, &
    & MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_all_energy_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write E vs T to file        
			!******************************
      DO i=1,Num_Hist
        DO j=1,Num_Hist
          WRITE(41,*) i, j, ROOT_ALL(i,j)
        END DO ! j
      END DO ! i      CLOSE(41)
      DEALLOCATE( ROOT_ALL )    
		END IF ! root
	END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_SC_X_Y .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( CS_X_Y(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_small_X_vs_Y.dat",ACCESS="APPEND")
      !******************************
      ! write X vs Y to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
	IF ( W_SC_X_Z .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( CS_X_Z(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_small_X_vs_Z.dat",ACCESS="APPEND")
      !******************************
      ! write X vs Z to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF ( W_SC_Y_Z .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( CS_Y_Z(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_small_Y_vs_Z.dat",ACCESS="APPEND")
      !******************************
      ! write X vs Z to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_C_X_Y .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_X_Y(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_X_vs_Y.dat",ACCESS="APPEND")
      !******************************
      ! write X vs Y to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
  
	IF ( W_C_X_Z .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_X_Z(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_X_vs_Z.dat",ACCESS="APPEND")
      !******************************
      ! write X vs Z to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

	IF ( W_C_Y_Z .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_Y_Z(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_Y_vs_Z.dat",ACCESS="APPEND")
      !******************************
      ! write X vs Z to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_C_E_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_E_T(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_energy_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write E vs T to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_C_X_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_X_T(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_X_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write X vs T to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_C_Y_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_Y_T(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_Y_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write Y vs T to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_C_Z_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_Z_T(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_Z_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write Z vs T to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_C_R_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( C_R_T(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_R_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write R vs T to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )

  IF ( W_SC_R_T .EQ. 1) THEN
    IF ( myid .EQ. 0 ) ALLOCATE( ROOT_CLICK( Tot_Click, Num_Hist, Num_Hist) )
    CALL MPI_REDUCE( CS_R_T(:,:,:), ROOT_CLICK(:,:,:), Tot_Click*Num_Hist*Num_Hist, &
    & 								MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr )
    CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )
    IF ( myid .EQ. 0 ) THEN
      OPEN(UNIT=41,FILE="../Data/lism_small_R_vs_time.dat",ACCESS="APPEND")
      !******************************
      ! write R vs T to file  
      !******************************
      DO ck=1,Tot_Click
        DO i=1,Num_Hist
          DO j=1,Num_Hist
            WRITE(41,*) ck*Click_Skip, i, j, ROOT_CLICK(ck,i,j)
          END DO ! j
        END DO ! i
      END DO ! ck
      CLOSE(41)
      DEALLOCATE( ROOT_CLICK )
    END IF ! root
  END IF
  CALL MPI_BARRIER( MPI_COMM_WORLD, ierr )


END SUBROUTINE write_lism_phase_spaces


