
module elastic_atom_ion

	real(kind=8),allocatable,dimension(:)	:: He_Hp_E
	real(kind=8),allocatable,dimension(:)	:: He_Hp_T
	real(kind=8),allocatable,dimension(:)	:: He_Hep_E
	real(kind=8),allocatable,dimension(:)	:: He_Hep_T

	integer																:: N_He_Hp
	integer																:: N_He_Hep

	real(kind=8),parameter								:: lam_He_Hp = 1.5D0
	real(kind=8),parameter								:: lam_He_Hep = 1.5D0
	real(kind=8),parameter								:: C0_S = 4.2D0
	real(kind=8),parameter								:: C1_S = -0.059D0
	real(kind=8),parameter								:: C0_B = 5.9D0
	real(kind=8),parameter								:: C1_B = -0.76D0
	real(kind=8),parameter								:: CX0  = 2.4251D0

contains

!!!!!!!!!!!!!

subroutine test_elastic_atom_ion

	implicit none

	real(kind=8)	:: Ein, Efn, E, dE, tcs, ang, t_ang
	integer				:: i, j, N, N_MC

	N    = 100
	N_MC = 10000	

	Ein  = 1.0D0
	Efn  = 10000.0D0
	dE   = (Efn-Ein)/REAL(N-1)
	
	call read_elastic_atom_ion_tables

	open(unit=10,file='../Data/He_Hp_test.dat',access='append')
	open(unit=20,file='../Data/He_Hep_test.dat',access='append')

	WRITE(*,*) 'He+Hp TEST'
	WRITE(*,*) '----------'
	do i=1,N
		E = Ein + (i-1)*dE
		call elastic_He_Hp_tcs( E, tcs )
		write(*,*) 'E: tcs: ', E, tcs
		t_ang= 0.0D0
		do j=1,N_MC
			call elastic_He_Hp_rand_angle( E, ang )				
!			write(*,*) 'ang: ', ang
			t_ang = t_ang + ang
		end do
		WRITE(10,*) E, tcs, t_ang/REAL(N_MC)
		WRITE(*,*) E, tcs, t_ang/REAL(N_MC)
	end do	

	WRITE(*,*) 'He+Hep TEST'
	WRITE(*,*) '----------'
	do i=1,N
		E = Ein + (i-1)*dE
		call elastic_He_Hep_tcs( E, tcs )
		t_ang = 0.0D0
		do j=1,N_MC
			call elastic_He_Hep_rand_angle( E, ang )				
			t_ang = t_ang + ang
		end do
		WRITE(20,*) E, tcs, t_ang/REAL(N_MC)
		WRITE(*,*) E, tcs, t_ang/REAL(N_MC)
	end do	

	close(10)
	close(20)
	call clean_elastic_atom_ion_tables

end subroutine

!!!!!!!!!!!!!

subroutine read_elastic_atom_ion_tables

	use mpi_info

	implicit none

	include 'mpif.h'

	integer :: i

	if (myid .eq. 0) then
    OPEN(UNIT=10, FILE="../Tables/Elastic_Hp_He_TCS.dat", STATUS="old",ACTION="read")
    OPEN(UNIT=20, FILE="../Tables/Elastic_Hep_He_TCS.dat", STATUS="old",ACTION="read")
		read(10,*) N_He_Hp
		read(20,*) N_He_Hep
		write(*,*) N_He_Hp, N_He_Hep
		allocate( He_Hp_E(N_He_Hp), He_Hp_T(N_He_Hp) )
		allocate( He_Hep_E(N_He_Hp), He_Hep_T(N_He_Hp) )
		do i=1,N_He_Hp
			read(10,*) He_Hp_E(i), He_Hp_T(i)
		end do
		do i=1,N_He_Hep
			read(20,*) He_Hep_E(i), He_Hep_T(i)
		end do
		close(10)
		close(20)
	end if	

	call MPI_BARRIER( MPI_COMM_WORLD, ierr )
	call MPI_BCAST( N_He_Hp,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	call MPI_BCAST( N_He_Hep,  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr )
	if (myid .ne. 0) then
		allocate( He_Hp_E(N_He_Hp), He_Hp_T(N_He_Hp) )
		allocate( He_Hep_E(N_He_Hp), He_Hep_T(N_He_Hp) )
	end if	
	call MPI_BARRIER( MPI_COMM_WORLD, ierr )
	call MPI_BCAST( He_Hp_E, N_He_Hp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	call MPI_BCAST( He_Hp_T, N_He_Hp, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	call MPI_BCAST( He_Hep_E, N_He_Hep, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	call MPI_BCAST( He_Hep_T, N_He_Hep, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr )
	call MPI_BARRIER( MPI_COMM_WORLD, ierr )

end subroutine read_elastic_atom_ion_tables

!!!!!!!!!!!!!

subroutine clean_elastic_atom_ion_tables

	implicit none

	deallocate( He_Hp_E, He_Hp_T )
	deallocate( He_Hep_E, He_Hep_T )

end subroutine clean_elastic_atom_ion_tables

!!!!!!!!!!!!!

subroutine elastic_He_Hp_tcs( E, tcs )

	implicit none

	!! inputs
	real(kind=8)	:: E 		! center of mass [eV]

	!! outputs
	real(kind=8)	:: tcs 	! [m^2]

	!! internal
	integer				:: i, i0, i1
	real(kind=8)	:: E_now, x0, x1, y0, y1, m, b

	i = 1
	E_now = He_Hp_E(i)

	do while ( (i .le. N_He_Hp) .and. (E_now .le. E) )
		i = i + 1
		E_now = He_Hp_E(i)
	end do	

	if (i .eq. 1) then
		i1 = 2
		i0 = 1
	else
		i1 = i
		i0 = i-1
	end if

	x0 = He_Hp_E(i0)	
	x1 = He_Hp_E(i1)	
	y0 = He_Hp_T(i0)	
	y1 = He_Hp_T(i1)	

	tcs = y0 + (y1-y0)*((E-x0)/(x1-x0))

end subroutine elastic_He_Hp_tcs

!!!!!!!!!!!!!

subroutine elastic_He_Hep_tcs( E, tcs )

	implicit none

	!! inputs
	real(kind=8)	:: E 		! center of mass [eV]

	!! outputs
	real(kind=8)	:: tcs 	! [m^2]

	!! internal
	integer				:: i, i0, i1
	real(kind=8)	:: E_now, x0, x1, y0, y1, m, b

	i = 1
	E_now = He_Hep_E(i)

	do while ( (i .le. N_He_Hep) .and. (E_now .le. E) )
		i = i + 1
		E_now = He_Hep_E(i)
	end do	

  if (i .eq. 1) then
    i1 = 2
    i0 = 1
  else
    i1 = i
    i0 = i-1
  end if

	x0 = He_Hep_E(i0)	
	x1 = He_Hep_E(i1)	
	y0 = He_Hep_T(i0)	
	y1 = He_Hep_T(i1)	

	tcs = y0 + (y1-y0)*((E-x0)/(x1-x0))

end subroutine elastic_He_Hep_tcs

!!!!!!!!!!!!!

subroutine elastic_He_Hp_PD( E, theta, tcs, PD )

	use physics_constants, only : pi

	implicit none

	real(kind=8)	:: E, theta, tcs, PD
	real(kind=8)	:: Ti, Tf, dT, C, p0, T, mu, Tau, CC, amp
	integer				:: i, N

	N  = 100000
	Ti = 0.01D0
	Tf = 170.0D0
	dT = (Tf-Ti)/REAL(N-1)
	C  = pi/180.0D0
	p0 = 0.0D0
	i  = 1
	T  = Ti + (i-1)*dT
	mu = 4.0D0/5.0D0

	do while ( T .le. theta )
		Tau = E*T/mu
    CC  = lam_He_Hp/(T*SIN(T*C))
    if (LOG(Tau) .le. CX0) then
      amp = CC*EXP(C0_S + C1_S*LOG(Tau))
    else
      amp = CC*EXP(C0_B + C1_B*LOG(Tau))
    end if
		p0 = p0 + 2.0D0*pi*SIN(T*C)*amp*dT*C/tcs
		i  = i+1	
		T  = Ti + (i-1)*dT
	end do

	PD = p0

end subroutine

!!!!!!!!!!!!!

subroutine elastic_He_Hp_rand_angle( E, theta )

	use physics_constants, only : pi, BOHRTOM

	implicit none

	!! input
	real(kind=8)	:: E ! center of mass energy [eV]
	
	!! output
	real(kind=8)	:: theta ! center of mass scattering angle [deg]

	!! internal
	real(kind=8)	:: lfg
	real(kind=8)	:: tcs, ti, tf, dt, t, tot1, tot0, randy, amp, mu, Tau, CC
	real(kind=8)	:: t0, t1, p0, p1, x0, x1
	integer				:: N, i, Z

	Z  	= 1
	N  	= 100000
	i  	= 1
	ti 	= 0.01D0
	tf 	= 170.0D0
	dt  = (tf-ti)/REAL(N-1)
	tot1= 0.0D0
	tot0= 0.0D0
	mu  = 4.0D0/5.0D0

	randy = lfg()

	call elastic_He_Hp_tcs( E, tcs )	
	tcs = tcs/(BOHRTOM*BOHRTOM)
	do while (tot1 .LE. randy)
		t   = ti + (i-1)*dt
		Tau = E*t/mu
    CC  = lam_He_Hep/(t*SIN(t*pi/180.0D0))
    if (LOG(Tau) .le. CX0) then
      amp = CC*EXP(C0_S + C1_S*LOG(Tau))
    else
      amp = CC*EXP(C0_B + C1_B*LOG(Tau))
    end if
    tot1 = tot1 + dt*(pi/180.0D0)*SIN(t*pi/180.0D0)*amp*2.0D0*pi/tcs
		i    = i+1
	end do

	if ( i .eq. 2) then
		p0 = 0.0D0
		p1 = tot1
		x0 = 0.0D0
		x1 = t
		theta = p0 + (p1-p0)*((randy-x0)/(x1-x0))
	else
		theta = t
	end if
	
!	write(*,*) 'rand: theta: ', randy, theta

end subroutine elastic_He_Hp_rand_angle

!!!!!!!!!!!!!

subroutine elastic_He_Hep_rand_angle( E, theta )

	use physics_constants, only : pi, BOHRTOM

	implicit none

	!! input
	real(kind=8)	:: E ! center of mass energy [eV]
	
	!! output
	real(kind=8)	:: theta ! center of mass scattering angle [deg]

	!! internal
	real(kind=8)	:: lfg
	real(kind=8)	:: tcs, ti, tf, dt, t, tot1, tot0, randy, amp, mu, Tau, CC
	real(kind=8)	:: t0, t1, p0, p1, x0, x1
	integer				:: N, i, Z

	Z  	= 1
	N  	= 100000
	i  	= 1
	ti 	= 0.01D0
	tf 	= 170.0D0
	dt  = (tf-ti)/REAL(N-1)
	tot1= 0.0D0
	tot0= 0.0D0
	mu  = 2.0D0

	randy = lfg()

	call elastic_He_Hep_tcs( E, tcs )	
	tcs = tcs/(BOHRTOM*BOHRTOM)
	do while (tot1 .LE. randy)
		t   = ti + (i-1)*dt
		Tau = E*t/mu
    CC  = lam_He_Hep/(t*SIN(t*pi/180.0D0))
    if (LOG(Tau) .le. CX0) then
      amp = CC*EXP(C0_S + C1_S*LOG(Tau))
    else
      amp = CC*EXP(C0_B + C1_B*LOG(Tau))
    end if
    tot1 = tot1 + dt*(pi/180.0D0)*SIN(t*pi/180.0D0)*amp*2.0D0*pi/tcs
		i    = i+1
	end do

	if ( i .eq. 2) then
		p0 = 0.0D0
		p1 = tot1
		x0 = 0.0D0
		x1 = t
		theta = p0 + (p1-p0)*((randy-x0)/(x1-x0))
	else
		theta = t
	end if

!	write(*,*) 'rand: theta: ', randy, theta

end subroutine elastic_He_Hep_rand_angle

!!!!!!!!!!!!!

end module



