program PlacaPlana

    	use rectangular_wing_module
	use uvlm
	use solvers
	use tecplot_wrappers

      	implicit none

      	type( RectWingType ) :: wing
	real(8), dimension(:,:), allocatable :: A_mat, A_copy, RHS
	real(8), dimension(:,:), allocatable :: cp_mat
	real(8), dimension(3,3) :: w_c_sys
      	real(8), dimension(3) :: V_B_dim, V_inf_dim, V_B, V_inf
      	real(8), dimension(3) :: w_origin, w_wvel, euler_ang
	real(8) :: V_char, L_char, delta_t
      	real(8) :: span, chord
      	real(8) :: span_dim, chord_dim
      	real(8) :: AR, alfa, core
      	integer   :: disc, ns, nc, globnp
	integer   :: t_steps, wake_size
	integer   :: t, i, an
	character(len=100) :: numchar

	!---------
	! Normalizing velocities

      	V_inf_dim  = 0.0d+0
      	V_B_dim    = [-100.0d+0, 0.d+0, 0.d+0]
	V_char = norm2( V_inf_dim - V_B_dim )
	
	!V_inf = V_inf_dim / V_char
	!V_B = V_B_dim / V_char
	V_inf = 0.d+0 !V_inf_dim / V_char
	V_B = 0.d+0 !V_B_dim / V_char

	!----------
	! Wing geometry params

	disc = 4
	AR = 20.d+0
      	chord_dim = 1.d+0
      	span_dim  = chord_dim * AR

      	ns = 68 !disc * int(AR)
      	nc = 10 !disc * 2

	allocate( cp_mat( ns*nc, 5 ))

	L_char = chord_dim / nc

	chord  = chord_dim / L_char
	span   = span_dim  / L_char

	core = 0.1d-1 * max( span/ns, chord/nc )

	!---------
	! Time stepping params

	delta_t   = 1.d+0
	t_steps   = 25 * nc
	wake_size = t_steps

	!---------
	! Angle of attack

	do an = 0,4

		alfa = an * 2.5d+0
		euler_ang = [0.d+0, 0.d+0, 0.d+0]
		V_inf = [cos(deg2rad(alfa)), 0.d+0, sin(deg2rad(alfa))]
		!------------
		! Wing position params

		w_c_sys  = reshape([1.d+0, 0.d+0, 0.d+0,&
				    0.d+0, 1.d+0, 0.d+0,&
				    0.d+0, 0.d+0, 1.d+0],[3,3])

 		!call euler_rot(2, alfa, w_c_sys)
		w_origin = [-2.0d+0, 0.d+0, 0.d+0]
		w_wvel   = 0.d+0

      		wing = initRectWing( span, chord, ns, nc, wake_size, w_c_sys, w_origin, V_B, w_wvel )

		!call printGeometry

		globnp = wing%NP

		!---------------
		! Setting up linear system por T = 0

		if( an == 0 ) then
			call allocateLinSystem
		endif

		A_mat = 0.0d+0
		RHS   = 0.0d+0

		call infMat( wing%NODESXYZ, wing%LOCMAT, wing%NBMAT, wing%CPXYZ, wing%NORMAL, A_mat )

		!--------------
		! Starting time advance

		write(numchar,*) int(alfa)
		call iniTec360File( 'placaplana'//trim(adjustl(numchar)), 'X Y Z Cp', 1 )

		time: do t = 0,t_steps
			
			!-------------------
			! convecting Wake

			if ( t > 0 ) then

				wing%wake%NODESUVW = 0.d+0
				
				call inducedVelocity( wing%wake%NODESXYZ, wing%wake%NODESUVW(:,1:wing%NEN*t), wing, t, core, V_inf)!, debug=.true. )
			
				call updateCoordinates( wing, euler_ang, delta_t ) 

				call convect( wing, delta_t, t )

				!call printWakePositions
				!call printWakeG

			endif
		
			!-----------
			! Setting up linear system
			
			A_copy = A_mat
			RHS = 0.d+0

			call RHSCALC( wing, wing, RHS, t, core, V_inf )
			!call printRHS

			!--------------
			! Solving Linear System

			call gesv( A_copy, RHS )
			!call printResult


			!--------------
			! Assingning Loop circulations on surface

			call updateG( wing, RHS )

			!-----------
			! Calculating induced velocities

			wing%MVEL = 0.d+0

			print*, 'calculating induced velocities on the wing'
			call inducedVelocity( wing%CPXYZ, wing%MVEL, wing, t, core, V_inf )

			call LOADCALC( wing, delta_t )

			if( t == t_steps ) then
				cp_mat(:,an+1) = wing%LOAD(:,1)
			endif

			SolTime = t * delta_t
			write(numchar,*) t

			StrandID = 1
			call newTec360Zone( 'placa'//trim(adjustl(numchar)), wing%NODESXYZ, wing%LOCMAT, wing%LOAD, [1, 1, 1, 0], 1 )

			StrandID = 2
			if( t > 0 ) then
				call newTec360Zone( 'estela'//trim(adjustl(numchar)), wing%wake%NODESXYZ(:,1:wing%NEN*(t+1)),&
							wing%wake%LOCMAT(:,1:wing%NEP*t), wing%wake%LOOPG(1:wing%NEP*t,:), [1, 1, 1, 0], 1 )
			endif

		end do time

		i = TecEnd112()

	        call wing%DESTROY()

	enddo

	call writeLoads

	contains

	subroutine writeLoads

		open(unit=11, file='placaplana.dat')
		write(11,*) '########################'
		write(11,*) '## Resultados CP      ##'
		write(11,*) '########################'
		write(11,'(5f18.12)') (cp_mat( i,: ), i=1,ns*nc)
		write(11,*) 
		close(11)

	end subroutine

	subroutine printLoads
		
		print*
		print*, '---------------'
		print*, ' Body Cp', t
		print'(1f18.12)', wing%LOAD
		print*

	end subroutine 

	subroutine printGammas

		print*
		print*, '---------------'
		print*, ' Body Gammas', t
		print'(4f18.12)', wing%GAMMAS
		print*
		print*, '---------------'
		print*, ' Wake Gammas', t
		print'(4f18.12)', wing%wake%GAMMAS
		print*
		
	end subroutine printGammas 

	subroutine printWakePositions

		print*
		print*, '---------------'
		print*, ' Wake positions for timestep', t
		print'(3f18.12)', wing%wake%NODESXYZ(:,1:wing%NEN * (t+1))
		print*

	end subroutine printWakePositions

	subroutine printRHS

		print*
		print*, '-----------------------------'
		print*, 'RHS for timestep', t
		print'(1f18.12)', RHS
		print*

	end subroutine printRHS

	subroutine printWakeG

		print*
		print*, '-----------------------------'
		print*, 'Wake G', t
		print'(1f18.12)', wing%wake%LOOPG(1:wing%NEP*t,1)
		print*

	end subroutine printWakeG 

	subroutine printResult

		print*
		print*, '-----------------------------'
		print*, 'Result for timestep', t
		print'(1f18.12)', RHS
		print*

	end subroutine printResult

	subroutine printGeometry

		print*
		print*, 'local nodes'
		print'(3f12.6)', wing%NODESLOC
		print*
		print*, 'global nodes'
		print'(3f12.6)', wing%NODESXYZ
		print*
		print*, 'global CPNODES'
		print'(3f12.6)', wing%CPXYZ
		print*
		print*, 'global cpvels'
		print'(3f12.6)', wing%CPUVW
		print*
		print*, 'Location Matrix'
		print'(4i8)', wing%LOCMAT
		print*
		print*, 'Neibour Matrix'
		print'(4i8)', wing%NBMAT
		print*
		print*, 'Wake Location Matrix'
		print'(4i8)', wing%wake%LOCMAT
		print*
		print*, 'Wake Neibour Matrix'
		print'(4i8)', wing%wake%NBMAT
		print*
		print*, 'Wake TELSNI'
		print'(1i8)', wing%TELSNI
		print*
		print*, 'Wake TELSSI'
		print'(1i8)', wing%TELSSI
		print*
		print*, 'Wake TEWAKEI'
		print'(1i8)', wing%wake%TEWAKEI
		print*

	end subroutine printGeometry

    	subroutine allocateLinSystem

		allocate( A_mat  ( globnp,globnp ) )
		allocate( A_copy ( globnp,globnp ) )
		allocate( RHS    ( globnp, 1     ) )

	end subroutine allocateLinSystem

	subroutine deallocateLinSystem

		deallocate( A_mat  )
		deallocate( A_copy )
		deallocate( RHS    )

	end subroutine deallocateLinSystem 

end program PlacaPlana
