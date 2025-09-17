program Esfera

	use geometry_uvlm_wrappers
	use geometry_module
	use math_routines
	use uvlm
	use solvers
	use tecplot_wrappers
	use ellipticFuselage

      	implicit none

	type( GeometryType ), pointer :: GEOM
	class( SurfaceType ), pointer :: bodyp

	integer, dimension(:,:), allocatable :: geom_index

      	real(8) :: major, minor, delta_t, core

	integer :: disc

	real(8), dimension(:,:), allocatable, target :: A_mat, A_copy, RHS

      	real(8), dimension(3) :: V_inf
      	real(8), dimension(3) :: origin
      	real(8), dimension(3,1) :: euler_angles

      	integer   :: globnp
	integer   :: t_steps
	integer   :: t, tecerr, ng, i

	character(len=100), dimension( 1 ) :: zone

	!----------
	! geometry params

	minor = 1.d+0
	major = minor * 2.d+0

	disc = 20

	!---------
	! Normalization

	core = 0.10d-4
	V_inf = [1.0d+0, 0.0d+0, 0.0d+0]

	origin = [0.0d+0, 0.0d+0, 0.0d+0]
	euler_angles(:,1) = [0.0d+0, 0.0d+0, 0.0d+0]

	!--------------
	! Generating geometry

	call geom_create( geom, initEllipticGeom( major, disc, minor, disc*2, 1, euler_angles(:,1), origin,&
								[0.d+0,0.d+0,0.d+0],[0.d+0,0.d+0,0.d+0] ) )

	!!---------------
	!! Setting up linear system por T = 0

	globnp = geom_total_np( geom )
	ng = geom_count( geom )
	allocate( geom_index( 2, ng ) )
	geom_index = generateGeomIndex( geom )

	print'(2i6)', geom_index
	print*, globnp
	!call printGeometry

	call allocateLinSystem

	A_mat = 0.0d+0
	RHS   = 0.0d+0

	!!--------------
	!! Starting time advance

	call iniTec360File( 'esfera', 'X Y Z Cp', 1 )
	
	zone(1) = 'esfera'

	t_steps = 1
	delta_t = 1.d+0

	time: do t = 0,t_steps

		SolTime = t * delta_t

		!-------------------
		! convecting Wake and moving geometry

   		call GEOMCONVECT( geom, euler_angles, V_inf, delta_t, t, core )

		!call printWakeG

		!!-----------
		! Setting up linear system

		call GEOMAI( geom, geom_index, A_mat, t, core, .false. ) 

		A_copy = A_mat

		RHS = 0.d+0

		call GEOMRHS( geom, geom_index, V_inf, t, core, RHS )

		!call printRHS

		!--------------
		! Solving Linear System

		call gesv( A_copy, RHS )

		!--------------
		! Assingning Loop circulations on surface

		call GEOMG( geom, geom_index, RHS )

		!call printResult

		!-----------
		! Calculating induced mean velocities

		call GEOMIVEL( geom, V_inf, t, core )

		!---------------
		! Calculating Cp con control points

		call GEOMLOADS( geom, delta_t )

		!---------------
		! Writing results to TECPLOT .plt file
		
		call GEOMTECIO( geom, zone, t )  
		
	end do time

	tecerr = TecEnd112()

	contains

	!subroutine printAI

	!	print*, 'influence matrix for timestep', t
	!	write(numchar,*) globnp
	!	print '('//trim(adjustl(numchar))//'f18.12)', transpose( A_copy )
	!	print*

	!end subroutine

	
	!subroutine writePenetrationError
	!	integer :: i, j

	!	write(11,*) 'TIMESTEP', t

	!	do i=1,ng
	!		bodyp => geom_get_data( geom, i )
	!		write(11,*)
	!		write(11,*) 'body', i
	!		write(11,*)
	!		do j = 1,bodyp%NP
	!		write(11,'(1f18.12)') dot_product( bodyp%NORMAL(:,j), bodyp%MVEL(:,j) - bodyp%CPUVW(:,j) )
	!		enddo
	!		write(11,*)
	!		write(11,*) 'MEAN VEL'
	!		do j = 1,bodyp%NP
	!		write(11,'(3f18.12)') bodyp%MVEL(:,j)
	!		enddo
	!		write(11,*)
	!		write(11,*) 'C_p'
	!		do j = 1,bodyp%NP
	!		write(11,'(1f18.12)') bodyp%LOAD(j,1)
	!		enddo
	!		write(11,*)
	!	enddo
	!		
	!end subroutine

	!subroutine printGammas
	!	integer :: i

	!	do i=1,ng
	!	bodyp => geom_get_data( geom, i )
	!	print*
	!	print'(4f12.6)', bodyp%GAMMAS
	!	print*
	!	enddo

	!end subroutine

	subroutine printGeometry
		integer :: i

		print*, 'total bodies', ng

		do i=1,ng
		
			print*, 'body', i
			bodyp => geom_get_data( geom, i )
			print*, 'nodes'
			print'(3f18.12)', bodyp%NODESXYZ
			print*
			print*, 'panels'
			print'(4i6)', bodyp%LOCMAT
			print*
			print*, 'neighbour matrix'
			print'(4i6)', bodyp%NBMAT
			print*
			print*, 'control points'
			print'(3f18.12)', bodyp%CPXYZ
			print*
			print*, 'control points'
			print'(3f18.12)', bodyp%NORMAL
			print*
			print*, 'Velocities'
			print'(3f18.12)', bodyp%CPUVW
			print*

			select type ( bodyp )
			class is ( LiftingSurfaceType )
				print*, 'TELSSI'
				print'(i4)', bodyp%TELSNI
				print*, 'TELSPI'
				print'(i4)', bodyp%TELSSI
				print*, 'TEWAKEI'
				print'(i4)', bodyp%wake%TEWAKEI
				print*, 'Wake nodes'
				print'(3f18.12)', bodyp%wake%NODESXYZ
				print*
				print*, 'Wake Loc mat'
				print'(4i6)', bodyp%wake%LOCMAT
				print*
				print*, 'Wake NBMAT'
				print'(4i6)', bodyp%wake%NBMAT
				print*
			end select
		enddo

	end subroutine

	!subroutine printWake
	!	integer :: i

	!	do i=1,ng
	!		bodyp => geom_get_data( geom, i )
	!		print*, 'Wake for timestep', t
	!		print*
	!		select type (bodyp)
	!		class is ( LiftingSurfaceType )
	!		print'(3f12.6)', bodyp%wake%NODESXYZ( :,1:bodyp%NEN*(t+1) )
	!		end select
	!		print*
	!	enddo
	!end subroutine

	subroutine printRHS

		print*, 'RHS FOR TIMESTEP', t
		print*
		print'(1f18.12)', RHS
		print*
	end subroutine

	!subroutine printWakeG
	!	integer :: i
	!	print*, ' WAKE CIRCULATIONS', t
	!	print*
	!	do i=1,ng
	!	bodyp => geom_get_data( geom, i )
	!	select type( bodyp )
	!	class is ( LiftingSurfaceType )
	!	print'(1f18.12)', bodyp%wake%LOOPG
	!	print*
	!	end select
	!	enddo

	!end subroutine

	subroutine printResult
		integer :: i
		print*, ' RESULT FOR TIMESTEP', t
		print*
		do i=1,ng
		bodyp => geom_get_data( geom, i )
		print'(1f18.12)', bodyp%NEWG
		print*
		enddo

	end subroutine

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

end program Esfera
