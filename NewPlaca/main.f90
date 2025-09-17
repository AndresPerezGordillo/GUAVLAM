program scorpion

	use uvlm_wrappers
	use solvers
	use case_module

      	implicit none

	type( SurfaceType ), dimension(:), allocatable, target :: GEOM
	integer, dimension(:,:), allocatable :: geom_index

	real(8), dimension(:,:), allocatable, target :: A_mat, A_copy, RHS
	integer, dimension(:), allocatable :: ipiv
	integer :: info

      	real(8), dimension(3) :: V_inf
	real(8), dimension(:,:), allocatable :: force, moment
	real(8)   :: delta_t, core
      	integer   :: globnp
	integer   :: t_steps, t_start
	integer   :: t, ng, i
	integer   :: stat_var
	integer   :: num_zones

	character(len=100) :: casefile
	character(len=100) :: tecfile
	character(len=100) :: numchar
	character(len=100) :: date, time

	character(len=100), dimension(:), allocatable :: zones

	real(8) :: t0, t1

	write(*,*) 'Load case file:'
	casefile = 'CaseControl'
	write(*,*) casefile

	call date_and_time(date=date, time=time)
	write(numchar,*) time(1:6)//'_'//date(1:4)//'-'//date(5:6)//'-'//date(7:8)

	call initCase( casefile, ng, num_zones, delta_t, t_steps, t_start, V_inf, core, tecfile )

	allocate( geom(ng) )
	allocate( zones( num_zones ) )

	call initGeom( ng, geom, zones, t_steps, t_start )

	!!---------------
	!! Setting up linear system por T = 0

	allocate( geom_index( 2, ng ) )
	call generateGeomIndex( geom, geom_index )
	globnp = geom_total_np( geom )

	call allocateLinSystem

	A_mat = 0.0d+0
	RHS   = 0.0d+0

	!!--------------
	!! Starting time advance

	call iniTec360File( trim(adjustl(tecfile))//'_'//trim(adjustl(numchar)),&
		'X Y Z Nx Ny Nz Cp', 1 )

	allocate( force( 3, t_steps-t_start + 1 ) )
	allocate( moment( 3, t_steps-t_start + 1 ) )
	
        t0 = omp_get_wtime() 
	!call cpu_time( t0 )

	open( unit=11, file='LOADS.DAT')
	open( unit=12, file='BODYG.DAT')
	open( unit=13, file='GAMMA.DAT')

	timeAdvance: do t = t_start,t_steps

		SolTime = t * delta_t

		print*, '--------------------'
		print*, ' TIMESTEP', t
		print*, '--------------------'
		print*

		!-------------------
		! convecting Wake and moving geometry

   		call GEOMCONVECT( geom, V_inf, delta_t, t, core )

		!!-----------
		! Setting up linear system

		call GEOMAI( geom, geom_index, A_mat, t, t==t_start ) 

		A_copy = A_mat

		RHS = 0.d+0
		
		call GEOMRHS( geom, geom_index, V_inf, t, core, RHS )

		!--------------
		! Solving Linear System

		call DGESV( globnp, 1, A_copy, globnp, ipiv, RHS, globnp, info )

		!--------------
		! Assingning Loop circulations on surface

		call GEOMG( geom, geom_index, RHS )

		!-----------
		! Calculating induced mean velocities

		call GEOMMVEL( geom, V_inf, t, core )

		!---------------
		! Calculating Cp con control points

		call GEOMLOADS( geom, delta_t )

		write(11,*) 'Timestep', t
		write(11,*) '--------------------'
		write(11,'(1f18.12)') geom(1)%LOAD
		write(11,*)

		write(12,*) 'Timestep', t
		write(12,*) '--------------------'
		write(12,'(1f18.12)') geom(1)%NEWG
		write(12,*)
		write(12,'(1f18.12)') geom(1)%OLDG
		write(12,*)

		write(13,*) 'Timestep', t
		write(13,*) '--------------------'
		write(13,'(4f18.12)') geom(1)%GAMMAS
		write(13,*)

		!---------------
		! Writing results to TECPLOT .plt file
		
		call GEOMTECIO( geom, zones, t )  

		call saveForceAndMoment( geom, ng, force(:,t-t_start+1), moment(:,t-t_start+1) )
		
	end do timeAdvance

	close( 11 )
	close( 12 )
	close( 13 )
	
	t1 = omp_get_wtime()
	!call cpu_time( t1 )

	print '("Elapsed time:",f20.2," seconds.")', t1 - t0

	i = TecEnd112()
	print*, 'Tecplot file created with status:', i

	print*, 'Writing forces and moments'
	call writeForceAndMoment( force, moment, numchar )

	print*, 'Writing Case...'
	call writeCase( casefile, ng, geom, zones, delta_t, t_steps, V_inf, core, tecfile )
	print*, 'Case written'

	do i=1,ng
		print*,'Destroying geom', i
		call destroy_surface( geom(i), stat_var )
	enddo
	
	stat_var = 0
	!print*,1
	if( stat_var == 0 ) deallocate( GEOM, STAT=stat_var )
	!print*,2, stat_Var
	if( stat_var == 0 )deallocate( geom_index, STAT=stat_var)
	!print*,3, stat_Var
	if( stat_var == 0 )deallocate( A_mat, STAT=stat_var)
	!print*,4, stat_Var
	if( stat_var == 0 )deallocate( A_copy, STAT=stat_var)
	!print*,5, stat_Var
	if( stat_var == 0 )deallocate( RHS, STAT=stat_var)
	!print*,6, stat_Var
	if( stat_var == 0 )deallocate( ipiv, STAT=stat_var)
	!print*, 'fin'

	contains

	subroutine writeForceAndMoment( force_array, moment_array, num_char )
		real(8), dimension(:,:), intent(inout) :: force_array, moment_array
		character(len=*), intent(in) :: num_char

		open( unit=15, file='FORCE_'//trim(adjustl(num_char))//'.DAT' )
		write(15,'(3f20.12)') force_array
		close(15)

		open( unit=15, file='MOMENT_'//trim(adjustl(num_char))//'.DAT' )
		write(15,'(3f20.12)') moment_array
		close(15)

	end subroutine writeForceAndMoment

	subroutine saveForceAndMoment( geom, ng, force_array, moment_array )
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: bodyp
		integer, intent(in) :: ng
		real(8), dimension(3), intent(inout) :: force_array, moment_array
		real(8), dimension(3) :: force_loc, force_temp, moment_loc, pos
		integer :: i, ii

		force_loc = 0.0d+0
		moment_loc = 0.0d+0
		
		do i = 1,ng

			bodyp => geom(i)

			if( .not. bodyp%BLUFF )then

			do ii = 1,bodyp%NP
				pos = bodyp%CPXYZ(:,ii)
				force_temp = bodyp%NORMAL(:,ii) * bodyp%LOAD(ii,1) * bodyp%AREA(ii,1)
				force_loc = force_temp + force_loc
				moment_loc = cross( pos, force_temp ) &
						+ moment_loc
			enddo

			endif

		enddo

		force_array = force_loc
		moment_array = moment_loc

	end subroutine

    	subroutine allocateLinSystem

		allocate( A_mat  ( globnp,globnp ) )
		allocate( A_copy ( globnp,globnp ) )
		allocate( RHS    ( globnp, 1     ) )
		allocate( ipiv   ( globnp	 ) )

	end subroutine allocateLinSystem

	!subroutine deallocateLinSystem

	!	deallocate( A_mat  )
	!	deallocate( A_copy )
	!	deallocate( RHS    )

	!end subroutine deallocateLinSystem 

end program scorpion
