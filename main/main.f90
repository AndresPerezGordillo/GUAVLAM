program GUAVLAM
!---------------------------------------------------------------------
!
!	General Unsteady Aerodynamics Vortex LAttice Method
!------------------------  GUAVLAM 1.0  ------------------------------
!	
!	Developed  by : Juan D. Colmenares
!	Supervised by : Dr. Sergio Preidikman
!			Dr. Omar Lopez
!
!	Department of Mechanical Engineering
!	Universidad de los Andes
!
!---------------------------------------------------------------------

	!!--------------------------------------------------------------
	!! Modules

	use uvlm_wrappers
	use case_module

      	implicit none

	!!--------------------------------------------------------------
	!! Declaration of Variables

	type( SurfaceType ), dimension(:), allocatable, target :: geom
	integer, dimension(:,:), allocatable :: geom_index

	real(8), dimension(:,:), allocatable, target :: a_mat, a_copy, rhs
	integer, dimension(:), allocatable :: ipiv
	integer :: info

      	real(8), dimension(3) :: v_inf
	real(8), dimension(:,:,:), allocatable :: force, moment
	real(8)   :: delta_t, core
	real(8)	:: KnecViscosity, Alpha, a1, coreVis! VariableS de Correccion Viscosa
      	integer   :: globnp
	integer   :: t_steps, t_start
	integer   :: t, ng, i
	integer   :: stat_var
	integer   :: num_zones
	integer   :: num_groups
	integer :: Viscous! VariableS de Correccion Viscosa

	character(len=100) :: casefile
	character(len=100) :: tecfile
	character(len=100) :: numchar
	character(len=100) :: ViscousControl!Variable de texto para cargar archivo de correccion viscosa
	character(len=100) :: ViscousCorrection!Variable de texto para Escribir si o no correccion viscosa

	character(len=100), dimension(:), allocatable :: zones

	real(8) :: t0, t1

	external :: dgesv

	!!--------------------------------------------------------------
	!! Loading Case File & Geometry

	write(*,*) 'Load case file:'
	casefile = 'CaseControl'! Default Case file
	! read(*,*) casefile ! Enter alternate name for a case file.
	write(*,*) casefile

	call initCase( casefile, ng, num_zones, num_groups, delta_t, t_steps, t_start, V_inf, core, tecfile )

	allocate( geom(ng) )
	allocate( zones( num_zones ) )

	call initGeom( ng, geom, zones, t_steps, t_start )

	!!--------------------------------------------------------------
	!! Loading Viscous Parameters

	ViscousControl = 'ViscousControl.txt'
	open(unit=17, file=ViscousControl, status='OLD', action='READ')
	read(17,*) Viscous
	read(17,*) KnecViscosity
	read(17,*) Alpha
	read(17,*) a1
	read(17,*) coreVis
	close(17)
	if( Viscous==1 )then
			ViscousCorrection = 'si'
			!core = coreVis
	else
	ViscousCorrection = 'no'
	end if

	!!--------------------------------------------------------------
	!! Allocating linear system por T = 0

	allocate( geom_index( 2, ng ) )
	call generateGeomIndex( geom, geom_index )
	globnp = geomtotalnp( geom )

	call allocateLinSystem

	A_mat = 0.0d+0
	RHS   = 0.0d+0

	!!--------------------------------------------------------------
	!! Initiating Tecplot File and Force and Moment Arrays

	call iniTec360File( trim(adjustl(tecfile)),&
		'X Y Z Nx Ny Nz Cp', 1 )

	allocate( force( 3, num_groups+1, t_steps-t_start + 1 ) )
	allocate( moment( 3, num_groups+1, t_steps-t_start + 1 ) )
	
        t0 = omp_get_wtime()  ! Use when linked with OMP Library
	!call cpu_time( t0 )   ! Use for Sequential Code 

	!!--------------------------------------------------------------
	!! Starting time advance

	timeAdvance: do t = t_start,t_steps

		SolTime = t * delta_t

		print*, '--------------------'
		print*, ' TIMESTEP', t
		print*, '--------------------'
		print*, ' Viscosity = ', ViscousCorrection
		print*

		!!-------------------
		!! convecting Wake and moving geometry

   		call GEOMCONVECT( geom, V_inf, delta_t, t, core, KnecViscosity, Alpha, a1, Viscous )

		!!-----------
		!! Setting up linear system

		call GEOMAI( geom, geom_index, A_mat, t, t==t_start, core ) 

		A_copy = A_mat

		RHS = 0.d+0
		
		call GEOMRHS( geom, geom_index, V_inf, t, core, RHS, KnecViscosity, Alpha, a1, delta_t, Viscous )

		!!--------------
		!1 Solving Linear System

		call DGESV( globnp, 1, A_copy, globnp, ipiv, RHS, globnp, info )

		!!--------------
		!! Assingning Loop circulations on surface

		call GEOMG( geom, geom_index, RHS )

		!!-----------
		!! Calculating induced mean velocities

		call GEOMMVEL( geom, V_inf, t, core )

		!!---------------
		!! Calculating Cp con control points

		call GEOMLOADS( geom, delta_t )

		!!---------------
		!! Writing results to TECPLOT .plt file
		
		call GEOMTECIO( geom, zones, t )  

		call saveForceAndMoment( geom, ng, force, moment, t-t_start+1 )
		
	end do timeAdvance
	
	t1 = omp_get_wtime()
	!call cpu_time( t1 )

	print '("Elapsed time:",f20.2," seconds.")', t1 - t0

	i = TecEnd112()
	print*, 'Tecplot file created with status:', i

	print*, 'Writing forces and moments'
	call writeForceAndMoment( force, moment, tecfile )

	print*, 'Writing Case...'
	call writeCase( casefile, ng, num_groups, geom, zones, delta_t, t_steps, V_inf, core, tecfile )
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

	subroutine writeForceAndMoment( force_array, moment_array, tecfile )
		real(8), dimension(:,:,:), intent(in) :: force_array, moment_array
		character(len=*), intent(in) :: tecfile
		character(len=100) :: numchar
		integer :: ts, ti

		ts = size(force_array,3)
		write(numchar,*) 3 * size(force_array,2)

		open( unit=15, file='FORCE_'//trim(adjustl(tecfile))//'.DAT' )
		write(15,*)
		write(15,'('//trim(adjustl(numchar))//'f20.12)') (force_array(:,:,ti), ti=1,ts)
		close(15)

		open( unit=15, file='MOMENT_'//trim(adjustl(tecfile))//'.DAT' )
		write(15,*)
		write(15,'('//trim(adjustl(numchar))//'f20.12)') (moment_array(:,:,ti), ti=1,ts)
		close(15)

	end subroutine writeForceAndMoment

	subroutine saveForceAndMoment( geom, ng, force_array, moment_array, t )
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: bodyp
		integer, intent(in) :: ng,t
		real(8), dimension(:,:,:), intent(inout) :: force_array, moment_array
		real(8), dimension(:,:),allocatable :: force_loc, moment_loc
		real(8), dimension(3) :: force_temp, moment_temp, pos
		integer :: i, ii

		allocate( force_loc(3,size(force_array,2)))
		allocate( moment_loc(3,size(force_array,2)))

		force_loc = 0.0d+0
		moment_loc = 0.0d+0

		
		do i = 1,ng

			bodyp => geom(i)

			!if( .not. bodyp%BLUFF )then

			do ii = 1,bodyp%NP
				pos = bodyp%CPXYZ(:,ii)
				force_temp = bodyp%NORMAL(:,ii) * bodyp%LOAD(ii,1) * bodyp%AREA(ii,1)
				force_loc(:,1) = force_temp + force_loc(:,1)
				if( bodyp%group > 0 .and. size(force_array,2) > 1 )then
					force_loc(:,bodyp%group+1) = force_temp + force_loc(:,bodyp%group+1)
				endif
						
				moment_temp = cross( pos, force_temp )
				moment_loc(:,1) = moment_temp + moment_loc(:,1)
				if( bodyp%group > 0 .and. size(force_array,2) > 1)then
					pos = bodyp%CPXYZ(:,ii) - bodyp%ORIGIN
					moment_temp = cross( pos, force_temp )
					moment_loc(:,bodyp%group+1) = moment_temp + moment_loc(:,bodyp%group+1)
				endif
			enddo

			!endif

		enddo

		force_array(:,:,t) = force_loc
		moment_array(:,:,t) = moment_loc

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

end program GUAVLAM
