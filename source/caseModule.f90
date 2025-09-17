module case_module
!---------------------------------------------------------------------
!	caseModule.f90
!-------------------------------------------------------------
!	
!       This module contains the subroutines used to initiate and
!	terminate a case run by the GUAVLAM code.
!
!	Developed  by : Juan D. Colmenares
!	Department of Mechanical Engineering
!	Universidad de los Andes
!
!---------------------------------------------------------------------

	use surface_module

	implicit none

	contains

	subroutine initCase( casefile, ng, num_zones, num_groups, delta_t, timesteps, t_start, V_inf, core, tecfile )
	!
	! Load Case file
	!
		real(8), dimension(3), intent(out) :: V_inf
		real(8), intent(out) :: delta_T, core
		integer, intent(out) :: t_start, timesteps
		integer, intent(out) :: ng, num_zones, num_groups
		integer :: io, i
		character(len=*), intent(out):: tecfile
		character(len=*), intent(in) :: casefile

		open(unit=10, file=trim(adjustl(casefile)), status='OLD', action='READ')

		read(10,*,iostat=io) ng
		read(10,*,iostat=io) num_zones
		read(10,*,iostat=io) num_groups
		read(10,*,iostat=io) timesteps
		read(10,*,iostat=io) t_start
		read(10,'(D)',iostat=io) delta_t
		read(10,'(D)',iostat=io) core
		read(10,'(D)',iostat=io) V_inf
		read(10,'(A)',iostat=io) tecfile

		if( timesteps < t_start )then
			write(*,*) 'Incorrect timesteps'
			stop
		endif

	end subroutine initCase


	subroutine initGeom( ng, geom, zones, timesteps, t_start )
	!
	! Load data for a single geometry (surface) from a .geom file.
	!
		type( SurfaceType ), dimension(:), intent(inout) :: geom
		integer, intent(in) :: t_start, timesteps, ng
		integer :: i
		integer :: io, zone_counter
		character(len=100) :: geomfile
		character(len=100), dimension(:), intent(inout) :: zones

		zone_counter = 1
		io = 0

		do i=1,ng
			read(10,'(A)',iostat=io) geomfile
			if( io == 0 ) geom(i) = loadGeom(geomfile, zones, zone_counter, t_start, timesteps)
		enddo

		close(10)

	end subroutine initGeom


	function loadGeom( geomfile, zones, zone_counter, t_start, timesteps ) result( geom )
	!
	!  Allocate data from the .geom file to the components of a SurfaceType
	!  variable.
	!
		character(len=*), intent(in) :: geomfile
		character(len=*), dimension(:), intent(inout) :: zones
		character(len=100) :: variable, wakeFile
		integer, intent(inout) :: zone_counter
		integer, intent(in) :: t_start, timesteps
		integer :: gstat, num_zones
		integer :: i
		type( SurfaceType ) :: geom

		print*, 'Opening geometry file: ', trim(adjustl(geomfile))
		open(unit=11, file=trim(adjustl(geomfile)), status='OLD', action='READ')

		gstat = 0

		!This data is read in a predefined order
		read(11,*) geom%NP
		read(11,*) geom%NN
		read(11,*) geom%NVS

		allocate( geom%NODESXYZ ( 3,geom%NN  ) )
		allocate( geom%NODESLOC ( 3,geom%NN  ) )
		allocate( geom%CPXYZ    ( 3,geom%NP  ) )
		allocate( geom%CPUVW    ( 3,geom%NP  ) )
		allocate( geom%NORMAL   ( 3,geom%NP  ) )
		allocate( geom%MVEL     ( 3,geom%NP  ) )
		allocate( geom%IVEL     ( 3,geom%NP  ) )
		allocate( geom%GAMMAS   ( 4,geom%NP  ) )
		allocate( geom%VORTEX   ( 7,geom%NVS ) )
		allocate( geom%NEWG     ( geom%NP,1  ) )
		allocate( geom%OLDG     ( geom%NP,1  ) )
		allocate( geom%LOAD     ( geom%NP,1  ) )
		allocate( geom%AREA     ( geom%NP,1  ) )
		allocate( geom%LOCMAT   ( 4,geom%NP  ) )
		allocate( geom%NBMAT    ( 4,geom%NP  ) )

		geom%MVEL  = 0.d+0
		geom%IVEL  = 0.d+0
		geom%GAMMAS= 0.d+0
		geom%VORTEX= 0.d+0
		geom%NEWG  = 0.d+0
		geom%OLDG  = 0.d+0
		geom%LOAD  = 0.d+0

		! This data may be read in a random order.
		readGeom: do while( gstat == 0 )
			read(11,'(A)',iostat=gstat) variable
			if( gstat == 0 )then
				if( index(variable,'!') /= 0 ) cycle readGeom
				print*, 'Reading variable ', trim(adjustl(variable))
				select case( trim(adjustl(variable)) )
				case('ZONES')
					read(11,'(i)',iostat=gstat)&
						num_zones
					do i=1,num_zones
						read(11,'(A)',iostat=gstat) zones(zone_counter)
						zone_counter = zone_counter + 1
					enddo
				case('COORDSYS')
					read(11,'(D)',iostat=gstat)	geom%COORDSYS
				case('ORIGIN')
					read(11,'(D)',iostat=gstat) geom%ORIGIN
				case('VEL')
					read(11,'(D)',iostat=gstat) geom%VEL
				case('WVEL')
					read(11,'(D)',iostat=gstat) geom%WVEL
				case('EANGLES')
					read(11,'(D)',iostat=gstat) geom%EANGLES
				case('NODESXYZ')
					read(11,'(D)',iostat=gstat) geom%NODESXYZ
				case('NODESLOC')
					read(11,'(D)',iostat=gstat) geom%NODESLOC
				case('CPXYZ')
					read(11,'(D)',iostat=gstat) geom%CPXYZ
				case('CPUVW')
					read(11,'(D)',iostat=gstat) geom%CPUVW
				case('NORMAL')
					read(11,'(D)',iostat=gstat) geom%NORMAL
				case('AREA')
					read(11,'(D)',iostat=gstat) geom%AREA
				case('MVEL')
					read(11,'(D)',iostat=gstat) geom%MVEL
				case('IVEL')
					read(11,'(D)',iostat=gstat) geom%IVEL
				case('GAMMAS')
					read(11,'(D)',iostat=gstat) geom%GAMMAS
				case('VORTEX')
					read(11,'(D)',iostat=gstat) geom%VORTEX
				case('NEWG')
					read(11,'(D)',iostat=gstat) geom%NEWG
				case('OLDG')
					read(11,'(D)',iostat=gstat) geom%OLDG
				case('LOAD')
					read(11,'(D)',iostat=gstat) geom%LOAD
				case('LOCMAT')
					read(11,*,iostat=gstat) geom%LOCMAT
				case('NBMAT')
					read(11,*,iostat=gstat) geom%NBMAT
				case('SLOWS')
					read(11,*,iostat=gstat) geom%SLOWS
				case('BLUFF')
					read(11,'(L)',iostat=gstat) geom%BLUFF
				case('GROUP')
					read(11,*,iostat=gstat) geom%GROUP
				case('WAKE')
					allocate( geom%wake )
					read(11,'(A)',iostat=gstat) wakeFile
					if( gstat == 0 ) call loadWake( geom%wake, wakeFile, t_start, timesteps )
					geom%wake%NODESXYZ( :,1:geom%wake%NEN ) = geom%NODESXYZ( :, geom%wake%TELSNI )
				end select
			endif
		enddo readGeom

		close(11)

	end function

	subroutine loadWake( wake, wakeFile, t_start_1, timesteps_1 )
	!
	! Load data for a WakeType variable from a given .wake input file.
	!
		type( WakeType ), intent(inout) :: wake
		character(len=*), intent(in) :: wakeFile
		integer, intent(in) :: timesteps_1, t_start_1
		integer :: timesteps, t_start
		integer :: wstat, i, wakeNodes, wakePanelsN, num_corners
		character(len=100) :: variable
        
		open(unit=21,file=trim(adjustl(wakeFile)),status='OLD',action='READ')

		read(21,*,iostat=wstat) wake%START
		read(21,*,iostat=wstat) wake%NEN
		read(21,*,iostat=wstat) wake%NEP 
		read(21,*,iostat=wstat) wake%NVS

		timesteps = timesteps_1 - wake%START
		t_start = t_start_1 - wake%START

		if( wake%NVS == -1 )then
			wake%NVS = wake%NEP * (timesteps+1) + wake%NEN * (timesteps)
			print*, 'wakeNVS', wake%NVS
		endif

		read(21,*,iostat=wstat) num_corners

		wakeNodes = wake%NEN * (timesteps + 1)
		wakePanelsN = wake%NEP * (timesteps)

		allocate( wake%NODESXYZ ( 3,wakeNodes ) )
		allocate( wake%NODESUVW ( 3,wakeNodes - wake%NEN ) )
		allocate( wake%LOOPG    ( wakePanelsN,1 ) )
		allocate( wake%VORTEXG  ( wake%NVS,1 ) )
		allocate( wake%VORTEXP  ( 6,wake%NVS ) )
		allocate( wake%VORTEXC  ( wake%NVS,1 ) )! Espacio para vector de vortex core
		allocate( wake%LOCMAT   ( 4,wakePanelsN ) )
		allocate( wake%VORTEXPI ( 2,wake%NVS ) )
		allocate( wake%TEWAKEI  ( wake%NEP ) )
		allocate( wake%TELSNI   ( wake%NEN ) )
		allocate( wake%TELSSI   ( wake%NEP ) )
		allocate( wake%CORNERS  ( num_corners ) )
		allocate( wake%CORMAP   ( num_corners ) )
		
		wake%NODESXYZ = 0.d+0 
		wake%NODESUVW = 0.d+0 
		wake%LOOPG    = 0.d+0
		wake%VORTEXG = 0.d+0

		call wakeVortexI( wake, timesteps )
		call wakePanels( wake, timesteps )

		readingWake: do while( wstat == 0 )
			read(21,'(A)',iostat=wstat) variable
			if( wstat == 0 )then
				if( index(variable,'!') /= 0 ) cycle readingWake
				print*, 'Reading wake variable ', trim(adjustl(variable))
				select case( trim(adjustl(variable)) )
				case('NODESXYZ')
					read(21,'(D)',iostat=wstat) &
						wake%NODESXYZ(:,1:wake%NEN*(t_start))
				case('LOOPG')
					read(21,'(D)',iostat=wstat)&
						wake%LOOPG(1:wake%NEP*(t_start-1),1)
				case('VORTEXG')
					read(21,'(D)',iostat=wstat)&
						wake%VORTEXG(1:(wake%NEN*(t_start-1)+wake%NEP*(t_start)),1)
				case('VORTEXP')
					read(21,'(D)',iostat=wstat)&
						wake%VORTEXP(:,1:(wake%NEN*(t_start-1)+wake%NEP*(t_start)))
				case('VORTEXC')! Vector VortexC
					read(21,'(D)',iostat=wstat)&
						wake%VORTEXC(1:(wake%NEN*(t_start-1)+wake%NEP*(t_start)),1)! Vector VortexC
				case('TEWAKEI')
					read(21,*,iostat=wstat)&
						wake%TEWAKEI
				case('CORNERS')
					read(21,*,iostat=wstat)&
						wake%CORNERS
				case('CORMAP')
					read(21,*,iostat=wstat)&
						wake%CORMAP
				case('TELSSI')
					read(21,*,iostat=wstat)&
						wake%TELSSI
				case('TELSNI')
					read(21,*,iostat=wstat)&
						wake%TELSNI
				case('ROOT')
					read(21,'(L)',iostat=wstat)&
						wake%ROOT
				end select
			endif
		enddo readingWake

		close(21)

	end subroutine


	subroutine writeCase( casefile, ng, num_groups, geom, zones, delta_t, timesteps, V_inf, core, tecfile )
	!
	! Writes saved data to the corresponding Case file, .geom files and
	! .wake files.
	!
	! NOTE: All pre-existing files that have the same names as the new files
	! will be overwritten. 
	!
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: body
		real(8), dimension(:), intent(in) :: V_inf
		real(8), intent(in) :: delta_T, core
		integer, intent(in) :: timesteps, ng, num_groups
		integer :: i, num_zones
		integer :: io
		integer :: zone_counter
		character(len=*), intent(in):: tecfile
		character(len=*), intent(in)   :: casefile
		character(len=100) :: geomfile
		character(len=100) :: numchar, id
		character(len=100), dimension(:), intent(in) :: zones

		id = tecfile
		
		open(unit=30, file='CaseControl_'//trim(adjustl(id)))

		num_zones = size(zones)
		zone_counter = 1

		write(30,*) ng
		write(30,*) num_zones
		write(30,*) num_groups
		write(30,*) -1
		write(30,*) timesteps+1
		write(30,'(D)') delta_t
		write(30,'(D)') core
		write(30,'(D)') V_inf
		write(30,*) trim(adjustl(tecfile))

		do i=1,ng
			write(numchar,*) i
			geomfile = 'geom'//trim(adjustl(numchar))//'_'//trim(adjustl(id))
			write(30,*) trim(adjustl(geomfile))//'.geom'
			write(*,*) 'Writing to file '//trim(adjustl(geomfile))//'.geom'
			body => geom(i)
			call writeGeom(body, geomfile, zones, zone_counter)
		enddo

		close(30)

	end subroutine writeCase

	subroutine writeGeom( geom, geomfile, zones, zone_counter )
	!
	! Writes data from a SurfaceType variable to a .geom output file.
	!
		type( SurfaceType ), intent(in), pointer :: geom
		character(len=*), intent(in) :: geomfile
		character(len=*), dimension(:), intent(in) :: zones
		integer, intent(inout) :: zone_counter

		open(unit=31, file=trim(adjustl(geomfile))//'.geom')
		
		write(31,*) geom%NP
		write(31,*) geom%NN
		write(31,*) geom%NVS
		
		write(31,*) 'COORDSYS'
		write(31,'(D)') geom%COORDSYS
		write(31,*) 'ORIGIN'
		write(31,'(D)') geom%ORIGIN
		write(31,*) 'VEL'
		write(31,'(D)') geom%VEL
		write(31,*) 'WVEL'
		write(31,'(D)') geom%WVEL
		write(31,*) 'EANGLES'
		write(31,'(D)') geom%EANGLES
		write(31,*) 'NODESXYZ'
		write(31,'(D)') geom%NODESXYZ
		write(31,*) 'NODESLOC'
		write(31,'(D)') geom%NODESLOC
		write(31,*) 'CPXYZ'
		write(31,'(D)') geom%CPXYZ
		write(31,*) 'CPUVW'
		write(31,'(D)') geom%CPUVW
		write(31,*) 'NORMAL'
		write(31,'(D)') geom%NORMAL
		write(31,*) 'AREA'
		write(31,'(D)') geom%AREA
		write(31,*) 'MVEL'
		write(31,'(D)') geom%MVEL
		write(31,*) 'IVEL'
		write(31,'(D)') geom%IVEL
		write(31,*) 'GAMMAS'
		write(31,'(D)') geom%GAMMAS
		write(31,*) 'VORTEX'
		write(31,'(D)') geom%VORTEX
		write(31,*) 'NEWG'
		write(31,'(D)') geom%NEWG
		write(31,*) 'OLDG'
		write(31,'(D)') geom%OLDG
		write(31,*) 'LOAD'
		write(31,'(D)') geom%LOAD
		write(31,*) 'LOCMAT'
		write(31,*) geom%LOCMAT
		write(31,*) 'NBMAT'
		write(31,*) geom%NBMAT
		write(31,*) 'SLOWS'
		write(31,*) geom%SLOWS
		write(31,*) 'GROUP'
		write(31,*) geom%GROUP
		write(31,*) 'BLUFF'
		write(31,'(L)') geom%BLUFF

		if( associated(geom%wake) )then
			call writeWake( geom%wake, geomfile )
			write(31,*) 'ZONES'
			write(31,*) 2
			write(31,*) trim(adjustl(zones(zone_counter)))
			zone_counter= zone_counter + 1
			write(31,*) trim(adjustl(zones(zone_counter)))
			zone_counter= zone_counter + 1
			write(31,*) 'WAKE'
			write(31,*) trim(adjustl(geomfile))//'.wake'
		else
			write(31,*) 'ZONES'
			write(31,*) 1
			write(31,*) trim(adjustl(zones(zone_counter)))
			zone_counter= zone_counter + 1
		endif

		close(31)

	end subroutine writeGeom

	subroutine writeWake( wake, wakeFile )
	!
	! Writes data from a WakeType variable to a .wake output file.
	!
		type( WakeType ), intent(in) :: wake
		character(len=*), intent(in) :: wakeFile
		
		open(unit=32, file=trim(adjustl(wakeFile))//'.wake')
		print*, 'Writing wake to file '//trim(adjustl(wakeFile))//'.wake'
	
		write(32,*) wake%START
		write(32,*) wake%NEN
		write(32,*) wake%NEP 
		write(32,*) -1
		write(32,*) size(wake%CORNERS)

		write(32,*) 'NODESXYZ'
		write(32,'(D)') wake%NODESXYZ
		write(32,*) 'LOOPG'
		write(32,'(D)') wake%LOOPG
		write(32,*) 'VORTEXG'
		write(32,'(D)') wake%VORTEXG
		write(32,*) 'VORTEXP'
		write(32,'(D)') wake%VORTEXP
		write(32,*) 'VORTEXC'
		write(32,'(D)') wake%VORTEXC
		write(32,*) 'TEWAKEI'
		write(32,*) wake%TEWAKEI
		write(32,*) 'CORNERS'
		write(32,*) wake%CORNERS
		write(32,*) 'CORMAP'
		write(32,*) wake%CORMAP
		write(32,*) 'TELSSI'
		write(32,*) wake%TELSSI
		write(32,*) 'TELSNI'
		write(32,*) wake%TELSNI
		write(32,*) 'ROOT'
		write(32,*) wake%ROOT

		close(32)
		
	end subroutine writeWake

end module
