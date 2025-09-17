program crearRotor

	use surface_module
	use rectangular_blade_module
	use elliptic_module
	
	implicit none

	type( SurfaceType ) :: geom
	integer :: t_start, timesteps
	real(8), dimension(3) :: zeros
	real(8), dimension(3) :: origin
	real(8), dimension(3,1) :: origin2
	real(8), dimension(3) :: eangles
	real(8), dimension(3) :: ea
	integer :: direction
	integer :: bla
	integer :: Viscous! Variable de seleccion para hacer o no correccion Viscosa 1 si
	real(8) :: delta_t
	real(8) :: azimuth
	real(8) :: KnecViscosity! Variable de Viscosidad Cinematica
	real(8) :: Alpha! Variable de constante alfa para correccion Viscosa
	real(8) :: a1! Variable de constante de cambio para calculo de Turbulent Viscosity Coefficient
	real(8) :: coreVis! Variable de constante de cambio para calculo de Turbulent Viscosity Coefficient

	integer :: num_props
	integer :: num_blades
	integer :: num_hub
	integer :: num_fuselage
	real(8) :: rotor_radius
	real(8) :: blade_span
	real(8) :: blade_chord
	real(8) :: blade_tip_chord
	real(8) :: blade_twist
	real(8) :: hub_radius
	real(8) :: gap_frac
	real(8) :: major_rad
	real(8) :: min_rad
	real(8) :: length
	real(8) :: width
	real    :: revs
	real(8) :: delta_psi
	integer :: slow_start
	real(8) :: alfa
	real(8) :: precone
	real(8) :: core
	!real(8) :: body_core
	integer :: blade_ns
	integer :: blade_nc
	integer :: nmaj 
	integer :: nmin 

	real(8) :: dt_aux, ds_aux

	character(len=100) :: geomFile, tecFile
	character(len=100) :: numchar, numaux
	integer :: io, i, j

	zeros=0.d+0

	open(unit=20, file='PARAMETERS.PAR', status='OLD', action='READ')

	read(20,*) num_props
	read(20,*) num_blades
	read(20,*) num_hub
	read(20,*) num_fuselage
	read(20,*) rotor_radius
	read(20,*) blade_chord
	read(20,*) blade_tip_chord
	read(20,*) blade_twist
	read(20,*) alfa
	read(20,*) precone
	read(20,*) ea(1)
	read(20,*) ea(2)
	read(20,*) ea(3)
	read(20,*) hub_radius
	read(20,*) gap_frac
	read(20,*) major_rad
	read(20,*) min_rad
	read(20,*) length
	read(20,*) width
	read(20,*) revs
	read(20,*) delta_psi
	read(20,*) slow_start
	read(20,*) core
	!read(20,*) body_core
	read(20,*) blade_ns
	read(20,*) blade_nc
	read(20,*) nmaj 
	read(20,*) nmin 

	close(20)

	t_start = 0
	timesteps = int((360.0 / delta_psi)*revs) + slow_start
	delta_t = 2.d+0 * pi * (delta_psi/3.60d+02)
	blade_span = rotor_radius - hub_radius
	dt_aux = pi / (blade_ns+2)
	ds_aux = 0.25d+0 * blade_span * (1.d+0 - cos( dt_aux ))
	dt_aux = 0.25d+0 * blade_chord / blade_nc
	core = min(core,ds_aux,dt_aux)

	Viscous = 1
	KnecViscosity = (7.4679d-7)
	Alpha = 1.25643d+0
	a1 = 2d-4
	coreVis = 0.05 * blade_tip_chord
	
	print*, core
	write(numchar,'(3i1)') int(ea(1))/10,int(abs(ea(2)))/10,int(ea(3))/10
	write(numaux,*) num_blades

! Inicio Genera el archivo de control de Viscosidad
	open(unit=17, file='ViscousControl.txt')
	write(17,*) Viscous
	write(17,*) KnecViscosity
	write(17,*) Alpha
	write(17,*) a1
	write(17,*) coreVis
	close(17)
! Fin Genera el archivo de control de Viscosidad

	open(unit=10, file='CaseControl')

	write(10,*) num_props * (num_blades + num_hub) + num_fuselage
	write(10,*) num_props * (2*num_blades + num_hub ) + num_fuselage
	write(10,*) num_props
	write(10,*) timesteps
	write(10,*) t_start
	write(10,'(D)') delta_t
	write(10,'(D)') core
	write(10,'(D)') zeros
	write(10,*) 'rotor_'//trim(adjustl(numaux))//'b'//trim(adjustl(numchar))

	createprops: do j=1,num_props

		direction = (-1)**(j+1)
		
		origin = [0.d+0, direction * min(num_props-1,1) * 1.5d+0 , 0.d+0]
		origin2(:,1) = origin
		eangles = [direction*ea(1), direction*ea(2), ea(3)]

		call eulerRot( 3,1,2, eangles(1),eangles(2),eangles(3),origin2 )

		origin = origin2(:,1)

		createBlades: do i=1,num_blades

	
			write(numchar,*) j*1000 + i
			geomFile = 'blade'//trim(adjustl(numchar))
			write(10,*) trim(adjustl(geomFile))//'.geom'
	
			azimuth = dble(i)* (360.d+0 / dble(num_blades))
		
			geom = initRectBlade( blade_span, blade_chord, blade_tip_chord, core, blade_ns, blade_nc, &
					direction, alfa, blade_twist, precone, azimuth, &
					hub_radius, timesteps, eangles, origin, zeros, direction * [0.d+0, 0.d+0, 1.d+0], slow_start )
	
			open(unit=11, file=trim(adjustl(geomFile))//'.geom' )
	
			write(11,*) geom%NP
			write(11,*) geom%NN
			write(11,*) geom%NVS
			
			write(11,*)'COORDSYS'
			write(11,'(D26.18)') geom%COORDSYS
			write(11,*)'ORIGIN'
			write(11,'(D26.18)') geom%ORIGIN
			write(11,*)'VEL'
			write(11,'(D26.18)') geom%VEL
			write(11,*)'WVEL'
			write(11,'(D26.18)') geom%WVEL
			write(11,*)'EANGLES'
			write(11,'(D26.18)') geom%EANGLES
			write(11,*)'SLOWS'
			write(11,*) geom%SLOWS
			write(11,*)'GROUP'
			write(11,*) j
			write(11,*)'BLUFF'
			write(11,*) geom%BLUFF
	
			write(11,*) 'ZONES'
			write(11,*) 2
			write(11,*) trim(adjustl(geomFile))
			write(11,*) trim(adjustl(geomFile))//'_wake'
			write(11,*) 'NODESXYZ'
			write(11,'(D26.18)') geom%NODESXYZ 
			write(11,*) 'NODESLOC'
			write(11,'(D26.18)') geom%NODESLOC 
			write(11,*) 'CPXYZ'
			write(11,'(D26.18)')geom%CPXYZ    
			write(11,*) 'CPUVW'
			write(11,'(D26.18)')geom%CPUVW    
			write(11,*) 'NORMAL'
			write(11,'(D26.18)')geom%NORMAL   
			write(11,*) 'MVEL'
			write(11,'(D26.18)')geom%MVEL     
			write(11,*) 'IVEL'
			write(11,'(D26.18)')geom%IVEL     
			write(11,*) 'GAMMAS'
			write(11,'(D26.18)')geom%GAMMAS   
			write(11,*) 'VORTEX'
			write(11,'(D26.18)')geom%VORTEX   
			write(11,*) 'NEWG'
			write(11,'(D26.18)')geom%NEWG     
			write(11,*) 'OLDG'
			write(11,'(D26.18)')geom%OLDG     
			write(11,*) 'LOAD'
			write(11,'(D26.18)')geom%LOAD     
			write(11,*) 'AREA'
			write(11,'(D26.18)')geom%AREA     
			write(11,*) 'LOCMAT'
			write(11,*)geom%LOCMAT   
			write(11,*) 'NBMAT'
			write(11,*) geom%NBMAT    
			write(11,*) 'WAKE'
			write(11,*) trim(adjustl(geomFile))//'.wake'
	
			close(11)
	
			open(unit=21, file=trim(adjustl(geomFile))//'.wake')
			write(21,*) geom%wake%START
			write(21,*) geom%wake%NEN
			write(21,*) geom%wake%NEP 
			write(21,*) -1
			write(21,*) size(geom%wake%CORNERS)
			write(21,*)'TEWAKEI'
			write(21,'(1i6)') geom%wake%TEWAKEI
			write(21,*)'CORNERS'
			write(21,'(1i6)') geom%wake%CORNERS
			write(21,*)'CORMAP'
			write(21,'(1i6)') geom%wake%CORMAP
			write(21,*)'TELSSI'
			write(21,'(1i6)') geom%wake%TELSSI
			write(21,*)'TELSNI'
			write(21,'(1i6)') geom%wake%TELSNI
			write(21,*)'ROOT'
			write(21,'(L)') geom%wake%ROOT
			close(21)
	
		enddo createBlades
	
	!-------------------------------------------------------------------------------------------
	
	!	createDuct: do i=1,num_shroud
	!
	!		write(numchar,*) j*1000 + i
	!		geomFile = 'shroud'//trim(adjustl(numchar))
	!		write(10,*) trim(adjustl(geomFile))//'.geom'
	!
	!		shroud_radius = rotor_radius * ( 1.d+0 + gap_frac )
	!		shroud_chord = shroud_radius / shroud_AR
	!
	!		geom = initShroud( origin + [0.d+0,0.d+0,shroud_pos * shroud_chord],&
	!			zeros, zeros, zeros, shroud_chord, shroud_radius, shroud_nc, shroud_nt, 0 )
	!
	!		open(unit=11, file=trim(adjustl(geomFile))//'.geom' )
	!
	!		write(11,*) geom%NP
	!		write(11,*) geom%NN
	!		write(11,*) geom%NVS
	!		
	!		write(11,*)'COORDSYS'
	!		write(11,'(D26.18)') geom%COORDSYS
	!		write(11,*)'ORIGIN'
	!		write(11,'(D26.18)') geom%ORIGIN
	!		write(11,*)'VEL'
	!		write(11,'(D26.18)') geom%VEL
	!		write(11,*)'WVEL'
	!		write(11,'(D26.18)') geom%WVEL
	!		write(11,*)'EANGLES'
	!		write(11,'(D26.18)') geom%EANGLES
	!		write(11,*)'SLOWS'
	!		write(11,*) geom%SLOWS
	!		write(11,*)'BLUFF'
	!		write(11,*) geom%BLUFF
	!
	!		write(11,*) 'ZONES'
	!		!write(11,*) 1
	!		write(11,*) 2
	!		write(11,*) trim(adjustl(geomFile))
	!		write(11,*) trim(adjustl(geomFile))//'_wake'
	!		write(11,*) 'NODESXYZ'
	!		write(11,'(D26.18)') geom%NODESXYZ 
	!		write(11,*) 'NODESLOC'
	!		write(11,'(D26.18)') geom%NODESLOC 
	!		write(11,*) 'CPXYZ'
	!		write(11,'(D26.18)')geom%CPXYZ    
	!		write(11,*) 'CPUVW'
	!		write(11,'(D26.18)')geom%CPUVW    
	!		write(11,*) 'NORMAL'
	!		write(11,'(D26.18)')geom%NORMAL   
	!		write(11,*) 'MVEL'
	!		write(11,'(D26.18)')geom%MVEL     
	!		write(11,*) 'IVEL'
	!		write(11,'(D26.18)')geom%IVEL     
	!		write(11,*) 'GAMMAS'
	!		write(11,'(D26.18)')geom%GAMMAS   
	!		write(11,*) 'VORTEX'
	!		write(11,'(D26.18)')geom%VORTEX   
	!		write(11,*) 'NEWG'
	!		write(11,'(D26.18)')geom%NEWG     
	!		write(11,*) 'OLDG'
	!		write(11,'(D26.18)')geom%OLDG     
	!		write(11,*) 'LOAD'
	!		write(11,'(D26.18)')geom%LOAD     
	!		write(11,*) 'AREA'
	!		write(11,'(D26.18)')geom%AREA     
	!		write(11,*) 'LOCMAT'
	!		write(11,*) geom%LOCMAT   
	!		write(11,*) 'NBMAT'
	!		write(11,*) geom%NBMAT    
	!		write(11,*) 'WAKE'
	!		write(11,*) trim(adjustl(geomFile))//'.wake'
	!
	!		close(11)
	!
	!		open(unit=21, file=trim(adjustl(geomFile))//'.wake')
	!		write(21,*) geom%wake%START
	!		write(21,*) geom%wake%NEN
	!		write(21,*) geom%wake%NEP 
	!		write(21,*) -1
	!		write(21,*) size(geom%wake%CORNERS)
	!		write(21,*)'TEWAKEI'
	!		write(21,'(1i6)') geom%wake%TEWAKEI
	!		write(21,*)'CORNERS'
	!		write(21,'(1i6)') geom%wake%CORNERS
	!		write(21,*)'CORMAP'
	!		write(21,'(1i6)') geom%wake%CORMAP
	!		write(21,*)'TELSSI'
	!		write(21,'(1i6)') geom%wake%TELSSI
	!		write(21,*)'TELSNI'
	!		write(21,'(1i6)') geom%wake%TELSNI
	!		write(21,*)'ROOT'
	!		write(21,'(L)') geom%wake%ROOT
	!		close(21)
	!
	!	enddo createDuct
	
	!-------------------------------------------------------------------------------------------
	
		createHub: do i=1,num_hub
	
			write(numchar,*) j*1000 + i
			geomFile = 'hub'//trim(adjustl(numchar))
			write(10,*) trim(adjustl(geomFile))//'.geom'
	
			geom = initEllipticGeom( major_rad, nmaj, min_rad, nmin, 3, zeros, origin, zeros, zeros )
	
			open(unit=11, file=trim(adjustl(geomFile))//'.geom' )
	
			write(11,*) geom%NP
			write(11,*) geom%NN
			write(11,*) geom%NVS
			
			write(11,*)'COORDSYS'
			write(11,'(D26.18)') geom%COORDSYS
			write(11,*)'ORIGIN'
			write(11,'(D26.18)') geom%ORIGIN
			write(11,*)'VEL'
			write(11,'(D26.18)') geom%VEL
			write(11,*)'WVEL'
			write(11,'(D26.18)') geom%WVEL
			write(11,*)'EANGLES'
			write(11,'(D26.18)') geom%EANGLES
			write(11,*)'SLOWS'
			write(11,*) geom%SLOWS
			write(11,*)'BLUFF'
			write(11,*) geom%BLUFF
	
			write(11,*) 'ZONES'
			write(11,*) 1
			write(11,*) trim(adjustl(geomFile))
			write(11,*) 'NODESXYZ'
			write(11,'(D26.18)') geom%NODESXYZ 
			write(11,*) 'NODESLOC'
			write(11,'(D26.18)') geom%NODESLOC 
			write(11,*) 'CPXYZ'
			write(11,'(D26.18)')geom%CPXYZ    
			write(11,*) 'CPUVW'
			write(11,'(D26.18)')geom%CPUVW    
			write(11,*) 'NORMAL'
			write(11,'(D26.18)')geom%NORMAL   
			write(11,*) 'MVEL'
			write(11,'(D26.18)')geom%MVEL     
			write(11,*) 'IVEL'
			write(11,'(D26.18)')geom%IVEL     
			write(11,*) 'GAMMAS'
			write(11,'(D26.18)')geom%GAMMAS   
			write(11,*) 'VORTEX'
			write(11,'(D26.18)')geom%VORTEX   
			write(11,*) 'NEWG'
			write(11,'(D26.18)')geom%NEWG     
			write(11,*) 'OLDG'
			write(11,'(D26.18)')geom%OLDG     
			write(11,*) 'LOAD'
			write(11,'(D26.18)')geom%LOAD     
			write(11,*) 'AREA'
			write(11,'(D26.18)')geom%AREA     
			write(11,*) 'LOCMAT'
			write(11,*)geom%LOCMAT   
			write(11,*) 'NBMAT'
			write(11,*) geom%NBMAT    
	
			close(11)
	
		enddo createHub
	
	enddo createprops

	if( num_fuselage == 1 ) then
		
			geomFile = 'fuselage'
			write(10,*) trim(adjustl(geomFile))//'.geom'
	
			geom = initEllipticGeom( length, nmaj, width, nmin, 1, zeros, &
						[0.d+0,0.d+0,-0.5d+0], zeros, zeros )
	
			open(unit=11, file=trim(adjustl(geomFile))//'.geom' )
	
			write(11,*) geom%NP
			write(11,*) geom%NN
			write(11,*) geom%NVS
			
			write(11,*)'COORDSYS'
			write(11,'(D26.18)') geom%COORDSYS
			write(11,*)'ORIGIN'
			write(11,'(D26.18)') geom%ORIGIN
			write(11,*)'VEL'
			write(11,'(D26.18)') geom%VEL
			write(11,*)'WVEL'
			write(11,'(D26.18)') geom%WVEL
			write(11,*)'EANGLES'
			write(11,'(D26.18)') geom%EANGLES
			write(11,*)'SLOWS'
			write(11,*) geom%SLOWS
			write(11,*)'BLUFF'
			write(11,*) geom%BLUFF
	
			write(11,*) 'ZONES'
			write(11,*) 1
			write(11,*) trim(adjustl(geomFile))
			write(11,*) 'NODESXYZ'
			write(11,'(D26.18)') geom%NODESXYZ 
			write(11,*) 'NODESLOC'
			write(11,'(D26.18)') geom%NODESLOC 
			write(11,*) 'CPXYZ'
			write(11,'(D26.18)')geom%CPXYZ    
			write(11,*) 'CPUVW'
			write(11,'(D26.18)')geom%CPUVW    
			write(11,*) 'NORMAL'
			write(11,'(D26.18)')geom%NORMAL   
			write(11,*) 'MVEL'
			write(11,'(D26.18)')geom%MVEL     
			write(11,*) 'IVEL'
			write(11,'(D26.18)')geom%IVEL     
			write(11,*) 'GAMMAS'
			write(11,'(D26.18)')geom%GAMMAS   
			write(11,*) 'VORTEX'
			write(11,'(D26.18)')geom%VORTEX   
			write(11,*) 'NEWG'
			write(11,'(D26.18)')geom%NEWG     
			write(11,*) 'OLDG'
			write(11,'(D26.18)')geom%OLDG     
			write(11,*) 'LOAD'
			write(11,'(D26.18)')geom%LOAD     
			write(11,*) 'AREA'
			write(11,'(D26.18)')geom%AREA     
			write(11,*) 'LOCMAT'
			write(11,*)geom%LOCMAT   
			write(11,*) 'NBMAT'
			write(11,*) geom%NBMAT    

			close(11)
	endif

end program
