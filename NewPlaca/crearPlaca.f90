program crearPlaca

    	use rectangular_wing_module
	use surface_module

	implicit none

	type( SurfaceType ) :: geom
	real(8), dimension(3) :: zeros
	real(8) :: span, chord, alfa
	integer :: t_start, timesteps, wake_size
	integer :: ns, nc
	integer :: bla
	real(8) :: bla2

	character(len=100) :: geomFile, tecFile

	open(unit=10, file='CaseControl', action='READ', status='OLD')

	read(10,*) bla
	read(10,*) bla
	read(10,*) timesteps
	read(10,*) t_start
	read(10,*) bla2
	read(10,*) bla2
	read(10,*) bla2
	read(10,*) bla2
	read(10,*) bla2
	read(10,*) tecFile
	read(10,*) geomFile
	
	close(10)

	zeros=0.d+0
	wake_size = timesteps
	ns = 4
	nc = 4
	span = 4.d+0
	chord = 4.d+0
	alfa = 00.d+0

	geom = initRectWing( span, chord, ns, nc, alfa, timesteps, zeros, zeros, zeros, zeros )

	open(unit=11, file=trim(adjustl(geomFile)))

	write(11,*) geom%NP
	write(11,*) geom%NN
	write(11,*) geom%NVS
	write(11,*) geom%CORE
	
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

end program
