module shroud_module
	! 
	! This module declares the shroud type. The readius of the shroud is proportional
	! to the radius of the propeller (center body radius + shroud span). The shape
	! is assumed to be cylindrical (thin, symmetrical airfoil).
	!
	use math_routines
	use surface_module
	use wake_module

  	implicit none
	private

	public initShroud


 	contains

	function initShroud( origin, vel, wvel, euler_angs, chord, radius, nc, nt, wake_size ) result( shroud )
		real(8), dimension(3), intent(in) :: origin, vel, wvel, euler_angs
		real(8), intent(in)		  :: chord, radius
		real(8), dimension(3,1) :: wvel_aux
		integer, intent(in) :: nc, nt, wake_size
		integer :: i
		type( SurfaceType )	  :: shroud

		shroud%ORIGIN = origin
		shroud%VEL    = vel
		shroud%WVEL   = wvel

		shroud%bluff = .false.
		
		shroud%COORDSYS = reshape([1.d+0, 0.d+0, 0.d+0,&
			    		   0.d+0, 1.d+0, 0.d+0,&
			    		   0.d+0, 0.d+0, 1.d+0],[3,3])

		call eulerRot( 3, 1, 2, euler_angs(1), euler_angs(2), euler_angs(3), shroud%COORDSYS )

		shroud%NN = nt * ( nc + 1 )
		shroud%NP = nt * nc
		shroud%NVS = nt*nc + nt*(nc+1)


		allocate( shroud%NODESXYZ ( 3,shroud%NN  ) )
		allocate( shroud%NODESLOC ( 3,shroud%NN  ) )
		allocate( shroud%CPXYZ    ( 3,shroud%NP  ) )
		allocate( shroud%CPUVW    ( 3,shroud%NP  ) )
		allocate( shroud%NORMAL   ( 3,shroud%NP  ) )
		allocate( shroud%MVEL     ( 3,shroud%NP  ) )
		allocate( shroud%IVEL     ( 3,shroud%NP  ) )
		allocate( shroud%GAMMAS   ( 4,shroud%NP  ) )
		allocate( shroud%VORTEX   ( 7,shroud%NVS ) )
		allocate( shroud%NEWG     ( shroud%NP,1  ) )
		allocate( shroud%OLDG     ( shroud%NP,1  ) )
		allocate( shroud%LOAD     ( shroud%NP,1  ) )
		allocate( shroud%AREA     ( shroud%NP,1  ) )
		allocate( shroud%LOCMAT   ( 4,shroud%NP  ) )
		allocate( shroud%NBMAT    ( 4,shroud%NP  ) )


		shroud%NODESLOC = shroudNodes ( chord, radius, nc, nt, shroud%NN )
		shroud%LOCMAT   = shroudPanels( nc, nt, shroud%NP )
		shroud%NBMAT    = snbmatrix   ( nc, nt, shroud%NP )
		shroud%VORTEX   = 0.d+0

		shroud%NODESXYZ = matmul( shroud%COORDSYS, shroud%NODESLOC ) 

		forall( i = 1:shroud%NN )
			shroud%NODESXYZ( :,i ) = shroud%NODESXYZ( :,i ) + shroud%ORIGIN
		end forall

		wvel_aux(:,1) = wvel

		wvel_aux = matmul( shroud%COORDSYS, wvel_aux )

		call cPoints( shroud )
		call panelAreas( shroud )

		forall( i = 1:shroud%NP )
			shroud%CPUVW( :,i ) = shroud%VEL + cross( wvel_aux(:,1), shroud%CPXYZ( :,i ) - shroud%ORIGIN )
		end forall

		shroud%MVEL   = 0.d+0
		shroud%IVEL   = 0.d+0
		shroud%GAMMAS = 0.d+0
		shroud%NEWG   = 0.d+0
		shroud%OLDG   = 0.d+0
		shroud%LOAD   = 0.d+0
		
		allocate( shroud%wake )

		shroud%wake%NEP = nt
		shroud%wake%NEN = nt + 1
		shroud%wake%NVS = wake_size * (shroud%wake%NEN) &
				+ ( wake_size + 1) * shroud%wake%NEP

		allocate( shroud%wake%NODESXYZ ( 3,shroud%wake%NEN * (wake_size+1) ) )
		allocate( shroud%wake%NODESUVW ( 3,shroud%wake%NEN * wake_size     ) )
		allocate( shroud%wake%LOOPG    (   shroud%wake%NEP * wake_size,1   ) )
		allocate( shroud%wake%VORTEXG  ( shroud%wake%NVS,1 ) )
		allocate( shroud%wake%VORTEXP  ( 6,shroud%wake%NVS ) )
		allocate( shroud%wake%LOCMAT   ( 4,shroud%wake%NEP * wake_size     ) )
		allocate( shroud%wake%VORTEXPI ( 2,shroud%wake%NVS ) )
		allocate( shroud%wake%TEWAKEI  (   shroud%wake%NEP ) )
		allocate( shroud%wake%TELSNI   (   shroud%wake%NEN ) )
		allocate( shroud%wake%TELSSI   (   shroud%wake%NEP ) )

		shroud%wake%TELSNI = shroudtEdgeIndex ( nc, shroud%wake%NEN )
		shroud%wake%TELSSI = shroudtEdgePanels( shroud%wake%NEP )

		shroud%wake%TELSNI = shroudtEdgeIndex ( nc, shroud%wake%NEN )
		shroud%wake%TELSSI = shroudtEdgePanels( shroud%wake%NEP )
		shroud%wake%TEWAKEI = shroudWakePanelsMap( shroud%wake%NEP )

		call shroudWakeCorners( shroud )
		call wakeVortexP( shroud%wake, wake_size )

	end function initShroud 

  	function shroudNodes ( chord, radius, nc, nt, nn ) result(SNM)
		real(8), intent(in) :: chord, radius
  	  	real(8), dimension(:,:), allocatable :: SNM
		integer, intent(in) :: nc, nt, nn
  	  	integer :: node_count
		integer :: i, j
  	  	real(8) :: half_chord, theta, dc

  	  	theta = 2.d+0 * pi / nt
  	  	half_chord = chord * 0.5d+0   
		dc = chord / NC

  	  	allocate( SNM( 3,nn ) )
		
  	  	node_count = 1
  		do i=1,NT 
  		  do j=1,NC + 1
  		    SNM(1,node_count) = radius*cos(theta*(i-1))
  		    SNM(2,node_count) = radius*sin(theta*(i-1))
  		    SNM(3,node_count) = dc*(j-1)-half_chord
  		    node_count=node_count+1
  		  end do
  		end do

  	end function shroudNodes

	subroutine shroudWakeCorners( self )
		type( SurfaceType ), intent(inout) :: self

		allocate( self%wake%CORNERS( 1 ) )
		allocate( self%wake%CORMAP ( 1 ) )

		self%wake%CORNERS(1) = self%wake%NEN

		self%wake%CORMAP(1) = 1

	end subroutine

  	function shroudPanels( nc, nt, np ) result(SPM)
  	  	integer, dimension(:,:), allocatable :: SPM
		integer, intent(in) :: nc, nt, np
  	  	integer :: panel_count
		integer :: i, j

  		allocate(SPM(4,np))

  		panel_count=1

  		do j=1,NC
  			do i=1,NT
  				if (i == NT) then
  					SPM(1,panel_count) = j+(i-1) * (nc+1)
  					SPM(2,panel_count) = j
  					SPM(3,panel_count) = SPM(2,panel_count)+1
  					SPM(4,panel_count) = SPM(1,panel_count)+1
  					panel_count=panel_count + 1
  				else
  					SPM(1,panel_count) = j+(i-1)*(NC+1)
  					SPM(2,panel_count) = SPM(1,panel_count)+NC+1
  					SPM(3,panel_count) = SPM(2,panel_count)+1
  					SPM(4,panel_count) = SPM(1,panel_count)+1
  					panel_count=panel_count+1
  				end if
  			end do
  		end do

  	end function shroudPanels

	function shroudWakePanelsMap( nep ) result(wm)
    		!-------------------------------------------------
    		!   Maps the vortex segment on the trailing edge
    		!   of index 'i' with the wake panel of index 'wm(i)'
    		!
    		!   It is used by the subroutine 'gammas' of the
    		!   to obtain the circulations
		integer, intent(in) :: nep
		integer, dimension(:), allocatable :: wm
		integer :: i

		allocate( wm( nep ) )

		forall( i=1:NEP )
			wm( i ) = i
		endforall

    	end function shroudWakePanelsMap

	
	function shroudtEdgePanels( nep ) result(P)
        	!
        	! Maps the circulation of wake panel 'i' with trailing edge
        	! with bounded panel TE(i)
        	!
		integer, intent(in) :: nep
        	integer, dimension(:), allocatable :: P
		integer :: i

        	allocate( P( NEP ) )

		forall( i=1:NEP )
			P( i ) = i
		endforall

    	end function shroudtEdgePanels


    	function shroudtEdgeIndex( nc, nen ) result(TENI)
		integer, intent(in) :: nc, nen
        	integer, dimension(:), allocatable :: TENI
		integer :: i

        	allocate ( TENI ( nen ) )

		forall( i=1:nen - 1 )
			TENI( i ) = ( i-1 )*(NC+1) + 1
		endforall

		TENI( nen ) = TENI( 1 )

    	end function shroudtEdgeIndex


    	function snbmatrix( nc, nt, np ) result(NB)
        	integer, intent(in) :: nt, nc, np
        	integer, dimension(:,:), allocatable :: NB ! neighbour Matrix
		integer :: i, j, panel_counter

	        allocate(NB(4,NP))

		panel_counter = 1

		do i = 1,nc
			do j = 1,nt
				! First Row
				if( i == 1 )then
					NB(1,panel_counter)= - 1
					NB(2,panel_counter)= (i-1) * NT + mod( j,nt ) + 1 
					NB(3,panel_counter)= NT + j
					NB(4,panel_counter)= (i-1) * NT + mod( j-1,nt ) + nt*(1/j)
					panel_counter = panel_counter + 1
				elseif( i == nc )then
					NB(1,panel_counter)= (i-2)*nt + j
					NB(2,panel_counter)= (i-1) * NT + mod( j,nt ) + 1 
					NB(3,panel_counter)= -2
					NB(4,panel_counter)= (i-1) * NT + mod( j-1,nt ) + nt*(1/j)
					panel_counter = panel_counter + 1
				else
					NB(1,panel_counter)= (i-2)*nt + j
					NB(2,panel_counter)= (i-1) * NT + mod( j,nt ) + 1 
					NB(3,panel_counter)= i * nt + j
					NB(4,panel_counter)= (i-1) * NT + mod( j-1,nt ) + nt*(1/j)
					panel_counter = panel_counter + 1
				endif
			enddo
		enddo
	
		
    	end function snbmatrix



 end module shroud_module

