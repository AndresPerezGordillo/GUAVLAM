module rectangular_wing_module

    	use surface_module
    	use wake_module
	use math_routines

    	implicit none
    	private

	integer, private :: i,j,k

	public initRectWing



    contains

	function initRectWing( span, chord, ns, nc, alfa, wake_size, euler_angs, orig, bvel, avel ) result( wing )
		real(8), intent(in) :: span, chord, alfa
		real(8), dimension(3), intent(in):: orig, bvel, avel, euler_angs
		integer, intent(in) :: ns, nc, wake_size
		type(SurfaceType)  :: wing

		wing%NN    = ( nc + 1 ) * ( ns + 1 )
		wing%NP    = nc * ns
		wing%NVS   = nc*(ns+1) + ns*(nc+1)
		wing%bluff = .false.

		wing%EANGLES = euler_angs

		allocate( wing%wake )
		
		wing%wake%NEN   = nc * 2 + ns + 1
		wing%wake%NEP   = nc * 2 + ns
		wing%wake%NVS   = (wake_size+1) * (wing%wake%NEP) + wake_size*(wing%wake%NEN)

		wing%COORDSYS = reshape([1.0d+0, 0.d+0, 0.d+0,&
					0.0d+0, 1.d+0, 0.d+0,&
					0.0d+0, 0.d+0, 1.d+0],[3,3])

		call euler_rot( 3, 1, 2, euler_angs(1), euler_angs(1), euler_angs(1), wing%COORDSYS )

		wing%ORIGIN   = orig
		wing%VEL      = bvel
		wing%WVEL     = avel

		allocate( wing%NODESXYZ ( 3,wing%NN  ) )
		allocate( wing%NODESLOC ( 3,wing%NN  ) )
		allocate( wing%CPXYZ    ( 3,wing%NP  ) )
		allocate( wing%CPUVW    ( 3,wing%NP  ) )
		allocate( wing%NORMAL   ( 3,wing%NP  ) )
		allocate( wing%MVEL     ( 3,wing%NP  ) )
		allocate( wing%IVEL     ( 3,wing%NP  ) )
		allocate( wing%GAMMAS   ( 4,wing%NP  ) )
		allocate( wing%VORTEX   ( 7,wing%NVS ) )
		allocate( wing%NEWG     ( wing%NP,1  ) )
		allocate( wing%OLDG     ( wing%NP,1  ) )
		allocate( wing%LOAD     ( wing%NP,1  ) )
		allocate( wing%AREA     ( wing%NP,1  ) )
		allocate( wing%LOCMAT   ( 4,wing%NP  ) )
		allocate( wing%NBMAT    ( 4,wing%NP  ) )
		
		allocate( wing%wake%NODESXYZ ( 3,wing%wake%NEN * (wake_size+1) ) )
		allocate( wing%wake%NODESUVW ( 3,wing%wake%NEN * wake_size     ) )
		allocate( wing%wake%LOOPG    (   wing%wake%NEP * wake_size,1   ) )
		allocate( wing%wake%VORTEXG  ( wing%wake%NVS,1 ) )
		allocate( wing%wake%VORTEXP  ( 6,wing%wake%NVS ) )
		allocate( wing%wake%LOCMAT   ( 4,wing%wake%NEP * wake_size     ) )
		allocate( wing%wake%VORTEXPI ( 2,wing%wake%NVS ) )
		allocate( wing%wake%TEWAKEI  (   wing%wake%NEP ) )
		allocate( wing%wake%TELSNI   (   wing%wake%NEN ) )
		allocate( wing%wake%TELSSI   (   wing%wake%NEP ) )

		call wakeCorners( wing, ns, nc )
		call wakeVortexP( wing%wake, wake_size )

		wing%NODESLOC = nodes ( span, chord, ns, nc, wing%nn )
		wing%LOCMAT   = panels( ns, nc, wing%np )
		wing%NBMAT    = rwnbmatrix( ns, nc, wing%np )

		call euler_rot( 2, alfa, wing%NODESLOC )

		wing%NODESXYZ = matmul( wing%COORDSYS, wing%NODESLOC ) 

		forall( i = 1:wing%NN )
			wing%NODESXYZ( :,i ) = wing%NODESXYZ( :,i ) + wing%ORIGIN
		end forall

		call cPoints( wing )
		call panelAreas( wing )

		forall( i = 1:wing%NP )
			wing%CPUVW( :,i ) = wing%VEL + cross( wing%WVEL, wing%CPXYZ( :,i ) - wing%ORIGIN )
		end forall

		wing%MVEL   = 0.d+0
		wing%GAMMAS = 0.d+0
		wing%NEWG   = 0.d+0
		wing%OLDG   = 0.d+0
		wing%LOAD   = 0.d+0
		
		wing%wake%TELSNI = tEdgeIndex ( ns, nc, wing%wake%NEN )
		wing%wake%TELSSI = tEdgePanels( ns, nc, wing%wake%NEP )

		wing%wake%NODESXYZ = 0.d+0
		wing%wake%NODESXYZ( :,1:wing%wake%NEN ) = wing%NODESXYZ( :,wing%wake%TELSNI )
		!call wakeVortexPos( wing%wake )

		wing%wake%NODESUVW = 0.d+0
		wing%wake%LOOPG    = 0.d+0
		wing%wake%VORTEXG  = 0.d+0

		call wakePanels( wing%wake, wake_size )
		wing%wake%TEWAKEI = wakePanelsMap( ns, nc, wing%wake%NEP )

	end function initRectWing


	subroutine wakeCorners( self, ns, nc )
		type( SurfaceType ), intent(inout) :: self
		integer, intent(in) :: ns, nc
		
		allocate( self%wake%CORNERS( 2 ) )
		allocate( self%wake%CORMAP ( 2 ) )

		self%wake%CORNERS(1) = self%wake%NEN + NC + 1
		self%wake%CORNERS(2) = self%wake%NEN + NC + 1 + NS

		self%wake%CORMAP(1) = self%wake%NEN + NC + 1 - 1
		self%wake%CORMAP(2) = self%wake%NEN + NC + 1 + NS + 1
		
	end subroutine


	function wakePanelsMap( ns, nc, n ) result(wm)
    		!-------------------------------------------------
    		!   Maps the vortex segment on the trailing edge
    		!   of index 'i' with the wake panel of index 'wm(i)'
    		!
    		!   It is used by the subroutine 'gammas' of the
    		!   'vortex_sheets' module to obtain the circulations
    		!   of the wake when calculating the circulation on
    		!   the trailing edge segments
    		!-------------------------------------------------
        	integer, dimension(:), allocatable :: wm
        	integer, intent(in) :: n, ns, nc

        	allocate ( wm ( n ) )

        	forall (i=1:nc + ns - 1)
        	    wm(i) = i
        	end forall

        	forall( i = 1:nc + 1 )
        	    wm( nc + ns - 1 + i ) = n + 1 - i
        	end forall

    	end function wakePanelsMap


	
	function tEdgePanels( ns, nc, np ) result(P)
        	!
        	! Maps the circulation of wake panel 'i' with trailing edge
        	! with bounded panel TE(i)
        	!
        	integer, dimension(:), allocatable :: P
        	integer, intent(in) :: nc, ns, np
		integer :: ind

        	allocate (P(np))

        	do i = 1, nc
        	    P ( i ) = i
        	end do

        	P ( nc + 1 ) = P ( nc )

        	do i = 2,ns
        	    ind = i + nc
        	    P ( ind ) = P ( ind - 1 ) + nc
        	end do

        	P ( nc + ns + 1 ) = P ( nc + ns )

        	do i = 2, nc
        	    ind = nc + ns + i
        	    P ( ind ) = P ( ind - 1 ) - 1
        	end do

    	end function tEdgePanels


    	function tEdgeIndex( ns, nc, nn ) result(TENI)
        	integer, dimension(:), allocatable :: TENI
        	integer, intent(in) :: nn, nc, ns
		integer :: ind, nci, nsi

		nsi = ns + 1
		nci = nc + 1

        	allocate ( TENI ( nn ) )

        	do i = 1,nci
        	    TENI ( i ) = ( i - 1 ) * nsi + 1
        	end do

        	do i = 1,nsi - 1
        	    ind = nci + i
        	    TENI ( ind ) = TENI ( ind - 1 ) + 1
        	end do

        	do i=1,nci - 1
        	    ind = nci + nsi - 1 + i
        	    TENI ( ind ) = TENI ( ind - 1 ) - nsi
        	end do

    	end function tEdgeIndex


    	function nodes( span, chord, ns, nc, nn ) result(NM)
		real(8), intent(in) :: span, chord
        	real(8), dimension(:,:), allocatable :: NM ! Node Matrix
        	real(8) :: center
        	real(8) :: ds, dc         ! DeltaS and DeltaC
        	integer :: node_count
		integer, intent(in) :: ns, nc, nn

        	ds = span / ns
        	dc = chord / nc

        	node_count = 1
        	center = span * 0.5d+0

        	allocate( NM( 3,NN ) )

        	do i = 1,nc + 1
        	    do j = 1,ns + 1
        	        NM(1,node_count) = (i-1)*dc
        	        NM(2,node_count) = (j-1)*ds - center
        	        NM(3,node_count) = 0.0d+0
        	        node_count = node_count + 1
        	    end do
        	end do

    	end function nodes

	
	function panels( ns, nc, np ) result(PM)
	!
	! This function creates the node connectivity or location matrix of the panels.
	! The panel nodes are defined in a counter-clockwise direction.
	!
        	integer, dimension(:,:), allocatable :: PM ! Blade Panel Matrix
        	integer :: panel_count
		integer, intent(in) :: ns, nc, np

        	allocate (PM(4,NP))
        	panel_count = 1

 		! Do construct for the location matrix of the panels
        	do i=1,ns
        	    do j=1,nc
        	        PM(1,panel_count)=i+j*(ns + 1)
        	        PM(2,panel_count)=PM(1,panel_count)+1
        	        PM(3,panel_count)=PM(2,panel_count)-ns-1
        	        PM(4,panel_count)=PM(3,panel_count)-1
        	        panel_count=panel_count+1
        	    end do
        	end do

    	end function panels


    	function rwnbmatrix( ns, nc, np ) result(NB)
        	integer, dimension(:,:), allocatable :: NB ! neighbour Matrix
        	integer, intent(in) :: ns, nc, np

	        allocate(NB(4,np))
	
		!
		! Corners
		!
	        ! UL
	        NB(1,1) = 2
	        NB(2,1) = nc + 1
	        NB(3,1) = -2
	        NB(4,1) = -1
	        ! UR
	        NB(1, (ns-1)*nc + 1 ) = nc * (ns - 1) + 2
	        NB(2, (ns-1)*nc + 1 ) = -1
	        NB(3, (ns-1)*nc + 1 ) = -2
	        NB(4, (ns-1)*nc + 1 ) = nc * (ns - 2) + 1
	        ! LL
	        NB(1, nc) = -1
	        NB(2, nc) = 2*nc
	        NB(3, nc) = nc - 1
	        NB(4, nc) = -1
	        ! LR
	        NB(1,ns*nc) = -1
	        NB(2,ns*nc) = -1
	        NB(3,ns*nc) = nc*ns - 1
	        NB(4,ns*nc) = (ns - 1) * nc
		!
		! Interior
		!
	        forall ( i = 2:(ns-1), j = 2:(nc-1) )
	            NB(1,(i-1)*nc + j) = (i-1) * nc + j + 1
	            NB(2,(i-1)*nc + j) = (i)   * nc + j
	            NB(3,(i-1)*nc + j) = (i-1) * nc + j - 1
	            NB(4,(i-1)*nc + j) = (i-2) * nc + j
	        end forall
		!
		! Upper Row
		!
	        forall ( k = nc+1 : (ns-2)*nc + 1 : nc )
	            NB(1,k) = k + 1
	            NB(2,k) = k + nc
	            NB(3,k) = -2
	            NB(4,k) = k - nc
	        end forall
		!
		! Left Row
		!
	        forall ( k = 2:nc-1 )
	            NB(1,k) = k + 1
	            NB(2,k) = k + nc
	            NB(3,k) = k - 1
	            NB(4,k) = -1
	        end forall
		!
		! Right row
		!
	        forall ( k = (ns-1)*nc + 2 : nc*ns - 1 )
	            NB(1,k) = k + 1
	            NB(2,k) = -1
	            NB(3,k) = k - 1
	            NB(4,k) = k - nc
	        end forall
		!
		! Lower row
		!
	        forall ( k = 2*nc : (ns-1)*nc : nc )
	            NB(1,k) = -1
	            NB(2,k) = k + nc
	            NB(3,k) = k - 1
	            NB(4,k) = k - nc
	        end forall
		
    	end function rwnbmatrix
	
end module rectangular_wing_module
