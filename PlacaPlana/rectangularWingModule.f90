module rectangular_wing_module

    	use surface_module
    	use wake_module
	use math_routines

    	implicit none
    	public

	integer, private :: i,j,k

	private panels, rwnbmatrix

    	type, extends( LiftingSurfaceType ) :: RectWingType

    	    real(8) :: span, chord
    	    integer :: ns, nc

	    contains

		procedure, pass( wing ) :: DESTROY

    	end type RectWingType

    contains

	function initRectWing( span, chord, ns, nc, wake_size, c_sys, orig, bvel, avel ) result( wing )
		real(8), intent(in) :: span, chord
		real(8), dimension(3,3), intent(in) :: c_sys
		real(8), dimension(3), intent(in):: orig, bvel, avel
		integer, intent(in) :: ns, nc, wake_size
		type(RectWingType)  :: wing

		wing%span  = span
		wing%chord = chord
		wing%ns    = ns
		wing%nc    = nc

		wing%NN    = ( nc + 1 ) * ( ns + 1 )
		wing%NP    = nc * ns
		wing%NEN   = nc * 2 + ns + 1
		wing%NEP   = nc * 2 + ns

		wing%COORDSYS = c_sys
		wing%ORIGIN   = orig
		wing%VEL      = bvel
		wing%WVEL     = avel

		allocate( wing%NODESXYZ ( 3,wing%NN  ) )
		allocate( wing%NODESLOC ( 3,wing%NN  ) )
		allocate( wing%CPXYZ    ( 3,wing%NP  ) )
		allocate( wing%CPUVW    ( 3,wing%NP  ) )
		allocate( wing%NORMAL   ( 3,wing%NP  ) )
		allocate( wing%MVEL     ( 3,wing%NP  ) )
		allocate( wing%GAMMAS   ( 4,wing%NP  ) )
		allocate( wing%NEWG     ( wing%NP,1  ) )
		allocate( wing%OLDG     ( wing%NP,1  ) )
		allocate( wing%LOAD     ( wing%NP,1  ) )
		allocate( wing%LOCMAT   ( 4,wing%NP  ) )
		allocate( wing%NBMAT    ( 4,wing%NP  ) )
		allocate( wing%TELSNI   (   wing%NEN ) )
		allocate( wing%TELSSI   (   wing%NEP ) )
		
		allocate( wing%wake%NODESXYZ ( 3,wing%NEN * (wake_size+1) ) )
		allocate( wing%wake%NODESUVW ( 3,wing%NEN * wake_size     ) )
		allocate( wing%wake%GAMMAS   ( 4,wing%NEP * wake_size     ) )
		allocate( wing%wake%LOOPG    (   wing%NEP * wake_size,1   ) )
		allocate( wing%wake%LOCMAT   ( 4,wing%NEP * wake_size     ) )
		allocate( wing%wake%NBMAT    ( 4,wing%NEP * wake_size     ) )
		allocate( wing%wake%TEWAKEI  (   wing%NEP                 ) )

		call wakeCorners( wing )

		wing%NODESLOC = nodes ( wing )
		wing%LOCMAT   = panels( wing )
		wing%NBMAT  = rwnbmatrix ( wing )

		wing%NODESXYZ = matmul( wing%COORDSYS, wing%NODESLOC ) 

		forall( i = 1:wing%NN )
			wing%NODESXYZ( :,i ) = wing%NODESXYZ( :,i ) + wing%ORIGIN
		end forall

		call cPoints( wing )

		forall( i = 1:wing%NP )
			wing%CPUVW( :,i ) = wing%VEL + cross( wing%WVEL, wing%CPXYZ( :,i ) - wing%ORIGIN )
		end forall

		wing%MVEL   = 0.d+0
		wing%GAMMAS = 0.d+0
		wing%NEWG   = 0.d+0
		wing%OLDG   = 0.d+0
		wing%LOAD   = 0.d+0
		
		wing%TELSNI = tEdgeIndex ( wing )
		wing%TELSSI = tEdgePanels( wing )

		wing%wake%NODESXYZ( :,1:wing%NEN ) = wing%NODESXYZ( :,wing%TELSNI )

		wing%wake%NODESUVW = 0.d+0
		wing%wake%GAMMAS   = 0.d+0
		wing%wake%LOOPG    = 0.d+0

		wing%wake%LOCMAT  = wakePanels( wing%NEP, wake_size )
		wing%wake%NBMAT   = wakeNB( wing%NEP, wake_size )
		wing%wake%TEWAKEI = wakePanelsMap( wing )

	end function initRectWing

	subroutine destroy( wing )
		class( RectWingType ), intent(inout) :: wing

		deallocate( wing%NODESXYZ )
		deallocate( wing%NODESLOC )
		deallocate( wing%CPXYZ    )
		deallocate( wing%CPUVW    )
		deallocate( wing%NORMAL   )
		deallocate( wing%MVEL     )
		deallocate( wing%GAMMAS   )
		deallocate( wing%NEWG     )
		deallocate( wing%OLDG     )
		deallocate( wing%LOAD     )
		deallocate( wing%LOCMAT   )
		deallocate( wing%NBMAT    )
		deallocate( wing%TELSNI   )
		deallocate( wing%TELSSI   )
		deallocate( wing%wake%NODESXYZ )
		deallocate( wing%wake%NODESUVW )
		deallocate( wing%wake%GAMMAS   )
		deallocate( wing%wake%LOOPG    )
		deallocate( wing%wake%LOCMAT   )
		deallocate( wing%wake%NBMAT    )
		deallocate( wing%wake%TEWAKEI  )
		deallocate( wing%wake%CORNERS  )
		deallocate( wing%wake%CORMAP   )

	end subroutine destroy

	subroutine wakeCorners(self)
		class( RectWingType ), intent(inout) :: self
		
		allocate( self%wake%CORNERS( 2 ) )
		allocate( self%wake%CORMAP ( 2 ) )

		self%wake%CORNERS(1) = self%NEN + self%NC + 1
		self%wake%CORNERS(2) = self%NEN + self%NC + 1 + self%NS

		self%wake%CORMAP(1) = self%NEN + self%NC + 1 - 1
		self%wake%CORMAP(2) = self%NEN + self%NC + 1 + self%NS + 1
		
	end subroutine


	function wakePanelsMap(self) result(wm)
    		!-------------------------------------------------
    		!   Maps the vortex segment on the trailing edge
    		!   of index 'i' with the wake panel of index 'wm(i)'
    		!
    		!   It is used by the subroutine 'gammas' of the
    		!   'vortex_sheets' module to obtain the circulations
    		!   of the wake when calculating the circulation on
    		!   the trailing edge segments
    		!-------------------------------------------------
        	class(RectWingType), intent(in) :: self
        	integer, dimension(:), allocatable :: wm
        	integer :: n, ns, nc

        	n = self%NEP
        	ns = self%ns
        	nc = self%nc

        	allocate ( wm ( n ) )

        	forall (i=1:nc + ns - 1)
        	    wm(i) = i
        	end forall

        	forall( i = 1:nc + 1 )
        	    wm( nc + ns - 1 + i ) = n + 1 - i
        	end forall

    	end function wakePanelsMap


	
	function tEdgePanels(self) result(P)
        	!
        	! Maps the circulation of wake panel 'i' with trailing edge
        	! with bounded panel TE(i)
        	!
        	class(RectWingType), intent(in) :: self
        	integer, dimension(:), allocatable :: P
        	integer :: nc, ns, np, ind

        	np = self%NEP
        	nc = self%nc
        	ns = self%ns

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


    	function tEdgeIndex(self) result(TENI)
	
        	class(RectWingType), intent(in) :: self
        	integer, dimension(:), allocatable :: TENI
        	integer :: nn, nc, ns, ind

        	nn = self%NEN
        	nc = self%nc + 1
        	ns = self%ns + 1

        	allocate ( TENI ( nn ) )

        	do i = 1,nc
        	    TENI ( i ) = ( i - 1 ) * ns + 1
        	end do

        	do i = 1,ns - 1
        	    ind = nc + i
        	    TENI ( ind ) = TENI ( ind - 1 ) + 1
        	end do

        	do i=1,nc - 1
        	    ind = nc + ns - 1 + i
        	    TENI ( ind ) = TENI ( ind - 1 ) - ns
        	end do

    	end function tEdgeIndex


    	function nodes(self) result(NM)
		type(RectWingType), intent(in) :: self
        	real(8), dimension(:,:), allocatable :: NM ! Node Matrix
        	real(8) :: center
        	real(8) :: ds, dc         ! DeltaS and DeltaC
        	integer :: node_count, ns, nc

        	ds = self%span / self%ns
        	dc = self%chord / self%nc

        	ns = self%ns
        	nc = self%nc

        	node_count = 1
        	center = self%span * 0.5d+0

        	allocate( NM( 3,self%NN ) )

        	do i = 1,nc + 1
        	    do j = 1,ns + 1
        	        NM(1,node_count) = (i-1)*dc
        	        NM(2,node_count) = (j-1)*ds - center
        	        NM(3,node_count) = 0.0d+0
        	        node_count = node_count + 1
        	    end do
        	end do

    	end function nodes

	
	function panels(self) result(PM)
	!
	! This function creates the node connectivity or location matrix of the panels.
	! The panel nodes are defined in a counter-clockwise direction.
	!
        	class(RectWingType), intent(in) :: self
        	integer, dimension(:,:), allocatable :: PM ! Blade Panel Matrix
        	integer :: panel_count, ns, nc

        	ns = self%ns
        	nc = self%nc

        	allocate (PM(4,self%NP))
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


    	function rwnbmatrix(self) result(NB)
        	class(RectWingType), intent(in) :: self
        	integer, dimension(:,:), allocatable :: NB ! neighbour Matrix
        	integer :: ns, nc

	        ns = self%ns
	        nc = self%nc
	
	        allocate(NB(4,nc*ns))
	
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
