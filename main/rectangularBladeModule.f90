module rectangular_blade_module

    	use surface_module

    	implicit none
    	private

	public :: initRectBlade

	integer, private :: i, j, k

    	contains

	function initRectBlade( span, chord, tip_chord, core, ns, nc, &
				direction, alfa, twist, precone, azimuth, &
				radius, wake_size, euler_angs, orig, bvel, avel, slow_start ) result(blade)

		real(8), intent(in) :: span, chord, tip_chord, alfa, precone, azimuth, twist
		real(8), dimension(3), intent(in):: orig, bvel, avel, euler_angs
		real(8), intent(in) :: radius
		real(8), intent(inout) :: core
		real(8), dimension(3,1) :: wvel_aux
		integer, intent(in) :: ns, nc, wake_size, slow_start, direction
		!real(8) :: dt_aux, ds_aux
		type( SurfaceType )  :: blade

		blade%NN    = ( nc + 1 ) * ( ns + 1 )
		blade%NP    = nc * ns
		blade%NVS   = nc * (ns+1) + ns*(nc+1)

		allocate( blade%wake )

		blade%wake%NEN   = ns + 1
		blade%wake%NEP   = ns
		blade%wake%NVS   = (wake_size + 1) * blade%wake%NEP + wake_size * blade%wake%NEN

		blade%bluff = .false.
		blade%SLOWS = slow_start

		blade%COORDSYS = reshape([1.d+0, 0.d+0, 0.d+0,&
			    		0.d+0, 1.d+0, 0.d+0,&
			    		0.d+0, 0.d+0, 1.d+0],[3,3])

		call eulerRot( 3, 1, 2, euler_angs(1), euler_angs(2), euler_angs(3), blade%COORDSYS )

		blade%ORIGIN   = orig
		blade%VEL      = bvel
		blade%WVEL     = avel
		blade%EANGLES  = euler_angs

		allocate( blade%NODESXYZ ( 3,blade%NN  ) )
		allocate( blade%NODESLOC ( 3,blade%NN  ) )
		allocate( blade%CPXYZ    ( 3,blade%NP  ) )
		allocate( blade%CPUVW    ( 3,blade%NP  ) )
		allocate( blade%NORMAL   ( 3,blade%NP  ) )
		allocate( blade%MVEL     ( 3,blade%NP  ) )
		allocate( blade%IVEL     ( 3,blade%NP  ) )
		allocate( blade%GAMMAS   ( 4,blade%NP  ) )
		allocate( blade%VORTEX   ( 7,blade%NVS ) )
		allocate( blade%NEWG     ( blade%NP,1  ) )
		allocate( blade%OLDG     ( blade%NP,1  ) )
		allocate( blade%LOAD     ( blade%NP,1  ) )
		allocate( blade%AREA     ( blade%NP,1  ) )
		allocate( blade%LOCMAT   ( 4,blade%NP  ) )
		allocate( blade%NBMAT    ( 4,blade%NP  ) )
		
		allocate( blade%wake%NODESXYZ ( 3,blade%wake%NEN * (wake_size+1) ) )
		allocate( blade%wake%NODESUVW ( 3,blade%wake%NEN * wake_size     ) )
		allocate( blade%wake%LOOPG    (   blade%wake%NEP * wake_size,1   ) )
		allocate( blade%wake%VORTEXG  ( blade%wake%NVS,1 ) )
		allocate( blade%wake%VORTEXP  ( 6,blade%wake%NVS ) )
		allocate( blade%wake%LOCMAT   ( 4,blade%wake%NEP * wake_size     ) )
		allocate( blade%wake%VORTEXPI ( 2,blade%wake%NVS ) )
		allocate( blade%wake%TEWAKEI  (   blade%wake%NEP ) )
		allocate( blade%wake%TELSNI   (   blade%wake%NEN ) )
		allocate( blade%wake%TELSSI   (   blade%wake%NEP ) )

		blade%NODESLOC = bladeNodes( span, chord, tip_chord, ns, nc, twist )
		blade%LOCMAT   = bladePanels( blade%NP, ns, nc ) 
		blade%NBMAT    = bnbmatrix( ns, nc )

		forall( i = 1:blade%NN )
			blade%NODESLOC( 2,i ) = blade%NODESLOC( 2,i ) + radius
		end forall

		call eulerRot( 3, 1, 2, azimuth, precone, alfa, blade%NODESLOC )
		blade%NODESLOC(2,:) = direction * blade%NODESLOC(2,:)

		blade%NODESXYZ = matmul( blade%COORDSYS, blade%NODESLOC ) 

		forall( i = 1:blade%NN )
			blade%NODESXYZ( :,i ) = blade%NODESXYZ( :,i ) + blade%ORIGIN
		end forall

		call panelAreas( blade )
		call cPoints( blade )

		if( slow_start > 0 )then
			wvel_aux(:,1) = avel * (1./slow_start)
		else
			wvel_aux(:,1) = avel
		endif

		wvel_aux = matmul( blade%COORDSYS, wvel_aux )

		do i = 1,blade%NP
			blade%CPUVW( :,i ) = blade%VEL + cross( wvel_aux(:,1), blade%CPXYZ( :,i ) - blade%ORIGIN )
		end do

		blade%MVEL   = 0.d+0
		blade%IVEL   = 0.d+0
		blade%GAMMAS = 0.d+0
		blade%NEWG   = 0.d+0
		blade%OLDG   = 0.d+0
		blade%LOAD   = 0.d+0
		
		blade%wake%TELSNI = bladetEdgeIndex( blade%wake%NEN, ns, nc)
		blade%wake%TELSSI = bladetEdgePanels( blade%wake%NEP, ns, nc )
   		call bladeWakeCorners( blade )
		call wakeVortexI( blade%wake, wake_size )

		blade%wake%NODESXYZ( :,1:blade%wake%NEN ) = blade%NODESXYZ( :,blade%wake%TELSNI )

		blade%wake%NODESUVW = 0.d+0
		blade%wake%LOOPG    = 0.d+0
		blade%wake%VORTEXG  = 0.d+0

		call wakePanels( blade%wake, wake_size )
		blade%wake%TEWAKEI = bladeWakePanelsMap( blade%wake%NEP )
		blade%wake%ROOT = .true.
		!blade%wake%ROOT = .false.

	end function initRectBlade


	subroutine bladeWakeCorners( self )
		type( SurfaceType ), intent(inout) :: self

		allocate( self%wake%CORNERS( 1 ) )
		allocate( self%wake%CORMAP ( 1 ) )

		self%wake%CORNERS(1) = 1!self%wake%NEN

		self%wake%CORMAP(1) = 1!self%wake%NEN

	end subroutine


	function bladeWakePanelsMap( nep ) result(wm)
	!function bladeWakePanelsMap( nep, ns, nc ) result(wm)
    		!-------------------------------------------------
    		!   Maps the vortex segment on the trailing edge
    		!   of index 'i' with the wake panel of index 'wm(i)'
    		!
    		!   It is used by the subroutine 'gammas' of the
    		!   to obtain the circulations
    		!   of the wake when calculating the circulation on
    		!   the trailing edge segments
    		!-------------------------------------------------
        	integer, dimension(:), allocatable :: wm
        	integer, intent(in) :: nep!, ns, nc

        	allocate ( wm ( nep ) )

        	!forall (i=1:nc + ns - 1)
        	!    wm(i) = i
        	!end forall

        	!forall( i = 1:nc + 1 )
        	!    wm( nc + ns - 1 + i ) = nep + 1 - i
        	!end forall


        	forall (i=1:nep)
        	    wm(i) = i
        	end forall

    	end function bladeWakePanelsMap

	function bladetEdgePanels( np, ns, nc ) result(P)
        	!
        	! Maps the circulation of wake panel 'i' with trailing edge
        	! with bounded panel TE(i)
        	!
        	integer, dimension(:), allocatable :: P
        	integer, intent(in) :: nc, ns, np
		!integer :: ind

        	allocate (P(np))

        !	do i = 1, nc
        !	    P ( i ) = i
        !	end do

        !	P ( nc + 1 ) = P ( nc )

        !	do i = 2,ns
        !	    ind = i + nc
        !	    P ( ind ) = P ( ind - 1 ) + nc
        !	end do

        !	P ( nc + ns + 1 ) = P ( nc + ns )

        !	do i = 2, nc
        !	    ind = nc + ns + i
        !	    P ( ind ) = P ( ind - 1 ) - 1
        !	end do


        	do i = 1,ns
        	    P ( i ) = i * nc
        	end do

    	end function bladetEdgePanels


    	!function bladetEdgeIndex( nn, nse, nce ) result(TENI)
        !	integer, dimension(:), allocatable :: TENI
        !	integer, intent(in) :: nn, nce, nse
	!	integer :: ind, nc, ns

        !	nc = nce + 1
        !	ns = nse + 1

        !	allocate ( TENI ( nn ) )

        !	do i = 1,nc
        !	    TENI ( i ) = ( i - 1 ) * ns + 1
        !	end do

        !	do i = 1,ns - 1
        !	    ind = nc + i
        !	    TENI ( ind ) = TENI ( ind - 1 ) + 1
        !	end do

        !	do i=1,nc - 1
        !	    ind = nc + ns - 1 + i
        !	    TENI ( ind ) = TENI ( ind - 1 ) - ns
        !	end do

    	!end function bladetEdgeIndex


    	function bladetEdgeIndex( nen, ns, nc) result(TENI)
        	integer, dimension(:), allocatable :: TENI
        	integer, intent(in) :: nen, nc, ns
		integer :: nsi, nci

        	nci = nc + 1
        	nsi = ns + 1

        	allocate ( TENI ( nen ) )

        	do i = 1,nsi
        	    TENI ( i ) = nsi * (nci - 1) + i
        	end do

    	end function bladetEdgeIndex

    	!function bladeNodes( span, chord, ns, nc ) result(NM)
        !	real(8), dimension(:,:), allocatable :: NM ! Node Matrix
	!	real(8), intent(in) :: span, chord
	!	integer, intent(in) :: ns, nc
        !	real(8) :: center
        !	real(8) :: ds, dc         ! DeltaS and DeltaC
        !	integer :: node_count, NN

	!	NN = (ns + 1)*(nc + 1)

        !	ds = span / ns
        !	dc = chord / nc

        !	node_count = 1
        !	center = chord * 0.25+0

        !	allocate( NM( 3,NN ) )

        !	do i = 1,nc + 1
        !	    do j = 1,ns + 1
        !	        NM(1,node_count) = center - 0.5d+0 * chord + ((i-1)*dc)
        !	        NM(2,node_count) = ((j-1)*ds)
        !	        NM(3,node_count) = 0.0d+0
        !	        node_count = node_count + 1
        !	    end do
        !	end do

    	!end function bladeNodes

    	function bladeNodes( span, chord, tip_chord, ns, nc, twist_deg ) result(NM)
        	real(8), dimension(:,:), allocatable :: NM ! Node Matrix
		real(8), intent(in) :: span, chord, twist_deg, tip_chord
		integer, intent(in) :: ns, nc
        	real(8) :: center_chord, half_span
		real(8), dimension(:), allocatable :: span_loc
		real(8), dimension(:), allocatable :: twist_factor
		real(8), dimension(:), allocatable :: chord_factor
		real(8), dimension(3) :: mean_chord
        	real(8) :: ds, dc, twist         ! DeltaS and DeltaC
        	integer :: node_count, NN

		NN = (ns + 1)*(nc + 1)

        	ds = (1.0d+0 * pi) / (ns+2)
        	!ds = (0.5d+0 * pi) / (ns)
        	dc = chord / nc
		twist = deg2rad(twist_deg)

        	node_count = 1
        	center_chord = chord * 0.0d+0
		half_span = 0.5d+0 * span

        	allocate( NM( 3,NN ) )
		allocate( span_loc( ns + 1 ) )
		allocate( twist_factor( ns + 1 ) )
		allocate( chord_factor( ns + 1 ) )

		span_loc(1) = 0.d+0
		do i=2,ns
			span_loc(i) = half_span - half_span * cos( i*ds )
			!span_loc(i) = span * sin( (i-1)*ds )
		enddo
		span_loc(ns+1) = span

		twist_factor(1) = 0.d+0
		do i=2,ns+1
			twist_factor(i) = (span_loc(i) / span)
		enddo

		do i=1,ns+1
			chord_factor(i)= (1.d+0 - ( (chord - tip_chord) / chord) * ( span_loc(i) / span ))
		enddo

        	do i = 1,nc + 1
        	    do j = 1,ns + 1
		    	mean_chord = [cos(twist_factor(j)*twist)*chord_factor(j)*(center_chord - 0.5d+0 * chord + (i-1)*dc),&
				0.d+0,&
				sin(twist_factor(j)*twist)*chord_factor(j)*(center_chord - 0.5d+0 * chord + (i-1)*dc)]
        	        NM(1,node_count) = mean_chord(1)
        	        NM(2,node_count) = span_loc(j)
        	        NM(3,node_count) = mean_chord(3)
        	        node_count = node_count + 1
        	    end do
        	end do

    	end function bladeNodes

	
	function bladePanels( np, ns, nc ) result(PM)
	!
	! This function creates the node connectivity or location matrix of the panels.
	! The panel nodes are defined in a counter-clockwise direction.
	!
        	integer, dimension(:,:), allocatable :: PM ! Blade Panel Matrix
        	integer :: panel_count
		integer, intent(in) :: np, ns, nc

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

    	end function bladePanels

    	function bnbmatrix( ns, nc ) result(NB)
        	integer, dimension(:,:), allocatable :: NB ! neighbour Matrix
        	integer, intent(in) :: ns, nc

	        allocate(NB( 4,nc*ns ))
	
		!!
		!! Corners
		!!
	        !! UL
	        !NB(1,1) = 2
	        !NB(2,1) = nc + 1
	        !NB(3,1) = -2
	        !NB(4,1) = -1
	        !! UR
	        !NB(1, (ns-1)*nc + 1 ) = nc * (ns - 1) + 2
	        !NB(2, (ns-1)*nc + 1 ) = -1
	        !NB(3, (ns-1)*nc + 1 ) = -2
	        !NB(4, (ns-1)*nc + 1 ) = nc * (ns - 2) + 1
	        !! LL
	        !NB(1, nc) = -1
	        !NB(2, nc) = 2*nc
	        !NB(3, nc) = nc - 1
	        !NB(4, nc) = -1
	        !! LR
	        !NB(1,ns*nc) = -1
	        !NB(2,ns*nc) = -1
	        !NB(3,ns*nc) = nc*ns - 1
	        !NB(4,ns*nc) = (ns - 1) * nc
		!!
		!! Interior
		!!
	        !forall ( i = 2:(ns-1), j = 2:(nc-1) )
	        !    NB(1,(i-1)*nc + j) = (i-1) * nc + j + 1
	        !    NB(2,(i-1)*nc + j) = (i)   * nc + j
	        !    NB(3,(i-1)*nc + j) = (i-1) * nc + j - 1
	        !    NB(4,(i-1)*nc + j) = (i-2) * nc + j
	        !end forall
		!!
		!! Upper Row
		!!
	        !forall ( k = nc+1 : (ns-2)*nc + 1 : nc )
	        !    NB(1,k) = k + 1
	        !    NB(2,k) = k + nc
	        !    NB(3,k) = -2
	        !    NB(4,k) = k - nc
	        !end forall
		!!
		!! Left Row
		!!
	        !forall ( k = 2:nc-1 )
	        !    NB(1,k) = k + 1
	        !    NB(2,k) = k + nc
	        !    NB(3,k) = k - 1
	        !    NB(4,k) = -1
	        !end forall
		!!
		!! Right row
		!!
	        !forall ( k = (ns-1)*nc + 2 : nc*ns - 1 )
	        !    NB(1,k) = k + 1
	        !    NB(2,k) = -1
	        !    NB(3,k) = k - 1
	        !    NB(4,k) = k - nc
	        !end forall
		!!
		!! Lower row
		!!
	        !forall ( k = 2*nc : (ns-1)*nc : nc )
	        !    NB(1,k) = -1
	        !    NB(2,k) = k + nc
	        !    NB(3,k) = k - 1
	        !    NB(4,k) = k - nc
	        !end forall
		
		!
		! Corners
		!
	        ! UL
	        NB(1,1) = 2
	        NB(2,1) = nc + 1
	        NB(3,1) = -2
	        NB(4,1) = -2
	        ! UR
	        NB(1, (ns-1)*nc + 1 ) = nc * (ns - 1) + 2
	        NB(2, (ns-1)*nc + 1 ) = -2
	        NB(3, (ns-1)*nc + 1 ) = -2
	        NB(4, (ns-1)*nc + 1 ) = nc * (ns - 2) + 1
	        ! LL
	        NB(1, nc) = -1
	        NB(2, nc) = 2*nc
	        NB(3, nc) = nc - 1
	        NB(4, nc) = -2
	        ! LR
	        NB(1,ns*nc) = -1
	        NB(2,ns*nc) = -2
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
	            NB(4,k) = -2
	        end forall
		!
		! Right row
		!
	        forall ( k = (ns-1)*nc + 2 : nc*ns - 1 )
	            NB(1,k) = k + 1
	            NB(2,k) = -2
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
		
    	end function bnbmatrix
	
end module rectangular_blade_module
