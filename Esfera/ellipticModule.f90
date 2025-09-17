!-----------------------------------
! Project: pruebaVS
! File:
!-----------------------------------
! Author: Juan Diego Colmenares
!-----------------------------------
module ellipticFuselage

	use math_routines
	use surface_module

    	implicit none
    	private

    	integer :: i, j

    	public initEllipticGeom

    contains

    	function initEllipticGeom( length, nx, max_radius, na, axis, euler_a, origin, bvel, avel ) result( body )
        !
        !   This subroutine generates the nodes and location matrix for the
        !   elliptic fuselage.
        !
        	real(8), intent(in) :: length, &! Fuselage length in the X axis
                               		max_radius! Minor radius of the ellipsis
		real(8), dimension(3,1) :: wvel_aux
		real(8), dimension(3), intent(in) :: euler_a, origin, bvel, avel
        	integer, intent(in) :: nx, na, axis
		type( SurfaceType ) :: body

    	    	body%COORDSYS = reshape([1.0d+0, 0.0d+0, 0.d+0,&
	    			     	0.0d+0, 1.0d+0, 0.d+0,&
				     	0.0d+0, 0.0d+0, 1.d+0],[3,3])

	    	call euler_rot( 3, 1, 2, euler_a(1), euler_a(2), euler_a(3), body%COORDSYS )

    	    	body%ORIGIN = origin
    	    	body%VEL = bvel
    	    	body%WVEL = avel
	
	    	body%NN = (nx - 1) * na + 2
	    	body%NP = nx * na
	    	body%bluff = .true.

	    	body%VELF  => esferaVel
	    	body%WVELF => esferaWVel

		allocate( body%NODESXYZ ( 3,body%NN  ) )
		allocate( body%NODESLOC ( 3,body%NN  ) )
		allocate( body%CPXYZ    ( 3,body%NP  ) )
		allocate( body%CPUVW    ( 3,body%NP  ) )
		allocate( body%NORMAL   ( 3,body%NP  ) )
		allocate( body%MVEL     ( 3,body%NP  ) )
		allocate( body%GAMMAS   ( 4,body%NP  ) )
		allocate( body%NEWG     ( body%NP,1  ) )
		allocate( body%OLDG     ( body%NP,1  ) )
		allocate( body%LOAD     ( body%NP,1  ) )
		allocate( body%LOCMAT   ( 4,body%NP  ) )
		allocate( body%NBMAT    ( 4,body%NP  ) )

    	    	body%NODESLOC = ellipticNodes( length, max_radius, nx, na, body%NN )
		
		select case( axis )
		case( 2 )
			call euler_rot( 3, -90.d+0, body%NODESLOC )
		case( 3 )
			call euler_rot( 2, -90.d+0, body%NODESLOC )
		end select

    	    	body%NODESXYZ = matmul( body%COORDSYS, body%NODESLOC )

	    	forall( i=1:size(body%NODESXYZ,2) )
	    		body%NODESXYZ(:,i) = body%ORIGIN + body%NODESXYZ(:,i)
		endforall

		call ellipticPanels( body, nx, na )
    	    	body%NBMAT = ellipticNB( body%NP, na )

		call cPoints( body )

		wvel_aux(:,1) = body%WVELF( 0, avel )
		wvel_aux = matmul( body%COORDSYS, wvel_aux )

		do i = 1,body%NP
			body%CPUVW( :,i ) = body%VELF( 0, bvel ) + cross( wvel_aux(:,1), body%CPXYZ( :,i ) - body%ORIGIN )
		end do

		body%MVEL   = 0.d+0
		body%GAMMAS = 0.d+0
		body%NEWG   = 0.d+0
		body%OLDG   = 0.d+0
		body%LOAD   = 0.d+0

	end function initEllipticGeom

	function esferaWVel( timestep, vel_char ) result( vel )
		integer, intent(in) :: timestep
		real(8), dimension(3), intent(in) :: vel_char
		real(8), dimension(3) :: vel
		real(8) :: factor

		if( timestep >= 0 )then
		       factor = 1.d+0
		else
		       factor = 0.d+0
		endif

		vel = vel_char * factor		

	end function esferaWVel

	function esferaVel( timestep, vel_char ) result( vel )
		integer, intent(in) :: timestep
		real(8), dimension(3), intent(in) :: vel_char
		real(8), dimension(3) :: vel
		real(8) :: factor

		if( timestep >= 0 )then
		       factor = 1.d+0
		else
		       factor = 0.d+0
		endif

		vel = vel_char * factor		

	end function esferaVel

	function ellipticNodes( length, max_radius, nx, na, nn ) result( nodes )
		real(8), intent(in) :: length, max_radius
        	real(8), dimension(:), allocatable :: XOR, radius ! Linear space in the X axis
		real(8), dimension(:,:), allocatable :: nodes
        	real(8) :: dx, da
        	real(8) :: y, z, twopi
		integer, intent(in) :: nx, na, nn
		integer :: node_count

        	dx = pi/nx

        	allocate(XOR(nx+1))
        	allocate(radius(nx+1))
		allocate(nodes( 3,nn ))

        	!
        	! Do loop locates X coordinates near the nose and the tail of the fuselage
        	!
        	do i=1,nx+1
        	        XOR(i) =  - 0.5d+0 * length * cos((i-1)*dx)
        	        radius(i) = max_radius * sin((i-1)*dx)
        	end do

        	twopi = 2.d+0 * pi
        	da = twopi/na

        	node_count = 1

        	xpos: do i = 1,nx+1   ! loop over the positions on the X axis
        	    do j = 1,na
        	        if (i==1) then
        	            nodes(:,node_count) = [XOR(1), 0.0d+0, 0.0d+0]
        	            node_count = node_count + 1
        	            cycle xpos
        	        elseif (i==nx+1) then
        	            nodes(:,node_count) = [XOR(nx+1), 0.0d+0, 0.0d+0]
        	            node_count = node_count + 1
        	            cycle xpos
        	        else
        	            y = -radius(i)*cos((j-1)*da)
        	            z = radius(i)*sin((j-1)*da)
        	            nodes(:,node_count) = [XOR(i), y, z]
        	            node_count = node_count + 1
        	        endif
        	    end do
        	end do xpos

	end function ellipticNodes

	subroutine ellipticPanels( body, nx, na )
		class( SurfaceType ), intent(inout) :: body
		integer, intent(in) :: nx, na
		integer :: panel_count, nn

		nn = body%NN

        	panel_count = 1

        	do i=1,nx
        	    do j=1,na
        	        if (i==1) then
        	            if (j<na) then
        	                body%LOCMAT(1,panel_count)=1
        	                body%LOCMAT(2,panel_count)=j+1
        	                body%LOCMAT(3,panel_count)=j+2
        	                body%LOCMAT(4,panel_count)=1       !Node 1 repeats itself in the panel
        	                panel_count=panel_count+1
        	            else
        	                body%LOCMAT(1,panel_count)=1
        	                body%LOCMAT(2,panel_count)=j+1
        	                body%LOCMAT(3,panel_count)=2
        	                body%LOCMAT(4,panel_count)=1       !Node 1 repeats itself in the panel
        	                panel_count=panel_count+1
        	            end if
        	        else if (i==nx) then
        	            if (j<na) then
        	                body%LOCMAT(1,panel_count)=2+(i-2)*na+(j-1)
        	                body%LOCMAT(2,panel_count)=nn
        	                body%LOCMAT(3,panel_count)=nn
        	                body%LOCMAT(4,panel_count)=body%LOCMAT(1,panel_count)+1
        	                panel_count=panel_count+1
        	            else
        	                body%LOCMAT(1,panel_count)=2+(i-2)*na+(j-1)
        	                body%LOCMAT(2,panel_count)=nn       !The Final node also repeats itself
        	                body%LOCMAT(3,panel_count)=nn
        	                body%LOCMAT(4,panel_count)=nn-na
        	                panel_count=panel_count+1
        	            end if
        	        else
        	            if (j<na) then
        	                body%LOCMAT(1,panel_count)=2+(i-2)*na+(j-1)
        	                body%LOCMAT(2,panel_count)=body%LOCMAT(1,panel_count)+na
        	                body%LOCMAT(3,panel_count)=body%LOCMAT(2,panel_count)+1
        	                body%LOCMAT(4,panel_count)=body%LOCMAT(1,panel_count)+1
        	                panel_count=panel_count+1
        	            else
        	                body%LOCMAT(1,panel_count)=2+(i-2)*na+(j-1)
        	                body%LOCMAT(2,panel_count)=body%LOCMAT(1,panel_count)+na
        	                body%LOCMAT(3,panel_count)=body%LOCMAT(2,panel_count)-na+1
        	                body%LOCMAT(4,panel_count)=body%LOCMAT(3,panel_count)-na
        	                panel_count=panel_count+1
        	            end if
        	        end if
        	    end do
        	end do

	end subroutine ellipticPanels

    function ellipticNB(n,na) result(NB)
        !
        ! Function that generates the neighbour matrix for the fuselage body%LOCMAT.
        !
        integer, intent(in) :: n, na
        integer, dimension(:,:), allocatable :: NB

        allocate(NB(4,n))

        do i = 1,n
            if (i <= na) then
                NB(1,i) = mod((na+i-2),na) + 1
                NB(2,i) = i + na
                NB(3,i) = mod(i,na) + 1
                NB(4,i) = 0
            elseif (i > n - na) then
                NB(1,i) = n - na + mod((na+i-2),na) + 1
                NB(2,i) = 0
                NB(3,i) = n - na + mod((na+i),na) + 1
                NB(4,i) = i - na
            else
                NB(1,i) = ((i-1)/na)*na + mod(i-2,na) + 1
                NB(2,i) = i + na
                NB(3,i) = ((i-1)/na)*na + mod(i,na) + 1
                NB(4,i) = i - na
            endif
        enddo

    end function ellipticNB

end module ellipticFuselage
