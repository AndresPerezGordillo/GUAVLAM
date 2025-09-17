module uvlm_module
!-----------------------------------
!  uvlm.f90
!-----------------------------------
! Author: Juan Diego Colmenares
!-----------------------------------
!
!	Contains basic subroutines and functions for the unsteady vortex lattice
!	method using the user defined types SurfaceType and WakeType.
!
!---------------------------------------------------------------------
	use surface_module
	use omp_lib

    	implicit none
   	public
	

    contains

	subroutine inducedVelocity( points, vels, surface, time, core, v_inf )
	! 
	! Calculate the velocity induced by lifting surface and its wake at a
	! given set of coordinates.
	!
	! Arguments: 
	! -> points: 	3xN array containing the X,Y and Z coordinates of N
	! 		points where induced velocities will be computed.
	! -> vels:	3xN array containing the U,V and W components of the
	! 		induced velocities at the specified coordinates.
	! -> surface: 	Variable of type SurfaceType representing the 'sending'
	! 		surface.
	! -> time: 	Current time step.
	! -> core:	Viscous core radius.
	! -> v_inf: 	Free stream velocity.
	!
		type( SurfaceType ), intent(in) :: surface
		real(8), dimension(:,:), intent(in   ) :: points
		real(8), dimension(:,:), intent(inout) :: vels
		real(8), dimension(3), intent(in), optional :: v_inf
		real(8), dimension(3) :: vi, vp
		real(8), intent(in) :: core
		integer, intent(in) :: time
		real(8), dimension(3) :: p
		integer :: num_points, i

		num_points = size( vels, 2 )

		if( present( v_inf ) ) then
			vi = v_inf
		else
			vi = 0.d+0
		endif

		do i=1,num_points
				
			p = points(:,i)
			vp = vels(:,i)
			vp = surfaceInfluence( p, surface, core ) + vi + vp
			
			if( associated( surface%wake ) )then
				
				vp = wakeInfluence( p, surface%wake, time, core ) + vp

			endif

			vels(:,i) = vp

		end do

	end subroutine
	

	pure function fakeBiot( p, n1, n2, core ) result ( ivel )
        !
        ! Biot-savart for a finite segment at a point 'p'
        !
        ! Arguments:
        !-> p: 		coordinates of the recieving point.
        !-> n1: 	coordinates at the start of the segment.
        !-> n2: 	coordinates at the end of the segment.
	!-> core: 	viscous core radius.
        !
        	real(8), dimension(:), intent(in) :: p
        	real(8), dimension(:), intent(in) :: n1, n2
        	real(8), intent(in) :: CORE
        	real(8), dimension(3) :: r0, r1, r2
        	real(8), dimension(3) :: direction
        	real(8) :: aux
        	real(8), dimension(3) :: ivel
        	real(8) :: num, den
		real(8), parameter :: CUTOFF = 1.0d-8

        	r1 = p - n1

        	r2 = p - n2

        	r0 = n2 - n1

        	direction = cross(r0, r1)

		num = norm2(direction)
		den = norm2(r0)

		if ( num < CUTOFF .or. den < CUTOFF )  then

			ivel = 0.d+0
		
		else
        		ivel = direction / ( ( num**4 + (den*core)**4 )**(0.5d+0) )
        		aux = dot_product(r0, r1/norm2(r1) - r2/norm2(r2))
        		ivel =  ivel * aux 
		endif

	end function fakeBiot

	function surfaceInfluence( point, surface, core ) result( vel )
	!
	! Function returning the velocity induced by a set of bound vortex
	! sheets of a lifting surface at a given point in space.
	!
		type( SurfaceType ), intent(in) :: surface
		real(8), dimension(3), intent(in) :: point
		real(8), dimension(:,:), allocatable :: vel_array
		real(8), dimension(6) :: n
		real(8), dimension(3) :: vel, velp
		real(8), intent(in) :: core
		real(8) :: cx, cy, cz
		real(8) :: gg
		integer, dimension(4) :: ind
		integer :: i

		cx = 0.d+0
		cy = 0.d+0
		cz = 0.d+0

		!$omp parallel do default(shared) private(gg, n, velp) reduction(+:cx,cy,cz)
		do i = 1,surface%NVS
			gg = surface%VORTEX(1,i)
			n = surface%VORTEX(2:7,i)
			velp = fakeBiot(point,n(1:3),n(4:6),CORE)*(0.25d+0*gg/pi)
			cx = velp(1) + cx
			cy = velp(2) + cy
			cz = velp(3) + cz
		enddo
		!$omp end parallel do

		vel(1) = (cx)
		vel(2) = (cy)
		vel(3) = (cz)

	end function surfaceInfluence


	function wakeInfluence( point, wake, time_in, core ) result( vel )
	!
	! Function returning the velocity induced by a set of free vortex
	! sheets of a wake structure at a given point in space.
	!
		type( WakeType ), intent(in) :: wake
		real(8), dimension(3), intent(in) :: point
		real(8), intent(in) :: core
		integer, intent(in) :: time_in

		real(8), dimension(:,:), allocatable :: vel_array
		real(8), dimension(6) :: n
		real(8), dimension(3) :: vel, velp
		real(8) :: gg, vx, vy, vz
		real(8) :: rc
		integer :: ii, nen, nep, nv
		integer :: time
		
		

		vx = 0.d+0
		vy = 0.d+0
		vz = 0.d+0

		nen = wake%NEN
		nep = wake%NEP

		time = time_in - wake%START

		nv = (time + 1)*nep + time*nen
		
		if ( time > 0 ) then
			
			!$OMP parallel do default(shared) private(gg,n,velp) reduction(+:vx,vy,vz)
			do ii = 1, nv
				gg = wake%VORTEXG(ii,1)

				rc = wake%VORTEXC(ii,1)! Vortex Core Corregido
				
				n = wake%VORTEXP(:,ii)

				velp =( 0.25d+0 * gg / pi ) * fakeBiot( point, n(1:3), n(4:6), rc )
				vx = velp(1) + vx
				vy = velp(2) + vy
				vz = velp(3) + vz
			enddo
			!$OMP end parallel do

			vel(1) = vx
			vel(2) = vy
			vel(3) = vz

		endif

	end function

	subroutine FSRHS( surface_r, v_inf, RHS )
	!
	! Computes the normal component of the free stream velocity at a surface
	! and adds this term to the right-hand side vector of the linear system
	! of algebraic equations.
	!
		type(SurfaceType), intent(in) :: surface_r
		real(8), dimension(3), intent(in) :: v_inf
		real(8), dimension(:,:), intent(inout) :: RHS
		real(8), dimension(3) :: rel_vel
		integer :: i

		do i= 1,size(RHS,1)
			rel_vel = v_inf - surface_r%CPUVW( :,i )
			RHS( i,1 ) = RHS( i,1 ) - dot_product( rel_vel, surface_r%NORMAL(:,i) )
		end do

	end subroutine FSRHS
		
	
	subroutine WAKERHS( surface_r , surface_s, RHS, timestep, core )
	!
	! Computes the normal component of the velocity induced by a wake
	! structure at the control points of a surface
	! and adds this term to the right-hand side vector of the linear system
	! of algebraic equations.
	!
		type(SurfaceType), intent(inout)  :: surface_r
		type(SurfaceType), intent(in)  :: surface_s
        	real(8), dimension(:,:), intent(inout) :: RHS
		real(8), dimension(3) 	:: VIND
        	real(8), intent(in)     :: core
	        integer, intent(in)     :: timestep
        	integer 		:: num_points, i

       		num_points = size( surface_r%CPXYZ, 2 )

		if( surface_r%bluff )then
			num_points = num_points - 1
		endif

        	do i = 1,num_points

			VIND = wakeInfluence( surface_r%CPXYZ( :,i ), surface_s%wake, timestep, core )
			surface_r%IVEL(:,i) = surface_r%IVEL(:,i) + VIND
            		RHS( i,1 ) = RHS( i,1 ) - dot_product( VIND, surface_r%NORMAL( :,i ) )

        	end do

    	end subroutine WAKERHS
 
   

    	function influenceMatrix(cpointsA, normalA, nodesB, panelsB, npA, npB, core) result(Mat)
    	    !
    	    ! Build influence matrix for a set of panels on a set of control points
    	    ! The panels and the control points don't need to be on the same body.
    	    !
    	    ! Arguments:
    	    !
    	    ! --> cpointA: control points at receiving panels
    	    ! --> normalA: normal vectors at receiving panels
    	    ! --> nodesB: nodes of the sending panels
    	    ! --> panelsB: location matrix for the sending panels
    	    ! --> NBB: neighbor matrix for the sending panels
    	    !
    	    real(8), dimension(:,:), intent(in) :: nodesB
    	    integer, dimension(:,:), intent(in) :: panelsB
    	    real(8), dimension(:,:), intent(in) :: cpointsA
    	    real(8), dimension(:,:), intent(in) :: normalA
    	    real(8), dimension(:,:), allocatable :: Mat
	    real(8), intent(in) :: core
    	    real(8), dimension(3) :: normal, point
    	    real(8), dimension(3,4) :: n, velp
    	    real(8), dimension(3) :: b
    	    real(8) :: c
	    integer, intent(in) :: npA, npB
    	    integer :: i, j

    	    allocate(Mat(npA,npB))

    	    ! initialize all the elements of Mat as 0.0
    	    Mat = 0.0

    	    do i = 1,npA ! Receiving Panels

    	        normal = normalA(:,i)
    	        point = cpointsA(:,i)

    	        do j = 1,npB ! Sending Panels

    	            	n(:,1) = nodesB(:,panelsB(1,j))
    	            	n(:,2) = nodesB(:,panelsB(2,j))
    	            	n(:,3) = nodesB(:,panelsB(3,j))
    	            	n(:,4) = nodesB(:,panelsB(4,j))

		    	velp(:,1) = fakeBiot( point, n(:,1), n(:,2), core ) 
			velp(:,2) = fakeBiot( point, n(:,2), n(:,3), core ) 
			velp(:,3) = fakeBiot( point, n(:,3), n(:,4), core ) 
			velp(:,4) = fakeBiot( point, n(:,4), n(:,1), core ) 

			b = sum( velp, 2 ) * ( 0.25d+0 / pi )

			c = dot_product( normal, b )

			Mat( i,j ) = c
    	            
    	        enddo

    	    enddo

    	end function influenceMatrix

    	subroutine infMat( core, bodyA, B_to_A, bodyB, A_to_B )
    	    !
    	    ! Subroutine to generate the influence matrices between panels from body 'A' and 'B'.
    	    !
	    type( SurfaceType ), pointer		:: bodyA
	    type( SurfaceType ), pointer, optional  	:: bodyB
	    real(8), dimension(:,:), pointer		:: B_to_A
	    real(8), dimension(:,:), pointer, optional	:: A_to_B
	    real(8), intent(in) :: core
	    integer :: npA, npB


    	    if (present( bodyB ) .and. present( A_to_B ))then

		if(bodyA%bluff)then
			npA = bodyA%NP - 1
		else
			npA = bodyA%NP
		endif

		if(bodyB%bluff)then
			npB = bodyB%NP - 1
		else
			npB = bodyB%NP
		endif

    	        B_to_A = influenceMatrix( bodyA%CPXYZ, bodyA%NORMAL,&
					  bodyB%NODESXYZ, bodyB%LOCMAT, npA, npB, CORE )
    	        A_to_B = influenceMatrix( bodyB%CPXYZ, bodyB%NORMAL,&
					  bodyA%NODESXYZ, bodyA%LOCMAT, npB, npA, CORE )
    	    else

		if(bodyA%bluff)then
			npA = bodyA%NP - 1
		else
			npA = bodyA%NP
		endif

    	        B_to_A = influenceMatrix( bodyA%CPXYZ, bodyA%NORMAL,&
					  bodyA%NODESXYZ, bodyA%LOCMAT, npA, npA, CORE )
    	    end if

    	end subroutine infMat

end module uvlm_module
