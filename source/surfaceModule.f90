module surface_module
!
!   	surfaceModule.f90
!------------------------------------
! Author: Juan D. Colmenares
!------------------------------------
!
!	This module declares the derived type SurfaceType for the definition of
!	lifting and non-lifting surfaces. Contains the subroutines and functions
!	used to modify the components of instances of this type.
!
!---------------------------------------------------------------------

	use wake_module
	use math_routines

    	implicit none
    	public

    	integer, private :: i,j

    	type :: SurfaceType

    	    real(8), dimension(3,3) :: COORDSYS
    	    real(8), dimension(3)   :: ORIGIN
    	    real(8), dimension(3)   :: VEL
    	    real(8), dimension(3)   :: WVEL
	    real(8), dimension(3)   :: EANGLES

    	    real(8), dimension(:,:), allocatable :: NODESXYZ
    	    real(8), dimension(:,:), allocatable :: NODESLOC
    	    real(8), dimension(:,:), allocatable :: CPXYZ
    	    real(8), dimension(:,:), allocatable :: CPUVW
    	    real(8), dimension(:,:), allocatable :: NORMAL
    	    real(8), dimension(:,:), allocatable :: AREA

    	    real(8), dimension(:,:), allocatable :: MVEL
    	    real(8), dimension(:,:), allocatable :: IVEL
    	    real(8), dimension(:,:), allocatable :: GAMMAS
    	    real(8), dimension(:,:), allocatable :: VORTEX
    	    real(8), dimension(:,:), allocatable :: NEWG
    	    real(8), dimension(:,:), allocatable :: OLDG
    	    real(8), dimension(:,:), allocatable :: LOAD

    	    integer, dimension(:,:), allocatable :: LOCMAT
    	    integer, dimension(:,:), allocatable :: NBMAT

	    type( wakeType ), pointer :: wake => null()

    	    integer :: NN, NP, NVS
	    integer :: SLOWS = 0
	    integer :: GROUP = 0

	    logical :: BLUFF = .false.

    	end type SurfaceType

    contains

	subroutine destroy_surface( geom, stat_var )
		type( SurfaceType ), intent(inout) :: geom
		integer, intent(out) :: stat_var
		

	    if( associated(geom%wake) )then
		print*, 'deallocating wake'
		deallocate( geom%wake, STAT=stat_var )
	    endif

	end subroutine

	subroutine LOADCALC( surface, delta_t )
		type( SurfaceType ), intent(inout) :: surface
		real(8), intent(in) 		    :: delta_t
		real(8), dimension(:,:), allocatable :: velJump
		real(8), dimension(3) :: velvec

		allocate( velJump( 3,surface%NP ) )

		call velocityJump( surface, velJump )

		if( surface%bluff )then
			!$omp parallel do private( velvec )
			do i=1,surface%NP
				velvec = surface%MVEL( :,i ) - surface%CPUVW( :,i ) - 0.5d+0 * velJump( :,i )
				surface%LOAD( i,1 ) = 1.d+0 - dot_product( velvec, velvec )
			end do
			!$omp end parallel do
		else
			!$omp parallel do
			do i=1,surface%NP
				surface%LOAD( i,1 ) = -2.0d+0 * ( ( surface%NEWG( i,1 ) - surface%OLDG( i,1 ) ) / delta_t &
						   + dot_product( surface%MVEL( :,i ) - surface%CPUVW( :,i ), velJump( :,i )))
			end do
			!$omp end parallel do
		endif

	end subroutine

	subroutine panelAreas( surface )
		type( SurfaceType ), intent(inout) :: surface
        	real(8), dimension(3,4) :: n, L
		real(8) :: area1, area2

        	do i = 1, surface%NP

        	    n(:,1) = surface%NODESXYZ( :,surface%LOCMAT(1,i) )
        	    n(:,2) = surface%NODESXYZ( :,surface%LOCMAT(2,i) )
        	    n(:,3) = surface%NODESXYZ( :,surface%LOCMAT(3,i) )
        	    n(:,4) = surface%NODESXYZ( :,surface%LOCMAT(4,i) )

        	    L(:,1) = n(:,2) - n(:,1)
        	    L(:,2) = n(:,3) - n(:,2)
        	    L(:,3) = n(:,4) - n(:,3)
        	    L(:,4) = n(:,1) - n(:,4)

        	    area1 = norm2( cross( L(:,1), L(:,2) ) )
        	    area2 = norm2( cross( L(:,3), L(:,4) ) )

        	    surface%AREA(i,1) = 0.5d+0 * (area1 + area2)

        	end do

	end subroutine

	subroutine velocityJump( surface, velJump )
		type( SurfaceType ), intent(in) :: surface
		real(8), dimension(:,:), intent(inout) :: velJump
        	real(8), dimension(3,4) :: n, L
        	real(8), dimension(3) :: gamma_arrow, deltaV, normal
		real(8), dimension(4) :: gam_loc
		real(8) :: area
		integer, dimension(4) :: NB, ind

		!$omp parallel do private(n, L, gamma_arrow, area, deltaV, gam_loc, normal, ind)
		do i = 1, surface%NP

		    ind = surface%LOCMAT(:,i)

        	    n(:,1) = surface%NODESXYZ( :,ind(1) )
        	    n(:,2) = surface%NODESXYZ( :,ind(2) )
        	    n(:,3) = surface%NODESXYZ( :,ind(3) )
        	    n(:,4) = surface%NODESXYZ( :,ind(4) )

        	    L(:,1) = n(:,2) - n(:,1)
        	    L(:,2) = n(:,3) - n(:,2)
        	    L(:,3) = n(:,4) - n(:,3)
        	    L(:,4) = n(:,1) - n(:,4)

		    gam_loc = surface%GAMMAS(:,i)

        	    L(:,1) = gam_loc(1) * L(:,1)
        	    L(:,2) = gam_loc(2) * L(:,2) 
        	    L(:,3) = gam_loc(3) * L(:,3) 
        	    L(:,4) = gam_loc(4) * L(:,4) 

		    area = surface%AREA(i,1)

        	    gamma_arrow = sum( L, 2 )

		    normal = surface%NORMAL(:,i)
        	    deltaV = 0.5d+0 * cross( normal, gamma_arrow ) / area

        	    velJump(:,i) = deltaV
		enddo
		!$omp end parallel do

   	end subroutine velocityJump


	subroutine convect( surface, deltaT, timestep_in, KnecViscosity, Alpha, a1, Viscous, core )
	    !
	    !  Subroutine to convect the wake nodes given their
	    !  local velocities. This also updates the loop circulations
	    !  of the wake and segment circulations
	    !
	    	type( SurfaceType ), intent(inout) :: surface
        	real(8), intent(in) 	:: deltaT
		real(8), intent(in) 	:: KnecViscosity
		real(8), intent(in) 	:: Alpha
		real(8), intent(in) 	:: a1
		real(8), intent(in) 	:: core
        	integer, intent(in) 	:: timestep_in
		integer, intent(in) 	:: Viscous
		real(8) :: TurbViscCoef
		integer			:: num_nodes, num_panels, nv_step
		integer :: timestep
		integer :: i, j, n, NV

		timestep = timestep_in - surface%wake%START

		num_nodes  = surface%wake%NEN * ( timestep )
		num_panels = surface%wake%NEP * ( timestep )
		nv_step = surface%wake%NEP + surface%wake%NEN

		surface%wake%NODESXYZ(:,1:num_nodes) = deltaT * ( surface%wake%NODESUVW( :,1:num_nodes ) ) &
				+ surface%wake%NODESXYZ( :,1:num_nodes )

		surface%wake%NODESXYZ = eoshift( surface%wake%NODESXYZ, -surface%wake%NEN, dim=2 )

		surface%wake%NODESXYZ( :,1:surface%wake%NEN ) = surface%NODESXYZ( :, surface%wake%TELSNI )

		surface%wake%NODESXYZ( :,surface%wake%CORNERS ) = surface%wake%NODESXYZ( :,surface%wake%CORMAP )

		surface%wake%LOOPG = eoshift( surface%wake%LOOPG, -surface%wake%NEP, dim=1 )

		surface%wake%LOOPG( 1:surface%wake%NEP,1 ) = surface%NEWG( surface%wake%TELSSI,1 )

		surface%wake%VORTEXG = eoshift( surface%wake%VORTEXG, -nv_step, dim=1 )

		call wakeVortexG( surface%wake )

!Inicio Correcion
		if ( Viscous==1 )then
			surface%wake%VORTEXC = core
			j = 0
			n = 1
			NV = INT ( nv_step*0.05/2 )
			if ( NV < 3 ) then
				NV = 3
				else
			end if
			do i = 1, nv_step * ( timestep )
				j = j+1
				if ( j > surface%wake%NEP - NV ) then
					if ( j <= surface%wake%NEP ) then
						TurbViscCoef = 1 + (a1 * abs(surface%wake%VORTEXG(i,1))/KnecViscosity)
						surface%wake%VORTEXC(i,1) = (core**(2d+0) + 4*Alpha*TurbViscCoef*KnecViscosity*n*deltaT)**(0.5d+0)
						else
					end if
					else
				end if
				if ( j > nv_step-NV ) then
					TurbViscCoef = 1 + (a1 * abs(surface%wake%VORTEXG(i,1))/KnecViscosity)
					surface%wake%VORTEXC(i,1) = (core**(2d+0) + 4*Alpha*TurbViscCoef*KnecViscosity*n*deltaT)**(0.5d+0)
					else
				end if
				if ( j >= nv_step )then
					j = 0
					n = n+1
					else
				end if
			end do
		else
			surface%wake%VORTEXC = core
		end if
!Fin Correccion

		

		call wakeVortexPos( surface%wake, timestep )

    	end subroutine convect

	subroutine gammas(NB, G, VTX, nodes, panels, gam, extra_i, extra_G)
        	!
        	! Calculate the circulations for each vortex segment
        	! --> gamma = G1 - G2
        	! --> G1: loop circulation at current vortex ring
        	! --> G2: loop circulation at neighboring vortex ring
        	!
        	!  Arguments:
        	! --> NB: Neighbor matrix
        	! --> G: loop circulations
        	!
        	real(8), dimension(:,:), intent(in) :: G, nodes
        	integer, dimension(:,:), intent(in) :: NB, panels
        	real(8), dimension(:,:), intent(in), optional :: extra_G
        	integer, dimension(:), intent(in), optional :: extra_i
        	integer :: nup
        	real(8), dimension(:,:), intent(inout) :: gam, VTX
        	integer, dimension(4) :: ne
        	integer :: nei, k
        	real(8) :: G1, G2
        	integer :: wake_index, counter

        	wake_index = 1
		counter = 1

        	nup = size(NB,2)

        	! Initialize the gamma values as 0.0
        	gam = 0.0

        	do i = 1,nup
        	    ! neighbors of the panel i
        	    ne = NB(:,i)
        	    G1 = G(i,1)
        	    ! loop over the neighbors
        	    do j = 1,4
        	        nei = ne(j)

        	        if (nei /= 0) then
        	            if ((nei < -1) .or. (nei > size(G,1))) then
        	                G2 = 0.0
        	            elseif (nei == -1) then

        	                if (present(extra_i) .and. present(extra_G)) then
        	                    G2 = extra_G(extra_i(wake_index),1)
        	                    wake_index = wake_index + 1
        	                else
        	                    G2 = 0.0
        	                end if

        	            else

        	                G2 = G(nei,1)

        	            end if
        	            ! circulation value of segment 'j' at panel 'i'.
        	            gam(j,i) = G1 - G2

			    if( nei < i )then
				    k = mod(j,4)+1
				if( nei < 0 )then
				    	VTX(1,counter)  = G1
				else
				    	VTX(1,counter)  = G1 - G2
				endif
				    VTX(2:4,counter)= nodes(:,panels(j,i))
				    VTX(5:7,counter)= nodes(:,panels(k,i))
				    counter = counter + 1
			    endif
		     	else
			    gam(j,i) = 0.d+0
		    	end if

			if(nei == -2)then
				select case(j)
				case(3)
					gam(j,i) = 2.d+0*gam(j,i)
				end select
			endif

        	    end do
        	end do

    	end subroutine gammas


	subroutine updateG( surface, GG )
		type(SurfaceType), intent(inout) :: surface
		real(8), dimension(:,:), intent(in) :: GG

		surface%OLDG = surface%NEWG
		surface%NEWG(1:size(GG,1),:) = GG

		if( associated( surface%wake ))then
			 call gammas( surface%NBMAT, surface%NEWG, &
						surface%VORTEX, surface%NODESXYZ,&
						surface%LOCMAT,surface%GAMMAS,&
						surface%wake%TEWAKEI, surface%wake%LOOPG )
		else
			call gammas( surface%NBMAT, surface%NEWG, &
						surface%VORTEX, surface%NODESXYZ,&
						surface%LOCMAT, surface%GAMMAS)
		endif
	
	end subroutine updateG


	subroutine updateCoordinates( surface, delta_t, timestep )
		type( SurfaceType ), intent(inout), target :: surface
		real(8), dimension(3) :: euler_angs
		real(8), dimension(3) :: k1, k2, k3, k4, euler_aux
		real(8), dimension(3,1) :: wvel_aux, vel_aux
		real(8), intent(in) :: delta_t
		integer, intent(in) :: timestep
		integer :: slow_start

		slow_start = surface%SLOWS

		if( timestep >= 0 .and. timestep + 1 < slow_start )then
			wvel_aux(:,1) = surface%WVEL * (timestep + 1.)/slow_start 
		else
			wvel_aux(:,1) = surface%WVEL
		endif

		euler_angs = surface%EANGLES

		k1 = eulerDerivatives( euler_angs, wvel_aux(:,1) )
		euler_aux = euler_angs + 0.5d+0 * delta_t * k1
		k2 = eulerderivatives( euler_aux, wvel_aux(:,1) )
		euler_aux = euler_angs + 0.5d+0 * delta_t * k2
		k3 = eulerderivatives( euler_aux, wvel_aux(:,1) )
		euler_aux = euler_angs + delta_t * k3
		k4 = eulerderivatives( euler_aux, wvel_aux(:,1) )

		euler_angs = euler_angs + delta_t * ( k1 + 2.d+0*k2 + 2.d+0*k3+ k4) / 6.d+0

		surface%COORDSYS = reshape([1.d+0, 0.d+0, 0.d+0,&
			    		0.d+0, 1.d+0, 0.d+0,&
			    		0.d+0, 0.d+0, 1.d+0],[3,3])

		call eulerRot( 3, 1, 2, euler_angs(1), euler_angs(2), euler_angs(3), surface%COORDSYS )

		surface%EANGLES = euler_angs

		surface%NODESXYZ = matmul( surface%COORDSYS, surface%NODESLOC ) 

		surface%ORIGIN = surface%ORIGIN + delta_t * surface%VEL

		forall( i = 1:surface%NN )
			surface%NODESXYZ( :,i ) = surface%NODESXYZ( :,i ) + surface%ORIGIN
		end forall

		call cPoints( surface )

		wvel_aux = matmul( surface%COORDSYS, wvel_aux )

		forall( i = 1:surface%NP )
			surface%CPUVW( :,i ) = surface%VEL + cross( wvel_aux(:,1), surface%CPXYZ( :,i ) - surface%origin )
		end forall

	end subroutine


    	subroutine cPoints( surface )
        	!
        	! Calculate control points and normal vectors on a set of panels.
        	!
		type(SurfaceType), intent(inout), target :: surface
        	real(8), dimension(3,2) 	:: diag
        	real(8), dimension(3,1) 	:: pc1, pc2, pc3, pc4
        	real(8), dimension(3) 		:: cross_product
        	integer, pointer 		:: p1, p2, p3, p4

        	do i = 1,surface%NP

            	p1 => surface%LOCMAT( 1,i )
            	p2 => surface%LOCMAT( 2,i )
            	p3 => surface%LOCMAT( 3,i )
            	p4 => surface%LOCMAT( 4,i )
            	pc1(:,1) = surface%NODESXYZ( :,p1 )
            	pc2(:,1) = surface%NODESXYZ( :,p2 )
            	pc3(:,1) = surface%NODESXYZ( :,p3 )
            	pc4(:,1) = surface%NODESXYZ( :,p4 )

            	surface%CPXYZ( :,i ) = quadPoint(pc1, pc2, pc3, pc4)
            	diag = diagonals(pc1, pc2, pc3, pc4)
            	cross_product = cross(diag(:,1),diag(:,2))
            	surface%NORMAL( :,i ) = cross_product / norm2( cross_product )

        	end do

    	end subroutine cPoints

    function triangCPoint(c1, c2, c3) result(cpcoord)
        !
        ! Returns an array with the coordinates of the control point
        ! at a triangular panel.
        !
        ! Arguments:
        ! --> c1, c2, c3: XYZ coordinates of the first, second and third
        !                 node respectively.
        !
        real(8), dimension(3,1), intent(in) :: c1, c2, c3
        real(8), dimension(3) :: cpcoord

        forall (i=1:3)
            cpcoord(i) = (c1(i,1)+c2(i,1)+c3(i,1))/3
        end forall

    end function triangCPoint

    function quadPoint(c1, c2, c3, c4) result(cpcoord)
        !
        !   Calculate control point at a quadrilateral panel. Distinguishes between
        !   a rectangular ring and a triangular one, taking the correct action to
        !   calculate its centroid.
        !
        ! Arguments:
        ! --> C1, c2, c3, c4: coordinates of node 1, 2, 3 and 4 respectively.
        !
        ! Output:
        ! --> cpcoord: coordinates of the control point.
        !
        real(8), dimension(3,1), intent(in) :: c1, c2, c3, c4
        real(8), dimension(3) :: cpcoord
        logical, dimension(3,1) :: equiv1, equiv2, equiv3, equiv4

        equiv1 =  c1 == c2
        equiv2 =  c1 == c4
        equiv3 =  c3 == c2
        equiv4 =  c3 == c4

        if ((all(equiv1(:,1))) .or. &
        (all(equiv2(:,1)))) then
            cpcoord = triangCPoint(c2, c3, c4)
        elseif ((all(equiv3(:,1))) .or. &
        (all(equiv4(:,1)))) then
            cpcoord = triangCPoint(c1, c2, c4)
        else
            forall (i=1:3)
                cpcoord(i) = (c1(i,1)+c2(i,1)+c3(i,1)+c4(i,1))*0.25d+0
            end forall
        end if

    end function quadPoint

    function diagonals(c1, c2, c3, c4) result(d)
        real(8), dimension(3,1), intent(in) :: c1, c2, c3, c4
        real(8), dimension(3,2) :: d

        d(:,1) = c3(:,1) - c1(:,1)
        d(:,2) = c4(:,1) - c2(:,1)
    end function diagonals

    
end module surface_module
