module uvlm_wrappers
!---------------------------------------------------------------------
!	uvlm_wrappers.f90
!-------------------------------------------------------------
!	
!       This module contains the subroutines used by the main program and serves
!       as an interface between the geometries and the basic functions and
!       subroutines from the uvlm_module.
!
!	Developed  by : Juan D. Colmenares
!	Department of Mechanical Engineering
!	Universidad de los Andes
!
!---------------------------------------------------------------------

	use uvlm_module
	use tecplot_wrappers

	implicit none
	public

	contains

	function geomTotalNp( geom ) result( tnp )
	!
	! This function returns the total number of panels which is the sum of
	! the 'panels' (vortex sheets) on each surface.
	!
		type( SurfaceType ), dimension(:), intent(in) :: geom
		integer :: tnp, i

		tnp = 0

		do i = 1, size(geom)
			if( geom(i)%bluff )then
				tnp = tnp + geom(i)%NP - 1
			else
				tnp = tnp + geom(i)%NP
			endif
		enddo

	end function

	subroutine generateGeomIndex( geom, geom_index )
	! 
	! Generates the geom_index array with the location indexes of the
	! 'panels' of each surface
	!
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: current
		integer, dimension(:,:), intent(inout) :: geom_index
		integer :: num_geoms, i, istart, ifinish

		num_geoms = size(geom)

		istart  = 1
		ifinish = 0

		do i=1,num_geoms
			current => geom(i)
			if( current%bluff )then
				ifinish = ifinish + current%NP - 1
			else
				ifinish = ifinish + current%NP
			endif
			geom_index( 1,i ) = istart
			geom_index( 2,i ) = ifinish
			istart = ifinish + 1
		enddo

	end subroutine generateGeomIndex

	subroutine GEOMAI( geom, geom_index, A_mat, timestep, deformable, core )
	!
	! Generates the influence matrix for the linear system of algebraic
	! equations. The influence coefficients corresponding to the influence
	! of a rigid surface on to itself are calculated only once and stay
	! constant throughout the simulation.  
	!
		type( SurfaceType ), dimension(:), intent(in), target  :: geom
		type( SurfaceType ), pointer      :: body1, body2
		integer, dimension(:,:), intent(in) :: geom_index
		real(8), dimension(:,:), intent(inout), target  :: A_mat
		real(8), dimension(:,:), pointer   :: A_p1, A_p2
		real(8), intent(in) :: core
		integer, intent(in) :: timestep
		integer :: num_geoms, i, j, istart1, istart2, ifinish1, ifinish2
		logical, intent(in) :: deformable

		num_geoms = size( geom )
		
		do i = 1,num_geoms 

			istart1  = geom_index( 1,i )
			ifinish1 = geom_index( 2,i )

			body1 => geom( i )

			do j = i,num_geoms 

				istart2  = geom_index( 1,j )
				ifinish2 = geom_index( 2,j )

				body2 => geom( j )

				A_p1 => A_mat( istart1:ifinish1 , istart2:ifinish2 )
				A_p2 => A_mat( istart2:ifinish2 , istart1:ifinish1 )

				if( i == j )then
					if( timestep < 1 .or. deformable )then
						call infMat( core, body1, A_p1 )
					endif
				else
					call infMat( core, body1, A_p1, body2, A_p2 )
				endif

			enddo

		enddo			

	end subroutine GEOMAI

	subroutine GEOMRHS( geom, geom_index, v_inf, timestep, core, RHS, KnecViscosity, Alpha, a1, delta_t, Viscous )
	!
	! Generates the right-hand side of the linear system of algebraic
	! equations.
	!
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: body_r, body_s
		real(8), dimension(3), intent(in) :: v_inf
		real(8), dimension(:,:), intent(inout), target 	:: RHS
		real(8), intent(in) 				:: core
		real(8), dimension(:,:), pointer 		:: RHS_p
		integer, dimension(:,:), intent(in) 	:: geom_index
		integer, intent(in) 			:: timestep
		integer, intent(in) 			:: Viscous
		integer :: i, j, istart, ifinish, num_geoms, nth
		real(8), intent(in) :: KnecViscosity
		real(8), intent(in) :: Alpha
		real(8), intent(in) :: a1
		real(8), intent(in) :: delta_t

		num_geoms = size( geom )
		nth = min(4,num_geoms)

		do i=1,num_geoms
			istart  = geom_index( 1,i )
			ifinish = geom_index( 2,i )

			body_r => geom( i )

			body_r%IVEL = 0.d+0

			RHS_p => RHS( istart:ifinish, : )

			call FSRHS( body_r, v_inf, RHS_p )

			do j=1,num_geoms

				body_s => geom( j )

				if( associated( body_s%wake ) )then

					if( timestep > 0 )then
					call WAKERHS( body_r, body_s, RHS_p, timestep, core )
					endif

				endif

			enddo

		enddo

	end subroutine GEOMRHS

	
	subroutine GEOMCONVECT( geom, V_inf, delta_t, timestep, core, KnecViscosity, Alpha, a1, Viscous )
	!
	! Calculates induced velocities on the wake, update the wake position 
	! and update the positions of the lifting surfaces.
	!
		type( SurfaceType ), dimension(:), intent(inout), target :: geom
		type( SurfaceType ), pointer :: body_r, body_s
		real(8), dimension(3), intent(in) :: V_inf
		real(8), dimension(3) :: V_infif
		real(8), intent(in) :: delta_t
		real(8), intent(in) :: core
		real(8), intent(in) :: KnecViscosity
		real(8), intent(in) :: Alpha
		real(8), intent(in) :: a1
		integer, intent(in) :: timestep
		integer, intent(in) :: Viscous
		integer :: wake_time
		integer :: ng, i, j, nth

		ng = size( geom )
		nth = min(4,ng)

		if( timestep > 0 ) then

		do i=1,ng
			
			body_r => geom( i )
			
			if( associated( body_r%wake ) )then

				body_r%wake%NODESUVW = 0.d+0
				wake_time = timestep - body_r%wake%START

				if( wake_time > 0 )then
					
				do j = 1,ng
					
					body_s => geom( j )

					if( i == j )then
					call inducedVelocity( body_r%wake%NODESXYZ, &
								body_r%wake%NODESUVW( :,1:body_r%wake%NEN*wake_time),&
								body_s, timestep-1, core, V_inf )
					else
					
					call inducedVelocity( body_r%wake%NODESXYZ, &
								body_r%wake%NODESUVW( :,1:body_r%wake%NEN*wake_time),&
								body_s, timestep-1, core )
					endif

				enddo

				endif
			endif

		enddo
		
		endif
		
		do i=1,ng

			body_r => geom( i )
			
			call updateCoordinates( body_r, delta_t, timestep )

			if( associated( body_r%wake ))then
				
				call convect( body_r, delta_t, timestep, KnecViscosity, Alpha, a1, Viscous, core )

			endif

		enddo


	end subroutine GEOMCONVECT

	subroutine GEOMG( geom, geom_index, GG )
	!
	! Pass circulations from the solution vector to the SurfaceType
	! variables.
	!
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: body
		integer, dimension(:,:), intent(in) :: geom_index
		real(8), dimension(:,:), intent(in), target :: GG
		real(8), dimension(:,:), pointer :: GI
		integer :: ng, i, istart, ifinish

		ng = size( geom_index, 2 )

		do i = 1,ng
			body => geom( i )
			istart = geom_index( 1,i )
			ifinish = geom_index( 2,i )

			GI => GG( istart:ifinish,: )

			call updateG( body, GI )
		enddo

	end subroutine GEOMG


	subroutine GEOMMVEL( geom, V_inf, timestep, core )
	!
	! Calculate mean induced velocity at each control point.
	!
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: body_r, body_s
		real(8), dimension(3), intent(in) :: V_inf
		real(8), intent(in) :: core
		integer, intent(in) :: timestep
		integer :: nb, i, j, k, m
		
		nb = size( geom )

		do i = 1,nb

			body_r => geom( i )

			forall( m = 1:body_r%NP )
				body_r%MVEL(:,m) = V_inf
			endforall

			do j=1,nb
				body_s => geom( j )

				do k=1,body_r%NP
					body_r%MVEL(:,k)= body_r%MVEL(:,k) + surfaceInfluence( body_r%CPXYZ(:,k), body_s, core )
				enddo
			enddo

			body_r%MVEL = body_r%MVEL + body_r%IVEL
		enddo

	end subroutine	GEOMMVEL


	subroutine GEOMLOADS( geom, delta_t )
	!
	! Compute pressure jump coefficients at each bound vortex sheet.
	!
		type( SurfaceType ), dimension(:), intent(inout), target :: geom
		type( SurfaceType ), pointer :: body
		real(8), intent(in) :: delta_t
		integer :: ng, i

		ng = size( geom )
		
		do i=1,ng
			body => geom( i )
			call LOADCALC( body, delta_t )
		enddo

	end subroutine GEOMLOADS


	subroutine GEOMTECIO( geom, zones, timestep_1 )
	!
	! Write data to ounput tecplot file.
	!
		type( SurfaceType ), dimension(:), intent(in), target :: geom
		type( SurfaceType ), pointer :: body
		character(len=*), dimension(:), intent(in) :: zones
		integer, intent(in) :: timestep_1
		integer :: ng, zone_counter, i, timestep

		zone_counter = 1

		ng = size(geom)


		do i=1,ng

			body => geom( i )

			StrandID = zone_counter
			call newSurfaceZone( trim(adjustl(zones(zone_counter))),&
					    body%NODESXYZ,&
					    body%NORMAL,&					
					    body%LOAD,&
					    body%LOCMAT,&
					    [1,1,1,0,0,0,0],&
					    1 )
			zone_counter = zone_counter + 1

			if( associated(body%wake) )then
				StrandID = zone_counter
				timestep = timestep_1 - body%wake%START

				if( timestep > 0 )then
				call newWakeZone( trim(adjustl(zones(zone_counter))),&
						    body%wake%NODESXYZ( :, 1:body%wake%NEN*(timestep + 1) ),&
						    body%wake%LOCMAT( :, 1:body%wake%NEP*timestep ),&
						    [1,1,1,0,0,0,0],&
						    1 )
				endif
				zone_counter = zone_counter + 1

			end if

		end do

	end subroutine

end module uvlm_wrappers
