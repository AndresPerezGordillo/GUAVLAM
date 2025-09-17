module geometry_uvlm_wrappers

	use geometry_module
	use uvlm
	use tecplot_wrappers

	implicit none
	public

	contains

	function generateGeomIndex( geom ) result( geom_index )
		type( GeometryType ), pointer        :: geom
		type( GeometryType ), pointer        :: current, old
		integer, dimension(:,:), allocatable :: geom_index
		integer :: num_geoms, i, istart, ifinish

		num_geoms = geom_count( geom )

		if( num_geoms > 0 )then
			allocate( geom_index( 2, num_geoms ) )
		else
			allocate( geom_index( 2, 0 ) )
			geom_index = 1
			return
		endif

		istart  = 1
		ifinish = 0

		current => geom
		
		do i=1,num_geoms
			if( current%body%bluff )then
				ifinish = ifinish + current%body%NP - 1
			else
				ifinish = ifinish + current%body%NP
			endif
			geom_index( 1,i ) = istart
			geom_index( 2,i ) = ifinish
			istart = ifinish + 1
			old     => current
			current => old%next
		enddo

	end function generateGeomIndex

	subroutine GEOMAI( geom, geom_index, A_mat, timestep, core, deformable )
		type( GeometryType ), pointer      :: geom
		class( SurfaceType ), pointer      :: body1, body2
		integer, dimension(:,:), intent(in) :: geom_index
		real(8), intent(in)                :: core
		real(8), dimension(:,:), intent(inout), target  :: A_mat
		real(8), dimension(:,:), pointer   :: A_p1, A_p2
		integer, intent(in) :: timestep
		integer :: num_geoms, i, j, istart1, istart2, ifinish1, ifinish2
		logical, intent(in) :: deformable

		num_geoms = geom_count( geom )
		
		do i = 1,num_geoms

			istart1  = geom_index( 1,i )
			ifinish1 = geom_index( 2,i )

			body1 => geom_get_data( geom, i )

			do j = i,num_geoms

				istart2  = geom_index( 1,j )
				ifinish2 = geom_index( 2,j )

				body2 => geom_get_data( geom, j )

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

	subroutine GEOMRHS( geom, geom_index, v_inf, timestep, core, RHS, debug )
		type( GeometryType ), pointer :: geom
		class( SurfaceType ), pointer :: body_r, body_s
		real(8), dimension(3), intent(in) :: v_inf
		real(8), dimension(:,:), intent(inout), target 	:: RHS
		real(8), intent(in) 				:: core
		real(8), dimension(:,:), pointer 		:: RHS_p
		integer, dimension(:,:), intent(in) 	:: geom_index
		integer, intent(in) 			:: timestep
		integer :: i, j, istart, ifinish, num_geoms
		logical, intent(in), optional :: debug
		logical :: debug_i

		if( present( debug ))then
			debug_i = debug
		else
			debug_i = .false.
		endif

		num_geoms = size( geom_index, 2 )

		do i=1,num_geoms
			istart  = geom_index( 1,i )
			ifinish = geom_index( 2,i )

			body_r => geom_get_data( geom, i )

			!print*, 'RHS ON BODY', i
			!print*

			RHS_p => RHS( istart:ifinish, : )

			call FSRHS( body_r, v_inf, RHS_p, debug_i )

			do j=1,num_geoms

				body_s => geom_get_data( geom, j )

				!print*, 'INFLUENCE OF BODY WAKE', J
				!print*

				select type( body_s )
				class is( LiftingSurfaceType )
					if( timestep > 0 )then
					call WAKERHS( body_r, body_s, RHS_p, timestep, core, debug_i )
					endif
				end select

			enddo

		enddo

	end subroutine GEOMRHS

	
	subroutine GEOMCONVECT( geom, euler_angles, V_inf, delta_t, timestep, core, forced )
		type( GeometryType ), pointer :: geom
		class( SurfaceType ), pointer :: body_r, body_s
		real(8), dimension(3), intent(in) :: V_inf
		real(8), dimension(:,:), intent(inout) :: euler_angles
		real(8), intent(in) :: delta_t
		real(8), intent(in) :: core
		real(8), dimension(3), intent(in), optional :: forced
		integer, intent(in) :: timestep
		integer :: ng, i, j

		ng = geom_count( geom )

		if( timestep > 0 ) then

		do i=1,ng
			
			body_r => geom_get_data( geom, i )

			!print*, 'CONVECTING WAKE', I

			select type( body_r )
			class is ( LiftingSurfaceType )

				body_r%wake%NODESUVW = 0.d+0

				do j = 1,ng
					
					body_s => geom_get_data( geom, j )

					!PRINT*, 'BY BODY', J

					if( i == j )then
					call inducedVelocity( body_r%wake%NODESXYZ, body_r%wake%NODESUVW( :,1:body_r%NEN*timestep),&
								body_s, timestep-1, core, V_inf )
					else
					call inducedVelocity( body_r%wake%NODESXYZ, body_r%wake%NODESUVW( :,1:body_r%NEN*timestep),&
								body_s, timestep-1, core )
					endif

				enddo

			end select
		enddo
		
		do i=1,ng

			body_r => geom_get_data( geom, i )
			
			select type( body_r )
			class is ( LiftingSurfaceType )
				call updateCoordinates( body_r, euler_angles( :,i ), delta_t, timestep )
				if( present( forced ) )then
					call convect( body_r, delta_t, timestep, forced=forced )
				else
					call convect( body_r, delta_t, timestep )
				endif
			type is ( SurfaceType )
			call updateCoordinates( body_r, euler_angles( :,i ), delta_t, timestep )
			end select

		enddo

		endif

	end subroutine GEOMCONVECT

	subroutine GEOMG( geom, geom_index, GG )
		type( GeometryType ), pointer :: geom
		class( SurfaceType ), pointer :: body
		integer, dimension(:,:), intent(in) :: geom_index
		real(8), dimension(:,:), intent(in), target :: GG
		real(8), dimension(:,:), pointer :: GI
		integer :: ng, i, istart, ifinish

		ng = size( geom_index, 2 )

		do i = 1,ng
			body => geom_get_data( geom, i )
			istart = geom_index( 1,i )
			ifinish = geom_index( 2,i )

			GI => GG( istart:ifinish,: )

			call updateG( body, GI )
		enddo

	end subroutine GEOMG


	subroutine GEOMIVEL( geom, V_inf, timestep, core )
		type( GeometryType ), pointer :: geom
		class( SurfaceType ), pointer :: body_r, body_s
		real(8), dimension(3), intent(in) :: V_inf
		real(8), intent(in) :: core
		integer, intent(in) :: timestep
		integer :: nb, i, j
		
		nb = geom_count( geom )

		do i = 1,nb
			body_r => geom_get_data( geom, i )

			body_r%MVEL = 0.d+0

			do j=1,nb
				body_s => geom_get_data( geom, j )
				if( i == j )then
				call inducedVelocity( body_r%CPXYZ, body_r%MVEL,&
							body_s, timestep, core, V_inf )
				else
				call inducedVelocity( body_r%CPXYZ, body_r%MVEL,&
							body_s, timestep, core )
				endif
			enddo
		enddo

	end subroutine	GEOMIVEL


	subroutine GEOMLOADS( geom, delta_t )
		type( GeometryType ), pointer :: geom
		class( SurfaceType ), pointer :: body
		real(8), intent(in) :: delta_t
		integer :: ng, i

		ng = geom_count( geom )
		
		do i=1,ng
			body => geom_get_data( geom, i )
			call LOADCALC( body, delta_t )
		enddo

	end subroutine GEOMLOADS


	subroutine GEOMTECIO( geom, zones, timestep )
		type( GeometryType ), pointer :: geom
		class( SurfaceType ), pointer :: body
		character(len=*), dimension(:), intent(in) :: zones
		character(len=100) :: timechar
		integer, intent(in) :: timestep
		integer :: ng, zone_counter, i

		zone_counter = 1

		ng = geom_count(geom)

		write(timechar,*) timestep
		timechar = adjustl(timechar)

		do i=1,ng

			body => geom_get_data( geom, i )

			select type( body )
			class is ( LiftingSurfaceType )

				StrandID = zone_counter
				call newTec360Zone( trim(adjustl(zones(zone_counter)))//trim(adjustl(timechar)),&
						    body%NODESXYZ,&
						    body%LOCMAT,&
						    body%LOAD,&
						    [1,1,1,0],&
						    1 )
				zone_counter = zone_counter + 1

				StrandID = zone_counter
				if( timestep > 0 )then
				call newTec360Zone( trim(adjustl(zones(zone_counter)))//trim(adjustl(timechar)),&
						    body%wake%NODESXYZ( :, 1:body%NEN*(timestep + 1) ),&
						    body%wake%LOCMAT( :, 1:body%NEP*timestep ),&
						    body%wake%LOOPG( 1:body%NEP*timestep,: ),&
						    [1,1,1,0],&
						    1 )
				endif
				zone_counter = zone_counter + 1

			type is ( SurfaceType )
				StrandID = zone_counter
				call newTec360Zone( trim(adjustl(zones(zone_counter)))//trim(adjustl(timechar)),&
						    body%NODESXYZ,&
						    body%LOCMAT,&
						    body%LOAD,&
						    [1,1,1,0],&
						    1 )
				zone_counter = zone_counter + 1
			end select

		end do

	end subroutine

end module geometry_uvlm_wrappers
