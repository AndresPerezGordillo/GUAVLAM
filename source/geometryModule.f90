module geometry_module
	!
	! This module contains the type definition of a 
	! linked list of instances of the SurfaceType
	! derived type.
	!
      	use surface_module

      	implicit none
      	public

      	type :: GeometryType
      	        class( SurfaceType ), allocatable :: body
      	        type( GeometryType ), pointer :: next => null()
      	end type

      	contains

	subroutine geom_create( self, body_data )
		type( GeometryType ), pointer :: self
		class( SurfaceType ), intent(in) :: body_data
		
		allocate( self )
		call geom_put_data( self, body_data )

	end subroutine geom_create

	subroutine geom_put_data( self, body_data )
	        type( GeometryType ), pointer :: self
	        class( SurfaceType ), intent(in)    :: body_data
	        
	        allocate( self%body, source=body_data )

      	end subroutine geom_put_data

	subroutine geom_add( self, body_data )
	        type( GeometryType ), pointer    :: self
		type( GeometryType ), pointer    :: old, current
		class( SurfaceType ), intent(in) :: body_data
		type( GeometryType ), pointer    :: element

		current => self
                do while( associated( current%next ) )
			old     => current
			current => current%next
		enddo

		call geom_create( element, body_data )
		element%next => current%next
		current%next => element

	end subroutine

	function geom_get_data( self, ind ) result( body_data )
		type( GeometryType ), intent(in), pointer :: self
		type( GeometryType ), pointer :: current, old
		class( SurfaceType ), pointer :: body_data
		integer, intent(in) :: ind
		integer :: counter

		counter = 1

		current => self

		do while( counter < ind )
			counter = counter + 1
			if( associated( current%next )) then
				old => current
				current => old%next
			else
				exit
			endif
		enddo

		body_data => current%body

	end function

	function geom_count( geom ) result( counter )
		type( GeometryType ), pointer :: geom
		type( GeometryType ), pointer :: current
		type( GeometryType ), pointer :: old
		integer :: counter

		counter = 0

		current => geom

		do while( associated( current ) )
			counter = counter + 1
			old => current
			current => old%next
			enddo

	end function geom_count

	subroutine geom_destroy( geom )
		type( GeometryType ), pointer :: geom
		type( GeometryType ), pointer :: current
		type( GeometryType ), pointer :: old

		current => geom

		do while ( associated( current ) )
			old => current
			current => old%next
			deallocate(old)
		enddo

	end subroutine

	function geom_total_np( geom ) result( globnp )
		type( GeometryType ), pointer :: geom
		type( GeometryType ), pointer :: current
		type( GeometryType ), pointer :: old
		integer :: globnp
		
		globnp = 0

		current => geom

		do while( associated( current ))

			if( current%body%bluff )then
				globnp = globnp + current%body%NP - 1
			else
				globnp = globnp + current%body%NP
			endif

			old  => current
			current => old%next
		enddo
		
	end function


end module geometry_module
