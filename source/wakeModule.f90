module wake_module

    implicit none

    integer, private :: i,j

    type :: WakeType

        real(8), dimension(:,:), allocatable :: NODESXYZ
	real(8), dimension(:,:), allocatable :: NODESUVW
	real(8), dimension(:,:), allocatable :: LOOPG
	real(8), dimension(:,:), allocatable :: VORTEXG
	real(8), dimension(:,:), allocatable :: VORTEXP
	real(8), dimension(:,:), allocatable :: VORTEXC! Vector de vortex core

        integer, dimension(:,:), allocatable :: LOCMAT
        integer, dimension(:,:), allocatable :: VORTEXPI
	integer, dimension(:)  , allocatable :: TEWAKEI
	integer, dimension(:)  , allocatable :: CORNERS
	integer, dimension(:)  , allocatable :: CORMAP

    	integer, dimension(:), allocatable :: TELSSI, TELSNI
        integer :: NEN, NEP, NVS
	integer :: START = 0

	logical :: ROOT = .false.

    end type WakeType

    contains

   	subroutine wakeVortexI( wake, wake_size )
		type( WakeType ), intent(inout) :: wake
		integer, intent(in) :: wake_size
		integer :: nep, nen, ns, t, counter

		nep = wake%NEP
		nen = wake%NEN
		ns = wake%NVS

		counter = 1

		do t = 1, wake_size
			do i = 1,nep
				wake%VORTEXPI(1,counter) = (t-1)*nen + i
				wake%VORTEXPI(2,counter) = (t-1)*nen + i + 1
				counter = counter + 1
			enddo
			do i = 1,nen
				wake%VORTEXPI(1,counter) = t*nen + i
				wake%VORTEXPI(2,counter) = (t-1)*nen + i
				counter = counter + 1
			enddo
		enddo

		do i = 1,nep
			wake%VORTEXPI(1,counter) = wake_size * nen + i
			wake%VORTEXPI(2,counter) = wake_size * nen + i + 1
			counter = counter + 1
		enddo

	end subroutine


	subroutine wakeVortexG( wake )
		type( WakeType ), intent(inout) :: wake
		integer :: nep, nen, counter

		nep = wake%NEP
		nen = wake%NEN

		counter = 1

		do i = 1,nep
			wake%VORTEXG(i,1) = -wake%LOOPG(i,1)
			counter = counter + 1
		end do

		if( wake%ROOT )then
			wake%VORTEXG(counter,1) = 0.d+0
		else
			wake%VORTEXG(counter,1) = -wake%LOOPG(1,1)
		endif

		counter = counter + 1

		do i = 2,nep
			wake%VORTEXG(counter,1) = wake%LOOPG(i-1,1) &
						- wake%LOOPG(i,1)
			counter = counter + 1
		enddo

		wake%VORTEXG(counter,1) = wake%LOOPG(nep,1)
		counter = counter + 1

		do i = 1,nep
			wake%VORTEXG(counter,1) = wake%LOOPG(i,1) &
						- wake%LOOPG(i+nep,1)
			counter = counter + 1
		enddo

	end subroutine

	subroutine wakeVortexPos( wake, t )
		type( WakeType ), intent(inout) :: wake
		integer, intent(in) :: t
		integer :: nvs_i

		nvs_i = wake%NEP * (t+1) + wake%NEN * t
		
		do i = 1,nvs_i
			wake%VORTEXP(1:3,i) = wake%NODESXYZ(:,wake%VORTEXPI(1,i))
			wake%VORTEXP(4:6,i) = wake%NODESXYZ(:,wake%VORTEXPI(2,i))
		enddo

	end subroutine

    subroutine wakePanels( wake, wake_size)
        type( wakeType ), intent(inout) :: wake
        integer, intent(in) ::  wake_size
        integer :: NEP

        NEP = wake%NEP

        do i = 1, wake_size
            do j = 1, wake%NEP
                wake%LOCMAT(1, (i-1)*NEP + j) = i * (NEP + 1) + j
                wake%LOCMAT(2, (i-1)*NEP + j) = i * (NEP + 1) + j + 1
                wake%LOCMAT(3, (i-1)*NEP + j) = (i-1) * (NEP + 1) + j + 1
                wake%LOCMAT(4, (i-1)*NEP + j) = (i-1) * (NEP + 1) + j
            end do
        end do

    end subroutine wakePanels

end module wake_module
