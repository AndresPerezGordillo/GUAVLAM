module math_routines
!
!   	math_routines.f90
!------------------------------------
! Author: Juan D. Colmenares
!------------------------------------
!
!	This module implements auxiliary general-purpose subroutines and
!	functions that are not specific to the unsteady vortex lattice method. 
!
!---------------------------------------------------------------------

	implicit none
	public
	
	external :: DGEMM ! General Matrix-Matrix multiplication from ACML.
	!
	! Global constants
	!
	real(8), parameter :: pi = &
		0.31415926535897932384626433832795028841972d+01
   	real(8), parameter :: e = &
         	2.7182818284590452353602874713526624977572d+0
   	real(8), parameter :: degrees = &
         	180.0d+0
   	
   	interface rad2deg
     		module procedure rad2deg1, rad2deg2
   	end interface

  	interface deg2rad
     		module procedure deg2rad1, deg2rad2
  	end interface

  	interface eulerRot
     		module procedure euler3rot, euler1rot
   	end interface

   	interface cross
     		module procedure cross1, cross2
   	end interface

   	interface addorigin
     		module procedure addorigin_real1, &
				addorigin_real2,  &
				addorigin_integer
   	end interface

    	interface extend2mat
        	module procedure extend2mat1, extend2mat2
    	end interface extend2mat

 	contains

      pure function norm2( X ) result( n )
      !
      ! Function returning the L2 norm of a vector. 
      ! 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!
      ! 
      ! THIS SHOULD BE IMPLEMENTED ON FORTRAN 95 OR EARLIER VERSIONS OF FORTRAN,
      ! AS THIS IS INCLUDED AS AN INTRINSIC FUNCTION IN FORTRAN 2003.
      !
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!!!
	      real(8), dimension(:), intent(in) :: X
	      real(8) :: n

	      n = sqrt( dot_product(X,X) )

      end function norm2

    
	function rad2deg1(ang) result(tmp)
	! Converts angle in radians to degrees
	!
      		real, intent(in) :: ang
      		real :: tmp
      		tmp = ang*degrees/pi

   	end function rad2deg1
    
	function deg2rad1(ang) result(tmp)
	! Converts angle in degrees to radians
	!
      		real, intent(in) :: ang
      		real :: tmp
      		tmp = ang*pi/degrees
   	end function deg2rad1

	function rad2deg2(ang) result(tmp)
	! Converts angle in radians to degrees
	!
      		real(8), intent(in) :: ang
	      	real(8) :: tmp
      		tmp = ang*degrees/pi
	end function rad2deg2

   	function deg2rad2(ang) result(tmp)
	! Converts angle in degrees to radians
	!
      		real(8), intent(in) :: ang
      		real(8) :: tmp
      		tmp = ang*pi/degrees
   	end function deg2rad2

   	subroutine euler1rot(b1,ang1,M)
	!
	! Makes a rotation of the coordinates in the argument matrix M
	! Input parameters: Number of rotation, angle and body coordinates
	!
     		real(8), dimension(:,:), intent(in out) :: M
     		real(8), dimension(size(M,1),size(M,2)) :: M2
     		integer, intent(in) :: b1                    ! Rotation number (1,2 or 3)
     		real(8), intent (in) :: ang1                    ! angle in degrees
     		real(8), dimension(3,3) :: erot                ! Rotation matrix
     		real(8) :: tmp1, tmp2, zero, one
     		M2 = M
     		zero = 0.000d+0
     		one = 1.000d+0
     		tmp1=cos(deg2rad(ang1))
     		tmp2=sin(deg2rad(ang1))
     		select case (b1)
     		case (1)
     		    erot=reshape([one,zero,zero,zero,tmp1,tmp2,zero,-tmp2,tmp1],[3,3])
     		case (2)
     		    erot=reshape([tmp1,zero,-tmp2,zero,one,zero,tmp2,zero,tmp1],[3,3])
     		case (3)
     		    erot=reshape([tmp1,tmp2,zero,-tmp2,tmp1,zero,zero,zero,one],[3,3])
     		end select
     		call dgemm( 'n', 'n', 3, size(M2,2), 3, 1.d+0, erot, 3, M2, 3, 0.d+0, M, 3 )

  	end subroutine euler1rot

  	subroutine euler3rot(a1,a2,a3,ang1,ang2,ang3,M)
	!
	! Makes a rotation of the coordinates in the argument matrix M
	! Input parameters: rotation sequence, angles and body coordinates
	!
     		real(8), dimension(:,:), intent(in out) :: M
     		real(8), dimension(size(M,1),size(M,2)) :: M2
     		integer, intent(in) :: a1, a2, a3            ! Rotation sequence number (1,2 or 3)
     		integer, dimension(3) :: rotation
     		real(8), intent (in) :: ang1, ang2, ang3        ! angle in degrees
     		real(8), dimension(3) :: angles
     		real(8), dimension(3,3) :: erot, tempM, rotM    ! Rotation matrix
     		real(8) :: tmp1, tmp2, zero, one
     		integer :: i

     		M2 = M
     		zero = 0.000d+0
     		one = 1.000d+0
     		rotation(1)=a1
     		rotation(2)=a2
     		rotation(3)=a3
     		angles(1)=ang1
     		angles(2)=ang2
     		angles(3)=ang3

  		do i=1,3
  		    tmp1=cos(deg2rad(angles(i)))
  		    tmp2=sin(deg2rad(angles(i)))
  		   select case (rotation(i))
  		   case (1)
  		       erot=reshape([one,zero,zero,zero,tmp1,tmp2,zero,-tmp2,tmp1],[3,3])
  		   case (2)
  		       erot=reshape([tmp1,zero,-tmp2,zero,one,zero,tmp2,zero,tmp1],[3,3])
  		   case (3)
  		       erot=reshape([tmp1,tmp2,zero,-tmp2,tmp1,zero,zero,zero,one],[3,3])
  		   end select
  		 if (i==1) then
  		   rotM=erot
  		 else
  		   tempM=matmul(rotM,erot)
  		   rotM = tempM
  		 end if
  		end do
     		call dgemm( 'n', 'n', 3, size(M2,2), 3, 1.d+0, rotM, 3, M2, 3, 0.d+0, M, 3 )

  	end subroutine euler3rot

  	subroutine addorigin_real2(M,V)
  	  real(8), dimension(:), intent(in) :: V
  	  real(8), dimension(:,:), intent(in out) :: M
  	  integer :: i
  	  do i = 1,size(M(1,:))
  	    M(:,i)=M(:,i)+V
  	  end do
  	end subroutine addorigin_real2

  	subroutine addorigin_real1(M,V)
  	  real, dimension(:), intent(in) :: V
  	  real, dimension(:,:), intent(in out) :: M
  	   integer :: i
  	  do i = 1,size(M(1,:))
  	    M(:,i)=M(:,i)+V
  	  end do
  	end subroutine addorigin_real1

  	subroutine addorigin_integer(M,V)
  	  integer, dimension(:), intent(in) :: V
  	  integer, dimension(:,:), intent(in out) :: M
  	   integer :: i
  	  do i = 1,size(M(1,:))
  	    M(:,i)=M(:,i)+V
  	  end do
  	end subroutine addorigin_integer

  	subroutine extend2mat2(M1,M2)
	!
	! This subroutine takes a matrix and changes it's number of columns in order to
	! add the columns of a second matrix.
	!
    		real(8), dimension(:,:), intent(in out), allocatable :: M1
    		real(8), dimension(:,:), intent(in) :: M2
    		real(8), dimension(:,:), allocatable :: temp
    		allocate(temp(3,size(M1(1,:))))
    		temp=M1
    		deallocate(M1)
    		allocate(M1(size(temp(:,1)),size(temp(1,:))+size(M2(1,:))))
    		M1(:,1:size(temp(1,:)))=temp
    		M1(:,size(temp(1,:))+1:size(M1(1,:)))=M2
  	end subroutine extend2mat2

  	subroutine extend2mat1(M1,M2)
	!
	! This subroutine takes a matrix and changes it's number of columns in order to
	! add the columns of a second matrix.
	!
    		real, dimension(:,:), intent(in out), allocatable :: M1
    		real, dimension(:,:), intent(in) :: M2
    		real, dimension(:,:), allocatable :: temp
    		allocate(temp(3,size(M1(1,:))))
    		temp=M1
    		deallocate(M1)
    		allocate(M1(size(temp(:,1)),size(temp(1,:))+size(M2(1,:))))
    		M1(:,1:size(temp(1,:)))=temp
    		M1(:,size(temp(1,:))+1:size(M1(1,:)))=M2
  	end subroutine extend2mat1

	function eulerDerivatives(vang, W) result(der_ang)
	!
	! This function returns the derivatives of the euler angles.
	! !!!!!!!!!!! ONLY for the 3-1-2 rotation sequence !!!!!!!!!!
	!
	  	real(8), dimension(3), intent(in) :: W
	  	real(8), dimension(3), intent(in) :: vang
	  	real(8), dimension(3,1) :: der_ang_mat
	  	real(8), dimension(3)   :: der_ang
	  	real(8), dimension(3,1) :: w_aux
	  	real(8) :: st1,st2,st3,ct1,ct2,ct3
	  	real(8), dimension(3,3) :: Brot
	  	real(8), parameter :: zero = 0.0000d+0
	
	  	w_aux(:,1) = W

	  	st1 = sin(deg2rad(vang(1)))
	  	st2 = sin(deg2rad(vang(2)))
	  	st3 = sin(deg2rad(vang(3)))

	  	ct1 = cos(deg2rad(vang(1)))
	  	ct2 = cos(deg2rad(vang(2)))
	  	ct3 = cos(deg2rad(vang(3)))

	  	Brot = reshape([-st3, ct2*ct3, st2*st3, &
				zero, zero, ct2, 	&
				ct3, ct2*st3, -st2*ct3],&
				[3,3])
	  	Brot = Brot*(1.d+0/ct2)
	  	der_ang_mat = matmul(Brot,w_aux)
	  	der_ang(1) = rad2deg(der_ang_mat(1,1))
	  	der_ang(2) = rad2deg(der_ang_mat(2,1))
	  	der_ang(3) = rad2deg(der_ang_mat(3,1))

	end function eulerDerivatives

 function euler_func(a1,a2,a3,ang1,ang2,ang3,M) result(M2)
!
! Makes a rotation of the coordinates in the argument matrix M
! Input parameters: Number of rotations, angles and body coordinates
!
     real(8), dimension(:,:), intent(in) :: M
     real(8), dimension(:,:), allocatable :: M2
     integer, intent(in) :: a1, a2, a3            ! Rotation number (1,2 or 3)
     integer, dimension(3) :: rotation
     real(8), intent (in) :: ang1, ang2, ang3        ! angle in degrees
     real(8), dimension(3) :: angles
     real(8), dimension(3,3) :: erot, tempM, rotM    ! Rotation matrix
     real(8) :: tmp1, tmp2, zero, one
     integer :: i
     zero = 0.000d+0
     one = 1.000d+0
     rotation(1)=a1
     rotation(2)=a2
     rotation(3)=a3
     angles(1)=ang1
     angles(2)=ang2
     angles(3)=ang3
  do i=1,3
      tmp1=cos(deg2rad(angles(i)))
      tmp2=sin(deg2rad(angles(i)))
     select case (rotation(i))
     case (1)
         erot=reshape([one,zero,zero,zero,tmp1,tmp2,zero,-tmp2,tmp1],[3,3])
     case (2)
         erot=reshape([tmp1,zero,-tmp2,zero,one,zero,tmp2,zero,tmp1],[3,3])
     case (3)
         erot=reshape([tmp1,tmp2,zero,-tmp2,tmp1,zero,zero,zero,one],[3,3])
     end select
   if (i==1) then
     rotM=erot
   else
     tempM=matmul(rotM,erot)
     rotM = tempM
   end if
  end do
     allocate(M2(3,size(M(1,:))))
     call dgemm( 'n', 'n', 3, size(M2,2), 3, 1.d+0, rotM, 3, M, 3, 0.d+0, M2, 3 )
  end function euler_func

  	pure function cross1(a, b) result(prod)
  	  real(8), dimension(3), intent(in) :: a, b
  	  real(8), dimension(3) :: prod

  	  prod(1) = a(2)*b(3) - a(3)*b(2)
  	  prod(2) = a(3)*b(1) - a(1)*b(3)
  	  prod(3) = a(1)*b(2) - a(2)*b(1)

  	end function cross1

  	pure function cross2(a, b) result(prod)
  	  real(8), dimension(3,1), intent(in) :: a, b
  	  real(8), dimension(3,1) :: prod

  	  prod(1,1) = a(2,1)*b(3,1) - a(3,1)*b(2,1)
  	  prod(2,1) = a(3,1)*b(1,1) - a(1,1)*b(3,1)
  	  prod(3,1) = a(1,1)*b(2,1) - a(2,1)*b(1,1)

  	end function cross2

end module math_routines
