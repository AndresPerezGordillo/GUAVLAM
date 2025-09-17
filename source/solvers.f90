!-----------------------------------
! Author: Juan Diego Colmenares
!-----------------------------------
module solvers

    !use lapack95, only: gesv
    use math_routines, only: norm => norm2

    implicit none
    public


 	external :: DGESV


    contains

    	function bi_cg(A, b, x0, tol, maxiter) result(x)
        	real(8), dimension(:,:), intent(in) :: A
        	real(8), dimension(:,:), allocatable :: At
        	real(8), dimension(:,:), intent(in) :: b
        	real(8), dimension(:,:), intent(in) :: x0
        	real(8), intent(in), optional :: tol
        	integer, intent(in), optional :: maxiter
        	real(8) :: tol_local, beta, alfa, residual
        	real(8), dimension(2) :: rho
        	integer :: maxiter_local, i
        	real(8), dimension(:,:), allocatable :: x, r, rs, p, ps, q, qs

        	if (present(tol)) then
        	    tol_local = tol
        	else
        	    tol_local = 0.1d-08
        	endif

        	if (present(maxiter)) then
        	    maxiter_local = maxiter
        	else
        	    maxiter_local = 100
        	endif

        	allocate(x(size(x0,1),1))
        	allocate(r(size(x0,1),1))
        	allocate(rs(size(x0,1),1))
        	allocate(p(size(x0,1),1))
        	allocate(ps(size(x0,1),1))
        	allocate(q(size(x0,1),1))
        	allocate(qs(size(x0,1),1))
        	allocate(At(size(A,2),size(A,1)))

        	x = x0
        	At = transpose(A)
        	r = b - matmul(A,x0)
        	rs = r
        	rho = 0.0d+0

        	iteration: do i=1,maxiter_local
        	    rho(2) = rho(1)
        	    rho(1) = dot_product(r(:,1),rs(:,1))
        	    !if (rho(1) < tol_local) then
        	        !print*, 'Method Failed'
        	        !exit iteration
        	    !endif
        	    if (i==1) then
        	        p = r
        	        ps = rs
        	    else
        	        beta = rho(1)/rho(2)
        	        p = r + beta*p
        	        ps = rs + beta*ps
        	    endif

        	    q = matmul(A, p)
        	    qs = matmul(At, ps)
        	    alfa = dot_product(rs(:,1),r(:,1))/dot_product(ps(:,1),q(:,1))
        	    x = x + alfa*p
        	    r = r - alfa*q
        	    rs = rs - alfa*qs

        	    residual = norm(r(:,1))
        	    if (residual < tol_local) then
        	        exit iteration
        	    endif
        	end do iteration

        	if ((i==maxiter_local+1) .and. (.not. residual < tol_local)) then
        	    print*, 'Did not acheive convergence after ', i-1, ' iterations.'
        	endif

    	end function bi_cg


end module solvers
