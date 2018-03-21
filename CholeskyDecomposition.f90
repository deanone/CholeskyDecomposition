program CholeskyDecomposition
implicit none
	integer, parameter :: read_unit = 99
	integer :: i, j, n, s, k
	real :: start, finish
	real, dimension(:, :), allocatable :: A
	real, dimension(:, :), allocatable :: L
	
	print*, "Enter the size of the array:"
	read*, n
	allocate(A(n, n))
	allocate(L(n, n))
	L ( : , : ) = 0.0
	
	! load the matrix from a .dat file
	open(read_unit, file = 'matrix.csv', status = 'old')
	do i = 1, n
		read(read_unit, *) A(i, :)
	end do
	close(read_unit)

	! Cholesky - –Banachiewicz algorithm
	call cpu_time(start)
	do i = 1, n
		do j = 1, n
			if (i == j) then
				s = 0
				do k = 1, j - 1
					s = s + L(j, k) ** 2
				end do
				L(i, j) = sqrt(A(i, j) - s)
			else if (i > j) then
				s = 0;
				do k = 1, j - 1
					s = s + L(i, k) * L(j, k)
				end do
				L(i, j) = (1.0 / L(j, j)) * (A(i, j) - s)
			end if
		end do
	end do
	call cpu_time(finish)
	
	deallocate(A)
	
	print *, finish - start
	! print the matrix
	do i = 1, n
		write(*, "(100g15.5)") (L(i, j), j = 1, n)
	end do
	
	deallocate(L)
end program CholeskyDecomposition