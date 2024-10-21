module lapack_module
    implicit none
    
contains

    subroutine use_lapack_example
        integer :: n, info
        integer, dimension(2) :: ipiv
        real(8), dimension(2,2) :: A
        real(8), dimension(2) :: b

        ! Initialize the matrix and right-hand side
        n = 2
        A = reshape((/3.0_8, 1.0_8, 1.0_8, 2.0_8/), (/2,2/))
        b = (/5.0_8, 5.0_8/)

        ! Perform LU factorization of A
        call dgetrf(n, n, A, n, ipiv, info)

        ! Check if factorization was successful
        if (info /= 0) then
            print *, 'LU factorization failed with info = ', info
            stop
        endif

        ! Solve the system A*x = b using the LU factorization
        call dgetrs('N', n, 1, A, n, ipiv, b, n, info)

        ! Check if solving was successful
        if (info /= 0) then
            print *, 'Solving failed with info = ', info
            stop
        endif

        ! Print the solution vector x
        print *, 'Solution x: ', b
    end subroutine use_lapack_example
end module lapack_module