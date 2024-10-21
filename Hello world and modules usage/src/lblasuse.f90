module blas_module
    implicit none
contains
    subroutine use_blas_example

        integer :: n
        real(8) :: alpha, beta
        real(8), dimension(2,2) :: A
        real(8), dimension(2) :: x, y

        ! Initialize values
        n = 2
        A = reshape((/1.0_8, 2.0_8, 3.0_8, 4.0_8/), (/2,2/))
        x = (/1.0_8, 1.0_8/)
        y = (/0.0_8, 0.0_8/)
        alpha = 1.0_8
        beta = 0.0_8

        ! Call DGEMV: y = alpha*A*x + beta*y
        call dgemv('N', n, n, alpha, A, n, x, 1, beta, y, 1)

        ! Print the result
        print *, 'Result y: ', y
    end subroutine use_blas_example
end module blas_module
