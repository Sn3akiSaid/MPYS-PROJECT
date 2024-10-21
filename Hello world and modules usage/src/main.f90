program main
    use vectors
    use blas_module
    use lapack_module
    implicit none
    
    type(vector) :: v
    real :: l

    v = vector(1.0, 2.0, 3.0)
    l = v%length()

    print *,'Hello world'
    print *,'Length', l
    call use_blas_example()
    call use_lapack_example()

end program main