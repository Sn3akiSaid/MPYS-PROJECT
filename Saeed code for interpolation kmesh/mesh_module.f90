module mesh_module
    implicit none
    integer, parameter :: n_points = 3
    real(8), dimension(3, n_points) :: kpath
    real(8), dimension(100, 100, 3) :: mesh  ! Adjust dimensions as needed

contains

    subroutine create_mesh()
        integer :: i, j

        ! Initialize kpath with the given data
        kpath(:, 1) = [-0.1d0, 0.0d0, 0.5d0]  ! L
        kpath(:, 2) = [0.0d0, 0.0d0, 0.5d0]   ! A
        kpath(:, 3) = [0.1d0, 0.0d0, 0.5d0]  ! L

        ! Create the mesh (example: simple interpolation between points)
        do i = 1, n_points - 1
            do j = 1, 100  ! Assuming 100 points in the interpolation
                mesh(j, i, 1) = kpath(1, i) + real(j - 1) * (kpath(1, i + 1) - kpath(1, i)) / 99.0d0
                mesh(j, i, 2) = kpath(2, i) + real(j - 1) * (kpath(2, i + 1) - kpath(2, i)) / 99.0d0
                mesh(j, i, 3) = kpath(3, i) + real(j - 1) * (kpath(3, i + 1) - kpath(3, i)) / 99.0d0
            end do
        end do

        ! Fill the last column of the mesh
        do j = 1, 100
            mesh(j, n_points, 1) = kpath(1, n_points)
            mesh(j, n_points, 2) = kpath(2, n_points)
            mesh(j, n_points, 3) = kpath(3, n_points)
        end do
        print*, mesh
    end subroutine create_mesh
end module mesh_module

