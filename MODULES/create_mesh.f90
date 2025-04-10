module create_mesh
  use mpi
  implicit none

  private

  public :: Lattice3D, Lattice2D

contains

  subroutine Lattice3D(np, dim, kbox_x, kbox_y, kbox_z, mesh, bvec)
    implicit none
    integer, intent(in) :: np, dim
    real*8, intent(in) :: kbox_x, kbox_y, kbox_z, bvec(3,3)
    real*8, intent(out) :: mesh(3, (2*np+1)**dim)

    integer j, i, n, k
    real*8 delkx, delky, delkz

    delkx = kbox_x/(2*np+1)
    delky = kbox_y/(2*np+1)
    delkz = kbox_z/(2*np+1)

    j=0
    do i = -np, np ! -np,np
        do n = -np, np
            do k = -np, np
                j = j+1  ! Direct index calculation
                ! Calculate coordinates directly
                mesh(1, j) = i*delkx
                mesh(2, j) = n*delky
                mesh(3, j) = k*delkz
            end do
        end do
    end do
    mesh(3,:)=mesh(3,:)+0.5d0*bvec(3,3)
  end subroutine Lattice3D

  subroutine Lattice2D(np, dim, kbox_x, kbox_y, mesh, bvec)
    implicit none
    integer, intent(in) :: np, dim
    real*8, intent(in) :: kbox_x, kbox_y, bvec(3,3)
    real*8, intent(out) :: mesh(3, (2*np+1)**dim)

    integer j, i, n
    real*8 delkx, delky

    delkx = kbox_x/(2*np+1)
    delky = kbox_y/(2*np+1)

    j=0
    do i = -np, np ! -np,np
        do n = -np, np
                j = j+1  ! Direct index calculation
                ! Calculate coordinates directly
                mesh(1, j) = i*delkx
                mesh(2, j) = n*delky
                mesh(3, j) = 0.5d0*bvec(3,3)
        end do
    end do

  end subroutine

end module create_mesh