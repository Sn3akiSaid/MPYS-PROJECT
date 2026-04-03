module create_mesh
  use mpi
  implicit none

  private

  public :: Lattice3D, berryLattice3D, SpecificLattice, Lattice2D

contains

  subroutine Lattice3D(np, dim, kbox_x, kbox_y, kbox_z, delkx, delky, delkz, mesh, bvec, j)
    implicit none
    integer, intent(in) :: np, dim
    real*8, intent(in) :: kbox_x, kbox_y, kbox_z, bvec(3,3)
    real*8, intent(out) :: delkx, delky, delkz, &
                           mesh(3, (2*np+1)**dim)

    integer i, n, k
    integer, intent(out) :: j

    delkx = kbox_x/(2*np+1)
    delky = kbox_y/(2*np+1)
    delkz = kbox_z/(2*np+1)

    j=0
    do i = -np, np ! -np,np
        do n = -np, np
            do k = -np, np
                j = j+1  ! Direct index calculation
                ! Calculate coordinates directly
                mesh(1, j) = i*delkx + 0.049d0
                mesh(2, j) = n*delky + 0.030d0
                mesh(3, j) = k*delkz + 0.5d0*bvec(3,3)
            end do
        end do
    end do
    ! mesh(3,:)=mesh(3,:)
  end subroutine Lattice3D
	subroutine berryLattice3D(np, dim, kbox_x, kbox_y, kbox_z, delkx, delky, delkz, mesh, bvec, j)
    implicit none
    integer, intent(in) :: np, dim
    real*8, intent(in) :: kbox_x, kbox_y, kbox_z, bvec(3,3)
		real*8, intent(inout) :: mesh(3, 2*np+1, 2*np+1, 2*np+1)
    real*8, intent(out) :: delkx, delky, delkz
                           

    integer i, n, k, kx, ky, kz
    integer, intent(out) :: j

    delkx = kbox_x/(2*np)
    delky = kbox_y/(2*np)
    delkz = kbox_z/(2*np)

    j=0
    do kx = 1, 2*np+1 !
        do ky = 1, 2*np+1
            do kz = 1, 2*np+1
                j = j+1  ! Direct index calculation
                ! Calculate coordinates directly
                ! offset(:,1,3) = (/ ,,/) !+ve
                ! ,,
                mesh(1, kx, ky, kz) = (kx-1)*delkx-kbox_x/2!-0.01765908461242609!+0.01766767357940942 ! + 0.046d0!
                mesh(2, kx, ky, kz) = (ky-1)*delky-kbox_y/2!-0.04663940372172213!+(-0.04650968222518360) ! + 0.03d0!
                mesh(3, kx, ky, kz) = (kz-1)*delkz-kbox_z/2+0.5d0*bvec(3,3)! + 0.43932749486813732!0.47482330199264761  !
            end do
        end do
    end do
  end subroutine berryLattice3D

  subroutine SpecificLattice(mesh, np, kx_min, kx_max, ky_min, ky_max, kz_min, kz_max, bvec, j)
    implicit none
    
    ! Parameters
    integer, intent(in) :: np           ! Number of points per dimension     
    real(8), intent(in) :: bvec(3,3)
    real(8), intent(in) :: kx_min, kx_max        ! X-bounds
    real(8), intent(in) :: ky_min, ky_max        ! Y-bounds
    real(8), intent(in) :: kz_min, kz_max        ! Z-bounds
    real(8), intent(out) :: mesh(3, (2*np+1)**3) ! Output mesh
    
    ! Local variables
    integer :: i, n, k
    integer, intent(out) :: j
    real(8) :: delkx, delky, delkz
    
    ! Calculate step sizes
    delkx = (kx_max - kx_min) / (2*np)
    delky = (ky_max - ky_min) / (2*np)
    delkz = (kz_max - kz_min) / (2*np)
    
    ! Create the mesh
    j = 0
    do i = -np, np 
        do n = -np, np
            do k = -np, np
                j = j + 1
                ! Calculate coordinates directly with the specified bounds
                mesh(1, j) = kx_min + (i + np) * delkx
                mesh(2, j) = ky_min + (n + np) * delky
                mesh(3, j) = kz_min + (k + np) * delkz
            end do
        end do
    end do
    
    ! Apply any additional transformations if needed
    ! mesh(3,:) = mesh(3,:) + 0.5d0*bvec(3,3)
  end subroutine SpecificLattice

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
                mesh(1, j) = i*delkx !+ 0.03d0
                mesh(2, j) = n*delky !+ 0.04d0
                mesh(3, j) = 0.5d0*bvec(3,3)
        end do
    end do

  end subroutine Lattice2D

end module create_mesh