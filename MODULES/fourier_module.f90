module fourier_module
  use mpi
  implicit none
  
  ! Make everything private by default
  private
  
  ! Public interfaces
  public :: fourier_transform_optimized, fourier_transform_general
  
  ! Module-level constants
  real*8, parameter, private :: twopi = 4.0d0*atan(1.0d0)*2.0d0
  
contains
  !--------------------------------------------------------------------
  ! Subroutine for the optimized Fourier transform (18x18 case)
  subroutine fourier_transform_optimized(np, dim, nr, nb, ndeg, mesh, rvec, &
                                        Hamr_trivial, Hamr_topological, Hmag, alpha, &
                                        enep, ene, work, lwork, rwork, rank, ierr)
      implicit none
      ! Inputs:
      integer, intent(in) :: np, dim, nr, nb, lwork, rank
      real*8, intent(in) :: alpha
      integer, intent(in) :: ndeg(nr)
      real*8, intent(in) :: mesh(3, (np+1)**dim)
      real*8, intent(in) :: rvec(3, nr)
      complex*16, intent(in) :: Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr)
      complex*16, intent(in) :: Hmag(nb, nb)
      ! Outputs:
      real*8, intent(out) :: enep(nb, (np+1)**dim), ene(nb, (np+1)**dim)
      ! Workspace:
      complex*16, intent(inout) :: work(lwork)
      real*8, intent(inout) :: rwork(*)
      integer, intent(out) :: ierr
      ! Local variables:
      integer :: k, j, kspace, info
      real*8 :: phase
      complex*16 :: phase_factor
      complex*16, allocatable :: phases(:,:)
      complex*16, allocatable :: Hk_trivial(:,:), Hk_topological(:,:), Hk(:,:), H(:,:)
      
      kspace = (np+1)**dim
      allocate(phases(nr, kspace))
      allocate(Hk_trivial(nb, nb), Hk_topological(nb, nb), Hk(nb, nb), H(nb, nb))
      
      ! Compute phase factors: for each partition (k) and contribution (j)
      do k = 1, kspace
          do j = 1, nr
              phase = dot_product(mesh(:, k), rvec(:, j))
              phases(j, k) = dcmplx(cos(phase), -sin(phase)) / dble(ndeg(j))
          end do
      end do
      
      ! Loop over partitions: accumulate the Fourier sums and compute eigenvalues
      do k = 1, kspace
          ! Initialize accumulation arrays to zero
          Hk_trivial = cmplx(0d0, 0d0)
          Hk_topological = cmplx(0d0, 0d0)
          
          do j = 1, nr
              phase_factor = phases(j, k)
              Hk_trivial = Hk_trivial + Hamr_trivial(:, :, j) * phase_factor
              Hk_topological = Hk_topological + Hamr_topological(:, :, j) * phase_factor
          end do
          
          ! Interpolate and add perturbation:
          Hk = Hk_trivial * (1.0d0 - alpha) + Hk_topological * alpha
          H = Hk + Hmag
          
          ! Compute eigenvalues/eigenvectors using LAPACK's ZHEEV
          call zheev('V', 'U', nb, H, nb, enep(:, k), work, lwork, rwork, info)
          call zheev('V', 'U', nb, Hk, nb, ene(:, k), work, lwork, rwork, info)
          
          if (info /= 0) then
              if (rank == 0) then
                  write(*,*) "ZHEEV failed with info =", info
              end if
              call MPI_ABORT(MPI_COMM_WORLD, info, ierr)
          end if
      end do
      
      deallocate(phases, Hk_trivial, Hk_topological, Hk, H)
  end subroutine fourier_transform_optimized
  
  !--------------------------------------------------------------------
  ! Subroutine for the general Fourier transform (4x4 case)
  subroutine fourier_transform_general(np, dim, nr, nb, ndeg, mesh, bvec, avec, &
                                      rvec_trivial, rvec_topological, Hamr_trivial, Hamr_topological, &
                                      Hmag, alpha, enep, ene, work, lwork, rwork, rank, ierr)
      implicit none
      ! Inputs:
      integer, intent(in) :: np, dim, nr, nb, lwork, rank
      real*8, intent(in) :: alpha
      integer, intent(in) :: ndeg(nr)
      real*8, intent(in) :: mesh(3, (np+1)**dim)
      real*8, intent(in) :: bvec(3,3), avec(3,3)
      real*8, intent(in) :: rvec_trivial(3, nr), rvec_topological(3, nr)
      complex*16, intent(in) :: Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr)
      complex*16, intent(in) :: Hmag(nb, nb)
      ! Outputs:
      real*8, intent(out) :: enep(nb, (np+1)**dim), ene(nb, (np+1)**dim)
      ! Workspace:
      complex*16, intent(inout) :: work(lwork)
      real*8, intent(inout) :: rwork(*)
      integer, intent(out) :: ierr
      ! Local variables:
      integer :: k, j, kspace, i, info
      real*8 :: phase_trivial, phase_topological
      complex*16 :: phase_factor_trivial, phase_factor_topological
      complex*16, allocatable :: Hk_trivial(:,:), Hk_topological(:,:), Hk(:,:), H(:,:)
      
      kspace = (np+1)**dim
      allocate(Hk_trivial(nb, nb), Hk_topological(nb, nb), Hk(nb, nb), H(nb, nb))
      
      ! Loop over partitions
      do k = 1, kspace
          ! Initialize accumulators
          Hk_trivial = cmplx(0d0, 0d0)
          Hk_topological = cmplx(0d0, 0d0)
          
          ! Compute contributions from each site
          do j = 1, nr
              ! Calculate phases
              phase_trivial = 0.0d0
              phase_topological = 0.0d0
              
              do i = 1, size(bvec, 2)
                  phase_trivial = phase_trivial + dot_product(mesh(:, k)*bvec(:, i), &
                                                         rvec_trivial(:, j)*avec(:, i))
                  phase_topological = phase_topological + dot_product(mesh(:, k)*bvec(:, i), &
                                                             rvec_topological(:, j)*avec(:, i))
              end do
              
              ! Convert to complex phase factors
              phase_factor_trivial = dcmplx(cos(phase_trivial), -sin(phase_trivial)) / dble(ndeg(j))
              phase_factor_topological = dcmplx(cos(phase_topological), -sin(phase_topological)) / dble(ndeg(j))
              
              ! Add contributions to Hamiltonian
              Hk_trivial = Hk_trivial + Hamr_trivial(:, :, j) * phase_factor_trivial
              Hk_topological = Hk_topological + Hamr_topological(:, :, j) * phase_factor_topological
          end do
          
          ! Interpolate and add perturbation
          Hk = Hk_trivial * (1.0d0 - alpha) + Hk_topological * alpha
          H = Hk + Hmag
          
          ! Compute eigenvalues/eigenvectors
          call zheev('V', 'U', nb, H, nb, enep(:, k), work, lwork, rwork, info)
          call zheev('V', 'U', nb, Hk, nb, ene(:, k), work, lwork, rwork, info)
          
          if (info /= 0) then
              if (rank == 0) then
                  write(*,*) "ZHEEV failed with info =", info
              end if
              call MPI_ABORT(MPI_COMM_WORLD, info, ierr)
          end if
      end do
      
      deallocate(Hk_trivial, Hk_topological, Hk, H)
  end subroutine fourier_transform_general
  
end module fourier_module