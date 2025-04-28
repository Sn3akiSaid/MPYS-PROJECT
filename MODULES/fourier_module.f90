module fourier_module
  
  use mpi
  implicit none
  
  ! Make everything private by default
  private
  
  ! Public interfaces
  public :: fourier_transform_optimized, fourier_transform_general, inner_ft_optimized, inner_ft_berry
  
  ! Module-level constants
  real*8, parameter, private :: twopi = 4.0d0*atan(1.0d0)*2.0d0
  
contains
  subroutine inner_ft_optimized(k, nr, nb, ndeg, mesh, rvec, &
                                Hamr_trivial, Hamr_topological, Hmag, alpha, &
                                Hk_trivial, Hk_topological, H, Hk, rank, ierr)
  implicit none

  integer, intent(in) :: k, nr, nb, rank
  integer, intent(in) :: ndeg(nr)
  real*8, intent(in) :: mesh(3, *), rvec(3, nr), alpha
  complex*16, intent(in) :: Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr)
  complex*16, intent(in) :: Hmag(nb, nb)
  ! Outputs
  complex*16, intent(out) :: Hk_trivial(nb, nb), Hk_topological(nb, nb), H(nb, nb), Hk(nb, nb)
  integer, intent(out) :: ierr
  ! Local variables
  integer :: j
  real*8 :: phase
  complex*16 :: phase_factor
  ! complex*16 :: Hk(nb, nb)

 ! Initialize accumulation arrays to zero
    Hk_trivial = (0d0, 0d0)
    Hk_topological = (0d0, 0d0)
    
    do j = 1, nr
        phase = dot_product(mesh(:, k), rvec(:, j))
        phase_factor = dcmplx(cos(phase), -sin(phase)) / float(ndeg(j))
      !   phase_factor = phases(j, k)!!! CONT FROM HERE
        Hk_trivial = Hk_trivial + Hamr_trivial(:, :, j) * phase_factor
        Hk_topological = Hk_topological + Hamr_topological(:, :, j) * phase_factor
    end do
    
    ! Interpolate and add perturbation:
    Hk = Hk_trivial * (1.0d0 - alpha) + Hk_topological * alpha
    H = Hk + Hmag
    
    ierr = 0
  end subroutine inner_ft_optimized
   subroutine inner_ft_berry(kx, ky, kz, nr, nb, ndeg, mesh, rvec, &
                            Hamr_trivial, Hamr_topological, Hmag, alpha, &
                            Hk_trivial, Hk_topological, H, Hk, rank, ierr)
  implicit none

  integer, intent(in) :: kx, ky, kz, nr, nb, rank
  integer, intent(in) :: ndeg(nr)
  real*8, intent(in) :: mesh(3, kx, ky, kz), rvec(3, nr), alpha
  complex*16, intent(in) :: Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr)
  complex*16, intent(in) :: Hmag(nb, nb)
  ! Outputs
  complex*16, intent(out) :: Hk_trivial(nb, nb), Hk_topological(nb, nb), H(nb, nb), Hk(nb, nb)
  integer, intent(out) :: ierr
  ! Local variables
  integer :: j
  real*8 :: phase
  complex*16 :: phase_factor
  ! complex*16 :: Hk(nb, nb)

 ! Initialize accumulation arrays to zero
    Hk=0
    Hk_trivial = (0d0, 0d0)
    Hk_topological = (0d0, 0d0)
    
    do j = 1, nr
        phase = dot_product(mesh(:, kx, ky, kz), rvec(:, j))
        phase_factor = dcmplx(cos(phase), -sin(phase)) / float(ndeg(j))
      !   phase_factor = phases(j, k)!!! CONT FROM HERE
        Hk_trivial = Hk_trivial + Hamr_trivial(:, :, j) * phase_factor
        Hk_topological = Hk_topological + Hamr_topological(:, :, j) * phase_factor
    end do
    
    ! Interpolate and add perturbation:
    Hk = Hk_trivial * (1.0d0 - alpha) + Hk_topological * alpha
    H = Hk + Hmag
    
    ierr = 0
  end subroutine inner_ft_berry
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
      integer :: k, kspace, info, inner_ierr
      ! real*8 :: phase
      ! complex*16 :: phase_factor
      ! complex*16, allocatable :: phases(:,:)
      complex*16, allocatable :: Hk_trivial(:,:), Hk_topological(:,:), Hk(:,:), H(:,:)
      
      kspace = (np+1)**dim
    !   allocate(phases(nr, kspace))
      allocate(Hk_trivial(nb, nb), Hk_topological(nb, nb), Hk(nb, nb), H(nb, nb))
      ! Loop over partitions: accumulate the Fourier sums and compute eigenvalues
      do k = 1, kspace

        call inner_ft_optimized(k, nr, nb, ndeg, mesh, rvec, &
                                Hamr_trivial, Hamr_topological, Hmag, alpha, &
                                Hk_trivial, Hk_topological, H, Hk, rank, inner_ierr)         
      
          if (inner_ierr /= 0) then
              ierr = inner_ierr
              deallocate(Hk_trivial, Hk_topological, H)
              return
          end if

      
      ! Compute eigenvalues/eigenvectors using LAPACK's ZHEEV
      !call zheev('V', 'U', nb, H, nb, enep(:, k), work, lwork, rwork, info)
      call zheev('V', 'U', nb, Hk, nb, ene(:, k), work, lwork, rwork, info)
      end do
      
      deallocate(Hk_trivial, Hk_topological, Hk, H)
  end subroutine fourier_transform_optimized
  
  !--------------------------------------------------------------------
  ! Subroutine for the general Fourier transform (4x4 case)
  subroutine fourier_transform_general(np, dim, nr_trivial, nr_topological, nb, ndeg_trivial, ndeg_topological, mesh,&
                                      rvec_trivial, rvec_topological, Hamr_trivial, Hamr_topological,&
                                      Hmag, alpha, enep, ene, work, lwork, rwork, rank, ierr)
      implicit none
      ! Inputs:
      integer, intent(in) :: np, dim, nr_trivial,nr_topological, nb, lwork, rank
      real*8, intent(in) :: alpha
      integer, intent(in) :: ndeg_topological(nr_topological),ndeg_trivial(nr_trivial)
      real*8, intent(in) :: mesh(3, (np+1)**dim)
      real*8, intent(in) :: rvec_trivial(3, nr_trivial), rvec_topological(3, nr_topological)
      complex*16, intent(in) :: Hamr_trivial(nb, nb, nr_trivial), Hamr_topological(nb, nb, nr_topological)
      complex*16, intent(in) :: Hmag(nb, nb)
      ! Outputs:
      real*8, intent(out) :: enep(nb, (np+1)**dim), ene(nb, (np+1)**dim)
      ! Workspace:
      complex*16, intent(inout) :: work(lwork)
      real*8, intent(inout) :: rwork(*)
      integer, intent(out) :: ierr
      ! Local variables:
      integer :: k, j, kspace, info
      real*8 :: phase_trivial, phase_topological
      complex*16 :: phase_factor_trivial, phase_factor_topological
    !   complex*16, allocatable :: phase_factor_trivial(:,:), phase_factor_topological(:,:)
      complex*16, allocatable :: Hk_trivial(:,:), Hk_topological(:,:), Hk(:,:), H(:,:)
      
      kspace = (np+1)**dim
    !   allocate(phase_factor_trivial(nr_trivial, kspace), phase_factor_topological(nr_topological, kspace))
      allocate(Hk_trivial(nb, nb), Hk_topological(nb, nb), Hk(nb, nb), H(nb, nb))
      
      ! Loop over partitions
      do k = 1, kspace
        ! Initialize accumulators
        Hk_trivial = (0d0, 0d0)
        Hk_topological = (0d0, 0d0)
              
        ! Compute contributions from each site
        ! Sum contributions from trivial sites:
        do j = 1, nr_trivial
            phase_trivial = dot_product(mesh(:, k), rvec_trivial(:, j))
            phase_factor_trivial = dcmplx(cos(phase_trivial), -sin(phase_trivial)) / float(ndeg_trivial(j))
            Hk_trivial = Hk_trivial + Hamr_trivial(:, :, j) * phase_factor_trivial!(j,k)
        end do
        ! Sum contributions from topological sites:
        do j = 1, nr_topological
            phase_topological = dot_product(mesh(:, k), rvec_topological(:, j))
            phase_factor_topological = dcmplx(cos(phase_topological), -sin(phase_topological)) / float(ndeg_topological(j))
            Hk_topological = Hk_topological + Hamr_topological(:, :, j) * phase_factor_topological!(j,k)
        end do
        ! Interpolate and add perturbation
        Hk = Hk_trivial * (1.0d0 - alpha) + Hk_topological * alpha
        H = Hk + Hmag
              
        ! Compute eigenvalues/eigenvectors
        call zheev('V', 'U', nb, H, nb, enep(:, k), work, lwork, rwork, info)
        !call zheev('V', 'U', nb, Hk, nb, ene(:, k), work, lwork, rwork, info)
              
        if (info /= 0) then
             if (rank == 0) then
                 write(*,*) "ZHEEV failed with info =", info
             end if
             call MPI_ABORT(MPI_COMM_WORLD, info, ierr)
        end if
    end do
    
      deallocate(Hk_trivial, Hk_topological, Hk, H)!, phase_factor_trivial, phase_factor_topological)
  end subroutine fourier_transform_general
  
end module fourier_module
