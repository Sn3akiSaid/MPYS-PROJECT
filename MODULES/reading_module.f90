module reading_module
    use mpi
    implicit none

    private

    public :: read_optimized_hamiltonians, read_general_hamiltonians
  contains
    !--------------------------------------------------------------------
    ! Subroutine to read optimized Hamiltonians (18x18 case)
    subroutine read_optimized_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr, &
                                           Hamr_trivial, Hamr_topological, rvec, ndeg,avec)
      implicit none

      ! Inputs:
      character(len=*), intent(in) :: hamil_file_trivial, hamil_file_topological
      integer, intent(in) :: nb, nr
      real*8, intent(in) :: avec(3,3)

      ! Outputs:
      complex*16, intent(out) :: Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr)
      real*8, intent(out) :: rvec(3, nr)
      integer, intent(out) :: ndeg(nr)
      

      ! Local variables:
      integer :: i, j, k, i1, i2
      real*8 :: a, b, a1, b1, rvecs(3)
      integer :: unit_trivial, unit_topological

!-------------------------------------------------------------------!      
      
      ! Open the trivial and topological files
      open(newunit=unit_trivial, file=trim(hamil_file_trivial), status='old', action='read')
      open(newunit=unit_topological, file=trim(hamil_file_topological), status='old', action='read')
  
      ! Read header information and dimensions
      read(unit_trivial, *) ndeg
      do i = 1, 80
         read(unit_trivial, *)
      end do
  
      do k = 1, nr
         do i = 1, nb
           do j = 1, nb
             read(unit_trivial, *) rvecs(1), rvecs(2), rvecs(3), i1, i2, a, b
             Hamr_trivial(i1, i2, k) = dcmplx(a, b)
             read(unit_topological, *) rvecs(1), rvecs(2), rvecs(3), i1, i2, a1, b1
             Hamr_topological(i1, i2, k) = dcmplx(a1, b1)
           end do
         end do
         ! Here you might combine rvecs with another vector (e.g. 'avec') to build rvec
         rvec(:,k) = rvecs(1)*avec(:,1) + rvecs(2)*avec(:,2) + rvecs(3)*avec(:,3)
        end do
  
      close(unit_trivial)
      close(unit_topological)
    end subroutine read_optimized_hamiltonians
  
    !--------------------------------------------------------------------
    ! Subroutine to read general Hamiltonians (4x4 case)
    subroutine read_general_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr, &
                                         Hamr_trivial, Hamr_topological, rvec_trivial, rvec_topological, ndeg, avec)
      implicit none
      ! Inputs:
      character(len=*), intent(in) :: hamil_file_trivial, hamil_file_topological
      integer, intent(in) :: nb, nr
      real*8, intent(in) :: avec(3,3)
      ! Outputs:
      complex*16, intent(out) :: Hamr_trivial(nb, nb, nr), Hamr_topological(nb, nb, nr)
      real*8, intent(out) :: rvec_trivial(3, nr), rvec_topological(3, nr)
      integer, intent(out) :: ndeg(nr)
  
      ! Local variables:
      integer :: i, j, k, i1, i2
      real*8 :: a, b, a1, b1
      real*8 :: rvecs_trivial(3), rvecs_topological(3)
      integer :: unit_trivial, unit_topological
  
      ! Open files (assume the 4x4 files have different header structure)
      open(newunit=unit_trivial, file=trim(hamil_file_trivial), status='old', action='read')
      open(newunit=unit_topological, file=trim(hamil_file_topological), status='old', action='read')
  
      ! Skip header lines as needed for 4x4 reading (these numbers come from your commented code)
      read(unit_trivial, *) ndeg
      read(unit_topological, *)
      do i = 1, 79
         read(unit_topological, *)
      end do
      do i = 1, 80
         read(unit_trivial, *)
      end do
  
      do k = 1, nr
         do i = 1, nb
           do j = 1, nb
            read(unit_trivial, *) rvecs_trivial(1:3), i1, i2, a, b
            Hamr_trivial(i1, i2, k) = dcmplx(a, b)
            rvec_trivial(:,k) = rvecs_trivial
   
            read(unit_topological, *) rvecs_topological(1:3), i1, i2, a1, b1
            Hamr_topological(i1, i2, k) = dcmplx(a1, b1)
            rvec_topological(:,k) = rvecs_topological
           end do
         end do
      end do
  
      close(unit_trivial)
      close(unit_topological)
    end subroutine read_general_hamiltonians
  
  end module reading_module
  
