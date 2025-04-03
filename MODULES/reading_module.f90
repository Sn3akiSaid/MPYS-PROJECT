module reading_module
  use mpi
  implicit none
  
  ! Make everything private by default
  private
  
  ! Public interfaces
  public :: read_lattice, read_optimized_hamiltonians, read_general_hamiltonians, read_header
  
contains
  !--------------------------------------------------------------------!
  subroutine read_lattice(nnkp, rank, avec, bvec, ierr)
    implicit none
    ! Inputs
    character(len=*), intent(in) :: nnkp
    integer, intent(in) :: rank
    ! Outputs
    real*8, intent(out) :: avec(3,3), bvec(3,3)
    integer, intent(out) :: ierr
    ! Local variables
    character(len=512) :: line
    
    ! Initialize error code
    ierr = 0
!---------------  Read the vectors
        open(98,file=trim(adjustl(nnkp)), status='old')
111     read(98,'(a)')line
        if(trim(adjustl(line)).ne."begin real_lattice") goto 111
        read(98,*)avec
        read(98,'(a)')line
        read(98,'(a)')line
        read(98,'(a)')line
        read(98,*)bvec
        close(98)

end subroutine read_lattice
  ! Subroutine to read optimized Hamiltonians (18x18 case)
  subroutine read_optimized_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr, &
                                        Hamr_trivial, Hamr_topological, rvec, ndeg, avec)
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
      
      ! Open the trivial and topological files
      open(newunit=unit_trivial, file=trim(hamil_file_trivial), status='old', action='read')
      open(newunit=unit_topological, file=trim(hamil_file_topological), status='old', action='read')
      read(unit_trivial, *)        ! Skip line 1
      read(unit_trivial, *)        ! Skip line 2
      read(unit_trivial, *)        ! Skip line 3
      read(unit_trivial, *) ndeg
      ! skip header information and dimensions
      do i = 1, 80
          read(unit_topological, *)
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
          ! Here you combine rvecs with avec to build rvec
          rvec(:,k) = rvecs(1)*avec(:,1) + rvecs(2)*avec(:,2) + rvecs(3)*avec(:,3)
      end do
      
      close(unit_trivial)
      close(unit_topological)
  end subroutine read_optimized_hamiltonians
  
  !--------------------------------------------------------------------
  ! Subroutine to read general Hamiltonians (4x4 case)
  subroutine read_general_hamiltonians(hamil_file_trivial, hamil_file_topological, nb, nr_trivial, nr_topological, &
                                      Hamr_trivial, Hamr_topological, rvec_trivial, rvec_topological, ndeg_trivial, ndeg_topological,  avec)
      implicit none
      ! Inputs:
      character(len=*), intent(in) :: hamil_file_trivial, hamil_file_topological
      integer, intent(in) :: nb, nr_trivial, nr_topological
      real*8, intent(in) :: avec(3,3)
      ! Outputs:
      complex*16, intent(out) :: Hamr_trivial(nb, nb, nr_trivial), Hamr_topological(nb, nb, nr_topological)
      real*8, intent(out) :: rvec_trivial(3, nr_trivial), rvec_topological(3, nr_topological)
      integer, intent(out) :: ndeg_trivial(nr_trivial),ndeg_topological(nr_topological)
      ! Local variables:
      integer :: i, j, k, i1, i2
      real*8 :: a, b, a1, b1
      real*8 :: rvecs_trivial(3), rvecs_topological(3)
      integer :: unit_trivial, unit_topological
      
      ! Open files
      open(newunit=unit_trivial, file=trim(hamil_file_trivial), status='old', action='read')
      open(newunit=unit_topological, file=trim(hamil_file_topological), status='old', action='read')
      read(unit_topological, *)    ! Skip line 1 ("written on 5Nov2024 at 13:56:48")
      read(unit_topological, *)    ! Skip line 2 (nb)
      read(unit_topological, *)    ! Skip line 3 (nr)
      read(unit_topological, *) ndeg_topological
      ! Similarly for the trivial file:
      read(unit_trivial, *)        ! Skip line 1
      read(unit_trivial, *)        ! Skip line 2
      read(unit_trivial, *)        ! Skip line 3
      read(unit_trivial, *) ndeg_trivial
      
      do k = 1, nr_trivial
          do i = 1, nb
              do j = 1, nb
                  read(unit_trivial, *) rvecs_trivial(1),rvecs_trivial(2),rvecs_trivial(3), i1, i2, a, b
                  Hamr_trivial(i1, i2, k) = dcmplx(a, b)
              end do
          end do
          
          ! Store the vectors
          rvec_trivial(:,k) = rvecs_trivial(1)*avec(:,1) + rvecs_trivial(2)*avec(:,2) + rvecs_trivial(3)*avec(:,3)
      end do

      do k = 1, nr_topological
        do i = 1, nb
            do j = 1, nb                
                read(unit_topological, *) rvecs_topological(1),rvecs_topological(2),rvecs_topological(3), i1, i2, a1, b1
                Hamr_topological(i1, i2, k) = dcmplx(a1, b1)
            end do
        end do
        rvec_topological(:,k) = rvecs_topological(1)*avec(:,1) + rvecs_topological(2)*avec(:,2) + rvecs_topological(3)*avec(:,3)
    end do
      
      close(unit_trivial)
      close(unit_topological)
  end subroutine read_general_hamiltonians
  
  subroutine read_header(hamil_file,nb,nr)
    implicit none
    character(len=*), intent(in) :: hamil_file
    integer, intent(out) :: nb,nr
    integer :: unit

    open(newunit=unit, file=trim(hamil_file), status='old', action='read')

    read(unit, *)
    read(unit, *) nb
    read(unit, *) nr
    close(unit)
  end subroutine read_header

end module reading_module
