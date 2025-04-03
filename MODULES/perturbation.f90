module perturbation
  use mpi
  implicit none
  private
  public :: magnetic_field

contains
  subroutine magnetic_field(nb, B_x, B_y, B_z, Hmag)
    implicit none
    ! Input
    integer, intent(in) :: nb
    real*8, intent(in) :: B_x, B_y, B_z
    ! Output
    complex*16, intent(out) :: Hmag(nb,nb)
    ! Local
    integer :: i
    complex*16 sigx(2, 2), sigy(2, 2), sigz(2, 2)
    complex*16 :: Hm(2,2) 

    ! Define Pauli Matrices
    data sigx /(0d0,0d0),(1d0,0d0),(1d0, 0d0),( 0d0, 0d0)/
    data sigy /(0d0,0d0),(0d0,1d0),(0d0,-1d0),( 0d0, 0d0)/
    data sigz /(1d0,0d0),(0d0,0d0),(0d0, 0d0),(-1d0, 0d0)/

    ! allocate(Hm(2,2),Hmag(nb,nb))
    Hm = B_x*sigx + B_y*sigy + B_z*sigz
  !------ Turn Hm into an 18x18 to match Hk      
    Hmag = (0d0, 0d0)
    do i=1, nb/2
          Hmag(i,i)=Hm(1,1)
          Hmag(i,i+nb/2)=Hm(1,2)
          Hmag(i+nb/2,i)=Hm(2,1)
          Hmag(i+nb/2,i+nb/2)=Hm(2,2)
    end do
  end subroutine magnetic_field
end module perturbation