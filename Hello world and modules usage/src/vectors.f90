module vectors
    implicit none
    
    type vector
        real x, y, z
        contains
        procedure, pass(self) :: length
    end type

    contains
    
    function length(self) result(l)
        class(vector) :: self
        real :: l
        l = sqrt(self%x**2 + self%y**2 + self%z**2)
    end function

end module vectors