! massive scalar field
! V(phi) = m^2/2 phi^2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module massive; implicit none

real, parameter :: m2 =   0.0
real, parameter :: phi0 = 1.0

contains

! derivative of scalar field potential
elemental function DV(phi)
        real DV, phi; intent(in) :: phi
        
        DV = m2 * phi
end function DV

! second derivative of scalar field potential
elemental function DDV(phi)
        real DDV, phi; intent(in) :: phi
        
        DDV = m2
end function DDV

end module
