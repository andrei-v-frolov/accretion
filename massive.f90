! massive scalar field
! V(phi) = m^2/2 phi^2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module massive; implicit none

real, parameter :: m2 = 100.0


contains

! derivative of scalar field potential
elemental function DV(phi)
        real DV, phi; intent(in) :: phi
        
        DV = m2 * phi
end function DV

end module
