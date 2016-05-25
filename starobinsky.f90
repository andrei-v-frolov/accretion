! Starobinsky's disappearing cosmological constant model
! f(R)/R0 = x + lambda*[1/(1+x^2)^n - 1], with x = R/R0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module starobinsky; implicit none

real, parameter :: R0 = 1.0e-2
real, parameter :: lambda = 2.0
integer, parameter :: n = 2

real, parameter :: phi0 = -1.0e-4

contains

! return dimensionless Ricci curvature r = R/R0 corresponding to phi
elemental function r(phi)
	integer i; real r, x, p, q, phi; intent(in) phi
	
	! initial guess
	x = sqrt((-(2.0*lambda*n)/phi)**(1.0/(n+0.5)) - 1.0)
	
	! Newton's method
	do i = 1,16
		p = 1.0 + x*x; q = 1.0 - (2*n+1)*x*x
		x = x - (phi * p**(n+1)/(2.0*lambda*n) + x) * p/q
	end do
	
	r = x
end function r

! derivative of scalar field potential
elemental function DV(phi)
        real DV, x, p, q, phi; intent(in) :: phi
        
        x = r(phi); p = 1.0 + x*x; q = 1.0 + (n+1)*x*x
        DV = R0/3.0 * (x + 2.0*lambda*(q/p**(n+1) - 1.0))
end function DV

! second derivative of scalar field potential
elemental function DDV(phi)
        real DDV, x, p, q, phi; intent(in) :: phi
        real, parameter :: w = (2.0*R0)/(3.0*lambda*n)
        
        x = r(phi); p = 1.0 + x*x; q = 1.0 - (2*n+1)*x*x
        DDV = w * (n*(n+1)*lambda*x**3 - p**(n+2)/4.0)/q
end function DDV

end module
