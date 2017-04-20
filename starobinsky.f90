! Starobinsky's disappearing cosmological constant model
! f(R)/R0 = x + lambda*[1/(1+x^2)^n - 1] + mu*x^2, with x = R/R0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module starobinsky; implicit none

real, parameter :: R0 = 1.0e-2

integer, parameter :: n = 2
real, parameter :: lambda = 2.0, mu = 1.0e-6
real, parameter :: upsilon = sqrt((n*lambda*mu**(2*n+1))**(1.0/(n+1)) - mu**2)

real, parameter :: phi0 = -1.0e-4

contains

! return dimensionless Ricci curvature r = R/R0 corresponding to phi
elemental function curvature(phi)
	integer i; real curvature, x, p, q, s, phi; intent(in) phi
	
	! initial guess
	if (mu > 0.0 .and. phi > 0.0) then
		x = ((phi/2.0)**(n+1) + upsilon**(n+1))**(1.0/(n+1))/mu
	else
		x = sqrt((-(2.0*lambda*n)/(phi-2.0*upsilon))**(1.0/(n+0.5)) - 1.0)
	end if
	
	! Newton's method
	do i = 1,16
		p = x*x + 1.0; q = (2*n+1)*x*x - 1.0; s = p**(n+2)
		x = x + ((phi/2.0 - mu*x)*s + (n*lambda)*p*x)/(mu*s + (n*lambda)*q)
	end do
	
	curvature = x
end function curvature

! derivative of scalar field potential
elemental function DV(phi)
        real DV, x, p, q, phi; intent(in) :: phi
        
        x = curvature(phi); p = x*x + 1.0; q = (n+1)*x*x + 1.0
        DV = (R0/3.0) * (x + 2.0*lambda*(q/p**(n+1) - 1.0))
end function DV

! second derivative of scalar field potential
elemental function DDV(phi)
        real DDV, x, p, q, s, phi; intent(in) :: phi
        
        x = curvature(phi); p = x*x + 1.0; q = (2*n+1)*x*x - 1.0; s = p**(n+2)
        DDV = (2.0/3.0*R0) * (s/4.0 - n*(n+1)*lambda*x**3)/(mu*s + (n*lambda)*q)
end function DDV

end module
