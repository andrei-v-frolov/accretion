! Hu-Sawicki's disappearing cosmological constant model
! f(R)/R0 = x - alpha*x^n/(beta*x^n + 1) + mu*x^2, with x = R/R0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module husawicki; implicit none

real, parameter :: R0 = 1.0e-2

integer, parameter :: n = 2
real, parameter :: alpha = 2.0, beta = 1.0, mu = 1.0e-6
real, parameter :: upsilon = mu**(1-1.0/n) * (mu**(2.0/(n+2))*(n*alpha/beta**2/2.0)**(n/(n+2.0))  - mu/beta)**(1.0/n)

real, parameter :: phi0 = -1.0e-4

contains

! return dimensionless Ricci curvature r = R/R0 corresponding to phi
elemental function curvature(phi)
	integer i; real curvature, x, p, q, s, phi; intent(in) phi
	
	! initial guess
	if (mu > 0.0 .and. phi > 0.0) then
		x = ((phi/2.0)**(n+1) + upsilon**(n+1))**(1.0/(n+1))/mu
	else
		x = ((-(alpha*n/beta**2)/(phi-2.0*upsilon))**(n/(n+1.0)) - 1.0/beta)**(1.0/n)
	end if
	
	! Newton's method
	do i = 1,16
		p = beta*x**n + 1.0; q = x**(n-2); s = p*p
		x = x + ((phi - 2.0*mu*x)*s + (n*alpha)*q*x)/(2.0*mu*s + (n*alpha)*q*(n+1 - 2*n/p))
	end do
	
	curvature = x
end function curvature

! derivative of scalar field potential
elemental function DV(phi)
        real DV, x, p, phi; intent(in) :: phi
        
        x = curvature(phi); p = beta*x**n + 1.0
        DV = (R0/3.0) * (x + alpha*(n/p-2.0) * x**n/p)
end function DV

! second derivative of scalar field potential
elemental function DDV(phi)
        real DDV, x, p, q, s, phi; intent(in) :: phi
        
        x = curvature(phi); p = beta*x**n + 1.0; q = x**(n-2); s = p*p
        DDV = (R0/3.0) * (s - (n*alpha)*q*x*(n+2-2*n/p))/(2.0*mu*s + (n*alpha)*q*(n+2-2*n/p))
end function DDV

end module
