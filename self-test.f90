! code testing implementation of curvaton-like models for internal consistency

program test; use starobinsky, only : curvaton, curvature, DV, DDV; implicit none

integer, parameter :: n = 1000 ! number of log-spaced curvature values to test

! curvature and curvaton values, potential and its derivatives for integration
real x(n), phi(n), error(n), y(3)

integer i

! initialize test curvaton grid
forall (i=1:n) x(i) = exp(log(1000.0)*(i-1)/(n-1))
phi = curvaton(x); error = curvature(phi)-x

! initialize potential integrator
y = [phi(1), 0.0, DV(phi(1))]

write (*,('(A)')) "# OUTPUT: x, phi, integral(V'), integral(V''), V'-integral(V''), curvature(phi)-x;"

do i = 2,n
	call gl10(y, phi(i)-phi(i-1))
	write (*,'(1F24.16,5G24.16)') x(i), y, DV(phi(i))-y(3), error(i)
end do

contains

! evaluate potential derivatives
subroutine evalf(y, dy)
        real, dimension(3) :: y, dy
        
        ! unmangle phase space vector contents into human-readable form
        associate (phi => y(1))
                dy = [1.0, DV(phi), DDV(phi)]
        end associate
end subroutine evalf

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 3
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 10th order Gauss-Legendre method
        real, parameter :: a(s,s) = reshape([ &
                  0.5923172126404727187856601017997934066Q-1, -1.9570364359076037492643214050884060018Q-2, &
                  1.1254400818642955552716244215090748773Q-2, -0.5593793660812184876817721964475928216Q-2, &
                  1.5881129678659985393652424705934162371Q-3,  1.2815100567004528349616684832951382219Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -2.4592114619642200389318251686004016630Q-2, &
                  1.0318280670683357408953945056355839486Q-2, -2.7689943987696030442826307588795957613Q-3, &
                  1.1377628800422460252874127381536557686Q-1,  2.6000465168064151859240589518757397939Q-1, &
                  1.4222222222222222222222222222222222222Q-1, -2.0690316430958284571760137769754882933Q-2, &
                  4.6871545238699412283907465445931044619Q-3,  1.2123243692686414680141465111883827708Q-1, &
                  2.2899605457899987661169181236146325697Q-1,  3.0903655906408664483376269613044846112Q-1, &
                  1.1965716762484161701032287870890954823Q-1, -0.9687563141950739739034827969555140871Q-2, &
                  1.1687532956022854521776677788936526508Q-1,  2.4490812891049541889746347938229502468Q-1, &
                  2.7319004362580148889172820022935369566Q-1,  2.5888469960875927151328897146870315648Q-1, &
                  0.5923172126404727187856601017997934066Q-1], [s,s])
        real, parameter ::   b(s) = [ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1]
        
        ! iterate trial steps
        g = 0.0; do k = 1,16
                g = matmul(g,a)
                do i = 1,s
                        call evalf(y + g(:,i)*dt, g(:,i))
                end do
        end do
        
        ! update the solution
        y = y + matmul(g,b)*dt
end subroutine gl10

end program
