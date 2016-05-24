! $Id$
! 
! tutorial code solving 1D wave equation Box[phi] = m^2 phi in Schwarzschield metric
!   ds^2 = -g(r) dt^2 + dr^2/g(r) + r^2 dOmega^2, g(r) = 1-2M/r
! in tortoise coordinates dx = dr/g(r) where d'Alembert operator becomes
!   g(r) Box[phi] = [-(d/dt)^2 + (d/dx)^2 + 2 g(r)/r (d/dx)] phi(x,t)
! which is a 1D wave equation with radial damping profile leading to back-scatter
! 
! as we are looking for smooth solutions, the method of choice to calculate
! spatial derivatives is pseudo-spectral, using Chebyshev basis on compactified
! coordinate y = x/sqrt(1+x^2), or inversely x = y/sqrt(1-y^2)
! 
! time integration is done using Gauss-Legendre method, which is A-stable
! and symplectic for Hamiltonian problems, as we have here
! 
! spherical collapse in scalar-tensor theories usually reduces to this type
! of PDE's given a proper field and coordinate choice, e.g. hep-th/0409117
! 
! Andrei Frolov <frolov@sfu.ca>; provided as is, modify or redistribute
! please cite arXiv:XXXX.XXXXX if you use this code for your research


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! basic how to:
! 
! compile with: ifort -O3 -ipo -r8 {-parallel|-openmp} -L /opt/intel/mkl/lib -lmkl_rt spectral.f90 <model>.f90
! or with GCC: gfortran -O3 -fdefault-real-8 {-fopenmp} -llapack spectral.f90 <model>.f90
! 
! run as: ./a.out > DATA; plot output using gnuplot: splot 'DATA' u 2:1:4 w l
! basic stop-frame animations can be done with gnuplot, for example using:
! set yrange [-1:1]; frame=0; plot 'DATA' index frame u 2:4 w l; load 'animate.gpl'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program wave; use starobinsky, only : DV; implicit none

! solver control parameters
integer, parameter :: nn = 2**9                 ! number of nodes to sample on (i.e. spectral order)
integer, parameter :: tt = 2**11                ! total number of time steps to take (i.e. runtime)
real,    parameter :: dt = 0.005                ! time step size (simulated timespan is tt*dt)

! collocation grid, metric functions, and spectral Laplacian operator
! phase space state vector v packs phi (1:nn) and dot(phi) (nn+1:nn+nn)
real theta(nn), x(nn), r(nn), g(nn), F(nn), L(nn,nn), v(2*nn)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer i

! initialize grid and Laplacian oprator
call initg(); call initl()

! sanity check on the time step value selected
if (x(nn/2+1)-x(nn/2) < dt) pause "Time step violates Courant condition, do you really want to run?"

! initial field and velocity profiles
v(1:nn) = -1.0e-4; v(nn+1:2*nn) = 0.0

! force term corresponding to matter distributuion truncated at 6M
F = DV(v(nn))*(tanh(5.0*(r-3.0)) + 1.0)/2.0; !call static(F, 100.0, v(1:nn))

! output initial conditions
call dump(0.0, v)

! main time evolution loop
do i = 1,tt
        call gl10(v, dt)
        call dump(i*dt, v)
end do

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! areal radius r from tortoise coordinate x (in units of 2M)
elemental function radius(x)
        real radius, x, q, r; intent(in) :: x; integer i
        
        ! closer to horizon, radius is 2M to all significant digits
        radius = 1.0; if (x < -(digits(x)-5)*log(2.0)) return
        
        ! refine q instead of r = 1.0 + log(1.0+exp(q)); initial guess
        q = x - 1.0
        
        ! Newton's iteration for x = r + log(r-1.0) as a function of q
        do i = 1,16
                r = 1.0 + log(1.0+exp(q))
                q = q - (r + log(r-1.0) - x) * (1.0 - 1.0/r) * (1.0 + exp(-q))
        end do
        
        ! find converged radius value
        radius = 1.0 + log(1.0+exp(q))
end function radius

! initialize Gauss-Lobatto collocation grid
subroutine initg()
        integer i; real, parameter :: pi = 3.1415926535897932384626433832795028841971694Q0
        
        ! Chebyshev collocation grid (extrema and roots-of varieties)
        !forall (i=1:nn) theta(i) = (nn-i)*pi/(nn-1) ! includes interval ends
        forall (i=1:nn) theta(i) = (nn-i+0.5)*pi/nn ! excludes interval ends
        
        x = cos(theta)/sin(theta); r = radius(x); g = 1.0 - 1.0/r
end subroutine initg

! evaluate rational Chebyshev basis on collocation grid theta
subroutine evalb(n, theta, Tn, Dn)
        integer n; real, dimension(nn) :: theta, Tn, Tnx, Tnxx, Dn
        intent(in) n, theta; intent(out) Tn, Dn; !optional Dn
        
        ! Chebyshev basis and its derivatives
        Tn = cos(n*theta)
        Tnx = n * sin(n*theta) * sin(theta)**2
        Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3
        
        ! Laplacian operator for Schwarzschild metric
        Dn = Tnxx + 2.0*g/r * Tnx
end subroutine evalb

! initialize Laplacian operator matrix
subroutine initl()
        integer i, pivot(nn), status; real A(nn,nn), B(nn,nn)
        
        ! evaluate basis and differential operator values on a grid
        do i = 1,nn; call evalb(i-1, theta, A(i,:), B(i,:)); end do
        
        ! find linear differentiation matrix
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, nn, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, nn, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! to evaluate Laplacian of function f, simply do matmul(L,f)
        L = transpose(B)
end subroutine initl

! evaluate equations of motion
subroutine evalf(v, dvdt)
        real v(2*nn), dvdt(2*nn)
        
        ! unmangle phase space state vector contents into human-readable form
        associate (phi => v(1:nn), pi => v(nn+1:2*nn), dphi => dvdt(1:nn), dpi => dvdt(nn+1:2*nn))
                dphi = pi; dpi = matmul(L,phi) - g*(DV(phi) - F)
        end associate
end subroutine evalf

! find static solution for a massive scalar field
subroutine static(F, m2, phi)
        real m2, F(nn), phi(nn), A(nn,nn), B(nn,1)
        integer i, pivot(nn), status
        
        ! set up Laplace equation for massive field
        A = L; do i = 1,nn
                A(i,i) = A(i,i) - m2*g(i)
                B(i,1) = -g(i)*F(i)
        end do
        
        ! find static solution by direct inversion
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, 1, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, 1, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! return static solution
        phi = B(:,1)
end subroutine static

! dump simulation data in plain text format
subroutine dump(t, v)
        real t, v(2*nn); integer i
        
        do i = 1,nn
                write (*,'(3F24.16,3G24.16)') t, x(i), r(i), v(i), v(nn+i)
        end do
        
        ! separate timesteps by empty lines (for gnuplot's benefit)
        write (*,'(A)') "", ""
end subroutine

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 2*nn
        real y(n), g(n,s), dt; integer i, k
        
        ! Butcher tableau for 8th order Gauss-Legendre method
        real, parameter :: a(s,s) = (/ &
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
                  0.5923172126404727187856601017997934066Q-1 /)
        real, parameter ::   b(s) = (/ &
                  1.1846344252809454375713202035995868132Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  2.8444444444444444444444444444444444444Q-1,  2.3931433524968323402064575741781909646Q-1, &
                  1.1846344252809454375713202035995868132Q-1 /)
        
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

end