! tutorial code solving 1D wave equation Box[phi] = m^2 phi in Schwarzschild metric
!   ds^2 = -g(r) dt^2 + dr^2/g(r) + r^2 dOmega^2, g(r) = 1-2M/r
! in tortoise coordinates dx = dr/g(r) where d'Alembert operator becomes
!   g(r) Box[phi] = [-(d/dt)^2 + (d/dx)^2 + 2 g(r)/r (d/dx)] phi(x,t)
! which is a 1D wave equation with radial damping profile leading to back-scatter
! 
! as we are looking for smooth solutions, the method of choice to calculate
! spatial derivatives is pseudo-spectral, using Chebyshev basis on compactified
! coordinate y = x/sqrt(l^2+x^2), or inversely x/l = y/sqrt(1-y^2)
! [http://www.isbnsearch.org/isbn/0486411834]
! 
! absorbing boundary conditions are implemented using perfectly matched layers,
! applied to flux-conservative form of the free wave equation along the lines
!   du/dt = dv/dx - gamma*u, dv/dt = du/dx - gamma*v
! with auxilliary variables u = dphi/dt, v = dphi/dx and damping factor gamma
! [http://math.mit.edu/~stevenj/18.369/pml.pdf]
! 
! time integration is done using Gauss-Legendre method, which is A-stable
! and symplectic for Hamiltonian problems, as we have here
! [http://www.jstor.org/stable/2003405]
! 
! spherical collapse in scalar-tensor theories usually reduces to this type
! of PDE's given a proper field and coordinate choice, e.g. hep-th/0409117
! 
! Andrei Frolov <frolov@sfu.ca>; provided as is, modify or redistribute
! please cite arXiv:1704.04114 if you use this code for your research


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! QUICK START:
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! compile and link with LAPACK and CFITSIO libraries (use make or examples below)
! ifort -O3 -ipo -heap-arrays 256 -r8 {-parallel|-qopenmp} -lmkl_rt -lcfitsio spectral.f90 fitsio.f90 <model>.f90
! gfortran -O3 -fdefault-real-8 {-fopenmp} -llapack -lcfitsio spectral.f90 fitsio.f90 <model>.f90
! 
! run as: ./a.out > DATA; plot output using gnuplot: splot 'DATA' u 2:1:4 w l
! basic stop-frame animations can be done with gnuplot, for example using:
! set yrange [-1:1]; frame=0; plot 'DATA' index frame u 2:4 w l; load 'animate.gpl'
! 
! full evolution history is resampled onto an uniform grid and output as a FITS file,
! which can be imported into Python for advanced plotting (using PyFITS/Astropy)
! see bundled Python script 'plot.py' for the glue code and matplotlib examples
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program wave; use fitsio; use massive, only : phi0, DV, DDV; implicit none

! solver control parameters
integer, parameter :: nn = 2**11                ! number of nodes to sample on (i.e. spectral order)
integer, parameter :: tt = 2**17                ! total number of time steps to take (i.e. runtime)
integer, parameter :: nt = 2**4                 ! log output every nt time steps (i.e. downsampling)
real,    parameter :: dt = 0.005                ! time step size (simulated timespan is tt*dt)
real,    parameter :: ell = 64.0                ! grid compactification scale (set to 2.0 if in doubt)

! output control parameters
integer, parameter :: pml = 2**7                ! perfectly matched layer (supported on 3*pml nodes)
integer, parameter :: pts = 2**11 + 1           ! number of points on an uniform-spaced output grid
real,    parameter :: x0 = (pts-1)*(nt*dt/2.0)  ! output spans the range of x in an interval [-x0,x0]

logical, parameter :: output$log = .true.       ! output plain text log to stdout while code is running
logical, parameter :: output$fit = .true.       ! output binary data to FITS file when simulation ends

! this is exactly what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028841971694Q0

! collocation grid, metric functions, force term, and damping profile
real theta(nn), x(nn), r(nn), g(nn), Veff(nn), F(nn), gamma(nn)

! phase space state vector packs [phi,u,v,w] into a single array
real state(4*nn)

! spectral derivative, Laplacian, and prolongation operators
real D(nn,nn), L(nn,nn), Q(pts,nn)

! field evolution history for binary output, resampled on an uniform grid
! this can be a rather large array, so it is best to allocate it dynamically
real, allocatable :: history(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer i

! allocate storage for field evolution
if (output$fit) allocate(history(3,pts,tt/nt+1), source=0.0)

! initialize grid and linear operators (derivative D, Laplacian L, and prolongation Q)
call initg(); call initl()

! sanity checks on the time step value selected
if (dt > x(nn/2+1)-x(nn/2)) pause "Time step violates Courant condition, do you really want to run?"
if (sqrt(DDV(phi0))*dt > 0.25) pause "Time step will not resolve mass scale, do you really want to run?"

! sanity checks on the output grid extent selected
if (x0 > maxval(x, gamma == 0.0)) pause "Output grid will include boundary layer, do you really want to run?"

! force term is absent in scattering problem
F = 0.0

! initial field and velocity profiles
associate (phi => state(1:nn), u => state(nn+1:2*nn), v => state(2*nn+1:3*nn), w => state(3*nn+1:4*nn))
        phi = exp(-(x/8.0-7.5)**2) * (60.0/r) ! ingoing Gaussian wavepacket
        u = matmul(D,phi) + (g/r) * phi; v = (r*r) * matmul(D,phi); w = 0.0
end associate

! output initial conditions
call dump(0.0, state, history(:,:,1))

! main time evolution loop
do i = 1,tt
        call gl10(state, dt); if (mod(i,nt) == 0) call dump(i*dt, state, history(:,:,i/nt+1))
end do

! write entire evolution history into a FITS file
if (output$fit) call write2fits('output.fits', history, [-x0,x0], [0.0,tt*dt], ['phi','pi','flux'], '(x,t)')

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! evaluate ramp(x) = log(1+exp(x)) in an accurate way
elemental function ramp(x)
        real ramp, x; intent(in) x
        
        if (x >= 0.0) then; ramp = x + log(1.0 + exp(-x))
        else; ramp = 2.0*atanh(exp(x)/(2.0+exp(x))); end if
end function ramp

! areal radius r from tortoise coordinate x (in units of 2M)
elemental function radius(x)
        real radius, x, q, dr; intent(in) x; integer i
        
        ! closer to horizon, radius is 2M to all significant digits
        radius = 1.0; if (x < -digits(x)*log(2.0)) return
        
        ! refine q instead of r = 1.0 + log(1.0+exp(q)); initial guess
        q = x - 1.0
        
        ! Newton's iteration for x = r + log(r-1.0) as a function of q
        do i = 1,16
                dr = ramp(q); q = q - (1.0 + dr + log(dr) - x) * (dr + exp(log(dr)-q))/(1.0+dr)
        end do
        
        ! find converged radius value
        radius = 1.0 + ramp(q)
end function radius

! initialize Gauss-Lobatto collocation grid
subroutine initg()
        integer i; integer, parameter :: l = 1
        
        ! Chebyshev collocation grid (extrema and roots-of varieties)
        !forall (i=1:nn) theta(i) = (nn-i)*pi/(nn-1) ! includes interval ends
        forall (i=1:nn) theta(i) = (nn-i+0.5)*pi/nn ! excludes interval ends
        
        ! tortoise coordinate and metric functions
        x = ell/tan(theta); r = radius(x); g = 1.0 - 1.0/r; Veff = l*(l+1)/(r*r)
        
        ! PML absorption profile is truncated Gaussian (compactly supported on about 3*pml nodes each)
        forall (i=1:nn) gamma(i) = exp(-real(i-1)**2/pml**2) + exp(-real(i-nn)**2/pml**2) - 1.25e-4
        where (gamma < 0.0) gamma = 0.0; gamma = (0.5/dt)/gamma(1) * gamma
end subroutine initg

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(n, pts, theta, Tn, Tnx, Tnxx)
        integer n, pts; real, dimension(pts), intent(in) :: theta
        real, dimension(pts), intent(out), optional :: Tn, Tnx, Tnxx
        
        ! Chebyshev basis and its derivatives
        if (present(Tn))   Tn = cos(n*theta)
        if (present(Tnx))  Tnx = n * sin(n*theta) * sin(theta)**2/ell
        if (present(Tnxx)) Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3/ell**2
end subroutine evalb

! initialize linear spectral operator matrices (derivative D, Laplacian L, and prolongation Q)
subroutine initl()
        integer, parameter :: ia = 1, ib = nn, la = ib+1, lb = ib+nn, xa = lb+1, xb = lb+pts, ops = xb
        integer i, pivot(nn), status; real x(pts), grid(pts), A(nn,nn), B(nn,ops)
        
        ! output is resampled onto a uniform grid in x
        forall (i=1:pts) x(i) = (2*i-1-pts)*x0/(pts-1); grid = acos(x/sqrt(ell**2 + x**2))
        
        ! evaluate basis and differential operator values on collocation and output grids
        do i = 1,nn; associate (basis => A(i,:), grad => B(i,ia:ib), laplacian => B(i,la:lb), resampled => B(i,xa:xb))
                call evalb(i-1, nn, theta, basis, Tnx=grad, Tnxx=laplacian); laplacian = laplacian + (2.0*g/r) * grad - (g*Veff) * basis
                call evalb(i-1, pts, grid, resampled)
        end associate; end do
        
        ! find linear operator matrices
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, ops, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, ops, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! to evaluate derivative of a function f, simply do matmul(D,f)
        D = transpose(B(:,ia:ib))
        
        ! to evaluate Laplacian of a function f, simply do matmul(L,f)
        L = transpose(B(:,la:lb))
        
        ! to resample a function f to output grid, simply do matmul(Q,f)
        Q = transpose(B(:,xa:xb))
end subroutine initl

! evaluate equations of motion
subroutine evalf(y, dydt)
        real, dimension(4*nn) :: y, dydt
        
        ! unmangle phase space vector contents into human-readable form
        associate (phi => y(1:nn), u => y(nn+1:2*nn), v => y(2*nn+1:3*nn), w => y(3*nn+1:4*nn), &
                  dphi => dydt(1:nn), dudt => dydt(nn+1:2*nn), dvdt => dydt(2*nn+1:3*nn), dwdt => dydt(3*nn+1:4*nn))
                dphi = u-w; dudt = matmul(D,v)/(r*r) - gamma*u; dvdt = matmul(D,u-w)*(r*r) - gamma*v; dwdt = g*(DV(phi) - F + Veff*phi)
        end associate
end subroutine evalf

! solve linear Laplace problem [L - m_eff^2] phi = RHS
! for massive scalar field, static phi = lsolve(g*m^2, -g*F)
function lsolve(m2eff, rhs)
        real lsolve(nn), m2eff(nn), rhs(nn), A(nn,nn), B(nn,1)
        integer i, pivot(nn), status
        
        ! set up Laplace equation for massive field
        A = L; do i = 1,nn
                A(i,i) = A(i,i) - m2eff(i)
                B(i,1) = rhs(i)
        end do
        
        ! find static solution by direct inversion
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, 1, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, 1, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! return solution
        lsolve = B(:,1)
end function lsolve

! find static solution for a non-linear scalar field
! (even linear problem benefits from two iterations)
subroutine static(phi)
        real phi(nn); integer i
        
        do i = 1,16
                phi = phi + lsolve(g*DDV(phi), g*(DV(phi) - F) - matmul(L,phi))
        end do
end subroutine static

! dump simulation data in plain text format
subroutine dump(t, state, output)
        integer i; real t, state(4*nn), output(3,pts), residual(nn); optional output
        
        ! store resampled output if requested
        if (output$fit .and. present(output)) then; associate (phi => state(1:nn), u => state(nn+1:2*nn), v => state(2*nn+1:3*nn), w => state(3*nn+1:4*nn))
                output(1,:) = matmul(Q,phi)
                output(2,:) = matmul(Q,u-w)
                output(3,:) = matmul(Q,v)*output(2,:)
        end associate; end if
        
        ! bail unless evolution log is requested
        if (.not. output$log) return
        
        ! this should vanish if integration is perfect
        residual = (r*r) * matmul(D,state(1:nn)) - state(2*nn+1:3*nn)
        
        ! dump solution on collocation nodes
        do i = 1,nn
                if (gamma(i) == 0.0) write (*,'(3F24.16,5G24.16)') t, x(i), r(i), state([0:3]*nn + i), residual(i)
        end do
        
        ! separate timesteps by empty lines (for gnuplot's benefit)
        write (*,'(A)') "", ""
end subroutine

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 4*nn
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

end