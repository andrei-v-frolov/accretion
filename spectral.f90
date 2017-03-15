! tutorial code solving 1D wave equation Box[phi] = m^2 phi in flat metric
!   ds^2 = -dt^2 + dx^2, 
! where d'Alembert operator is
!   Box[phi] = [-(d/dt)^2 + (d/dx)^2] phi(x,t)
! 
! as we are looking for smooth solutions, the method of choice to calculate
! spatial derivatives is pseudo-spectral, using Chebyshev basis on compactified
! coordinate y = x/sqrt(1+x^2), or inversely x = y/sqrt(1-y^2)
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
! please cite arXiv:XXXX.XXXXX if you use this code for your research


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
integer, parameter :: nn = 2**9                 ! number of nodes to sample on (i.e. spectral order)
integer, parameter :: tt = 2**13                ! total number of time steps to take (i.e. runtime)
real,    parameter :: dt = 0.005                ! time step size (simulated timespan is tt*dt)

! output control parameters
integer, parameter :: pml = 2**4                ! perfectly matched layer depth
integer, parameter :: pts = 2**11 + 1           ! number of points on an uniform-spaced output grid
real,    parameter :: x0 = (pts-1)*(dt/2.0)     ! output spans the range of x in an interval [-x0,x0]

! this is exactly what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028841971694Q0

! collocation grid, metric functions, force term, damping, and phase space vector [phi,u,v,w]
real theta(nn), x(nn), F(nn), gamma(nn), state(4*nn)

! spectral integration, differentiation, Laplacian, and prolongation operators
real S(nn,nn), D(nn,nn), L(nn,nn), Q(pts,nn)

! field evolution history for binary output, resampled on an uniform grid
! this can be a rather large array, so it is best to allocate it dynamically
real, allocatable :: history(:,:,:)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer i

! allocate storage for field evolution
allocate(history(3,pts,tt+1), source=0.0)

! initialize grid and linear operators (Laplacian L and prolongation Q)
call initg(); call initl()

! sanity check on the time step value selected
if (dt > x(nn/2+1)-x(nn/2)) pause "Time step violates Courant condition, do you really want to run?"

! force term corresponding to matter distributuion truncated at 6M
F = 0.0

! initial field and velocity profiles
associate (phi => state(1:nn), u => state(nn+1:2*nn), v => state(2*nn+1:3*nn), w => state(3*nn+1:4*nn))
        phi = phi0*exp(-(x-1.0)**2*8.0); !call static(phi)
        u = 0.0; v = matmul(D,phi); w = 0.0
end associate

! output initial conditions
call dump(0.0, state, history(:,:,1))

! main time evolution loop
do i = 1,tt
        call gl10(state, dt); call dump(i*dt, state, history(:,:,i+1))
end do

! write entire evolution history into a FITS file
call write2fits('output.fits', history, [-x0,x0], [0.0,tt*dt], ['phi','pi','flux'], '(x,t)')

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! initialize Gauss-Lobatto collocation grid
subroutine initg()
        integer i
        
        ! Chebyshev collocation grid (extrema and roots-of varieties)
        !forall (i=1:nn) theta(i) = (nn-i)*pi/(nn-1) ! includes interval ends
        forall (i=1:nn) theta(i) = (nn-i+0.5)*pi/nn ! excludes interval ends
        
        ! PML absorption 
        
        x = cos(theta)/sin(theta); gamma = 0.0; !gamma(1:nn/2) = 10.0
        
        !forall (i=1:nn) gamma(i) = 100.0/cosh(9.0*(i-1)/(pml-1))**2 + 100.0/cosh(9.0*(i-nn)/(pml-1))**2
        forall (i=1:nn) gamma(i) = 100.0*(exp(-real(i-1)**2/pml**2) + exp(-real(i-nn)**2/pml**2))
end subroutine initg

! evaluate rational Chebyshev basis on a grid theta
subroutine evalb(n, pts, theta, Tn, Tnx, Tnxx)
        integer n, pts; real, dimension(pts), intent(in) :: theta
        real, dimension(pts), intent(out), optional :: Tn, Tnx, Tnxx
        
        ! Chebyshev basis and its derivatives
        if (present(Tn))   Tn = cos(n*theta)
        if (present(Tnx))  Tnx = n * sin(n*theta) * sin(theta)**2
        if (present(Tnxx)) Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3
end subroutine evalb

! integrate rational Chebyshev basis on a grid theta
subroutine integrate(nn, pts, theta, In)
        integer n, nn, pts
        real, dimension(pts), intent(in) :: theta
        real, dimension(nn,pts), intent(out) :: In
        
        ! we are expecting quite a large number of modes, so this should never trigger
        if (nn < 4) pause "Integration of Chebyshev basis is coded up for at least four modes"
        
        ! lowest four integrals
        In(1,:) = 1.0/tan(theta)
        In(2,:) = 1.0/sin(theta)
        In(3,:) = 2.0*theta + In(1,:)
        In(4,:) = (3.0 - 2.0*cos(2*theta))*In(2,:)
        
        ! compute higher order ones using recursion
        do n = 2,nn-3
                In(n+3,:) = (4.0/n)*sin(n*theta) + 2.0*In(n+1,:) - In(n-1,:)
        end do
end subroutine integrate

! initialize Laplacian operator matrix
subroutine initl()
        integer, parameter :: ia = 1, ib = nn, ka = ib+1, kb = ib+nn, la = kb+1, lb = kb+nn, xa = lb+1, xb = lb+pts, ops = xb
        integer i, pivot(nn), status; real x(pts), grid(pts), A(nn,nn), B(nn,ops)
        
        ! output is resampled onto a uniform grid in x
        forall (i=1:pts) x(i) = (2*i-1-pts)*x0/(pts-1); grid = acos(x/sqrt(1.0 + x**2))
        
        ! evaluate basis and differential operator values on collocation and output grids
        do i = 1,nn; associate (basis => A(i,:), grad => B(i,ka:kb), laplacian => B(i,la:lb), resampled => B(i,xa:xb))
                call evalb(i-1, nn, theta, basis, Tnx=grad, Tnxx=laplacian)
                call evalb(i-1, pts, grid, resampled)
        end associate; end do
        
        ! integrate basis on collocation grid
        call integrate(nn, nn, theta, B(:,ia:ib))
        
        ! find linear operator matrices
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, ops, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, ops, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! to evaluate Laplacian of function f, simply do matmul(L,f)
        S = transpose(B(:,ia:ib))
        D = transpose(B(:,ka:kb))
        L = transpose(B(:,la:lb))
        
        ! to resample function f to output grid, simply do matmul(Q,f)
        Q = transpose(B(:,xa:xb))
end subroutine initl

! evaluate equations of motion
subroutine evalf(y, dydt)
        real, dimension(4*nn) :: y, dydt
        
        ! unmangle phase space vector contents into human-readable form
        associate (phi => y(1:nn), u => y(nn+1:2*nn), v => y(2*nn+1:3*nn), w => y(3*nn+1:4*nn), &
                  dphi => dydt(1:nn), dudt => dydt(nn+1:2*nn), dvdt => dydt(2*nn+1:3*nn), dwdt => dydt(3*nn+1:4*nn))
                dphi = u-w; dudt = matmul(D,v) - gamma*u; dvdt = matmul(D,u-w) - gamma*v; dwdt = (DV(phi)-F) 
        end associate
end subroutine evalf

! solve linear Laplace problem [L - m_eff^2] phi = RHS
! for massive scalar field, static phi = lsolve(m2, -F)
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
                phi = phi + lsolve(DDV(phi), (DV(phi) - F) - matmul(L,phi))
        end do
end subroutine static

! dump simulation data in plain text format
subroutine dump(t, state, output)
        integer i; real t, state(4*nn), output(3,pts); optional output
        
        ! store resampled output if requested
        if (present(output)) then
                output(1,:) = matmul(Q,state(1:nn))
                output(2,:) = matmul(Q,state(nn+1:2*nn))
                output(3,:) = matmul(Q,state(2*nn+1:3*nn))
        end if
        
        ! dump solution on collocation nodes
        do i = 1,nn
                write (*,'(2F24.16,4G24.16)') t, x(i), state([0:3]*nn + i)
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