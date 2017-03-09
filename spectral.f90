! tutorial code solving 1D wave equation Box[phi] = m^2 phi in flat metric
!   ds^2 = -dt^2 + dx^2, 
! where d'Alembert operator is
!   Box[phi] = [-(d/dt)^2 + (d/dx)^2] phi(x,t)
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
integer, parameter :: pts = 2**11 + 1           ! number of points on an uniform-spaced output grid
real,    parameter :: x0 = (pts-1)*(dt/2.0)     ! output spans the range of x in an interval [-x0,x0]

! this is exactly what you think it is...
real, parameter :: pi = 3.1415926535897932384626433832795028841971694Q0

! collocation grid, metric functions, and spectral Laplacian operator
! phase space state vector v packs phi (1:nn) and dot(phi) (nn+1:nn+nn)
real theta(nn), x(nn), F(nn), L(nn,nn), Q(pts,nn), U(pts,nn), v(2*nn)

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

! initial field and velocity profiles
v(1:nn) = phi0*exp(-(x-0.5)**2*2.0); v(nn+1:2*nn) = 0.0

! force term corresponding to matter distributuion truncated at 6M
F = 0.0; !call static(v(1:nn))

! output initial conditions
call dump(0.0, v, history(:,:,1))

! main time evolution loop
do i = 1,tt
        call gl10(v, dt); call dump(i*dt, v, history(:,:,i+1))
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
        
        x = cos(theta)/sin(theta)
end subroutine initg

! evaluate rational Chebyshev basis on collocation grid theta
subroutine evalb(n, pts, theta, Tn, Dn, p, q)
        integer n, m, pts, mode; real, dimension(pts) :: theta, p, q, Tn, Tnx, Tnxx, Dn
        intent(in) n, pts, theta, p, q; intent(out) Tn, Dn; optional p, q, Dn
        
        ! Chebyshev basis and its derivatives
        Tn = cos(n*theta)
        if (.not. present(Dn)) return
        Tnx = n * sin(n*theta) * sin(theta)**2
        !Tnxx = -n * (n*cos(n*theta)*sin(theta) + 2.0*sin(n*theta)*cos(theta)) * sin(theta)**3
        
        ! Galerkin inner product (Tnxx,Tn)
        Tnxx = -(3.0/8.0)*n*n * Tn
        m = n+2; if (m < nn) Tnxx = Tnxx + n*(n+1)*cos(m*theta)/4.0
                 if (m == 3) Tnxx = Tnxx + cos(m*theta)/16.0
        m = n-2; if (m >= 0) Tnxx = Tnxx + n*(n-1)*cos(m*theta)/4.0
                 if (m == 0) Tnxx = Tnxx + (1.0/2.0)
                 if (m == 1) Tnxx = Tnxx - 3*cos(m*theta)/16.0
        m = n+4; if (m < nn) Tnxx = Tnxx - n*(n+2)*cos(m*theta)/16.0
        m = n-4; if (m >= 0) Tnxx = Tnxx - n*(n-2)*cos(m*theta)/16.0
                 if (m == 0) Tnxx = Tnxx - (1.0/2.0)
        
        ! Dn is a linear combination of first and second derivatives with coefficients p and q
        mode = 0; if (present(p)) mode = mode+2; if (present(q)) mode = mode+1
        
        ! Laplacian operator for Schwarzschild metric has q = 2.0*g/r
        select case (mode)
        	case(3); Dn = p*Tnxx + q*Tnx    ! both p and q are specified
        	case(2); Dn = p*Tnx             ! q is omitted, assume q = 0
        	case(1); Dn = Tnxx + q*Tnx      ! p is omitted, assume p = 1
        	case(0); Dn = Tnxx              ! both p and q are omitted
        end select
end subroutine evalb

! initialize Laplacian operator matrix
subroutine initl()
        integer i, pivot(nn), status; real x(pts), grid(pts), area(pts), A(nn,nn), B(nn,nn+pts+pts)
        
        ! output is resampled onto a uniform grid in x
        forall (i=1:pts) x(i) = (2*i-1-pts)*x0/(pts-1)
        grid = acos(x/sqrt(1.0 + x**2)); area = 1.0
        
        ! evaluate basis and differential operator values on collocation and output grids
        do i = 1,nn; associate (basis => A(i,:), laplacian => B(i,1:nn), resampled => B(i,nn+1:nn+pts), gradient => B(i,nn+pts+1:nn+pts+pts))
                call evalb(i-1, nn, theta, basis, laplacian)
                call evalb(i-1, pts, grid, resampled, gradient, p=area)
        end associate; end do
        
        ! find linear operator matrices
        status = 0; select case (kind(A))
                case(4); call sgesv(nn, nn+pts+pts, A, nn, pivot, B, nn, status)
                case(8); call dgesv(nn, nn+pts+pts, A, nn, pivot, B, nn, status)
                case default; call abort
        end select
        
        ! bail at first sign of trouble
        if (status /= 0) call abort
        
        ! to evaluate Laplacian of function f, simply do matmul(L,f)
        L = transpose(B(:,1:nn))
        
        ! to resample function f to output grid, simply do matmul(Q,f)
        Q = transpose(B(:,nn+1:nn+pts))
        
        ! to resample gradient f to output grid, simply do matmul(U,f)
        U = transpose(B(:,nn+pts+1:nn+pts+pts))
end subroutine initl

! evaluate equations of motion
subroutine evalf(v, dvdt)
        real v(2*nn), dvdt(2*nn)
        
        ! unmangle phase space state vector contents into human-readable form
        associate (phi => v(1:nn), pi => v(nn+1:2*nn), dphi => dvdt(1:nn), dpi => dvdt(nn+1:2*nn))
                dphi = pi; dpi = matmul(L,phi) - (DV(phi) - F)
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
subroutine dump(t, v, output)
        integer i; real t, v(2*nn), output(3,pts); optional output
        
        ! store resampled output if requested
        if (present(output)) then
                output(1,:) = matmul(Q,v(1:nn))
                output(2,:) = matmul(Q,v(nn+1:nn+nn))
                output(3,:) = matmul(U,v(1:nn))*output(2,:)
        end if
        
        ! dump solution on collocation nodes
        do i = 1,nn
                write (*,'(2F24.16,2G24.16)') t, x(i), v(i), v(nn+i)
        end do
        
        ! separate timesteps by empty lines (for gnuplot's benefit)
        write (*,'(A)') "", ""
end subroutine

! 10th order implicit Gauss-Legendre integrator
subroutine gl10(y, dt)
        integer, parameter :: s = 5, n = 2*nn
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