### Scalar field accretion onto a static Schwarzschild-like black hole

Tutorial code solving 1D wave equation $$ \Box\phi = m^2 \phi $$ in Schwarzschild metric
$$ ds^2 = -g(r) dt^2 + dr^2/g(r) + r^2 d\Omega^2 $$, where $$ g(r) = 1 - \frac{2M}{r} $$.
In tortoise coordinates $$ dx = dr/g(r) $$ the d'Alembert operator becomes
$$ \Box\phi = \frac{1}{g(r)} \left[-\frac{d^2~}{dt^2} + \frac{d^2~}{dx^2} + 2\,\frac{g(r)}{r}\, \frac{d~}{dx}\right] \phi(x,t) $$,
which is a 1D wave equation with radial damping profile leading to back-scatter.

As we are looking for smooth solutions, the method of choice to calculate
spatial derivatives is pseudo-spectral, using Chebyshev basis on compactified
coordinate $$ y = x/\sqrt{1+x^2} $$, or inversely $$ x = y/\sqrt{1-y^2} $$

Time integration is done using [Gauss-Legendre method](http://www.jstor.org/stable/2003405), which is A-stable
and symplectic for Hamiltonian problems, as we have here.

Spherical collapse in scalar-tensor theories usually reduces to this type of
PDE's given a proper field and coordinate choice, e.g. [hep-th/0409117](https://arxiv.org/abs/hep-th/0409117).

Written by Andrei Frolov <frolov@sfu.ca>; provided as is, permission is granted to
modify or redistribute. Please cite [arXiv:1704.04114](https://arxiv.org/abs/1704.04114) if you use this code in
your research.
