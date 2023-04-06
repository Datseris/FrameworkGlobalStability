using DynamicalSystems # our framework implementation
using OrdinaryDiffEq   # high-accuracy ODE solvers

# initialize Lorenz84 as a `DynamicalSystem`
function lorenz84_rule(u, p, t)
    F, G, a, b = p
    x, y, z = u
    dx = -y^2 -z^2 -a*x + a*F
    dy = x*y - y - b*x*z + G
    dz = b*x*y + x*z - z
    return SVector(dx, dy, dz)
end
p0 = [6.886, 1.347, 0.255, 4.0]
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9)
ds = CoupledODEs(lorenz84_rule, ones(3), p0; diffeq)

# initialize recurrences-based algorithm
xg = yg = zg = range(-3, 3; length = 600)
grid = (xg, yg, zg)
mapper = AttractorsViaRecurrences(ds, grid;
    mx_chk_fnd_att = 1000, mx_chk_loc_att = 2000,
    mx_chk_lost = 100, mx_chk_safety = 1e8,
    Dt = 0.05, force_non_adaptive = true,
)

# find and continue attractors
prange = range(1.34, 1.37; length = 101)
pidx = 2
sampler = statespace_sampler(grid)[1]
rsc = RecurrencesSeededContinuation(mapper)
fractions_curves, attractors_info = continuation(
    rsc, prange, pidx, sampler
)