using DrWatson
@quickactivate "FrameworkGlobalStability"
using Attractors
using OrdinaryDiffEq: Vern9

# Neural mass model as in
# https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/tutorials/ode/tutorialsODE/#Neural-mass-equation-(Hopf-aBS)

function neural_mass_rule(z, p, t)
	@unpack J, α, E0, τ, τD, τF, U0 = p
	E, x, u = z
	SS0 = J * u * x * E + E0
	SS1 = α * log(1 + exp(SS0 / α))
	dz1 = (-E + SS1) / τ
	dz2 =	(1.0 - x) / τD - u * x * E
	dz3 = (U0 - u) / τF +  U0 * (1.0 - u) * E
	SVector(dz1, dz2, dz3)
end

p0 = (α = 1.5, τ = 0.013, J = 3.07, E0 = -2.0, τD = 0.200, U0 = 0.3, τF = 1.5, τS = 0.007)
p0 = ntuple2dict(p0) # make it mutable
z0 = [10.0, 0.982747, 0.367876]

# Use `verbose = false` option because there are unstable initial conditions
neural_mass = CoupledODEs(neural_mass_rule, z0, p0;
	diffeq = (verbose = false, alg = Vern9(), abstol = 1e-9, reltol = 1e-9)
)

# Compute attractors at a range with multstability
E0_multi = -1.75
set_parameter!(neural_mass, :E0, E0_multi)
# Define some arbitrarily large enough grid
density = 21
Eg = range(0, 30; length = density)
xg = range(0, 2; length = density)
ug = range(0, 2; length = density)
grid = (Eg, xg, ug)

# Default mapper
mapper = AttractorsViaRecurrences(
    neural_mass, grid;
    sparse = true,
	# this is set only to store the limit cycle with more detail,
	# which gives smoother continuation curves. Works perfectly fine with
	# commenting out these commands (and it is incredibly fast)
	mx_chk_fnd_att = 1000,
	mx_chk_loc_att = 1000,
	# store_once_per_cell = false,
)

using Random: Xoshiro
sampler, = statespace_sampler(Xoshiro(1234);
    min_bounds = minimum.(grid), max_bounds = maximum.(grid)
)
ics = StateSpaceSet([SVector{3}(sampler()) for _ in 1:1000])

set_parameter!(neural_mass, :E0, -1.472)

# fs, labels = basins_fractions(mapper, ics)

# Basins fractions
E0_range = range(-2, -0.9; length = 101)

rsc = RecurrencesFindAndMatch(mapper;
	threshold = 3.0,
	# Here we use a special distance: the distance between
	# maximum value of first dimension of attractors
	# (which is the same as the value plotted)
	distance = (A, B) -> abs(maximum(A[:, 1]) - maximum(B[:, 1]))
)

fractions_curves, attractors_info = Attractors.continuation(
    rsc, E0_range, :E0, sampler, samples_per_parameter = 1000,
)

# fig = basins_curves_plot(fractions_curves, E0_range)

# So easy... ;)

# Plot Attractors
attractor_to_real = A -> maximum(x[1] for x in A)
# attractors_curves_plot(attractors_info, attractor_to_real, E0_range)

# plot both!
fig = plot_basins_attractors_curves(fractions_curves, attractors_info, attractor_to_real, E0_range)
wsave(papersdir("figures", "neuralmass_attractors.png"), fig)
fig