using DrWatson
@quickactivate "FrameworkGlobalStability"
using ChaosTools, OrdinaryDiffEq
using ChaosTools.DynamicalSystemsBase
using Random
include(srcdir("vis", "basins_plotting.jl"))

ds = Systems.lorenz()
M = 200
xg = yg = range(-25, 25; length = M)
zg = range(0, 60; length = M)
grid = (xg, yg, zg)
diffeq = (alg = Vern9(), reltol = 1e-9, abstol = 1e-9, maxiters = 1e12)


sampler, = statespace_sampler(Random.MersenneTwister(1234);
    min_bounds = minimum.(grid), max_bounds = maximum.(grid)
)

mapper = AttractorsViaRecurrencesSparse(ds, grid;
    mx_chk_fnd_att = 2000,
    mx_chk_loc_att = 100,
    mx_chk_att = 4,
    Ttr = 1000,
    safety_counter_max = 1e8,
    Δt = 0.1,
    diffeq,
)

# coexistance of periodic and chaotic, and then the chaotic collapses
# into the fixed point via a "crisis" (aka global bifurcation).
# stable fixed point exists always throughout the parameter range,
# but after the collapse, a fixed point and periodic attractor exist
prange = range(22.0, 28.0; length = 101)
pidx = 2
# threshold = 0.01 is the ε value we give at the mapper test
continuation = RecurrencesSeedingContinuation(mapper; threshold = Inf)
fractions_curves, attractors_info = basins_fractions_continuation(
    continuation, prange, pidx, sampler;
    show_progress = true, samples_per_parameter = 100
)

ukeys = unique_keys(fractions_curves)

# So if we get fractions_curves[80:90]
# we see that just after the transition of 3 to 2 attractors,
# we get  Dict(3 => 0.1400198609731877, 1 => 0.23535253227408143).
# So the fractions do not sum to 1.

# %%
animate_attractors_continuation(
    ds, attractors_info, prange, pidx;
    savename = "test.mp4", access = [1,2],
    limits = (-25,25,-25,25),
)

# Basins plot
# %%

# Barplot of all basins, order by keys
fig = basins_curves_plot(fractions_curves, prange)

# Makie.save("lorenz84_fracs.png", fig; px_per_unit = 4)
# negate_remove_bg("lorenz84_fracs.png")
