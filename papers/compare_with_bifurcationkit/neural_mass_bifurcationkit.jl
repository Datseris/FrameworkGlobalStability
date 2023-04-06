using DrWatson
@quickactivate "FrameworkGlobalStability"
using BifurcationKit, Setfield

using LinearAlgebra: norm
norminf(x) = norm(x, Inf)

# vector field
function TMvf!(dz, z, p, t)
	@unpack J, α, E0, τ, τD, τF, U0 = p
	E, x, u = z
	SS0 = J * u * x * E + E0
	SS1 = α * log(1 + exp(SS0 / α))
	dz[1] = (-E + SS1) / τ
	dz[2] =	(1.0 - x) / τD - u * x * E
	dz[3] = (U0 - u) / τF +  U0 * (1.0 - u) * E
	dz
end

# out of place method
TMvf(z, p) = TMvf!(similar(z), z, p, 0)

# parameter values
par_tm = (α = 1.5, τ = 0.013, J = 3.07, E0 = -2.0, τD = 0.200, U0 = 0.3, τF = 1.5, τS = 0.007)

# initial condition
z0 = [0.238616, 0.982747, 0.367876]

# Bifurcation Problem
prob = BifurcationProblem(TMvf, z0, par_tm, (@lens _.E0);
    recordFromSolution = (x, p) -> (E = x[1], x = x[2], u = x[3]),
)


# continuation options
opts_br = ContinuationPar(pMin = -10.0, pMax = -0.9,
	# parameters to have a smooth result
	ds = 0.04, dsmax = 0.05,
	# this is to detect bifurcation points precisely with bisection
	detectBifurcation = 3,
	# Optional: bisection options for locating bifurcations
	nInversion = 8, maxBisectionSteps = 25, nev = 3,
    # for obtaining the curves to plot them:
    saveSolEveryStep = 1,
)

# continuation of equilibria
@time br = BifurcationKit.continuation(prob, PALC(tangent=Bordered()), opts_br;
	plot = false, normC = norminf
)

# Okay, now get the bifurcation curves to plot. This is really hard...
E = br.E
p = br.param
isstable = br.stable
# we segment curves into stable and unstable:
function segment_curves(p, E, isstable)
	changes = findall(i -> isstable[i+1] != isstable[i], 1:length(isstable)-1)
	pushfirst!(changes, 1)
	push!(changes, length(isstable)-1)
	curvesE = [E[changes[i]:changes[i+1]] for i in 1:length(changes)-1]
	curvesp = [p[changes[i]:changes[i+1]] for i in 1:length(changes)-1]
	stabilities = isstable[changes[2:end]]
	return curvesp, curvesE, stabilities
end

curvesp_fp, curvesE_fp, stabilities_fp = segment_curves(p, E, isstable)


# Continue now with the periodic orbit
# newton parameters
optn_po = NewtonPar(verbose = false, tol = 1e-8,  maxIter = 8)

# continuation parameters
opts_po_cont = ContinuationPar(
	dsmax = 0.1, ds= -0.0001, dsmin = 1e-4, pMax = 0., pMin=-5.,
	# The higher the `maxSteps` here, the
	maxSteps = 120, newtonOptions = (@set optn_po.tol = 1e-7),
	nev = 3, tolStability = 1e-8, detectBifurcation = 3, saveSolEveryStep=1
)

# arguments for periodic orbits
args_po = (
	recordFromSolution = (x, p) -> begin
	xtt = BifurcationKit.getPeriodicOrbit(p.prob, x, @set par_tm.E0 = p.p)
	return (max = maximum(xtt[1,:]),
			min = minimum(xtt[1,:]),
			period = getPeriod(p.prob, x, @set par_tm.E0 = p.p))
	end,
	normC = norminf
)

Mt = 200 # number of time sections

# Wow, this takes 7 seconds or so.
# Our entire algorithm takes 1 second, including finding the fractions...
@time br_potrap = continuation(
	# we want to branch form the 4th bif. point
	br, 4, opts_po_cont,
	# we want to use the Trapeze method to locate PO
	# this jacobian is specific to ODEs
	# it is computed using AD and
	# updated inplace
	PeriodicOrbitTrapProblem(M = Mt, jacobian = :Dense, updateSectionEveryStep = 0);
	# regular continuation options
	verbosity = 0, plot = false,
	args_po...,
	callbackN = BifurcationKit.cbMaxNorm(1000.0)
)

# okay now lets plot the new branch.
p = br_potrap.param
E = br_potrap.max
isstable = br_potrap.stable
curvesp_lc, curvesE_lc, stabilities_lc = segment_curves(p, E, isstable)
