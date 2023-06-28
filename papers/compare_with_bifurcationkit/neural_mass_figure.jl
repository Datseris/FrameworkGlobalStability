# This script plots the comparison between BifurcationKit.jl
# and Attractors.jl. It first runs the two dedicated scripts.
using DrWatson
@quickactivate "FrameworkGlobalStability"
using CairoMakie, Attractors
include(srcdir("vis", "theme.jl"))
include("neural_mass_bifurcationkit.jl")
include("neural_mass_attractors.jl")

# %%
fig, axs = subplotgrid(2, 2; sharex = true, xlabels = "parameter",
    # titles = ["BifurcationKit.jl", "Attractors.jl"],
    titles = ["Traditional continuation-based bifurcation analysis (CBA)",
    "Recurrences-based attractor find-and-match continuation (RAFM)"],
	# ylabels = ["max(E) - step 1", "max(E) - step 2"],
    resolution = (1600, 600)
)

axs[1, 1].ylabel = "max(state), step 1"
axs[2, 1].ylabel = "max(state), step 2"
ylims!(axs[2, 1], (0, 25))
axs[1, 2].ylabel = "fractions %"
axs[2, 2].ylabel = "max(attractor)"
ylims!(axs[2, 2], (0, 25))
# reduce fontsize of axis titles slightly (they are so long)
axs[1, 1].titlesize = axs[1, 1].titlesize[] - 2
axs[1, 2].titlesize = axs[1, 2].titlesize[] - 2

fig

# left comumn is BifurcationKit.jl

# plot the branches in appropriate style; in the future
# BifurcationKit.jl will support Makie.
for ax in axs[:, 1]
for (p, E, s) in zip(curvesp_fp, curvesE_fp, stabilities_fp)
    lines!(ax, p, E;
        linestyle = s ? :solid : :dot, color = COLORS[2],
        label = s ? "stable" : "unstable", linewidth = 4,
    )
end
end
# scatterplot the special points
point_to_color = Dict(
    :bp => COLORS[3],
    :hopf => COLORS[4],
    :nd => COLORS[5],
    :ns => COLORS[6],
    :pd => "red",
)
point_to_marker = Dict(
    :bp => :circle,
    :hopf => :rect,
    :nd => :diamond,
    :ns => :utriangle,
    :pd => :rect,
)

for ax in axs[:, 1]
for sp in br.specialpoint
    t = sp.type
    t == :endpoint && continue
    p = sp.param
    E = sp.x[1]
    scatter!(ax, [p], [E]; markersize = 20,
        label = string(t),
        marker = point_to_marker[t], color = "transparent",
        strokewidth = 2.5, strokecolor = point_to_color[t],

    )
end
end
axislegend(axs[1, 1]; position = :rb, unique = true, patchsize = (30.0f0, 20.0f0))

p = br_potrap.param
E = br_potrap.max
isstable = br_potrap.stable
curvesp, curvesE, stabilities = segment_curves(p, E, isstable)

ax = axs[2, 1]

legend_elements = []
legend_entries = String[]
k = 0
for (p, E, s) in zip(curvesp_lc, curvesE_lc, stabilities_lc)
	l = lines!(ax, p, E;
		linestyle = s ? :solid : :dot, color = COLORS[1],
        linewidth = 4,
	)
	if s && k < 1
		push!(legend_elements, l)
		push!(legend_entries, "periodic")
        k += 1
	end
end

# and also plot the newly found bifurcation points
used_ns = false
for sp in br_potrap.γ.specialpoint
	t = sp.type
	t == :endpoint && continue
	p = sp.param
	E = br_potrap[sp.idx].max
	l = scatter!(ax, [p], [E];
        markersize = 15, color = "transparent", label = string(t),
        marker = point_to_marker[t],
        strokewidth = 2.5, strokecolor = point_to_color[t],
	)
    if t ∉ (:bp, :hopf) && k < 3
        if t == :ns && !used_ns
            used_ns = true
        elseif t == :ns
            continue
        end
        @show t
        push!(legend_elements, l)
        push!(legend_entries, string(t))
        k += 1
	end
end

axislegend(ax, legend_elements, legend_entries;
    position = :rt, patchsize = (30.0f0, 20.0f0)
)

# Right column is Attractors.jl
plot_basins_attractors_curves!(axs[1, 2], axs[2, 2],
    fractions_curves, attractors_info, attractor_to_real, E0_range;
    colors = Dict(-1 => COLORS[3], 1 => COLORS[2], 2 => COLORS[1], 3 => COLORS[4]),
    labels = Dict(-1 => "diverge", 1 => "fixed point", 2 => "limit cycle", 3 => "fixed point"),
    axislegend_kwargs = (position = :rt,),
)


fig


# notice that the stable periodic orbit exists _before_ the point
# that Bifurcation kit finds it. I've checked extensively, and indeed,
# there is a really slow converging periodic orbit at a parameter
# _larger_ than the hopf bifurcation parameter.

# %%
wsave(papersdir("figures", "figure4_comparison.png"), fig)
