using DrWatson 
@quickactivate
using Attractors, Random, BenchmarkTools, DataFrames, StatsBase
include(srcdir("vis", "theme.jl"))

function benchmark_groups(features, feature_type)

    # reusable comparison function which is given to `produce_or_load`
    function group_benchmark(config)
        (; N, group_config, feature_type, group_method) = config
        labels = group_features(features, group_config)
        num_atts = length(labels)
        idxs_unique_labels = [findfirst(x->x==ulabel, labels) for ulabel in unique(labels)]
        atts = features[idxs_unique_labels]
        b = @benchmark(
            group_features($(features), $(group_config)),
            seconds = 10, # max benchmarking time (will take more if run costs more time)
        )
        b = median(b)
        time = b.time/1e9
        num_allocs = b.allocs
        memory = b.memory #bytes
        @info "$group_method in $feature_type has used $allocs and $memory and has found these $(length(atts)) attractors: $atts"
        outdata = @strdict(feature_type, group_method, N, time, num_allocs, memory)
        return outdata
    end


    N = length(features)
    configs = Dict()

    # optimal_radius_method = "silhouettes_optim"
    optimal_radius_method = 0.1
    config_dbscan = GroupViaClustering(; optimal_radius_method, rescale_features=false)
    configs["dbscan"] = config_dbscan

    distance_threshold = 0.1
    config_pairwise = GroupViaPairwiseComparison(; distance_threshold, rescale_features=false)
    configs["pairwise"] = config_pairwise 

    for (group_method, group_config) in configs 
        full_config = (; N, group_config, feature_type, group_method)
        @produce_or_load(
            group_benchmark, full_config, datadir("grouping_benchmarks", "groups");
            force = true, verbose = false, tag=false
        )
    end

end

function test_randomly_chosen_features()
    function generate_random_feats(num_total_feats=100, num_dims=10, num_unique_feats=5)
        unique_features = [rand(MersenneTwister(i), Float64, num_dims) for i in 1:num_unique_feats]
        features = [rand(MersenneTwister(Int64(i)), unique_features) for i in 1:num_total_feats]
    end

    for (idx, num_total_feats) in enumerate([1e2, 5e2, 1e3, 2e3, 5e3, 1e4, 2e4])
        features = generate_random_feats(num_total_feats)
        benchmark_groups(features, "random_features")
    end

end

function test_henon()
    henon_rule(x, p, n) = SVector{2}(1.0 - p[1]*x[1]^2 + x[2], p[2]*x[1])
    henon() = DeterministicIteratedMap(henon_rule, zeros(2), [1.4, 0.3])
    ds = henon()

    xg = yg = range(-2.0, 2.0; length=100)
    grid = (xg, yg)
    function featurizer(A, t)
        # Notice that unsupervised clustering cannot support "divergence to infinity",
        # which it identifies as another attractor (in fact, the first one).
        x = SVector(mean(A[:, 1]), mean(A[:, 2]))
        return any(isinf, x) ? SVector(200.0, 200.0) : x
    end
    config = GroupViaPairwiseComparison(; distance_threshold=0)
    mapper = AttractorsViaFeaturizing(ds, featurizer, config; T=500, Ttr = 500)
    
    num_total_feats_all = [100, 500, 1000, 5000, 10000]
    sampler, = statespace_sampler(grid, 1234)
    ics = StateSpaceSet([copy(sampler()) for i in 1:maximum(num_total_feats_all)])
    features = Attractors.extract_features(mapper, ics; show_progress=true)
    for (idx, num_total_feats) in enumerate(num_total_feats_all)
        features_reduced = features[1:num_total_feats]
        benchmark_groups(features_reduced, "henon")
    end
end

function plot_results_benchmark()
   # Load all benchmarks into a dataframe:
    results = collect_results(datadir("grouping_benchmarks", "groups"))
    sort!(results, ["feature_type", "group_method", "N"])
    results_by_type = groupby(results, "feature_type")
    feature_types = only.(keys(results_by_type)) # only gives us the elements in nice vector format
    methods = string.(sort!(unique(results[:, "group_method"])))

    fig, axs = subplotgrid(3, length(feature_types);
    sharex = true, titles = feature_types, resolution = (1000, 900))

    legend_elements = []
    display(fig)

    for (i, (group_method, benchmarks)) in enumerate(zip(feature_types, results_by_type))
        ax_times = axs[1,i]
        ax_alloc = axs[3,i]
        ax_mem = axs[2,i]
        grouped = groupby(benchmarks, "group_method")
        for (j, x) in enumerate(grouped)
            # notice that we drop the "method" column
            Ns, time_vs_N, allocs_vs_N, memory_vs_N = eachcol(x[:, ["N", "time", "num_allocs", "memory"]])
            # c = Cycled(j)
            # m = Cycled(j)
            # kwargs = (color = c, marker = m, markersize = 20, linewidth = 2,)
            kwargs = (markersize = 20, linewidth = 2,)
            scatterlines!(ax_alloc, log2.(Ns), log2.(allocs_vs_N); kwargs...)
            ele = scatterlines!(ax_times, log2.(Ns), log2.(time_vs_N); kwargs...)
            scatterlines!(ax_mem, log2.(Ns), log2.(memory_vs_N); kwargs...)
            if i == 1
                push!(legend_elements, ele)
                ax_times.ylabel = "log₂(time [sec])"
                ax_alloc.ylabel = "log₂(# allocs)"
                ax_mem.ylabel = "log₂(memory [byte])"
            end
            ax_alloc.xlabel = "log₂(# i.c.)"
        end
        # add titles or whatever
    end

    # add legend
    Legend(fig[end+1, :], legend_elements, methods;
    nbanks=length(methods)÷2, tellheight=true, tellwidth = false)
    display(fig)
    wsave(plotsdir("benchmarks", "groups.png"), fig) 
    return fig 
end

test_henon()
test_randomly_chosen_features()
plot_results_benchmark()
