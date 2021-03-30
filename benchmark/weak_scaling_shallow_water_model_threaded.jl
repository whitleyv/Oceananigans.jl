using JLD2
using BenchmarkTools
using Benchmarks

threads = (1, 2, 4, 8, 16, 32)

grid_size(T) = (4096, 256T)

# Run benchmarks

print_system_info()

for T in threads
    Nx, Ny = grid_size(T)
    @info "Benchmarking serial shallow water model weak scaling with threading [N=($Nx, $Ny), T=$T]..."
    julia = Base.julia_cmd()
    run(`$julia -t $T --project weak_scaling_shallow_water_model_serial.jl $Nx $Ny`)
end

# Collect and merge benchmarks from all ranks

suite = BenchmarkGroup(["size", "threads"])

for T in threads
    Nx, Ny = grid_size(T)
    case = ((Nx, Ny), T)

    filename = string("weak_scaling_shallow_water_model_threads$(T).jld2")
    file = jldopen(filename, "r")
    suite[case] = file["trial"]
end

# Summarize benchmarks

df = benchmarks_dataframe(suite)
sort!(df, :threads)
benchmarks_pretty_table(df, title="Shallow water model weak scaling with multithreading benchmark")

base_case = (grid_size(1), threads[1])
suite_Δ = speedups_suite(suite, base_case=base_case)
df_Δ = speedups_dataframe(suite_Δ, slowdown=true, efficiency=:weak, base_case=base_case)
sort!(df_Δ, :threads)
benchmarks_pretty_table(df_Δ, title="Shallow water model weak multithreading scaling speedup")
