using Logging
using JLD2
using BenchmarkTools
using Benchmarks

using Oceananigans
using Oceananigans.Models
using CUDA

Logging.global_logger(OceananigansLogger())

Nx = parse(Int, ARGS[1])
Ny = parse(Int, ARGS[2])

T = Threads.nthreads()

@info "Setting up serial shallow water model with N=($Nx, $Ny) grid points and $T threads..."

 topo = (Periodic, Periodic, Bounded)   # Use Flat
 grid = RegularRectilinearGrid(topology=topo, size=(Nx, Ny, 1), extent=(1, 1, 1))
model = ShallowWaterModel(architecture=CPU(), grid=grid, gravitational_acceleration=1.0)
set!(model, h=1)

@info "Warming up serial shallow water model..."

time_step!(model, 1) # warmup

@info "Benchmarking serial shallow water model..."

trial = @benchmark begin
    @sync_gpu time_step!($model, 1)
    CUDA.@sync blocking=true time_step!($model, 1)
end samples=10

t_median = BenchmarkTools.prettytime(median(trial).time)
@info "Done benchmarking. Median time: $t_median"

jldopen("weak_scaling_shallow_water_model_threads$T.jld2", "w") do file
    file["trial"] = trial
end
