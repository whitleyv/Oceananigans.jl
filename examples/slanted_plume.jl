using Statistics
using Plots
using JLD2
using Printf

using Oceananigans
using Oceananigans.Advection
using Oceananigans.Grids
using Oceananigans.Fields
using Oceananigans.Buoyancy
using Oceananigans.AbstractOperations
using Oceananigans.OutputWriters
using Oceananigans.Utils

ENV["GKSwstype"] = "100"

## Setup grid

N = 256
L = 1000
topology=(Bounded, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(1, N, N), x=(0, 1), y=(0, L), z=(-L, 0))

## Immersed boundary: ice shelf with constant 45° slope

inside_ice_shelf(x, y, z) = z >= -L + y

## Model setup

model = IncompressibleModel(
                   grid = grid,
            timestepper = :RungeKutta3,
              advection = WENO5(),
                closure = IsotropicDiffusivity(ν=1e-3, κ=1e-3),
      immersed_boundary = inside_ice_shelf
)

## Initial condition: Square blob warm anomaly

T₀(x, y, z) = 0.25L <= y <= 0.4L && -0.95L <= z <= -0.8L ? 20.0 : 19.0
set!(model, T=T₀)

## Simulation setup

progress(sim) = @info @sprintf("Iteration: % 4d, time: %.2f, extrema(w) = (%.4e, %.4e), extrema(T) = (%.4e, %.4e)",
                               sim.model.clock.iteration, sim.model.clock.time,
                               minimum(interior(sim.model.velocities.w)), maximum(interior(sim.model.velocities.w)),
                               minimum(interior(sim.model.tracers.T)), maximum(interior(sim.model.tracers.T)))

simulation = Simulation(model, Δt=5seconds, stop_time=1hour, iteration_interval=10, progress=progress)

## Output writing with diagnosed fields

u, v, w = model.velocities
T, S = model.tracers

diagnosed_fields = (
    ke = ComputedField(@at (Cell, Cell, Cell) (u^2 + v^2 + w^2) / 2),
    b  = BuoyancyField(model)
)

simulation.output_writers[:fields] = JLD2OutputWriter(model,
                                                      merge(model.velocities, model.tracers, diagnosed_fields),
                                                      schedule = TimeInterval(1minute),
                                                      prefix = "slanted_plume",
                                                      force = true)
## Run simulation!

run!(simulation)

## Animate simulation

file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

xc, yc, zc = nodes(T)

mask = inside_ice_shelf.(nodes(T, reshape=true)...)[1, :, :]

anim = @animate for (i, iteration) in enumerate(iterations)
    @info "Plotting frame $i from iteration $iteration..."

    t = file["timeseries/t/$iteration"]
    T_slice = file["timeseries/T/$iteration"][1, :, :]
    T_slice[mask] .= NaN

    T_plot = heatmap(yc, zc, T_slice';
                      linewidth = 0,
                          color = :thermal,
                    aspectratio = 1,
                          title = @sprintf("T(x, y, t = %s)", prettytime(t)),
                         xlabel = "y",
                         ylabel = "z",
                          xlims = (grid.yF[1], grid.yF[grid.Ny+1]),
                          ylims = (grid.zF[1], grid.zF[grid.Nz+1]),
                          clims = (19, 20),
                            dpi = 200
    )
end

gif(anim, "slanted_plume.gif", fps=15)
mp4(anim, "slanted_plume.mp4", fps=15)

## Plot kinetic energy statistics

ke_domain_mean = zeros(length(iterations))
ke_domain_max = zeros(length(iterations))
ke_immersed_mean = zeros(length(iterations))
ke_immersed_max = zeros(length(iterations))

for (i, iteration) in enumerate(iterations)
    ke = file["timeseries/ke/$iteration"][1, :, :]
    ke_domain_mean[i] = mean(ke[.~mask])
    ke_domain_max[i] = maximum(ke[.~mask])
    ke_immersed_mean[i] = mean(ke[mask])
    ke_immersed_max[i] = maximum(ke[mask])
end

times = [file["timeseries/t/$i"] for i in iterations] / minutes
ke_plot = plot(times[2:end], ke_domain_mean[2:end], linewidth=2, label="domain mean", yaxis=:log, color=:blue,
	       title="slanted plume kinetic energy statistics", xlabel="minutes", ylabel="kinetic energy",
	       xlims=extrema(times), grid=false, legend=:outertopright, framestyle=:box,
               foreground_color_legend=nothing, background_color_legend=nothing, dpi=200)

plot!(ke_plot, times[2:end], ke_domain_max[2:end], linewidth=2, label="domain max", color=:red)
plot!(ke_plot, times[2:end], ke_immersed_mean[2:end], linewidth=2, label="immersed mean", color=:blue, linestyle=:dot)
plot!(ke_plot, times[2:end], ke_immersed_max[2:end], linewidth=2, label="immersed max", color=:red, linestyle=:dot)

savefig(ke_plot, "slanted_plume_kinetic_energy.png")

## Plot volume-integrated buoyancy statistics

∫b_domain = zeros(length(iterations))
∫b_immersed = zeros(length(iterations))

for (i, iteration) in enumerate(iterations)
    b = file["timeseries/b/$iteration"][1, :, :]
    ∫b_domain[i] = sum(b[.~mask])
    ∫b_immersed[i] = sum(b[mask])
end

∫Δb_domain = (∫b_domain .- ∫b_domain[1]) / ∫b_domain[1]
∫Δb_immersed = (∫b_immersed .- ∫b_immersed[1]) / ∫b_immersed[1]

times = [file["timeseries/t/$i"] for i in iterations] / minutes
∫Δb_plot = plot(times, ∫Δb_domain, linewidth=2, label="domain", color=:blue,
		title="slanted plume ∫ Δb dV statistics", xlabel="minutes", ylabel="∫ Δb dV (relative change)",
	        xlims=extrema(times), grid=false, legend=:outertopright, framestyle=:box,
                foreground_color_legend=nothing, background_color_legend=nothing, dpi=200)

plot!(∫Δb_plot, times, ∫Δb_immersed, linewidth=2, label="immersed", color=:red)
plot!(∫Δb_plot, times, ∫Δb_domain .+ ∫Δb_immersed, linewidth=2, label="total", color=:green)

savefig(∫Δb_plot, "slanted_plume_mass.png")

