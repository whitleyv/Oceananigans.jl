# # Steady-state flow around a cylinder in 2D using immersed boundaries

using Statistics
using Plots
using JLD2
using Printf

using Oceananigans
using Oceananigans.Advection
using Oceananigans.Grids
using Oceananigans.Fields
using Oceananigans.AbstractOperations
using Oceananigans.OutputWriters
using Oceananigans.Utils

ENV["GKSwstype"] = "100"

N = 256
L = 1000
topology=(Bounded, Bounded, Bounded)
grid = RegularCartesianGrid(topology=topology, size=(1, N, N), x=(0, 1), y=(0, L), z=(-L, 0))

# Ice shelf with constant slope
θ_ice_shelf = 45.0
inside_ice_shelf(x, y, z) = z >= -L + y*sind(θ_ice_shelf)

localized_heat_flux(x, y, t) = 0.2L <= y <= 0.3L ? -1e-4 : 0.0
T_top_bc = FluxBoundaryCondition(localized_heat_flux)

T_bcs = TracerBoundaryConditions(grid, bottom=T_top_bc)

# setting up incompressible model with immersed boundary
model = IncompressibleModel(
                   grid = grid,
            timestepper = :RungeKutta3,
              advection = UpwindBiasedFifthOrder(),
                closure = IsotropicDiffusivity(ν=1e-3, κ=1e-3),
    boundary_conditions = (T=T_bcs,),
      immersed_boundary = inside_ice_shelf
)

# initial condition
# setting velocitiy to zero inside the cylinder and 1 everywhere else
T₀(x, y, z) = 0.3L <= y <= 0.4L && -0.95L <= z <= -0.85L ? 20.0 : 19.0
# T₀(x, y, z) = inside_ice_shelf(x, y, z) ? 20.0 : 20 # + 0.001z
set!(model, T=T₀)

progress(sim) = @info @sprintf("Iteration: % 4d, time: %.2f, extrema(w) = (%.4e, %.4e), extrema(T) = (%.4e, %.4e)",
                               sim.model.clock.iteration, sim.model.clock.time,
                               minimum(interior(sim.model.velocities.w)), maximum(interior(sim.model.velocities.w)),
                               minimum(interior(sim.model.tracers.T)), maximum(interior(sim.model.tracers.T)))

simulation = Simulation(model, Δt=5seconds, stop_time=1hour, iteration_interval=10, progress=progress)

# ## Output

simulation.output_writers[:fields] = JLD2OutputWriter(model,
                                                      merge(model.velocities, model.tracers),
                                                      schedule = TimeInterval(1minute),
                                                      prefix = "slanted_plume",
                                                      force = true)

# run it
run!(simulation)

# Analyze Results
file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

u, v, w = model.velocities
T, S = model.tracers

xc, yc, zc = nodes(T)

anim = @animate for (i, iteration) in enumerate(iterations)
    @info "Plotting frame $i from iteration $iteration..."

    t = file["timeseries/t/$iteration"]
    T_slice = file["timeseries/T/$iteration"][1, :, :]

    T_min, T_max = extrema(T_slice)
    T_levels = vcat([T_min], range(T_min, T_max, length=50), [T_max])

    v_plot = heatmap(yc, zc, T_slice';
                      linewidth = 0,
                          color = :thermal,
                    aspectratio = 1,
                          title = @sprintf("T(x, y, t = %s)", prettytime(t)),
                         xlabel = "x",
                         ylabel = "y",
                        #  levels = T_levels,
                          xlims = (grid.yF[1], grid.yF[grid.Ny]),
                          ylims = (grid.zF[1], grid.zF[grid.Nz]),
                          clims = (19, 20)
    )

    # plot!(circleShape(30,20,1),seriestype=[:shape,],linecolor=:black,
    #       legend=false,fillalpha=0, aspect_ratio=1)
end

gif(anim, "slanted_plume.gif", fps=15)
mp4(anim, "slanted_plume.mp4", fps=15)
