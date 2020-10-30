# # Flow around a cylinder

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

grid = RegularCartesianGrid(size=(64, 64, 1), x=(-4, 10), y=(-5, 5), z=(0, 1))

Re = 1000
const R = 1

inside_cylinder(x, y, z) = (x^2 + y^2) <= R

model = IncompressibleModel(timestepper = :RungeKutta3, 
                              advection = UpwindBiasedFifthOrder(),
                                   grid = grid,
                               buoyancy = nothing,
                                tracers = nothing,
                                closure = IsotropicDiffusivity(ν=1/Re),
                      immersed_boundary = inside_cylinder
                           )

# ## Random initial conditions
#
# Our initial condition randomizes `model.velocities.u` and `model.velocities.v`.
# We ensure that both have zero mean for aesthetic reasons.

set!(model, u=1)

progress(sim) = @info @sprintf("Iteration: % 4d, time: %.2f, max(u): %.2f, min(u): %.2f",
sim.model.clock.iteration,
sim.model.clock.time,
maximum(sim.model.velocities.u.data),
minimum(sim.model.velocities.u.data))

simulation = Simulation(model, Δt=2e-2, stop_time=20, iteration_interval=10, progress=progress)

# ## Output

simulation.output_writers[:fields] = JLD2OutputWriter(model, model.velocities,
                                                      schedule = TimeInterval(0.5),
                                                      prefix = "flow_around_cylinder",
                                                      force = true)

run!(simulation)

file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

u, v, w = model.velocities

xu, yu, zu = nodes(u)

@info "Making a neat movie of vorticity and speed..."

anim = @animate for (i, iteration) in enumerate(iterations)

    @info "Plotting frame $i from iteration $iteration..."
    
    t = file["timeseries/t/$iteration"]

    u_slice = file["timeseries/u/$iteration"][:, :, 1]
    @info maximum(u_slice) minimum(u_slice)
    u_max = maximum(abs, u_slice)
    u_lim = 0.8 * u_max

    u_levels = vcat([-u_max], range(-u_lim, stop=u_lim, length=20), [u_max])

    u_plot = contourf(xu, yu, u_slice';
                      linewidth = 0,
                          color = :balance,
                    aspectratio = 1,
                          title = @sprintf("u(x, y, t = %.1f) around a cylinder", t),
                         xlabel = "x",
                         ylabel = "y",
                         levels = u_levels,
                          xlims = (grid.xF[1], grid.xF[grid.Nx]),
                          ylims = (grid.yF[1], grid.yF[grid.Ny]),
                          clims = (-u_lim, u_lim))
end

gif(anim, "flow_around_cylinder.gif", fps = 8) # hide
