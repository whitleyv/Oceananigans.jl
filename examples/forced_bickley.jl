using JLD2
using Plots

using Oceananigans
using Oceananigans.Advection
using Oceananigans.AbstractOperations
using Oceananigans.OutputWriters
using Oceananigans.Grids
using Oceananigans.Fields
using Oceananigans.Forcings

grid = RegularCartesianGrid(size=(128, 128, 1),
                            x = (-4π, 4π), y=(-4π, 4π), z=(0, 1),
                            topology = (Periodic, Bounded, Bounded))

C(y) = sin(y)
U(y) = sech(y)

relax_u(x, y, z, t, u, τ) = 1/τ * (U(y) - u)
relax_c(x, y, z, t, c, τ) = 1/τ * (C(y) - c)

force_u = Forcing(relax_u, field_dependencies=:u, parameters=1e-2)
force_c = Forcing(relax_c, field_dependencies=:c, parameters=1e-2)

model = IncompressibleModel(timestepper = :RungeKutta3, 
                              advection = WENO5(),
                                   grid = grid,
                               buoyancy = nothing,
                                tracers = :c,
                                closure = IsotropicDiffusivity(ν=1e-5),
                                forcing = (u=force_u, c=force_c)
                           )

# ## Initial conditions
#
# Our initial condition randomizes `model.velocities.u` and `model.velocities.v`.
# We ensure that both have zero mean for aesthetic reasons.

set!(model,
     u = (x, y, z) -> sech(y) - 0.1 * sech(y) * sin(x) * (tanh(y) * sin(y) - cos(y)),
     v = (x, y, z) -> 0.1 * sech(y) * cos(x) * sin(y),
     c = (x, y, z) -> sin(y),
)

u, v, w = model.velocities

ω = ComputedField(∂x(v) - ∂y(u))

progress(sim) = @info "Iteration: $(sim.model.clock.iteration), time: $(round(Int, sim.model.clock.time))"

simulation = Simulation(model, Δt=0.2, stop_time=200, iteration_interval=100, progress=progress)

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, (ω=ω,)),
                                                      schedule = TimeInterval(2),
                                                      prefix = "forced_bickley",
                                                      force = true)

run!(simulation)

file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

xu, yu, zu = nodes(model.velocities.u)
xv, yv, zv = nodes(model.velocities.v)
xω, yω, zω = nodes(ω)

@info "Making a neat movie of vorticity and speed..."

anim = @animate for (i, iteration) in enumerate(iterations)

    local u
    local ω

    @info "Plotting frame $i from iteration $iteration..."
    
    t = file["timeseries/t/$iteration"]
    ω = file["timeseries/ω/$iteration"][:, :, 1]
    u = file["timeseries/u/$iteration"][:, :, 1]

    ωlim = 2.0
    ωlevels = range(-ωlim, stop=ωlim, length=20)

    ulim = 1.0
    ulevels = range(-ulim, stop=ulim, length=20)

    kwargs = (xlabel="x", ylabel="y", aspectratio=1, linewidth=0, colorbar=true,
              xlims=(-model.grid.Lx/2, model.grid.Lx/2), ylims=(-model.grid.Ly/2, model.grid.Ly/2))
 
    u_plot = contourf(xu, yu, clamp.(u, -ulim, ulim)';
                       color = :balance,
                      levels = ulevels,
                       clims = (-ulim, ulim),
                      kwargs...)
             
    ω_plot = contourf(xω, yω, clamp.(ω, -ωlim, ωlim)';
                       color = :balance,
                      levels = ωlevels,
                       clims = (-ωlim, ωlim),
                      kwargs...)

    u_title = @sprintf("u at t = %.1f", t)
    ω_title = @sprintf("ω at t = %.1f", t)

    plot(u_plot, ω_plot, title=[u_title ω_title], layout=(1, 2), size=(1200, 500))
end

gif(anim, "forced_bickley.gif", fps = 8) # hide
