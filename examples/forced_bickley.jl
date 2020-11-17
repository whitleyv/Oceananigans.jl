using JLD2
using Plots
using Printf

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

force_u = Forcing(relax_u, field_dependencies=:u, parameters=1e2)
force_c = Forcing(relax_c, field_dependencies=:c, parameters=1e2)

model = IncompressibleModel(timestepper = :RungeKutta3, 
                              advection = WENO5(),
                                   grid = grid,
                               buoyancy = nothing,
                                tracers = :c,
                                closure = IsotropicDiffusivity(ν=1e-4, κ=1e-4),
                                forcing = (u=force_u, c=force_c))

# Initial conditions

ϵ = 0.01

set!(model,
     u = (x, y, z) -> sech(y)^2 * (1 + ϵ * randn()),
     v = (x, y, z) -> ϵ * sech(y)^2 * randn(),
     c = (x, y, z) -> sin(y / 4))

u, v, w = model.velocities

ω = ComputedField(∂x(v) - ∂y(u))

progress(sim) = @info(@sprintf("Iter: %d, time: %.1f, max|u|: %.2f",
                               sim.model.clock.iteration, sim.model.clock.time, maximum(abs, interior(u))))

wizard = TimeStepWizard(cfl=1.0, Δt=1e-1, max_change=1.1)

simulation = Simulation(model, Δt=wizard, stop_time=100, iteration_interval=100, progress=progress)

simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers, (ω=ω,)),
                                                      schedule = TimeInterval(1.0),
                                                      prefix = "forced_bickley",
                                                      force = true)

@info "Running a simulation of a forced Bickley jet..."

run!(simulation)

file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

xu, yu, zu = nodes(u)
xv, yv, zv = nodes(v)
xω, yω, zω = nodes(ω)
xc, yc, zc = nodes(model.tracers.c)

@info "Making a neat movie..."

anim = @animate for (i, iteration) in enumerate(iterations)

    local u
    local ω

    @info "Plotting frame $i from iteration $iteration..."
    
    t = file["timeseries/t/$iteration"]
    ω = file["timeseries/ω/$iteration"][:, :, 1]
    u = file["timeseries/u/$iteration"][:, :, 1]
    c = file["timeseries/c/$iteration"][:, :, 1]

    ωlim = 1
    ωlevels = range(-ωlim, stop=ωlim, length=20)

    ulim = 1
    ulevels = range(-ulim, stop=ulim, length=20)

    clim = 1
    clevels = range(-clim, stop=clim, length=20)

    kwargs = (xlabel="x", ylabel="y", aspectratio=1, linewidth=0, colorbar=true,
              xlims=(-model.grid.Lx/2, model.grid.Lx/2), ylims=(-model.grid.Ly/2, model.grid.Ly/2))
 
    u_plot = contourf(xu, yu, clamp.(u, -ulim, ulim)';
                       color = :balance,
                      levels = ulevels,
                       clims = (-ulim, ulim),
                      kwargs...)

    c_plot = contourf(xc, yc, clamp.(c, -clim, clim)';
                       color = :balance,
                      levels = clevels,
                       clims = (-clim, clim),
                      kwargs...)
             
    ω_plot = contourf(xω, yω, clamp.(ω, -ωlim, ωlim)';
                       color = :balance,
                      levels = ωlevels,
                       clims = (-ωlim, ωlim),
                      kwargs...)

    u_title = @sprintf("u at t = %.1f", t)
    ω_title = @sprintf("ω at t = %.1f", t)
    c_title = @sprintf("c at t = %.1f", t)

    plot(u_plot, ω_plot, c_plot, title=[u_title ω_title c_title], layout=(1, 3), size=(1500, 400))
end

gif(anim, "forced_bickley.gif", fps = 8) # hide
