# # Double Gyre
#
# This example simulates a double gyre following:
# https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html

using Oceananigans
using Oceananigans.Utils
using Oceananigans.Grids
using Oceananigans.Fields
using Oceananigans.BoundaryConditions
using Oceananigans.Advection
using Oceananigans.Diagnostics
using Oceananigans.OutputWriters
using Oceananigans.AbstractOperations

using Printf

grid = RegularCartesianGrid(size=(64, 64, 1), x=(-1e6, 1e6), y=(-1e6, 1e6), z=(-1e3, 0),
                            topology=(Bounded, Bounded, Bounded))

# ## Turbulence closure
closure = AnisotropicBiharmonicDiffusivity(νh=1)

# ## Model building

model = IncompressibleModel(architecture = CPU(),
                            timestepper = :RungeKutta3, 
                            advection = UpwindBiasedFifthOrder(),
                            grid = grid,
                            coriolis = BetaPlane(latitude=45),
                            buoyancy = nothing,
                            tracers = nothing,
                            closure = closure)

# ## Initial conditions

eddy_u(x, y, z, x₀, y₀, U, R) = - U * (y - y₀) / R * exp(-((x - x₀)^2 + (y - y₀)^2) / 2R^2)
eddy_v(x, y, z, x₀, y₀, U, R) = + U * (x - x₀) / R * exp(-((x - x₀)^2 + (y - y₀)^2) / 2R^2)

U = 0.1
R = 1e5

uᵢ(x, y, z) = eddy_u(x, y, z, -5e5, -5e5, U, R) +
              eddy_u(x, y, z, +0.0,  0.0, U, R) + 
              eddy_u(x, y, z, +5e5, +5e5, U, R)

vᵢ(x, y, z) = eddy_v(x, y, z, -5e5, -5e5, U, R) +
              eddy_v(x, y, z, +0.0,  0.0, U, R) + 
              eddy_v(x, y, z, +5e5, +5e5, U, R)

set!(model, u=uᵢ, v=vᵢ)

# ## Simulation setup

max_Δt = 10day

wizard = TimeStepWizard(cfl=1.0, Δt=12hour, max_change=1.1, max_Δt=max_Δt)

# Finally, we set up and run the the simulation.

umax = FieldMaximum(abs, model.velocities.u)
vmax = FieldMaximum(abs, model.velocities.v)
wmax = FieldMaximum(abs, model.velocities.w)

wall_clock = time_ns()

function print_progress(simulation)
    model = simulation.model

    ## Print a progress message
    msg = @sprintf("i: %04d, t: %s, Δt: %s, umax = (%.1e, %.1e, %.1e) ms⁻¹, wall time: %s\n",
                   model.clock.iteration,
                   prettytime(model.clock.time),
                   prettytime(wizard.Δt),
                   umax(), vmax(), wmax(),
                   prettytime(1e-9 * (time_ns() - wall_clock))
                  )

    @info msg

    return nothing
end

simulation = Simulation(model, Δt=wizard, stop_time=480days, iteration_interval=10, progress=print_progress)

# ## Output

u, v, w = model.velocities

speed = ComputedField(u^2 + v^2)

outputs = merge(model.velocities, model.tracers, (speed=speed,))

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs,
                                                      schedule = TimeInterval(8day),
                                                      prefix = "radiating_eddies",
                                                      field_slicer = FieldSlicer(k=model.grid.Nz),
                                                      force = true)

run!(simulation)

# # A neat movie

x, y, z = nodes(model.velocities.u)

xlims = (-grid.Lx/2 * 1e-3, grid.Lx/2 * 1e-3)
ylims = (-grid.Ly/2 * 1e-3, grid.Ly/2 * 1e-3)

x_km = x * 1e-3
y_km = y * 1e-3

nothing # hide

# Next, we open the JLD2 file, and extract the iterations we ended up saving at,

using JLD2, Plots

file = jldopen("radiating_eddies.jld2")

iterations = parse.(Int, keys(file["timeseries/t"]))

# These utilities are handy for calculating nice contour intervals:

""" Returns colorbar levels equispaced from `(-clim, clim)` and encompassing the extrema of `c`. """
function divergent_levels(c, clim, nlevels=21)
    levels = range(-clim, stop=clim, length=nlevels)
    cmax = maximum(abs, c)
    return ((-clim, clim), clim > cmax ? levels : levels = vcat([-cmax], levels, [cmax]))
end

""" Returns colorbar levels equispaced between `clims` and encompassing the extrema of `c`."""
function sequential_levels(c, clims, nlevels=20)
    levels = range(clims[1], stop=clims[2], length=nlevels)
    cmin, cmax = minimum(c), maximum(c)
    cmin < clims[1] && (levels = vcat([cmin], levels))
    cmax > clims[2] && (levels = vcat(levels, [cmax]))
    return clims, levels
end

# Finally, we're ready to animate.

@info "Making an animation from the saved data..."

anim = @animate for (i, iter) in enumerate(iterations)
    
    @info "Drawing frame $i from iteration $iter \n"

    t = file["timeseries/t/$iter"]
    u = file["timeseries/u/$iter"][:, :, 1]
    v = file["timeseries/v/$iter"][:, :, 1]
    s = file["timeseries/speed/$iter"][:, :, 1]

    ulims, ulevels = divergent_levels(u, U/2)
    slims, slevels = sequential_levels(s, (0.0, U/2))

    kwargs = (aspectratio=:equal, linewidth=0, xlims=xlims,
              ylims=ylims, xlabel="x (km)", ylabel="y (km)")

    u_plot = contourf(x_km, y_km, u';
                      color = :balance,
                      clims = ulims,
                      levels = ulevels,
                      kwargs...)
                        
    s_plot = contourf(x_km, y_km, s';
                      color = :thermal,
                      clims = slims,
                      levels = slevels,
                      kwargs...)
                             
    plot(u_plot, s_plot, size=(800, 500),
         title = ["u(t="*string(round(t/day, digits=1))*" day)" "speed"])

    iter == iterations[end] && close(file)
end

gif(anim, "radiating_eddies.gif", fps = 8) # hide
