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

using Oceananigans.Fields: PressureField

using Printf

grid = RegularCartesianGrid(size=(128, 128, 1), x=(-8, 8), y=(-8, 8), z=(-1, 0),
                            topology=(Bounded, Bounded, Bounded))

# ## Turbulence closure
#closure = AnisotropicBiharmonicDiffusivity(νh=1)
closure = IsotropicDiffusivity(ν=1e-2)

# ## Model building

model = IncompressibleModel(architecture = CPU(),
                            timestepper = :RungeKutta3, 
                            advection = UpwindBiasedFifthOrder(),
                            grid = grid,
                            coriolis = BetaPlane(f₀=1, β=0.1),
                            buoyancy = nothing,
                            tracers = nothing,
                            closure = closure)

# ## Initial conditions

eddy_u(x, y, z, ϵ) = - ϵ * y * exp(-(x^2 + y^2))
eddy_v(x, y, z, ϵ) = + ϵ * x * exp(-(x^2 + y^2))

ϵ = 0.01

uᵢ(x, y, z) = eddy_u(x, y, z, ϵ)
vᵢ(x, y, z) = eddy_v(x, y, z, ϵ)

set!(model, u=uᵢ, v=vᵢ)

# ## Simulation setup

# Finally, we set up and run the the simulation.

umax = FieldMaximum(abs, model.velocities.u)
vmax = FieldMaximum(abs, model.velocities.v)
wmax = FieldMaximum(abs, model.velocities.w)

wall_clock = [time_ns()]

function print_progress(simulation)
    model = simulation.model

    ## Print a progress message
    msg = @sprintf("i: %04d, t: %.2f, Δt: %.2f, umax = (%.1e, %.1e, %.1e), wall time: %s\n",
                   model.clock.iteration,
                   model.clock.time,
                   simulation.Δt,
                   umax(), vmax(), wmax(),
                   prettytime(1e-9 * (time_ns() - wall_clock[1]))
                  )

    @info msg

    wall_clock[1] = time_ns()

    return nothing
end

simulation = Simulation(model, Δt=0.5, stop_time=60.0, iteration_interval=10, progress=print_progress)

# ## Output

u, v, w = model.velocities

p = PressureField(model)

outputs = merge(model.velocities, model.tracers, (p=p,))

simulation.output_writers[:fields] = JLD2OutputWriter(model, outputs,
                                                      schedule = TimeInterval(1.0),
                                                      prefix = "radiating_eddies",
                                                      force = true)

run!(simulation)

# # A neat movie

xu, yu, zu = nodes(model.velocities.u)
xv, yv, zv = nodes(model.velocities.v)
xp, yp, zp = nodes(p)

xlims = (-grid.Lx/2, grid.Lx/2)
ylims = (-grid.Ly/2, grid.Ly/2)

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
    p = file["timeseries/p/$iter"][:, :, 1]

    ulims, ulevels = divergent_levels(u, ϵ/10)
    plims, plevels = divergent_levels(p, maximum(abs, p) + 1e-9)

    kwargs = (aspectratio=:equal, linewidth=0, xlims=xlims,
              ylims=ylims, xlabel="x (km)", ylabel="y (km)")

    u_plot = contourf(xu, yu, clamp.(u, ulims[1], ulims[2])';
                      color = :balance,
                      clims = ulims,
                      levels = ulevels,
                      kwargs...)
                        
    v_plot = contourf(xv, yv, clamp.(v, ulims[1], ulims[2])';
                      color = :balance,
                      clims = ulims,
                      levels = ulevels,
                      kwargs...)

    p_plot = contourf(xp, yp, clamp.(p, plims[1], plims[2])';
                      color = :balance,
                      clims = plims,
                      levels = plevels,
                      kwargs...)
                             
    u_title = @sprintf("u at t = %.2f", t)
    v_title = @sprintf("v at t = %.2f", t)
    p_title = @sprintf("p at t = %.2f", t)

    plot(u_plot, v_plot, p_plot, size=(1600, 400), layout=(1, 3),
         title = [u_title v_title p_title])

    iter == iterations[end] && close(file)
end

gif(anim, "radiating_eddy.gif", fps = 8) # hide
