# # Double Gyre
#
# This example simulates a double gyre following:
# https://mitgcm.readthedocs.io/en/latest/examples/baroclinic_gyre/baroclinic_gyre.html

using Oceananigans
using Oceananigans.Grids
using Oceananigans.Utils
using Oceananigans.BoundaryConditions
using Oceananigans.Advection
using Oceananigans.Diagnostics

using Printf

grid = RegularCartesianGrid(size=(48, 48, 16), x=(-2e6, 2e6), y=(-3e6, 3e6), z=(-2e3, 0),
                            topology=(Bounded, Bounded, Bounded))

# ## Boundary conditions

@inline wind_stress(x, y, t, parameters) = - parameters.τ * cos(2π * y / parameters.L)

@inline u_bottom_stress(i, j, grid, clock, model_fields, parameters) =
    @inbounds - parameters.μ * parameters.H * model_fields.u[i, j, 1]

@inline v_bottom_stress(i, j, grid, clock, model_fields, parameters) =
    @inbounds - parameters.μ * parameters.H * model_fields.v[i, j, 1]

wind_stress_bc = BoundaryCondition(Flux, wind_stress, parameters = (τ=1e-4, L=grid.Ly))

u_bottom_stress_bc = BoundaryCondition(Flux, u_bottom_stress,
                                       discrete_form=true, parameters=(μ=1/30day, H=grid.Lz))

v_bottom_stress_bc = BoundaryCondition(Flux, v_bottom_stress,
                                       discrete_form=true, parameters=(μ=1/30day, H=grid.Lz))

u_bcs = UVelocityBoundaryConditions(grid,
                                     north = BoundaryCondition(Value, 0),
                                     south = BoundaryCondition(Value, 0),
                                       top = wind_stress_bc,
                                    bottom = u_bottom_stress_bc)

v_bcs = VVelocityBoundaryConditions(grid,
                                      east = BoundaryCondition(Value, 0),
                                      west = BoundaryCondition(Value, 0),
                                    bottom = v_bottom_stress_bc)

w_bcs = WVelocityBoundaryConditions(grid,
                                    north = BoundaryCondition(Value, 0),
                                    south = BoundaryCondition(Value, 0),
                                     east = BoundaryCondition(Value, 0),
                                     west = BoundaryCondition(Value, 0))

b_reference(y, parameters) = parameters.Δb / parameters.Ly * y

@inline buoyancy_flux(i, j, grid, clock, model_fields, parameters) =
    @inbounds - parameters.μ * (model_fields.b[i, j, grid.Nz] - b_reference(grid.yC[j], parameters))

buoyancy_flux_bc = BoundaryCondition(Flux, buoyancy_flux,
                                     discrete_form = true,
                                     parameters = (μ=1/day, Δb=0.06, Ly=grid.Ly))

b_bcs = TracerBoundaryConditions(grid, 
                                 bottom = BoundaryCondition(Value, 0),
                                 top = buoyancy_flux_bc)

closure = AnisotropicDiffusivity(νh=500, νz=1e-2, κh=100, κz=1e-2)

model = IncompressibleModel(architecture = CPU(),
                            timestepper = :RungeKutta3, 
                            advection = UpwindBiasedFifthOrder(),
                            grid = grid,
                            coriolis = BetaPlane(latitude=45),
                            buoyancy = BuoyancyTracer(),
                            tracers = :b,
                            closure = closure,
                            boundary_conditions = (u=u_bcs, v=v_bcs, w=w_bcs, b=b_bcs))

## Temperature initial condition: a stable density gradient with random noise superposed.
b₀(x, y, z) = b_bcs.top.condition.parameters.Δb * (1 + z / grid.Lz)

set!(model, b=b₀)

# ## Running the simulation
#
# To run the simulation, we instantiate a `TimeStepWizard` to ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 0.2.

max_Δt = min(hour, 0.5 * min(grid.Δz^2 / closure.κz, grid.Δx^2 / closure.νx))

wizard = TimeStepWizard(cfl=0.5, Δt=hour/2, max_change=1.1, max_Δt=max_Δt)

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

simulation = Simulation(model, Δt=wizard, stop_time=20*365day, iteration_interval=100, progress=print_progress)

# ## Set up output
#
# We set up an output writer that saves all velocity fields, tracer fields, and the subgrid
# turbulent diffusivity associated with `model.closure`. The `prefix` keyword argument
# to `JLD2OutputWriter` indicates that output will be saved in
# `double_gyre.jld2`.

using Oceananigans.OutputWriters

#=
## Instantiate a JLD2OutputWriter to write fields. We will add it to the simulation before
## running it.
simulation.output_writers[:fields] = JLD2OutputWriter(model, merge(model.velocities, model.tracers);
                                                      time_interval=2day,
                                                      prefix="double_gyre",
                                                      field_slicer=FieldSlicer(k=model.grid.Nz),
                                                      force=true)

run!(simulation)
=#

# # Making a neat movie
#
# We look at the results by plotting vertical slices of $u$ and $w$, and a horizontal
# slice of $w$ to look for Langmuir cells.

# Making the coordinate arrays takes a few lines of code,

x, y, z = nodes(model.tracers.b)
nothing # hide

# Next, we open the JLD2 file, and extract the iterations we ended up saving at,

using JLD2, Plots

file = jldopen("double_gyre.jld2")

iterations = parse.(Int, keys(file["timeseries/t"]))
nothing # hide

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

    ## Load 3D fields from file, omitting halo regions
    u = file["timeseries/u/$iter"]
    v = file["timeseries/v/$iter"]
    w = file["timeseries/w/$iter"]
    t = file["timeseries/t/$iter"]

    ## Extract slices
    uxy = 1/2 * (u[1:end-1, :, end] .+ u[2:end, :, end])
    vxy = 1/2 * (v[:, 1:end-1, end] .+ v[:, 2:end, end])
    wxy = w[:, :, 1]
    
    speed = @. sqrt(uxy^2 + vxy^2)
    
    ulims, ulevels = divergent_levels(u, ulim)
    slims, slevels = sequential_levels(speed, (0.0, 1.0))

    xlims = (-grid.Lx/2 * 1e-3, grid.Lx/2 * 1e-3)
    ylims = (-grid.Ly/2 * 1e-3, grid.Ly/2 * 1e-3)

    uxy_plot = contourf(x / 1e3, y / 1e3, uxy';
                          linewidth = 0,
                              color = :balance,
                        aspectratio = :equal,
                              clims = ulims,
                             levels = ulevels,
                              xlims = xlims,
                              ylims = ylims,
                             xlabel = "x (km)",
                             ylabel = "y (km)")
                        
     wxy_plot = contourf(x / 1e3, y / 1e3, wxy';
                           linewidth = 0,
                               color = :balance,
                         aspectratio = :equal,
                               clims = ulims,
                              levels = ulevels,
                               xlims = xlims,
                               ylims = ylims,
                              xlabel = "x (km)",
                              ylabel = "y (km)")
                         
    speed_plot = contourf(x / 1e3, y / 1e3 , speed';
                          linewidth = 0,
                              color = :thermal,
                        aspectratio = :equal,
                              clims = slims,
                             levels = slevels,
                              xlims = xlims,
                              ylims = ylims,
                             xlabel = "x (km)",
                             ylabel = "y (km)")
                             
    plot(uxy_plot, speed_plot, size=(1100, 500), title = ["u(t="*string(round(t/day, digits=1))*" day)" "speed"])

    iter == iterations[end] && close(file)
end

gif(anim, "double_gyre.gif", fps = 12) # hide
