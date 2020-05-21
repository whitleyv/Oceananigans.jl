# # Diurnal boundary layer with biogeochemistry
#
# This example simulates the diurnal variation of a boundary layer that is heated
# during the day, and cooled at night. Because warm ocean water tends to rise to the surface,
# heating applied to the boundary layer during the day acts to "stabilize" the fluid, producing
# density stratification with warm ocean water on top, and cold ocean water below, and
# supressing turbulence. At night, however, radiation out of the top of boundary layer
# acts to cool the surface and leads to convection. On top of all of this, we add a bit
# of forcing by atmospheric winds to sustain turbulence during the day time.
#
# In addition to these physical effects, we implement a simple model for a carbon cycle
# in the ocean surface boundary layer. We thus demonstrate
#
#   * how to implement time-varying internal forcing and a time-varying surface boundary condition
#   * how to implement a system of passive, reacting tracers
#
# In addition to `Oceananigans.jl` we need `Plots` for plotting, `Random` for
# generating random initial conditions, and `Printf` for printing progress messages.
# We also need `Oceananigans.OutputWriters` and `Oceananigans.Diagnostics` to access
# some nice features for writing output data to disk.

using Random, Statistics, Printf, Plots, JLD2

using Oceananigans, Oceananigans.OutputWriters, Oceananigans.Diagnostics, Oceananigans.Utils,
      Oceananigans.BoundaryConditions, Oceananigans.Grids, Oceananigans.Forcing

# ## The domain
#
# We use a square, two-dimensional domain with vertical and horizontal
# resolution

Nx = 256  # Number of grid points in x
Nz = 192  # Number of grid points in z
Δx = 0.5  # Grid spacing in x (m)
Δz = 0.5  # Grid spacing in z (m)

# We can then make the grid,

grid = RegularCartesianGrid(size = (Nx, 1, Nz), extent = (Δx*Nx, 1, Δz*Nz))

# which by default is horizontally-periodic and bounded in the vertical.

# ## Buoyancy and equation of state
#
# We use the default gravitational acceleration and thermal expansion associated
# with SeawaterBuoyancy and LinearEquationOfState. We prescribe salinity to be constant.

buoyancy = SeawaterBuoyancy(equation_of_state=LinearEquationOfState(), constant_salinity=35.0)

# ## Boundary conditions
#
# ### Bottom boundary condition
#
# The bottom boundary condition imposes a constant stable temperature stratification
# with buoyancy gradient

N² = 1e-5 # s⁻²

# We derive the temperature gradient from the buoyancy gradient using gravitational
# acceleraiton and the thermal expansion coefficient, 

α = buoyancy.equation_of_state.α
g = buoyancy.gravitational_acceleration

∂T∂z = N² / (α * g)

# ### Top boundary condition
#
# At the top of the domain, we impose an outgoing temperature flux to model 
# the effects of radiative cooling at night. The outgoing temperature flux 
# thus varies on a diurnal cycle, and is neutral during the day.

peak_outgoing_radiation = 300 # Watts / m²

# The temperature flux associated with outgoing radiation of heat depends on 
# the reference density and heat capacity,

reference_density = 1035 # kg m⁻³
    heat_capacity = 3991 # J / kg / ᵒC

peak_temperature_flux = peak_outgoing_radiation / (reference_density * heat_capacity)

# For reference, we also calculate the buoyancy flux at the top due to cooling.

@show Qᵇ = α * g * peak_temperature_flux

# The diurnal variation of outgoing radiation is modeled with a clipped cosine function,

@inline outgoing_flux(x, y, t, p) = max(0, p.peak * cos(2π * (t / p.day - 0.5)))

outgoing_flux_bc = TracerBoundaryCondition(Flux, :z, outgoing_flux,
                                           (day=day, peak=peak_temperature_flux))

# To summarize, the temperature boundary condition are:
#
# * Bottom: stable temperature gradient
# * Top: diurnally varying flux due to outgoing radiation at night

T_bcs = TracerBoundaryConditions(grid, bottom = BoundaryCondition(Gradient, ∂T∂z),
                                          top = outgoing_flux_bc)

# ## Interior forcing due to solar heating
#
# Next, we specify interior forcing due to solar heating. We set the decay scale of
# the forcing

insolation_decay_scale = 20 # m

# and the solar insolation at the surface, 

surface_solar_insolation = 300 # Watts / m²

# which determines the surface temperature flux due to solar insolation,

surface_solar_temperature_flux = surface_solar_insolation / (reference_density * heat_capacity)

# We model the diurnal variations in solar heating with a clipped cosine,

@inline solar_flux_divergence(z, t, Qᴵ, λ, day) = Qᴵ / λ * exp(z / λ) * cos(2π * (t / day - 0.0))

@inline diurnal_solar_flux_divergence(x, y, z, t, p) =
    max(0, solar_flux_divergence(z, t, p.surface_flux, p.decay_scale, p.day))

# Finally, the temperature forcing functions are wrapped in a `SimpleForcing` object,

interior_heating = SimpleForcing(diurnal_solar_flux_divergence,
                                 parameters = (surface_flux = surface_solar_temperature_flux,
                                                decay_scale = insolation_decay_scale,
                                                        day = day))

# ## Model instantiation

model = IncompressibleModel(
                                   architecture = CPU(),
                                           grid = grid,
                                       coriolis = nothing,
                                        tracers = (:T,),
                                       buoyancy = buoyancy,
                                        closure = ConstantIsotropicDiffusivity(ν=1e-4, κ=1e-4),
                            boundary_conditions = (T=T_bcs,),
                                        forcing = ModelForcing(T=interior_heating),
                           )
nothing # hide

## Random noise concentrated at the top
Ξ(z) = randn() * exp(z / (8 * Δz))

## Temperature initial condition: a stable density tradient with random noise superposed.
T₀(x, y, z) = 20 + ∂T∂z * z + ∂T∂z * model.grid.Lz * 1e-4 * Ξ(z)

set!(model, T=T₀)

# ## Running the simulation
#
# To run the simulation, we instantiate a `TimeStepWizard` to ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 0.2.

wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=5.0)
nothing # hide

# A diagnostic that returns the maximum absolute value of `w` by calling
# `wmax(model)`:

wmax = FieldMaximum(abs, model.velocities.w)
nothing # hide

# # Running the simulation
#
# We set up a function to print a progress message, and then run the simulation.

function progress_message(simulation)
    model = simulation.model

    @printf("i: %04d, t: %s, Δt: %s, wmax = %.1e ms⁻¹\n",
            model.clock.iteration, prettytime(model.clock.time), prettytime(wizard.Δt),
            wmax(model))

    return nothing
end

simulation = Simulation(model, progress_frequency = 100,
                                        stop_time = 48hour,
                                               Δt = wizard, 
                                         progress = progress_message)
                      
## Add field writer
field_writer = JLD2OutputWriter(model, FieldOutputs(merge(model.velocities, model.tracers)),
                                interval = 1minute,
                                  prefix = "diurnal_boundary_layer_with_biogeochemistry",
                                   force = true)

simulation.output_writers[:fields] = field_writer

## Run
run!(simulation)

# # Visualize the results
#
# There's little point to running a fluids simulation if you don't make a pretty movie.
# Here we use Plots.jl to visualize the evolution of vertical velocity and temperature.

@info "Making an animation from the saved data..."

file = jldopen(field_writer.filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

anim = @animate for (i, iter) in enumerate(iterations)

    @printf "frame: %d, iter: %d \n" i iter

    ## Load the data and make a plot
    x, zw, zθ = xnodes(Cell, grid)[:], znodes(Face, grid)[:], znodes(Cell, grid)[:]

    w = file["timeseries/w/$iter"][2:end-1, 2, 2:end-1]
    θ = file["timeseries/T/$iter"][2:end-1, 2, 2:end-1]

    kwargs = (aspectratio=:equal, xlabel="x", ylabel="z", linecolor=:transparent,
              legend=false, xlims=(0, grid.Lx), zlims=(-grid.Lz, 0))

    wlim = maximum(abs, w) / 2
    wmax = maximum(abs, w)
    wlevels = range(-wlim, stop=wlim, length=10)
    wlevels = wmax > wlim ? vcat([-wmax], wlevels, [wmax]) : wlevels

    θlim⁺ = 20.1
    θlim⁻ = 19.8
    θmax = maximum(θ)
    θlevels = range(θlim⁻, stop=θlim⁺, length=10)
    θlevels = θmax > θlim⁺ ? vcat([0], θlevels, [θmax]) : vcat([0], θlevels)

    w_plot = heatmap(x, zw, w'; c=:balance, clims=(-wlim, wlim),  kwargs... ) #levels=wlevels, kwargs...)
    θ_plot = heatmap(x, zθ, θ'; c=:thermal, clims=(θlim⁻, θlim⁺), kwargs... ) #levels=θlevels, kwargs...)

    plot(w_plot, θ_plot, layout=(1, 2), size=(2000, 600))
end

close(file)

mp4(anim, "diurnal_boundary_layer_with_biogeochemistry.mp4", fps = 15) # hide
