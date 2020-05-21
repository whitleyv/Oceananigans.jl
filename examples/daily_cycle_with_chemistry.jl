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

using SeawaterPolynomials

# ## The domain
#
# Our domain is square, two-dimensional, horizontally-periodic, and deep enough to support 
# a diurnal boundary layer relatively unaffected by internal wave propagation into the 
# stratified interior below,

grid = RegularCartesianGrid(
                                size = (128, 1, 96), 
                              extent = (128, 1, 96), # meters
                            topology = (Periodic, Periodic, Bounded),
                           )

# ## Buoyancy and equation of state
#
# We use the default gravitational acceleration and the TEOS10 equation of state.
# We prescribe salinity to be constant, which means we don't need to calculate its
# distribution as we run the model.

buoyancy = SeawaterBuoyancy(
                            equation_of_state = SeawaterPolynomials.TEOS10EquationOfState(),
                            constant_salinity = 35.0 # psu
                           )

# ## Boundary conditions
#
# ### Bottom boundary condition
#
# The bottom boundary condition imposes a constant stable temperature stratification
# with buoyancy gradient

N² = 1e-5 # s⁻²

# We derive the temperature gradient from the buoyancy gradient using gravitational
# acceleration and the thermal expansion coefficient at 20ᵒC, and 35 psu at atmospheric
# pressure and thus z=0,

 α = SeawaterPolynomials.thermal_expansion(20, 35, 0, buoyancy.equation_of_state)
 g = buoyancy.gravitational_acceleration
ρᵣ = buoyancy.equation_of_state.reference_density

## Note that b = α g T / ρᵣ, and N² ≡ ∂b∂z
@show ∂T∂z = ρᵣ * N² / (α * g)

bottom_temperature_boundary_condition = BoundaryCondition(Gradient, ∂T∂z)

# ### Top boundary condition
#
# At the top of the domain, we impose an outgoing temperature flux to model 
# the effects of radiative cooling at night. The outgoing temperature flux 
# thus varies on a diurnal cycle, and is neutral during the day.

peak_outgoing_radiation = 300 # Watts / m²

# The temperature flux associated with outgoing radiation of heat depends on 
# the reference density and heat capacity. For conservative temperature
# described by TEOS10, the heat_capacity is

heat_capacity = 3991 # J / kg / ᵒC

# The reference density is also defined by TEOS10,

reference_density = buoyancy.equation_of_state.reference_density # kg m⁻³

# We can then calculate the outgoing temperature flux using the definition
#
#   $ \rm{internal_energy} = \rho * c_P * \Theta $,
#
# where $\Theta$ is Oceananigans' `T`.

peak_outgoing_flux = peak_outgoing_radiation / (reference_density * heat_capacity)

# The buoyancy flux at the top during peak outward radiation is, roughly,

Qᵇ = α * g * peak_outgoing_flux

@printf "The outward buoyancy flux at midnight is Qᵇ = %.2e m² s⁻³ \n" Qᵇ

# ## The diurnal and nocturnal cycles
#
# We arbitrarily define t = 0 as night noon.

""" Returns a function which is a positive half cosine during the day and 0 at night. """
@inline diurnal_cycle(t, day) = max(0, cos(2π * t / day))

""" Returns a function which is 0 during the day, and a (positive) half cosine at night. """
@inline nocturnal_cycle(t, day) = max(0, - cos(2π * t / day))

# The diurnal variation of outgoing radiation is thus

@inline outgoing_flux(x, y, t, p) = p.peak * nocturnal_cycle(t, p.day)

# and then wrapped inside a `BoundaryFunction` object,

surface_temperature_boundary_condition = 
    TracerBoundaryCondition(Flux, :z, outgoing_flux, (day=day, peak=peak_outgoing_flux))

# To summarize, the temperature boundary condition are:
#
# * Bottom at z = -Lz: stable temperature gradient
# * Top at z = 0: diurnally varying flux due to outgoing radiation at night

T_bcs = TracerBoundaryConditions(grid, bottom = bottom_temperature_boundary_condition,
                                          top = surface_temperature_boundary_condition)

# ## Interior forcing due to solar heating
#
# Next, we specify interior forcing due to solar heating. We set the decay scale of
# the forcing

light_attenuation_scale = 20 # m

# and the solar insolation at the surface, 

surface_solar_insolation = 600 # Watts / m²

# which determines the surface temperature flux due to solar insolation,

surface_solar_temperature_flux = surface_solar_insolation / (reference_density * heat_capacity)

# We model the diurnal variations in solar heating with a clipped cosine. The phase
# of the solar heating is such that t=0 corresponds to high noon.

@inline daylight(z, t, λ, day) = exp(z / λ) * diurnal_cycle(t, day)

@inline solar_flux_divergence(z, t, Qᴵ, λ, day) = Qᴵ / λ * daylight(z, t, λ, day)

@inline diurnal_solar_flux_divergence(x, y, z, t, p) =
    max(0, solar_flux_divergence(z, t, p.surface_flux, p.attenuation, p.day))

# Finally, the temperature forcing functions are wrapped in a `SimpleForcing` object,

interior_heating = SimpleForcing(diurnal_solar_flux_divergence,
                                 parameters = (surface_flux = surface_solar_temperature_flux,
                                                attenuation = light_attenuation_scale,
                                                        day = day))

# ## Biology and chemistry model

struct TracerReaction{F, P, R}
             forcing :: F
          parameters :: P
    reactive_tracers :: R

    function TracerReaction(forcing; reactive_tracers, parameters=nothing) 
        return new{typeof(forcing), 
                   typeof(parameters), 
                   typeof(reactive_tracers)}(forcing, parameters, reactive_tracers)
    end
end

@inline function (r::TracerReaction)(i, j, k, grid, clock, state)

    @inbounds reactive_tracers = Tuple(get_property(state.tracers, c)[i, j, k] 
                                       for c in r.reactive_tracers)

    return r.forcing(
                     xnode(Cell, i, grid),
                     ynode(Cell, j, grid),
                     znode(Cell, k, grid),
                     clock.time,
                     reactive_tracers...,
                     r.parameters
                    )
end

@inline γ(c, K) = c / (c + K)

@inline function biomass_rate(x, y, z, t, biomass, NH₄, NO₃, PO₄, p)

    γ_NH₄ = γ(NH₄, p.K_NH₄)
    γ_NO₃ = γ(NO₃, p.K_NO₃)
    γ_PO₃ = γ(PO₃, p.K_PO₃)

    return biomass * (
              p.max_growth_rate * daylight(z, t, p.attenuation, p.day) * min(γ_NH₄ + γ_NO₃, γ_PO₃)
            - p.respiration_rate
            - p.mortality_rate
           )
end

biochemical_parameters = (
                          attenuation = light_attenuation_scale,
                          day = day,
                          respiration_rate = 1,
                          mortality_rate = 1,
                          K_NH₄ = 1,
                          K_NO₃ = 1,
                          K_PO₃ = 1,
                         )

biomass_forcing = TracerReaction(
                                 biomass_rate,
                                       parameters = biochemical_parameters,
                                 reactive_tracers = (:biomass, :NH₄, :NO₃, :PO₃)
                                )

# ## Model instantiation
# 
# We specify the diffusivity that, through trial and error, was determined to produce 
# a solution "without too much grid-scale noise".
#
tracer_names = (
                :T,                  # temperature
                :biomass,            # phytoplankton biomass
                :inorganic_carbon,   # dissolved inorganic carbon
                :organic_carbon,     # dissolved organic carbon
                :particulate_carbon, # particular (not dissolved) carbon
               )

model = IncompressibleModel(
                                   architecture = CPU(),
                                           grid = grid,
                                       coriolis = nothing,
                                        tracers = :T, #tracer_names,
                                       buoyancy = buoyancy,
                                        closure = ConstantIsotropicDiffusivity(ν=1e-3, κ=1e-3),
                            boundary_conditions = (T=T_bcs,),
                                        forcing = ModelForcing(T=interior_heating),
                           )

# ## The initial condition
#
# Our initial condition is a stably stratified temperature gradient with a bit
# of Gaussian-distributed random noise superimposed,

T₀(x, y, z) = (
               20 + ∂T∂z * z # constant plus linear gradient
                  + ∂T∂z * grid.Lz * 1e-4 * randn() * exp(z / (8 * grid.Δz)) # noise
               )

set!(model, T=T₀)

# ## Build the Simulation
#
# To run the simulation, we instantiate a `TimeStepWizard` to ensure stable time-stepping
# with a Courant-Freidrichs-Lewy (CFL) number of 0.2.

wizard = TimeStepWizard(cfl=0.2, Δt=1.0, max_change=1.1, max_Δt=5.0)

# We also make a diagnostic that returns the maximum absolute value of `w` by calling
# `wmax(model)`, for convenient logging:

wmax = FieldMaximum(abs, model.velocities.w)

# We set up a function to print a progress message, and then run the simulation.

wall_clock = time_ns()

function print_message(simulation)
    model = simulation.model

    msg = @sprintf("i: %04d, t: %s, Δt: %s, wmax = %.1e ms⁻¹, wall time: %s\n", 
                   model.clock.iteration, 
                   prettytime(model.clock.time),
                   prettytime(wizard.Δt), 
                   wmax(model), 
                   prettytime(1e-9 * (time_ns()-wall_clock))
                  )

    @info msg

    return nothing
end

# We instantiate the simulation, setting the stop time to 48 hours and asking for
# a progress update every 100 iterations:

simulation = Simulation(model, progress_frequency = 100,
                                        stop_time = 48hour,
                                               Δt = wizard, 
                                         progress = print_message)
                      
# We add a field writer so we can make a movie from the output afterwards

field_writer = JLD2OutputWriter(model, FieldOutputs(merge(model.velocities, model.tracers)),
                                interval = 2minute,
                                  prefix = "diurnal_boundary_layer_with_biogeochemistry",
                                   force = true)

simulation.output_writers[:fields] = field_writer

# ## Run the simulation (this takes some time)

run!(simulation)

# ## Visualize the results
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

    θlim⁺ = 20.00
    θlim⁻ = 19.75
    θmax = maximum(θ)
    θlevels = range(θlim⁻, stop=θlim⁺, length=10)
    θlevels = θmax > θlim⁺ ? vcat([0], θlevels, [θmax]) : vcat([0], θlevels)

    w_plot = heatmap(x, zw, w'; c=:balance, clims=(-wlim, wlim),  kwargs... ) #levels=wlevels, kwargs...)
    θ_plot = heatmap(x, zθ, θ'; c=:thermal, clims=(θlim⁻, θlim⁺), kwargs... ) #levels=θlevels, kwargs...)

    plot(w_plot, θ_plot, layout=(1, 2), size=(2000, 800))
end

close(file)

mp4(anim, "diurnal_boundary_layer_with_biogeochemistry.mp4", fps = 15) # hide
