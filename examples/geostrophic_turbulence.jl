# # Two dimensional turbulence example
#
# In this example, we initialize a random velocity field and observe its turbulent decay
# in a two-dimensional domain. This example demonstrates:
#
#   * How to run a model with no tracers and no buoyancy model.
#   * How to use `AbstractOperations`.
#   * How to use `ComputedField`s to generate output.

# ## Install dependencies
#
# First let's make sure we have all required packages installed.

using Pkg
# pkg"add Oceananigans, Plots"

# ## Model setup

# We instantiate the model with an isotropic diffusivity. We use a grid with 128² points,
# a fifth-order advection scheme, third-order Runge-Kutta time-stepping,
# and a small isotropic viscosity.

using Oceananigans, Oceananigans.Advection

grid = RegularCartesianGrid(size=(128, 128, 1), extent=(2π, 2π, 2π),
                            topology=(Periodic, Bounded, Bounded))

model = IncompressibleModel(timestepper = :RungeKutta3,
                              advection = UpwindBiasedFifthOrder(),
                                   grid = grid,
                               coriolis = BetaPlane(f₀=2e-2, β=40), # f and beta are scaled
                               buoyancy = nothing,
                                tracers = nothing,
                                closure = IsotropicDiffusivity(ν=1e-5)
                           )

# ## Random initial conditions
#
# Our initial condition randomizes `model.velocities.u` and `model.velocities.v`.
# We ensure that both have zero mean for aesthetic reasons.

using Statistics

u₀ = rand(size(model.velocities.u)...)
u₀ .-= mean(u₀)
v₀ = rand(size(model.velocities.v)...)
v₀ .-= mean(v₀)

set!(model, u=u₀, v=v₀)

# ## Computing vorticity and speed

using Oceananigans.Fields, Oceananigans.AbstractOperations

# To make our equations prettier, we unpack `u`, `v`, and `w` from
# the `NamedTuple` model.velocities:
u, v, w = model.velocities

# Next we create two objects called `ComputedField`s that calculate
# _(i)_ vorticity that measures the rate at which the fluid rotates
# and is defined as
#
# ```math
# ω = ∂_x v - ∂_y u \, ,
# ```

ω = ∂x(v) - ∂y(u)

ω_field = ComputedField(ω)

# We also calculate _(ii)_ the _speed_ of the flow,
#
# ```math
# s = \sqrt{u^2 + v^2} \, .
# ```

s = sqrt(u^2 + v^2)

s_field = ComputedField(s)

# We'll pass these `ComputedField`s to an output writer below to calculate them during the simulation.
# Now we construct a simulation that prints out the iteration and model time as it runs.

progress(sim) = @info "Iteration: $(sim.model.clock.iteration), time: $(round(Int, sim.model.clock.time))"

wizard = TimeStepWizard(cfl=0.5, Δt=0.01)
simulation = Simulation(model, Δt=wizard, stop_time=100, iteration_interval=1, progress=progress)

# ## Output
#
# We set up an output writer for the simulation that saves the vorticity every 20 iterations.

using Oceananigans.OutputWriters

simulation.output_writers[:fields] =
    NetCDFOutputWriter(model, (u=model.velocities.u, ω=ω_field, s=s_field),
                       schedule = TimeInterval(1),
                       filepath = "geostrophic_turbulence.nc", mode="c")

# ## Running the simulation
#
# Pretty much just

run!(simulation)
close(simulation.output_writers[:fields])

# ## Visualizing the results
#
# We load the output and make a movie.

using GeoData, NCDatasets
using GeoData: GeoXDim, GeoYDim, GeoZDim
@dim xC GeoXDim "x"
@dim xF GeoXDim "x"
@dim yC GeoYDim "y"
@dim yF GeoYDim "y"
@dim zC GeoZDim "z"
@dim zF GeoZDim "z"


ds = NCDstack(simulation.output_writers[:fields].filepath)
u, ω, s = ds[:u], ds[:ω], ds[:s]
times = dims(u)[4]
Nt = length(times)

# animate the vorticity and fluid speed.

using Plots
using Oceananigans.Utils

@info "Making a neat movie of vorticity and speed..."

anim = @animate for n in 1:Nt

    @info "Plotting frame $n/$Nt..."

    ω_lim = 5.0
    ω_plot = contourf(clamp.(ω[Ti=n], -ω_lim, ω_lim),
                      color=:balance, clims=(-ω_lim, ω_lim), aspect_ratio=:auto,
                      title="Vorticity: $(prettytime(times[n]))",
                      aspectratio=1, linewidth=0)

    s_lim = 0.5
    s_plot = contourf(clamp.(s[Ti=n], 0, s_lim),
                      color=:thermal, clims=(0, s_lim), aspect_ratio=:auto,
                      title="Speed: $(prettytime(times[n]))",
                      aspectratio=1, linewidth=0)

    plot(ω_plot, s_plot, layout=(1, 2), size=(1600, 900))
end

gif(anim, "geostrophic_turbulence.gif", fps = 8) # hide

@info "Plotting Hovmoller plot of u..."

um = mean(u,dims=xF)
u_plot = heatmap(um.data[1,:,1,:])

gif(u_plot, "geostrophic_turbulence_u.gif") # hide

@info "Plotting wavenumber spectrum of energy..."

ke = 0.5*s.^2
ky = abs2.(fft(ke.data,2))
kxky = abs2.(fft(ky,1))
