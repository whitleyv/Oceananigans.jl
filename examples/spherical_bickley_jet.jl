# # Geostrophic adjustment using Oceananigans.HydrostaticFreeSurfaceModel
#
# ## Install dependencies
#
# First let's make sure we have all required packages installed.

# ```julia
# using Pkg
# pkg"add Oceananigans, JLD2, Plots"
# ```

# ## A spherical domain
#
# We use a one-dimensional domain of geophysical proportions,

using Oceananigans
using Oceananigans.Grids: RegularLatitudeLongitudeGrid

grid = RegularLatitudeLongitudeGrid(size = (360, 160, 1), longitude = (-180, 180), latitude = (-80, 80), z = (-1, 0))

using Oceananigans.Coriolis: HydrostaticSphericalCoriolis

coriolis = HydrostaticSphericalCoriolis()

# ## Building a `HydrostaticFreeSurfaceModel`
#
# We use `grid` and `coriolis` to build a simple `HydrostaticFreeSurfaceModel`,

using Oceananigans.Models.HydrostaticFreeSurfaceModels: VectorInvariant
using Oceananigans.Models: HydrostaticFreeSurfaceModel
using Oceananigans.TurbulenceClosures: HorizontallyCurvilinearAnisotropicDiffusivity

closure = HorizontallyCurvilinearAnisotropicDiffusivity(νh=1, κh=1)

model = HydrostaticFreeSurfaceModel(grid = grid,
                                    momentum_advection = VectorInvariant(),
                                    tracers = (),
                                    buoyancy = nothing,
                                    coriolis = nothing,
                                    closure = closure)

# ## The Bickley jet on a sphere

Ψ(ϕ) = - tanh(ϕ)
U(ϕ) = sech(ϕ)^2

# A sinusoidal tracer
C(ϕ) = sind(ϕ)

# Slightly off-center vortical perturbations
ψ̃(λ, ϕ, ℓ, k) = exp(-(ϕ + ℓ/10)^2 / 2ℓ^2) * cosd(k * λ) * cosd(k * ϕ)

# Vortical velocity fields (ũ, ṽ) = (-∂_y, +∂_x) ψ̃
ũ(λ, ϕ, ℓ, k) = + ψ̃(λ, ϕ, ℓ, k) * (k * tan(k * ϕ) + ϕ / ℓ^2)
ṽ(λ, ϕ, ℓ, k) = - ψ̃(λ, ϕ, ℓ, k) * k * tan(k * λ)

# Parameters
ϵ = 0.1 # perturbation magnitude
ℓ = 4   # Gaussian width (degrees)
k = 4   # Sinusoidal wavenumber

# Total initial conditions
uᵢ(λ, ϕ, z) = U(ϕ) + ϵ * ũ(λ, ϕ, ℓ, k)
vᵢ(λ, ϕ, z) = ϵ * ṽ(λ, ϕ, ℓ, k)
cᵢ(λ, ϕ, z) = C(ϕ, grid.Lϕ)

set!(model, u=uᵢ, v=vᵢ) #, c=cᵢ)

# ## Running a `Simulation`
#
# We pick a time-step that resolves the surface dynamics,

using Oceananigans.Utils: prettytime

g = model.free_surface.gravitational_acceleration

gravity_wave_speed = sqrt(g * grid.Lz) # hydrostatic (shallow water) gravity wave speed

wave_propagation_time_scale = 100e3 * model.grid.Δλ / gravity_wave_speed

simulation = Simulation(model, Δt = 0.05wave_propagation_time_scale, stop_time = 100wave_propagation_time_scale,
                        progress = s -> @info "Time = $(prettytime(s.model.clock.time)) / $(prettytime(s.stop_time))")

# ## Output
#
# We output the velocity field and free surface displacement,

output_fields = merge(model.velocities, (η=model.free_surface.η,))

using Oceananigans.OutputWriters: JLD2OutputWriter, TimeInterval

simulation.output_writers[:fields] = JLD2OutputWriter(model, output_fields,
                                                      schedule = TimeInterval(0.05wave_propagation_time_scale),
                                                      prefix = "spherical_bickley",
                                                      force = true)

run!(simulation)

# ## Visualizing the results

using JLD2, Printf, Oceananigans.Grids, GLMakie
using Oceananigans.Utils: hours

λ, ϕ, r = nodes(model.free_surface.η, reshape=true)

λ = λ .+ 180  # Convert to λ ∈ [0°, 360°]
ϕ = 90 .- ϕ   # Convert to ϕ ∈ [0°, 180°] (0° at north pole)

file = jldopen(simulation.output_writers[:fields].filepath)

iterations = parse.(Int, keys(file["timeseries/t"]))

iter = Node(0)
plot_title = @lift @sprintf("Oceananigans.jl on the sphere! Rossby splash: u, v, η @ time = %s", prettytime(file["timeseries/t/" * string($iter)]))
u = @lift file["timeseries/u/" * string($iter)][:, :, 1]
v = @lift file["timeseries/v/" * string($iter)][:, :, 1]
η = @lift file["timeseries/η/" * string($iter)][:, :, 1]

# Plot on the unit sphere to align with the spherical wireframe.
# Multiply by 1.01 so the η field is a bit above the wireframe.
x = @. 1.01 * cosd(λ) * sind(ϕ)
y = @. 1.01 * sind(λ) * sind(ϕ)
z = @. 1.01 * cosd(ϕ) * λ ./ λ

x = x[:, :, 1]
y = y[:, :, 1]
z = z[:, :, 1]

fig = Figure(resolution = (1920, 1080))

clims = [(-0.003, 0.003), (-0.003, 0.003), (-0.01, 0.01)]

for (n, var) in enumerate([u, v, η])
    ax = fig[1, n] = LScene(fig, title="$n")
    wireframe!(ax, Sphere(Point3f0(0), 1f0), show_axis=false)
    surface!(ax, x, y, z, color=var, colormap=:balance, colorrange=clims[n])
    rotate_cam!(ax.scene, (2π/3, 0, 0))
    zoom!(ax.scene, (0, 0, 0), 5, false)
end

supertitle = fig[0, :] = Label(fig, plot_title, textsize=30)

record(fig, "rossby_splash.mp4", iterations, framerate=30) do i
    @info "Animating iteration $i/$(iterations[end])..."
    iter[] = i
end
