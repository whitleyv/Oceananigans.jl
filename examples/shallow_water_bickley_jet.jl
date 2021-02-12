using Oceananigans.Architectures: CPU, architecture
using Oceananigans.Models: ShallowWaterModel
using Oceananigans.Grids: Periodic, Bounded, RegularCartesianGrid
using Oceananigans.Grids: xnodes, ynodes, interior
using Oceananigans.Simulations: Simulation, set!, run!, TimeStepWizard
using Oceananigans.Coriolis: FPlane
using Oceananigans.Advection: WENO5
using Oceananigans.OutputWriters: JLD2OutputWriter, IterationInterval, TimeInterval
using Oceananigans.Fields, Oceananigans.AbstractOperations

using Plots
using Printf
using JLD2
using LinearAlgebra
using IJulia

import Oceananigans.Utils: cell_advection_timescale

### Parameters

Lx = 2π       # Geometry
Ly = 20
Nx = 128
Ny = 128

f = 1            # Physics
g = 10

# Jet parameters
Uj = 1.0
Lj = 1.0
Δη = Uj * Lj * f /g

# Perturbation parameters
ϵ = 1e-1 # Perturbation amplitude
ℓ = 0.5  # Perturbation width
k = 1.0  # Perturbation wavenumber

grid = RegularCartesianGrid(
    size=(Nx, Ny, 1),
    x=(-Lx, Lx),
    y=(-Ly, Ly),
    z=( -1,  1),
    topology=(Periodic, Bounded, Bounded)
    )

#uh_bcs = BoundaryConditions(grid)

model = ShallowWaterModel(
    architecture=CPU(),
    advection=WENO5(),
    grid=grid,
    gravitational_acceleration=g,
    coriolis=FPlane(f=f)
    )

### Basic State
H₀(x, y, z) =   10.0
H(x, y, z)  =   H₀(x, y, z) - Δη * tanh(y)
U(x, y, z)  =         g / f * Δη * sech(y)^2
UH(x, y, z) = U(x, y, z) * H(x, y, z)
V(x, y, z)  =    0*y
VH(x, y, z) = V(x, y, z) * H(x, y, z)
Ω( x, y, z) =  2 * g / f * Δη * sech(y)^2 * tanh(y)

ψ(x, y, z, ℓ , k) = exp(- y / 2ℓ^2) * cos(k * x) * cos(k * y)

uhᵢ(x, y, z) = UH(x, y, z) + 0* ϵ * ψ(x, y, z, ℓ, k)
vhᵢ(x, y, z) = VH(x, y, z)
hᵢ(x, y, z) =  H(x, y, z)

set!(model, uh = uhᵢ , vh = vhᵢ, h = hᵢ)

u_op   = model.solution.uh / model.solution.h
v_op   = model.solution.vh / model.solution.h
η_op   = model.solution.h - H₀
ω_op   = @at (Center, Center, Center) ∂x(v_op) - ∂y(u_op)
ω_pert = @at (Center, Center, Center) ω_op - Ω

u_field = ComputedField(u_op)
v_field = ComputedField(v_op)
η_field = ComputedField(η_op)
ω_field = ComputedField(ω_op)
ω_pert  = ComputedField(ω_pert)

wizard = TimeStepWizard(cfl=0.2, Δt=1e-3, max_change=1.1, min_Δt=1e-4, max_Δt=1e-3)

function progress(simulation)
    @show simulation.Δt.Δt 
    @show simulation.model.clock.time
end

simulation = Simulation(model, Δt=wizard, progress=progress, iteration_interval=1, stop_time=200.0)

simulation.output_writers[:fields] =
    JLD2OutputWriter(
        model,
        (u = u_field, v = v_field, η = η_field, ω = ω_field, ωp = ω_pert),
        prefix = "Bickley_Jet",
        schedule=TimeInterval(0.5),
        force = true)


run!(simulation)

print(" after run \n ")

uhmax = maximum(abs, interior(model.solution.uh))
vhmax = maximum(abs, interior(model.solution.vh))
hmin  = minimum(abs, interior(model.solution.h))

print("uhmax = ", uhmax, " vhmax = ", vhmax, " hmin = ", hmin, "\n")

file = jldopen(simulation.output_writers[:fields].filepath)

kwargs = (
         xlabel = "y",
         ylabel = "x",
           fill = true,
         levels = 20,
      linewidth = 0,
          color = :balance,
       colorbar = true,
           ylim = (-Lx, Lx),
           xlim = (-Ly, Ly)
)

xc = xnodes(Center, grid)
yc = ynodes(Center, grid)

iterations = parse.(Int, keys(file["timeseries/t"]))

@info "Making a movie of the total and perturbation vorticity fields..."

anim = @animate for (i, iteration) in enumerate(iterations)

    @info "Plotting frame $i from iteration $iteration..."

    t = file["timeseries/t/$iteration"]
    η_snapshot  = file["timeseries/η/$iteration"][:, :, 1]
    ω_snapshot  = file["timeseries/ω/$iteration"][:, :, 1]
    ωp_snapshot = file["timeseries/ωp/$iteration"][:, :, 1]

    #η_plot = contour(yc, xc, η_snapshot, title="free-surface"; kwargs...)
    ω_plot  = contour(yc, xc, ω_snapshot,  title="Total Vorticity"; kwargs...)
    ωp_plot = contour(yc, xc, ωp_snapshot, title="Perturbation vorticity";  kwargs...)

    print("t = ", t, "\n")
    plot(
        ω_plot, 
        ωp_plot, 
        layout = (1,2), 
        size=(1200,500),
        title = @sprintf("Vorticities at t = %6.3f", t) 
        )
    
    #print(@sprintf("Vorticities at t = %6.3f", t),  "\n")
    print("Norm of perturbation = ", norm(ωp_snapshot), " with N = ", model.grid.Nx, "\n")

end

gif(anim, "swm_bickley_jet.gif", fps=8)

# To-Do-List
# 1) Fix halo regions
# 2) Fix how Δt is computed
# 3) Speed up

