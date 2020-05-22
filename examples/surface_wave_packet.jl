# # Propagating surface wave packet

using Oceananigans, Oceananigans.Grids, Oceananigans.SurfaceWaves, Plots, Printf

Nx = 256 # resolution
Lx = 100  # domain extent
Lz = 1   # domain extent

struct Envelope{T}
                 a :: T
                 k :: T
                 g :: T
    group_velocity :: T
             scale :: T
                x₀ :: T
end

function Envelope(; a, k, g, scale, x₀=0)
    group_velocity = 1/2 * sqrt(g / k)
    return Envelope{Float64}(a, k, g, group_velocity, scale, x₀)
end

@inline    a(x, t, e) = e.a * exp(- (x - e.x₀ - e.group_velocity * t)^2 / (2 * e.scale^2)) 
@inline ∂t_a(x, t, e) = e.group_velocity * (x - e.x₀ - e.group_velocity * t) / e.scale^2 * a(x, t, e) 

struct PacketStokesDrift{D, A, T}
    envelope :: A
           k :: T
           g :: T

    function PacketStokesDrift{D}(env) where D
        return new{D, typeof(env), Float64}(env, env.k, env.g)
    end
end

PacketStokesDrift(env) = PacketStokesDrift{nothing}(env)

@inline (p::PacketStokesDrift{nothing})(x, y, z, t) =
    exp(2 * p.k * z) * a(x, t, p.envelope)^2 * p.k * sqrt(p.g * p.k)

@inline (p::PacketStokesDrift{:z})(x, y, z, t) =
    exp(2 * p.k * z) * 2 * a(x, t, p.envelope)^2 * p.k^2 * sqrt(p.g * p.k)

@inline (p::PacketStokesDrift{:t})(x, y, z, t) =
    exp(2 * p.k * z) * 2 * a(x, t, p.envelope) * ∂t_a(x, t, p.envelope) * p.k * sqrt(p.g * p.k)

env = Envelope(a=0.01, k=10, g=1, scale=1, x₀=-75)

   uˢ = PacketStokesDrift(env)
∂z_uˢ = PacketStokesDrift{:z}(env)
∂t_uˢ = PacketStokesDrift{:t}(env)

model = IncompressibleModel(         grid = RegularCartesianGrid(size=(Nx, 1, Nx), x=(-Lx/2, Lx/2), y=(0, 1), z=(-Lz, 0)),
                                  closure = ConstantIsotropicDiffusivity(ν=1e-6, κ=1e-6),
                                 coriolis = FPlane(f=0.1),
                                  tracers = :b, buoyancy = BuoyancyTracer(),
                                 # tracers = nothing, buoyancy = nothing,
                            surface_waves = SurfaceWaves.StokesDrift(∂z_uˢ=∂z_uˢ, ∂t_uˢ=∂t_uˢ)
                           )
nothing # hide

# We initialize the velocity and buoyancy fields
# with our internal wave initial condition.

uᵢ(x, y, z) = uˢ(x, y, z, 0)
bᵢ(x, y, z) = 0.1 * z

set!(model, u=uᵢ, b=bᵢ)
#set!(model, u=uᵢ)

# ## A wave packet on the loose
#
# Finally, we release the packet and watch it go!
advection_time_scale = env.scale / env.group_velocity

simulation = Simulation(model, Δt = 0.005 * advection_time_scale, stop_iteration = 0)

Uˢ = uˢ(0, 0, 0, 0)
Wˢ = Uˢ * env.k * env.scale

anim = @animate for i=0:1000

    x, zw, zu = xnodes(Cell, model.grid)[:], znodes(Face, model.grid)[:], znodes(Cell, model.grid)[:]

    u = interior(model.velocities.u)[:, 1, :]
    w = interior(model.velocities.w)[:, 1, :]

    packet = plot(x, ∂t_uˢ.(x, 0, 0, model.clock.time), legend=false)
    
    u_flow = contourf(x, zu, u',
                   title = @sprintf("t = %.2f", model.clock.time),
                  levels = range(-Uˢ, stop=Uˢ, length=10),
                   clims = (-Uˢ, Uˢ),
                  xlabel = "x",
                  ylabel = "z",
                   xlims = (-Lx/2, Lx/2),
                   ylims = (-Lz, 0),
                   color = :balance,
                  legend = false)
             #aspectratio = :equal)
             
    w_flow = contourf(x, zw, w',
                   title = @sprintf("t = %.2f", model.clock.time),
                  #levels = range(-Wˢ, stop=Wˢ, length=10),
                  # clims = (-Wˢ, Wˢ),
                  xlabel = "x",
                  ylabel = "z",
                   xlims = (-Lx/2, Lx/2),
                   ylims = (-Lz, 0),
                   color = :balance,
                  legend = false)
             #aspectratio = :equal)

            
    plot(packet, u_flow, w_flow, layout=(3, 1), size = (2000, 800))

    simulation.stop_iteration += 20

    run!(simulation)
end

mp4(anim, "weak_rotating_strong_stratified_flow_beneath_a_surface_wave_packet.mp4", fps = 15) # hide
