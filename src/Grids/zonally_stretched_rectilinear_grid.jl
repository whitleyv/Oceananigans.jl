"""
    ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ, R, A} <: AbstractRectilinearGrid{FT, TX, TY, TZ}

A rectilinear grid with with constant grid spacings `Δy` and `Δz`, and
non-uniform or stretched zonal grid spacing `Δx` between cell centers and cell faces,
topology `{TX, TY, TZ}`, and coordinate ranges of type `R` (where a range can be used) and
`A` (where an array is needed).
"""
struct ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ, R, A} <: AbstractRectilinearGrid{FT, TX, TY, TZ}

    # Number of grid points in (x,y,z).
    Nx :: Int
    Ny :: Int
    Nz :: Int

    # Halo size in (x,y,z).
    Hx :: Int
    Hy :: Int
    Hz :: Int

    # Domain size [m].
    Lx :: FT
    Ly :: FT
    Lz :: FT

    # Grid spacing [m].
    Δxᶜᵃᵃ :: A
    Δxᶠᵃᵃ :: A
       Δy :: FT
       Δz :: FT

    # Range of coordinates at the centers of the cells.
    xᶜᵃᵃ :: A
    yᵃᶜᵃ :: R
    zᵃᵃᶜ :: R

    # Range of grid coordinates at the faces of the cells.
    # Note: there are Nx+1 faces in the x-dimension, Ny+1 in the y, and Nz+1 in the z.
    xᶠᵃᵃ :: A
    yᵃᶠᵃ :: R
    zᵃᵃᶠ :: R
end

# FJP: since this is for ShallowWater, should the default in z be Flat?
function ZonallyStretchedRectilinearGrid(FT=Float64; 
                                   architecture = CPU(),
                                           size, 
                                           xF, 
                                              y = nothing, 
                                              z = nothing,
                                         extent = nothing,
                                       topology = (Bounded, Periodic, Periodic),
                                           halo = nothing
                                       )

      TX, TY, TZ = validate_topology(topology)
      print(TX, TY, TZ, size, "\n")
            size = validate_size(TX, TY, TZ, size)
            halo = validate_halo(TX, TY, TZ, halo)
    Ly, Lz, y, z = validate_zonally_stretched_grid_yz(TY, TZ, FT, extent, y, z)

    # Unpacking
    Nx, Ny, Nz = size
    Hx, Hy, Hz = halo
             L = (xF[end], Ly, Lz)

    # Initialize zonally-stretched arrays on CPU
    Lx, xᶠᵃᵃ, xᶜᵃᵃ, Δxᶜᵃᵃ, Δxᶠᵃᵃ = generate_stretched_zonal_grid(FT, topology[3], Nx, Hx, xF)

    # Construct uniform horizontal grid
    # FJP: don't want to force the grid to be uniform in the vertical slice!
    Lyz, Nyz, Hyz, X₁ = (Ly, Lz), size[2:3], halo[2:3], (y[1], z[1])
    Δy = Ly/size[2]
    Δz = Lz/size[3]        # what if Flat?

    # Face-node limits in x, y, z
    yF₋ = y[1] - Hyz[1] * Δy
    zF₋ = z[1] - Hyz[2] * Δz
    yF₊ = yF₋  + Ly + (2*Hyz[1] - 1) * Δy 
    zF₊ = zF₋  + Lz + (2*Hyz[2] - 1) * Δz 

    YF₋ = yF₋, zF₋     # are these needed?
    YF₊ = yF₊, zF₊ 

    # Center-node limits in x, y, z
    yC₋ = yF₋ + Δy/2
    zC₋ = zF₋ + Δz/2
    yC₊ = yC₋ + Ly * Δy * (2*Hyz[1] - 1)
    zC₊ = zC₋ + Lz * Δz * (2*Hyz[2] - 1)

    YC₋ = yC₋, zC₋     # are these needed?
    YC₊ = yC₊, zC₊

    # Total length of Center and Face quantities
    TFx, TFy, TFz = total_length.(Face, topology, size, halo)
    TCx, TCy, TCz = total_length.(Center, topology, size, halo)

    # Include halo points in coordinate arrays
    Δxᶠᵃᵃ = OffsetArray(Δxᶠᵃᵃ, -Hx)
    Δxᶜᵃᵃ = OffsetArray(Δxᶜᵃᵃ, -Hx)

    yᵃᶠᵃ = range(yF₋, yF₊; length = TFy)
    zᵃᵃᶠ = range(zF₋, zF₊; length = TFz)

    yᵃᶜᵃ = range(yC₋, yC₊; length = TCy)
    zᵃᵃᶜ = range(zC₋, zC₊; length = TCz)

    xᶜᵃᵃ = OffsetArray(xᶜᵃᵃ,  -Hx)
    yᵃᶜᵃ = OffsetArray(yᵃᶜᵃ,  -Hy)
    zᵃᵃᶜ = OffsetArray(zᵃᵃᶜ,  -Hz)

    xᶠᵃᵃ = OffsetArray(xᶠᵃᵃ,  -Hx)
    yᵃᶠᵃ = OffsetArray(yᵃᶠᵃ,  -Hy)
    zᵃᵃᶠ = OffsetArray(zᵃᵃᶠ,  -Hz)

    # Needed for pressure solver solution to be divergence-free.
    # Will figure out why later...
    # FJP: Is this still needed?   Not for ShallowWater but others?
    Δxᶠᵃᵃ[Nz] = Δxᶠᵃᵃ[Nx-1]

    # Seems needed to avoid out-of-bounds error in viscous dissipation
    # operators wanting to access Δxᶠᵃᵃ[Nx+2].
    Δxᶠᵃᵃ = OffsetArray(cat(Δxᶠᵃᵃ[0], Δxᶠᵃᵃ..., Δxᶠᵃᵃ[Nx], dims=1), -Hx-1)

    # Convert to appropriate array type for arch
    xᶠᵃᵃ  = OffsetArray(arch_array(architecture,  xᶠᵃᵃ.parent),  xᶠᵃᵃ.offsets...)
    xᶜᵃᵃ  = OffsetArray(arch_array(architecture,  xᶜᵃᵃ.parent),  xᶜᵃᵃ.offsets...)
    Δxᶜᵃᵃ = OffsetArray(arch_array(architecture, Δxᶜᵃᵃ.parent), Δxᶜᵃᵃ.offsets...)
    Δxᶠᵃᵃ = OffsetArray(arch_array(architecture, Δxᶠᵃᵃ.parent), Δxᶠᵃᵃ.offsets...)
    
    return ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ, typeof(xᶠᵃᵃ), typeof(zᵃᵃᶠ)}(
        Nx, Ny, Nz, Hx, Hy, Hz, Lx, Ly, Lz, Δxᶜᵃᵃ, Δxᶠᵃᵃ , Δy, Δz, xᶜᵃᵃ, yᵃᶜᵃ, zᵃᵃᶜ, xᶠᵃᵃ, yᵃᶠᵃ, zᵃᵃᶠ)
end

#####
##### Vertically stretched grid utilities
#####

get_x_face(x::Function, i) = x(i)
get_x_face(x::AbstractVector, i) = x[i]

lower_exterior_Δxᶜᵃᵃ(x_topo,          xFi, Hx) = [xFi[end - Hx + i] - xFi[end - Hx + i - 1] for i = 1:Hx]
lower_exterior_Δxᶜᵃᵃ(::Type{Bounded}, xFi, Hx) = [xFi[2]  - xFi[1] for i = 1:Hx]

upper_exterior_Δxᶜᵃᵃ(x_topo,          xFi, Hx) = [xFi[i + 1] - xFi[i] for i = 1:Hx]
upper_exterior_Δxᶜᵃᵃ(::Type{Bounded}, xFi, Hx) = [xFi[end]   - xFi[end - 1] for i = 1:Hx]

function generate_stretched_zonal_grid(FT, x_topo, Nx, Hx, xF_generator)

    # Ensure correct type for xF and derived quantities
    interior_xF = zeros(FT, Nx+1)

    for i = 1:Nx+1
        interior_xF[i] = get_x_face(xF_generator, i)
    end

    Lx = interior_xF[Nx+1] - interior_xF[1]

    # Build halo regions
    ΔxF₋ = lower_exterior_Δxᶜᵃᵃ(x_topo, interior_xF, Hx)
    ΔxF₊ = lower_exterior_Δxᶜᵃᵃ(x_topo, interior_xF, Hx)

    x¹, xᴺ⁺¹ = interior_xF[1], interior_xF[Nx+1]

    xF₋ = [x¹   - sum(ΔxF₋[i:Hx]) for i = 1:Hx] # locations of faces in lower halo
    xF₊ = [xᴺ⁺¹ + ΔxF₊[i]         for i = 1:Hx] # locations of faces in width of top halo region

    # FJP: why not use OffsetArrays to put halos beyond physical limits?
    xF = vcat(xF₋, interior_xF, xF₊)

    # Build cell centers, cell center spacings, and cell interface spacings
    TCx = total_length(Center, x_topo, Nx, Hx)
     xC = [ (xF[i + 1] + xF[i]) / 2 for i = 1:TCx ]
    ΔxC = [  xC[i] - xC[i - 1]      for i = 2:TCx ]

    # Trim face locations for periodic domains
    TFx = total_length(Face, x_topo, Nx, Hx)
    # FJP: this cuts off the halo value at the end.  Seems wrong
    #xF = xF[1:TFx]

    ΔxF = [xF[i + 1] - xF[i] for i = 1:TFx-1]

    return Lx, xF, xC, ΔxF, ΔxC
end

# We cannot reconstruct a ZonallyStretchedRectilinearGrid without the xF_generator.   FJP: why?
# So the best we can do is tell the user what they should have done.
function with_halo(new_halo, old_grid::ZonallyStretchedRectilinearGrid)
    new_halo != halo_size(old_grid) &&
        @error "You need to construct your ZonallyStretchedRectilinearGrid with the keyword argument halo=$new_halo"
    return old_grid
end

@inline x_domain(grid::ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ} = domain(TX, grid.Nx, grid.xᶠᵃᵃ)
@inline y_domain(grid::ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ} = domain(TY, grid.Ny, grid.yᵃᶠᵃ)
@inline z_domain(grid::ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ} = domain(TZ, grid.Nz, grid.zᵃᵃᶠ)

short_show(grid::ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ} =
    "ZonallyStretchedRectilinearGrid{$FT, $TX, $TY, $TZ}(Nx=$(grid.Nx), Ny=$(grid.Ny), Nz=$(grid.Nz))"

function show(io::IO, g::ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ}
    Δx_min = minimum(view(g.Δxᶜᵃᵃ, 1:g.Nx))
    Δx_max = maximum(view(g.Δxᶜᵃᵃ, 1:g.Nx))
    print(io, "ZonallyStretchedRectilinearGrid{$FT, $TX, $TY, $TZ}\n",
              "                   domain: $(domain_string(g))\n",
              "                 topology: ", (TX, TY, TZ), '\n',
              "  resolution (Nx, Ny, Nz): ", (g.Nx, g.Ny, g.Nz), '\n',
              "   halo size (Hx, Hy, Hz): ", (g.Hx, g.Hy, g.Hz), '\n',
              "grid spacing (Δx, Δy, Δz): , [min=", Δx_min, ", max=", Δx_max,"])", g.Δy, ", ", g.Δz,)
end

Adapt.adapt_structure(to, grid::ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ}) where {FT, TX, TY, TZ} =
    ZonallyStretchedRectilinearGrid{FT, TX, TY, TZ, typeof(Adapt.adapt(to, grid.xᶠᵃᵃ)), typeof(grid.zᵃᵃᶠ)}(
        grid.Nx, grid.Ny, grid.Nz,
        grid.Hx, grid.Hy, grid.Hz,
        grid.Lx, grid.Ly, grid.Lz,
        Adapt.adapt(to, grid.Δxᶜᵃᵃ),
        Adapt.adapt(to, grid.Δxᶠᵃᵃ),
        grid.Δy, grid.Δz,
        Adapt.adapt(to, grid.xᶜᵃᵃ),
        Adapt.adapt(to, grid.xᶠᵃᵃ),
        grid.yᵃᶜᵃ, grid.zᵃᵃᶜ,
        grid.yᵃᶠᵃ, grid.zᵃᵃᶠ)

#####
##### Should merge with grid_utils.jl at some point
#####

@inline xnode(::Type{Center}, i, grid::ZonallyStretchedRectilinearGrid) = @inbounds grid.xᶜᵃᵃ[i]
@inline xnode(::Type{Face},   i, grid::ZonallyStretchedRectilinearGrid) = @inbounds grid.xᶠᵃᵃ[i]

@inline ynode(::Type{Center}, j, grid::ZonallyStretchedRectilinearGrid) = @inbounds grid.yᵃᶜᵃ[j]
@inline ynode(::Type{Face},   j, grid::ZonallyStretchedRectilinearGrid) = @inbounds grid.yᵃᶠᵃ[j]

@inline znode(::Type{Center}, k, grid::ZonallyStretchedRectilinearGrid) = @inbounds grid.zᵃᵃᶜ[k]
@inline znode(::Type{Face},   k, grid::ZonallyStretchedRectilinearGrid) = @inbounds grid.zᵃᵃᶠ[k]


all_x_nodes(::Type{Center}, grid::ZonallyStretchedRectilinearGrid) = grid.xᶜᵃᵃ
all_x_nodes(::Type{Face},   grid::ZonallyStretchedRectilinearGrid) = grid.xᶠᵃᵃ
all_y_nodes(::Type{Center}, grid::ZonallyStretchedRectilinearGrid) = grid.yᵃᶜᵃ
all_y_nodes(::Type{Face},   grid::ZonallyStretchedRectilinearGrid) = grid.yᵃᶠᵃ
all_z_nodes(::Type{Center}, grid::ZonallyStretchedRectilinearGrid) = grid.zᵃᵃᶜ
all_z_nodes(::Type{Face},   grid::ZonallyStretchedRectilinearGrid) = grid.zᵃᵃᶠ



#
# Get minima of grid
#

min_Δx(grid::ZonallyStretchedRectilinearGrid) = minimum(view(grid.Δxᶜᵃᵃ, 1:grid.Nx))
min_Δy(grid::ZonallyStretchedRectilinearGrid) = grid.Δy
min_Δz(grid::ZonallyStretchedRectilinearGrid) = grid.Δz

