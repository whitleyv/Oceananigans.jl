#####
##### Convinience functions and utils
#####

Base.length(loc, topo, N) = N
Base.length(::Type{Face}, ::Type{Bounded}, N) = N+1
Base.length(::Type{Nothing}, topo, N) = 1

function Base.size(loc, grid::AbstractGrid)
    N = (grid.Nx, grid.Ny, grid.Nz)
    return Tuple(length(loc[d], topology(grid, d), N[d]) for d in 1:3)
end

Base.size(loc, grid, d) = size(loc, grid)[d]

"""
    size(loc, grid)

Returns the size of a field at `loc` on `grid`, not including halos.
This is a 3-tuple of integers corresponding to the number of interior nodes
of `f` along `x, y, z`.
"""
@inline size(loc, grid) = (length(loc[1], topology(grid, 1), grid.Nx),
                           length(loc[2], topology(grid, 2), grid.Ny),
                           length(loc[3], topology(grid, 3), grid.Nz))

total_size(a) = size(a) # fallback

"""
    total_size(loc, grid)

Returns the "total" size of a field at `loc` on `grid`.
This is a 3-tuple of integers corresponding to the number of grid points
contained by `f` along `x, y, z`.
"""
@inline total_size(loc, grid) = (total_length(loc[1], topology(grid, 1), grid.Nx, grid.Hx),
                                 total_length(loc[2], topology(grid, 2), grid.Ny, grid.Hy),
                                 total_length(loc[3], topology(grid, 3), grid.Nz, grid.Hz))

"""
    total_extent(topology, H, Δ, L)

Returns the total extent, including halo regions, of constant-spaced
`Periodic` and `Flat` dimensions with number of halo points `H`,
constant grid spacing `Δ`, and interior extent `L`.
"""
@inline total_extent(topology, H, Δ, L) = L + (2H - 1) * Δ

"""
    total_extent(::Type{Bounded}, H, Δ, L)

Returns the total extent of, including halo regions, of constant-spaced
`Bounded` and `Flat` dimensions with number of halo points `H`,
constant grid spacing `Δ`, and interior extent `L`.
"""
@inline total_extent(::Type{Bounded}, H, Δ, L) = L + 2H * Δ

"""
    total_length(loc, topo, N, H=0)

Returns the total length (number of nodes), including halo points, of a field
located at `Cell` centers along a grid dimension of length `N` and with halo points `H`.
"""
@inline total_length(loc, topo, N, H=0) = N + 2H

"""
    total_length(::Type{Face}, ::Type{Bounded}, N, H=0)

Returns the total length, including halo points, of a field located at
cell `Face`s along a grid dimension of length `N` and with halo points `H`.
"""
@inline total_length(::Type{Face}, ::Type{Bounded}, N, H=0) = N + 1 + 2H

"""
    total_length(::Type{Nothing}, topo, N, H=0)

Returns 1, which is the 'length' of a field along a reduced dimension.
"""
@inline total_length(::Type{Nothing}, topo, N, H=0) = 1

# Grid domains
@inline domain(topo, N, ξ) = ξ[1], ξ[N+1]
@inline domain(::Type{Flat}, N, ξ) = ξ[1], ξ[1]

@inline x_domain(grid) = domain(topology(grid, 1), grid.Nx, grid.xF)
@inline y_domain(grid) = domain(topology(grid, 2), grid.Ny, grid.yF)
@inline z_domain(grid) = domain(topology(grid, 3), grid.Nz, grid.zF)

unpack_grid(grid) = grid.Nx, grid.Ny, grid.Nz, grid.Lx, grid.Ly, grid.Lz

flatten_halo(TX, TY, TZ, halo) = Tuple(T === Flat ? 0 : halo[i] for (i, T) in enumerate((TX, TY, TZ)))
flatten_size(TX, TY, TZ, halo) = Tuple(T === Flat ? 0 : halo[i] for (i, T) in enumerate((TX, TY, TZ)))

"""
    pop_flat_elements(tup, topo)

Returns a new tuple that contains the elements of `tup`,
except for those elements corresponding to the `Flat` directions
in `topo`.
"""
function pop_flat_elements(tup, topo)
    new_tup = []
    for i = 1:3
        topo[i] != Flat && push!(new_tup, tup[i])
    end
    return Tuple(new_tup)
end
