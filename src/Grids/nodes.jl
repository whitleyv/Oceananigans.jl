#####
##### << Nodes >>
#####

# Node by node
@inline xnode(::Type{Cell}, i, grid) = @inbounds grid.xC[i]
@inline xnode(::Type{Face}, i, grid) = @inbounds grid.xF[i]

@inline ynode(::Type{Cell}, j, grid) = @inbounds grid.yC[j]
@inline ynode(::Type{Face}, j, grid) = @inbounds grid.yF[j]

@inline znode(::Type{Cell}, k, grid) = @inbounds grid.zC[k]
@inline znode(::Type{Face}, k, grid) = @inbounds grid.zF[k]

# Convenience is king
@inline xC(i, grid) = xnode(Cell, i, grid)
@inline xF(i, grid) = xnode(Face, i, grid)

@inline yC(j, grid) = ynode(Cell, j, grid)
@inline yF(j, grid) = ynode(Face, j, grid)

@inline zC(k, grid) = znode(Cell, k, grid)
@inline zF(k, grid) = znode(Face, k, grid)

all_x_nodes(::Type{Cell}, grid) = grid.xC
all_x_nodes(::Type{Face}, grid) = grid.xF
all_y_nodes(::Type{Cell}, grid) = grid.yC
all_y_nodes(::Type{Face}, grid) = grid.yF
all_z_nodes(::Type{Cell}, grid) = grid.zC
all_z_nodes(::Type{Face}, grid) = grid.zF

"""
    xnodes(loc, grid, reshape=false)

Returns a view over the interior `loc=Cell` or `loc=Face` nodes
on `grid` in the x-direction. For `Bounded` directions,
`Face` nodes include the boundary points. `reshape=false` will
return a 1D array while `reshape=true` will return a 3D array
with size Nx×1×1.

See `znodes` for examples.
"""
function xnodes(loc, grid; reshape=false)

    x = view(all_x_nodes(loc, grid),
             interior_indices(loc, topology(grid, 1), grid.Nx))

    return reshape ? Base.reshape(x, length(x), 1, 1) : x
end

"""
    ynodes(loc, grid, reshape=false)

Returns a view over the interior `loc=Cell` or `loc=Face` nodes
on `grid` in the y-direction. For `Bounded` directions,
`Face` nodes include the boundary points. `reshape=false` will
return a 1D array while `reshape=true` will return a 3D array
with size 1×Ny×1.


See `znodes` for examples.
"""
function ynodes(loc, grid; reshape=false)

    y = view(all_y_nodes(loc, grid),
             interior_indices(loc, topology(grid, 2), grid.Ny))

    return reshape ? Base.reshape(y, 1, length(y), 1) : y
end

"""
    znodes(loc, grid, reshape=false)

Returns a view over the interior `loc=Cell` or `loc=Face` nodes
on `grid` in the z-direction. For `Bounded` directions,
`Face` nodes include the boundary points. `reshape=false` will
return a 1D array while `reshape=true` will return a 3D array
with size 1×1×Nz.


Examples
========

```jldoctest znodes
julia> using Oceananigans, Oceananigans.Grids

julia> horz_periodic_grid = RegularCartesianGrid(size=(3, 3, 3), extent=(2π, 2π, 1),
                                                 topology=(Periodic, Periodic, Bounded));

julia> zC = znodes(Cell, horz_periodic_grid)
3-element view(OffsetArray(::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, 0:4), 1:3) with eltype Float64:
 -0.8333333333333331
 -0.4999999999999999
 -0.16666666666666652
```

``` jldoctest znodes
julia> zF = znodes(Face, horz_periodic_grid)
4-element view(OffsetArray(::StepRangeLen{Float64,Base.TwicePrecision{Float64},Base.TwicePrecision{Float64}}, 0:5), 1:4) with eltype Float64:
 -1.0
 -0.6666666666666666
 -0.33333333333333337
 -4.44089209850063e-17
```
"""
function znodes(loc, grid; reshape=false)

    z = view(all_z_nodes(loc, grid),
             interior_indices(loc, topology(grid, 3), grid.Nz))

    return reshape ? Base.reshape(z, 1, 1, length(z)) : z
end

"""
    nodes(loc, grid; reshape=false)

Returns a 3-tuple of views over the interior nodes
at the locations in `loc` in `x, y, z`.

If `reshape=true`, the views are reshaped to 3D arrays
with non-singleton dimensions 1, 2, 3 for `x, y, z`, respectively.
These reshaped arrays can then be used in broadcast operations with 3D fields
or arrays.

See `xnodes`, `ynodes`, and `znodes`.
"""
function nodes(loc, grid::AbstractGrid; reshape=false)
    if reshape
        x, y, z = nodes(loc, grid; reshape=false)

        N = (length(x), length(y), length(z))

        x = Base.reshape(x, N[1], 1, 1)
        y = Base.reshape(y, 1, N[2], 1)
        z = Base.reshape(z, 1, 1, N[3])

        return (x, y, z)
    else
        return (xnodes(loc[1], grid),
                ynodes(loc[2], grid),
                znodes(loc[3], grid))
    end
end

