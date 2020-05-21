module SurfaceWaves

export
    ∂t_uˢ,
    ∂t_vˢ,
    ∂t_wˢ,
    x_curl_Uˢ_cross_U,
    y_curl_Uˢ_cross_U,
    z_curl_Uˢ_cross_U

using Oceananigans.Grids: AbstractGrid

using Oceananigans.Fields
using Oceananigans.Operators

"""
    abstract type AbstractStokesDrift end

Parent type for parameter structs for Stokes drift fields
associated with surface waves
"""
abstract type AbstractStokesDrift end

#####
##### Functions for "no surface waves"
#####

@inline ∂t_uˢ(i, j, k, grid, ::Nothing, time) = zero(eltype(grid))
@inline ∂t_vˢ(i, j, k, grid, ::Nothing, time) = zero(eltype(grid))
@inline ∂t_wˢ(i, j, k, grid, ::Nothing, time) = zero(eltype(grid))

@inline x_curl_Uˢ_cross_U(i, j, k, grid, ::Nothing, U, time) = zero(eltype(grid))
@inline y_curl_Uˢ_cross_U(i, j, k, grid, ::Nothing, U, time) = zero(eltype(grid))
@inline z_curl_Uˢ_cross_U(i, j, k, grid, ::Nothing, U, time) = zero(eltype(grid))

#####
##### Uniform surface waves
#####

"""
    UniformStokesDrift{UZ, VZ, UT, VT} <: AbstractStokesDrift

Parameter struct for Stokes drift fields associated with surface waves.
"""
struct UniformStokesDrift{UZ, VZ, UT, VT} <: AbstractStokesDrift
    ∂z_uˢ :: UZ
    ∂z_vˢ :: VZ
    ∂t_uˢ :: UT
    ∂t_vˢ :: VT
end

addzero(args...) = 0

"""
    UniformStokesDrift(; ∂z_uˢ=addzero, ∂z_vˢ=addzero,
                                    ∂t_uˢ=addzero, ∂t_vˢ=addzero)

Construct a set of functions that describes the Stokes drift field beneath
a uniform surface gravity wave field.
"""
UniformStokesDrift(; ∂z_uˢ=addzero, ∂z_vˢ=addzero, ∂t_uˢ=addzero, ∂t_vˢ=addzero) =
    UniformStokesDrift(∂z_uˢ, ∂z_vˢ, ∂t_uˢ, ∂t_vˢ)

const USD = UniformStokesDrift

@inline ∂t_uˢ(i, j, k, grid, sw::USD, time) = sw.∂t_uˢ(znode(Cell, k, grid), time)
@inline ∂t_vˢ(i, j, k, grid, sw::USD, time) = sw.∂t_vˢ(znode(Cell, k, grid), time)
@inline ∂t_wˢ(i, j, k, grid, sw::USD, time) = zero(eltype(grid))

@inline x_curl_Uˢ_cross_U(i, j, k, grid, sw::USD, U, time) =
    @inbounds ℑxzᶠᵃᶜ(i, j, k, grid, U.w) * sw.∂z_uˢ(znode(Cell, k, grid), time)

@inline y_curl_Uˢ_cross_U(i, j, k, grid, sw::USD, U, time) =
    @inbounds ℑyzᵃᶠᶜ(i, j, k, grid, U.w) * sw.∂z_vˢ(znode(Cell, k, grid), time)

@inline z_curl_Uˢ_cross_U(i, j, k, grid, sw::USD, U, time) = @inbounds begin (
    - ℑxzᶜᵃᶠ(i, j, k, grid, U.u) * sw.∂z_uˢ(znode(Face, k, grid), time)
    - ℑyzᵃᶜᶠ(i, j, k, grid, U.v) * sw.∂z_vˢ(znode(Face, k, grid), time) )
end

#####
##### Three-dimensional, horizontally-modulated surface waves
#####

"""
    UniformStokesDrift{UZ, VZ, UT, VT} <: AbstractStokesDrift

Parameter struct for Stokes drift fields associated with surface waves.
"""
struct StokesDrift{
                   UY, 
                   UZ, 
                   VX, 
                   VZ, 
                   WX, 
                   WY, 
                   UT, 
                   VT
                   WT
                   } <: AbstractStokesDrift

    ∂y_uˢ :: UY
    ∂z_uˢ :: UZ

    ∂x_vˢ :: VX
    ∂z_vˢ :: VZ

    ∂x_wˢ :: WX
    ∂y_wˢ :: WY

    ∂t_uˢ :: UT
    ∂t_vˢ :: VT
    ∂t_wˢ :: WT

end

addzero(args...) = 0

"""
    StokesDrift(; 
                ∂y_uˢ=addzero,
                ∂z_uˢ=addzero,
                ∂x_vˢ=addzero,
                ∂z_vˢ=addzero,
                ∂x_wˢ=addzero,
                ∂y_wˢ=addzero,
                ∂t_uˢ=addzero, 
                ∂t_vˢ=addzero
                ∂t_wˢ=addzero
                )

Construct a set of functions that describes the Stokes drift field beneath
a uniform surface gravity wave field.
"""
function StokesDrift(; 
                       ∂y_uˢ=addzero,
                       ∂z_uˢ=addzero,

                       ∂x_vˢ=addzero,
                       ∂z_vˢ=addzero,

                       ∂x_wˢ=addzero,
                       ∂y_wˢ=addzero,

                       ∂t_uˢ=addzero,
                       ∂t_vˢ=addzero
                       ∂t_wˢ=addzero
                    )

    return StokesDrift(
                       ∂y_uˢ,
                       ∂z_uˢ,
                       ∂x_vˢ,
                       ∂z_vˢ,
                       ∂x_wˢ,
                       ∂y_wˢ,
                       ∂t_uˢ,
                       ∂t_vˢ
                       ∂t_wˢ
                      )

const SD = StokesDrift

@inline ∂t_uˢ(i, j, k, grid, sw::SD, time) = sw.∂t_uˢ(xnode(Face, i, grid), ynode(Cell, j, grid), znode(Cell, k, grid), time)
@inline ∂t_vˢ(i, j, k, grid, sw::SD, time) = sw.∂t_vˢ(xnode(Cell, i, grid), ynode(Face, j, grid), znode(Cell, k, grid), time)
@inline ∂t_wˢ(i, j, k, grid, sw::SD, time) = sw.∂t_wˢ(xnode(Cell, i, grid), ynode(Cell, j, grid), znode(Face, k, grid), time)

@inline function x_curl_Uˢ_cross_U(i, j, k, grid, sw::SD, U, time)

    x = xnode(Face, i, grid)
    y = ynode(Cell, j, grid)
    z = znode(Cell, k, grid)

    return (
              ℑxzᶠᵃᶜ(i, j, k, grid, U.w) * (sw.∂z_uˢ(x, y, z, time) - sw.∂x_wˢ(x, y, z, time))
            + ℑxyᶠᶜᵃ(i, j, k, grid, U.v) * (sw.∂y_uˢ(x, y, z, time) - sw.∂x_vˢ(x, y, z, time))
           )
end

@inline function y_curl_Uˢ_cross_U(i, j, k, grid, sw::SD, U, time)

    x = xnode(Cell, i, grid)
    y = ynode(Face, j, grid)
    z = znode(Cell, k, grid)

    return (
              ℑyzᵃᶠᶜ(i, j, k, grid, U.w) * (sw.∂z_vˢ(x, y, z, time) - sw.∂x_wˢ(x, y, z, time))
            + ℑxyᶜᶠᵃ(i, j, k, grid, U.u) * (sw.∂z_vˢ(x, y, z, time) - sw.∂x_wˢ(x, y, z, time))
           )
end

@inline function z_curl_Uˢ_cross_U(i, j, k, grid, sw::SD, U, time)

    x = xnode(Cell, i, grid)
    y = ynode(Cell, j, grid)
    z = znode(Face, k, grid)

    return (
              ℑxzᶜᵃᶠ(i, j, k, grid, U.u) * (sw.∂x_wˢ(x, y, z, time) - sw.∂z_uˢ(x, y, z, time))
            + ℑyzᵃᶜᶠ(i, j, k, grid, U.v) * (sw.∂y_wˢ(x, y, z, time) - sw.∂z_vˢ(x, y, z, time))
           )
end

end # module
