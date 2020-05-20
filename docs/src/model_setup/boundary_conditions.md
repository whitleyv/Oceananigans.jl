# Boundary conditions

Boundary conditions are conditions that tracers and velocity components satisfy
at boundaries. Boudnary conditions are specified by users before constructing a
model, such as as `IncompressibleModel`.

## Overview

A boundary condition is applied to each field at the "right" and "left" ends of each dimension.
We use the cardinal directions to distinguish between the end points in each dimension; in other words:

1. x has west and east boundaries
2. y has south and north boundaries
3. z has bottom and top boundaries

Each field therefore has 6 boundary conditions, each of which may be specified.
Boudnary conditions can be specified as numeric constants, arrays, or functions that
depend on space, time, and other variables in the problem.

See [Numerical implementation of boundary conditions](@ref numerical_bcs) for more details.

## Examples

We start off with a few examples that build in complexity:

1. [Constant value on a boundary](@ref)
2. [Constant flux across a boundary](@ref)
3. [Constant gradient on a boundary](@ref)
4. [Flux that is randomly varying in space but constant in time](@ref)
5. [Flux that is Gaussian in space and sinusoidal in time](@ref)

#### Constant value on a boundary

To specify that a field has the value "20" on the boundary, we write

```@example
using Oceananigans # hide

constant_T_bc = ValueBoundaryCondition(20)
```

#### Constant flux across a boundary

For example, to specify a constant surface stress on a boundary,

```@example
using Oceananigans # hide

reference_density = 1027  # Reference density [kg/m³]
surface_stress = 0.08  # Wind stress [N/m²]

wind_stress_bc = FluxBoundaryCondition(surface_stress / reference_density)
```

This can also be used to represent heating, cooling, evaporation, or the flux of 
chemical species across a boundary.

#### Constant gradient on a boundary

```@example
using Oceananigans # hide

buoyancy_gradient = 1e-5 # Buoyancy gradient, the square of the buoyancy frequency [s⁻²]

bottom_buoyancy_bc = GradientBoundaryCondition(buoyancy_gradient)
```

#### Random, spatially-varying, constant-in-time flux

For example, to specify a spatially-varying cooling using an array,

```@example
using Oceananigans # hide

Nx = Ny = 16  # Number of grid points.

reference_density = 1027  # Reference density [kg/m³]
heat_capacity = 4000  # Heat capacity of water at constant pressure [J/kg/K]

Q  = randn(Nx, Ny) ./ (reference_density * heat_capacity)

white_noise_T_bc = FluxBoundaryCondition(Q)
```

When running on the GPU, `Q` must be converted to a `CuArray`.

#### Gaussian in space, sinusoidal in time flux

```@example
using Oceananigans, Oceananigans.Grids # hide

@inline Q(i, j, grid, clock, state) = @inbounds exp(-(xC(i, grid)^2 + yC(j, grid)^2)) * sin(2π * model.clock.t)

localized_heating_bc = FluxBoundaryCondition(Q)
```

This may also be written

```@example
using Oceananigans, Oceananigans.Grids # hide

@inline Q(x, y, t) = @inbounds exp(-(x^2 + y^2)) * sin(2π * t)

localized_heating_bc = TracerBoundaryCondition(Flux, :z, Q)
```

Constructors with the same syntax as `TracerBoundaryCondition`, for use with `BoundaryFunction`s,
but for velocity components, are:

1. `UVelocityBoundaryCondition`
2. `VVelocityBoundaryCondition`
3. `WVelocityBoundaryCondition`

!!! info "Performance of functions in boundary conditions"
    For performance reasons, you should define all functions used in boundary conditions as inline functions via the
    `@inline` macro. If any arrays are accessed within the function, disabling bounds-checking with `@inbounds` will
    also speed things up.

## Types of boundary conditions

The _type_ of boundary condition that may be specified depends on the "topology" of the domain.
For more information about specifying the `topology` of a domain, see [Specifying the grid's topology](@ref).

A `Periodic` direction, for example, only accepts `Periodic` boundary conditions (technically,
such a domain is unbounded and therefore does not have boundary conditions.) 
In a `Bounded` direction, tracer fields and tangential velocity components accept three boundary conditions types:

1. [`Value`](@ref Value)
2. [`Gradient`](@ref Gradient)
3. [`Flux`](@ref Flux)

which determine the `Value` of the field on the boundary, the `Gradient` of the field on the boundary,
or specify the `Flux` of the field across the boundary.
For example, the velocity components `u` and `v` are tangent to the top boundary (in the `z`-direction),
and therfore accept one of the three boundary conditions above when `z` is `Bounded`.

Velocity components that are normal to a boundary (for example, the vertical velocity `w` with respect to the
`z` direction) accept only one boundary condition:

4. [`No-penetration`](@ref NoPenetration)

Finally, directions that are `Periodic` have the boundary condition

5. [`Periodic`](@ref Periodic)

Neither `NoPenetration` nor `Periodic` boundary conditions are set by the user. Instead, these boundary conditions
are inferred during model construction from `topology(grid)`.

## Default boundary conditions

As mentioned above, `Periodic` boundary conditions apply in periodic directions, while `NoPenetration` boundary
conditions are applied to a velocity component normal to a `Bounded` direction.

In addition to these, a zero flux or "no flux" boundary condition is applied to tracers and tangential velocity 
components in `Bounded` directions by default.

## Boundary condition structures

Oceananigans uses a hierarchical structure to expressing boundary conditions.

1. A [`BoundaryCondition`](@ref) is associated with every field, dimension, and endpoint.
2. Boundary conditions specifying the condition at the left and right endpoints are
   grouped into [`CoordinateBoundaryConditions`](@ref).
3. A set of three `CoordinateBoundaryConditions` specifying the boundary conditions along the x, y, and z dimensions
   for a single field are grouped into a [`FieldBoundaryConditions`](@ref) named tuple.
4. A set of `FieldBoundaryConditions`, one for each field, are grouped together into a named tuple and passed to the
   `Model` constructor.

Boundary conditions are defined at model construction time by passing a named tuple of `FieldBoundaryConditions`
specifying boundary conditions on every field: velocities ($u$, $v$, $w$) and all tracers.

See the sections below for more details. The examples and verification experiments also provide examples for setting up
many different kinds of boundary conditions.

### Boundary condition functions

Oceananigans supports boundary condition numeric constants, arrays, and functions.
The function signature for a z boundary condition is

```
z_bc(i, j, grid, clock, state)
```

where `i, j` is the grid index, `grid` is `model.grid`, `clock` is `model.clock`, which contains `model.clock.time`, 
and `model.clock.iteration`, and `state` contains the `state.velocities` components and `state.tracers`.
The signature is similar for x and y boundary conditions expect that `i, j` is replaced with `j, k` and `i, k` respectively.

## Specifying boundary conditions on a field

To, for example, create a set of horizontally periodic field boundary conditions

```@example
using Oceananigans # hide

topology = (Periodic, Periodic, Bounded)

grid = RegularCartesianGrid(size=(16, 16, 16), extent=(1, 1, 1), topology=topology)

T_bcs = TracerBoundaryConditions(grid,    top = ValueBoundaryCondition(20),
                                       bottom = GradientBoundaryCondition(0.01))
```
which will create a [`FieldBoundaryConditions`](@ref) object for temperature T appropriate for horizontally periodic
model configurations where the x and y boundary conditions are all periodic.

## Specifying model boundary conditions

A named tuple of [`FieldBoundaryConditions`](@ref) objects must be passed to the Model constructor specifying boundary
conditions on all fields. To, for example, impose non-default boundary conditions on the u-velocity and temperature

```@example
using Oceananigans # hide

topology = (Periodic, Periodic, Bounded)

grid = RegularCartesianGrid(size=(16, 16, 16), extent=(1, 1, 1), topology=topology)

u_bcs = UVelocityBoundaryConditions(grid,    top = ValueBoundaryCondition(+0.1),
                                          bottom = ValueBoundaryCondition(-0.1))

T_bcs = TracerBoundaryConditions(grid,    top = ValueBoundaryCondition(20),
                                       bottom = GradientBoundaryCondition(0.01))

model = IncompressibleModel(grid=grid, boundary_conditions=(u=u_bcs, T=T_bcs))

nothing # hide
```
