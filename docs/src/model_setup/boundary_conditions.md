# Boundary conditions

Boundary conditions are conditions that tracers and velocity components satisfy
at boundaries. Boudnary conditions are specified by users before constructing a
model, such as as `IncompressibleModel`.

## Overview

A boundary condition is applied to each field at the "right" and "left" ends of each dimension.
We use the cardinal directions to distinguish between the end points in each dimension:

1. x has west and east boundaries,
2. y has south and north boundaries,
3. z has bottom and top boundaries.

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

value_bc = BoundaryCondition(Value, 20)
```

This is useful for setting a constant temperature on a boundary, for example.

#### Constant flux across a boundary

To specify a constant stress on a boundary,

```@example
using Oceananigans # hide

reference_density = 1027  # kg / m³
surface_stress = 0.08     # N / m²

stress_bc = BoundaryCondition(Flux, surface_stress / reference_density)
```

This is a model for the transfer of atmospheric momentum to the ocean.
A constant flux boundary condition can be used to represent heating, cooling,
evaporation, or the flux of chemical species across a boundary.

#### Constant gradient on a boundary

To specify a constant buoyancy gradient, write

```@example
using Oceananigans # hide

buoyancy_gradient = 1e-5 # aka the square of the buoyancy frequency [s⁻²]

buoyancy_bc = BoundaryCondition(Gradient, buoyancy_gradient)
```

#### Flux that is randomly varying in space but constant in time

To specify spatially-varying cooling using an array, write

```@example
using Oceananigans # hide

Nx = Ny = 16  # Number of grid points.

reference_density = 1027  # kg / m³
heat_capacity = 4000  # (at constant pressure) [J / kg / K]

cooling = randn(Nx, Ny) ./ (reference_density * heat_capacity)

noisy_bc = BoundaryCondition(Flux, cooling)
```

When running on the GPU, `heat_flux` must be converted to a `CuArray`.

#### Flux that is Gaussian in space but sinusoidal in time

```@example
using Oceananigans, Oceananigans.Grids # hide

@inline Q(i, j, grid, clock, state) = @inbounds exp(-(xC(i, grid)^2 + yC(j, grid)^2)) * sin(2π * model.clock.t)

localized_cooling_bc = BoundaryCondition(Flux, Q)
```

`Oceananigans` also provides a wrapper for using simpler function signatures that are closer to the
mathematical specification of the problem,

```@example
using Oceananigans, Oceananigans.BoundaryConditions # hide

@inline Q(x, y, t) = @inbounds exp(-(x^2 + y^2)) * sin(2π * t)

localized_cooling_bc = TracerBoundaryCondition(Flux, :z, Q)
```

Constructors with the same syntax as `TracerBoundaryCondition` for use with `BoundaryFunction`s,
but for velocity components, are:

1. `UVelocityBoundaryCondition`
2. `VVelocityBoundaryCondition`
3. `WVelocityBoundaryCondition`

!!! info "Performance of functions in boundary conditions"
    For performance reasons, you should define all functions used in boundary conditions as inline functions via the
    `@inline` macro. If any arrays are accessed within the function, disabling bounds-checking with `@inbounds` will
    also speed things up.

### Types of boundary conditions

The _type_ of boundary condition that may be specified depends on the "topology" of the domain.
For more information about specifying the `topology` of a domain, see [Specifying the grid's topology](@ref).

A `Periodic` direction, for example, only accepts `Periodic` boundary conditions (technically,
such a domain is unbounded and therefore does not have boundary conditions.) 
In a `Bounded` direction, tracer fields and tangential velocity components accept three boundary conditions types:

* [`Value`](@ref Value)
* [`Gradient`](@ref Gradient)
* [`Flux`](@ref Flux)

which determine the `Value` of the field on the boundary, the `Gradient` of the field on the boundary,
or specify the `Flux` of the field across the boundary.
For example, the velocity components `u` and `v` are tangent to the top boundary (in the `z`-direction),
and therfore accept one of the three boundary conditions above when `z` is `Bounded`.

Velocity components that are normal to a boundary (for example, the vertical velocity `w` with respect to the
`z` direction) accept only one boundary condition:

* [`NoPenetration`](@ref NoPenetration)

Finally, directions that are `Periodic` have the boundary condition

* [`Periodic`](@ref Periodic)

Neither `NoPenetration` nor `Periodic` boundary conditions are set by the user. Instead, these boundary conditions
are inferred during model construction from `topology(grid)`.

### Boundary condition functions

Oceananigans supports boundary conditions that are numeric constants, two-dimensinoal arrays, and functions.
The function signature for a z boundary condition is

```julia
func(i, j, grid, clock, state)
```

where `i, j` are the `x, y` indices, `grid` is `model.grid`, `clock` is `model.clock`, which contains `model.clock.time`, 
and `model.clock.iteration`, and `state` contains the `state.velocities` components and `state.tracers`.
The expected signature for `x` and `y` boundary conditions is identical, except that `i, j`
correspond to `j, k` (`y, z` for `x` boundaries) and `i, k` (`x, z` for `y` boundaries).

## Building boundary conditions for a field

The above examples show how to create individual boundary condition that apply at a particular boundary.
Once all of the non-default boundary conditions for a field are determined, we must create an object
to collects all six of the field's boundary conditions. Under the hood, this object is called `FieldBoundaryConditions`.

To create a set of horizontally periodic boundary conditions with non-default top and bottom boundary conditions for
a tracer field, write

```@example
using Oceananigans # hide

grid = RegularCartesianGrid(topology = (Periodic, Periodic, Bounded),
                                size = (16, 16, 16),
                              extent = (1, 1, 1))

temperature_boundary_conditions = TracerBoundaryConditions(grid, bottom = BoundaryCondition(Gradient, 0.01),
                                                                    top = BoundaryCondition(Value, 20))
```

The display above shows a [`FieldBoundaryConditions`](@ref) object.
The `Periodic` horizontal boundary conditions are inferred from `topology(grid)`.

## Specifying boundary conditions for an `IncompressibleModel`

Boundary conditions on the fields of a model are specified through the `boundary_conditions` keyword
argument in the `IncompressibleModel` constructor. For example,

```@example
using Oceananigans # hide

grid = RegularCartesianGrid(topology = (Periodic, Periodic, Bounded),
                                size = (16, 16, 16),
                              extent = (1, 1, 1))

u_bcs = UVelocityBoundaryConditions(grid, bottom = BoundaryCondition(Value, -0.1),
                                             top = BoundaryCondition(Value,  0.1))

T_bcs = TracerBoundaryConditions(grid, bottom = BoundaryCondition(Gradient, 0.01),
                                          top = BoundaryCondition(Value, 20))

model = IncompressibleModel(grid=grid, tracers=(:T, :S), boundary_conditions=(u=u_bcs, T=T_bcs))
```

imposes non-default boundary condition on the `u` velocity component and tracer `T`.

## Default boundary conditions

### Periodic directions

`Periodic` boundary conditions apply in periodic directions.

```@example
using Oceananigans # hide

grid = RegularCartesianGrid(topology = (Periodic, Periodic, Bounded),
                                size = (16, 16, 16),
                              extent = (1, 1, 1))
                            
w_bcs = WVelocityBoundaryConditions(grid)

@show w_bcs.x.left
```

### Bounded directions

#### Normal component of the velocity field 

A `NoPenetration` boundary condition is applied to a velocity component normal to a `Bounded` direction.

```@example
using Oceananigans # hide

grid = RegularCartesianGrid(topology = (Periodic, Periodic, Bounded),
                                size = (16, 16, 16),
                              extent = (1, 1, 1))
                            
w_bcs = WVelocityBoundaryConditions(grid)

@show w_bcs.z.top
```

#### Tracers and tangential components of the velocity field

A zero flux or "no flux" boundary condition is applied to tracers and tangential velocity 
components in `Bounded` directions by default.

For example,

```@example
using Oceananigans # hide

grid = RegularCartesianGrid(topology = (Periodic, Periodic, Bounded),
                                size = (16, 16, 16),
                              extent = (1, 1, 1))
                            
u_bcs = UVelocityBoundaryConditions(grid)
```

Defaults are applied within the model constructor.
Adapting the first example in the README, we find,

```@example
using Oceananigans # hide

model = IncompressibleModel(grid=RegularCartesianGrid(size=(16, 16, 16), extent=(1, 1, 1)))

@show model.velocities.u.boundary_conditions.z.bottom
```
