import Adapt

using Oceananigans: short_show
using Oceananigans.Operators: interpolation_operator
using Oceananigans.Fields: assumed_field_location, show_location
using Oceananigans.Utils: tupleit

"""
    ContinuousForcing{X, Y, Z, P, F, D, I}

A callable object that implements a "continuous form" forcing function
on a field at the location `X, Y, Z` with optional parameters.
"""
struct ContinuousForcing{X, Y, Z, D, N, F, P, TD, Tℑ}
                          func :: F
                    parameters :: P
            field_dependencies :: TD
                             ℑ :: Tℑ

    # Non-public "temporary" constructor that stores func, parameters, and field_dependencies
    # for later regularization
    function ContinuousForcing(func, parameters, field_dependencies)
        field_dependencies = tupleit(field_dependencies)

        return new{Nothing, Nothing, Nothing,
                   field_dependencies,
                   length(field_dependencies),
                   typeof(func), 
                   typeof(parameters), 
                   typeof(field_dependencies),
                   Nothing}(func, parameters, field_dependencies, nothing)
    end

    # Non-public "final" constructor.
    function ContinuousForcing{X, Y, Z}(func, parameters=nothing, field_dependencies=(),
                                        field_dependencies_indices=(), field_dependencies_interp=()) where {X, Y, Z}
        return new{X, Y, Z,
                   field_dependencies,
                   length(field_dependencies),
                   typeof(func),
                   typeof(parameters),
                   typeof(field_dependencies),
                   typeof(g)}(func, parameters, field_dependencies, field_dependencies_interp)
    end
end

"""
    ContinuousForcing(func; parameters=nothing, field_dependencies=())

Construct a "continuous form" forcing with optional `parameters` and optional
`field_dependencies` on other fields in a model.

If neither `parameters` nor `field_dependencies` are provided, then `func` must be 
callable with the signature

    `func(x, y, z, t)`

where `x, y, z` are the east-west, north-south, and vertical spatial coordinates, and `t` is time.

If `field_dependencies` are provided, the signature of `func` must include them.
For example, if `field_dependencies=(:u, :S)` (and `parameters` are _not_ provided), then
`func` must be callable with the signature

    `func(x, y, z, t, u, S)`

where `u` is assumed to be the `u`-velocity component, and `S` is a tracer. Note that any field
which does not have the name `u`, `v`, or `w` is assumed to be a tracer and must be present
in `model.tracers`.

If `parameters` are provided, then the _last_ argument to `func` must be `parameters`.
For example, if `func` has no `field_dependencies` but does depend on `parameters`, then
it must be callable with the signature

    `func(x, y, z, t, parameters)`

With `field_dependencies=(:u, :v, :w, :c)` and `parameters`, then `func` must be
callable with the signature

    `func(x, y, z, t, u, v, w, c, parameters)`

"""
ContinuousForcing(func; parameters=nothing, field_dependencies=()) =
    ContinuousForcing(func, parameters, field_dependencies)

""" 
    regularize_forcing(forcing::ContinuousForcing, field_name, model_field_names)

Regularize `forcing::ContinuousForcing` for use during time-stepping in an `IncompressibleModel`.

To do this, we

    * Obtain the location `X, Y, Z` at which `ContinuousForcing` is applied
    * Use `field_dependencies` tuple of `Symbols` to infer the tuple `field_dependencies_indices`
      that maps `field_dependencies` to `model_fields` (for use on the GPU)
    * Diagnose the interpolation operators `g` needed to interpolate
      `field_dependencies` to `X, Y, Z`
"""
function regularize_forcing(forcing::ContinuousForcing, field_name, model_field_names)

    X, Y, Z = assumed_field_location(field_name)

    g = Tuple(interpolation_operator(assumed_field_location(name), (X, Y, Z))
                                      for name in forcing.field_dependencies)

    field_dependencies_indices = ntuple(length(forcing.field_dependencies)) do i
        name = forcing.field_dependencies[i]
        findfirst(isequal(name), model_field_names)
    end

    return ContinuousForcing{X, Y, Z}(forcing.func, forcing.parameters, forcing.field_dependencies,
                                      field_dependencies_indices, g)
end

#####
##### Functions for calling ContinuousForcing in a time-stepping kernel
#####

@inline field_arguments(i, j, k, grid, forcing::ContinuousForcing{X, Y, Z, D, 1}, model_fields) where {X, Y, Z, D} =
    @inbounds (forcing.ℑ[1](i, j, k, grid, getproperty(model_fields, D[1])),)

@inline field_arguments(i, j, k, grid, forcing::ContinuousForcing{X, Y, Z, D, 2}, model_fields) where {X, Y, Z, D} =
    @inbounds (forcing.ℑ[1](i, j, k, grid, getproperty(model_fields, D[1])),
               forcing.ℑ[2](i, j, k, grid, getproperty(model_fields, D[2])))

@inline field_arguments(i, j, k, grid, forcing::ContinuousForcing{X, Y, Z, D, N}, model_fields) where {X, Y, Z, D, N} =
    @inbounds ntuple(n -> forcing.ℑ[n](i, j, k, grid, getproperty(model_fields, D[n])), Val(N))

""" Returns the arguments that follow `x, y, z, t` in a `ContinuousForcing` object without parameters. """
@inline forcing_func_arguments(i, j, k, grid, ::Nothing, forcing, model_fields) =
    field_arguments(i, j, k, grid, forcing, model_fields)

@inline function forcing_func_arguments(i, j, k, grid, parameters, forcing, model_fields)

    field_args = field_arguments(i, j, k, grid, forcing, model_fields)

    return tuple(field_args..., parameters)
end

@inline function (forcing::ContinuousForcing{X, Y, Z, F})(i, j, k, grid, clock, model_fields) where {X, Y, Z, F}

    args = forcing_func_arguments(i, j, k, grid, forcing.parameters, forcing, model_fields)

    return @inbounds forcing.func(xnode(X, i, grid), ynode(Y, j, grid), znode(Z, k, grid), clock.time, args...)
end

"""Show the innards of a `ContinuousForcing` in the REPL."""
Base.show(io::IO, forcing::ContinuousForcing{X, Y, Z, P}) where {X, Y, Z, P} =
    print(io, "ContinuousForcing{$P} at ", show_location(X, Y, Z), '\n',
        "├── func: $(short_show(forcing.func))", '\n',
        "├── parameters: $(forcing.parameters)", '\n',
        "└── field dependencies: $(forcing.field_dependencies)")

"""Show the innards of an "non-regularized" `ContinuousForcing` in the REPL."""
Base.show(io::IO, forcing::ContinuousForcing{Nothing, Nothing, Nothing, P}) where P =
    print(io, "ContinuousForcing{$P}", '\n',
        "├── func: $(short_show(forcing.func))", '\n',
        "├── parameters: $(forcing.parameters)", '\n',
        "└── field dependencies: $(forcing.field_dependencies)")

Adapt.adapt_structure(to, forcing::ContinuousForcing{X, Y, Z}) where {X, Y, Z} =
    ContinuousForcing{X, Y, Z}(Adapt.adapt(to, forcing.func),
                               Adapt.adapt(to, forcing.parameters),
                               nothing,
                               nothing,
                               Adapt.adapt(to, forcing.ℑ))
