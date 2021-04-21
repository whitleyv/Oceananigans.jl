using Oceananigans.Grids: interior_parent_indices
using Statistics
using LinearAlgebra

import Oceananigans.Solvers: solve!, iterate!

function Statistics.norm(a::CubedSphereField)
    ii = interior_parent_indices(location(a, 1), topology(a.grid, 1), a.grid.Nx, a.grid.Hx)
    ji = interior_parent_indices(location(a, 2), topology(a.grid, 2), a.grid.Ny, a.grid.Hy)
    ki = interior_parent_indices(location(a, 3), topology(a.grid, 3), a.grid.Nz, a.grid.Hz)
    return sqrt(mapreduce(x -> x * x, +, view(parent(a), ii, ji, ki)))
end

function Statistics.dot(a::CubedSphereField, b::CubedSphereField)
    ii = interior_parent_indices(location(a, 1), topology(a.grid, 1), a.grid.Nx, a.grid.Hx)
    ji = interior_parent_indices(location(a, 2), topology(a.grid, 2), a.grid.Ny, a.grid.Hy)
    ki = interior_parent_indices(location(a, 3), topology(a.grid, 3), a.grid.Nz, a.grid.Hz)
    return mapreduce((x, y) -> x * y, +,
                     view(parent(a), ii, ji, ki),
                     view(parent(b), ii, ji, ki))
end

function solve!(x::CubedSphereField, solver::PreconditionedConjugateGradientSolver, b, args...)

    # Initialize
    solver.iteration = 0

    # q = A*x
    q = solver.linear_operator_product
    solver.linear_operation!(q, x, args...)

    # r = b - A*x
    # REWRITE FOR CUBEDSPHEREFIELD
    for (residual, b, q) in zip(faces(solver.residual), faces(b), faces(q)) # maybe
        parent(residual) .= parent(b) .- parent(q)
    end

    @debug "PreconditionedConjugateGradientSolver, |b|: $(norm(b))"
    @debug "PreconditionedConjugateGradientSolver, |A(x)|: $(norm(q))"

    while iterating(solver)
        iterate!(x, solver, b, args...)
    end

    fill_halo_regions!(x, solver.architecture)

    return nothing
end

function iterate!(x::CubedSphereField, solver, b, args...)
    r = solver.residual
    p = solver.search_direction
    q = solver.linear_operator_product

    @debug "PreconditionedConjugateGradientSolver $(solver.iteration), |r|: $(norm(r))"

    # Preconditioned:   z = P * r
    # Unpreconditioned: z = r
    z = maybe_precondition(solver.precondition!, solver.preconditioner_product, r, args...) 
    ρ = dot(z, r)

    @debug "PreconditionedConjugateGradientSolver $(solver.iteration), ρ: $ρ"
    @debug "PreconditionedConjugateGradientSolver $(solver.iteration), |z|: $(norm(z))"

    if solver.iteration == 0
        # REWRITE FOR CUBEDSPHEREFIELD
        parent(p) .= parent(z)
    else
        # REWRITE FOR CUBEDSPHEREFIELD
        β = ρ / solver.ρⁱ⁻¹
        parent(p) .= parent(z) .+ β .* parent(p)

        @debug "PreconditionedConjugateGradientSolver $(solver.iteration), β: $β"
    end

    # q = A * p
    solver.linear_operation!(q, p, args...)
    α = ρ / dot(p, q)

    @debug "PreconditionedConjugateGradientSolver $(solver.iteration), |q|: $(norm(q))"
    @debug "PreconditionedConjugateGradientSolver $(solver.iteration), α: $α"
        
    # REWRITE FOR CUBEDSPHEREFIELD
    parent(x) .+= α .* parent(p)
    parent(r) .-= α .* parent(q)

    solver.iteration += 1
    solver.ρⁱ⁻¹ = ρ

    return nothing
end

function iterating(solver)
    # End conditions
    solver.iteration >= solver.maximum_iterations && return false
    norm(solver.residual) <= solver.tolerance && return false
    return true
end

function Base.show(io::IO, solver::PreconditionedConjugateGradientSolver)
    print(io, "Oceananigans-compatible preconditioned conjugate gradient solver.\n")
    print(io, " Problem size = "  , size(solver.q), '\n')
    print(io, " Grid = "  , solver.grid)
    return nothing
end
