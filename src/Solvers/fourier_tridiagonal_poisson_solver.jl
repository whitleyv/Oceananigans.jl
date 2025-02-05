using Oceananigans.Operators: Δzᵃᵃᶜ, Δzᵃᵃᶠ

struct FourierTridiagonalPoissonSolver{A, G, B, R, S, β, T}
                  architecture :: A
                          grid :: G
    batched_tridiagonal_solver :: B
                   source_term :: R
                       storage :: S
                        buffer :: β
                    transforms :: T
end

@kernel function compute_diagonals!(D, grid, λx, λy)
    i, j = @index(Global, NTuple)
    Nz = grid.Nz

    # Using a homogeneous Neumann (zero Gradient) boundary condition:
    D[i, j, 1] = -1 / Δzᵃᵃᶠ(i, j, 2, grid) - Δzᵃᵃᶜ(i, j, 1, grid) * (λx[i] + λy[j])

    @unroll for k in 2:Nz-1
        D[i, j, k] = - (1 / Δzᵃᵃᶠ(i, j, k+1, grid) + 1 / Δzᵃᵃᶠ(i, j, k, grid)) - Δzᵃᵃᶜ(i, j, k, grid) * (λx[i] + λy[j])
    end

    D[i, j, Nz] = -1 / Δzᵃᵃᶠ(i, j, Nz, grid) - Δzᵃᵃᶜ(i, j, Nz, grid) * (λx[i] + λy[j])
end

function FourierTridiagonalPoissonSolver(arch, grid, planner_flag=FFTW.PATIENT)
    TX, TY, TZ = topology(grid)
    TZ != Bounded && error("FourierTridiagonalPoissonSolver can only be used with a Bounded z topology.")

    if grid isa VerticallyStretchedRectilinearGrid && any([T() isa Flat for T in (TX, TY)])
        @warn "FourierTridiagonalPoissonSolver is probably wrong for topologies that contain " *
              "Flat dimensions."
    end

    Nx, Ny, Nz = size(grid)

    # Compute discrete Poisson eigenvalues
    λx = poisson_eigenvalues(grid.Nx, grid.Lx, 1, TX())
    λy = poisson_eigenvalues(grid.Ny, grid.Ly, 2, TY())

    λx = arch_array(arch, λx)
    λy = arch_array(arch, λy)

    # Plan required transforms for x and y
    sol_storage = arch_array(arch, zeros(complex(eltype(grid)), size(grid)...))
    transforms = plan_transforms(arch, grid, sol_storage, planner_flag)

    # Lower and upper diagonals are the same
    lower_diagonal = CUDA.@allowscalar [1 / Δzᵃᵃᶠ(1, 1, k, grid) for k in 2:Nz]
    lower_diagonal = arch_array(arch, lower_diagonal)
    upper_diagonal = lower_diagonal

    # Compute diagonal coefficients for each grid point
    diagonal = arch_array(arch, zeros(Nx, Ny, Nz))
    event = launch!(arch, grid, :xy, compute_diagonals!, diagonal, grid, λx, λy,
                    dependencies=Event(device(arch)))
    wait(device(arch), event)

    # Set up batched tridiagonal solver
    btsolver = BatchedTridiagonalSolver(arch, grid;
                                         lower_diagonal = lower_diagonal,
                                               diagonal = diagonal,
                                         upper_diagonal = upper_diagonal)

    # Need buffer for index permutations and transposes.
    buffer_needed = arch isa GPU && Bounded in (TX, TY) ? true : false
    buffer = buffer_needed ? similar(sol_storage) : nothing

    # Storage space for right hand side of Poisson equation
    rhs = arch_array(arch, zeros(complex(eltype(grid)), size(grid)...))

    return FourierTridiagonalPoissonSolver(arch, grid, btsolver, rhs, sol_storage, buffer, transforms)
end

function solve_poisson_equation!(solver::FourierTridiagonalPoissonSolver)
    ϕ = solver.storage
    RHS = solver.source_term

    # Apply forward transforms in order
    [transform!(RHS, solver.buffer) for transform! in solver.transforms.forward]

    # Solve tridiagonal system of linear equations in z at every column.
    solve!(ϕ, solver.batched_tridiagonal_solver, RHS)

    # Apply backward transforms in order
    [transform!(ϕ, solver.buffer) for transform! in solver.transforms.backward]

    ϕ .= real.(ϕ)

    # Set the volume mean of the solution to be zero.
    # Solutions to Poisson's equation are only unique up to a constant (the global mean
    # of the solution), so we need to pick a constant. We choose the constant to be zero
    # so that the solution has zero-mean.
    ϕ .= ϕ .- mean(ϕ)

    return nothing
end

@kernel function calculate_pressure_source_term_fourier_tridiagonal_solver!(RHS, grid, Δt, U★)
    i, j, k = @index(Global, NTuple)

    @inbounds RHS[i, j, k] = Δzᵃᵃᶜ(i, j, k, grid) * divᶜᶜᶜ(i, j, k, grid, U★.u, U★.v, U★.w) / Δt
end

"""
    set_source_term!(solver, source_term)

Sets the source term in the discrete Poisson equation `solver`
to `source_term` by multiplying it by the vertical grid spacing at z cell centers.
"""
function set_source_term!(solver::FourierTridiagonalPoissonSolver, source_term)
    grid = solver.grid
    arch = solver.architecture
    solver.source_term .= source_term

    event = launch!(arch, grid, :xyz, multiply_by_Δzᵃᵃᶜ!, solver.source_term, grid, dependencies=Event(device(arch)))
    wait(device(arch), event)
                    
    return nothing
end

@kernel function multiply_by_Δzᵃᵃᶜ!(a, grid)
    i, j, k = @index(Global, NTuple)
    a[i, j, k] *= Δzᵃᵃᶜ(i, j, k, grid)
end
