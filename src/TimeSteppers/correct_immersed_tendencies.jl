
correct_immersed_tendencies!(model, Δt, γⁿ, ζⁿ) =
    correct_immersed_tendencies!(model, model.immersed_boundary, Δt, γⁿ, ζⁿ)

correct_immersed_tendencies!(model, ::Nothing, Δt, γⁿ, ζⁿ) = nothing

function correct_immersed_tendencies!(model, immersed_boundary, Δt, γⁿ, ζⁿ)

    workgroup, worksize = work_layout(model.grid, :xyz)

    barrier = Event(device(model.architecture))

    correct_immersed_tendencies_kernel! = _correct_immersed_tendencies!(device(model.architecture), workgroup, worksize)

    correct_tendencies_event =
        correct_immersed_tendencies_kernel!(model.timestepper.Gⁿ,
                                            model.grid,
                                            immersed_boundary,
                                            model.timestepper.G⁻,
                                            model.velocities,
                                            Δt, γⁿ, ζⁿ,
                                            dependencies=barrier)

    wait(device(model.architecture), correct_tendencies_event)

    return nothing
end

@kernel function _correct_immersed_tendencies!(Gⁿ, grid::AbstractGrid{FT}, immersed, G⁻, velocities, Δt, γⁿ, ζⁿ) where FT
    i, j, k = @index(Global, NTuple)

    x = xnode(Cell, i, grid)
    y = ynode(Cell, j, grid)
    z = znode(Cell, k, grid)

    @inbounds begin
        Gⁿ.u[i, j, k] = ifelse(immersed(x, y, z),
                               - (velocities.u[i, j, k] + ζⁿ * Δt * G⁻.u[i, j, k]) / (γⁿ * Δt),
                               Gⁿ.u[i, j, k])

        Gⁿ.v[i, j, k] = ifelse(immersed(x, y, z),
                               - (velocities.v[i, j, k] + ζⁿ * Δt * G⁻.v[i, j, k]) / (γⁿ * Δt),
                               Gⁿ.v[i, j, k])

        Gⁿ.w[i, j, k] = ifelse(immersed(x, y, z),
                               - (velocities.w[i, j, k] + ζⁿ * Δt * G⁻.w[i, j, k]) / (γⁿ * Δt),
                               Gⁿ.w[i, j, k])
    end
end

