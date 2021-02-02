struct ImplicitFreeSurface end

function pressure_correct_velocities!(model::IncompressibleModel{T, E, <:ImplicitFreeSurface}, Δt) where {T, E}
    return nothing
end
