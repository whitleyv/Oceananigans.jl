function instantiate_uniform_stokes_drift()
    ∂t_uˢ(z, t) = exp(z/20) * cos(t)
    ∂t_vˢ(z, t) = exp(z/20) * cos(t)
    ∂z_uˢ(z, t) = exp(z/20) * cos(t)
    ∂z_vˢ(z, t) = exp(z/20) * cos(t)

    surface_waves = SurfaceWaves.UniformStokesDrift(∂t_uˢ=∂t_uˢ, ∂t_vˢ=∂t_vˢ,
                                                    ∂z_uˢ=∂z_uˢ, ∂z_vˢ=∂z_vˢ)
                                                    
    return true
end

envelope(x, y) = exp(-(x^2 + y^2) / 100)

function instantiate_stokes_drift()

    ∂t_uˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)
    ∂t_vˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)
    ∂t_wˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)

    ∂y_uˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)
    ∂z_uˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)

    ∂x_vˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)
    ∂z_vˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)

    ∂x_wˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)
    ∂y_wˢ(x, y, z, t) = envelope(x, y) * exp(z/20) * cos(t)

    surface_waves = SurfaceWaves.StokesDrift(
                                             ∂t_uˢ=∂t_uˢ, ∂t_vˢ=∂t_vˢ, ∂t_wˢ=∂t_wˢ,
                                             ∂z_uˢ=∂z_uˢ, ∂y_uˢ=∂y_uˢ, 
                                             ∂z_vˢ=∂z_vˢ, ∂x_vˢ=∂x_vˢ, 
                                             ∂x_wˢ=∂x_wˢ, ∂y_wˢ=∂y_wˢ, 
                                            )
                                                    
    return true
end



@testset "Surface waves" begin
    @info "Testing surface waves..."

    @testset "Surface waves" begin
        @test instantiate_uniform_stokes_drift()
        @test instantiate_stokes_drift()
    end
end
