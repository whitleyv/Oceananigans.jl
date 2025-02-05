#####
##### NetCDFOutputWriter tests
#####

function test_DateTime_netcdf_output(arch)
    grid = RegularRectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
    clock = Clock(time=DateTime(2021, 1, 1))
    model = IncompressibleModel(architecture=arch, grid=grid, clock=clock)

    Δt = 5days + 3hours + 44.123seconds
    simulation = Simulation(model, Δt=Δt, stop_time=DateTime(2021, 2, 1))

    filepath = "test_DateTime.nc"
    simulation.output_writers[:cal] = NetCDFOutputWriter(model, fields(model), filepath=filepath, schedule=IterationInterval(1))

    run!(simulation)

    ds = NCDataset(filepath)
    @test ds["time"].attrib["units"] == "seconds since 2000-01-01 00:00:00"

    Nt = length(ds["time"])
    @test Nt == 8

    for n in 1:Nt-1
        @test ds["time"][n] == DateTime(2021, 1, 1) + (n-1) * Millisecond(1000Δt)
    end

    @test ds["time"][Nt] == DateTime(2021, 2, 1)

    close(ds)
    rm(filepath)

    return nothing
end

function test_TimeDate_netcdf_output(arch)
    grid = RegularRectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
    clock = Clock(time=TimeDate(2021, 1, 1))
    model = IncompressibleModel(architecture=arch, grid=grid, clock=clock)

    Δt = 5days + 3hours + 44.123seconds
    simulation = Simulation(model, Δt=Δt, stop_time=TimeDate(2021, 2, 1))

    filepath = "test_TimeDate.nc"
    simulation.output_writers[:cal] = NetCDFOutputWriter(model, fields(model), filepath=filepath, schedule=IterationInterval(1))

    run!(simulation)

    ds = NCDataset(filepath)
    @test ds["time"].attrib["units"] == "seconds since 2000-01-01 00:00:00"

    Nt = length(ds["time"])
    @test Nt == 8

    for n in 1:Nt-1
        @test ds["time"][n] == DateTime(2021, 1, 1) + (n-1) * Millisecond(1000Δt)
    end

    @test ds["time"][Nt] == DateTime(2021, 2, 1)

    close(ds)
    rm(filepath)

    return nothing
end

function test_thermal_bubble_netcdf_output(arch)
    Nx, Ny, Nz = 16, 16, 16
    Lx, Ly, Lz = 100, 100, 100

    topo = (Periodic, Periodic, Bounded)
    grid = RegularRectilinearGrid(topology=topo, size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz))
    closure = IsotropicDiffusivity(ν=4e-2, κ=4e-2)
    model = IncompressibleModel(architecture=arch, grid=grid, closure=closure)
    simulation = Simulation(model, Δt=6, stop_iteration=10)

    # Add a cube-shaped warm temperature anomaly that takes up the middle 50%
    # of the domain volume.
    i1, i2 = round(Int, Nx/4), round(Int, 3Nx/4)
    j1, j2 = round(Int, Ny/4), round(Int, 3Ny/4)
    k1, k2 = round(Int, Nz/4), round(Int, 3Nz/4)
    CUDA.@allowscalar model.tracers.T.data[i1:i2, j1:j2, k1:k2] .+= 0.01

    outputs = Dict("v" => model.velocities.v,
                   "u" => model.velocities.u,
                   "w" => model.velocities.w,
                   "T" => model.tracers.T,
                   "S" => model.tracers.S)

    nc_filepath = "test_dump_$(typeof(arch)).nc"
    nc_writer = NetCDFOutputWriter(model, outputs, filepath=nc_filepath, schedule=IterationInterval(10), verbose=true)
    push!(simulation.output_writers, nc_writer)

    i_slice = 1:2:10
    j_slice = 13
    k_slice = 9:11
    field_slicer = FieldSlicer(i=i_slice, j=j_slice, k=k_slice)
    j_slice = j_slice:j_slice  # So we can correctly index with it for later tests.

    nc_sliced_filepath = "test_dump_sliced_$(typeof(arch)).nc"
    nc_sliced_writer = NetCDFOutputWriter(model, outputs, filepath=nc_sliced_filepath, schedule=IterationInterval(10),
                                          field_slicer=field_slicer, verbose=true)

    push!(simulation.output_writers, nc_sliced_writer)

    run!(simulation)

    ds3 = Dataset(nc_filepath)

    @test haskey(ds3.attrib, "date") && !isnothing(ds3.attrib["date"])
    @test haskey(ds3.attrib, "Julia") && !isnothing(ds3.attrib["Julia"])
    @test haskey(ds3.attrib, "Oceananigans") && !isnothing(ds3.attrib["Oceananigans"])
    @test haskey(ds3.attrib, "schedule") && ds3.attrib["schedule"] == "IterationInterval"
    @test haskey(ds3.attrib, "interval") && ds3.attrib["interval"] == 10
    @test haskey(ds3.attrib, "output iteration interval") && !isnothing(ds3.attrib["output iteration interval"])

    @test eltype(ds3["time"]) == eltype(model.clock.time)

    @test eltype(ds3["xC"]) == Float64
    @test eltype(ds3["xF"]) == Float64
    @test eltype(ds3["yC"]) == Float64
    @test eltype(ds3["yF"]) == Float64
    @test eltype(ds3["zC"]) == Float64
    @test eltype(ds3["zF"]) == Float64

    @test length(ds3["xC"]) == Nx
    @test length(ds3["yC"]) == Ny
    @test length(ds3["zC"]) == Nz
    @test length(ds3["xF"]) == Nx
    @test length(ds3["yF"]) == Ny
    @test length(ds3["zF"]) == Nz+1  # z is Bounded

    @test ds3["xC"][1] == grid.xC[1]
    @test ds3["xF"][1] == grid.xF[1]
    @test ds3["yC"][1] == grid.yC[1]
    @test ds3["yF"][1] == grid.yF[1]
    @test ds3["zC"][1] == grid.zC[1]
    @test ds3["zF"][1] == grid.zF[1]

    @test ds3["xC"][end] == grid.xC[Nx]
    @test ds3["xF"][end] == grid.xF[Nx]
    @test ds3["yC"][end] == grid.yC[Ny]
    @test ds3["yF"][end] == grid.yF[Ny]
    @test ds3["zC"][end] == grid.zC[Nz]
    @test ds3["zF"][end] == grid.zF[Nz+1]  # z is Bounded

    @test eltype(ds3["u"]) == Float32
    @test eltype(ds3["v"]) == Float32
    @test eltype(ds3["w"]) == Float32
    @test eltype(ds3["T"]) == Float32
    @test eltype(ds3["S"]) == Float32

    u = ds3["u"][:, :, :, end]
    v = ds3["v"][:, :, :, end]
    w = ds3["w"][:, :, :, end]
    T = ds3["T"][:, :, :, end]
    S = ds3["S"][:, :, :, end]

    close(ds3)

    @test all(u .≈ Array(interior(model.velocities.u)))
    @test all(v .≈ Array(interior(model.velocities.v)))
    @test all(w .≈ Array(interior(model.velocities.w)))
    @test all(T .≈ Array(interior(model.tracers.T)))
    @test all(S .≈ Array(interior(model.tracers.S)))

    ds2 = Dataset(nc_sliced_filepath)

    @test haskey(ds2.attrib, "date") && !isnothing(ds2.attrib["date"])
    @test haskey(ds2.attrib, "Julia") && !isnothing(ds2.attrib["Julia"])
    @test haskey(ds2.attrib, "Oceananigans") && !isnothing(ds2.attrib["Oceananigans"])
    @test haskey(ds2.attrib, "schedule") && ds2.attrib["schedule"] == "IterationInterval"
    @test haskey(ds2.attrib, "interval") && ds2.attrib["interval"] == 10
    @test haskey(ds2.attrib, "output iteration interval") && !isnothing(ds2.attrib["output iteration interval"])

    @test eltype(ds2["time"]) == eltype(model.clock.time)

    @test eltype(ds2["xC"]) == Float64
    @test eltype(ds2["xF"]) == Float64
    @test eltype(ds2["yC"]) == Float64
    @test eltype(ds2["yF"]) == Float64
    @test eltype(ds2["zC"]) == Float64
    @test eltype(ds2["zF"]) == Float64

    @test length(ds2["xC"]) == length(i_slice)
    @test length(ds2["xF"]) == length(i_slice)
    @test length(ds2["yC"]) == length(j_slice)
    @test length(ds2["yF"]) == length(j_slice)
    @test length(ds2["zC"]) == length(k_slice)
    @test length(ds2["zF"]) == length(k_slice)

    @test ds2["xC"][1] == grid.xC[i_slice[1]]
    @test ds2["xF"][1] == grid.xF[i_slice[1]]
    @test ds2["yC"][1] == grid.yC[j_slice[1]]
    @test ds2["yF"][1] == grid.yF[j_slice[1]]
    @test ds2["zC"][1] == grid.zC[k_slice[1]]
    @test ds2["zF"][1] == grid.zF[k_slice[1]]

    @test ds2["xC"][end] == grid.xC[i_slice[end]]
    @test ds2["xF"][end] == grid.xF[i_slice[end]]
    @test ds2["yC"][end] == grid.yC[j_slice[end]]
    @test ds2["yF"][end] == grid.yF[j_slice[end]]
    @test ds2["zC"][end] == grid.zC[k_slice[end]]
    @test ds2["zF"][end] == grid.zF[k_slice[end]]

    @test eltype(ds2["u"]) == Float32
    @test eltype(ds2["v"]) == Float32
    @test eltype(ds2["w"]) == Float32
    @test eltype(ds2["T"]) == Float32
    @test eltype(ds2["S"]) == Float32

    u_sliced = ds2["u"][:, :, :, end]
    v_sliced = ds2["v"][:, :, :, end]
    w_sliced = ds2["w"][:, :, :, end]
    T_sliced = ds2["T"][:, :, :, end]
    S_sliced = ds2["S"][:, :, :, end]

    close(ds2)

    @test all(u_sliced .≈ Array(interior(model.velocities.u))[i_slice, j_slice, k_slice])
    @test all(v_sliced .≈ Array(interior(model.velocities.v))[i_slice, j_slice, k_slice])
    @test all(w_sliced .≈ Array(interior(model.velocities.w))[i_slice, j_slice, k_slice])
    @test all(T_sliced .≈ Array(interior(model.tracers.T))[i_slice, j_slice, k_slice])
    @test all(S_sliced .≈ Array(interior(model.tracers.S))[i_slice, j_slice, k_slice])

    rm(nc_filepath)
    rm(nc_sliced_filepath)

    return nothing
end

function test_thermal_bubble_netcdf_output_with_halos(arch)
    Nx, Ny, Nz = 16, 16, 16
    Lx, Ly, Lz = 100, 100, 100

    topo = (Periodic, Periodic, Bounded)
    grid = RegularRectilinearGrid(topology=topo, size=(Nx, Ny, Nz), extent=(Lx, Ly, Lz))
    closure = IsotropicDiffusivity(ν=4e-2, κ=4e-2)
    model = IncompressibleModel(architecture=arch, grid=grid, closure=closure)
    simulation = Simulation(model, Δt=6, stop_iteration=10)

    # Add a cube-shaped warm temperature anomaly that takes up the middle 50%
    # of the domain volume.
    i1, i2 = round(Int, Nx/4), round(Int, 3Nx/4)
    j1, j2 = round(Int, Ny/4), round(Int, 3Ny/4)
    k1, k2 = round(Int, Nz/4), round(Int, 3Nz/4)
    CUDA.@allowscalar model.tracers.T.data[i1:i2, j1:j2, k1:k2] .+= 0.01

    nc_filepath = "test_dump_with_halos_$(typeof(arch)).nc"
    nc_writer = NetCDFOutputWriter(model, merge(model.velocities, model.tracers),
                                   filepath=nc_filepath,
                                   schedule=IterationInterval(10),
                                   field_slicer=FieldSlicer(with_halos=true))

    push!(simulation.output_writers, nc_writer)

    run!(simulation)

    ds = Dataset(nc_filepath)

    @test haskey(ds.attrib, "date") && !isnothing(ds.attrib["date"])
    @test haskey(ds.attrib, "Julia") && !isnothing(ds.attrib["Julia"])
    @test haskey(ds.attrib, "Oceananigans") && !isnothing(ds.attrib["Oceananigans"])
    @test haskey(ds.attrib, "schedule") && ds.attrib["schedule"] == "IterationInterval"
    @test haskey(ds.attrib, "interval") && ds.attrib["interval"] == 10
    @test haskey(ds.attrib, "output iteration interval") && !isnothing(ds.attrib["output iteration interval"])

    @test eltype(ds["time"]) == eltype(model.clock.time)

    @test eltype(ds["xC"]) == Float64
    @test eltype(ds["xF"]) == Float64
    @test eltype(ds["yC"]) == Float64
    @test eltype(ds["yF"]) == Float64
    @test eltype(ds["zC"]) == Float64
    @test eltype(ds["zF"]) == Float64

    Hx, Hy, Hz = grid.Hx, grid.Hy, grid.Hz
    @test length(ds["xC"]) == Nx+2Hx
    @test length(ds["yC"]) == Ny+2Hy
    @test length(ds["zC"]) == Nz+2Hz
    @test length(ds["xF"]) == Nx+2Hx
    @test length(ds["yF"]) == Ny+2Hy
    @test length(ds["zF"]) == Nz+2Hz+1  # z is Bounded

    @test ds["xC"][1] == grid.xC[1-Hx]
    @test ds["xF"][1] == grid.xF[1-Hx]
    @test ds["yC"][1] == grid.yC[1-Hy]
    @test ds["yF"][1] == grid.yF[1-Hy]
    @test ds["zC"][1] == grid.zC[1-Hz]
    @test ds["zF"][1] == grid.zF[1-Hz]

    @test ds["xC"][end] == grid.xC[Nx+Hx]
    @test ds["xF"][end] == grid.xF[Nx+Hx]
    @test ds["yC"][end] == grid.yC[Ny+Hy]
    @test ds["yF"][end] == grid.yF[Ny+Hy]
    @test ds["zC"][end] == grid.zC[Nz+Hz]
    @test ds["zF"][end] == grid.zF[Nz+Hz+1]  # z is Bounded

    @test eltype(ds["u"]) == Float32
    @test eltype(ds["v"]) == Float32
    @test eltype(ds["w"]) == Float32
    @test eltype(ds["T"]) == Float32
    @test eltype(ds["S"]) == Float32

    u = ds["u"][:, :, :, end]
    v = ds["v"][:, :, :, end]
    w = ds["w"][:, :, :, end]
    T = ds["T"][:, :, :, end]
    S = ds["S"][:, :, :, end]

    close(ds)

    @test all(u .≈ Array(model.velocities.u.data.parent))
    @test all(v .≈ Array(model.velocities.v.data.parent))
    @test all(w .≈ Array(model.velocities.w.data.parent))
    @test all(T .≈ Array(model.tracers.T.data.parent))
    @test all(S .≈ Array(model.tracers.S.data.parent))

    rm(nc_filepath)

    return nothing
end

function test_netcdf_function_output(arch)
    N = 16
    L = 1
    Δt = 1.25
    iters = 3

    grid = RegularRectilinearGrid(size=(N, N, N), extent=(L, 2L, 3L))
    model = IncompressibleModel(architecture=arch, grid=grid)
    simulation = Simulation(model, Δt=Δt, stop_iteration=iters)
    grid = model.grid

    # Define scalar, vector, and 2D slice outputs
    f(model) = model.clock.time^2

    g(model) = model.clock.time .* exp.(znodes(Center, grid))

    h(model) = model.clock.time .* (   sin.(xnodes(Center, grid, reshape=true)[:, :, 1])
                                    .* cos.(ynodes(Face, grid, reshape=true)[:, :, 1]))

    outputs = (scalar=f, profile=g, slice=h)
    dims = (scalar=(), profile=("zC",), slice=("xC", "yC"))

    output_attributes = (
        scalar = (longname="Some scalar", units="bananas"),
        profile = (longname="Some vertical profile", units="watermelons"),
        slice = (longname="Some slice", units="mushrooms")
    )

    global_attributes = (location="Bay of Fundy", onions=7)

    nc_filepath = "test_function_outputs_$(typeof(arch)).nc"

    simulation.output_writers[:food] =
        NetCDFOutputWriter(model, outputs; filepath=nc_filepath,
                           schedule=TimeInterval(Δt), dimensions=dims, array_type=Array{Float64}, verbose=true,
                           global_attributes=global_attributes, output_attributes=output_attributes)

    run!(simulation)

    ds = Dataset(nc_filepath, "r")

    @test haskey(ds.attrib, "date") && !isnothing(ds.attrib["date"])
    @test haskey(ds.attrib, "Julia") && !isnothing(ds.attrib["Julia"])
    @test haskey(ds.attrib, "Oceananigans") && !isnothing(ds.attrib["Oceananigans"])
    @test haskey(ds.attrib, "schedule") && !isnothing(ds.attrib["schedule"])
    @test haskey(ds.attrib, "interval") && !isnothing(ds.attrib["interval"])
    @test haskey(ds.attrib, "output time interval") && !isnothing(ds.attrib["output time interval"])

    @test eltype(ds["time"]) == eltype(model.clock.time)

    @test eltype(ds["xC"]) == Float64
    @test eltype(ds["xF"]) == Float64
    @test eltype(ds["yC"]) == Float64
    @test eltype(ds["yF"]) == Float64
    @test eltype(ds["zC"]) == Float64
    @test eltype(ds["zF"]) == Float64

    @test length(ds["xC"]) == N
    @test length(ds["yC"]) == N
    @test length(ds["zC"]) == N
    @test length(ds["xF"]) == N
    @test length(ds["yF"]) == N
    @test length(ds["zF"]) == N+1  # z is Bounded

    @test ds["xC"][1] == grid.xC[1]
    @test ds["xF"][1] == grid.xF[1]
    @test ds["yC"][1] == grid.yC[1]
    @test ds["yF"][1] == grid.yF[1]
    @test ds["zC"][1] == grid.zC[1]
    @test ds["zF"][1] == grid.zF[1]

    @test ds["xC"][end] == grid.xC[N]
    @test ds["yC"][end] == grid.yC[N]
    @test ds["zC"][end] == grid.zC[N]
    @test ds["xF"][end] == grid.xF[N]
    @test ds["yF"][end] == grid.yF[N]
    @test ds["zF"][end] == grid.zF[N+1]  # z is Bounded

    @test ds.attrib["location"] == "Bay of Fundy"
    @test ds.attrib["onions"] == 7

    @test eltype(ds["scalar"]) == Float64
    @test eltype(ds["profile"]) == Float64
    @test eltype(ds["slice"]) == Float64

    @test length(ds["time"]) == iters+1
    @test ds["time"][:] == [n*Δt for n in 0:iters]

    @test length(ds["scalar"]) == iters+1
    @test ds["scalar"].attrib["longname"] == "Some scalar"
    @test ds["scalar"].attrib["units"] == "bananas"
    @test ds["scalar"][:] == [(n*Δt)^2 for n in 0:iters]
    @test dimnames(ds["scalar"]) == ("time",)

    @test ds["profile"].attrib["longname"] == "Some vertical profile"
    @test ds["profile"].attrib["units"] == "watermelons"
    @test size(ds["profile"]) == (N, iters+1)
    @test dimnames(ds["profile"]) == ("zC", "time")

    for n in 0:iters
        @test ds["profile"][:, n+1] == n*Δt .* exp.(znodes(Center, grid))
    end

    @test ds["slice"].attrib["longname"] == "Some slice"
    @test ds["slice"].attrib["units"] == "mushrooms"
    @test size(ds["slice"]) == (N, N, iters+1)
    @test dimnames(ds["slice"]) == ("xC", "yC", "time")

    for n in 0:iters
        @test ds["slice"][:, :, n+1] == n*Δt .* (   sin.(xnodes(Center, grid, reshape=true)[:, :, 1])
                                                 .* cos.(ynodes(Face, grid, reshape=true)[:, :, 1]))
    end

    close(ds)

    #####
    ##### Take 1 more time step and test that appending to a NetCDF file works
    #####

    iters += 1
    simulation = Simulation(model, Δt=Δt, stop_iteration=iters)

    simulation.output_writers[:food] =
        NetCDFOutputWriter(model, outputs; filepath=nc_filepath, mode="a",
                           schedule=IterationInterval(1), array_type=Array{Float64}, dimensions=dims, verbose=true,
                           global_attributes=global_attributes, output_attributes=output_attributes)

    run!(simulation)

    ds = Dataset(nc_filepath, "r")

    @test length(ds["time"]) == iters+1
    @test length(ds["scalar"]) == iters+1
    @test size(ds["profile"]) == (N, iters+1)
    @test size(ds["slice"]) == (N, N, iters+1)

    @test ds["time"][:] == [n*Δt for n in 0:iters]
    @test ds["scalar"][:] == [(n*Δt)^2 for n in 0:iters]

    for n in 0:iters
        @test ds["profile"][:, n+1] == n*Δt .* exp.(znodes(Center, grid))
        @test ds["slice"][:, :, n+1] == n*Δt .* (   sin.(xnodes(Center, grid, reshape=true)[:, :, 1])
                                                 .* cos.(ynodes(Face, grid, reshape=true)[:, :, 1]))
    end

    close(ds)

    rm(nc_filepath)

    return nothing
end

function test_netcdf_time_averaging(arch)
    topo = (Periodic, Periodic, Periodic)
    domain = (x=(0, 1), y=(0, 1), z=(0, 1))
    grid = RegularRectilinearGrid(topology=topo, size=(4, 4, 4); domain...)

    λ(x, y, z) = x + (1 - y)^2 + tanh(z)
    Fc(x, y, z, t, c) = - λ(x, y, z) * c

    c_forcing = Forcing(Fc, field_dependencies=(:c,))

    model = IncompressibleModel(
                grid = grid,
        architecture = arch,
         timestepper = :RungeKutta3,
             tracers = :c,
             forcing = (c=c_forcing,),
            coriolis = nothing,
            buoyancy = nothing,
             closure = nothing
    )

    set!(model, c=1)

    Δt = 1/512  # Nice floating-point number
    simulation = Simulation(model, Δt=Δt, stop_time=50Δt)

    ∫c_dxdy = AveragedField(model.tracers.c, dims=(1, 2))

    nc_outputs = Dict("c" => ∫c_dxdy)
    nc_dimensions = Dict("c" => ("zC",))

    horizontal_average_nc_filepath = "decay_averaged_field_test.nc"
    simulation.output_writers[:horizontal_average] =
        NetCDFOutputWriter(model, nc_outputs, filepath=horizontal_average_nc_filepath, schedule=TimeInterval(10Δt),
                           dimensions=nc_dimensions, array_type=Array{Float64}, verbose=true)

    time_average_nc_filepath = "decay_windowed_time_average_test.nc"
    window = 6Δt
    stride = 2
    simulation.output_writers[:time_average] =
        NetCDFOutputWriter(model, nc_outputs, filepath=time_average_nc_filepath, array_type=Array{Float64},
                           schedule=AveragedTimeInterval(10Δt, window=window, stride=stride),
                           dimensions=nc_dimensions, verbose=true)

    run!(simulation)

    ##### Horizontal average should evaluate to
    #####
    #####     c̄(z, t) = ∫₀¹ ∫₀¹ exp{- λ(x, y, z) * t} dx dy
    #####             = 1 / (Nx*Ny) * Σᵢ₌₁ᴺˣ Σⱼ₌₁ᴺʸ exp{- λ(i, j, k) * t}
    #####
    ##### which we can compute analytically.

    ds = NCDataset(horizontal_average_nc_filepath)

    Nx, Ny, Nz = size(grid)
    xs, ys, zs = nodes(model.tracers.c)

    c̄(z, t) = 1 / (Nx * Ny) * sum(exp(-λ(x, y, z) * t) for x in xs for y in ys)

    for (n, t) in enumerate(ds["time"])
        @test ds["c"][:, n] ≈ c̄.(zs, t)
    end

    close(ds)

    #####
    ##### Test strided windowed time average against analytic solution
    #####

    ds = NCDataset(time_average_nc_filepath)

    attribute_names = ("schedule", "interval", "output time interval",
                       "time_averaging_window", "time averaging window",
                       "time_averaging_stride", "time averaging stride")

    for name in attribute_names
        @test haskey(ds.attrib, name) && !isnothing(ds.attrib[name])
    end

    c̄(ts) = 1/length(ts) * sum(c̄.(zs, t) for t in ts)

    window_size = Int(window/Δt)
    for (n, t) in enumerate(ds["time"][2:end])
        averaging_times = [t - n*Δt for n in 0:stride:window_size-1]
        @test ds["c"][:, n+1] ≈ c̄(averaging_times)
    end

    close(ds)

    rm(horizontal_average_nc_filepath)
    rm(time_average_nc_filepath)

    return nothing
end

function test_netcdf_output_alignment(arch)
    grid = RegularRectilinearGrid(size=(1, 1, 1), extent=(1, 1, 1))
    model = IncompressibleModel(architecture=arch, grid=grid)
    simulation = Simulation(model, Δt=0.2, stop_time=40)

    test_filename1 = "test_output_alignment1.nc"
    simulation.output_writers[:stuff] =
        NetCDFOutputWriter(model, model.velocities, filepath=test_filename1,
                           schedule=TimeInterval(7.3))

    test_filename2 = "test_output_alignment2.nc"
    simulation.output_writers[:something] =
        NetCDFOutputWriter(model, model.tracers, filepath=test_filename2,
                           schedule=TimeInterval(3.0))

    run!(simulation)

    Dataset(test_filename1, "r") do ds
        @test all(ds["time"] .== 0:7.3:40)
    end

    Dataset(test_filename2, "r") do ds
        @test all(ds["time"] .== 0:3.0:40)
    end

    rm(test_filename1)
    rm(test_filename2)

    return nothing
end

function test_netcdf_vertically_stretched_grid_output(arch)
    Nx = Ny = 8
    Nz = 16
    zF = [k^2 for k in 0:Nz]
    grid = VerticallyStretchedRectilinearGrid(architecture=arch, size=(Nx, Ny, Nz), x=(0, 1), y=(-π, π), z_faces=zF)

    model = IncompressibleModel(architecture=arch, grid=grid)

    Δt = 1.25
    iters = 3
    simulation = Simulation(model, Δt=Δt, stop_iteration=iters)

    nc_filepath = "test_netcdf_vertically_stretched_grid_output_$(typeof(arch)).nc"

    simulation.output_writers[:fields] =
        NetCDFOutputWriter(model, merge(model.velocities, model.tracers),
                             filepath = nc_filepath,
                             schedule = IterationInterval(1),
                           array_type = Array{Float64},
                              verbose = true)

    run!(simulation)

    grid = model.grid
    ds = NCDataset(nc_filepath)

    @test length(ds["xC"]) == Nx
    @test length(ds["yC"]) == Ny
    @test length(ds["zC"]) == Nz
    @test length(ds["xF"]) == Nx
    @test length(ds["yF"]) == Ny
    @test length(ds["zF"]) == Nz+1  # z is Bounded

    @test ds["xC"][1] == grid.xᶜᵃᵃ[1]
    @test ds["xF"][1] == grid.xᶠᵃᵃ[1]
    @test ds["yC"][1] == grid.yᵃᶜᵃ[1]
    @test ds["yF"][1] == grid.yᵃᶠᵃ[1]
    @test ds["zC"][1] == grid.zᵃᵃᶜ[1]
    @test ds["zF"][1] == grid.zᵃᵃᶠ[1]

    @test ds["xC"][end] == grid.xᶜᵃᵃ[Nx]
    @test ds["xF"][end] == grid.xᶠᵃᵃ[Nx]
    @test ds["yC"][end] == grid.yᵃᶜᵃ[Ny]
    @test ds["yF"][end] == grid.yᵃᶠᵃ[Ny]
    @test ds["zC"][end] == grid.zᵃᵃᶜ[Nz]
    @test ds["zF"][end] == grid.zᵃᵃᶠ[Nz+1]  # z is Bounded

    close(ds)
    rm(nc_filepath)

    return nothing
end

using Oceananigans.Models.HydrostaticFreeSurfaceModels: VectorInvariant

function test_netcdf_regular_lat_lon_grid_output(arch)
    Nx = Ny = Nz = 16
    grid = RegularLatitudeLongitudeGrid(size=(Nx, Ny, Nz), longitude=(-180, 180), latitude=(-80, 80), z=(-100, 0))
    model = HydrostaticFreeSurfaceModel(architecture=arch, momentum_advection = VectorInvariant(), grid=grid)

    Δt = 1.25
    iters = 3
    simulation = Simulation(model, Δt=Δt, stop_iteration=iters)

    nc_filepath = "test_netcdf_regular_lat_lon_grid_output_$(typeof(arch)).nc"

    simulation.output_writers[:fields] =
        NetCDFOutputWriter(model, merge(model.velocities, model.tracers),
                             filepath = nc_filepath,
                             schedule = IterationInterval(1),
                           array_type = Array{Float64},
                              verbose = true)

    run!(simulation)

    grid = model.grid
    ds = NCDataset(nc_filepath)

    @test length(ds["xC"]) == Nx
    @test length(ds["yC"]) == Ny
    @test length(ds["zC"]) == Nz
    @test length(ds["xF"]) == Nx
    @test length(ds["yF"]) == Ny+1  # y is Bounded
    @test length(ds["zF"]) == Nz+1  # z is Bounded

    @test ds["xC"][1] == grid.λᶜᵃᵃ[1]
    @test ds["xF"][1] == grid.λᶠᵃᵃ[1]
    @test ds["yC"][1] == grid.φᵃᶜᵃ[1]
    @test ds["yF"][1] == grid.φᵃᶠᵃ[1]
    @test ds["zC"][1] == grid.zᵃᵃᶜ[1]
    @test ds["zF"][1] == grid.zᵃᵃᶠ[1]

    @test ds["xC"][end] == grid.λᶜᵃᵃ[Nx]
    @test ds["xF"][end] == grid.λᶠᵃᵃ[Nx]
    @test ds["yC"][end] == grid.φᵃᶜᵃ[Ny]
    @test ds["yF"][end] == grid.φᵃᶠᵃ[Ny+1]  # y is Bounded
    @test ds["zC"][end] == grid.zᵃᵃᶜ[Nz]
    @test ds["zF"][end] == grid.zᵃᵃᶠ[Nz+1]  # z is Bounded

    close(ds)
    rm(nc_filepath)

    return nothing
end

for arch in archs
    @testset "NetCDF output writer [$(typeof(arch))]" begin
        @info "  Testing NetCDF output writer [$(typeof(arch))]..."
        test_DateTime_netcdf_output(arch)
        test_TimeDate_netcdf_output(arch)
        test_thermal_bubble_netcdf_output(arch)
        test_thermal_bubble_netcdf_output_with_halos(arch)
        test_netcdf_function_output(arch)
        test_netcdf_output_alignment(arch)
        test_netcdf_time_averaging(arch)
        test_netcdf_vertically_stretched_grid_output(arch)
        test_netcdf_regular_lat_lon_grid_output(arch)
    end
end
