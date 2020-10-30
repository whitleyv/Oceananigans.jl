#####
##### << Indexing >>
#####

@inline left_halo_indices(loc, topo, N, H) = 1-H:0
@inline left_halo_indices(::Type{Nothing}, topo, N, H) = 1:0 # empty

@inline right_halo_indices(loc, topo, N, H) = N+1:N+H
@inline right_halo_indices(::Type{Face}, ::Type{Bounded}, N, H) = N+2:N+1+H
@inline right_halo_indices(::Type{Nothing}, topo, N, H) = 1:0 # empty

@inline underlying_left_halo_indices(loc, topo, N, H) = 1:H
@inline underlying_left_halo_indices(::Type{Nothing}, topo, N, H) = 1:0 # empty

@inline underlying_right_halo_indices(loc, topo, N, H) = N+1+H:N+2H
@inline underlying_right_halo_indices(::Type{Face}, ::Type{Bounded}, N, H) = N+2+H:N+1+2H
@inline underlying_right_halo_indices(::Type{Nothing}, topo, N, H) = 1:0 # empty

@inline interior_indices(loc, topo, N) = 1:N
@inline interior_indices(::Type{Face}, ::Type{Bounded}, N) = 1:N+1
@inline interior_indices(::Type{Nothing}, topo, N) = 1:1

@inline interior_x_indices(loc, grid) = interior_indices(loc, topology(grid, 1), grid.Nx)
@inline interior_y_indices(loc, grid) = interior_indices(loc, topology(grid, 2), grid.Ny)
@inline interior_z_indices(loc, grid) = interior_indices(loc, topology(grid, 3), grid.Nz)

@inline interior_parent_indices(loc, topo, N, H) = 1+H:N+H
@inline interior_parent_indices(::Type{Face}, ::Type{Bounded}, N, H) = 1+H:N+1+H
@inline interior_parent_indices(::Type{Nothing}, topo, N, H) = 1:1

# All indices including halos.
@inline all_indices(loc, topo, N, H) = 1-H:N+H
@inline all_indices(::Type{Face}, ::Type{Bounded}, N, H) = 1-H:N+1+H
@inline all_indices(::Type{Nothing}, topo, N, H) = 1:1

@inline all_x_indices(loc, grid) = all_indices(loc, topology(grid, 1), grid.Nx, grid.Hx)
@inline all_y_indices(loc, grid) = all_indices(loc, topology(grid, 2), grid.Ny, grid.Hy)
@inline all_z_indices(loc, grid) = all_indices(loc, topology(grid, 3), grid.Nz, grid.Hz)

@inline all_parent_indices(loc, topo, N, H) = 1:N+2H
@inline all_parent_indices(::Type{Face}, ::Type{Bounded}, N, H) = 1:N+1+2H
@inline all_parent_indices(::Type{Nothing}, topo, N, H) = 1:1

@inline all_parent_x_indices(loc, grid) = all_parent_indices(loc, topology(grid, 1), grid.Nx, grid.Hx)
@inline all_parent_y_indices(loc, grid) = all_parent_indices(loc, topology(grid, 2), grid.Ny, grid.Hy)
@inline all_parent_z_indices(loc, grid) = all_parent_indices(loc, topology(grid, 3), grid.Nz, grid.Hz)
