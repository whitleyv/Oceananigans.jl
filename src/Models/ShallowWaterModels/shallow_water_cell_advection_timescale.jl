"Returns the time-scale for advection on a regular grid across a single grid cell 
 for ShallowWaterModel."

import Oceananigans.Utils: cell_advection_timescale
using Oceananigans.Fields: interior

function shallow_water_cell_advection_timescale(uh, vh, h, grid)
    uhmax = maximum(abs, interior(uh))
    vhmax = maximum(abs, interior(vh))
    hmin  = minimum(abs, interior(h))

    return min(grid.Δx / uhmax, grid.Δy / vhmax) * hmin
end

cell_advection_timescale(model::ShallowWaterModel) =
    shallow_water_cell_advection_timescale(
        model.solution.uh, 
        model.solution.vh,
        model.solution.h,
        model.grid
        )
