using Oceananigans.BoundaryConditions: PressureBoundaryConditions
using Oceananigans.Fields: CenterField

import Oceananigans.Fields: PressureFields

struct HydrostaticPressure end

function PressureFields(::HydrostaticPressure, arch, grid, bcs)
    pHY′_bcs = :pHY′ ∈ keys(bcs) ? bcs[:pHY′] : PressureBoundaryConditions(grid)
    pHY′ = CenterField(arch, grid, pHY′_bcs)
    return (pHY′=pHY′, pNHS=nothing)
end

extract_boundary_conditions(::HydrostaticPressure) = NamedTuple()

solve_for_pressure!(::Nothing, args...) = nothing


