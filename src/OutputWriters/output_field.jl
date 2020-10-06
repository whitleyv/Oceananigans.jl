using Oceananigans.Fields: offset_data

struct OutputField{X, Y, Z, G, B, N, A}
    filepath :: String
    grid :: G
    boundary_conditions :: B
    name :: N
    architecture :: A

    function OutputField{X, Y, Z}(filepath, architecture, grid, name, boundary_conditions) where {X, Y, Z}

        # Extract the grid and location from file
        jldopen(output.filepath, "r") do file
            field_group = file["output/$name"]
            if boundary_conditions != nothing && !ismissing(field_group["boundary_conditions"])
                boundary_conditions = field_group["boundary_conditions"]
            end
        end

        return new{X, Y, Z, typeof(grid), typeof(boundary_conditions),
                   typeof(name), typeof(architecture)}(filepath, grid, boundary_conditions, name, architecture)
    end
end

function OutputField(filepath, name; boundary_conditions=nothing, architecture=CPU())
    jldopen(output.filepath, "r") do file
        grid = file["grid"]
        X, Y, Z = location = file["output/$name/location"]
    end

    return OutputField{X, Y, Z}(filepath, architecture, grid, name, boundary_conditions)
end

function Base.getindex(output::OutputField, time_index)

    jldopen(output.filepath, "r") do file

        field_group = file["output/$(output.name)"]
        location = field_group["location"]
        underlying_data = field_group[string(time_index)]

        data = offset_data(underlying_data, output.grid, location)

        return Field(location, output.architecture, output.grid, output.boundary_conditions, data)
    end

end

@inline function Base.getindex(output::OutputField, time_index, i, j, k)
    jldopen(output.filepath, "r") do file
        return file["output/$(output.name)"][string(time_index)][i, j, k]
    end
end
