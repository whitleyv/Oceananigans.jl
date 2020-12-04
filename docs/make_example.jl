using Documenter
using Literate
using Plots  # to avoid capturing precompilation output by Literate
using Oceananigans

formatted_name = Dict(
    "one_dimensional_diffusion" => "One-dimensional diffusion"
)

# Gotta set this environment variable when using the GR run-time on Travis CI.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "examples")

example = ARGS[1]
example_filepath = joinpath(EXAMPLES_DIR, "$example.jl")

@info "Making example Literate: $example_filepath"
Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)

@info "Building docs..."
makedocs(
    sitename = "Oceananigans.jl",
    source = "examples",
    pages = [formatted_name[example] => "$example.md"]
)
