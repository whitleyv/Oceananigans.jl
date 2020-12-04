push!(LOAD_PATH, "..")

using Documenter
using DocumenterCitations
using Literate
using Plots  # to avoid capturing precompilation output by Literate
using Oceananigans

#####
##### Build example
#####

# Gotta set this environment variable when using the GR run-time on Travis CI.
# This happens as examples will use Plots.jl to make plots and movies.
# See: https://github.com/jheinen/GR.jl/issues/278
ENV["GKSwstype"] = "100"

bib_filepath = joinpath(dirname(@__FILE__), "oceananigans.bib")
bib = CitationBibliography(bib_filepath)

const EXAMPLES_DIR = joinpath(@__DIR__, "..", "examples")
const OUTPUT_DIR   = joinpath(@__DIR__, "src", "generated")

example = ARGS[1]
example_filepath = joinpath(EXAMPLES_DIR, example * ".jl")

@info "Making example Literate: $example_filepath"
Literate.markdown(example_filepath, OUTPUT_DIR, documenter=true)

@info "Building docs..."
makedocs(bib,
	 sitename = "Oceananigans.jl",
          authors = "Ali Ramadhan, Gregory Wagner, John Marshall, Jean-Michel Campin, Chris Hill",
            pages = [example => "generated/$example.md"]
)
