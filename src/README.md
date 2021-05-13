# AMSC 664 Final Code Report

The actual updates are being done through the [vw/arbitrary_immersedboundary](https://github.com/CliMA/Oceananigans.jl/tree/vw/arbitrary_immersedboundary) branch in the actual Oceananigans repository. For presentation purposes, I have forked my work here. We are currently in the src folder since there is already an Oceananigans "README" file in the outermost directory.

## Implementation
The incompressible model is defined in [here](Models/IncompressibleModels/incompressible_model.jl) where ```immersed_boundary``` is included as an option for the user. Most of the changes have been implemented here, though some are in the [TimeSteppers](TimeSteppers) directory.

- Within this directory, the script that runs the Runge Kutta time-stepper [runge_kutta_3.jl](TimeSteppers/runge_kutta_3.jl), has been altered to include the velocity and tracer corrections within the stepping.

- There is a outer script function [correct_immersed_tendencies.jl](TimeSteppers/correct_immersed_tendencies.jl) within this same time-stepper folder to make "nothing" be the default action for the corrections

- The main fucntion script is the [correct_incompressible_immersed_tendencies.jl](Models/IncompressibleModels/correct_incompressible_immersed_tendencies.jl), which contains the main functions that implement the forcing. The Runge-Kutta script calls on this file at each stage of the stepper through the function ```correct_immersed_tendencies!```.

- There is some interpolation work within [interpolate.jl](Fields/interpolate.jl), where the trilinear interpolation is called with ```interpolate```. This file has been included in other interpolation work by other developers, so is not entirely work done by me. 

## Packages
The required packages to use for Oceananigans are included within the [Project.toml](../Project.toml), such that if you enter into a ```julia --project``` environment within your local Oceananigans directory, you have all of the dependencies and packages needed. Different from those that the master branch of Oceananigans started with include ```IJulia, LinearAlgebra, ForwardDiff, PyPlot, and Revise```

- [```IJulia```](https://github.com/JuliaLang/IJulia.jl) allows you to access Julia language within Jupyter interactive notebooks. These were used extensively in analysis, and some tests as well.
- [```LinearAlgebra```](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) gives you access to common norms, which was used in the error anylsis for this project.
- [```ForwardDiff```](https://github.com/JuliaDiff/ForwardDiff.jl) is part of the automatic differentiation packaging that julia has. Distinct from symbolic differnetiation or finite differencing, the function is decomposed into subdifferentials through chainrule. So in forward differencing, one first fixes the independent variable with respect to which differentiation is performed and computes the derivative of each sub-expression recursively. This is especially effective when gradients are involved, as are necessary to the geometric distancing aspect of this project.
- [```PyPlot```](https://github.com/JuliaPy/PyPlot.jl) allows julia to use the ```Matplotlib``` plotting library from Python. This was used because the usual visualization julia package, ```Plots``` did not have a good method of plotting streamlines.
- [```Revise```](https://timholy.github.io/Revise.jl/stable/) is not something necssary to work with this product, but it was used in its creation. This package keeps track of your files and edits that you make so that you do not have to restart a julia terminal to track new changes in a script.

## Examples

We have created a simple example, for a general immersed boundary method user to follow as they create their own simulations. This is in the validation folder, [flow_around_cylinder.jl](../validation/immersed_boundaries/flow_around_cylinder.jl). This test is still steady state flow, but it runs on a small domain for a short amount of time, just long enough for the user to get an idea of the differents steps.
There is another example that includes tracers in the implementation. This is given in, [flow_around_cylinder_withTemp.jl](../validation/immersed_boundaries/flow_around_cylinder_withTemp.jl)
