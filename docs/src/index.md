```@meta
CurrentModule = PiecewiseLinearOpt
```

# PiecewiseLinearOpt.jl

[PiecewiseLinearOpt.jl](https://github.com/jump-dev/PiecewiseLinearOpt.jl) is a
[JuMP](https://github.com/jump-dev/JuMP.jl) extension for modeling optimization
problems containing piecewise linear functions. It provides MIP (mixed-integer
programming) formulations for embedding piecewise linear functions into
optimization models.

This package is an accompaniment to the paper:

> J. Huchette and J.P. Vielma, [_Nonconvex piecewise linear functions: Advanced
> formulations and simple modeling tools_](https://arxiv.org/abs/1708.00050),
> Operations Research, 2019.

## Features

- **Univariate piecewise linear functions**: Approximate any univariate function
  with a piecewise linear function defined on a set of breakpoints.
- **Bivariate piecewise linear functions**: Approximate bivariate functions on
  triangulated grids with multiple triangulation patterns.
- **Multiple MIP formulations**: Choose from a variety of formulations with
  different trade-offs in terms of the number of binary variables, continuous
  variables, and constraints.
- **Graph, epigraph, and hypograph modes**: Model the graph of the function
  exactly, or relax to epigraph/hypograph constraints.

## Installation

Install PiecewiseLinearOpt using the Julia package manager:

```julia
import Pkg
Pkg.add("PiecewiseLinearOpt")
```

## Getting help

If you need help, please ask a question on the
[JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please open a
[GitHub issue](https://github.com/jump-dev/PiecewiseLinearOpt.jl/issues/new).

## License

`PiecewiseLinearOpt.jl` is licensed under the
[MIT license](https://github.com/jump-dev/PiecewiseLinearOpt.jl/blob/master/LICENSE.md).

## Quick start

```julia
using JuMP, PiecewiseLinearOpt, HiGHS

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
# Approximate sin(x) on [0, 2π] with breakpoints every 0.5
z = piecewiselinear(model, x, 0:0.5:2pi, sin)
@objective(model, Max, z)
optimize!(model)
value(x)  # ≈ π/2
value(z)  # ≈ 1.0
```
