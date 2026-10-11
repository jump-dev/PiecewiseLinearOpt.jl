# PiecewiseLinearOpt.jl

[![Build Status](https://github.com/jump-dev/PiecewiseLinearOpt.jl/workflows/CI/badge.svg)](https://github.com/jump-dev/PiecewiseLinearOpt.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/PiecewiseLinearOpt.jl/stable)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/PiecewiseLinearOpt.jl/dev)

[PiecewiseLinearOpt.jl](https://github.com/jump-dev/PiecewiseLinearOpt.jl) is a
[JuMP](https://jump.dev) extension for modeling optimization problems
containing piecewise linear functions, possibly nonconvex and of several
variables.

The package accompanies the paper
[_Nonconvex piecewise linear functions: Advanced formulations and simple modeling tools_](https://arxiv.org/abs/1708.00050)
by Joey Huchette and Juan Pablo Vielma, and implements the formulations from
it together with other formulations from the literature.

## Installation

```julia
import Pkg
Pkg.add("PiecewiseLinearOpt")
```

## Example

```julia
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
@variable(model, x)
@variable(model, y)
# A piecewise linear approximation of sin(x) with breakpoints 0, 0.5, ..., 6
f = piecewiselinear(model, x, 0:0.5:2π, sin)
# A piecewise linear approximation of a bivariate function on a grid, with
# a formulation chosen with the `method` keyword
g = piecewiselinear(model, x, y, 0:0.25:1, 0:0.25:1, (a, b) -> (a - 1)^2 + (b - 0.5)^2; method = SixStencil(ZigZagBinary()))
@objective(model, Min, f + g)
optimize!(model)
```

The [documentation](https://jump.dev/PiecewiseLinearOpt.jl/stable) describes
how to construct functions (univariate, bivariate on grid triangulations, or
of several variables), the available
[formulation methods](https://jump.dev/PiecewiseLinearOpt.jl/stable/methods/),
how to [extend the package](https://jump.dev/PiecewiseLinearOpt.jl/stable/extending/)
with new formulations, and how to
[migrate from version 0.4](https://jump.dev/PiecewiseLinearOpt.jl/stable/migration/).

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please open a [GitHub issue](https://github.com/jump-dev/PiecewiseLinearOpt.jl/issues/new).

## License

`PiecewiseLinearOpt.jl` is licensed under the [MIT license](https://github.com/jump-dev/PiecewiseLinearOpt.jl/blob/master/LICENSE.md).
