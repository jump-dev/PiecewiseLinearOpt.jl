# PiecewiseLinearOpt.jl

[![Build Status](https://github.com/jump-dev/PiecewiseLinearOpt.jl/workflows/CI/badge.svg)](https://github.com/jump-dev/PiecewiseLinearOpt.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl)

[PiecewiseLinearOpt.jl](https://github.com/jump-dev/PiecewiseLinearOpt.jl) is a
JuMP extension for modeling optimization problems containing piecewise linear
functions.

This package is an accompaniment to a paper entitled
[_Nonconvex piecewise linear functions: Advanced formulations and simple modeling tools_](https://arxiv.org/abs/1708.00050),
by Joey Huchette and Juan Pablo Vielma.

## Getting help

If you need help, please ask a question on the [JuMP community forum](https://jump.dev/forum).

If you have a reproducible example of a bug, please open a [GitHub issue](https://github.com/jump-dev/PiecewiseLinearOpt.jl/issues/new).

## License

`PiecewiseLinearOpt.jl` is licensed under the [MIT license](https://github.com/jump-dev/PiecewiseLinearOpt.jl/blob/master/LICENSE.md).

## Installation

Install PiecewiseLinearOpt using `Pkg.add`:

```julia
import Pkg
Pkg.add("PiecewiseLinearOpt")
```

## Use with JuMP

Current support is limited to modeling the graph of a continuous piecewise
linear function, either univariate or bivariate, with the goal of adding support
for the epigraphs of lower semicontinuous piecewise linear functions.

### Univariate

Consider a piecewise linear function `f`. The function is described a domain `d`,
which is a set of breakpoints between pieces, and the function value `fd` at
those breakpoints:

```julia
julia> f(x) = sin(x)
f (generic function with 1 method)

julia> d = 0:0.5:2pi
0.0:0.5:6.0

julia> fd = f.(d)
13-element Vector{Float64}:
  0.0
  0.479425538604203
  0.8414709848078965
  0.9974949866040544
  0.9092974268256817
  0.5984721441039564
  0.1411200080598672
 -0.35078322768961984
 -0.7568024953079282
 -0.977530117665097
 -0.9589242746631385
 -0.7055403255703919
 -0.27941549819892586
```

To represent this function in a JuMP model, do:

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x)
z = PiecewiseLinearOpt.piecewiselinear(model, x, d, fd; method = :CC)
@objective(model, Min, z) # minimize f(x)
```

### Bivariate

Consider piecewise linear approximation for the function $f(x, y) = exp(x + y)$:

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x)
@variable(model, y)
z = PiecewiseLinearOpt.piecewiselinear(
    model,
    x,
    y,
    0:0.1:1,
    0:0.1:1,
    (u, v) -> exp(u + v);
    method = :DisaggLogarithmic,
)
@objective(model, Min, z)
```

## Methods

Supported univariate formulations:

* Convex combination (`:CC`)
* Multiple choice (`:MC`)
* Native SOS2 branching (`:SOS2`)
* Incremental (`:Incremental`)
* Logarithmic (`:Logarithmic`; default)
* Disaggregated Logarithmic (`:DisaggLogarithmic`)
* Binary zig-zag (`:ZigZag`)
* General integer zig-zag (`:ZigZagInteger`)

Supported bivariate formulations for entire constraint:

* Convex combination (`:CC`)
* Multiple choice (`:MC`)
* Disaggregated Logarithmic (`:DisaggLogarithmic`)

Also, you can use any univariate formulation for bivariate functions as well.
They will be used to impose two axis-aligned SOS2 constraints, along with the
"6-stencil" formulation for the triangle selection portion of the constraint.
See the associated paper for more details. In particular, the following are also
acceptable bivariate formulation choices:

* Native SOS2 branching (`:SOS2`)
* Incremental (`:Incremental`)
* Logarithmic (`:Logarithmic`)
* Binary zig-zag (`:ZigZag`)
* General integer zig-zag (`:ZigZagInteger`)
