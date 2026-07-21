# PiecewiseLinearOpt.jl

[![Build Status](https://github.com/jump-dev/PiecewiseLinearOpt.jl/workflows/CI/badge.svg)](https://github.com/jump-dev/PiecewiseLinearOpt.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl)
[![Aqua QA](https://juliatesting.github.io/Aqua.jl/dev/assets/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

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
linear function, with a primary focus on univariate or bivariate functions.
There are also methods for more general multivariate problems.

### Univariate

Consider a piecewise linear function described by a domain `d`,
which is a set of breakpoints between pieces, and the function value at
those breakpoints given by the function `f` at those points:

```julia
julia> d = 0:0.5:2pi
0.0:0.5:6.0

julia> f(x) = sin(x)
f (generic function with 1 method)
```

To represent this function in a JuMP model, do:

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x)
z = PiecewiseLinearOpt.piecewiselinear(model, x, d, f; method = Logarithmic())
@objective(model, Min, z) # minimize f(x)
```

### Bivariate

Consider a piecewise linear approximation for the function $f(x, y) = exp(x + y)$
on a triangular grid with a best fit pattern: 

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
    method = SixStencil(),
    pattern = :BestFit
)
@objective(model, Min, z)
```

## Methods

The following formualations are available in the package and is provided through the 
`method` argument:

Supported multivariate formulations:
* `ConvexCombination()`
* `DisaggregatedLogarithmic()`
* `MultipleChoice()`: Limited support as it currently needs an explicit formulations with hyperplanes

Supported univariate formulations:
* `Incremental()`
* `Logarithmic()`
* `LogarithmicIndependentBranching()`
* `NativeSOS2()`
* `ZigZagBinary()`
* `ZigZagInteger()`

The following bivariate formulations are available and can be combined with most univariate 
formulations to impose two axis-aligned SOS2 constraints. See the associated paper for more details.
* `K1(sos2_method)`: requires a K1 grid triangulation 
* `UnionJack(sos2_method)`: requires a UnionJack grid triangulation 
* `SixStencil(sos2_method)`: requires a grid triangulation
* `NineStencil(sos2_method)`: requires a grid triangulation

