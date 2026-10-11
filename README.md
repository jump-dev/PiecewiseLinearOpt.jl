# PiecewiseLinearOpt.jl

[![Build Status](https://github.com/jump-dev/PiecewiseLinearOpt.jl/workflows/CI/badge.svg)](https://github.com/jump-dev/PiecewiseLinearOpt.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jump-dev/PiecewiseLinearOpt.jl)
[![Docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://jump.dev/PiecewiseLinearOpt.jl/stable)
[![Docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://jump.dev/PiecewiseLinearOpt.jl/dev)

[PiecewiseLinearOpt.jl](https://github.com/jump-dev/PiecewiseLinearOpt.jl) is a
JuMP extension for modeling optimization problems containing piecewise linear
functions.

This package is an accompaniment to a paper entitled
[_Nonconvex piecewise linear functions: Advanced formulations and simple modeling tools_](https://arxiv.org/abs/1708.00050),
by Joey Huchette and Juan Pablo Vielma.

See the [documentation](https://jump.dev/PiecewiseLinearOpt.jl/stable) for
details, including an overview of the formulation methods and a guide for
migrating from version 0.4.

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
on a grid where each cell is split into two triangles along the diagonal that
best fits the function:

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
    triangulation = BestFitTriangulation(),
)
@objective(model, Min, z)
```

The function can also be given by a matrix of values at the grid points,
`z[i, j] = f(x[i], y[j])`. The available triangulations are
`K1Triangulation()` (the default), `UnionJackTriangulation()`,
`UpperTriangulation()`, `LowerTriangulation()`, `BestFitTriangulation()` (not
for matrix input) and `RandomTriangulation(seed)`. A `PWLFunction` can be
evaluated at a point, e.g. `BivariatePWLFunction(x, y, f)(0.5, 0.5)`.

### General interface

The convenience methods above return the output only. The general method takes
a `PWLFunction` and a vector of inputs, and also returns a `PWLFormulation` with
the variables and constraints that were added:

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x[1:2])
pwl = BivariatePWLFunction(0:0.1:1, 0:0.1:1, (u, v) -> exp(u + v))
y, formulation = piecewiselinear(model, x, pwl; name = "cost")
formulation.variables  # the variables added, named cost_...
delete(model, formulation)  # removes the function from the model again
```

## Methods

The following formualations are available in the package and is provided through the 
`method` argument:

Supported multivariate formulations:
* `ConvexCombination()`
* `DisaggregatedLogarithmic()`
* `MultipleChoice()`: Requires each segment to be a simplex (e.g. a triangulation in 2D); supports discontinuous functions

Supported univariate formulations:
* `Incremental()`
* `IndependentBranching(scheme)`: Logarithmically sized ideal formulations from
  [Lyu, Hicks and Huchette (2023)](https://arxiv.org/abs/2304.14542), with
  `scheme` one of `:biclique` (default), `:balanced` or `:brgc`; supports
  discontinuous functions
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

