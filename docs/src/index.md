# PiecewiseLinearOpt.jl

[PiecewiseLinearOpt.jl](https://github.com/jump-dev/PiecewiseLinearOpt.jl) is a
[JuMP](https://jump.dev) extension for modeling optimization problems
containing piecewise linear functions, possibly nonconvex and of several
variables.

The package accompanies the paper
[_Nonconvex piecewise linear functions: Advanced formulations and simple modeling tools_](https://arxiv.org/abs/1708.00050)
by Joey Huchette and Juan Pablo Vielma, and implements the formulations from it
together with other formulations from the literature, see
[Formulation methods](@ref) and [References](@ref).

Version 0.5 is a rewrite of the code introduced with the paper, with a new
interface. The original code can be found in the
[v0.4.2 release](https://github.com/jump-dev/PiecewiseLinearOpt.jl/tree/v0.4.2),
the last of version 0.4.

## Installation

```julia
import Pkg
Pkg.add("PiecewiseLinearOpt")
```

## Quick start

To minimize a piecewise linear approximation of ``\sin(x)`` with breakpoints
at `0:0.5:2π`:

```@example quickstart
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
f = PiecewiseLinear(sin, 0:0.5:2π)
@expression(model, y, f(x))
@objective(model, Min, y)
optimize!(model)
value(x), value(y)
```

The output `y` equals the piecewise linear function of `x`, which is
formulated with auxiliary variables and constraints when `f(x)` is called. Functions of two
variables work the same way, with breakpoints along each axis:

```@example quickstart
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, y)
g = PiecewiseLinear(0:0.25:1, 0:0.25:1) do a, b
    return (a - 0.3)^2 + (b - 0.6)^2
end
@expression(model, z, g(x, y))
@objective(model, Min, z)
optimize!(model)
value(x), value(y), value(z)
```

## Contents

The manual:

- [Piecewise linear functions](@ref) describes how to construct functions.
- [Adding functions to a model](@ref) describes [`PiecewiseLinear`](@ref),
  [`add_formulation`](@ref) and the formulations they add.
- [Formulation methods](@ref) gives an overview of the available methods and
  advice on choosing one.

The tutorials show the package at work on complete models:

- [Transportation with economies of scale](@ref) approximates concave and
  fixed-charge costs of a single variable, and compares with solving the
  nonlinear model directly.
- [Hydropower scheduling](@ref) approximates a production function of two
  variables in a model with on/off decisions.

Further pages:

- [References](@ref) lists the papers behind the formulations.
- [Extending the package](@ref) describes how to add a method.
- [Release notes](@ref) and [Migrating from version 0.4](@ref) list the
  changes from previous versions.

## Getting help

If you need help, please ask a question on the
[JuMP community forum](https://jump.dev/forum). If you have a reproducible
example of a bug, please open a
[GitHub issue](https://github.com/jump-dev/PiecewiseLinearOpt.jl/issues/new).
