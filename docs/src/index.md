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
y = piecewiselinear(model, x, 0:0.5:2π, sin)
@objective(model, Min, y)
optimize!(model)
value(x), value(y)
```

The output `y` equals the piecewise linear function of `x`, which is
formulated with auxiliary variables and constraints. Functions of two
variables work the same way, with breakpoints along each axis:

```@example quickstart
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, y)
z = piecewiselinear(model, x, y, 0:0.25:1, 0:0.25:1, (a, b) -> (a - 0.3)^2 + (b - 0.6)^2)
@objective(model, Min, z)
optimize!(model)
value(x), value(y), value(z)
```

## Contents

- [Piecewise linear functions](@ref) describes how to construct functions.
- [Adding functions to a model](@ref) describes [`piecewiselinear`](@ref) and
  the formulations it returns.
- [Formulation methods](@ref) gives an overview of the available methods.
- [Migrating from version 0.4](@ref) lists the changes from previous versions.

## Getting help

If you need help, please ask a question on the
[JuMP community forum](https://jump.dev/forum). If you have a reproducible
example of a bug, please open a
[GitHub issue](https://github.com/jump-dev/PiecewiseLinearOpt.jl/issues/new).
