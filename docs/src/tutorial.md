# Tutorial

This tutorial demonstrates how to use PiecewiseLinearOpt.jl to model piecewise
linear functions in JuMP optimization models.

## Univariate piecewise linear functions

The simplest use case is approximating a univariate function with a piecewise
linear function. You specify a set of breakpoints and either a function or
explicit function values at those breakpoints.

### Using a function

```@example univariate
using JuMP, PiecewiseLinearOpt, HiGHS

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)

# Approximate x^2 on [0, 2] with breakpoints at 0, 0.5, 1.0, 1.5, 2.0
z = piecewiselinear(model, x, 0:0.5:2, xi -> xi^2)

@objective(model, Min, z)
optimize!(model)
println("x = $(value(x)), z = $(value(z))")
```

### Using explicit values

You can also pass the function values directly instead of a function:

```@example univariate
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)

breakpoints = [0.0, 1.0, 2.0, 3.0]
values = [0.0, 1.0, 4.0, 9.0]  # x^2 at each breakpoint
z = piecewiselinear(model, x, breakpoints, values)

@objective(model, Max, z)
@constraint(model, x <= 2.5)
optimize!(model)
println("x = $(value(x)), z = $(value(z))")
```

### Choosing a formulation method

PiecewiseLinearOpt.jl supports several MIP formulation methods. The default for
univariate functions is [`Logarithmic`](@ref) (an alias for
[`LogarithmicEmbedding`](@ref)), which uses ``\log_2(n)`` binary variables for
``n`` breakpoints. You can choose a different method using the `method` keyword:

```@example univariate
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)

# Use the Incremental formulation (n-1 binary variables)
z = piecewiselinear(model, x, 0:0.5:2, xi -> xi^2; method = Incremental())

@objective(model, Min, z)
optimize!(model)
println("x = $(value(x)), z = $(value(z))")
```

See [Formulation Methods](@ref) for a full list of available methods.

### Epigraph and hypograph constraints

By default, `piecewiselinear` models the graph of the function, i.e., it
enforces `z == f(x)`. You can relax this to an epigraph (`z >= f(x)`) or
hypograph (`z <= f(x)`) constraint using the `direction` keyword:

```@example univariate
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)

# Epigraph: z >= f(x)
z = piecewiselinear(
    model, x, 0:0.5:2, xi -> xi^2;
    direction = PiecewiseLinearOpt.Epigraph,
)

@objective(model, Min, z)  # minimize z subject to z >= f(x)
optimize!(model)
println("x = $(value(x)), z = $(value(z))")
```

The three options are:
- `PiecewiseLinearOpt.Graph` (default): `z == f(x)`
- `PiecewiseLinearOpt.Epigraph`: `z >= f(x)`
- `PiecewiseLinearOpt.Hypograph`: `z <= f(x)`

## Bivariate piecewise linear functions

PiecewiseLinearOpt.jl also supports bivariate piecewise linear functions defined
on triangulated rectangular grids.

```@example bivariate
using JuMP, PiecewiseLinearOpt, HiGHS

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, y)

# Approximate exp(x + y) on [0,1] × [0,1]
z = piecewiselinear(
    model, x, y,
    0:0.25:1, 0:0.25:1,
    (xi, yi) -> exp(xi + yi),
)

@objective(model, Min, z)
optimize!(model)
println("x = $(value(x)), y = $(value(y)), z = $(value(z))")
```

### Triangulation patterns

The `pattern` keyword controls how each rectangular cell in the grid is split
into two triangles. The available patterns are:

- `:K1` (default): A standard triangulation compatible with the [`K1`](@ref) method.
- `:UnionJack`: A triangulation compatible with the [`UnionJack`](@ref) method
  where the diagonal direction alternates.
- `:BestFit`: Chooses the diagonal that best approximates the function at the
  midpoint.
- `:Upper`: Chooses the diagonal that overestimates the function.
- `:Lower`: Chooses the diagonal that underestimates the function.
- `:Random`: Randomly chooses the diagonal direction.

```@example bivariate
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, y)

z = piecewiselinear(
    model, x, y,
    0:0.25:1, 0:0.25:1,
    (xi, yi) -> xi^2 + yi^2;
    pattern = :BestFit,
    method = SixStencil(),
)

@objective(model, Min, z)
optimize!(model)
println("x = $(value(x)), y = $(value(y)), z = $(value(z))")
```

## Using the general API

For more control, you can construct a [`PWLFunction`](@ref) object directly and
pass it to [`piecewiselinear`](@ref) with a tuple of input variables.

### Univariate

```@example general
using JuMP, PiecewiseLinearOpt, HiGHS

# Construct a UnivariatePWLFunction explicitly
pwl = UnivariatePWLFunction([0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 0.5, 1.5])

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)

# Pass input as a tuple
output = piecewiselinear(model, (x,), pwl)
# output is a 1-tuple; extract the variable
z = output[1]

@objective(model, Max, z)
optimize!(model)
println("x = $(value(x)), z = $(value(z))")
```

### Bivariate

```@example general
pwl = BivariatePWLFunction(
    0:0.5:1, 0:0.5:1,
    (xi, yi) -> xi * yi;
    pattern = :K1,
)

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, y)

output = piecewiselinear(model, (x, y), pwl)
z = output[1]

@objective(model, Max, z)
optimize!(model)
println("x = $(value(x)), y = $(value(y)), z = $(value(z))")
```

## Providing an existing output variable

If you already have a variable you want to use as the output, you can pass it
using the `output_var` (univariate/bivariate) or `output_vars` (general)
keyword:

```@example univariate
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, z)

piecewiselinear(model, x, 0:0.5:2, xi -> xi^2; output_var = z)

@objective(model, Max, z)
@constraint(model, x <= 1.5)
optimize!(model)
println("x = $(value(x)), z = $(value(z))")
```
