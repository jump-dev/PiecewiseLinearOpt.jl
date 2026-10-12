# Adding functions to a model

## Convenience methods

For univariate and bivariate functions, [`piecewiselinear`](@ref) takes the
breakpoints and the function, or its values at the breakpoints, and returns the
output:

```@example convenience
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 ≤ x ≤ 4)
@variable(model, 0 ≤ y ≤ 4)
cost_x = piecewiselinear(model, x, [0, 1, 2, 4], [0, 3, 4, 5])
cost_y = piecewiselinear(model, y, 0:4, t -> sqrt(t); method = ZigZagBinary())
@constraint(model, x + y == 5)
@objective(model, Min, cost_x + cost_y)
optimize!(model)
value(x), value(y), objective_value(model)
```

The inputs may be variables or affine expressions, and keyword arguments choose
the method, the direction and the triangulation of bivariate functions.

## The general method

The general method takes a [`PWLFunction`](@ref) and a vector of inputs, and
returns a vector of outputs together with a [`PWLFormulation`](@ref):

```@example general
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:2])
pwl = BivariatePWLFunction(0:0.5:2, 0:0.5:2, (a, b) -> a * b)
z, formulation = piecewiselinear(model, x, pwl; name = "product")
formulation
```

The formulation holds the variables and constraints that were added. Their
names start with the `name` given, or `pwl1`, `pwl2`, … by default:

```@example general
formulation.variables[1:3]
```

Deleting the formulation removes the function from the model again:

```@example general
delete(model, formulation)
num_variables(model)
```

## Outputs and directions

By default, new output variables are created, bounded by the smallest and
largest function values. Existing variables or expressions can be given with
`output_vars` in the general method, or `output_var` in the convenience
methods.

The `direction` keyword relaxes the equality between output and function
value, see [`PiecewiseLinearOpt.Direction`](@ref), and the new output
variables are then only bounded from the side of the function. For example,
with [`Epigraph`](@ref) the output is only bounded from below, which is
enough when minimizing:

```@example direction
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
y = piecewiselinear(
    model,
    x,
    [0, 1, 2, 3],
    [2, 0, 1, 0];
    direction = Epigraph(),
)
@objective(model, Min, y + x / 10)
optimize!(model)
value(x), value(y)
```

The relaxed direction gives the same optimum whenever a smaller output is never
worse with [`Epigraph`](@ref), or a larger output with [`Hypograph`](@ref).
This is the case when the output only appears in the objective with the right
sign, e.g. a cost that is minimized or a revenue that is maximized, or in a
constraint that a larger output (smaller for [`Hypograph`](@ref)) only makes
harder to satisfy, e.g. a cost in `cost ≤ budget`. The relaxed direction is then recommended, as
it is typically about 10% faster to solve, and up to 30% for some methods.

## Several outputs

A function with several outputs shares the formulation between them, which
is cheaper than formulating each output separately:

```@example outputs
using JuMP, PiecewiseLinearOpt, HiGHS
xs = [0.0, 1.0, 2.0]
segments = [
    SegmentPointRep{1,2}(
        [(xs[i],), (xs[i+1],)],
        [(sin(xs[i]), cos(xs[i])), (sin(xs[i+1]), cos(xs[i+1]))],
    ) for i in 1:2
]
pwl = PWLFunction(segments, Intervals())
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
(s, c), _ = piecewiselinear(model, [x], pwl)
@objective(model, Max, s + c)
optimize!(model)
value(x), value(s), value(c)
```

## Separable reformulation of products

A bivariate function that is built from products and univariate terms can
instead be written as a sum of univariate functions of linear expressions
(D'Ambrosio, Lodi and Martello, 2010), using

```math
x y = \frac{1}{4} \left((x + y)^2 - (x - y)^2\right).
```

As the inputs of [`piecewiselinear`](@ref) may be affine expressions, each
square is a univariate function of `x + y` or `x - y`. The univariate functions
are cheaper than a bivariate function on a grid, so they can use many more
breakpoints, which matters most when the optimum is not near the grid points.
Each term can also use the direction that suits its sign. For example, to
maximize `x y` subject to `x + 2y ≤ 4`:

```@example separable
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 ≤ x ≤ 4)
@variable(model, 0 ≤ y ≤ 2)
@constraint(model, x + 2y ≤ 4)
square(t) = t^2 / 4
# x + y is in [0, 6] and x - y in [-2, 4]
u = piecewiselinear(model, x + y, 0:0.25:6, square; direction = Hypograph())
v = piecewiselinear(model, x - y, -2:0.25:4, square; direction = Epigraph())
@objective(model, Max, u - v)
optimize!(model)
value(x), value(y), objective_value(model)
```

Whether this is faster than a bivariate grid formulation depends on the problem
and the solver, so it is worth comparing the two on representative instances.
