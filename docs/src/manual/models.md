# Adding functions to a model

## Function objects

A [`PiecewiseLinear`](@ref) holds a piecewise linear function together with
the method and direction to formulate it with. Univariate and bivariate
functions are built from the function, or its values, and the breakpoints;
other functions from a [`PWLFunction`](@ref). Keyword arguments choose the
method, the direction and the triangulation of bivariate functions.

The object belongs to no model, and is called like a function: with variables
or affine expressions as inputs, the call adds the function to their model and
returns the output.

```@example objects
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 0 ≤ x ≤ 4)
@variable(model, 0 ≤ y ≤ 4)
cost = PiecewiseLinear([0, 1, 2, 4], [0, 1, 1.5, 2])
g = PiecewiseLinear(0:4, 0:4; method = SixStencil(ZigZagBinary())) do a, b
    return 3 * sqrt(a * b)
end
@expression(model, cost_x, cost(x))
@constraint(model, cost_x ≤ 1.5)
@objective(model, Max, g(x, y) - cost_x - cost(y))
optimize!(model)
value(x), value(y), objective_value(model)
```

Each call adds a new formulation with its own variables, also for the same
inputs, so a function that is used several times, as `cost(x)` above, should be
added once and stored, e.g. with `@expression`. Above, `cost(x)` and `cost(y)`
have different inputs and need a formulation each.

With constant inputs, the call evaluates the function instead:

```@example objects
cost(1.5), g(1, 4)
```

## Formulations

[`add_formulation`](@ref) adds a function like a call does, but returns a
vector of outputs together with a [`PWLFormulation`](@ref):

```@example general
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:2])
pwl = BivariatePWLFunction(0:0.5:2, 0:0.5:2, (a, b) -> a * b)
z, formulation = add_formulation(model, PiecewiseLinear(pwl), x...; name = "product")
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
`output_var` in the call, or `output_vars` in [`add_formulation`](@ref).

The `direction` keyword relaxes the equality between output and function
value, see [`Direction`](@ref), and the new output
variables are then only bounded from the side of the function. For example,
with [`Epigraph`](@ref) the output is only bounded from below, which is
enough when minimizing:

```@example direction
using JuMP, PiecewiseLinearOpt, HiGHS
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
y = PiecewiseLinear([0, 1, 2, 3], [2, 0, 1, 0]; direction = Epigraph())(x)
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
s, c = PiecewiseLinear(pwl)(x)
@objective(model, Max, s + c)
optimize!(model)
value(x), value(s), value(c)
```

## Discontinuous functions

A univariate function with jumps is given to [`PiecewiseLinear`](@ref) in
pieces, each with its own function, or values, and breakpoints. For example,
``f(x) = -x`` for ``x < 0`` and ``f(x) = 1 + x^2`` for ``x \ge 0``:

```@example pieces
using JuMP, PiecewiseLinearOpt, HiGHS
f = PiecewiseLinear([t -> -t, t -> 1 + t^2], [[-2, 0], 0:0.5:2])
```

The default method for a discontinuous function is
[`MultipleChoice`](@ref), see [Formulation methods](@ref). A strict
inequality such as ``x < 0`` cannot be modeled, so at the jump the model may
take the value of either piece:

```@example pieces
function value_range(x0)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x == x0)
    @expression(model, fx, f(x))
    return map((MIN_SENSE, MAX_SENSE)) do sense
        @objective(model, sense, fx)
        optimize!(model)
        return value(fx)
    end
end
value_range(-1), value_range(0), value_range(1.5)
```

If the side matters, end the left piece slightly before the jump, e.g. at
`-1e-4`. The pieces may also leave a gap, which is then excluded from the
domain.

## Separable reformulation of products

A bivariate function that is built from products and univariate terms can
instead be written as a sum of univariate functions of linear expressions
(D'Ambrosio, Lodi and Martello, 2010), using

```math
x y = \frac{1}{4} \left((x + y)^2 - (x - y)^2\right).
```

As the inputs of a [`PiecewiseLinear`](@ref) may be affine expressions, each
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
u = PiecewiseLinear(square, 0:0.25:6; direction = Hypograph())(x + y)
# The square is convex, so its epigraph needs no binary variables
v = PiecewiseLinear(square, -2:0.25:4; method = ConvexPieces(), direction = Epigraph())(x - y)
@objective(model, Max, u - v)
optimize!(model)
value(x), value(y), objective_value(model)
```

Whether this is faster than a bivariate grid formulation depends on the problem
and the solver, so it is worth comparing the two on representative instances.
