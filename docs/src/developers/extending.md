```@meta
CurrentModule = PiecewiseLinearOpt
```

# Extending the package

A formulation method is a subtype of [`PWLMethod`](@ref) that
implements how a function is added to a model. Three levels of extension are
possible, from the least to the most work:

1. A [`SOS2Method`](@ref) implements
   [`formulate_sos2!`](@ref), a formulation of the SOS2 constraint on a vector
   of convex multipliers. The package then formulates continuous univariate
   functions with it, and it can be used as the axis method of the bivariate
   grid methods.
2. A [`TriangleSelectionMethod`](@ref) implements
   [`formulate_triangle_selection!`](@ref), which selects a triangle in a
   cell of a grid, after the package has selected the cell with SOS2
   constraints along the axes.
3. Any other [`PWLMethod`](@ref) implements
   [`formulate_pwl!`](@ref), the whole formulation for the dimensions it
   supports.

In all cases, every variable and constraint added to the model must be passed
through [`record!`](@ref), so that the [`PWLFormulation`](@ref) can report and
delete them, and variables should be named with [`variable_name`](@ref). The
helpers [`constrain_output!`](@ref) and [`add_biclique_cover!`](@ref) cover
common parts of a formulation. The interface is exported, but adding a method
to an interface function requires importing it or writing it with the module
prefix.

The methods below are examples of the interface. They are tested, but they are
not meant for production use: `LinearBranching` has a quadratic number of
nonzeros, and `CellBinaries` one binary variable per cell, which makes them
slower than the methods of the package.

## An SOS2 method

The SOS2 constraint on multipliers ``\lambda_1, \ldots, \lambda_n`` requires
that at most two, and those adjacent, are nonzero. In the language of
independent branching, the multipliers ``\lambda_u`` and ``\lambda_v`` with
``v \ge u + 2`` must not be nonzero together. One way to enforce this uses a
binary variable ``z_k`` for each ``k = 1, \ldots, n - 2``, which separates
``\lambda_1, \ldots, \lambda_k`` from ``\lambda_{k+2}, \ldots, \lambda_n``:

```@example sos2
using JuMP, PiecewiseLinearOpt, HiGHS
import PiecewiseLinearOpt: SOS2Method, formulate_sos2!, add_biclique_cover!

struct LinearBranching <: SOS2Method end

function formulate_sos2!(formulation, λ::Vector, ::LinearBranching)
    n = length(λ)
    cover = [(collect(1:k), collect(k+2:n)) for k in 1:n-2]
    add_biclique_cover!(formulation, λ, cover)
    return
end
```

The new method works everywhere an SOS2 method does:

```@example sos2
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
@variable(model, y)
f = PiecewiseLinear(sin, 0:0.5:2π; method = LinearBranching())
g = PiecewiseLinear((a, b) -> (a - 1)^2 + (b - 2)^2, 0:0.5:3, 0:0.5:3; method = SixStencil(LinearBranching()))
@objective(model, Min, f(x) + g(x, y))
optimize!(model)
value(x), value(y), objective_value(model)
```

## A triangle selection method

A [`TriangleSelectionMethod`](@ref)
receives the matrix of multipliers of the grid points, of which only the four
corners of one cell can be nonzero, and must restrict them to one of the two
triangles of that cell. The simplest way is a binary variable for each cell,
which separates the two corners that are not joined by the diagonal of the
cell, as given by [`cell_diagonals`](@ref):

```@example triangle
using JuMP, PiecewiseLinearOpt, HiGHS
import PiecewiseLinearOpt:
    TriangleSelectionMethod,
    formulate_triangle_selection!,
    add_biclique_cover!,
    cell_diagonals

struct CellBinaries{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
end
CellBinaries() = CellBinaries(Logarithmic())

function formulate_triangle_selection!(
    formulation,
    λ::Matrix,
    nw_se::Matrix{Bool},
    ::CellBinaries,
)
    cover = Tuple{Vector{Tuple{Int,Int}},Vector{Tuple{Int,Int}}}[]
    for i in 1:size(nw_se, 1), j in 1:size(nw_se, 2)
        _, (a, b) = cell_diagonals(nw_se, i, j)
        push!(cover, ([a], [b]))
    end
    add_biclique_cover!(formulation, λ, cover)
    return
end
```

The cell is selected by the package, with the `axis_method` field that
[`axis_method`](@ref) returns by default. With a `3 × 3` grid, the method adds
four binary variables for the cells to those of the axis method:

```@example triangle
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x[1:2])
pwl = BivariatePWLFunction(0:2, 0:2, (a, b) -> a^2 - b; triangulation = UpperTriangulation())
(z,), formulation = add_formulation(model, PiecewiseLinear(pwl; method = CellBinaries()), x...)
count(is_binary, formulation.variables)
```

and formulates the function correctly, here checked at the centroid of a
triangle:

```@example triangle
fix(x[1], 0.5)
fix(x[2], 1 / 3)
@objective(model, Max, z)
optimize!(model)
value(z), pwl(0.5, 1 / 3)
```

## A general method

For functions that the two kinds above do not cover, e.g. of three or more
variables, implement [`formulate_pwl!`](@ref) for the dimensions `D` and `F`
the method supports. It receives the inputs and outputs as tuples of
[`Scalar`](@ref)s, i.e. variables, expressions or constants, and relates the
outputs to the function values with [`constrain_output!`](@ref), which takes
the [`Direction`](@ref) into account. The pieces of the function are read with
[`segments`](@ref), [`input_values`](@ref) and [`output_values`](@ref).
Functions that the method does not support should be rejected with an error;
the package then removes anything the method added to the model before the
error.

As an example, the epigraph of a convex univariate function is the
intersection of the epigraphs of its segments, which needs no binary
variables:

```@example general
using JuMP, PiecewiseLinearOpt, HiGHS
import PiecewiseLinearOpt:
    Direction,
    Scalar,
    formulate_pwl!,
    constrain_output!,
    record!,
    segments,
    input_values,
    output_values

struct ConvexEpigraph <: PWLMethod end

function formulate_pwl!(
    formulation,
    (x,)::Tuple{Scalar},
    (y,)::Tuple{Scalar},
    pwl::UnivariatePWLFunction,
    ::ConvexEpigraph,
    direction::Direction,
)
    if !(direction isa Epigraph)
        error("ConvexEpigraph only formulates the epigraph of a function.")
    end
    slope = -Inf
    for segment in segments(pwl)
        (a,), (b,) = input_values(segment)
        (fa,), (fb,) = output_values(segment)
        if (fb - fa) / (b - a) < slope
            error("ConvexEpigraph requires a convex function.")
        end
        slope = (fb - fa) / (b - a)
        constrain_output!(formulation, y, fa + slope * (x - a), direction)
    end
    # The input is restricted to the domain of the function
    breakpoints = [only(p) for s in segments(pwl) for p in input_values(s)]
    model = formulation.model
    record!(formulation, @constraint(model, x >= minimum(breakpoints)))
    record!(formulation, @constraint(model, x <= maximum(breakpoints)))
    return
end
```

The segments of a [`UnivariatePWLFunction`](@ref) are ordered, so checking
that the slopes increase verifies convexity. The method only formulates the
epigraph, which is enough when the function is minimized:

```@example general
model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, x)
pwl = UnivariatePWLFunction(-2:2, t -> t^2)
(y,), formulation = add_formulation(model, PiecewiseLinear(pwl; method = ConvexEpigraph(), direction = Epigraph()), x)
@objective(model, Min, y - 1.5x)
optimize!(model)
value(x), value(y), count(is_binary, all_variables(model))
```
