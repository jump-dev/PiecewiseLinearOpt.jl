# Formulation methods

The `method` keyword of [`piecewiselinear`](@ref) chooses how a function is
formulated with auxiliary variables and constraints. All methods give a
correct model; they differ in size and in how fast a solver finds the optimum.
Which method is fastest depends on the problem and the solver. The
logarithmically sized methods often scale best to functions with many
segments, see [Huchette and Vielma (2023)](https://arxiv.org/abs/1708.00050)
for a computational comparison, while [`MultipleChoice`](@ref) and
[`Incremental`](@ref) can be faster for functions with few segments, and
[`MultipleChoice`](@ref) for discontinuous functions. It is worth trying a few
methods on representative instances.

Below, `d` is the number of segments of a univariate function.

## Univariate functions

The default is [`Logarithmic`](@ref LogarithmicEmbedding).

| Method | Integer variables | Discontinuous |
|:-------|:------------------|:--------------|
| [`LogarithmicEmbedding`](@ref), also `Logarithmic()` | `ceil(log2(d))` binary | no |
| [`LogarithmicIndependentBranching`](@ref) | `ceil(log2(d))` binary | no |
| [`IndependentBranching`](@ref) | `ceil(log2(d))` binary | yes |
| [`ZigZagBinary`](@ref) | `ceil(log2(d))` binary | no |
| [`ZigZagInteger`](@ref) | `ceil(log2(d))` general integer | no |
| [`DisaggregatedLogarithmic`](@ref) | `ceil(log2(d))` binary | yes |
| [`Incremental`](@ref) | `d - 1` binary | no |
| [`ConvexCombination`](@ref) | `d` binary | no |
| [`MultipleChoice`](@ref) | `d` binary | yes |
| [`NativeSOS2`](@ref) | none, branching by the solver | no |

[`NativeSOS2`](@ref) leaves the branching to solvers that support SOS2
constraints, e.g. Gurobi; for other solvers, such as HiGHS, MathOptInterface
reformulates the constraint with binary variables.

## Bivariate functions

Bivariate functions are defined on a grid, where each cell is split into two
triangles. The default is [`SixStencil`](@ref).

The grid methods select a cell with an SOS2 constraint along each axis,
formulated with any univariate method in the table above (the `axis_method`,
by default [`Logarithmic`](@ref LogarithmicEmbedding)), and then select a
triangle in the cell:

| Method | Triangle selection | Triangulation |
|:-------|:-------------------|:--------------|
| [`SixStencil`](@ref) | 6 binary | any |
| [`NineStencil`](@ref) | 9 binary | any |
| [`OptimalTriangleSelection`](@ref) | smallest number, typically 2 to 4, computed with a solver | any |
| [`K1`](@ref) | 2 binary | [`K1Triangulation`](@ref) |
| [`UnionJack`](@ref) | 1 binary | [`UnionJackTriangulation`](@ref) |

With `n₁ × n₂` grid points and a logarithmic axis method, the cell selection
uses `ceil(log2(n₁ - 1)) + ceil(log2(n₂ - 1))` binary variables. For solvers
with support for SOS2 constraints, [`NativeSOS2`](@ref) is also a good axis
method.

Computing the optimal triangle selection of
[`OptimalTriangleSelection`](@ref) takes about a second for regular
triangulations of grids up to 8×8, but can take much longer for irregular
triangulations, e.g. [`BestFitTriangulation`](@ref) or
[`RandomTriangulation`](@ref), while saving only a few binary variables
compared with [`SixStencil`](@ref). If it takes longer than its time limit, it
falls back to [`SixStencil`](@ref).

[`OptimalIndependentBranching`](@ref) instead selects a triangle directly with
the smallest possible number of binary variables, which is only practical to
compute for small grids, up to about 4×4 points.

## Functions of more variables

[`ConvexCombination`](@ref), [`DisaggregatedLogarithmic`](@ref) and
[`MultipleChoice`](@ref) formulate functions of any number of variables. The
segments must be simplices for [`MultipleChoice`](@ref).
[`DisaggregatedLogarithmic`](@ref), the default, uses the fewest integer
variables, `ceil(log2(d))` for `d` segments. Both support discontinuous
functions.

## Example

Different methods give the same optimum, but different model sizes:

```@example methods
using JuMP, PiecewiseLinearOpt, HiGHS
f(a, b) = sin(3a) * cos(2b)
pwl = BivariatePWLFunction(0:0.1:1, 0:0.1:1, f)
for method in (
    SixStencil(),
    NineStencil(),
    K1(),
    OptimalTriangleSelection(HiGHS.Optimizer),
    SixStencil(IndependentBranching()),
    ConvexCombination(),
)
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, x[1:2])
    (z,), formulation = piecewiselinear(model, x, pwl; method = method)
    @objective(model, Min, z)
    optimize!(model)
    binaries = count(is_binary, formulation.variables)
    println(rpad(method, 50), "binaries: ", lpad(binaries, 3), "  minimum: ", round(objective_value(model); digits = 4))
end
```
