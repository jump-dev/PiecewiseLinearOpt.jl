# Migrating from version 0.4

Version 0.5 is a rewrite with a new interface. The main changes are:

- Methods are values of their own types instead of symbols, e.g.
  `method = Logarithmic()` instead of `method = :Logarithmic`.
- `piecewiselinear` is replaced by [`PiecewiseLinear`](@ref) objects, which
  hold the function, the method and the direction, and are called with the
  inputs to add the function to their model, e.g. `f(x)`. The function is given
  before the breakpoints, so that a `do` block can be used.
  [`add_formulation`](@ref) also returns the [`PWLFormulation`](@ref) with the
  variables and constraints that were added.
- The triangulation of bivariate functions is given as a value with the
  `triangulation` keyword instead of the `pattern` keyword, and the default is
  now [`K1Triangulation`](@ref) instead of best fit. The default method of
  bivariate functions depends on the triangulation, see
  [Formulation methods](@ref), instead of always being the six-stencil
  formulation; pass `method = SixStencil()` to keep it.
- Functions may have several outputs, more than two inputs and, for some
  methods, discontinuities.
- The `direction` keyword takes `Epigraph()` or `Hypograph()` for a function
  that is only bounded from one side, instead of modeling its graph.

## Methods

| Version 0.4 | Version 0.5 |
|:------------|:------------|
| `:Logarithmic`, `:Log` | [`Logarithmic()`](@ref LogarithmicEmbedding) |
| `:LogarithmicIB`, `:LogIB` | [`LogarithmicIndependentBranching()`](@ref LogarithmicIndependentBranching) |
| `:ZigZag`, `:ZZB` | [`ZigZagBinary()`](@ref ZigZagBinary) |
| `:ZigZagInteger`, `:ZZI` | [`ZigZagInteger()`](@ref ZigZagInteger) |
| `:GeneralizedCelaya`, `:SymmetricCelaya` | removed, use e.g. [`ZigZagInteger()`](@ref ZigZagInteger) for a general integer formulation |
| `:SOS2` | [`NativeSOS2()`](@ref NativeSOS2) |
| `:CC` | [`ConvexCombination()`](@ref ConvexCombination) |
| `:MC` | [`MultipleChoice()`](@ref MultipleChoice) |
| `:DisaggLogarithmic`, `:DLog` | [`DisaggregatedLogarithmic()`](@ref DisaggregatedLogarithmic) |
| `:Incremental` | [`Incremental()`](@ref Incremental) |
| `:OptimalIB` with `subsolver = s` | [`OptimalIndependentBranching(s)`](@ref OptimalIndependentBranching) |
| bivariate with an SOS2 method `m`, e.g. `:Logarithmic` | [`SixStencil(m)`](@ref SixStencil), [`NineStencil(m)`](@ref NineStencil), [`K1(m)`](@ref K1), [`UnionJack(m)`](@ref UnionJack) or [`OptimalTriangleSelection(s, m)`](@ref OptimalTriangleSelection) |

## Functions and triangulations

| Version 0.4 | Version 0.5 |
|:------------|:------------|
| `UnivariatePWLFunction(x, z)` | [`UnivariatePWLFunction(x, z)`](@ref UnivariatePWLFunction), unchanged |
| `BivariatePWLFunction(x, y, f; pattern = :BestFit)` | [`BivariatePWLFunction(x, y, f; triangulation = BestFitTriangulation())`](@ref BivariatePWLFunction) |
| `pattern = :K1` | `triangulation = K1Triangulation()` |
| `pattern = :UnionJack` | `triangulation = UnionJackTriangulation()` |
| `pattern = :Upper` | `triangulation = UpperTriangulation()` |
| `pattern = :Lower` | `triangulation = LowerTriangulation()` |
| `pattern = :Random` | `triangulation = RandomTriangulation(seed)` |
| fields `pwl.x`, `pwl.z`, `pwl.T` | `segments(pwl)`, and evaluation with `pwl(x...)` |

## Adding functions

| Version 0.4 | Version 0.5 |
|:------------|:------------|
| `z = piecewiselinear(m, x, d, f; method = :CC)` with a function `f` | `z = PiecewiseLinear(f, d; method = ConvexCombination())(x)` |
| `z = piecewiselinear(m, x, d, fd)` with values `fd` | `z = PiecewiseLinear(d, fd)(x)` |
| `z = piecewiselinear(m, x, y, dx, dy, f; method = :Log)` | `z = PiecewiseLinear(f, dx, dy; method = SixStencil())(x, y)` |
| `z = piecewiselinear(m, x, pwl)` | `z = PiecewiseLinear(pwl)(x)`, or `(z,), formulation = add_formulation(m, PiecewiseLinear(pwl), x)` |
| `z = piecewiselinear(m, x, y, pwl)` | `z = PiecewiseLinear(pwl)(x, y)` |
