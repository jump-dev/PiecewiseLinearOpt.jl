# Migrating from version 0.4

Version 0.5 is a rewrite with a new interface. The main changes are:

- Methods are values of their own types instead of symbols, e.g.
  `method = Logarithmic()` instead of `method = :Logarithmic`.
- The general [`piecewiselinear`](@ref) method takes a vector of inputs and
  returns a vector of outputs together with a [`PWLFormulation`](@ref). The
  convenience methods for univariate and bivariate functions still return the
  output only.
- The triangulation of bivariate functions is given as a value with the
  `triangulation` keyword instead of the `pattern` keyword, and the default is
  now [`K1Triangulation`](@ref) instead of best fit.
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
| fields `pwl.x`, `pwl.z`, `pwl.T` | `pwl.segments`, and evaluation with `pwl(x...)` |

## Adding functions

| Version 0.4 | Version 0.5 |
|:------------|:------------|
| `z = piecewiselinear(m, x, d, f; method = :CC)` | `z = piecewiselinear(m, x, d, f; method = ConvexCombination())` |
| `z = piecewiselinear(m, x, y, dx, dy, f; method = :Log)` | `z = piecewiselinear(m, x, y, dx, dy, f; method = SixStencil())` |
| `z = piecewiselinear(m, x, pwl)` | `(z,), formulation = piecewiselinear(m, [x], pwl)` |
| `z = piecewiselinear(m, x, y, pwl)` | `(z,), formulation = piecewiselinear(m, [x, y], pwl)` |
