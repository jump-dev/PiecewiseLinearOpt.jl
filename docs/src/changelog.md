# Release notes

The format is based on [Keep a Changelog](https://keepachangelog.com).

## v0.5.0 (unreleased)

Version 0.5 is a rewrite with a new interface. See
[Migrating from version 0.4](@ref) for the correspondence between the old and
the new names.

### Added

- Functions with several outputs, with three or more inputs, and with
  discontinuities, built from their segments with [`PWLFunction`](@ref).
- Piecewise linear relaxations, i.e. bands between a lower and an upper bound,
  with [`UnivariatePWLRelaxation`](@ref), formulated by
  [`IndependentBranching`](@ref) and [`DisaggregatedLogarithmic`](@ref).
- The `direction` keyword with [`Epigraph`](@ref) and [`Hypograph`](@ref),
  which only bound the output from one side.
- The methods [`IndependentBranching`](@ref) (Lyu, Hicks and Huchette),
  [`ConvexPieces`](@ref), [`GreedyTriangleSelection`](@ref) and
  [`MultipleChoice`](@ref) for functions of any dimension.
- [`piecewiselinear`](@ref) returns a [`PWLFormulation`](@ref) with the
  variables and constraints it added, which can be deleted from the model
  again.
- Evaluation of functions with `pwl(x...)`, and bivariate functions from a
  matrix of values.
- An extension interface with [`SOS2Method`](@ref),
  [`TriangleSelectionMethod`](@ref) and [`formulate_pwl!`](@ref), see
  [Extending the package](@ref).
- Documentation with a manual, tutorials and an API reference.

### Changed

- Methods are values of their own types instead of symbols, e.g.
  `Logarithmic()` instead of `:Logarithmic`.
- The default triangulation of bivariate functions is
  [`K1Triangulation`](@ref), and the default method depends on the
  triangulation: [`K1`](@ref), [`UnionJack`](@ref) or
  [`GreedyTriangleSelection`](@ref).
- The triangle selections of [`SixStencil`](@ref) and [`NineStencil`](@ref)
  leave out levels that separate nothing, so they use fewer binary variables on
  small grids and on triangulations with a single diagonal direction.
- [`OptimalTriangleSelection`](@ref) and [`OptimalIndependentBranching`](@ref)
  take a time limit, cache their results, and fall back to the greedy
  selection.
- The inputs of a function may be affine expressions or constants.

### Removed

- The Celaya formulations. Use [`ZigZagInteger`](@ref) for a general integer
  formulation.
