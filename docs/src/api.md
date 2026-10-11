# API reference

## Adding functions

```@docs
piecewiselinear
PWLFormulation
PiecewiseLinearOpt.Direction
PiecewiseLinearOpt.Graph
Epigraph
Hypograph
```

## Functions

```@docs
PWLFunction
UnivariatePWLFunction
BivariatePWLFunction
SegmentPointRep
SegmentStructure
Intervals
Unstructured
```

## Triangulations

```@docs
GridTriangulation
K1Triangulation
UnionJackTriangulation
UpperTriangulation
LowerTriangulation
BestFitTriangulation
RandomTriangulation
UnstructuredTriangulation
```

## Univariate methods

```@docs
LogarithmicEmbedding
LogarithmicIndependentBranching
IndependentBranching
ZigZagBinary
ZigZagInteger
Incremental
NativeSOS2
```

## Bivariate methods

```@docs
SixStencil
NineStencil
OptimalTriangleSelection
K1
UnionJack
OptimalIndependentBranching
```

## Multivariate methods

```@docs
ConvexCombination
DisaggregatedLogarithmic
MultipleChoice
```

## Extending

See [Extending the package](@ref).

```@docs
PiecewiseLinearOpt.Method
PiecewiseLinearOpt.SOS2Method
PiecewiseLinearOpt.TriangleSelectionMethod
PiecewiseLinearOpt.formulate_pwl!
PiecewiseLinearOpt.formulate_sos2!
PiecewiseLinearOpt.formulate_triangle_selection!
PiecewiseLinearOpt.axis_method
PiecewiseLinearOpt.record!
PiecewiseLinearOpt.variable_name
PiecewiseLinearOpt.constrain_output!
PiecewiseLinearOpt.add_biclique_cover!
```
