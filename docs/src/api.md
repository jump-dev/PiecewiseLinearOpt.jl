# API reference

## Adding functions

```@docs
piecewiselinear
PWLFormulation
Direction
Graph
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
ConvexPieces
```

## Bivariate methods

```@docs
SixStencil
NineStencil
GreedyTriangleSelection
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
PWLMethod
SOS2Method
TriangleSelectionMethod
formulate_pwl!
formulate_sos2!
formulate_triangle_selection!
axis_method
record!
variable_name
constrain_output!
add_biclique_cover!
cell_diagonals
Scalar
segments
input_values
output_values
```
