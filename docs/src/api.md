# API Reference

## Main function

```@docs
piecewiselinear
```

## PWL function types

```@docs
PWLFunction
UnivariatePWLFunction
BivariatePWLFunction
```

## Direction enum

```@docs
PiecewiseLinearOpt.DIRECTION
```

## Univariate methods

```@docs
Logarithmic
LogarithmicEmbedding
LogarithmicIndependentBranching
Incremental
NativeSOS2
ZigZagBinary
ZigZagInteger
```

## Bivariate methods

```@docs
K1
SixStencil
NineStencil
UnionJack
```

## Multivariate methods

```@docs
ConvexCombination
DisaggregatedLogarithmic
MultipleChoice
```

## Internal types

```@docs
PiecewiseLinearOpt.SegmentPointRep
PiecewiseLinearOpt.SegmentHyperplaneRep
PiecewiseLinearOpt.AffineFunction
```
