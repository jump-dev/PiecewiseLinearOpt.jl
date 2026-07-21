# Formulation Methods

PiecewiseLinearOpt.jl provides a variety of MIP (mixed-integer programming)
formulations for embedding piecewise linear functions into optimization models.
Each formulation makes different trade-offs in terms of the number of binary
variables, continuous variables, and constraints. The choice of formulation can
significantly affect solver performance.

For a detailed theoretical treatment, see:

> J. Huchette and J.P. Vielma, [_Nonconvex piecewise linear functions: Advanced
> formulations and simple modeling tools_](https://arxiv.org/abs/1708.00050),
> Operations Research, 2019.

## Univariate methods

These methods are designed for one-dimensional piecewise linear functions. Given
a function with ``n`` breakpoints (i.e., ``n - 1`` linear segments), the methods
differ in how they enforce the SOS2 (Special Ordered Set of Type 2) condition on
the convex combination weights.

| Method | Binary variables | Type |
|:-------|:-----------------|:-----|
| [`Logarithmic`](@ref) / [`LogarithmicEmbedding`](@ref) | ``\lceil \log_2(n-1) \rceil`` | Embedding |
| [`LogarithmicIndependentBranching`](@ref) | ``\lceil \log_2(n-1) \rceil`` | Independent branching |
| [`ZigZagBinary`](@ref) | ``\lceil \log_2(n-1) \rceil`` | Encoding |
| [`ZigZagInteger`](@ref) | ``\lceil \log_2(n-1) \rceil`` | Encoding (integer) |
| [`Incremental`](@ref) | ``n - 1`` | Incremental |
| [`NativeSOS2`](@ref) | 0 (uses solver SOS2) | Native |
| [`ConvexCombination`](@ref) | ``n - 1`` | Big-M |

### Logarithmic / LogarithmicEmbedding

```julia
Logarithmic()  # or equivalently, LogarithmicEmbedding()
```

The default method for univariate functions. Uses a logarithmic number of binary
variables based on reflected Gray codes. This is often the best choice for
general-purpose use.

### LogarithmicIndependentBranching

```julia
LogarithmicIndependentBranching()
```

Similar to `Logarithmic`, but uses an independent branching scheme. Also uses
a logarithmic number of binary variables.

### ZigZagBinary

```julia
ZigZagBinary()
```

Uses a zig-zag encoding with binary variables. The number of binary variables is
logarithmic in the number of breakpoints.

### ZigZagInteger

```julia
ZigZagInteger()
```

A variant of the zig-zag formulation that uses general integer variables instead
of purely binary variables. The number of integer variables is logarithmic.

### Incremental

```julia
Incremental()
```

Uses one binary variable per segment. The formulation is based on an incremental
representation of the piecewise linear function. Simple but uses more binary
variables than the logarithmic methods.

### NativeSOS2

```julia
NativeSOS2()
```

Delegates the SOS2 constraint directly to the solver using JuMP's built-in SOS2
support. This method does not add any binary variables itself, but relies on the
solver to handle the SOS2 constraint natively. Not all solvers support SOS2
constraints.

## Bivariate methods

These methods are designed for two-dimensional piecewise linear functions defined
on triangulated rectangular grids. They combine SOS2 constraints along each axis
with a triangle selection mechanism.

Each bivariate method accepts an optional `axis_method` argument that specifies
which univariate formulation to use for the axis-aligned SOS2 constraints. The
default is `Logarithmic()`.

| Method | Triangle selection variables | Grid requirement |
|:-------|:---------------------------|:-----------------|
| [`SixStencil`](@ref) | 6 binary | Any grid triangulation |
| [`NineStencil`](@ref) | 9 binary | Any grid triangulation |
| [`K1`](@ref) | 2 binary | K1 triangulation |
| [`UnionJack`](@ref) | 1 binary | Union Jack triangulation |

### SixStencil

```julia
SixStencil()              # uses Logarithmic() for axis SOS2
SixStencil(Incremental()) # uses Incremental() for axis SOS2
```

The default method for bivariate functions. Works with any grid triangulation
and uses 6 binary variables for triangle selection. A good general-purpose
choice.

### NineStencil

```julia
NineStencil()
NineStencil(Incremental())
```

Works with any grid triangulation. Uses 9 binary variables for triangle
selection.

### K1

```julia
K1()
K1(Incremental())
```

An efficient method that uses only 2 binary variables for triangle selection, but
requires the grid to use a K1 triangulation (set `pattern = :K1` when
constructing the bivariate function, which is the default).

### UnionJack

```julia
UnionJack()
UnionJack(Incremental())
```

Uses only 1 binary variable for triangle selection, but requires a Union Jack
triangulation (set `pattern = :UnionJack` when constructing the bivariate
function).

## Multivariate methods

These methods work for piecewise linear functions of any dimension.

### ConvexCombination

```julia
ConvexCombination()
```

A general-purpose method based on the convex combination (or "lambda") formulation.
Uses one binary variable per segment. Works for any dimensionality and can also
be used as a univariate method.

### DisaggregatedLogarithmic

```julia
DisaggregatedLogarithmic()
```

A disaggregated version of the logarithmic formulation that works for general
piecewise linear functions of any dimension. Uses a logarithmic number of binary
variables in the number of segments.

### MultipleChoice

```julia
MultipleChoice()
```

A multiple choice formulation for piecewise linear functions specified in
hyperplane representation (i.e., each segment is described by affine functions
and linear constraints defining its domain). This is the only method that
supports the [`PWLFunction`](@ref) with
[`SegmentHyperplaneRep`](@ref PiecewiseLinearOpt.SegmentHyperplaneRep) segments.

## Choosing a method

For **univariate** problems, the default [`Logarithmic`](@ref) method is
generally a good choice. It provides a compact formulation with a logarithmic
number of binary variables. If the solver supports native SOS2 constraints,
[`NativeSOS2`](@ref) may be faster.

For **bivariate** problems, the default [`SixStencil`](@ref) works with any
triangulation. If you can use a specific triangulation pattern, [`K1`](@ref)
or [`UnionJack`](@ref) provide more compact formulations.

For **general multivariate** problems, [`DisaggregatedLogarithmic`](@ref) or
[`ConvexCombination`](@ref) are the main options, with
[`DisaggregatedLogarithmic`](@ref) being more compact.
