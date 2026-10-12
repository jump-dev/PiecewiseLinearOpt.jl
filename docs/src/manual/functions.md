# Piecewise linear functions

A [`PWLFunction`](@ref) from ``\mathbb{R}^D`` to ``\mathbb{R}^F`` consists of
segments. Each segment is a [`SegmentPointRep`](@ref), given by its vertices
and the function values at them; the function is affine on the convex hull of
the vertices. Most users construct functions from breakpoints with
[`UnivariatePWLFunction`](@ref) and [`BivariatePWLFunction`](@ref).

## Univariate functions

A univariate function is given by strictly increasing breakpoints and either
the function values at them or the function itself:

```@example univariate
using PiecewiseLinearOpt
pwl = UnivariatePWLFunction([0, 1, 2, 4], [0, 2, 1, 3])
pwl(1.5)
```

```@example univariate
pwl = UnivariatePWLFunction(0:0.5:2π, sin)
pwl(1.25), sin(1.25)
```

The function is linear between neighboring breakpoints, and is only defined
between the first and the last breakpoint. The figure shows the first function
above, and a coarser approximation of `sin` with breakpoints `0:1.0:2π`, which
ends at 6, together with `sin` itself (dashed):

![Two univariate piecewise linear functions with their breakpoints](../assets/univariate.svg)

## Bivariate functions

A bivariate function is given by breakpoints along each axis, which form a
grid, and the function or a matrix `z` with `z[i, j] = f(x[i], y[j])`:

```@example bivariate
using PiecewiseLinearOpt
f(a, b) = (a + b)^2
x, y = 0:2, 0:3
pwl = BivariatePWLFunction(x, y, f)
pwl_z = BivariatePWLFunction(x, y, [f(a, b) for a in x, b in y])
pwl(0.5, 0.5), pwl_z(0.5, 0.5)
```

Each cell of the grid is split into two triangles along one of its diagonals,
as given by the `triangulation` keyword. The figure shows the triangulations of
``f(x, y) = \sin(3x) \cos(3y)`` on a ``5 \times 5`` grid, with ``x`` to the
right and ``y`` upwards, and each triangle shaded by the mean function value at
its corners, from dark (low) to light (high):

![The grid triangulations of a function on a 5 × 5 grid](../assets/triangulations.svg)

[`K1Triangulation`](@ref) and [`UnionJackTriangulation`](@ref) have a fixed
pattern, independent of the function. [`UpperTriangulation`](@ref),
[`LowerTriangulation`](@ref) and [`BestFitTriangulation`](@ref) choose the
diagonal of each cell from the function, and [`RandomTriangulation`](@ref) at
random. The choice affects the function inside the cells:

```@example bivariate
[
    t => BivariatePWLFunction(x, y, f; triangulation = t)(0.5, 0.5) for t in (
        K1Triangulation(),
        UnionJackTriangulation(),
        UpperTriangulation(),
        LowerTriangulation(),
        BestFitTriangulation(),
    )
]
```

The true value at the center of the cell is `f(0.5, 0.5) = 1`.
[`BestFitTriangulation`](@ref) chooses the diagonal that best matches the
function at the center of each cell, while [`UpperTriangulation`](@ref) and
[`LowerTriangulation`](@ref) give the largest and smallest interpolant.
[`K1Triangulation`](@ref), the default, and [`UnionJackTriangulation`](@ref)
have a fixed structure that the [`K1`](@ref) and [`UnionJack`](@ref) methods
exploit.

## General functions

Functions of any dimension, possibly with several outputs or discontinuities,
can be built from their segments. For example, a discontinuous univariate
function:

```@example general
using PiecewiseLinearOpt
segments = [
    SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)]),
    SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)]),
]
pwl = PWLFunction(segments, Intervals())
pwl(0.5), pwl(1.5)
```

and a function of three variables on a single tetrahedron, with two outputs:

```@example general
vertices = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
values = [(0.0, 1.0), (1.0, 1.0), (2.0, 0.0), (3.0, 2.0)]
pwl = PWLFunction([SegmentPointRep{3,2}(vertices, values)], Unstructured{3}())
pwl(0.25, 0.25, 0.25)
```

Not every method supports every function, see [Formulation methods](@ref).
