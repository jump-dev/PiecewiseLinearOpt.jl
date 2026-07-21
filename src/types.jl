# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

abstract type Method end
abstract type UnivariateMethod <: Method end

"""
    DIRECTION

Enum controlling how the output variable relates to the piecewise linear function value.

- `Graph`: The output equals the function value, i.e., `z == f(x)`.
- `Epigraph`: The output is an upper bound, i.e., `z >= f(x)`.
- `Hypograph`: The output is a lower bound, i.e., `z <= f(x)`.
"""
@enum DIRECTION Graph Epigraph Hypograph

# TODO: Make eltypes of input_vals and output_vals a type parameter
abstract type Segment{D,F} end

"""
    SegmentPointRep{D,F}

A segment of a piecewise linear function represented by its vertices (point representation).
Each segment is a simplex defined by `D + 1` vertices in `D`-dimensional input space.

# Fields
- `input_vals::Vector{NTuple{D,Float64}}`: The input coordinates of the segment vertices.
- `output_vals::Vector{NTuple{F,Float64}}`: The output values at each vertex.
"""
struct SegmentPointRep{D,F} <: Segment{D,F}
    input_vals::Vector{NTuple{D,Float64}}
    output_vals::Vector{NTuple{F,Float64}}

    function SegmentPointRep{D,F}(
        input_vals::Vector{NTuple{D,Float64}},
        output_vals::Vector{NTuple{F,Float64}},
    ) where {D,F}
        if length(input_vals) != length(output_vals)
            error("Must specify the same number of input and output values.")
        end
        # TODO: Run verifier to ensure this is actually a PWL function
        return new{D,F}(input_vals, output_vals)
    end
end

"""
    AffineFunction{D}

An affine function `f(x) = coeffs ⋅ x + offset` in `D` dimensions.
"""
struct AffineFunction{D}
    coeffs::NTuple{D,Float64}
    offset::Float64
end

"""
    SegmentHyperplaneRep{D,F}

A segment of a piecewise linear function represented by affine constraints (hyperplane
representation). The domain of the segment is defined by `f_i(x) >= 0` for each
constraint, and the function values are given by affine functions.

# Fields
- `constraints::Vector{AffineFunction{D}}`: Linear constraints defining the segment domain.
- `funcs::NTuple{F,AffineFunction{D}}`: Affine functions giving the output on this segment.
"""
struct SegmentHyperplaneRep{D,F} <: Segment{D,F}
    # Domain given by f_i(x) >= 0 where f_i is i-th constraint in constraints
    constraints::Vector{AffineFunction{D}}
    funcs::NTuple{F,AffineFunction{D}}
end

abstract type SegmentStructure{D} end

struct Intervals <: SegmentStructure{1} end

abstract type GridTriangulation <: SegmentStructure{2} end
struct UnstructuredTriangulation <: GridTriangulation end
struct K1Triangulation <: GridTriangulation end
struct UnionJackTriangulation <: GridTriangulation end

"""
    PWLFunction{D,F,T}

A piecewise linear function from ℝ^D to ℝ^F, composed of a collection of
segments of type `T`.

`D` is the input dimension, `F` is the output dimension, and `T` is the
segment representation type (either [`SegmentPointRep`](@ref PiecewiseLinearOpt.SegmentPointRep)
or [`SegmentHyperplaneRep`](@ref PiecewiseLinearOpt.SegmentHyperplaneRep)).

# Fields
- `segments::Vector{T}`: The segments composing the piecewise linear function.
- `structure::SegmentStructure{D}`: Metadata about the structure of the segments (e.g., grid triangulation type).

# Type aliases
- `UnivariatePWLFunction = PWLFunction{1,1,SegmentPointRep{1,1}}`: A univariate PWL function.
- `BivariatePWLFunction = PWLFunction{2,1,SegmentPointRep{2,1}}`: A bivariate PWL function.
"""
struct PWLFunction{D,F,T<:Segment{D,F}}
    segments::Vector{T}
    structure::SegmentStructure{D}
end

const PWLFunctionPointRep{D,F} = PWLFunction{D,F,SegmentPointRep{D,F}}
const PWLFunctionHyperplaneRep{D,F} = PWLFunction{D,F,SegmentHyperplaneRep{D,F}}

#const UnivariatePWLFunction{F} = PWLFunctionPointRep{1, F}
#const BivariatePWLFunction{F} = PWLFunctionPointRep{2, F}

"""
    UnivariatePWLFunction

A univariate piecewise linear function (ℝ → ℝ) in point representation.

# Constructors

    UnivariatePWLFunction(x::Vector, z::Vector)

Construct from breakpoints `x` and corresponding function values `z`.

    UnivariatePWLFunction(x, f::Function)

Construct from breakpoints `x` and a function `f` evaluated at those points.

# Examples
```jldoctest
julia> pwl = UnivariatePWLFunction([0.0, 1.0, 2.0], [0.0, 1.0, 0.0]);

julia> pwl = UnivariatePWLFunction(0:0.5:2, x -> x^2);
```
"""
const UnivariatePWLFunction = PWLFunctionPointRep{1,1}

"""
    BivariatePWLFunction

A bivariate piecewise linear function (ℝ² → ℝ) in point representation,
defined on a triangulated rectangular grid.

# Constructor

    BivariatePWLFunction(x, y, fz::Function; pattern=:K1, seed=...)

Construct from grid breakpoints `x` and `y`, and a function `fz(xi, yi)`.
The `pattern` keyword controls how grid rectangles are triangulated:
- `:K1` (default), `:UnionJack`, `:BestFit`, `:Upper`, `:Lower`, `:Random`

# Examples
```jldoctest
julia> pwl = BivariatePWLFunction(0:0.5:1, 0:0.5:1, (x, y) -> x + y);
```
"""
const BivariatePWLFunction = PWLFunctionPointRep{2,1}

function PWLFunctionPointRep{1,1}(x::Vector, z::Vector)
    if length(x) != length(z)
        error("Mismatch in the number of points and function values")
    end
    xs = [convert(Float64, xi) for xi in x]
    zs = [convert(Float64, zi) for zi in z]
    segments = [
        PiecewiseLinearOpt.SegmentPointRep{1,1}(
            [(xs[i],), (xs[i+1],)],
            [(zs[i],), (zs[i+1],)],
        ) for i in 1:length(x)-1
    ]

    return PWLFunction(segments, Intervals())
end

function PWLFunctionPointRep{1,1}(x, f::Function)
    d = collect(x)
    fd = [f(xi) for xi in d]
    return PWLFunctionPointRep{1,1}(d, fd)
end

function PWLFunctionPointRep{1,1}(x, z)
    return PWLFunctionPointRep{1,1}(collect(x), collect(z))
end

function PWLFunctionPointRep{2,1}(
    x,
    y,
    fz::Function;
    pattern = :K1,
    seed = hash((length(x), length(y))),
)
    xs = [convert(Float64, xi) for xi in x]
    ys = [convert(Float64, yi) for yi in y]

    segments = SegmentPointRep{2,1}[]
    structure = UnstructuredTriangulation()
    if pattern == :K1
        structure = K1Triangulation()
    elseif pattern == :UnionJack
        structure = UnionJackTriangulation()
    end

    mt = Random.MersenneTwister(seed)

    # run for each square on [x[i],x[i+1]] × [y[i],y[i+1]]
    for i in 1:length(xs)-1, j in 1:length(ys)-1
        xL, xU, yL, yU = xs[i], xs[i+1], ys[j], ys[j+1]
        mid1 = 0.5 * (fz(xL, yL) + fz(xU, yU))
        mid2 = 0.5 * (fz(xL, yU) + fz(xU, yL))
        mid3 = fz(0.5 * (xL + xU), 0.5 * (yL + yU))
        diagonal_nw_se = true
        if pattern == :Upper
            diagonal_nw_se = (mid1 > mid2)
        elseif pattern == :Lower
            diagonal_nw_se = (mid1 < mid2)
        elseif pattern == :BestFit
            diagonal_nw_se = (abs(mid1 - mid3) < abs(mid2 - mid3))
        elseif pattern == :K1
            diagonal_nw_se = false
        elseif pattern == :UnionJack
            diagonal_nw_se = isodd(i + j)
        elseif pattern == :Random
            diagonal_nw_se = rand(mt, Bool)
        end

        if diagonal_nw_se
            corners1 = [(xL, yL), (xL, yU), (xU, yL)] # SW, NW, SE
            corners2 = [(xU, yL), (xL, yU), (xU, yU)] # SE, NW, NE
        else
            corners1 = [(xL, yL), (xU, yU), (xU, yL)] # SW, NE, SE
            corners2 = [(xL, yL), (xL, yU), (xU, yU)] # SW, NW, NE
        end

        push!(
            segments,
            SegmentPointRep{2,1}(corners1, [(fz(c...),) for c in corners1]),
        )
        push!(
            segments,
            SegmentPointRep{2,1}(corners2, [(fz(c...),) for c in corners2]),
        )
    end

    return PWLFunction(segments, structure)
end
