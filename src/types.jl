# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

abstract type Method end
abstract type UnivariateMethod <: Method end

@enum DIRECTION Graph Epigraph Hypograph

# TODO: Make eltypes of input_vals and output_vals a type parameter
abstract type Segment{D,F} end

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

struct AffineFunction{D}
    coeffs::NTuple{D,Float64}
    offset::Float64
end

struct SegmentHyperplaneRep{D,F} <: Segment{D,F}
    # Domain given by f_i(x) >= 0 where f_i is i-th constraint in constraints
    constraints::Vector{AffineFunction{D}}
    funcs::NTuple{F,AffineFunction{D}}
end

abstract type SegmentStructure{D} end

# Segments without any known structure, e.g. a triangulation in higher
# dimensions
struct Unstructured{D} <: SegmentStructure{D} end

struct Intervals <: SegmentStructure{1} end

"""
    GridTriangulation

A triangulation of a rectangular grid, where each cell is split into two
triangles along one of its diagonals. The subtypes describe how the diagonals
are chosen when constructing a [`BivariatePWLFunction`](@ref).
"""
abstract type GridTriangulation <: SegmentStructure{2} end

"""
    UnstructuredTriangulation()

A grid triangulation without any known structure, e.g. one built segment by
segment.
"""
struct UnstructuredTriangulation <: GridTriangulation end

"""
    K1Triangulation()

Every cell is split along the diagonal from its upper left to its lower right
corner.
"""
struct K1Triangulation <: GridTriangulation end

"""
    UnionJackTriangulation()

The diagonals alternate between neighboring cells, forming a Union Jack pattern.
"""
struct UnionJackTriangulation <: GridTriangulation end

"""
    UpperTriangulation()

Every cell is split along the diagonal with the larger average of the function
values at its end points, which gives the largest interpolant.
"""
struct UpperTriangulation <: GridTriangulation end

"""
    LowerTriangulation()

Every cell is split along the diagonal with the smaller average of the function
values at its end points, which gives the smallest interpolant.
"""
struct LowerTriangulation <: GridTriangulation end

"""
    BestFitTriangulation()

Every cell is split along the diagonal whose midpoint value is closest to the
function value at the center of the cell. This requires the function, not only
its values at the grid points.
"""
struct BestFitTriangulation <: GridTriangulation end

"""
    RandomTriangulation(seed = nothing)

Every cell is split along a random diagonal. If `seed` is `nothing`, the seed is
derived from the size of the grid.
"""
struct RandomTriangulation <: GridTriangulation
    seed::Union{Nothing,Integer}
end
RandomTriangulation() = RandomTriangulation(nothing)

struct PWLFunction{D,F,T<:Segment{D,F}}
    segments::Vector{T}
    structure::SegmentStructure{D}
end

const PWLFunctionPointRep{D,F} = PWLFunction{D,F,SegmentPointRep{D,F}}
const PWLFunctionHyperplaneRep{D,F} = PWLFunction{D,F,SegmentHyperplaneRep{D,F}}

#const UnivariatePWLFunction{F} = PWLFunctionPointRep{1, F}
#const BivariatePWLFunction{F} = PWLFunctionPointRep{2, F}

const UnivariatePWLFunction = PWLFunctionPointRep{1,1}
const BivariatePWLFunction = PWLFunctionPointRep{2,1}

function _check_strictly_increasing(xs::Vector{Float64}, name::String)
    for i in 1:length(xs)-1
        if !(xs[i] < xs[i+1])
            error("The breakpoints in $name must be strictly increasing.")
        end
    end
    return
end

function PWLFunctionPointRep{1,1}(x::Vector, z::Vector)
    if length(x) != length(z)
        error("Mismatch in the number of points and function values")
    end
    xs = [convert(Float64, xi) for xi in x]
    zs = [convert(Float64, zi) for zi in z]
    _check_strictly_increasing(xs, "x")
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

"""
    BivariatePWLFunction(x, y, f::Function; triangulation = K1Triangulation())
    BivariatePWLFunction(x, y, z::AbstractMatrix; triangulation = K1Triangulation())

The piecewise linear function on the grid with breakpoints `x` and `y`, which
interpolates `f(x[i], y[j])`, or `z[i, j]`, at the grid points. Each cell of
the grid is split into two triangles as given by `triangulation`, see
[`GridTriangulation`](@ref).
"""
BivariatePWLFunction

function PWLFunctionPointRep{2,1}(
    x,
    y,
    f::Function;
    triangulation::GridTriangulation = K1Triangulation(),
)
    xs = [convert(Float64, xi) for xi in x]
    ys = [convert(Float64, yi) for yi in y]
    z = [convert(Float64, f(xi, yj)) for xi in xs, yj in ys]
    center(i, j) = f(0.5 * (xs[i] + xs[i+1]), 0.5 * (ys[j] + ys[j+1]))
    return _grid_pwl_function(xs, ys, z, triangulation, center)
end

function PWLFunctionPointRep{2,1}(
    x,
    y,
    z::AbstractMatrix;
    triangulation::GridTriangulation = K1Triangulation(),
)
    if size(z) != (length(x), length(y))
        error(
            "The matrix of function values must have size $((length(x), length(y))), got $(size(z)).",
        )
    end
    xs = [convert(Float64, xi) for xi in x]
    ys = [convert(Float64, yi) for yi in y]
    zs = [convert(Float64, zij) for zij in z]
    return _grid_pwl_function(xs, ys, zs, triangulation, nothing)
end

function _grid_pwl_function(
    xs::Vector{Float64},
    ys::Vector{Float64},
    z::Matrix{Float64},
    triangulation::GridTriangulation,
    center::Union{Nothing,Function},
)
    _check_strictly_increasing(xs, "x")
    _check_strictly_increasing(ys, "y")
    nw_se = _nw_se_diagonals(triangulation, xs, ys, z, center)
    segments = SegmentPointRep{2,1}[]
    for i in 1:length(xs)-1, j in 1:length(ys)-1
        sw, se, nw, ne = (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)
        if nw_se[i, j]
            triangles = ((sw, nw, se), (se, nw, ne))
        else
            triangles = ((sw, ne, se), (sw, nw, ne))
        end
        for t in triangles
            push!(
                segments,
                SegmentPointRep{2,1}(
                    [(xs[k], ys[l]) for (k, l) in t],
                    [(z[k, l],) for (k, l) in t],
                ),
            )
        end
    end
    return PWLFunction(segments, triangulation)
end

# Returns whether each cell (i, j) is split along its NW-SE diagonal, in which
# case the value at the center of the cell is the average of the NW and SE
# values, instead of the SW and NE values.
function _nw_se_diagonals(::K1Triangulation, xs, ys, z, center)
    return fill(true, length(xs) - 1, length(ys) - 1)
end

function _nw_se_diagonals(::UnionJackTriangulation, xs, ys, z, center)
    return [isodd(i + j) for i in 1:length(xs)-1, j in 1:length(ys)-1]
end

function _nw_se_diagonals(::UpperTriangulation, xs, ys, z, center)
    return [
        z[i, j+1] + z[i+1, j] > z[i, j] + z[i+1, j+1] for i in 1:length(xs)-1,
        j in 1:length(ys)-1
    ]
end

function _nw_se_diagonals(::LowerTriangulation, xs, ys, z, center)
    return [
        z[i, j+1] + z[i+1, j] < z[i, j] + z[i+1, j+1] for i in 1:length(xs)-1,
        j in 1:length(ys)-1
    ]
end

function _nw_se_diagonals(::BestFitTriangulation, xs, ys, z, center)
    if center === nothing
        error(
            "BestFitTriangulation requires the function itself, not only its values at the grid points.",
        )
    end
    return [
        abs(0.5 * (z[i, j+1] + z[i+1, j]) - center(i, j)) <
        abs(0.5 * (z[i, j] + z[i+1, j+1]) - center(i, j)) for
        i in 1:length(xs)-1, j in 1:length(ys)-1
    ]
end

function _nw_se_diagonals(t::RandomTriangulation, xs, ys, z, center)
    seed = something(t.seed, hash((length(xs), length(ys))))
    return rand(
        Random.MersenneTwister(seed),
        Bool,
        length(xs) - 1,
        length(ys) - 1,
    )
end

function _nw_se_diagonals(t::GridTriangulation, xs, ys, z, center)
    return error("Cannot construct a grid triangulation of type $(typeof(t)).")
end

"""
    (pwl::PWLFunction)(x::Real...)

Evaluates the piecewise linear function at `x`, returning a number for a
function with one output and a tuple otherwise. Every segment must be a simplex.
Where the function is discontinuous, the value of the first segment containing
`x` is returned.
"""
function (pwl::PWLFunctionPointRep{D,F})(x::Vararg{Real,D}) where {D,F}
    for seg in pwl.segments
        h = _hyperplane_rep(seg)
        if all(c -> _evaluate(c, x) ≥ -1e-9, h.constraints)
            y = map(g -> _evaluate(g, x), h.funcs)
            return F == 1 ? y[1] : y
        end
    end
    return error("The point $x is outside the domain of the function.")
end

function _evaluate(f::AffineFunction{D}, x) where {D}
    return sum(f.coeffs[j] * x[j] for j in 1:D) + f.offset
end
