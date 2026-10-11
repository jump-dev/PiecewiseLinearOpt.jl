# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

abstract type Method end

# Print methods as they are constructed, e.g. SixStencil(LogarithmicEmbedding())
# instead of SixStencil{LogarithmicEmbedding}(LogarithmicEmbedding()).
function Base.show(io::IO, method::Method)
    print(io, nameof(typeof(method)), "(")
    fields = (getfield(method, f) for f in fieldnames(typeof(method)))
    join(io, (_show_field(io, x) for x in fields), ", ")
    print(io, ")")
    return
end

_show_field(io::IO, x) = repr(x; context = io)

# A sub-solver given with attributes is shown by its constructor only, as the
# attributes make the output very long
function _show_field(io::IO, x::JuMP.MOI.OptimizerWithAttributes)
    return repr(x.optimizer_constructor; context = io)
end

"""
    DIRECTION

How the outputs `y` of [`piecewiselinear`](@ref) relate to the function value
`f(x)`:

- `PiecewiseLinearOpt.Graph`: `y == f(x)`
- `PiecewiseLinearOpt.Epigraph`: `y >= f(x)`
- `PiecewiseLinearOpt.Hypograph`: `y <= f(x)`
"""
@enum DIRECTION Graph Epigraph Hypograph

"""
    SegmentPointRep{D,F}(input_vals, output_vals)

A piece of a piecewise linear function from R^D to R^F, given by the points
`input_vals` and the function values `output_vals` at them. The piece is the
convex hull of the points, and the function is affine on it.
"""
struct SegmentPointRep{D,F}
    input_vals::Vector{NTuple{D,Float64}}
    output_vals::Vector{NTuple{F,Float64}}

    function SegmentPointRep{D,F}(
        input_vals::Vector{NTuple{D,Float64}},
        output_vals::Vector{NTuple{F,Float64}},
    ) where {D,F}
        if length(input_vals) != length(output_vals)
            error("Must specify the same number of input and output values.")
        end
        # What each method requires of the segments, e.g. that they are
        # simplices or aligned with a grid, is checked when formulating.
        return new{D,F}(input_vals, output_vals)
    end
end

struct _AffineFunction{D}
    coeffs::NTuple{D,Float64}
    offset::Float64
end

# The hyperplane representation of a simplex segment, used internally by the
# MultipleChoice method and for evaluation: the domain is given by f(x) ≥ 0 for
# each f in constraints, and the function value by funcs.
struct _HyperplaneSegment{D,F}
    constraints::Vector{_AffineFunction{D}}
    funcs::NTuple{F,_AffineFunction{D}}
end

"""
    SegmentStructure{D}

How the segments of a piecewise linear function of `D` inputs are arranged.
"""
abstract type SegmentStructure{D} end

"""
    Unstructured{D}()

Segments without any known structure, e.g. a triangulation in three or more
dimensions.
"""
struct Unstructured{D} <: SegmentStructure{D} end

"""
    Intervals()

The segments of a univariate function, which are intervals.
"""
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

"""
    PWLFunction(segments::Vector{SegmentPointRep{D,F}}, structure::SegmentStructure{D})

A piecewise linear function from R^D to R^F, given by its `segments`. The
`structure` describes how the segments are arranged, e.g. `Intervals()` for a
univariate function, a [`GridTriangulation`](@ref) for a bivariate function on
a grid, or `Unstructured{D}()`.

See also [`UnivariatePWLFunction`](@ref) and [`BivariatePWLFunction`](@ref) for
constructing functions from breakpoints and function values.
"""
struct PWLFunction{D,F,S<:SegmentStructure{D}}
    segments::Vector{SegmentPointRep{D,F}}
    structure::S
    # The hyperplane representation of the segments, computed on first
    # evaluation of the function
    hyperplanes::Base.RefValue{Union{Nothing,Vector{_HyperplaneSegment{D,F}}}}

    function PWLFunction(
        segments::Vector{SegmentPointRep{D,F}},
        structure::S,
    ) where {D,F,S<:SegmentStructure{D}}
        return new{D,F,S}(
            segments,
            structure,
            Ref{Union{Nothing,Vector{_HyperplaneSegment{D,F}}}}(nothing),
        )
    end
end

function PWLFunction{D,F}(
    segments::Vector{SegmentPointRep{D,F}},
    structure::SegmentStructure{D},
) where {D,F}
    return PWLFunction(segments, structure)
end

const UnivariatePWLFunction = PWLFunction{1,1}
const BivariatePWLFunction = PWLFunction{2,1}

function _check_strictly_increasing(xs::Vector{Float64}, name::String)
    for i in 1:length(xs)-1
        if !(xs[i] < xs[i+1])
            error("The breakpoints in $name must be strictly increasing.")
        end
    end
    return
end

"""
    UnivariatePWLFunction(x, f::Function)
    UnivariatePWLFunction(x, z::AbstractVector)

The piecewise linear function with the strictly increasing breakpoints `x`,
which interpolates `f(x[i])`, or `z[i]`, at the breakpoints.
"""
UnivariatePWLFunction

function PWLFunction{1,1}(x::Vector, z::Vector)
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

function PWLFunction{1,1}(x::AbstractVector, f::Function)
    d = collect(x)
    fd = [f(xi) for xi in d]
    return PWLFunction{1,1}(d, fd)
end

function PWLFunction{1,1}(x::AbstractVector, z::AbstractVector)
    return PWLFunction{1,1}(collect(x), collect(z))
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

function PWLFunction{2,1}(
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

function PWLFunction{2,1}(
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
function (pwl::PWLFunction{D,F})(x::Vararg{Real,D}) where {D,F}
    hyperplanes = pwl.hyperplanes[]
    if hyperplanes === nothing
        hyperplanes = [_hyperplane_rep(seg, Float64) for seg in pwl.segments]
        pwl.hyperplanes[] = hyperplanes
    end
    for h in hyperplanes
        if all(c -> _evaluate(c, x) ≥ -1e-9, h.constraints)
            y = map(g -> _evaluate(g, x), h.funcs)
            return F == 1 ? y[1] : y
        end
    end
    return error("The point $x is outside the domain of the function.")
end

function _evaluate(f::_AffineFunction{D}, x) where {D}
    return sum(f.coeffs[j] * x[j] for j in 1:D) + f.offset
end
