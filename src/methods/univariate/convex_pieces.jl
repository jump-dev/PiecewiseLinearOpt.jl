# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ConvexPieces(encoding = :mc)

Formulation of univariate piecewise linear functions with a binary variable for
each convex piece.

It is meant for the [`Epigraph`](@ref) or [`Hypograph`](@ref) direction, where
only the convex pieces of the function need integer variables, rather than its
segments. With [`Epigraph`](@ref), consecutive segments that join continuously with
nondecreasing slopes form a convex piece, whose epigraph is a convex
polyhedron. The epigraph of the function is the union of the epigraphs of its
pieces, which is formulated with a copy of the input and output for each piece,
as [`MultipleChoice`](@ref) does for segments (Balas, 1998). With
[`Hypograph`](@ref), the pieces are concave instead. A function that is convex
in the direction needs no integer variables, and e.g. an S-shaped utility that
is maximized needs one for the concave part and one for each segment of the
convex part. With the default direction, every segment is a piece, which gives
the [`MultipleChoice`](@ref) formulation.

With `encoding = :mc`, there is a binary variable for each piece, and the
formulation is ideal. With `encoding = :log`, the active piece is encoded with
`ceil(log2(p))` binary variables for `p` pieces, as in
[`DisaggregatedLogarithmic`](@ref).

The function may be discontinuous, and the pieces then end at the jumps.
"""
struct ConvexPieces <: PWLMethod
    encoding::Symbol

    function ConvexPieces(encoding::Symbol = :mc)
        if !(encoding in (:mc, :log))
            error(
                "Unknown encoding $(repr(encoding)) for ConvexPieces, expected :mc or :log.",
            )
        end
        return new(encoding)
    end
end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::Tuple{Scalar},
    output_vars::Tuple{Scalar},
    pwl::PWLFunction{1,1},
    method::ConvexPieces,
    direction::Direction,
)
    model = formulation.model
    x, = input_vars
    y, = output_vars
    pieces = _convex_pieces(pwl, direction)
    if length(pieces) == 1
        # The function is convex in the direction
        piece = only(pieces)
        record!(formulation, JuMP.@constraint(model, x ≥ first(piece).lower))
        record!(formulation, JuMP.@constraint(model, x ≤ last(piece).upper))
        for f in piece
            constrain_output!(formulation, y, f.slope * x + f.offset, direction)
        end
        return nothing
    end
    n = length(pieces)
    x_hat = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:n],
            base_name = variable_name(formulation, "x_hat")
        ),
    )
    y_hat = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:n],
            base_name = variable_name(formulation, "y_hat")
        ),
    )
    z = _piece_selection!(formulation, n, method.encoding)
    record!(formulation, JuMP.@constraint(model, sum(x_hat) == x))
    for (p, piece) in enumerate(pieces)
        record!(
            formulation,
            JuMP.@constraint(model, x_hat[p] ≥ first(piece).lower * z[p]),
        )
        record!(
            formulation,
            JuMP.@constraint(model, x_hat[p] ≤ last(piece).upper * z[p]),
        )
        for f in piece
            record!(
                formulation,
                _output_constraint(
                    model,
                    y_hat[p],
                    f.slope * x_hat[p] + f.offset * z[p],
                    direction,
                ),
            )
        end
    end
    constrain_output!(formulation, y, sum(y_hat), direction)
    return nothing
end

# The variables z[p] that select piece p of n
function _piece_selection!(formulation::PWLFormulation, n::Int, encoding)
    model = formulation.model
    if encoding == :mc
        z = record!(
            formulation,
            JuMP.@variable(
                model,
                [1:n],
                Bin,
                base_name = variable_name(formulation, "z")
            ),
        )
    else
        z = record!(
            formulation,
            JuMP.@variable(
                model,
                [1:n],
                lower_bound = 0,
                upper_bound = 1,
                base_name = variable_name(formulation, "z")
            ),
        )
        r = ceil(Int, log2(n))
        codes = _reflected_gray_codes(r)
        w = record!(
            formulation,
            JuMP.@variable(
                model,
                [1:r],
                Bin,
                base_name = variable_name(formulation, "w")
            ),
        )
        for k in 1:r
            record!(
                formulation,
                JuMP.@constraint(
                    model,
                    sum(codes[p][k] * z[p] for p in 1:n) == w[k]
                ),
            )
        end
    end
    record!(formulation, JuMP.@constraint(model, sum(z) == 1))
    return z
end

# A segment of a univariate function, as the affine function on an interval
struct _IntervalSegment
    lower::Float64
    upper::Float64
    slope::Float64
    offset::Float64
end

_value(f::_IntervalSegment, x) = f.slope * x + f.offset

# The segments of `pwl` in order, grouped into the pieces that are convex in
# `direction`: consecutive segments that join continuously, with nondecreasing
# slopes for Epigraph and nonincreasing slopes for Hypograph. With Graph, every
# segment is a piece.
function _convex_pieces(pwl::PWLFunction{1,1}, direction::Direction)
    intervals = _IntervalSegment[]
    for segment in segments(pwl)
        if length(input_values(segment)) != 2
            error(
                "Every segment of a univariate piecewise linear function must have two points.",
            )
        end
        ((a,), (b,)), ((fa,), (fb,)) =
            input_values(segment), output_values(segment)
        if a == b
            error(
                "ConvexPieces requires segments of positive length, but a segment has the points $(input_values(segment)).",
            )
        elseif a > b
            a, b, fa, fb = b, a, fb, fa
        end
        slope = (fb - fa) / (b - a)
        push!(intervals, _IntervalSegment(a, b, slope, fa - slope * a))
    end
    sort!(intervals; by = f -> (f.lower, f.upper))
    pieces = Vector{_IntervalSegment}[]
    for (k, f) in enumerate(intervals)
        if k > 1
            g = intervals[k-1]
            continuous =
                g.upper == f.lower && isapprox(
                    _value(g, f.lower),
                    _value(f, f.lower);
                    atol = 1e-9 * (1 + abs(_value(f, f.lower))),
                )
            if continuous && _is_convex_join(g.slope, f.slope, direction)
                push!(pieces[end], f)
                continue
            end
        end
        push!(pieces, [f])
    end
    return pieces
end

_is_convex_join(left, right, ::Graph) = false

function _is_convex_join(left, right, ::Epigraph)
    return right ≥ left - 1e-9 * (1 + abs(left))
end

function _is_convex_join(left, right, ::Hypograph)
    return right ≤ left + 1e-9 * (1 + abs(left))
end
