# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    MultipleChoice()

The multiple choice formulation of piecewise linear functions of any dimension
(see Vielma, Ahmed and Nemhauser, 2010), with a copy of the input and a binary
variable for each segment. The segments must be simplices, but the function
may be discontinuous.
"""
struct MultipleChoice <: Method end

# Each segment gets its own copy of the input, so the function need not be
# continuous.
function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{D,_Scalar},
    output_vars::NTuple{F,_Scalar},
    pwl::PWLFunction{D,F},
    method::MultipleChoice,
    direction::DIRECTION,
) where {D,F}
    segments = [_hyperplane_rep(seg, Rational{BigInt}) for seg in pwl.segments]
    _formulate_multiple_choice!(
        formulation,
        input_vars,
        output_vars,
        segments,
        direction,
    )
    return nothing
end

function _formulate_multiple_choice!(
    formulation::PWLFormulation,
    input_vars::NTuple{D,_Scalar},
    output_vars::NTuple{F,_Scalar},
    segments::Vector{_HyperplaneSegment{D,F}},
    direction::DIRECTION,
) where {D,F}
    model = formulation.model
    n = length(segments)

    x_hat = _record!(
        formulation,
        JuMP.@variable(
            model,
            [1:n, 1:D],
            base_name = _name(formulation, "x_hat")
        ),
    )
    z = _record!(
        formulation,
        JuMP.@variable(model, [1:n], Bin, base_name = _name(formulation, "z")),
    )
    _record!(formulation, JuMP.@constraint(model, sum(z) == 1))
    for j in 1:D
        _record!(
            formulation,
            JuMP.@constraint(model, sum(x_hat[:, j]) == input_vars[j]),
        )
    end
    for (i, seg) in enumerate(segments), constraint in seg.constraints
        _record!(
            formulation,
            JuMP.@constraint(
                model,
                _homogenize(constraint, x_hat[i, :], z[i]) ≥ 0
            ),
        )
    end
    for k in 1:F
        rhs = sum(
            _homogenize(seg.funcs[k], x_hat[i, :], z[i]) for
            (i, seg) in enumerate(segments)
        )
        _constrain_output_var(formulation, output_vars[k], rhs, direction)
    end
    return nothing
end

# Returns f.coeffs ⋅ x + f.offset * z, which equals f(x) when z = 1 and is
# zero when z = 0 (forcing x = 0, as segments are bounded).
function _homogenize(f::_AffineFunction{D}, x, z) where {D}
    return sum(f.coeffs[j] * x[j] for j in 1:D) + f.offset * z
end

# Converts a simplex with vertices v_0, ..., v_D into the hyperplane
# representation given by its barycentric coordinates μ_0, ..., μ_D ≥ 0,
# computing in the number type T.
#
# With T = Rational{BigInt} the computation is exact, as every Float64 is a
# rational number, so it adds no round-off errors; this is used for the
# coefficients of the MultipleChoice formulation. With T = Float64 it is fast,
# which suits evaluation. As the inputs themselves usually carry round-off
# errors, e.g. 0.1 or products like 0.1 * 0.3, the result can have tiny
# entries where zero is intended, so entries below 1e-12 times the largest
# entry of each affine function are set to zero.
function _hyperplane_rep(seg::SegmentPointRep{D,F}, ::Type{T}) where {D,F,T}
    if length(seg.input_vals) != D + 1
        error(
            "Every segment of a R^$D -> R^$F function must be a simplex with $(D + 1) points, but a segment has $(length(seg.input_vals)).",
        )
    end
    exact(x) = convert(T, x)
    v_0 = exact.(collect(seg.input_vals[1]))
    M = Matrix{T}(undef, D, D)
    for i in 1:D
        M[:, i] = exact.(collect(seg.input_vals[i+1])) .- v_0
    end
    degenerate =
        T <: AbstractFloat ? LinearAlgebra.rank(M) < D :
        iszero(LinearAlgebra.det(M))
    if degenerate
        error("Degenerate segment with the points $(seg.input_vals).")
    end
    M_inv = inv(M)
    function affine(coeffs, offset)
        a = Float64[coeffs; offset]
        tol = 1e-12 * maximum(abs, a)
        a[abs.(a).<tol] .= 0.0
        return _AffineFunction{D}(ntuple(j -> a[j], D), a[D+1])
    end

    # μ_i(x) = M_inv[i, :] ⋅ (x - v_0) for i = 1, ..., D
    constraints =
        [affine(M_inv[i, :], -LinearAlgebra.dot(M_inv[i, :], v_0)) for i in 1:D]
    # μ_0(x) = 1 - sum of the others
    c = vec(sum(M_inv; dims = 1))
    push!(constraints, affine(-c, 1 + LinearAlgebra.dot(c, v_0)))

    # f_k(x) = y_0 + sum((y_i - y_0) * μ_i(x) for i in 1:D)
    funcs = ntuple(F) do k
        y_0 = exact(seg.output_vals[1][k])
        Δ = [exact(seg.output_vals[i+1][k]) - y_0 for i in 1:D]
        g = transpose(M_inv) * Δ
        return affine(g, y_0 - LinearAlgebra.dot(g, v_0))
    end
    return _HyperplaneSegment{D,F}(constraints, funcs)
end
