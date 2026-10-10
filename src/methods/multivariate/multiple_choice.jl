# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct MultipleChoice <: Method end

# Each segment gets its own copy of the input, so the function need not be
# continuous.
function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{D,_Scalar},
    output_vars::NTuple{F,_Scalar},
    pwl::PWLFunctionPointRep{D,F},
    method::MultipleChoice,
    direction::DIRECTION,
) where {D,F}
    segments = [_hyperplane_rep(seg) for seg in pwl.segments]
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
    segments::Vector{SegmentHyperplaneRep{D,F}},
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
function _homogenize(f::AffineFunction{D}, x, z) where {D}
    return sum(f.coeffs[j] * x[j] for j in 1:D) + f.offset * z
end

# Converts a simplex with vertices v_0, ..., v_D into the hyperplane
# representation given by its barycentric coordinates μ_0, ..., μ_D ≥ 0.
function _hyperplane_rep(seg::SegmentPointRep{D,F}) where {D,F}
    if length(seg.input_vals) != D + 1
        error(
            "The MultipleChoice method requires every segment of a R^$D -> R^$F function to be a simplex with $(D + 1) points.",
        )
    end
    v_0 = collect(seg.input_vals[1])
    M = Matrix{Float64}(undef, D, D)
    for i in 1:D
        M[:, i] = collect(seg.input_vals[i+1]) .- v_0
    end
    if LinearAlgebra.rank(M) < D
        error(
            "The MultipleChoice method does not support degenerate segments: $(seg.input_vals)",
        )
    end
    M_inv = inv(M)

    # μ_i(x) = M_inv[i, :] ⋅ (x - v_0) for i = 1, ..., D
    constraints = [
        AffineFunction{D}(
            ntuple(j -> M_inv[i, j], D),
            -LinearAlgebra.dot(M_inv[i, :], v_0),
        ) for i in 1:D
    ]
    # μ_0(x) = 1 - sum of the others
    c = vec(sum(M_inv; dims = 1))
    push!(
        constraints,
        AffineFunction{D}(ntuple(j -> -c[j], D), 1 + LinearAlgebra.dot(c, v_0)),
    )

    # f_k(x) = y_0 + sum((y_i - y_0) * μ_i(x) for i in 1:D)
    funcs = ntuple(F) do k
        y_0 = seg.output_vals[1][k]
        Δ = [seg.output_vals[i+1][k] - y_0 for i in 1:D]
        g = transpose(M_inv) * Δ
        return AffineFunction{D}(
            ntuple(j -> g[j], D),
            y_0 - LinearAlgebra.dot(g, v_0),
        )
    end
    return SegmentHyperplaneRep{D,F}(constraints, funcs)
end
