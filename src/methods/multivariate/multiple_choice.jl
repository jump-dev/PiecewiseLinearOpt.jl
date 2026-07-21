# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    MultipleChoice()

The multiple choice formulation for piecewise linear functions in hyperplane
representation.

This is the only method that supports `PWLFunction` objects with
[`SegmentHyperplaneRep`](@ref PiecewiseLinearOpt.SegmentHyperplaneRep) segments,
where each segment is defined by affine constraints and affine output functions.
"""
struct MultipleChoice <: Method end

function formulate_pwl!(
    model::JuMP.Model,
    input_vars::NTuple{D,VarOrAff},
    output_vars::NTuple{F,VarOrAff},
    pwl::PWLFunctionHyperplaneRep{D,F},
    method::MultipleChoice,
    direction::DIRECTION,
) where {D,F}
    segments = pwl.segments
    S = 1:length(segments)
    x_hat =
        JuMP.@variable(model, [S, 1:D], base_name = _pwl_name(model, "x_hat"))
    y_hat =
        JuMP.@variable(model, [S, 1:F], base_name = _pwl_name(model, "y_hat"))
    z = JuMP.@variable(model, [S], Bin, base_name = _pwl_name(model, "z"))
    JuMP.@constraint(model, sum(z) == 1)
    for i in 1:D
        JuMP.@constraint(model, sum(x_hat[:, i]) == input_vars[i])
    end
    for i in 1:F
        _constrain_output_var(
            model,
            output_vars[i],
            sum(y_hat[:, i]),
            direction,
        )
    end
    for (s, seg) in enumerate(segments)
        for constraint in seg.constraints
            coeffs, offset = constraint.coeffs, constraint.offset
            JuMP.@constraint(
                model,
                LinearAlgebra.dot(coeffs, x_hat[s, :]) + offset * z[s] ≥ 0
            )
        end
        for i in 1:F
            output_func = seg.funcs[i]
            coeffs, offset = output_func.coeffs, output_func.offset
            JuMP.@constraint(
                model,
                y_hat[s, i] ==
                LinearAlgebra.dot(coeffs, x_hat[s, :]) + offset * z[s]
            )
        end
    end
    return nothing
end
