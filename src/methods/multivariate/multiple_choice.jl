struct MultipleChoice <: Method end

function formulate_pwl!(model::JuMP.Model, input_vars::Tuple{VarOrAff}, output_vars::NTuple{F, VarOrAff}, pwl::PWLFunctionHyperplaneRep{D, F}, method::MultipleChoice, direction::DIRECTION) where {D, F}
    x_hat = JuMP.@variable(model, [segments, 1:D], base_name="x_hat_$counter")
    y_hat = JuMP.@variable(model, [segments, 1:F], base_name="y_hat_$counter")
    z = JuMP.@variable(model, [segments], Bin, base_name="z_$counter")
    JuMP.@constraint(model, sum(z) == 1)
    for i in 1:D
        JuMP.@constraint(model, sum(x_hat[:, i]) == x[i])
    end
    for i in 1:F
        JuMP.@constraint(model, sum(y_hat[:, i]) == y[i])
    end
    for seg in segments
        for constraint in seg.constraints
            coeffs, offset = constraint.coeffs, constraint.offset
            JuMP.@constraint(model, dot(coeffs, x_hat[seg, :]) + offset * z[seg] ≥ 0)
        end
        for i in 1:F
            output_func = seg.output_funcs[i]
            coeffs, offset = output_func.coeffs, output_func.offset
            JuMP.@constraint(model, y_hat[seg, i] == dot(coeffs, x_hat[seg, :]) + offset * z[seg])
        end
    end
    return nothing
end
