struct Incremental <: Method end

function _formulate!(model::JuMP.Model, input_vars::Tuple{VarOrAff}, output_vars::NTuple{F,VarOrAff}, pwl::UnivariatePWLFunction, method::Incremental, direction::DIRECTION) where {F}
    xs, ys = _continuous_breakpoints_or_die(pwl)

    counter = model.ext[:PWL].counter
    n = length(pwl.segments) + 1
    @assert length(xs) == length(ys) == n

    δ = JuMP.@variable(model, [1:n], lower_bound=0, upper_bound=1, base_name="δ_$counter")
    z = JuMP.@variable(model, [1:n-1], Bin, base_name="z_$counter")
    JuMP.@constraint(model, input_vars[1] == xs[1][1] + sum(δ[i]*(xs[i+1][1] - xs[i][1]) for i in 1:n-1))
    for j in 1:F
        rhs = ys[1][j] + sum(δ[i]*(ys[i+1][j] - ys[i][j]) for i in 1:n-1)
        _constrain_output_var(model, output_vars[j], rhs, direction)
    end
    for i in 1:n-1
        JuMP.@constraint(model, δ[i+1] ≤ z[i])
        JuMP.@constraint(model, z[i] ≤ δ[i])
    end
    return
end
