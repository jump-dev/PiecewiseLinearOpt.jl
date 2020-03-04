struct Incremental <: Method end

function _formulate!(model::JuMP.Model, input_var::VarOrAff, output_vars::NTuple{F,VarOrAff}, pwl::UnivariatePWLFunction, method::Incremental, direction::DIRECTION) where {F}

    d, fds = _continuous_breakpoints_or_die(pwl)

    counter = model.ext[:PWL].counter
    n = length(pwl.segments)

    δ = JuMP.@variable(model, [1:n], lower_bound=0, upper_bound=1, base_name="δ_$counter")
    z = JuMP.@variable(model, [1:n-1], Bin, base_name="z_$counter")
    JuMP.@constraint(model, input_var == d[1] + sum(δ[i]*( d[i+1]- d[i]) for i in 1:n-1))
    for j in 1:F
        rhs = fds[j][1] + sum(δ[i]*(fd[j][i+1]-fd[j][i]) for i in 1:n-1)
        _constrain_output_var(model, output_vars[j], rhs, direction)
    end
    for i in 1:n-1
        JuMP.@constraint(model, δ[i+1] ≤ z[i])
        JuMP.@constraint(model, z[i] ≤ δ[i])
    end
    return
end
