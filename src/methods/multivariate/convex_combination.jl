struct ConvexCombination <: Method end

function formulate_pwl!(model::JuMP.Model, input_vars::NTuple{D, VarOrAff}, output_vars::NTuple{F, VarOrAff}, pwl::PWLFunctionPointRep{D, F}, method::ConvexCombination, direction::DIRECTION) where {D, F}
    # TODO: assert PWL function is continuous
    counter = model.ext[:PWL].counter
    segments = pwl.segments
    z = JuMP.@variable(model, [segments], Bin, base_name = "z_$counter")
    JuMP.@constraint(model, sum(z) == 1)

    all_input_vals = Set{NTuple{D, Float64}}()
    output_val_map = Dict{NTuple{D, Float64}, NTuple{F, Float64}}()
    for seg in segments
        for (it, input_val) in enumerate(seg.input_vals)
            push!(all_input_vals, input_val)
            output_val = seg.output_vals[it]
            if haskey(output_val_map, input_val)
                if output_val_map[input_val] != output_val
                    error("The ConvexCombination method does not currently support discontinuous piecewise linear functions.")
                end
            else
                output_val_map[input_val] = output_val
            end
        end
    end
    λ = JuMP.@variable(model, [all_input_vals], lower_bound = 0, upper_bound = 1, base_name = "λ_$counter")
    JuMP.@constraint(model, sum(λ) == 1)
    for j in 1:D
        JuMP.@constraint(model, input_vars[j] == sum(λ[input_val] * input_val[j] for input_val in all_input_vals))
    end
    for j in 1:F
        rhs = sum(λ[input_val] * output_val_map[input_val][j] for input_val in all_input_vals)
        _constrain_output_var(model, output_vars[j], rhs, direction)
    end
    for input_val in all_input_vals
        JuMP.@constraint(model, λ[input_val] ≤ sum(z[seg] for seg in segments if input_val in seg.input_vals))
    end
    return nothing
end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::ConvexCombination) where {T <: VarOrAff}
    counter = model.ext[:PWL].counter
    n = length(λ)
    z = JuMP.@variable(model, [1:n-1], Bin, base_name="z_$counter")
    JuMP.@constraint(model, sum(z) == 1)
    JuMP.@constraint(model, λ[1] ≤ z[1])
    for i in 2:n-1
        JuMP.@constraint(model, λ[i] ≤ z[i-1] + z[i])
    end
    JuMP.@constraint(model, λ[n] ≤ z[n-1])
    return nothing
end
