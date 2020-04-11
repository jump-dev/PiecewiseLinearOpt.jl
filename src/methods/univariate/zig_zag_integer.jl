struct ZigZagInteger <: Method end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::ZigZagInteger) where {T <: VarOrAff}
    counter = model.ext[:PWL].counter
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    if k == 0
        return nothing
    end
    y = JuMP.@variable(model, [i in 1:k], Int, lower_bound = 0, upper_bound = 2^(k - i), base_name = "y_$counter")
    _sos2_encoding_constraints!(model, λ, y, _integer_zigzag_codes(k), _unit_vector_hyperplanes(k))
    return nothing
end
