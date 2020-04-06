struct ZigZagInteger <: Method end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::ZigZagInteger) where {T <: VarOrAff}
    counter = model.ext[:PWL].counter
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    codes = _integer_zigzag_codes(k)
    lb = [minimum(t[i] for t in codes) for i in 1:k]
    ub = [maximum(t[i] for t in codes) for i in 1:k]
    y = JuMP.@variable(model, [i=1:k], Int, lower_bound = lb[i], upper_bound = ub[i], base_name = "y_$counter")
    _sos2_encoding_constraints!(model, λ, y, codes, _unit_vector_hyperplanes(k))
    return nothing
end
