struct ZigZagBinary <: Method end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::ZigZagBinary) where {T <: VarOrAff}
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
    y = JuMP.@variable(model, [1:k], Bin, base_name = "y_$counter")
    _sos2_encoding_constraints!(model, λ, y, _zigzag_codes(k), _zigzag_hyperplanes(k))
    return nothing
end
