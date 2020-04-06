struct LogarithmicIndependentBranching <: Method end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::LogarithmicIndependentBranching) where {T <: VarOrAff}
    counter = model.ext[:PWL].counter
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    z = JuMP.@variable(model, [1:k], Bin, base_name = "y_$counter")
    _H = _reflected_gray_codes(k)
    H = Dict(i => _H[i] for i in 1:d)
    H[0] = H[1]
    H[d+1] = H[d]
    for j in 1:k
        JuMP.@constraints(model, begin
            sum(λ[i] for i in 1:n if H[i-1][j] == H[i][j] == 1) ≤     z[j]
            sum(λ[i] for i in 1:n if H[i-1][j] == H[i][j] == 0) ≤ 1 - z[j]
        end)
    end
    return nothing
end
