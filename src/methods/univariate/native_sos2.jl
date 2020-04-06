struct NativeSOS2 <: Method end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::NativeSOS2) where {T <: VarOrAff}
    JuMP.@constraint(model, λ in MOI.SOS2([k for k in 1:length(λ)]))
    return nothing
end
