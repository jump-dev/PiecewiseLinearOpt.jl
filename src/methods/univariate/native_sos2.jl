# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct NativeSOS2 <: Method end

function formulate_sos2!(
    model::JuMP.Model,
    λ::Vector{T},
    method::NativeSOS2,
) where {T<:VarOrAff}
    JuMP.@constraint(model, λ in JuMP.SOS2([k for k in 1:length(λ)]))
    return nothing
end
