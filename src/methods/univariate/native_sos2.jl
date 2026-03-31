# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    NativeSOS2()

Delegates the SOS2 constraint directly to the solver using JuMP's built-in
SOS2 support. No additional binary variables are introduced.

Requires a solver that supports SOS2 constraints natively.
"""
struct NativeSOS2 <: Method end

function formulate_sos2!(
    model::JuMP.Model,
    λ::Vector{T},
    method::NativeSOS2,
) where {T<:VarOrAff}
    JuMP.@constraint(model, λ in JuMP.SOS2([k for k in 1:length(λ)]))
    return nothing
end
