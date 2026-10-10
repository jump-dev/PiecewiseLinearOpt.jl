# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct SymmetricCelaya <: Method end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::SymmetricCelaya,
) where {T<:_Scalar}
    _celaya_sos2!(formulation, λ, _symmetric_celaya_codes)
    return nothing
end
