# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    GeneralizedCelaya()

The generalized Celaya formulation of univariate continuous piecewise linear
functions (Huchette and Vielma, 2023).

For `d` segments, it uses `ceil(log2(d))` general integer variables with values
in `{-1, 0, 1}`, in addition to the `d + 1` convex multipliers.
"""
struct GeneralizedCelaya <: Method end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::GeneralizedCelaya,
) where {T<:_Scalar}
    _celaya_sos2!(formulation, λ, _generalized_celaya_codes)
    return nothing
end
