# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    LogarithmicIndependentBranching()

The logarithmic independent branching formulation of univariate continuous
piecewise linear functions (Vielma and Nemhauser, 2011), based on the binary
reflected Gray code.

For `d` segments, it uses `ceil(log2(d))` binary variables and
`2 * ceil(log2(d))` constraints in addition to the `d + 1` convex multipliers.
See also [`IndependentBranching`](@ref).
"""
struct LogarithmicIndependentBranching <: SOS2Method end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::LogarithmicIndependentBranching,
) where {T<:Scalar}
    model = formulation.model
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    if k == 0
        return nothing
    end
    z = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:k],
            Bin,
            base_name = variable_name(formulation, "z")
        ),
    )
    _H = _reflected_gray_codes(k)
    H = Dict(i => _H[i] for i in 1:d)
    H[0] = H[1]
    H[d+1] = H[d]
    for j in 1:k
        record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(λ[i] for i in 1:n if H[i-1][j] == H[i][j] == 1) ≤ z[j]
            ),
        )
        record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(λ[i] for i in 1:n if H[i-1][j] == H[i][j] == 0) ≤ 1 - z[j]
            ),
        )
    end
    return nothing
end
