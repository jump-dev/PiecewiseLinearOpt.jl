# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    LogarithmicEmbedding()

The logarithmic embedding formulation for SOS2 constraints.

Uses `⌈log₂(n-1)⌉` binary variables for `n` breakpoints, based on reflected
Gray codes. This is the default method for univariate piecewise linear functions.

`Logarithmic` is an alias for `LogarithmicEmbedding`.
"""
struct LogarithmicEmbedding <: Method end

function formulate_sos2!(
    model::JuMP.Model,
    λ::Vector{T},
    method::LogarithmicEmbedding,
) where {T<:VarOrAff}
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    if k == 0
        return nothing
    end
    y = JuMP.@variable(model, [1:k], Bin, base_name = _pwl_name(model, "y"))
    _sos2_encoding_constraints!(
        model,
        λ,
        y,
        _reflected_gray_codes(k),
        _unit_vector_hyperplanes(k),
    )
    return nothing
end
