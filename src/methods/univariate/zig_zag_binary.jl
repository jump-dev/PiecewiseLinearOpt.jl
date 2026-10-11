# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ZigZagBinary()

The binary zig-zag formulation of univariate continuous piecewise linear
functions (Huchette and Vielma, 2023).

For `d` segments, it uses `ceil(log2(d))` binary variables and
`2 * ceil(log2(d))` constraints in addition to the `d + 1` convex multipliers.
"""
struct ZigZagBinary <: SOS2Method end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::ZigZagBinary,
) where {T<:_Scalar}
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    if k == 0
        return nothing
    end
    y = record!(
        formulation,
        JuMP.@variable(
            formulation.model,
            [1:k],
            Bin,
            base_name = variable_name(formulation, "z")
        ),
    )
    _sos2_encoding_constraints!(
        formulation,
        λ,
        y,
        _zigzag_codes(k),
        _zigzag_hyperplanes(k),
    )
    return nothing
end
