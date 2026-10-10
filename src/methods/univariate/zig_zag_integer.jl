# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ZigZagInteger()

The general integer zig-zag formulation of univariate continuous piecewise
linear functions (Huchette and Vielma, 2023).

For `d` segments, it uses `ceil(log2(d))` bounded general integer variables and
`2 * ceil(log2(d))` constraints in addition to the `d + 1` convex multipliers.
"""
struct ZigZagInteger <: Method end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::ZigZagInteger,
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
    y = _record!(
        formulation,
        JuMP.@variable(
            formulation.model,
            [i in 1:k],
            Int,
            lower_bound = 0,
            upper_bound = 2^(k - i),
            base_name = _name(formulation, "z")
        ),
    )
    _sos2_encoding_constraints!(
        formulation,
        λ,
        y,
        _integer_zigzag_codes(k),
        _unit_vector_hyperplanes(k),
    )
    return nothing
end
