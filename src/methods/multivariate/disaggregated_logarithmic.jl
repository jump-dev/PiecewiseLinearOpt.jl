# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    DisaggregatedLogarithmic()

The disaggregated logarithmic formulation for piecewise linear functions.

Uses `⌈log₂(S)⌉` binary variables where `S` is the number of segments.
Works for any input dimension. More compact than [`ConvexCombination`](@ref)
for functions with many segments.
"""
struct DisaggregatedLogarithmic <: Method end

function formulate_pwl!(
    model::JuMP.Model,
    input_vars::NTuple{D,VarOrAff},
    output_vars::NTuple{F,VarOrAff},
    pwl::PWLFunctionPointRep{D,F},
    method::DisaggregatedLogarithmic,
    direction::DIRECTION,
) where {D,F}
    segments = pwl.segments
    num_bps = Dict(seg => length(seg.input_vals) for seg in segments)

    γ = JuMP.@variable(
        model,
        [seg in segments, i in 1:num_bps[seg]],
        lower_bound = 0,
        upper_bound = 1,
        base_name = _pwl_name(model, "γ")
    )
    JuMP.@constraint(model, sum(γ) == 1)

    for j in 1:D
        JuMP.@constraint(
            model,
            input_vars[j] == sum(
                sum(γ[seg, i] * seg.input_vals[i][j] for i in 1:num_bps[seg])
                for seg in segments
            )
        )
    end
    for j in 1:F
        rhs = sum(
            sum(γ[seg, i] * seg.output_vals[i][j] for i in 1:num_bps[seg])
            for seg in segments
        )
        _constrain_output_var(model, output_vars[j], rhs, direction)
    end

    r = ceil(Int, log2(length(segments)))
    if r == 0
        return nothing
    end
    _H = _reflected_gray_codes(r)
    H = Dict(segments[i] => _H[i] for i in 1:length(segments))
    z = JuMP.@variable(model, [1:r], Bin, base_name = _pwl_name(model, "z"))
    for j in 1:r
        JuMP.@constraint(
            model,
            sum(
                sum(γ[seg, i] * H[seg][j] for i in 1:num_bps[seg]) for
                seg in segments
            ) == z[j]
        )
    end
end
