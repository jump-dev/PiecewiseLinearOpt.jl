# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    DisaggregatedLogarithmic()

The logarithmic disaggregated convex combination formulation.

It works for piecewise linear functions of any dimension (see Vielma, Ahmed and
Nemhauser, 2010), with a convex multiplier for each vertex of each segment. As
the segments do not share multipliers, the function may be discontinuous.

For `d` segments, it uses `ceil(log2(d))` binary variables.
"""
struct DisaggregatedLogarithmic <: Method end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{D,_Scalar},
    output_vars::NTuple{F,_Scalar},
    pwl::PWLFunction{D,F},
    method::DisaggregatedLogarithmic,
    direction::Direction,
) where {D,F}
    model = formulation.model
    segments = pwl.segments
    num_bps = Dict(seg => length(seg.input_vals) for seg in segments)

    γ = record!(
        formulation,
        JuMP.@variable(
            model,
            [seg in segments, i in 1:num_bps[seg]],
            lower_bound = 0,
            upper_bound = 1,
            base_name = variable_name(formulation, "γ")
        ),
    )
    record!(formulation, JuMP.@constraint(model, sum(γ) == 1))

    for j in 1:D
        record!(
            formulation,
            JuMP.@constraint(
                model,
                input_vars[j] == sum(
                    sum(
                        γ[seg, i] * seg.input_vals[i][j] for i in 1:num_bps[seg]
                    ) for seg in segments
                )
            ),
        )
    end
    for j in 1:F
        rhs = sum(
            sum(γ[seg, i] * seg.output_vals[i][j] for i in 1:num_bps[seg])
            for seg in segments
        )
        constrain_output!(formulation, output_vars[j], rhs, direction)
    end

    r = ceil(Int, log2(length(segments)))
    if r == 0
        return nothing
    end
    _H = _reflected_gray_codes(r)
    H = Dict(segments[i] => _H[i] for i in 1:length(segments))
    z = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:r],
            Bin,
            base_name = variable_name(formulation, "z")
        ),
    )
    for j in 1:r
        record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(
                    sum(γ[seg, i] * H[seg][j] for i in 1:num_bps[seg]) for
                    seg in segments
                ) == z[j]
            ),
        )
    end
    return nothing
end
