# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    ConvexCombination()

The convex combination formulation of continuous piecewise linear functions.

It works for functions of any dimension (see Vielma, Ahmed and Nemhauser,
2010), with a convex multiplier for each breakpoint and a binary variable for
each segment.

It can also be used as the axis method of the bivariate grid formulations.
"""
struct ConvexCombination <: PWLMethod end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{D,Scalar},
    output_vars::NTuple{F,Scalar},
    pwl::PWLFunction{D,F},
    method::ConvexCombination,
    direction::Direction,
) where {D,F}
    model = formulation.model
    segments = pwl.segments
    z = record!(
        formulation,
        JuMP.@variable(
            model,
            [segments],
            Bin,
            base_name = variable_name(formulation, "z")
        ),
    )
    record!(formulation, JuMP.@constraint(model, sum(z) == 1))

    all_input_vals = Set{NTuple{D,Float64}}()
    output_val_map = Dict{NTuple{D,Float64},NTuple{F,Float64}}()
    for seg in segments
        for (it, input_val) in enumerate(seg.input_vals)
            push!(all_input_vals, input_val)
            output_val = seg.output_vals[it]
            if haskey(output_val_map, input_val)
                if output_val_map[input_val] != output_val
                    error(
                        "The ConvexCombination method does not currently support discontinuous piecewise linear functions.",
                    )
                end
            else
                output_val_map[input_val] = output_val
            end
        end
    end
    λ = record!(
        formulation,
        JuMP.@variable(
            model,
            [all_input_vals],
            lower_bound = 0,
            upper_bound = 1,
            base_name = variable_name(formulation, "λ")
        ),
    )
    record!(formulation, JuMP.@constraint(model, sum(λ) == 1))
    for j in 1:D
        record!(
            formulation,
            JuMP.@constraint(
                model,
                input_vars[j] == sum(
                    λ[input_val] * input_val[j] for input_val in all_input_vals
                )
            ),
        )
    end
    for j in 1:F
        rhs = sum(
            λ[input_val] * output_val_map[input_val][j] for
            input_val in all_input_vals
        )
        constrain_output!(formulation, output_vars[j], rhs, direction)
    end
    for input_val in all_input_vals
        record!(
            formulation,
            JuMP.@constraint(
                model,
                λ[input_val] ≤
                sum(z[seg] for seg in segments if input_val in seg.input_vals)
            ),
        )
    end
    return nothing
end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::ConvexCombination,
) where {T<:Scalar}
    model = formulation.model
    n = length(λ)
    z = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:n-1],
            Bin,
            base_name = variable_name(formulation, "z")
        ),
    )
    record!(formulation, JuMP.@constraint(model, sum(z) == 1))
    record!(formulation, JuMP.@constraint(model, λ[1] ≤ z[1]))
    for i in 2:n-1
        record!(formulation, JuMP.@constraint(model, λ[i] ≤ z[i-1] + z[i]))
    end
    record!(formulation, JuMP.@constraint(model, λ[n] ≤ z[n-1]))
    return nothing
end
