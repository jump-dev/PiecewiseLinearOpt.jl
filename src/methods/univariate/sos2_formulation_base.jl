# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

function _sos2_encoding_constraints!(
    formulation::PWLFormulation,
    λ::Vector{T},
    y::Vector{<:JuMP.AbstractVariableRef},
    h::Vector{Vector{Float64}},
    B::Vector{Vector{Float64}},
) where {T<:_Scalar}
    model = formulation.model
    n = length(λ) - 1
    for b in B
        record!(
            formulation,
            JuMP.@constraint(
                model,
                LinearAlgebra.dot(b, h[1]) * λ[1] +
                sum(
                    min(
                        LinearAlgebra.dot(b, h[v]),
                        LinearAlgebra.dot(b, h[v-1]),
                    ) * λ[v] for v in 2:n
                ) +
                LinearAlgebra.dot(b, h[n]) * λ[n+1] ≤ LinearAlgebra.dot(b, y)
            ),
        )
        record!(
            formulation,
            JuMP.@constraint(
                model,
                LinearAlgebra.dot(b, h[1]) * λ[1] +
                sum(
                    max(
                        LinearAlgebra.dot(b, h[v]),
                        LinearAlgebra.dot(b, h[v-1]),
                    ) * λ[v] for v in 2:n
                ) +
                LinearAlgebra.dot(b, h[n]) * λ[n+1] ≥ LinearAlgebra.dot(b, y)
            ),
        )
    end
    return nothing
end

function _reflected_gray_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _reflected_gray_codes(k - 1)
        return vcat(
            [vcat(code, 0.0) for code in codes],
            [vcat(code, 1.0) for code in reverse(codes)],
        )
    end
end

function _zigzag_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _zigzag_codes(k - 1)
        return vcat(
            [vcat(code, 0.0) for code in codes],
            [vcat(code, 1.0) for code in codes],
        )
    end
end

function _integer_zigzag_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _integer_zigzag_codes(k - 1)
        offset = [2^(j - 2) for j in k:-1:2]
        return vcat(
            [vcat(code, 0.0) for code in codes],
            [vcat(code .+ offset, 1.0) for code in codes],
        )
    end
end

function _zigzag_hyperplanes(k::Int)::Vector{Vector{Float64}}
    hps = Vector{Int}[]
    for i in 1:k
        hp = zeros(Int, k)
        hp[i] = 1
        for j in (i+1):k
            hp[j] = 2^(j - i - 1)
        end
        push!(hps, hp)
    end
    return hps
end

function _unit_vector_hyperplanes(k::Int)::Vector{Vector{Float64}}
    hps = Vector{Float64}[]
    for i in 1:k
        hp = zeros(Int, k)
        hp[i] = 1
        push!(hps, hp)
    end
    return hps
end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::Tuple{_Scalar},
    output_vars::NTuple{F,_Scalar},
    pwl::PWLFunction{1,F},
    method::SOS2Method,
    direction::Direction,
) where {F}
    grid = _continuous_gridpoints_or_die(pwl)
    λ = _create_convex_multiplier_vars(
        formulation,
        grid,
        input_vars,
        output_vars,
        direction,
    )
    formulate_sos2!(formulation, λ, method)
    return nothing
end

# Any univariate method gives an SOS2 formulation, by formulating the function
# whose outputs are the unit vectors at the breakpoints.
# The arguments other than the method are untyped, so that any implementation
# for a specific method is more specific than this fallback.
function formulate_sos2!(formulation, λ, method::Method)
    n = length(λ)
    D = 1
    F = n
    d = range(0, 1; length = n)
    segments = SegmentPointRep{D,F}[]
    for i in 1:(n-1)
        output_left = ntuple(t -> t == i ? 1.0 : 0.0, n)
        output_right = ntuple(t -> t == i + 1 ? 1.0 : 0.0, n)
        push!(
            segments,
            SegmentPointRep{D,F}(
                [(d[i],), (d[i+1],)],
                [output_left, output_right],
            ),
        )
    end
    dummy_input_var = record!(
        formulation,
        JuMP.@variable(
            formulation.model,
            lower_bound = 0,
            upper_bound = 1,
            base_name = variable_name(formulation, "t")
        ),
    )
    dummy_pwl = PWLFunction(segments, Intervals())
    formulate_pwl!(
        formulation,
        (dummy_input_var,),
        tuple(λ...),
        dummy_pwl,
        method,
        Graph(),
    )
    return nothing
end
