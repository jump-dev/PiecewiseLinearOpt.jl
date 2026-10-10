# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

function _sos2_encoding_constraints!(
    m::JuMP.Model,
    λ::Vector{T},
    y::Vector{JuMP.VariableRef},
    h::Vector{Vector{Float64}},
    B::Vector{Vector{Float64}},
) where {T<:VarOrAff}
    n = length(λ) - 1
    for b in B
        JuMP.@constraints(
            m,
            begin
                LinearAlgebra.dot(b, h[1]) * λ[1] +
                sum(
                    min(
                        LinearAlgebra.dot(b, h[v]),
                        LinearAlgebra.dot(b, h[v-1]),
                    ) * λ[v] for v in 2:n
                ) +
                LinearAlgebra.dot(b, h[n]) * λ[n+1] ≤ LinearAlgebra.dot(b, y)
                LinearAlgebra.dot(b, h[1]) * λ[1] +
                sum(
                    max(
                        LinearAlgebra.dot(b, h[v]),
                        LinearAlgebra.dot(b, h[v-1]),
                    ) * λ[v] for v in 2:n
                ) +
                LinearAlgebra.dot(b, h[n]) * λ[n+1] ≥ LinearAlgebra.dot(b, y)
            end
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

function _generalized_celaya_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    elseif k == 2
        return [[0.0, 0.0], [1.0, 1.0], [1.0, 0.0], [0.0, 1.0]]
    else
        codes = _generalized_celaya_codes(k - 1)
        n = length(codes)
        hp = div(n, 2)
        return vcat(
            [vcat(code, 0.0) for code in codes],
            [vcat(codes[i], 1.0) for i in 1:hp],
            [vcat(codes[i], -1.0) for i in (hp+1):n],
        )
    end
end

function _symmetric_celaya_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _generalized_celaya_codes(k - 1)
        n = length(codes)
        hp = div(n, 2)
        return vcat(
            [vcat(code, 0.0) for code in codes],
            [vcat(codes[i], 1.0) for i in 1:hp],
            [vcat(codes[i], -1.0) for i in (hp+1):n],
        )
    end
end

function _celaya_sos2!(
    model::JuMP.Model,
    λ::Vector{T},
    codes::Function,
) where {T<:VarOrAff}
    counter = model.ext[:PWL].counter
    n = length(λ)
    d = n - 1
    if 0 <= d <= 1
        return nothing
    end
    k = ceil(Int, log2(d))
    C = codes(k)
    # Hyperplanes can only be computed from codes of length at least 2; with a
    # single bit the codes are binary, as in the logarithmic formulation.
    H = k == 1 ? _unit_vector_hyperplanes(k) : _compute_hyperplanes(C)
    lb = [minimum(c[i] for c in C) for i in 1:k]
    ub = [maximum(c[i] for c in C) for i in 1:k]
    y = JuMP.@variable(
        model,
        [i in 1:k],
        Int,
        lower_bound = lb[i],
        upper_bound = ub[i],
        base_name = "y_$counter"
    )
    _sos2_encoding_constraints!(model, λ, y, C, H)
    return nothing
end

# Computes the hyperplanes spanned by the differences between consecutive codes
function _compute_hyperplanes(C::Vector{Vector{T}}) where {T<:Number}
    n = length(C)
    k = length(C[1])

    d = zeros(n - 1, k)
    for i in 1:n-1
        d[i, :] = C[i+1] - C[i]
    end

    if k <= 1
        error("Cannot process codes of length $k")
    elseif k == 2
        @assert n == 4
        spanners = Vector{Float64}[]
        for i in 1:n-1
            v = _canonical!(d[i, :])
            v = [v[2], -v[1]]
            push!(spanners, v)
        end
        indices = [1]
        approx12 = isapprox(spanners[1], spanners[2])
        approx13 = isapprox(spanners[1], spanners[3])
        approx23 = isapprox(spanners[2], spanners[3])
        if !approx12
            push!(indices, 2)
        end
        if !approx13 && !(!approx12 && approx23)
            push!(indices, 3)
        end
        return spanners[indices]
    end

    indices = [1, 2]
    spanners = Vector{Float64}[]
    while !isempty(indices)
        if indices == [n - 1]
            break
        end
        if LinearAlgebra.rank(d[indices, :]) == length(indices) &&
           length(indices) <= k - 1
            if length(indices) == k - 1
                nullsp = LinearAlgebra.nullspace(d[indices, :])
                @assert size(nullsp, 2) == 1
                v = vec(nullsp)
                push!(spanners, _canonical!(v))
            end
            if indices[end] != n - 1
                push!(indices, indices[end] + 1)
            else
                pop!(indices)
                indices[end] += 1
            end
        else
            if indices[end] != n - 1
                indices[end] += 1
            else
                pop!(indices)
                indices[end] += 1
            end
        end
    end

    keepers = [1]
    for i in 2:length(spanners)
        if !any(isapprox(spanners[i], spanners[j]) for j in keepers)
            push!(keepers, i)
        end
    end

    return spanners[keepers]
end

# Normalizes `v`, rounds near-zero entries to zero, and makes the first nonzero
# entry positive.
function _canonical!(v::Vector{Float64})
    LinearAlgebra.normalize!(v)
    for j in 1:length(v)
        if abs(v[j]) < 1e-8
            v[j] = 0
        end
    end
    sgn = sign(v[findfirst(!iszero, v)])
    v .*= sgn
    return v
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

# ConvexCombination could also be added to this list, but the multivariate
# implementation works just fine.
const SOS2Method = Union{
    GeneralizedCelaya,
    LogarithmicEmbedding,
    LogarithmicIndependentBranching,
    NativeSOS2,
    SymmetricCelaya,
    ZigZagBinary,
    ZigZagInteger,
}

function formulate_pwl!(
    model::JuMP.Model,
    input_vars::Tuple{VarOrAff},
    output_vars::NTuple{F,VarOrAff},
    pwl::PWLFunctionPointRep{1,F},
    method::SOS2Method,
    direction::DIRECTION,
) where {F}
    grid = _continuous_gridpoints_or_die(pwl)
    λ = _create_convex_multiplier_vars(
        model,
        grid,
        input_vars,
        output_vars,
        direction,
    )
    formulate_sos2!(model, λ, method)
    return nothing
end

function formulate_sos2!(
    model::JuMP.Model,
    λ::Vector{T},
    method::Method,
) where {T<:VarOrAff}
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
    dummy_input_var = JuMP.@variable(model, lower_bound = 0, upper_bound = 1)
    dummy_pwl = PWLFunctionPointRep{1,F}(segments, Intervals())
    formulate_pwl!(
        model,
        (dummy_input_var,),
        tuple(λ...),
        dummy_pwl,
        method,
        Graph,
    )
    return nothing
end
