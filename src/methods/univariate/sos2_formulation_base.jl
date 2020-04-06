const SOS2Method = Union{LogarithmicEmbedding, LogarithmicIndependentBranching, NativeSOS2, ZigZagBinary, ZigZagInteger}

function _sos2_encoding_constraints!(m::JuMP.Model, λ::Vector{T}, y::Vector{JuMP.VariableRef}, h::Vector{Vector{Float64}}, B::Vector{Vector{Float64}}) where {T <: VarOrAff}
    n = length(λ) - 1
    for b in B
        JuMP.@constraints(m, begin
            dot(b, h[1]) * λ[1] + sum(min(dot(b, h[v]), dot(b, h[v-1])) * λ[v] for v in 2:n) + dot(b, h[n]) * λ[n+1] ≤ dot(b, y)
            dot(b, h[1]) * λ[1] + sum(max(dot(b, h[v]), dot(b, h[v-1])) * λ[v] for v in 2:n) + dot(b, h[n]) * λ[n+1] ≥ dot(b, y)
        end)
    end
    return nothing
end

function _reflected_gray_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _reflected_gray_codes(k-1)
        return vcat([vcat(code, 0.0) for code in codes],
                    [vcat(code, 1.0) for code in reverse(codes)])
    end
end

function _zigzag_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _zigzag_codes(k-1)
        return vcat([vcat(code, 0.0) for code in codes],
                    [vcat(code, 1.0) for code in codes])
    end
end

function _integer_zigzag_codes(k::Int)::Vector{Vector{Float64}}
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes = _integer_zigzag_codes(k-1)
        offset = [2^(j - 2) for j in k:-1:2]
        return vcat([vcat(code,           0.0) for code in codes],
                    [vcat(code .+ offset, 1.0) for code in codes])
    end
end

function _zigzag_hyperplanes(k::Int)::Vector{Vector{Float64}}
    hps = Vector{Int}[]
    for i in 1:k
        hp = zeros(Int, k)
        hp[i] = 1
        for j in (i+1):k
            hp[j] = 2^(j-i-1)
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

function formulate_pwl!(model::JuMP.Model, input_vars::Tuple{VarOrAff}, output_vars::NTuple{F,VarOrAff}, pwl::UnivariatePWLFunction{F}, method::SOS2Method, direction::DIRECTION) where {F}
    grid = _continuous_gridpoints_or_die(pwl)
    λ = _create_convex_multiplier_vars(model, grid, input_vars, output_vars, direction)
    formulate_sos2!(model, λ, method)
    return nothing
end
