function _continuous_breakpoints_or_die(pwl::UnivariatePWLFunction{F}) where {F}
    xs = Set{Tuple{Float64}}()
    for segment in pwl.segments
        for val in segment.input_vals
            push!(xs, val)
        end
    end
    ys = Dict{Tuple{Float64},NTuple{F,Float64}}()
    for segment in pwl.segments
        for i in 1:length(segment.input_vals)
            x = segment.input_vals[i]
            y = segment.output_vals[i]
            if haskey(ys, x)
                if y != ys[x]
                    error("Expected piecewise linear function to be continuous.")
                end
            else
                ys[x] = y
            end
        end
    end
    X = sort!(collect(xs))
    Y = [ys[x] for x in X]
    return X, Y
end

function _constrain_output_var(model::JuMP.Model, output_var::VarOrAff, func_val::VarOrAff, direction::DIRECTION)
    if direction == Graph
        JuMP.@constraint(model, output_var == func_val)
    elseif direction == Hypograph
        JuMP.@constraint(model, output_var <= func_val)
    else
        @assert direction == Epigraph
        JuMP.@constraint(model, output_var >= func_val)
    end
    return
end

function _create_convex_multiplier_vars(model::JuMP.Model, xs::Vector{NTuple{D,Float64}}, ys::Vector{NTuple{F,Float64}}, input_vars::NTuple{D,JuMP.VariableRef}, output_vars::NTuple{F,JuMP.VariableRef}, direction::DIRECTION) where {D,F}
    counter = model.ext[:PWL].counter
    n = length(xs)
    @assert n == length(ys)

    λ = JuMP.@variable(model, [1:n], lower_bound=0, upper_bound=1, base_name="λ_$counter")
    JuMP.@constraint(model, sum(λ) == 1)
    for j in 1:D
        JuMP.@constraint(model, sum(λ[i]* xs[i][j] for i in 1:n) == input_vars[j])
    end
    for j in 1:F
        _constrain_output_var(model, output_vars[j], sum(λ[i]*ys[i][j] for i in 1:n), direction)
    end
    return λ
end

function _sos2_encoding_constraints!(m::JuMP.Model, λ::Vector{JuMP.VariableRef}, y::Vector{JuMP.VariableRef}, h::Vector{Vector{Float64}}, B::Vector{Vector{Float64}})
    n = length(λ) - 1
    for b in B
        JuMP.@constraints(m, begin
            dot(b,h[1])*λ[1] + sum(min(dot(b,h[v]),dot(b,h[v-1]))*λ[v] for v in 2:n) + dot(b,h[n])*λ[n+1] ≤ dot(b,y)
            dot(b,h[1])*λ[1] + sum(max(dot(b,h[v]),dot(b,h[v-1]))*λ[v] for v in 2:n) + dot(b,h[n])*λ[n+1] ≥ dot(b,y)
        end)
    end
    return nothing
end

function _zigzag_codes(k::Int)
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0],[1]]
    else
        codes′ = _zigzag_codes(k-1)
        return vcat([vcat(code,0) for code in codes′],
                    [vcat(code,1) for code in codes′])
    end
end

function _integer_zigzag_codes(k::Int)
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0],[1]]
    else
        codes′ = _integer_zigzag_codes(k-1)
        offset = [2^(j-2) for j in k:-1:2]
        return vcat([vcat(code,        0) for code in codes′],
                    [vcat(code.+offset,1) for code in codes′])
    end
end

function _zigzag_hyperplanes(k::Int)
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

function _unit_vector_hyperplanes(k::Int)
    hps = Vector{Float64}[]
    for i in 1:k
        hp = zeros(Int,k)
        hp[i] = 1
        push!(hps, hp)
    end
    return hps
end

_generalized_celaya_hyperplanes(k::Int) = _compute_hyperplanes(_generalized_celaya_codes(k))

_symmetric_celaya_hyperplanes(k::Int) = _compute_hyperplanes(_symmetric_celaya_codes(k))

# TODO: Figure out assumptions at play here
function _compute_hyperplanes(C::Vector{Vector{T}}) where T <: Number
    n = length(C)
    k = length(C[1])

    d = zeros(n-1,k)
    for i in 1:n-1
        d[i,:] = C[i+1] - C[i]
    end

    if k <= 1
        error("Cannot process codes of length $k")
    elseif k == 2
        @assert n == 4
        spanners = Vector{Float64}[]
        for i in 1:n-1
            v = canonical!(d[i,:])
            v = [v[2], -v[1]]
            push!(spanners, v)
        end
        indices = [1]
        approx12 = isapprox(spanners[1],spanners[2])
        approx13 = isapprox(spanners[1],spanners[3])
        approx23 = isapprox(spanners[2],spanners[3])
        if !approx12
            push!(indices,2)
        end
        if !approx13 && !(!approx12 && approx23)
            push!(indices,3)
        end
        return spanners[indices]
    end

    indices = [1,2]
    spanners = Vector{Float64}[]
    while !isempty(indices)
        if indices == [n-1]
            break
        end
        if rank(d[indices,:]) == length(indices) && length(indices) <= k-1
            if length(indices) == k-1
                nullsp = nullspace(d[indices,:])
                @assert size(nullsp,2) == 1
                v = vec(nullsp)
                push!(spanners, canonical!(v))
            end
            if indices[end] != n-1
                push!(indices, indices[end]+1)
            else
                pop!(indices)
                indices[end] += 1
            end
        else
            if indices[end] != n-1
                indices[end] += 1
            else
                pop!(indices)
                indices[end] += 1
            end
        end
    end

    keepers = [1]
    for i in 2:length(spanners)
        alreadyin = false
        for j in keepers
            if isapprox(spanners[i], spanners[j])
                alreadyin = true
                break
            end
        end
        if !alreadyin
            push!(keepers, i)
        end
    end

    return spanners[keepers]
end

function _canonical!(v::Vector{Float64})
    normalize!(v)
    for j in 1:length(v)
        if abs(v[j]) < 1e-8
            v[j] = 0
        end
    end
    sgn = sign(v[findfirst([x!=0 for x in v])])
    for j in 1:length(v)
        if abs(v[j]) < 1e-8
            v[j] = 0
        else
            v[j] *= sgn
        end
    end
    return v
end
