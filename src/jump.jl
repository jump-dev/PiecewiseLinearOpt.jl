# TODO: choose method based on problem size
defaultmethod() = :Logarithmic

type PWLData
    counter::Int
    PWLData() = new(0)
end

function initPWL!(m::JuMP.Model)
    if !haskey(m.ext, :PWL)
        m.ext[:PWL] = PWLData()
    end
    nothing
end

function piecewiselinear(m::JuMP.Model, x::JuMP.Variable, d, f::Function; method=defaultmethod())
    initPWL!(m)
    fd = [f(xx) for xx in d]
    piecewiselinear(m, x, d, fd; method=method)
end

piecewiselinear(m::JuMP.Model, x::JuMP.Variable, d, fd; method=defaultmethod()) =
    piecewiselinear(m, x, UnivariatePWLFunction(d, fd); method=method)

function piecewiselinear(m::JuMP.Model, x::JuMP.Variable, pwl::UnivariatePWLFunction; method=defaultmethod())
    initPWL!(m)
    counter = m.ext[:PWL].counter
    counter += 1
    m.ext[:PWL].counter = counter
    d = [_x[1] for _x in pwl.x]
    fd = pwl.z
    n = length(d)
    z = JuMP.@variable(m, lowerbound=minimum(fd), upperbound=maximum(fd), basename="z_$counter")
    if method == :Incremental
        δ = JuMP.@variable(m, [1:n], lowerbound=0, upperbound=1, basename="δ_$counter")
        y = JuMP.@variable(m, [1:n-1], Bin, basename="y_$counter")
        JuMP.@constraint(m, x ==  d[1] + sum(δ[i]*( d[i+1]- d[i]) for i in 1:n-1))
        JuMP.@constraint(m, z == fd[1] + sum(δ[i]*(fd[i+1]-fd[i]) for i in 1:n-1))
        for i in 1:n-1
            JuMP.@constraint(m, δ[i+1] ≤ y[i])
            JuMP.@constraint(m, y[i] ≤ δ[i])
        end
    elseif method == :MC
        x̂ = JuMP.@variable(m, [1:n-1],      basename="x̂_$counter")
        ẑ = JuMP.@variable(m, [1:n-1],      basename="ẑ_$counter")
        y = JuMP.@variable(m, [1:n-1], Bin, basename="y_$counter")
        JuMP.@constraint(m, sum(y) == 1)
        JuMP.@constraint(m, sum(x̂) == x)
        JuMP.@constraint(m, sum(ẑ) == z)
        Δ = [(fd[i+1]-fd[i])/(d[i+1]-d[i]) for i in 1:n-1]
        for i in 1:n-1
            JuMP.@constraints(m, begin
                x̂[i] ≥ d[i]  *y[i]
                x̂[i] ≤ d[i+1]*y[i]
                ẑ[i] == fd[i]*y[i] + Δ[i]*(x̂[i]-d[i]*y[i])
            end)
        end
    else # V-formulation method
        λ = JuMP.@variable(m, [1:n], lowerbound=0, upperbound=1, basename="λ_$counter")
        JuMP.@constraint(m, sum(λ) == 1)
        JuMP.@constraint(m, sum(λ[i]* d[i] for i in 1:n) == x)
        JuMP.@constraint(m, sum(λ[i]*fd[i] for i in 1:n) == z)
        if method == :Logarithmic
            sos2_logarthmic_formulation!(m, λ)
        elseif method == :CC
            sos2_cc_formulation!(m, λ)
        elseif method == :ZigZag
            sos2_zigzag_formulation!(m, λ)
        elseif method == :ZigZagInteger
            sos2_zigzag_general_integer_formulation!(m, λ)
        elseif method == :GeneralizedCelaya
            sos2_generalized_celaya_formulation!(m, λ)
        elseif method == :SymmetricCelaya
            sos2_symmetric_celaya_formulation!(m, λ)
        elseif method == :SOS2
            JuMP.addSOS2(m, [λ[i] for i in 1:n])
        else
            error("Unrecognized method $method")
        end
    end
    z
end

function sos2_cc_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)
    y = JuMP.@variable(m, [1:n-1], Bin, basename="y_$counter")
    JuMP.@constraint(m, sum(y) == 1)
    JuMP.@constraint(m, λ[1] ≤ y[1])
    for i in 2:n-1
        JuMP.@constraint(m, λ[i] ≤ y[i-1] + y[i])
    end
    JuMP.@constraint(m, λ[n] ≤ y[n-1])
    nothing
end

function sos2_mc_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)
    γ = JuMP.@variable(m, [1:n-1, 1:n],     basename="γ_$counter")
    y = JuMP.@variable(m, [1:n-1],     Bin, basename="y_$counter")
    JuMP.@constraint(m, sum(y) == 1)
    JuMP.@constraint(m, sum(γ[i,:] for i in 1:n-1) .== λ)
    for i in 1:n-1
        JuMP.@constraint(m, γ[i,i] + γ[i,i+1] ≥ y[i])
    end
    nothing
end

function sos2_logarthmic_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))
    y = JuMP.@variable(m, [1:k], Bin, basename="y_$counter")

    sos2_encoding_constraints!(m, λ, y, reflected_gray_codes(k), unit_vector_hyperplanes(k))
    nothing
end

function sos2_zigzag_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))
    y = JuMP.@variable(m, [1:k], Bin, basename="y_$counter")

    sos2_encoding_constraints!(m, λ, y, zigzag_codes(k), zigzag_hyperplanes(k))
    nothing
end

function sos2_zigzag_general_integer_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))
    # TODO: tighter upperbounds
    y = JuMP.@variable(m, [i=1:k], Int, lowerbound=0, upperbound=2^(k-i), basename="y_$counter")

    sos2_encoding_constraints!(m, λ, y, integer_zigzag_codes(k), unit_vector_hyperplanes(k))
    nothing
end

function sos2_generalized_celaya_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))

    codes = generalized_celaya_codes(k)
    lb = [minimum(t[i] for t in codes) for i in 1:k]
    ub = [maximum(t[i] for t in codes) for i in 1:k]
    y = JuMP.@variable(m, [i=1:k], Int, lowerbound=lb[i], upperbound=ub[i], basename="y_$counter")

    sos2_encoding_constraints!(m, λ, y, codes, generalized_celaya_hyperplanes(k))
    nothing
end

function sos2_symmetric_celaya_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))

    codes = symmetric_celaya_codes(k)
    lb = [minimum(t[i] for t in codes) for i in 1:k]
    ub = [maximum(t[i] for t in codes) for i in 1:k]
    y = JuMP.@variable(m, [i=1:k], Int, lowerbound=lb[i], upperbound=ub[i], basename="y_$counter")

    sos2_encoding_constraints!(m, λ, y, codes, symmetric_celaya_hyperplanes(k))
    nothing
end

function sos2_encoding_constraints!(m, λ, y, h, B)
    n = length(λ)-1
    for b in B
        JuMP.@constraints(m, begin
            dot(b,h[1])*λ[1] + sum(min(dot(b,h[v]),dot(b,h[v-1]))*λ[v] for v in 2:n) + dot(b,h[n])*λ[n+1] ≤ dot(b,y)
            dot(b,h[1])*λ[1] + sum(max(dot(b,h[v]),dot(b,h[v-1]))*λ[v] for v in 2:n) + dot(b,h[n])*λ[n+1] ≥ dot(b,y)
        end)
    end
    nothing
end

function reflected_gray_codes(k::Int)
    if k == 0
        codes = Vector{Int}[]
    elseif k == 1
        codes = [[0],[1]]
    else
        codes′ = reflected_gray_codes(k-1)
        codes = vcat([vcat(code,0) for code in codes′],
                     [vcat(code,1) for code in reverse(codes′)])
        return codes
    end
    codes
end

function zigzag_codes(k::Int)
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        codes = [[0],[1]]
    else
        codes′ = zigzag_codes(k-1)
        codes = vcat([vcat(code,0) for code in codes′],
                     [vcat(code,1) for code in codes′])
    end
    codes
end

function integer_zigzag_codes(k::Int)
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        codes = [[0],[1]]
    else
        codes′ = integer_zigzag_codes(k-1)
        offset = [2^(j-2) for j in k:-1:2]
        codes = vcat([vcat(code,        0) for code in codes′],
                     [vcat(code.+offset,1) for code in codes′])
    end
    codes
end

function generalized_celaya_codes(k::Int)
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0],[1]]
    elseif k == 2
        return [[0,0],[1,1],[1,0],[0,1]]
    else
        codes′ = generalized_celaya_codes(k-1)
        n = length(codes′)
        hp = Int(n/2)
        firstcodes  = [codes′[i] for i in 1:hp]
        secondcodes = [codes′[i] for i in (hp+1):n]
        return vcat([vcat(codes′[i], 0) for i in     1 :  n],
                    [vcat(codes′[i], 1) for i in     1 : hp],
                    [vcat(codes′[i],-1) for i in (hp+1):  n])
    end
end

function symmetric_celaya_codes(k::Int)
    if k <= 0
        error("Invalid code length $k")
    elseif k == 1
        return [[0],[1]]
    else
        codes′ = generalized_celaya_codes(k-1)
        n = length(codes′)
        hp = Int(n/2)
        firstcodes  = [codes′[i] for i in 1:hp]
        secondcodes = [codes′[i] for i in (hp+1):n]
        return vcat([vcat(codes′[i], 0) for i in     1 :  n],
                    [vcat(codes′[i], 1) for i in     1 : hp],
                    [vcat(codes′[i],-1) for i in (hp+1):  n])
    end
end

function zigzag_hyperplanes(k::Int)
    hps = Vector{Int}[]
    for i in 1:k
        hp = zeros(Int, k)
        hp[i] = 1
        for j in (i+1):k
            hp[j] = 2^(j-i-1)
        end
        push!(hps, hp)
    end
    hps
end

function unit_vector_hyperplanes(k::Int)
    hps = Vector{Int}[]
    for i in 1:k
        hp = zeros(Int,k)
        hp[i] = 1
        push!(hps, hp)
    end
    hps
end

function generalized_celaya_hyperplanes(k::Int)
    C = generalized_celaya_codes(k)
    compute_hyperplanes(C)
end

function symmetric_celaya_hyperplanes(k::Int)
    C = symmetric_celaya_codes(k)
    compute_hyperplanes(C)
end

function compute_hyperplanes{T}(C::Vector{Vector{T}})
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
        approx12 || push!(indices,2)
        approx13 || (!approx12 && approx23) || push!(indices,3)
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
        alreadyin || push!(keepers, i)
    end
    spanners[keepers]
end

function canonical!(v::Vector{Float64})
    normalize!(v)
    for j in 1:length(v)
        if abs(v[j]) < 1e-8
            v[j] = 0
        end
    end
    sgn = sign(v[findfirst(v)])
    for j in 1:length(v)
        if abs(v[j]) < 1e-8
            v[j] = 0
        else
            v[j] *= sgn
        end
    end
    v
end

piecewiselinear(m::JuMP.Model, x::JuMP.Variable, y::JuMP.Variable, dˣ, dʸ, f::Function; method=defaultmethod()) =
    piecewiselinear(m, x, y, BivariatePWLFunction(dˣ, dʸ, f); method=method)

function piecewiselinear(m::JuMP.Model, x₁::JuMP.Variable, x₂::JuMP.Variable, pwl::BivariatePWLFunction; method=defaultmethod())
    initPWL!(m)
    counter = m.ext[:PWL].counter
    counter += 1
    m.ext[:PWL].counter = counter
    dˣ = [_x[1] for _x in pwl.x]
    dʸ = [_x[2] for _x in pwl.x]
    uˣ, uʸ = unique(dˣ), unique(dʸ)
    @assert issorted(uˣ)
    @assert issorted(uʸ)

    T = pwl.T

    nˣ, nʸ = length(uˣ), length(uʸ)

    ˣtoⁱ = Dict(uˣ[i] => i for i in 1:nˣ)
    ʸtoʲ = Dict(uʸ[i] => i for i in 1:nʸ)

    fd = Array(Float64, nˣ, nʸ)
    for (v,fv) in zip(pwl.x, pwl.z)
        # i is the linear index into pwl.x...really want (i,j) pair
        fd[ˣtoⁱ[v[1]],ʸtoʲ[v[2]]] = fv
    end

    z = JuMP.@variable(m, lowerbound=minimum(fd), upperbound=maximum(fd), basename="z_$counter")

    if method == :MC
        x̂₁ = JuMP.@variable(m, [T],      basename="x̂₁_$counter")
        x̂₂ = JuMP.@variable(m, [T],      basename="x̂₂_$counter")
        ẑ  = JuMP.@variable(m, [T],      basename="ẑ_$counter")
        y  = JuMP.@variable(m, [T], Bin, basename="y_$counter")
        JuMP.@constraint(m, sum(y)  == 1)
        JuMP.@constraint(m, sum(x̂₁) == x₁)
        JuMP.@constraint(m, sum(x̂₂) == x₂)
        JuMP.@constraint(m, sum(ẑ)  == z)
        for t in T
            @assert length(t) == 3
            r¹,  r²,  r³  = pwl.x[t[1]], pwl.x[t[2]], pwl.x[t[3]]
            fz¹, fz², fz³ = pwl.z[t[1]], pwl.z[t[2]], pwl.z[t[3]]
            for P in [[1,2,3],[2,3,1],[3,1,2]]
                p¹, p², p³ = [r¹, r², r³][P]
                p¹, p², p³ = pwl.x[t[P[1]]], pwl.x[t[P[2]]], pwl.x[t[P[3]]]

                A = [p¹[1] p¹[2] 1
                     p²[1] p²[2] 1
                     p³[1] p³[2] 1]
                b = [0,0,1]
                q = A \ b

                @assert isapprox(q[1]*p¹[1] + q[2]*p¹[2] + q[3], 0, atol=1e-4)
                @assert isapprox(q[1]*p²[1] + q[2]*p²[2] + q[3], 0, atol=1e-4)
                @assert isapprox(q[1]*p³[1] + q[2]*p³[2] + q[3], 1, atol=1e-4)

                JuMP.@constraint(m, q[1]*x̂₁[t] + q[2]*x̂₂[t] + q[3]*y[t] ≥ 0)
            end
            A = [r¹[1] r¹[2] 1
                 r²[1] r²[2] 1
                 r³[1] r³[2] 1]
            b = [fz¹, fz², fz³]
            q = A \ b

            @assert isapprox(q[1]*r¹[1] + q[2]*r¹[2] + q[3], fz¹, atol=1e-4)
            @assert isapprox(q[1]*r²[1] + q[2]*r²[2] + q[3], fz², atol=1e-4)
            @assert isapprox(q[1]*r³[1] + q[2]*r³[2] + q[3], fz³, atol=1e-4)
            JuMP.@constraint(m, ẑ[t] == q[1]*x̂₁[t] + q[2]*x̂₂[t] + q[3]*y[t])
        end
        return z
    end

    λ = JuMP.@variable(m, [1:nˣ,1:nʸ], lowerbound=0, upperbound=1, basename="λ_$counter")
    JuMP.@constraint(m, sum(λ) == 1)
    JuMP.@constraint(m, sum(λ[i,j]*uˣ[i]   for i in 1:nˣ, j in 1:nʸ) == x₁)
    JuMP.@constraint(m, sum(λ[i,j]*uʸ[j]   for i in 1:nˣ, j in 1:nʸ) == x₂)
    JuMP.@constraint(m, sum(λ[i,j]*fd[i,j] for i in 1:nˣ, j in 1:nʸ) == z)

    if method == :CC
        T = pwl.T
        y = JuMP.@variable(m, [T], Bin, basename="y_$counter")
        JuMP.@constraint(m, sum(y) == 1)

        Ts = Dict((i,j) => Vector{Int}[] for i in 1:nˣ, j in 1:nʸ)
        for t in T, ind in t
            rˣ, rʸ = pwl.x[ind]
            i, j = ˣtoⁱ[rˣ], ʸtoʲ[rʸ]
            push!(Ts[(i,j)], t)
        end
        for i in 1:nˣ, j in 1:nʸ
            JuMP.@constraint(m, λ[i,j] ≤ sum(y[t] for t in Ts[(i,j)]))
        end
    elseif method == :Incremental
        error()
    else # formulations with SOS2 along each dimension
        if method == :Logarithmic
            sos2_logarthmic_formulation!(m, [sum(λ[tˣ,tʸ] for tˣ in 1:nˣ) for tʸ in 1:nʸ])
            sos2_logarthmic_formulation!(m, [sum(λ[tˣ,tʸ] for tʸ in 1:nʸ) for tˣ in 1:nˣ])
        elseif method == :ZigZag
            sos2_zigzag_formulation!(m, [sum(λ[tˣ,tʸ] for tˣ in 1:nˣ) for tʸ in 1:nʸ])
            sos2_zigzag_formulation!(m, [sum(λ[tˣ,tʸ] for tʸ in 1:nʸ) for tˣ in 1:nˣ])
        elseif method == :ZigZagInteger
            sos2_zigzag_general_integer_formulation!(m, [sum(λ[tˣ,tʸ] for tˣ in 1:nˣ) for tʸ in 1:nʸ])
            sos2_zigzag_general_integer_formulation!(m, [sum(λ[tˣ,tʸ] for tʸ in 1:nʸ) for tˣ in 1:nˣ])
        elseif method == :GeneralizedCelaya
            sos2_generalized_celaya_formulation!(m, [sum(λ[tˣ,tʸ] for tˣ in 1:nˣ) for tʸ in 1:nʸ])
            sos2_generalized_celaya_formulation!(m, [sum(λ[tˣ,tʸ] for tʸ in 1:nʸ) for tˣ in 1:nˣ])
        elseif method == :SymmetricCelaya
            sos2_symmetric_celaya_formulation!(m, [sum(λ[tˣ,tʸ] for tˣ in 1:nˣ) for tʸ in 1:nʸ])
            sos2_symmetric_celaya_formulation!(m, [sum(λ[tˣ,tʸ] for tʸ in 1:nʸ) for tˣ in 1:nˣ])
        elseif method == :SOS2
            γˣ = JuMP.@variable(m, [1:nˣ], lowerbound=0, upperbound=1, basename="γˣ_$counter")
            γʸ = JuMP.@variable(m, [1:nʸ], lowerbound=0, upperbound=1, basename="γʸ_$counter")
            JuMP.@constraint(m, [tˣ in 1:nˣ], γˣ[tˣ] == sum(λ[tˣ,tʸ] for tʸ in 1:nʸ))
            JuMP.@constraint(m, [tʸ in 1:nʸ], γʸ[tʸ] == sum(λ[tˣ,tʸ] for tˣ in 1:nˣ))
            JuMP.addSOS2(m, γˣ)
            JuMP.addSOS2(m, γʸ)
        else
            error("Unrecognized method $method")
        end

        pattern = pwl.meta[:structure]
        if pattern == :UnionJack
            numT = 0
            minˣ, minʸ = uˣ[1], uʸ[1]
            # find the index of the bottom-left point
            idx = findfirst(pwl.x) do w; w == (minˣ, minʸ); end
            for t in pwl.T
                if idx in t
                    numT += 1
                end
            end
            @assert 1 <= numT <= 2

            w = JuMP.@variable(m, category=:Bin, basename="w_$counter")
            # If numT==1, then bottom-left is contained in only one point, and so needs separating; otherwise numT==2, and need to offset by one
            JuMP.@constraints(m, begin
                sum(λ[tx,ty] for tx in 1:2:nˣ, ty in    numT :2:nʸ) ≤     w
                sum(λ[tx,ty] for tx in 2:2:nˣ, ty in (3-numT):2:nʸ) ≤ 1 - w
            end)
        elseif pattern == :K1
            # Assumption is that the triangulation is uniform, as defined in types.jl
            w = JuMP.@variable(m, [1:2], Bin, basename="w_$counter")
            JuMP.@constraints(m, begin
                sum(λ[tx,ty] for tx in 1:nˣ, ty in 1:nʸ if mod(tx,2) == mod(ty,2) && (tx+ty) in 2:4:(nˣ+nʸ)) ≤ w[1]
                sum(λ[tx,ty] for tx in 1:nˣ, ty in 1:nʸ if mod(tx,2) == mod(ty,2) && (tx+ty) in 4:4:(nˣ+nʸ)) ≤ 1 - w[1]
                sum(λ[tx,ty] for tx in 1:nˣ, ty in 1:nʸ if mod(tx,2) != mod(ty,2) && (tx+ty) in 3:4:(nˣ+nʸ)) ≤ w[2]
                sum(λ[tx,ty] for tx in 1:nˣ, ty in 1:nʸ if mod(tx,2) != mod(ty,2) && (tx+ty) in 5:4:(nˣ+nʸ)) ≤ 1 - w[2]
            end)
        else
            w = JuMP.@variable(m, [1:3,1:3], Bin, basename="w_$counter")
            for oˣ in 1:3, oʸ in 1:3
                innoT = fill(true, nˣ, nʸ)
                for (i,j,k) in pwl.T
                    xⁱ, xʲ, xᵏ = pwl.x[i], pwl.x[j], pwl.x[k]
                    iiˣ, iiʸ = ˣtoⁱ[xⁱ[1]], ʸtoʲ[xⁱ[2]]
                    jjˣ, jjʸ = ˣtoⁱ[xʲ[1]], ʸtoʲ[xʲ[2]]
                    kkˣ, kkʸ = ˣtoⁱ[xᵏ[1]], ʸtoʲ[xᵏ[2]]
                    # check to see if one of the points in the triangle falls on the grid
                    if (mod1(iiˣ,3) == oˣ && mod1(iiʸ,3) == oʸ) || (mod1(jjˣ,3) == oˣ && mod1(jjʸ,3) == oʸ) || (mod1(kkˣ,3) == oˣ && mod1(kkʸ,3) == oʸ)
                        innoT[iiˣ,iiʸ] = false
                        innoT[jjˣ,jjʸ] = false
                        innoT[kkˣ,kkʸ] = false
                    end
                end
                JuMP.@constraints(m, begin
                    sum(λ[i,j] for i in oˣ:3:nˣ, j in oʸ:3:nʸ) ≤  1 - w[oˣ,oʸ]
                    sum(λ[i,j] for i in 1:nˣ, j in 1:nʸ if innoT[i,j]) ≤ w[oˣ,oʸ]
                end)
            end
        end
    end
    z
end
