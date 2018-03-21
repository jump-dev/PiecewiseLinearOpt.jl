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

const VarOrAff = Union{JuMP.Variable,JuMP.AffExpr}

function piecewiselinear(m::JuMP.Model, x::VarOrAff, d, f::Function; method=defaultmethod())
    initPWL!(m)
    fd = [f(xx) for xx in d]
    piecewiselinear(m, x, d, fd; method=method)
end

piecewiselinear(m::JuMP.Model, x::VarOrAff, d, fd; method=defaultmethod()) =
    piecewiselinear(m, x, UnivariatePWLFunction(d, fd); method=method)

function piecewiselinear(m::JuMP.Model, x::VarOrAff, pwl::UnivariatePWLFunction; method=defaultmethod())
    initPWL!(m)
    counter = m.ext[:PWL].counter
    counter += 1
    m.ext[:PWL].counter = counter
    d = [_x[1] for _x in pwl.x]
    fd = pwl.z
    n = length(d)

    if n != length(fd)
        error("You provided a different number of breakpoints ($n) and function values at the breakpoints ($(length(fd)))")
    end

    n == 0 && error("I don't know how to handle a piecewise linear function with no breakpoints")

    z = JuMP.@variable(m, lowerbound=minimum(fd), upperbound=maximum(fd), basename="z_$counter")

    if n == 1
        JuMP.@constraint(m, x ==  d[1])
        JuMP.@constraint(m, z == fd[1])
        return z
    end

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
        ẑ = JuMP.@variable(m, [1:n-1],      basename="ẑ_$counter")
        y = JuMP.@variable(m, [1:n-1], Bin, basename="y_$counter")
        JuMP.@constraint(m, sum(y) == 1)
        JuMP.@constraint(m, sum(x̂) == x)
        JuMP.@constraint(m, sum(ẑ) == z)
        Δ = [(fd[i+1]-fd[i])/(d[i+1]-d[i]) for i in 1:n-1]
        for i in 1:n-1
            JuMP.@constraints(m, begin
                x̂[i] ≥ d[i]  *y[i]
                x̂[i] ≤ d[i+1]*y[i]
                ẑ[i] == fd[i]*y[i] + Δ[i]*(x̂[i]-d[i]*y[i])
            end)
        end
    elseif method in (:DisaggLogarithmic,:DLog)
        γ = JuMP.@variable(m, [i=1:n,j=max(1,i-1):min(n-1,i)], lowerbound=0, upperbound=1)
        JuMP.@constraint(m, sum(γ) == 1)
        JuMP.@constraint(m, γ[1,1]* d[1] + sum((γ[i,i-1]+γ[i,i])* d[i] for i in 2:n-1) + γ[n,n-1]* d[n] == x)
        JuMP.@constraint(m, γ[1,1]*fd[1] + sum((γ[i,i-1]+γ[i,i])*fd[i] for i in 2:n-1) + γ[n,n-1]*fd[n] == z)
        r = ceil(Int, log2(n-1))
        H = reflected_gray_codes(r)
        y = JuMP.@variable(m, [1:r], Bin)
        for j in 1:r
            JuMP.@constraint(m, sum((γ[i,i]+γ[i+1,i])*H[i][j] for i in 1:(n-1)) == y[j])
        end
    else # V-formulation method
        λ = JuMP.@variable(m, [1:n], lowerbound=0, upperbound=1, basename="λ_$counter")
        JuMP.@constraint(m, sum(λ) == 1)
        JuMP.@constraint(m, sum(λ[i]* d[i] for i in 1:n) == x)
        JuMP.@constraint(m, sum(λ[i]*fd[i] for i in 1:n) == z)
        if method in (:Logarithmic,:Log)
            sos2_logarthmic_formulation!(m, λ)
        elseif method in (:LogarithmicIB,:LogIB)
            sos2_logarthmic_IB_formulation!(m, λ)
        elseif method == :CC
            sos2_cc_formulation!(m, λ)
        elseif method in (:ZigZag,:ZZB)
            sos2_zigzag_formulation!(m, λ)
        elseif method in (:ZigZagInteger,:ZZI)
            sos2_zigzag_general_integer_formulation!(m, λ)
        elseif method == :GeneralizedCelaya
            sos2_generalized_celaya_formulation!(m, λ)
        elseif method == :SymmetricCelaya
            sos2_symmetric_celaya_formulation!(m, λ)
        elseif method == :SOS2
            JuMP.addSOS2(m, [i*λ[i] for i in 1:n])
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

function sos2_logarthmic_IB_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))
    y = JuMP.@variable(m, [1:k], Bin, basename="y_$counter")
    _H = reflected_gray_codes(k)
    d = length(_H)
    H = Dict(i => _H[i] for i in 1:d)
    H[0] = H[1]
    H[d+1] = H[d]
    for j in 1:k
        JuMP.@constraints(m, begin
            sum(λ[i] for i in 1:(n+1) if H[i-1][j] == H[i][j] == 1) ≤     y[j]
            sum(λ[i] for i in 1:(n+1) if H[i-1][j] == H[i][j] == 0) ≤ 1 - y[j]
        end)
    end
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

function optimal_IB_scheme!(m::JuMP.Model, λ, pwl, subsolver)
    m.ext[:OptimalIB] = Int[]

    if !haskey(m.ext, :OptimalIBCache)
        m.ext[:OptimalIBCache] = Dict()
    end

    T = pwl.T
    n = maximum(maximum(t) for t in T)
    J = 1:n
    E = fill(false, n, n)
    for t in T
        for i in t, j in t
            E[i,j] = true
        end
    end

    if haskey(m.ext[:OptimalIBCache], pwl.T)
        xx, yy = m.ext[:OptimalIBCache][pwl.T]
        t = size(xx,1)
    else
        t = ceil(Int, log2(length(T)))

        xx, yy = Array(Float64,0,0), Array(Float64,0,0)
        while true
            model = JuMP.Model(solver=subsolver)
            JuMP.@variable(model, x[1:t,1:n], Bin)
            JuMP.@variable(model, y[1:t,1:n], Bin)
            JuMP.@variable(model, z[1:t,1:n,1:n],Bin)

            for j in 1:t
                for r in J, s in J
                    r >= s && continue
                    JuMP.@constraints(model, begin
                        z[j,r,s] <= x[j,r] + x[j,s]
                        z[j,r,s] <= x[j,r] + y[j,r]
                        z[j,r,s] <= x[j,s] + y[j,s]
                        z[j,r,s] <= y[j,r] + y[j,s]
                        z[j,r,s] >= x[j,r] + y[j,s] - 1
                        z[j,r,s] >= x[j,s] + y[j,r] - 1
                    end)
                end
                for r in J
                    JuMP.@constraint(model, x[j,r] + y[j,r] <= 1)
                end
            end

            for r in J, s in J
                r >= s && continue
                if E[r,s]
                    JuMP.@constraint(model, sum(z[j,r,s] for j in 1:t) == 0)
                else
                    JuMP.@constraint(model, sum(z[j,r,s] for j in 1:t) >= 1)
                end
            end

            JuMP.@objective(model, Min, sum(x) + sum(y))
            stat = JuMP.solve(model)
            xx = JuMP.getvalue(x)
            yy = JuMP.getvalue(y)
            if any(isnan, xx) || any(isnan, yy)
                t += 1
            else
                break
            end
        end
        m.ext[:OptimalIBCache][pwl.T] = (xx,yy)
    end

    y = JuMP.@variable(m, [1:t], Bin)

    dˣ = [_x[1] for _x in pwl.x]
    dʸ = [_x[2] for _x in pwl.x]
    uˣ, uʸ = unique(dˣ), unique(dʸ)
    @assert issorted(uˣ)
    @assert issorted(uʸ)
    nˣ, nʸ = length(uˣ), length(uʸ)
    ˣtoⁱ = Dict(uˣ[i] => i for i in 1:nˣ)
    ʸtoʲ = Dict(uʸ[i] => i for i in 1:nʸ)

    for i in 1:t
        JuMP.@constraints(m, begin
            sum(λ[ˣtoⁱ[pwl.x[j][1]],ʸtoʲ[pwl.x[j][2]]] for j in J if xx[i,j] == 1) ≤     y[i]
            sum(λ[ˣtoⁱ[pwl.x[j][1]],ʸtoʲ[pwl.x[j][2]]] for j in J if yy[i,j] == 1) ≤ 1 - y[i]
        end)
    end
    push!(m.ext[:OptimalIB], t)
    nothing
end

piecewiselinear(m::JuMP.Model, x::JuMP.Variable, y::VarOrAff, dˣ, dʸ, f::Function; method=defaultmethod()) =
    piecewiselinear(m, x, y, BivariatePWLFunction(dˣ, dʸ, f); method=method)

function piecewiselinear(m::JuMP.Model, x₁::VarOrAff, x₂::JuMP.Variable, pwl::BivariatePWLFunction; method=defaultmethod(), subsolver=nothing)
    if (method == :OptimalIB) && (subsolver === nothing)
        error("No MIP solver provided to construct optimal IB scheme. Pass a solver object to the piecewiselinear function, e.g. piecewiselinear(m, x₁, x₂, bivariatefunc, method=:OptimalIB, subsolver=GurobiSolver())")
    end

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

    if nˣ == 0 || nʸ == 0
        error("I don't know how to handle a piecewise linear function with zero breakpoints")
    elseif nˣ == 1
        @assert length(dʸ) == nʸ == length(pwl.z)
        return piecewiselinear(m, x₂, UnivariatePWLFunction(dʸ, [pwl.z[i] for i in 1:nʸ]))
    elseif nʸ == 1
        @assert length(dˣ) == nˣ == length(pwl.z)
        return piecewiselinear(m, x₁, UnivariatePWLFunction(dˣ, [pwl.z[i] for i in 1:nˣ]))
    end

    ˣtoⁱ = Dict(uˣ[i] => i for i in 1:nˣ)
    ʸtoʲ = Dict(uʸ[i] => i for i in 1:nʸ)

    fd = Array{Float64}(nˣ, nʸ)
    for (v,fv) in zip(pwl.x, pwl.z)
        # i is the linear index into pwl.x...really want (i,j) pair
        fd[ˣtoⁱ[v[1]],ʸtoʲ[v[2]]] = fv
    end

    z = JuMP.@variable(m, lowerbound=minimum(fd), upperbound=maximum(fd), basename="z_$counter")

    if method == :MC
        x̂₁ = JuMP.@variable(m, [T],      basename="x̂₁_$counter")
        x̂₂ = JuMP.@variable(m, [T],      basename="x̂₂_$counter")
        ẑ  = JuMP.@variable(m, [T],      basename="ẑ_$counter")
        y  = JuMP.@variable(m, [T], Bin, basename="y_$counter")
        JuMP.@constraint(m, sum(y)  == 1)
        JuMP.@constraint(m, sum(x̂₁) == x₁)
        JuMP.@constraint(m, sum(x̂₂) == x₂)
        JuMP.@constraint(m, sum(ẑ)  == z)
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
            JuMP.@constraint(m, ẑ[t] == q[1]*x̂₁[t] + q[2]*x̂₂[t] + q[3]*y[t])
        end
        return z
    elseif method in (:DisaggLogarithmic,:DLog)
        T = pwl.T
        X = pwl.x
        Z = pwl.z
        n = length(X)
        γ = JuMP.@variable(m, [t=T,v=t], lowerbound=0, upperbound=1, basename="γ_$counter")
        Tv = Dict(v => Any[] for v in 1:n)
        for t in T, v in t
            push!(Tv[v], t)
        end
        JuMP.@constraint(m, sum(γ) == 1)
        JuMP.@constraint(m, sum(sum(γ[t,i] for t in Tv[i]) * X[i][1] for i in 1:n) == x₁)
        JuMP.@constraint(m, sum(sum(γ[t,i] for t in Tv[i]) * X[i][2] for i in 1:n) == x₂)
        JuMP.@constraint(m, sum(sum(γ[t,i] for t in Tv[i]) * Z[i]    for i in 1:n) == z)

        r = ceil(Int, log2(length(T)))
        H = reflected_gray_codes(r)
        y = JuMP.@variable(m, [1:r], Bin, basename="y_$counter")
        for j in 1:r
            JuMP.@constraint(m, sum(sum(γ[T[i],v] for v in T[i])*H[i][j] for i in 1:length(T)) == y[j])
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
        return z
    elseif method == :Incremental
        error("Incremental formulation for bivariate functions is not currently implemented.")
        # TODO: Implement algorithm of Geissler et al. (2012)
    elseif method == :OptimalIB
        optimal_IB_scheme!(m, λ, pwl, subsolver)
        return z
    else # formulations with SOS2 along each dimension
        Tx = [sum(λ[tˣ,tʸ] for tˣ in 1:nˣ) for tʸ in 1:nʸ]
        Ty = [sum(λ[tˣ,tʸ] for tʸ in 1:nʸ) for tˣ in 1:nˣ]
        if method in (:Logarithmic,:Log)
            sos2_logarthmic_formulation!(m, Tx)
            sos2_logarthmic_formulation!(m, Ty)
        elseif method in (:LogarithmicIB,:LogIB)
            sos2_logarthmic_IB_formulation!(m, Tx)
            sos2_logarthmic_IB_formulation!(m, Ty)
        elseif method in (:ZigZag,:ZZB)
            sos2_zigzag_formulation!(m, Tx)
            sos2_zigzag_formulation!(m, Ty)
        elseif method in (:ZigZagInteger,:ZZI)
            sos2_zigzag_general_integer_formulation!(m, Tx)
            sos2_zigzag_general_integer_formulation!(m, Ty)
        elseif method == :GeneralizedCelaya
            sos2_generalized_celaya_formulation!(m, Tx)
            sos2_generalized_celaya_formulation!(m, Ty)
        elseif method == :SymmetricCelaya
            sos2_symmetric_celaya_formulation!(m, Tx)
            sos2_symmetric_celaya_formulation!(m, Ty)
        elseif method == :SOS2
            γˣ = JuMP.@variable(m, [1:nˣ], lowerbound=0, upperbound=1, basename="γˣ_$counter")
            γʸ = JuMP.@variable(m, [1:nʸ], lowerbound=0, upperbound=1, basename="γʸ_$counter")
            JuMP.@constraint(m, [tˣ in 1:nˣ], γˣ[tˣ] == sum(λ[tˣ,tʸ] for tʸ in 1:nʸ))
            JuMP.@constraint(m, [tʸ in 1:nʸ], γʸ[tʸ] == sum(λ[tˣ,tʸ] for tˣ in 1:nˣ))
            JuMP.addSOS2(m, [i*γˣ[i] for i in 1:nˣ])
            JuMP.addSOS2(m, [i*γʸ[i] for i in 1:nʸ])
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
        elseif pattern == :OptimalTriangleSelection
            m.ext[:OptimalTriSelect] = Int[]

            if !haskey(m.ext, :OptimalTriSelectCache)
                m.ext[:OptimalTriSelectCache] = Dict()
            end

            J = [(i,j) for i in 1:nˣ, j in 1:nʸ]
            E = Set{Tuple{Tuple{Int,Int},Tuple{Int,Int}}}()
            for t in T
                @assert length(t) == 3
                IJ = [(ˣtoⁱ[pwl.x[i][1]],ʸtoʲ[pwl.x[i][2]]) for i in t]
                im = minimum(ij[1] for ij in IJ)
                iM = maximum(ij[1] for ij in IJ)
                jm = minimum(ij[2] for ij in IJ)
                jM = maximum(ij[2] for ij in IJ)
                @assert im < iM
                @assert im < iM
                if ((im,jm) in IJ) && ((iM,jM) in IJ)
                    push!(E, ((im,jM),(iM,jm)))
                elseif ((im,jM) in IJ) && ((iM,jm) in IJ)
                    push!(E, ((im,jm),(iM,jM)))
                else
                    error()
                end
            end

            if haskey(m.ext[:OptimalTriSelectCache], E)
                xx, yy = m.ext[:OptimalTriSelectCache][E]
                t = JuMP.size(xx,1)
            else
                if subsolver === nothing
                    error("No MIP solver provided to construct optimal triangle selection. Pass a solver object to the piecewiselinear function, e.g. piecewiselinear(m, x₁, x₂, bivariatefunc, method=:Logarithmic, subsolver=GurobiSolver())")
                end
                t = 1
                xx, yy = Array(Float64,0,0), Array(Float64,0,0)
                while true
                    subm = JuMP.Model(solver=subsolver)
                    JuMP.@variable(subm, xˢᵘᵇ[1:t,J], Bin)
                    JuMP.@variable(subm, yˢᵘᵇ[1:t,J], Bin)
                    JuMP.@variable(subm, zˢᵘᵇ[1:t,J,J], Bin)

                    for j in 1:t
                        for r in J, s in J
                            # lexiographic ordering on points on grid
                            (r[1] > s[1] || (r[1] == s[1] && r[2] >= s[2])) && continue
                            JuMP.@constraints(subm, begin
                                zˢᵘᵇ[j,r,s] <= xˢᵘᵇ[j,r] + xˢᵘᵇ[j,s]
                                zˢᵘᵇ[j,r,s] <= xˢᵘᵇ[j,r] + yˢᵘᵇ[j,r]
                                zˢᵘᵇ[j,r,s] <= xˢᵘᵇ[j,s] + yˢᵘᵇ[j,s]
                                zˢᵘᵇ[j,r,s] <= yˢᵘᵇ[j,r] + yˢᵘᵇ[j,s]
                                zˢᵘᵇ[j,r,s] >= xˢᵘᵇ[j,r] + yˢᵘᵇ[j,s] - 1
                                zˢᵘᵇ[j,r,s] >= xˢᵘᵇ[j,s] + yˢᵘᵇ[j,r] - 1
                            end)
                        end
                        for r in J
                            JuMP.@constraint(subm, xˢᵘᵇ[j,r] + yˢᵘᵇ[j,r] <= 1)
                        end
                    end

                    for r in J, s in J
                        # lexiographic ordering on points on grid
                        (r[1] > s[1] || (r[1] == s[1] && r[2] >= s[2])) && continue
                        if (r,s) in E
                            JuMP.@constraint(subm, sum(zˢᵘᵇ[j,r,s] for j in 1:t) >= 1)
                        elseif max(abs(r[1]-s[1]), abs(r[2]-s[2])) == 1
                            JuMP.@constraint(subm, sum(zˢᵘᵇ[j,r,s] for j in 1:t) == 0)
                        end
                    end

                    JuMP.@objective(subm, Min, sum(xˢᵘᵇ) + sum(yˢᵘᵇ))
                    stat = JuMP.solve(subm)
                    if any(isnan, subm.colVal)
                        t += 1
                    else
                        xx = JuMP.getvalue(xˢᵘᵇ)
                        yy = JuMP.getvalue(yˢᵘᵇ)
                        m.ext[:OptimalTriSelectCache][E] = (xx,yy)
                        break
                    end
                end
            end
            y = JuMP.@variable(m, [1:t], Bin, basename="Δselect_$counter")

            for i in 1:t
                JuMP.@constraints(m, begin
                    sum(λ[v[1],v[2]] for v in J if xx[i,v] ≈ 1) ≤     y[i]
                    sum(λ[v[1],v[2]] for v in J if yy[i,v] ≈ 1) ≤ 1 - y[i]
                end)
            end
            push!(m.ext[:OptimalTriSelect], t)

            z
        elseif pattern in (:Stencil,:Stencil9)
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
        else
            # Eⁿᵉ[i,j] = true means that we must cover the edge {(i,j),(i+1,j+1)}
            Eⁿᵉ = fill(false, nˣ-1, nʸ-1)
            for (i,j,k) in pwl.T
                xⁱ, xʲ, xᵏ = pwl.x[i], pwl.x[j], pwl.x[k]
                iiˣ, iiʸ = ˣtoⁱ[xⁱ[1]], ʸtoʲ[xⁱ[2]]
                jjˣ, jjʸ = ˣtoⁱ[xʲ[1]], ʸtoʲ[xʲ[2]]
                kkˣ, kkʸ = ˣtoⁱ[xᵏ[1]], ʸtoʲ[xᵏ[2]]
                IJ = [(iiˣ,iiʸ), (jjˣ,jjʸ), (kkˣ,kkʸ)]
                im = min(iiˣ, jjˣ, kkˣ)
                iM = max(iiˣ, jjˣ, kkˣ)
                jm = min(iiʸ, jjʸ, kkʸ)
                jM = max(iiʸ, jjʸ, kkʸ)
                if ((im,jM) in IJ) && ((iM,jm) in IJ)
                    Eⁿᵉ[im,jm] = true
                else
                    @assert (im,jm) in IJ && (iM,jM) in IJ
                end
            end

            # Diagonal lines running from SW to NE. Grouped with an offset of 3.
            wⁿᵉ = JuMP.@variable(m, [0:2], Bin, basename="wⁿᵉ_$counter")
            for o in 0:2
                Aᵒ = Set{Tuple{Int,Int}}()
                Bᵒ = Set{Tuple{Int,Int}}()
                for offˣ in o:3:(nˣ-2)
                    SWinA = true # Whether we put the SW corner of the next
                                 # triangle to cover in set A
                    for i in (1+offˣ):(nˣ-1)
                        j = i - offˣ
                        (1 ≤ i ≤ nˣ-1) || continue
                        (1 ≤ j ≤ nʸ-1) || continue # should never happen
                        if Eⁿᵉ[i,j] # if we need to cover the edge...
                            if SWinA # figure out which set we need to put it in.
                                     # This depends on previous triangle covered
                                     # in our current line
                                push!(Aᵒ, (i  ,j  ))
                                push!(Bᵒ, (i+1,j+1))
                            else
                                push!(Aᵒ, (i+1,j+1))
                                push!(Bᵒ, (i  ,j  ))
                            end
                            SWinA = !SWinA
                        end
                    end
                end
                for offʸ in (3-o):3:(nʸ-1)
                    SWinA = true
                    for j in (offʸ+1):(nʸ-1)
                        i = j - offʸ
                        (1 ≤ i ≤ nˣ-1) || continue
                        if Eⁿᵉ[i,j]
                            if SWinA
                                push!(Aᵒ, (i  ,j  ))
                                push!(Bᵒ, (i+1,j+1))
                            else
                                push!(Aᵒ, (i+1,j+1))
                                push!(Bᵒ, (i  ,j  ))
                            end
                            SWinA = !SWinA
                        end
                    end
                end
                JuMP.@constraints(m, begin
                    sum(λ[i,j] for (i,j) in Aᵒ) ≤     wⁿᵉ[o]
                    sum(λ[i,j] for (i,j) in Bᵒ) ≤ 1 - wⁿᵉ[o]
                end)
            end

            wˢᵉ = JuMP.@variable(m, [0:2], Bin, basename="wˢᵉ_$counter")
            for o in 0:2
                Aᵒ = Set{Tuple{Int,Int}}()
                Bᵒ = Set{Tuple{Int,Int}}()
                for offˣ in o:3:(nˣ-2)
                    SEinA = true
                    # for i in (1+offˣ):-1:1
                        # j = offˣ - i + 2
                    for j in 1:(nʸ-1)
                        i = nˣ - j - offˣ
                        (1 ≤ i ≤ nˣ-1) || continue
                        if !Eⁿᵉ[i,j]
                            if SEinA
                                push!(Aᵒ, (i+1,j  ))
                                push!(Bᵒ, (i  ,j+1))
                            else
                                push!(Aᵒ, (i  ,j+1))
                                push!(Bᵒ, (i+1,j  ))
                            end
                            SEinA = !SEinA
                        end
                    end
                end
                for offʸ in (3-o):3:(nʸ-1)
                    SEinA = true
                    for j in (offʸ+1):(nʸ-1)
                        i = nˣ - j + offʸ
                        (1 ≤ i ≤ nˣ-1) || continue
                        if !Eⁿᵉ[i,j]
                            if SEinA
                                push!(Aᵒ, (i+1,j  ))
                                push!(Bᵒ, (i  ,j+1))
                            else
                                push!(Aᵒ, (i  ,j+1))
                                push!(Bᵒ, (i+1,j  ))
                            end
                            SEinA = !SEinA
                        end
                    end
                end
                JuMP.@constraints(m, begin
                    sum(λ[i,j] for (i,j) in Aᵒ) ≤     wˢᵉ[o]
                    sum(λ[i,j] for (i,j) in Bᵒ) ≤ 1 - wˢᵉ[o]
                end)
            end
        end
    end
    z
end
