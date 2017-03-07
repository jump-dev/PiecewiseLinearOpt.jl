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
    elseif k == 2
        codes = [[0,0],[0,1],[1,1],[1,0]]
    elseif k < 0
        error()
    else
        codes′ = reflected_gray_codes(k-1)
        codes = vcat([vcat(code,0) for code in codes′],
                     [vcat(code,1) for code in reverse(codes′)])
        return codes
    end
    codes
end

function zigzag_codes(k::Int)
    if k == 0
        codes = Vector{Int}[]
    elseif k == 1
        codes = [[0],[1]]
    elseif k == 2
        codes = [[0,0],[1,0],[0,1],[1,1]]
    elseif k < 0
        error()
    else
        codes′ = zigzag_codes(k-1)
        codes = vcat([vcat(code,0) for code in codes′],
                     [vcat(code,1) for code in codes′])
    end
    codes
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

function integer_zigzag_codes(k::Int)
    if k == 0
        codes = Vector{Int}[]
    elseif k == 1
        codes = [[0],[1]]
    elseif k == 2
        codes = [[0,0],[1,0],[1,1],[2,1]]
    elseif k < 0
        error()
    else
        codes′ = integer_zigzag_codes(k-1)
        offset = [2^(j-2) for j in k:-1:2]
        codes = vcat([vcat(code,        0) for code in codes′],
                     [vcat(code.+offset,1) for code in codes′])
    end
    codes
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
        #foobar
    else # formulations with SOS2 along each dimension
        if method == :Logarithmic
            for tx in 1:nˣ
                sos2_logarthmic_formulation!(m, [λ[tx,ty] for ty in 1:nʸ])
            end
            for ty in 1:nʸ
                sos2_logarthmic_formulation!(m, [λ[tx,ty] for tx in 1:nˣ])
            end
        elseif method == :ZigZag
            for tx in 1:nˣ
                sos2_zigzag_formulation!(m, [λ[tx,ty] for ty in 1:nʸ])
            end
            for ty in 1:nʸ
                sos2_zigzag_formulation!(m, [λ[tx,ty] for tx in 1:nˣ])
            end
        elseif method == :ZigZagInteger
            for tx in 1:nˣ
                sos2_zigzag_general_integer_formulation!(m, [λ[tx,ty] for ty in 1:nʸ])
            end
            for ty in 1:nʸ
                sos2_zigzag_general_integer_formulation!(m, [λ[tx,ty] for tx in 1:nˣ])
            end
        elseif method == :SOS2
            for tx in 1:nˣ
                JuMP.addSOS2(m, [λ[tx,ty] for ty in 1:nʸ])
            end
            for ty in 1:nʸ
                JuMP.addSOS2(m, [λ[tx,ty] for tx in 1:nˣ])
            end
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

            w = JuMP.@variable(m, category=:Bin)
            # If numT==1, then bottom-left is contained in only one point, and so needs separating; otherwise numT==2, and need to offset by one
            JuMP.@constraints(m, begin
                sum(λ[tx,ty] for tx in 1:2:nˣ, ty in    numT :2:nʸ) ≤     w
                sum(λ[tx,ty] for tx in 2:2:nˣ, ty in (3-numT):2:nʸ) ≤ 1 - w
            end)
        else
            w = JuMP.@variable(m, [0:2,0:2], Bin)
            for oˣ in 0:2, oʸ in 0:2
                innoT = fill(true, nˣ, nʸ)
                for (i,j,k) in pwl.T
                    xⁱ, xʲ, xᵏ = pwl.x[i], pwl.x[j], pwl.x[k]
                    iiˣ, iiʸ = ˣtoⁱ[xⁱ[1]], ʸtoʲ[xⁱ[2]]
                    jjˣ, jjʸ = ˣtoⁱ[xʲ[1]], ʸtoʲ[xʲ[2]]
                    kkˣ, kkʸ = ˣtoⁱ[xᵏ[1]], ʸtoʲ[xᵏ[2]]
                    # check to see if one of the points in the triangle falls on the grid
                    if (mod(iiˣ,3) == oˣ && mod(iiʸ,3) == oʸ) || (mod(jjˣ,3) == oˣ && mod(jjʸ,3) == oʸ) || (mod(kkˣ,3) == oˣ && mod(kkʸ,3) == oʸ)
                        innoT[iiˣ,iiʸ] = false
                        innoT[jjˣ,jjʸ] = false
                        innoT[kkˣ,kkʸ] = false
                    end
                end
                JuMP.@constraints(m, begin
                    sum(λ[i,j] for i in 1+oˣ:3:nˣ, j in 1+oʸ:3:nʸ) ≤  1 - w[oˣ,oʸ]
                    sum(λ[i,j] for i in 1:nˣ, j in 1:nʸ if !innoT[i,j]) ≤ w[oˣ,oʸ]
                end)
            end
        end
    end
    z
end
