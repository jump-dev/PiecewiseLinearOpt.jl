# TODO: choose method based on problem size
defaultmethod(x) = :Logarithmic

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

function piecewiselinear(m::JuMP.Model, x::JuMP.Variable, d, f::Function; method=defaultmethod(d))
    initPWL!(m)
    fd = [f(xx) for xx in d]
    piecewiselinear(m, x, d, fd; method=method)
end

piecewiselinear(m::JuMP.Model, x::JuMP.Variable, d, fd; method=defaultmethod(d)) =
    piecewiselinear(m, x, UnivariatePWLFunction(d, fd); method=method)

function piecewiselinear(m::JuMP.Model, x::JuMP.Variable, pwl::UnivariatePWLFunction; method=defaultmethod(pwl))
    initPWL!(m)
    counter = m.ext[:PWL].counter + 1
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
    else # V-formulation method
        λ = JuMP.@variable(m, [1:n], lowerbound=0, upperbound=1, basename="λ_$counter")
        JuMP.@constraint(m, sum(λ) == 1)
        JuMP.@constraint(m, sum(λ[i]* d[i] for i in 1:n) == x)
        JuMP.@constraint(m, sum(λ[i]*fd[i] for i in 1:n) == z)
        if method == :Logarithmic
            sos2_logarthmic_formulation!(m, λ)
        elseif method == :CC
            sos2_cc_formulation!(m, λ)
        elseif method == :MC
            sos2_mc_formulation!(m, λ)
        end
    end
    m.ext[:PWL].counter = counter + 1
    z
end

function piecewiselinear(m::JuMP.Model, x::JuMP.Variable, y::JuMP.Variable, pwl::BivariatePWLFunction; method=defaultmethod(pwl))
    d = pwl.x
    fd = pwl.z
    n = length(d)
    JuMP.@variable(m, minimum(fd) ≤ z ≤ maximum(fd))
    error("Bivariate PWL functions currently not implemented")
end

function sos2_cc_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    y = JuMP.@variable(m, [1:n], Bin, basename="y_$counter")
    JuMP.@constraint(m, sum(y) == 1)
    for i in 1:n-1
        JuMP.@constraint(m, λ[i] ≤ y[i] + y[i+1])
    end
    nothing
end

function sos2_mc_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    γ = JuMP.@variable(m, [1:n,1:n+1], basename="γ_$counter")
    y = JuMP.@variable(m, [1:n], Bin,  basename="y_$counter")
    JuMP.@constraint(m, sum(γ[i,:] for i in 1:n) .== λ)
    for i in 1:n
        JuMP.@constraint(m, γ[i,i] + γ[i,i+1] ≥ y[i])
    end
    nothing
end

function reflected_gray(k::Int)
    if k == 0
        codes = Vector{Int}[]
    elseif k == 1
        codes = [[0],[1]]
    elseif k == 2
        codes = [[0,0],[0,1],[1,1],[1,0]]
    elseif k < 0
        error()
    else
        codes′ = reflected_gray(k-1)
        codes = vcat([vcat(code,0) for code in codes′],
                     [vcat(code,1) for code in reverse(codes′)])
        return codes
    end
    codes
end

function sos2_logarthmic_formulation!(m::JuMP.Model, λ)
    counter = m.ext[:PWL].counter
    n = length(λ)-1
    k = ceil(Int,log2(n))
    H = reflected_gray(k)
    y = JuMP.@variable(m, [1:k], Bin, basename="y_$counter")
    for i in 1:k
        JuMP.@constraints(m, begin
            H[1][i]*λ[1] + sum(min(H[v][i],H[v-1][i])*λ[v] for v in 2:n) + H[n][i]*λ[n+1] ≤ y[i]
            H[1][i]*λ[1] + sum(max(H[v][i],H[v-1][i])*λ[v] for v in 2:n) + H[n][i]*λ[n+1] ≥ y[i]
        end)
    end
    nothing
end
