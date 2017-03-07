for instance in ["10104_1_concave_1"]
    folder = joinpath(dirname(@__FILE__),"1D-pwl-instances",instance)
    objs = Dict()

    demand = readdlm(joinpath(folder, "dem.dat"))
    supply = readdlm(joinpath(folder, "sup.dat"))
    numdem = size(demand, 1)
    numsup = size(supply, 1)

    d  = readdlm(joinpath(folder, "mat.dat"))
    fd = readdlm(joinpath(folder, "obj.dat"))
    K = size(d, 2)

    for method in (:Incremental,:MC,:CC,:Logarithmic,:ZigZag,:ZigZagInteger)
        model = Model(solver=solver)
        @variable(model, x[1:numsup,1:numdem] ≥ 0)
        for j in 1:numdem
            # demand constraint
            @constraint(model, sum(x[i,j] for i in 1:numsup) == demand[j])
        end
        for i in 1:numsup
            # supply constraint
            @constraint(model, sum(x[i,j] for j in 1:numdem) == supply[i])
        end

        idx = 1
        obj = AffExpr()
        for i in 1:numsup, j in 1:numdem
            z = piecewiselinear(model, x[i,j], d[idx,:], fd[idx,:], method=method)
            obj += z
        end
        @objective(model, Min, obj)

        stat = solve(model)
        objs[method] = getobjectivevalue(model)
    end
    vals = collect(values(objs))
    for i in 2:length(vals)
        @test isapprox(vals[i-1], vals[i], rtol=1e-4)
    end
end

# for numpieces in [4,8,16,32], variety in 1:5, objective in 1:20
for numpieces in [4], variety in 1:5, objective in 1:20
    folder = joinpath(dirname(@__FILE__),"2D-pwl-instances",string("55",numpieces,"_",variety,"_1"))
    objs = Dict()

    demand = readdlm(joinpath(folder, "dem.dat"))
    supply = readdlm(joinpath(folder, "sup.dat"))
    numdem = size(demand, 1)
    numsup = size(supply, 1)

    d = readdlm(joinpath(folder, "mat.dat"))
    fd = readdlm(joinpath(dirname(@__FILE__),"2D-pwl-instances","objectives",string("55",numpieces,"_",variety,"_",objective)))
    K = size(d, 2)

    for method in (:MC,:CC,:Logarithmic,:ZigZag,:ZigZagInteger)
        model = Model(solver=solver)
        @variable(model, x[1:numsup,1:numdem] ≥ 0)
        @variable(model, y[1:numsup,1:numdem] ≥ 0)

        for j in 1:numdem
            # demand constraint
            @constraint(model, sum(x[i,j] for i in 1:numsup) == demand[j])
            @constraint(model, sum(y[i,j] for i in 1:numsup) == demand[j])
        end
        for i in 1:numsup
            # supply constraint
            @constraint(model, sum(x[i,j] for j in 1:numdem) == supply[i])
            @constraint(model, sum(y[i,j] for j in 1:numdem) == supply[i])
        end

        for i in 1:numdem, j in 1:numsup
            @constraint(model, x[i,j] + y[i,j] ≤ 1.5d[end])
        end

        idx = 1
        obj = AffExpr()
        for i in 1:numsup, j in 1:numdem
            dˣ =  d[idx,:]
            fˣ = reshape(fd[idx,:], length(dˣ), length(dˣ))
            ˣtoⁱ = Dict(dˣ[i] => i for i in 1:length(dˣ))
            fp = (pˣ,pʸ) -> fˣ[ˣtoⁱ[pˣ],ˣtoⁱ[pʸ]]
            z = piecewiselinear(model, x[i,j], y[i,j], BivariatePWLFunction(dˣ, dˣ, fp, pattern=:UnionJack), method=method)
            obj += z
        end
        @objective(model, Min, obj)

        stat = solve(model)
        objs[method] = getobjectivevalue(model)
    end
    vals = collect(values(objs))
    for i in 2:length(vals)
        @test isapprox(vals[i-1], vals[i], rtol=1e-4)
    end
end
