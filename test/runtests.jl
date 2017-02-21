using PiecewiseLinearOpt
using Base.Test

using JuMP, Cbc

const solver = CbcSolver()

let d = linspace(1,2π,8), f = sin
    mCC = Model()
    @variable(mCC, xCC)
    zCC = piecewiselinear(mCC, xCC, d, sin, method=:CC)
    @objective(mCC, Max, zCC)
    @test solve(mCC) == :Optimal
    @test isapprox(getvalue(xCC), 1.75474, rtol=1e-4)
    @test isapprox(getvalue(zCC), 0.98313, rtol=1e-4)

    @constraint(mCC, xCC ≤ 1.5zCC)
    @test solve(mCC) == :Optimal
    @test isapprox(getvalue(xCC), 1.36495, rtol=1e-4)
    @test isapprox(getvalue(zCC), 0.90997, rtol=1e-4)
end

let d = linspace(1,2π,8), f = sin
    mMC = Model()
    @variable(mMC, xMC)
    zMC = piecewiselinear(mMC, xMC, d, sin, method=:MC)
    @objective(mMC, Max, zMC)
    @test solve(mMC) == :Optimal
    @test isapprox(getvalue(xMC), 1.75474, rtol=1e-4)
    @test isapprox(getvalue(zMC), 0.98313, rtol=1e-4)

    @constraint(mMC, xMC ≤ 1.5zMC)
    @test solve(mMC) == :Optimal
    @test isapprox(getvalue(xMC), 1.36495, rtol=1e-4)
    @test isapprox(getvalue(zMC), 0.90997, rtol=1e-4)
end

let d = linspace(1,2π,8), f = sin
    mInc = Model()
    @variable(mInc, xInc)
    zInc = piecewiselinear(mInc, xInc, d, sin, method=:Incremental)
    @objective(mInc, Max, zInc)
    @test solve(mInc) == :Optimal
    @test isapprox(getvalue(xInc), 1.75474, rtol=1e-4)
    @test isapprox(getvalue(zInc), 0.98313, rtol=1e-4)

    @constraint(mInc, xInc ≤ 1.5zInc)
    @test solve(mInc) == :Optimal
    @test isapprox(getvalue(xInc), 1.36495, rtol=1e-4)
    @test isapprox(getvalue(zInc), 0.90997, rtol=1e-4)
end

let d = linspace(1,2π,8), f = sin
    mLog = Model()
    @variable(mLog, xLog)
    zLog = piecewiselinear(mLog, xLog, d, sin, method=:Logarithmic)
    @objective(mLog, Max, zLog)
    @test solve(mLog) == :Optimal
    @test isapprox(getvalue(xLog), 1.75474, rtol=1e-4)
    @test isapprox(getvalue(zLog), 0.98313, rtol=1e-4)

    @constraint(mLog, xLog ≤ 1.5zLog)
    @test solve(mLog) == :Optimal
    @test isapprox(getvalue(xLog), 1.36495, rtol=1e-4)
    @test isapprox(getvalue(zLog), 0.90997, rtol=1e-4)
end

for instance in ["10104_1_concave_1"]
    folder = joinpath(Pkg.dir("PiecewiseLinearOpt"),"test","1D-pwl-instances",instance)
    objs = Dict()

    demand = readdlm(joinpath(folder, "dem.dat"))
    supply = readdlm(joinpath(folder, "sup.dat"))
    numdem = size(demand, 1)
    numsup = size(supply, 1)

    d  = readdlm(joinpath(folder, "mat.dat"))
    fd = readdlm(joinpath(folder, "obj.dat"))
    K = size(d, 2)

    for method in (:Incremental,:MC,:CC,:Logarithmic)
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
