using PiecewiseLinearOpt
using Base.Test

using JuMP, CPLEX

const solver = CplexSolver(CPX_PARAM_TILIM=30*60.0, CPX_PARAM_MIPCBREDLP=0)

fp = open("1D-pwl-results.csv", "w+")
fp2 = open("1D-pwl-objective-value.csv", "w+")

methods = (:MomentCurve,:Incremental,:MC,:CC,:Logarithmic)

println(fp,  "instance, ", join(methods, ", "))
println(fp2, "instance, ", join(methods, ", "))

for instance in readdir(joinpath(Pkg.dir("PiecewiseLinearOpt"),"test","1D-pwl-instances"))
    print(fp, "$instance")
    print(fp2, "$instance")

    folder = joinpath(Pkg.dir("PiecewiseLinearOpt"),"test","1D-pwl-instances",instance)

    demand = readdlm(joinpath(folder, "dem.dat"))
    supply = readdlm(joinpath(folder, "sup.dat"))
    numdem = size(demand, 1)
    numsup = size(supply, 1)

    d  = readdlm(joinpath(folder, "mat.dat"))
    fd = readdlm(joinpath(folder, "obj.dat"))
    K = size(d, 2)

    for method in methods
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

        tm = @elapsed solve(model)
	print(fp, ", $tm")
	flush(fp)
	print(fp2, ", $(getobjectivevalue(model))")
	flush(fp2)

	error()
    end
    println(fp)
    println(fp2)
end
