using PiecewiseLinearOpt
using Base.Test

using JuMP, Cbc
const solver = CbcSolver()

methods_1D = (:CC,:MC,:Logarithmic,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:Incremental)
methods_2D = (:CC,:MC,:Logarithmic,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya)

using HDF5

const instance_data_1D = joinpath(dirname(@__FILE__),"1D-pwl-instances.h5")
const instance_data_2D = joinpath(dirname(@__FILE__),"2D-pwl-instances.h5")


for instance in ["10104_1_concave_1"]
    objs = Dict()

    demand = h5read(instance_data_1D, string(instance,"/demand"))
    supply = h5read(instance_data_1D, string(instance,"/supply"))
    numdem = size(demand, 1)
    numsup = size(supply, 1)

    d  = h5read(instance_data_1D, string(instance,"/d"))
    fd = h5read(instance_data_1D, string(instance,"/fd"))
    K = size(d, 2)

    for method in methods_1D
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
    instance = string(numpieces,"_",variety,"_",objective)
    objs = Dict()

    demand = h5read(instance_data_2D, string(instance,"/demand"))
    supply = h5read(instance_data_2D, string(instance,"/supply"))
    numdem = size(demand, 1)
    numsup = size(supply, 1)

    d  = h5read(instance_data_2D, string(instance,"/d"))
    fd = h5read(instance_data_2D, string(instance,"/fd"))
    K = size(d, 2)

    for method in methods_2D
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
