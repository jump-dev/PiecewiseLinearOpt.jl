
# using Cbc # slow, not recommended
# solver = CbcSolver(logLevel=0, integerTolerance=1e-9, primalTolerance=1e-9, ratioGap=1e-8)

using Gurobi
solver = GurobiSolver(OutputFlag=0)

# using CPLEX
# solver = CplexSolver(CPX_PARAM_SCRIND=0)


using JuMP
using PiecewiseLinearOpt
using Base.Test
using HDF5

methods_1D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:Incremental,:DisaggLogarithmic)
methods_2D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:DisaggLogarithmic)
# methods_2D = (:OptimalIB,) # very slow on all solvers
patterns_2D = (:Upper,:Lower,:UnionJack,:K1,:Random) # not :BestFit because requires function values at midpoints, :OptimalTriangleSelection and :Stencil not supported currently

# tests on network flow model with piecewise-linear objective
# instance data loaded from .h5 files
instance_data_1D = joinpath(dirname(@__FILE__),"1D-pwl-instances.h5")
instance_data_2D = joinpath(dirname(@__FILE__),"2D-pwl-instances.h5")

println("\nunivariate tests")
@testset "instance $instance" for instance in ["10104_1_concave_1"]
    demand = h5read(instance_data_1D, string(instance,"/demand"))
    supply = h5read(instance_data_1D, string(instance,"/supply"))
    ndem = size(demand, 1)
    nsup = size(supply, 1)

    rawd = h5read(instance_data_1D, string(instance,"/d"))
    d = rawd[1,:]
    rawfd = h5read(instance_data_1D, string(instance,"/fd"))
    fd = rawfd[1,:]

    objval1 = NaN
    @testset "1D: $method" for method in methods_1D
        model = Model(solver=solver)
        @variable(model, x[1:nsup,1:ndem] ≥ 0)
        @constraint(model, [j in 1:ndem], sum(x[i,j] for i in 1:nsup) == demand[j])
        @constraint(model, [i in 1:nsup], sum(x[i,j] for j in 1:ndem) == supply[i])
        @objective(model, Min, sum(piecewiselinear(model, x[i,j], d, fd, method=method) for i in 1:nsup, j in 1:ndem))

        @test solve(model) == :Optimal
        if isnan(objval1)
            objval1 = getobjectivevalue(model)
        else
            @test getobjectivevalue(model) ≈ objval1 rtol=1e-4
        end
    end
end

println("\nbivariate tests")
# @testset "numpieces $numpieces, variety $variety, objective $objective" for numpieces in [4,8], variety in 1:5, objective in 1:20
@testset "numpieces $numpieces, variety $variety, objective $objective" for numpieces in [4], variety in 1:1, objective in 1:1
    instance = string(numpieces,"_",variety,"_",objective)

    demand = h5read(instance_data_2D, string(instance,"/demand"))
    supply = h5read(instance_data_2D, string(instance,"/supply"))
    ndem = size(demand, 1)
    nsup = size(supply, 1)

    rawd = h5read(instance_data_2D, string(instance,"/d"))
    d = rawd[1,:]
    rawfd = h5read(instance_data_2D, string(instance,"/fd"))
    fˣ = reshape(rawfd[1,:], length(d), length(d))
    ˣtoⁱ = Dict(d[p] => p for p in 1:length(d))
    f = (pˣ,pʸ) -> fˣ[ˣtoⁱ[pˣ],ˣtoⁱ[pʸ]]

    objval1 = NaN
    @testset "2D: $method, $pattern" for method in methods_2D, pattern in patterns_2D
        model = Model(solver=solver)
        @variable(model, x[1:nsup,1:ndem] ≥ 0)
        @constraint(model, [j in 1:ndem], sum(x[i,j] for i in 1:nsup) == demand[j])
        @constraint(model, [i in 1:nsup], sum(x[i,j] for j in 1:ndem) == supply[i])
        @variable(model, y[1:nsup,1:ndem] ≥ 0)
        @constraint(model, [j in 1:ndem], sum(y[i,j] for i in 1:nsup) == demand[j])
        @constraint(model, [i in 1:nsup], sum(y[i,j] for j in 1:ndem) == supply[i])

        @objective(model, Min, sum(piecewiselinear(model, x[i,j], y[i,j], BivariatePWLFunction(d, d, f, pattern=pattern), method=method, subsolver=solver) for i in 1:nsup, j in 1:ndem))

        @test solve(model) == :Optimal
        if isnan(objval1)
            objval1 = getobjectivevalue(model)
        else
            @test getobjectivevalue(model) ≈ objval1 rtol=1e-4
        end
    end
end
