
using Cbc
using Test
using LinearAlgebra
solver = CbcSolver(logLevel=0, integerTolerance=1e-9, primalTolerance=1e-9, ratioGap=1e-8)

# using Gurobi
# solver = GurobiSolver(OutputFlag=0)

# using CPLEX
# solver = CplexSolver(CPX_PARAM_SCRIND=0)


using JuMP
using PiecewiseLinearOpt


methods_1D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:Incremental,:DisaggLogarithmic)
methods_2D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:DisaggLogarithmic)
patterns_2D = (:Upper,:Lower,:BestFit,:UnionJack,:K1,:Random) # :OptimalTriangleSelection and :Stencil not supported currently

println("\nunivariate tests")
@testset "1D: $method" for method in methods_1D
    model = Model(solver=solver)
    @variable(model, x)
    z = piecewiselinear(model, x, linspace(1,2π,8), sin, method=method)
    @objective(model, Max, z)

    @test solve(model) == :Optimal
    @test getvalue(x) ≈ 1.75474 rtol=1e-4
    @test getvalue(z) ≈ 0.98313 rtol=1e-4

    @constraint(model, x ≤ 1.5z)

    @test solve(model) == :Optimal
    @test getvalue(x) ≈ 1.36495 rtol=1e-4
    @test getvalue(z) ≈ 0.90997 rtol=1e-4
end

println("\nbivariate tests")
@testset "2D: $method, $pattern" for method in methods_2D, pattern in patterns_2D
    model = Model(solver=solver)
    @variable(model, x)
    @variable(model, y)
    d = linspace(0,1,8)
    f = (x,y) -> 2*(x-1/3)^2 + 3*(y-4/7)^4
    z = piecewiselinear(model, x, y, BivariatePWLFunction(d, d, f, pattern=pattern), method=method, subsolver=solver)
    @objective(model, Min, z)

    @test solve(model) == :Optimal
    @test getvalue(x) ≈ 0.285714 rtol=1e-4
    @test getvalue(y) ≈ 0.571429 rtol=1e-4
    @test getvalue(z) ≈ 0.004535 rtol=1e-3

    @constraint(model, x ≥ 0.6)

    @test solve(model) == :Optimal
    @test getvalue(x) ≈ 0.6 rtol=1e-4
    @test getvalue(y) ≈ 0.571428 rtol=1e-4
    @test getvalue(z) ≈ 0.148753 rtol=1e-4
end

println("\nbivariate optimal IB scheme tests")
@testset "2D: optimal IB, UnionJack" begin
    model = Model(solver=solver)
    @variable(model, x)
    @variable(model, y)
    d = linspace(0,1,3)
    f = (x,y) -> 2*(x-1/3)^2 + 3*(y-4/7)^4
    z = piecewiselinear(model, x, y, BivariatePWLFunction(d, d, f, pattern=:UnionJack), method=:OptimalIB, subsolver=solver)
    @objective(model, Min, z)

    @test solve(model) == :Optimal
    @test getvalue(x) ≈ 0.5 rtol=1e-4
    @test getvalue(y) ≈ 0.5 rtol=1e-4
    @test getvalue(z) ≈ 0.055634 rtol=1e-3

    @constraint(model, x ≥ 0.6)

    @test solve(model) == :Optimal
    @test getvalue(x) ≈ 0.6 rtol=1e-4
    @test getvalue(y) ≈ 0.5 rtol=1e-4
    @test getvalue(z) ≈ 0.222300 rtol=1e-3
end
