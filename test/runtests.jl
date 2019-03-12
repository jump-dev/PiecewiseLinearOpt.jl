using Cbc
using Gurobi
using Test
using LinearAlgebra

# using Gurobi
# solver = GurobiSolver(OutputFlag=0)

# using CPLEX
# solver = CplexSolver(CPX_PARAM_SCRIND=0)


import JuMP
import MathOptInterface
const MOI = MathOptInterface

using PiecewiseLinearOpt

# TODO: Add :SOS2 to this list, once Cbc supports it
methods_1D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:GeneralizedCelaya,:SymmetricCelaya,:Incremental,:DisaggLogarithmic)
# TODO: Add :SOS2 to this list, once Cbc supports it
# TODO: Add :MC to this list, Cbc (but not Gurobi) gives a different answer below, only for :MC (maybe a bug in Cbc?)
methods_2D = (:CC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:GeneralizedCelaya,:SymmetricCelaya,:DisaggLogarithmic)
patterns_2D = (:Upper,:Lower,:BestFit,:UnionJack,:K1,:Random) # :OptimalTriangleSelection and :Stencil not supported currently

println("\nunivariate tests")
@testset "1D: $method" for method in methods_1D
    model = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer))
    JuMP.@variable(model, x)
    z = piecewiselinear(model, x, range(1,stop=2π, length=8), sin, method=method)
    JuMP.@objective(model, Max, z)

    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.value(x) ≈ 1.75474 rtol=1e-4
    @test JuMP.value(z) ≈ 0.98313 rtol=1e-4

    JuMP.@constraint(model, x ≤ 1.5z)

    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.value(x) ≈ 1.36495 rtol=1e-4
    @test JuMP.value(z) ≈ 0.90997 rtol=1e-4
    @test JuMP.objective_value(model) ≈ 0.90997 rtol=1e-4
    @test JuMP.objective_value(model) ≈ JuMP.value(z) rtol=1e-4
end

println("\nbivariate tests")
@testset "2D: $method, $pattern" for method in methods_2D, pattern in patterns_2D
    model = JuMP.Model(JuMP.with_optimizer(Gurobi.Optimizer))
    JuMP.@variable(model, x[1:2])
    d = range(0,stop=1,length=8)
    f = (x1,x2) -> 2*(x1-1/3)^2 + 3*(x2-4/7)^4
    z = piecewiselinear(model, x[1], x[2], BivariatePWLFunction(d, d, f, pattern=pattern), method=method)
    JuMP.@objective(model, Min, z)

    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.value(x[1]) ≈ 0.285714 rtol=1e-4
    @test JuMP.value(x[2]) ≈ 0.571429 rtol=1e-4
    @test JuMP.value(z) ≈ 0.004535 rtol=1e-3
    @test JuMP.objective_value(model) ≈ 0.004535 rtol=1e-3
    @test JuMP.objective_value(model) ≈ JuMP.value(z) rtol=1e-3

    JuMP.@constraint(model, x[1] ≥ 0.6)

    JuMP.optimize!(model)

    @test JuMP.termination_status(model) == MOI.OPTIMAL
    @test JuMP.value(x[1]) ≈ 0.6 rtol=1e-4
    @test JuMP.value(x[2]) ≈ 0.571428 rtol=1e-4
    @test JuMP.value(z) ≈ 0.148753 rtol=1e-4
    @test JuMP.objective_value(model) ≈ 0.148753 rtol=1e-3
    @test JuMP.objective_value(model) ≈ JuMP.value(z) rtol=1e-3
end
#
# println("\nbivariate optimal IB scheme tests")
# @testset "2D: optimal IB, UnionJack" begin
#     model = JuMP.Model(JuMP.with_optimizer(Cbc.Optimizer))
#     JuMP.@variable(model, x)
#     JuMP.@variable(model, y)
#     d = range(0,stop=1,length=3)
#     f = (x,y) -> 2*(x-1/3)^2 + 3*(y-4/7)^4
#     z = piecewiselinear(model, x, y, BivariatePWLFunction(d, d, f, pattern=:UnionJack), method=:OptimalIB, subsolver=solver)
#     JuMP.@objective(model, Min, z)
#
#     JuMP.optimize!(model)
#
#     @test JuMP.termination_status(model) == MOI.OPTIMAL
#     @test JuMP.value(x) ≈ 0.5 rtol=1e-4
#     @test JuMP.value(y) ≈ 0.5 rtol=1e-4
#     @test JuMP.value(z) ≈ 0.055634 rtol=1e-3
#     @test getobjectivevalue(model) ≈ 0.055634 rtol=1e-3
#     @test getobjectivevalue(model) ≈ JuMP.value(z) rtol=1e-3
#
#     JuMP.@constraint(model, x ≥ 0.6)
#
#     JuMP.optimize!(model)
#
#     @test JuMP.termination_status(model) == MOI.OPTIMAL
#     @test JuMP.value(x) ≈ 0.6 rtol=1e-4
#     @test JuMP.value(y) ≈ 0.5 rtol=1e-4
#     @test JuMP.value(z) ≈ 0.222300 rtol=1e-3
#     @test JuMP.objective_value(model) ≈ 0.222300 rtol=1e-3
#     @test JuMP.objective_value(model) ≈ JuMP.value(z) rtol=1e-3
# end
