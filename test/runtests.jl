using Cbc, Gurobi
using Test
using LinearAlgebra

import JuMP
import MathOptInterface
const MOI = MathOptInterface

using PiecewiseLinearOpt
const PLO = PiecewiseLinearOpt

# TODO: Add :SOS2 to this list, once Cbc supports it
methods_1D = (Incremental(), Logarithmic())
@testset "Simple univariate" begin
    m = JuMP.Model(Gurobi.Optimizer)
    JuMP.@variable(m, x)

    s1 = PLO.SegmentPointRep{1,1}([(1.,), (2.,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.,), (3.,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction{1,1,PLO.SegmentPointRep{1,1}}([s1, s2])

    piecewiselinear(m, x, pwl, method=PLO.Incremental())
end

@testset "1D: $method" for method in methods_1D
    model = JuMP.Model(Gurobi.Optimizer)
    JuMP.@variable(model, x)
    segments = PLO.SegmentPointRep{1,1}[]
    d = 7
    xs = range(1,stop=2π, length=(d + 1))
    for i in 1:d
        x_l = xs[i]
        x_r = xs[i+1]
        push!(segments, PLO.SegmentPointRep{1,1}([(x_l,), (x_r,)], [(sin(x_l),), (sin(x_r),)]))
    end
    pwl = PLO.PWLFunction{1,1,PLO.SegmentPointRep{1,1}}(segments)
    z = piecewiselinear(model, x, pwl, method=method)
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

# println("\nbivariate tests")
# @testset "2D: $method, $pattern" for method in methods_2D, pattern in patterns_2D
#     model = JuMP.Model(JuMP.with_optimizer(Cbc.Optimizer))
#     JuMP.@variable(model, x[1:2])
#     d = range(0,stop=1,length=8)
#     f = (x1,x2) -> 2*(x1-1/3)^2 + 3*(x2-4/7)^4
#     z = piecewiselinear(model, x[1], x[2], BivariatePWLFunction(d, d, f, pattern=pattern), method=method)
#     JuMP.@objective(model, Min, z)
#
#     JuMP.optimize!(model)
#
#     @test JuMP.termination_status(model) == MOI.OPTIMAL
#     @test JuMP.value(x[1]) ≈ 0.285714 rtol=1e-4
#     @test JuMP.value(x[2]) ≈ 0.571429 rtol=1e-4
#     @test JuMP.value(z) ≈ 0.004535 rtol=1e-3
#     @test JuMP.objective_value(model) ≈ 0.004535 rtol=1e-3
#     @test JuMP.objective_value(model) ≈ JuMP.value(z) rtol=1e-3
#
#     JuMP.@constraint(model, x[1] ≥ 0.6)
#
#     JuMP.optimize!(model)
#
#     @test JuMP.termination_status(model) == MOI.OPTIMAL
#     @test JuMP.value(x[1]) ≈ 0.6 rtol=1e-4
#     @test JuMP.value(x[2]) ≈ 0.571428 rtol=1e-4
#     @test JuMP.value(z) ≈ 0.148753 rtol=1e-4
#     @test JuMP.objective_value(model) ≈ 0.148753 rtol=1e-3
#     @test JuMP.objective_value(model) ≈ JuMP.value(z) rtol=1e-3
# end
