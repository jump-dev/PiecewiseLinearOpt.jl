# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using PiecewiseLinearOpt
using HiGHS
using JuMP
using LinearAlgebra
using Test

const PLO = PiecewiseLinearOpt

optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

const methods_1D = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    Incremental(),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    NativeSOS2(),
    ZigZagBinary(),
    ZigZagInteger()
]

@testset "Simple univariate" for method in methods_1D
    model = Model(optimizer)
    @variable(model, x)

    s1 = PLO.SegmentPointRep{1,1}([(1.,), (2.,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.,), (3.,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())

    y = piecewiselinear(model, (x,), pwl, method=method)
    @objective(model, Min, y[1])

    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 3.0 rtol=1e-4
    @test value(y[1]) ≈ 1.0 rtol=1e-4
end
const sos2_methods = [
    ConvexCombination(),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    NativeSOS2(),
    ZigZagBinary(),
    ZigZagInteger()
]
const methods_2D_gen = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    #OptimalIndependentBranching(optimizer),
    [NineStencil(sos2_method) for sos2_method in methods_1D]...,
    [OptimalTriangleSelection(optimizer, sos2_method) for sos2_method in methods_1D]...,
    [SixStencil(sos2_method) for sos2_method in methods_1D]...,
]
@testset "Simple bivariate" for method in methods_2D_gen
    model = Model(optimizer)
    @variable(model, x[1:2])
    s1 = PLO.SegmentPointRep{2,1}([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)], [(0.0,), (1.0,), (2.0,)])
    s2 = PLO.SegmentPointRep{2,1}([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)], [(0.0,), (3.0,), (2.0,)])
    pwl = PLO.PWLFunction{2,1,PLO.SegmentPointRep{2,1}}([s1, s2], PLO.UnstructuredTriangulation())

    y = piecewiselinear(model, (x[1], x[2]), pwl, method = method)
    @objective(model, Min, y[1])

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.0 rtol=1e-4
    @test value(x[2]) ≈ 0.0 rtol=1e-4
    @test value(y[1]) ≈ 0.0 rtol=1e-4
end

@testset "1D: $method" for method in methods_1D
    model = Model(optimizer)
    @variable(model, x)
    d = 7
    xs = collect(range(1,stop=2π, length=(d + 1)))
    zs = sin.(xs)
    pwl = PLO.UnivariatePWLFunction(xs, zs)
    y = piecewiselinear(model, x, pwl, method=method)
    @objective(model, Max, y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.75474 rtol=1e-4
    @test value(y) ≈ 0.98313 rtol=1e-4
    @test objective_value(model) ≈ 0.98313 rtol=1e-4
    @test objective_value(model) ≈ value(y) rtol=1e-4
    @constraint(model, x ≤ 1.5y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.36495 rtol=1e-4
    @test value(y) ≈ 0.90997 rtol=1e-4
    @test objective_value(model) ≈ 0.90997 rtol=1e-4
    @test objective_value(model) ≈ value(y) rtol=1e-4
end

patterns = [:Upper, :Lower, :BestFit, :K1, :UnionJack, :Random]
method_pattern = vec(collect(Iterators.product(methods_2D_gen, patterns)))
k1_methods = [(method, :K1) for method in [K1(sos2_method) for sos2_method in methods_1D]]
uj_methods = [(method, :UnionJack) for method in [UnionJack(sos2_method) for sos2_method in methods_1D]]
append!(method_pattern, k1_methods)
append!(method_pattern, uj_methods)

@testset "2D: $method, $pattern" for (method, pattern) in method_pattern
    model = Model(optimizer)
    @variable(model, x[1:2])
    d = range(0,stop=1,length=8)
    f = (x1,x2) -> 2*(x1-1/3)^2 + 3*(x2-4/7)^4
    pwl = PLO.BivariatePWLFunction(d, d, f, pattern=pattern)
    z = piecewiselinear(model, x[1], x[2], pwl, method=method)
    @objective(model, Min, z)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.285714 rtol=1e-4
    @test value(x[2]) ≈ 0.571429 rtol=1e-4
    @test value(z) ≈ 0.004535 rtol=1e-3
    @test objective_value(model) ≈ 0.004535 rtol=1e-3
    @test objective_value(model) ≈ value(z) rtol=1e-3

    @constraint(model, x[1] ≥ 0.6)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.6 rtol=1e-4
    @test value(x[2]) ≈ 0.571428 rtol=1e-4
    @test value(z) ≈ 0.148753 rtol=1e-4
    @test objective_value(model) ≈ 0.148753 rtol=1e-3
    @test objective_value(model) ≈ value(z) rtol=1e-3
end
