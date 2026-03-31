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
    ZigZagInteger(),
]

@testset "Simple univariate" for method in methods_1D
    model = Model(optimizer)
    @variable(model, x)

    s1 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())

    y = piecewiselinear(model, (x,), pwl; method = method)
    @objective(model, Min, y[1])

    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 3.0 rtol = 1e-4
    @test value(y[1]) ≈ 1.0 rtol = 1e-4
end

@testset "Univariate pwlinear" begin
    d = 0:0.01:1
    f = (xi -> xi^2)
    fd = [f(xi) for xi in d]
    pwl = UnivariatePWLFunction(d, f)

    model = Model(optimizer)
    @variable(model, x)
    y1 = piecewiselinear(model, (x,), pwl)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y1[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y1[1]) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(optimizer)
    @variable(model, x)
    y2 = piecewiselinear(model, x, d, f)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y2)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y2) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(optimizer)
    @variable(model, x)
    y3 = piecewiselinear(model, x, d, fd)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y3)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y3) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4
end

@testset "Bivariate pwlinear" begin
    d = 0:0.05:1
    f = (xi, yi) -> xi^2 + yi^2
    pwl = BivariatePWLFunction(d, d, f)

    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    z1 = piecewiselinear(model, (x, y), pwl)
    @constraint(model, x ≤ 0.75)
    @constraint(model, y ≤ 0.75)
    @objective(model, Max, z1[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 0.75 rtol = 1e-4
    @test value(y) ≈ 0.75 rtol = 1e-4
    @test value(z1[1]) ≈ 1.125 rtol = 1e-4

    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    z2 = piecewiselinear(model, x, y, d, d, f)
    @constraint(model, x ≤ 0.75)
    @constraint(model, y ≤ 0.75)
    @objective(model, Max, z2)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 0.75 rtol = 1e-4
    @test value(y) ≈ 0.75 rtol = 1e-4
    @test value(z2) ≈ 1.125 rtol = 1e-4
end

const sos2_methods = [
    ConvexCombination(),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    NativeSOS2(),
    ZigZagBinary(),
    ZigZagInteger(),
]
const methods_2D_gen = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    #OptimalIndependentBranching(optimizer),
    [NineStencil(sos2_method) for sos2_method in methods_1D]...,
    [
        OptimalTriangleSelection(optimizer, sos2_method) for
        sos2_method in methods_1D
    ]...,
    [SixStencil(sos2_method) for sos2_method in methods_1D]...,
]
@testset "Simple bivariate" for method in methods_2D_gen
    model = Model(optimizer)
    @variable(model, x[1:2])
    s1 = PLO.SegmentPointRep{2,1}(
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0,), (1.0,), (2.0,)],
    )
    s2 = PLO.SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(0.0,), (3.0,), (2.0,)],
    )
    pwl = PLO.PWLFunction{2,1,PLO.SegmentPointRep{2,1}}(
        [s1, s2],
        PLO.UnstructuredTriangulation(),
    )

    y = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    @objective(model, Min, y[1])

    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.0 rtol = 1e-4
    @test value(x[2]) ≈ 0.0 rtol = 1e-4
    @test value(y[1]) ≈ 0.0 rtol = 1e-4
end

@testset "1D: $method" for method in methods_1D
    model = Model(optimizer)
    @variable(model, x)
    d = 7
    xs = collect(range(1; stop = 2π, length = (d + 1)))
    zs = sin.(xs)
    y = piecewiselinear(model, x, xs, zs; method = method)
    @objective(model, Max, y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.75474 rtol = 1e-4
    @test value(y) ≈ 0.98313 rtol = 1e-4
    @test objective_value(model) ≈ 0.98313 rtol = 1e-4
    @test objective_value(model) ≈ value(y) rtol = 1e-4
    @constraint(model, x ≤ 1.5y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.36495 rtol = 1e-4
    @test value(y) ≈ 0.90997 rtol = 1e-4
    @test objective_value(model) ≈ 0.90997 rtol = 1e-4
    @test objective_value(model) ≈ value(y) rtol = 1e-4
end

patterns = [:Upper, :Lower, :BestFit, :K1, :UnionJack, :Random]
method_pattern = vec(collect(Iterators.product(methods_2D_gen, patterns)))
k1_methods = [
    (method, :K1) for method in [K1(sos2_method) for sos2_method in methods_1D]
]
uj_methods = [
    (method, :UnionJack) for
    method in [UnionJack(sos2_method) for sos2_method in methods_1D]
]
append!(method_pattern, k1_methods)
append!(method_pattern, uj_methods)

@testset "2D: $method, $pattern" for (method, pattern) in method_pattern
    model = Model(optimizer)
    @variable(model, x[1:2])
    d = range(0; stop = 1, length = 8)
    f = (x1, x2) -> 2 * (x1 - 1 / 3)^2 + 3 * (x2 - 4 / 7)^4
    z = piecewiselinear(
        model,
        x[1],
        x[2],
        d,
        d,
        f;
        method = method,
        pattern = pattern,
    )
    @objective(model, Min, z)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.285714 rtol = 1e-4
    @test value(x[2]) ≈ 0.571429 rtol = 1e-4
    @test value(z) ≈ 0.004535 rtol = 1e-3
    @test objective_value(model) ≈ 0.004535 rtol = 1e-3
    @test objective_value(model) ≈ value(z) rtol = 1e-3

    @constraint(model, x[1] ≥ 0.6)
    optimize!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.6 rtol = 1e-4
    @test value(x[2]) ≈ 0.571428 rtol = 1e-4
    @test value(z) ≈ 0.148753 rtol = 1e-4
    @test objective_value(model) ≈ 0.148753 rtol = 1e-3
    @test objective_value(model) ≈ value(z) rtol = 1e-3
end

# Test OptimalTriangleSelection with one method for coverage
@testset "OptimalTriangleSelection" begin
    method = OptimalTriangleSelection(optimizer, ConvexCombination())
    model = Model(optimizer)
    @variable(model, x[1:2])
    s1 = PLO.SegmentPointRep{2,1}([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)], [(0.0,), (1.0,), (2.0,)])
    s2 = PLO.SegmentPointRep{2,1}([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)], [(0.0,), (3.0,), (2.0,)])
    pwl = PLO.PWLFunction{2,1,PLO.SegmentPointRep{2,1}}([s1, s2], PLO.UnstructuredTriangulation())
    y = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    @objective(model, Min, y[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.0 rtol = 1e-4
    @test value(x[2]) ≈ 0.0 rtol = 1e-4
    @test value(y[1]) ≈ 0.0 rtol = 1e-4
end

# Test OptimalIndependentBranching (small problem only, the MIP can be very slow)
@testset "OptimalIndependentBranching" begin
    method = OptimalIndependentBranching(optimizer)
    model = Model(optimizer)
    @variable(model, x[1:2])
    s1 = PLO.SegmentPointRep{2,1}([(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)], [(0.0,), (1.0,), (2.0,)])
    s2 = PLO.SegmentPointRep{2,1}([(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)], [(0.0,), (3.0,), (2.0,)])
    pwl = PLO.PWLFunction{2,1,PLO.SegmentPointRep{2,1}}([s1, s2], PLO.UnstructuredTriangulation())
    y = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    @objective(model, Min, y[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
end

# Test MultipleChoice (requires SegmentHyperplaneRep)
@testset "MultipleChoice" begin
    method = MultipleChoice()

    # 1D: Two segments: f(x) = x + 1.5 on [1,2], f(x) = -2.5x + 8.5 on [2,3]
    # Segment 1: domain x >= 1 and x <= 2, output y = x + 1.5
    s1 = PLO.SegmentHyperplaneRep{1,1}(
        [PLO.AffineFunction{1}((1.0,), -1.0), PLO.AffineFunction{1}((-1.0,), 2.0)],
        (PLO.AffineFunction{1}((1.0,), 1.5),),
    )
    # Segment 2: domain x >= 2 and x <= 3, output y = -2.5x + 8.5
    s2 = PLO.SegmentHyperplaneRep{1,1}(
        [PLO.AffineFunction{1}((1.0,), -2.0), PLO.AffineFunction{1}((-1.0,), 3.0)],
        (PLO.AffineFunction{1}((-2.5,), 8.5),),
    )
    pwl = PLO.PWLFunction{1,1,PLO.SegmentHyperplaneRep{1,1}}(
        [s1, s2], PLO.Intervals()
    )

    model = Model(optimizer)
    @variable(model, 1 <= x <= 3)
    y = piecewiselinear(model, (x,), pwl; method = method)
    @objective(model, Min, y[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    # Min of y = x + 1.5 on [1,2] is 2.5 at x=1
    # Min of y = -2.5x + 8.5 on [2,3] is 1.0 at x=3
    @test value(y[1]) ≈ 1.0 rtol = 1e-4
    @test value(x) ≈ 3.0 rtol = 1e-4
end

# Test direction parameter (Epigraph and Hypograph)
@testset "Direction tests" begin
    s1 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())

    model = Model(optimizer)
    @variable(model, x)
    y_epi = piecewiselinear(model, (x,), pwl; direction = PLO.Epigraph)
    @objective(model, Min, y_epi[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL

    model = Model(optimizer)
    @variable(model, x)
    y_hypo = piecewiselinear(model, (x,), pwl; direction = PLO.Hypograph)
    @objective(model, Max, y_hypo[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
end

# Test error cases
@testset "Error cases" begin
    # Mismatched input/output in SegmentPointRep
    @test_throws ErrorException PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,)])

    # Empty segments error
    empty_pwl = PLO.PWLFunction(PLO.SegmentPointRep{1,1}[], PLO.Intervals())
    model = Model(optimizer)
    @variable(model, x)
    @test_throws ErrorException piecewiselinear(model, (x,), empty_pwl)

    # Test output_vars kwarg
    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y_out)
    s1 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    y = piecewiselinear(model, (x,), pwl; output_vars = (y_out,))
    @test y[1] === y_out

    # Test univariate output_var kwarg with function
    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y_out2)
    d = 0:0.1:1
    f = xi -> xi^2
    y = piecewiselinear(model, x, d, f; output_var = y_out2)
    @test y === y_out2

    # Test univariate output_var with fd (vector) version
    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y_out3)
    fd = [xi^2 for xi in d]
    y = piecewiselinear(model, x, d, fd; output_var = y_out3)
    @test y === y_out3

    # Test bivariate output_var kwarg
    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    @variable(model, z_out)
    d2 = 0:0.5:1
    f2 = (xi, yi) -> xi + yi
    z = piecewiselinear(model, x, y, d2, d2, f2; output_var = z_out)
    @test z === z_out

    # Test formulate_pwl! fallback error
    struct FakeMethod <: PLO.Method end
    model = Model(optimizer)
    @variable(model, x)
    s1 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    @test_throws MethodError piecewiselinear(model, (x,), pwl; method = FakeMethod())
end

# Test internal utility functions
@testset "Utility functions" begin
    # Test _canonical!
    v = [3.0, 4.0]
    result = PLO._canonical!(v)
    @test result[1] ≈ 0.6 atol = 1e-6
    @test result[2] ≈ 0.8 atol = 1e-6

    # Test with negative leading element
    v2 = [-3.0, 4.0]
    result2 = PLO._canonical!(v2)
    @test result2[1] > 0  # should be positive after sign flip

    # Test with near-zero elements
    v3 = [1e-10, 1.0]
    result3 = PLO._canonical!(v3)
    @test result3[1] == 0.0
    @test result3[2] ≈ 1.0 atol = 1e-6

    # Test _compute_hyperplanes with k=2
    C = [[0, 0], [1, 0], [1, 1], [0, 1]]
    hyperplanes = PLO._compute_hyperplanes(C)
    @test length(hyperplanes) >= 1

    # Test _compute_hyperplanes with k > 2
    C3 = [[0, 0, 0], [1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 1, 1]]
    hyperplanes3 = PLO._compute_hyperplanes(C3)
    @test length(hyperplanes3) >= 1

    # Test _compute_hyperplanes error with k <= 1
    @test_throws ErrorException PLO._compute_hyperplanes([[1], [2], [3]])

    # Test _continuous_gridpoints_or_die error: discontinuous function
    s1 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(999.0,), (1.0,)])
    pwl_disc = PLO.PWLFunction([s1, s2], PLO.Intervals())
    @test_throws ErrorException PLO._continuous_gridpoints_or_die(pwl_disc)

    # Test Grid size mismatch error
    input_vals = Array{NTuple{1,Float64},1}([(1.0,), (2.0,)])
    output_vals = Array{NTuple{1,Float64},1}([(1.0,)])
    @test_throws ErrorException PLO.Grid(input_vals, output_vals)
end

include("aqua.jl")