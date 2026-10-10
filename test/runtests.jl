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

# With the default tolerances, HiGHS occasionally reports a solve error because
# rounding a binary within its integrality tolerance (1e-6) leaves a constraint
# violated by more than the primal feasibility tolerance.
optimizer = optimizer_with_attributes(
    HiGHS.Optimizer,
    MOI.Silent() => true,
    "mip_feasibility_tolerance" => 1e-8,
    "primal_feasibility_tolerance" => 1e-8,
)

const methods_1D = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    GeneralizedCelaya(),
    Incremental(),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    MultipleChoice(),
    NativeSOS2(),
    SymmetricCelaya(),
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

@testset "Univariate sizes: $method" for method in methods_1D
    # The encodings of the SOS2 formulations depend on the number of segments.
    # Fixing x inside a segment, y must be pinned to the interpolated value
    # (a weak formulation would let y range over the convex hull).
    for n in 2:10
        xs = collect(range(0; stop = 2π, length = n))
        zs = sin.(xs)
        model = Model(optimizer)
        @variable(model, x)
        y = piecewiselinear(model, x, xs, zs; method = method)
        for i in 1:n-1
            x0 = 0.3 * xs[i] + 0.7 * xs[i+1]
            y0 = 0.3 * zs[i] + 0.7 * zs[i+1]
            fix(x, x0)
            for sense in (MAX_SENSE, MIN_SENSE)
                @objective(model, sense, y)
                optimize!(model)
                @test termination_status(model) == MOI.OPTIMAL
                @test value(y) ≈ y0 atol = 1e-5
            end
        end
    end
end

@testset "Discontinuous univariate" begin
    model = Model(optimizer)
    @variable(model, x)
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    y = piecewiselinear(model, (x,), pwl; method = MultipleChoice())
    @objective(model, Max, y[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.0 atol = 1e-6
    @test value(y[1]) ≈ 3.0 atol = 1e-6
end

@testset "Discontinuous bivariate" begin
    # Two triangles sharing the diagonal from (0, 0) to (1, 1), with a jump
    # across it.
    s1 = PLO.SegmentPointRep{2,1}(
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0,), (1.0,), (2.0,)],
    )
    s2 = PLO.SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(5.0,), (3.0,), (7.0,)],
    )
    pwl = PLO.PWLFunction{2,1,PLO.SegmentPointRep{2,1}}(
        [s1, s2],
        PLO.UnstructuredTriangulation(),
    )
    model = Model(optimizer)
    @variable(model, x[1:2])
    y = piecewiselinear(model, (x[1], x[2]), pwl; method = MultipleChoice())
    @objective(model, Max, y[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 1.0 atol = 1e-6
    @test value(x[2]) ≈ 1.0 atol = 1e-6
    @test value(y[1]) ≈ 7.0 atol = 1e-6
    # On the shared edge the lower triangle (s2) gives the larger value
    @constraint(model, x[1] == 0.5)
    @constraint(model, x[2] == 0.5)
    optimize!(model)
    @test value(y[1]) ≈ 6.0 atol = 1e-6
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
    GeneralizedCelaya(),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    NativeSOS2(),
    SymmetricCelaya(),
    ZigZagBinary(),
    ZigZagInteger(),
]
const methods_2D_gen = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
    [NineStencil(sos2_method) for sos2_method in methods_1D]...,
    [
        OptimalTriangleSelection(optimizer, sos2_method) for
        sos2_method in methods_1D
    ]...,
    [SixStencil(sos2_method) for sos2_method in methods_1D]...,
]
@testset "Simple bivariate" for method in [
    methods_2D_gen
    OptimalIndependentBranching(optimizer)
]
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

# Fixing the input at the centroid of each triangle, the output must be pinned
# to the interpolated value. A formulation that selects the wrong triangle in
# some cell lets the output deviate from it.
function test_triangle_selection(method, pattern, n)
    d = range(0; stop = 1, length = n)
    f = (a, b) -> a^2 + 3b^2 - a * b + sin(3a)
    pwl = BivariatePWLFunction(d, d, f; pattern = pattern)
    model = Model(optimizer)
    @variable(model, x[1:2])
    z = piecewiselinear(model, (x[1], x[2]), pwl; method = method)[1]
    for seg in pwl.segments
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        z0 = sum(v[1] for v in seg.output_vals) / 3
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, z)
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(z) ≈ z0 atol = 1e-5
        end
    end
    return
end

const methods_2D_selection = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
    NineStencil(),
    OptimalTriangleSelection(optimizer),
    SixStencil(),
]

@testset "2D triangle selection: $method, $pattern" for method in
                                                        methods_2D_selection,
    pattern in patterns

    test_triangle_selection(method, pattern, 4)
end

@testset "2D triangle selection: $method" for (method, pattern) in [
        (K1(), :K1),
        (UnionJack(), :UnionJack),
    ],
    n in (4, 5)

    test_triangle_selection(method, pattern, n)
end

# The cover computed for OptimalIndependentBranching gets expensive quickly
@testset "2D triangle selection: OptimalIndependentBranching, $pattern" for pattern in
                                                                            patterns
    test_triangle_selection(OptimalIndependentBranching(optimizer), pattern, 3)
end

@testset "Triangulation checks" begin
    d = range(0; stop = 1, length = 4)
    f = (a, b) -> a + b
    for (method, pattern, structure) in (
        (K1(), :UnionJack, PLO.K1Triangulation()),
        (UnionJack(), :K1, PLO.UnionJackTriangulation()),
    )
        pwl = BivariatePWLFunction(d, d, f; pattern = pattern)
        pwl = PLO.PWLFunction(pwl.segments, structure)
        model = Model()
        @variable(model, x[1:2])
        @test_throws ErrorException piecewiselinear(
            model,
            (x[1], x[2]),
            pwl;
            method = method,
        )
    end
    # K1 and UnionJack require the matching structure
    model = Model()
    @variable(model, x[1:2])
    pwl = BivariatePWLFunction(d, d, f; pattern = :Upper)
    @test_throws ErrorException piecewiselinear(
        model,
        (x[1], x[2]),
        pwl;
        method = K1(),
    )
end

@testset "Input validation" begin
    @test_throws ErrorException UnivariatePWLFunction([0, 2, 1], [0, 1, 2])
    @test_throws ErrorException UnivariatePWLFunction([0, 1, 1], [0, 1, 2])
    @test_throws ErrorException BivariatePWLFunction(
        [0, 1],
        [1, 0],
        (a, b) -> a + b,
    )
    @test_throws ErrorException BivariatePWLFunction(
        [0, 1],
        [0, 1],
        (a, b) -> a + b;
        pattern = :Unknown,
    )
    # Segments that are not simplices cannot be used with MultipleChoice
    s = PLO.SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (1.0,)],
    )
    pwl = PLO.PWLFunction([s], PLO.UnstructuredTriangulation())
    model = Model()
    @variable(model, x[1:2])
    @test_throws ErrorException piecewiselinear(
        model,
        (x[1], x[2]),
        pwl;
        method = MultipleChoice(),
    )
end

@testset "Affine inputs: $method" for method in methods_1D
    # y = f(2x + 1), with breakpoints of f at 1:5
    xs = 1:5
    zs = [0.0, 2.0, 1.0, 3.0, -1.0]
    model = Model(optimizer)
    @variable(model, x)
    y = piecewiselinear(model, 2x + 1, xs, zs; method = method)
    for (x0, y0) in ((0.25, 0.5 * 0.0 + 0.5 * 2.0), (1.75, 0.5 * 3.0 - 0.5))
        fix(x, x0)
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y)
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ y0 atol = 1e-5
        end
    end
end

@testset "Affine inputs: $method" for method in
                                      [MultipleChoice(), SixStencil(), K1()]
    d = 0:2
    f = (a, b) -> a^2 - b
    model = Model(optimizer)
    @variable(model, x[1:2])
    z = piecewiselinear(model, 2x[1], x[2] + 1, d, d, f; method = method)
    fix(x[1], 0.5)
    fix(x[2], 0.0)
    for sense in (MAX_SENSE, MIN_SENSE)
        @objective(model, sense, z)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(z) ≈ f(1, 1) atol = 1e-5
    end
end
