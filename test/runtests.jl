# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

using PiecewiseLinearOpt
using HiGHS
using JuMP
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
    IndependentBranching(:biclique),
    IndependentBranching(:balanced),
    IndependentBranching(:brgc),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    MultipleChoice(),
    NativeSOS2(),
    SymmetricCelaya(),
    ZigZagBinary(),
    ZigZagInteger(),
]

# Checks that the levels (L, R) form a biclique cover of the conflict graph of
# the CDC given by S: no level separates two points of the same set, and any
# two points that are not in a common set are separated by some level.
function is_valid_cover(S, cover)
    J = reduce(union, S)
    for (L, R) in cover
        if !isdisjoint(L, R) || !issubset(union(L, R), J)
            return false
        end
        for s in S
            if !isdisjoint(s, L) && !isdisjoint(s, R)
                return false
            end
        end
    end
    for u in J, v in J
        if u < v && !any(s -> u in s && v in s, S)
            if !any(
                ((L, R),) -> (u in L && v in R) || (v in L && u in R),
                cover,
            )
                return false
            end
        end
    end
    return true
end

@testset "Ordered CDC covers: $scheme" for scheme in
                                           (:biclique, :balanced, :brgc)
    for d in 1:20
        # SOS2
        S = [[i, i + 1] for i in 1:d]
        cover = PLO._ordered_cdc_cover(S, scheme)
        @test is_valid_cover(S, cover)
        @test length(cover) == (d == 1 ? 0 : ceil(Int, log2(d)))
        # Polytopes sharing none, one or two points with their neighbors
        for shared in (0, 1, 2)
            S, next = Vector{Int}[], 1
            for i in 1:d
                first = i == 1 ? next : next - shared
                push!(S, collect(first:first+2*shared+mod(i, 3)))
                next = S[end][end] + 1
            end
            cover = PLO._ordered_cdc_cover(S, scheme)
            @test is_valid_cover(S, cover)
            @test length(cover) == (d == 1 ? 0 : ceil(Int, log2(d)))
        end
    end
    # Remark 4: the first part of each cut gets floor(n / 2) of the n sets
    @test PLO._balanced_reversed_edge_ranking(6) == [2, 3, 1, 2, 3]
end

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

const methods_discontinuous = [
    IndependentBranching(:biclique),
    IndependentBranching(:balanced),
    IndependentBranching(:brgc),
    MultipleChoice(),
]

@testset "Discontinuous univariate: $method" for method in methods_discontinuous
    model = Model(optimizer)
    @variable(model, x)
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    y = piecewiselinear(model, (x,), pwl; method = method)
    @objective(model, Max, y[1])
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.0 atol = 1e-6
    @test value(y[1]) ≈ 3.0 atol = 1e-6

    # A function with jumps, gaps in the domain and continuous stretches.
    # Inside each segment, y must equal the interpolated value.
    xs = [(0, 1), (1, 2), (2, 3), (3.5, 4), (4, 5), (5, 5.5), (6, 7), (7, 8)]
    ys = [(0, 2), (1, 0), (0, 3), (3, -1), (2, 2), (2, 4), (4, 1), (1, 0)]
    segments = [
        PLO.SegmentPointRep{1,1}(
            [(float(a),), (float(b),)],
            [(float(c),), (float(e),)],
        ) for ((a, b), (c, e)) in zip(xs, ys)
    ]
    pwl = PLO.PWLFunction(segments, PLO.Intervals())
    model = Model(optimizer)
    @variable(model, x)
    y = piecewiselinear(model, (x,), pwl; method = method)[1]
    for ((a, b), (c, e)) in zip(xs, ys)
        fix(x, 0.3 * a + 0.7 * b)
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y)
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ 0.3 * c + 0.7 * e atol = 1e-5
        end
    end
    # Points in the gaps of the domain are infeasible
    fix(x, 3.2)
    optimize!(model)
    @test termination_status(model) == MOI.INFEASIBLE
end

@testset "Independent branching: $scheme" for scheme in
                                              (:biclique, :balanced, :brgc)
    # ceil(log2(d)) binaries for d segments, also for discontinuous functions
    for d in (2, 3, 5, 8, 13), jump in (0.0, 1.0)
        segments = [
            PLO.SegmentPointRep{1,1}(
                [(float(i),), (i + 1.0,)],
                [(float(i),), (i + 1.0 + jump * isodd(i),)],
            ) for i in 1:d
        ]
        pwl = PLO.PWLFunction(segments, PLO.Intervals())
        model = Model()
        @variable(model, x)
        piecewiselinear(model, (x,), pwl; method = IndependentBranching(scheme))
        @test num_constraints(model, VariableRef, MOI.ZeroOne) ==
              ceil(Int, log2(d))
    end
    # Segments must be ordered so that only neighbors share points
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (2.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (3.0,)], [(1.0,), (0.0,)])
    s3 = PLO.SegmentPointRep{1,1}([(2.0,), (4.0,)], [(1.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2, s3], PLO.Intervals())
    model = Model()
    @variable(model, x)
    @test_throws ErrorException piecewiselinear(
        model,
        (x,),
        pwl;
        method = IndependentBranching(scheme),
    )
    @test_throws ErrorException IndependentBranching(:unknown)
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

const patterns = [:Upper, :Lower, :BestFit, :K1, :UnionJack, :Random]
const method_pattern = [
    vec(collect(Iterators.product(methods_2D_gen, patterns)))
    [(K1(sos2_method), :K1) for sos2_method in methods_1D]
    [(UnionJack(sos2_method), :UnionJack) for sos2_method in methods_1D]
]

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

@testset "Direction: $method" for method in methods_1D
    xs = [0.0, 1.0, 2.0, 3.0]
    zs = [0.0, 2.0, 1.0, 3.0]
    for (direction, sense, other_bound) in (
        (PLO.Epigraph, MIN_SENSE, maximum(zs)),
        (PLO.Hypograph, MAX_SENSE, minimum(zs)),
    )
        model = Model(optimizer)
        @variable(model, x)
        y = piecewiselinear(
            model,
            x,
            xs,
            zs;
            method = method,
            direction = direction,
        )
        fix(x, 1.5)
        @objective(model, sense, y)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ 1.5 atol = 1e-5
        # In the other direction, y is only limited by its bounds
        @objective(model, sense == MIN_SENSE ? MAX_SENSE : MIN_SENSE, y)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ other_bound atol = 1e-5
    end
end

@testset "Direction: $method" for method in [
    ConvexCombination(),
    MultipleChoice(),
    OptimalTriangleSelection(optimizer),
    SixStencil(),
]
    d = 0:2
    f = (a, b) -> a^2 - b
    pwl = BivariatePWLFunction(d, d, f)
    seg = pwl.segments[3]
    z0 = sum(v[1] for v in seg.output_vals) / 3
    for (direction, sense) in
        ((PLO.Epigraph, MIN_SENSE), (PLO.Hypograph, MAX_SENSE))
        model = Model(optimizer)
        @variable(model, x[1:2])
        z = piecewiselinear(
            model,
            (x[1], x[2]),
            pwl;
            method = method,
            direction = direction,
        )[1]
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        @objective(model, sense, z)
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(z) ≈ z0 atol = 1e-5
    end
end

@testset "Multiple outputs: $method" for method in methods_1D
    xs = [0.0, 1.0, 2.0, 3.0]
    fs = ([0.0, 2.0, 1.0, 3.0], [5.0, -1.0, 4.0, 0.0])
    segments = [
        PLO.SegmentPointRep{1,2}(
            [(xs[i],), (xs[i+1],)],
            [(fs[1][i], fs[2][i]), (fs[1][i+1], fs[2][i+1])],
        ) for i in 1:3
    ]
    pwl = PLO.PWLFunction(segments, PLO.Intervals())
    model = Model(optimizer)
    @variable(model, x)
    y = piecewiselinear(model, (x,), pwl; method = method)
    @test length(y) == 2
    for x0 in (0.5, 1.25, 2.9)
        fix(x, x0)
        i = floor(Int, x0) + 1
        t = x0 - xs[i]
        for k in 1:2, sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y[k])
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y[k]) ≈ (1 - t) * fs[k][i] + t * fs[k][i+1] atol = 1e-5
        end
    end
end

@testset "Multiple outputs: $method" for method in [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    K1(),
    MultipleChoice(),
    OptimalTriangleSelection(optimizer),
    SixStencil(),
]
    d = range(0; stop = 1, length = 3)
    f = (a, b) -> a^2 + 3b^2 - a * b
    g = (a, b) -> sin(3a) - b
    base = BivariatePWLFunction(d, d, f)
    segments = [
        PLO.SegmentPointRep{2,2}(
            s.input_vals,
            [(f(v...), g(v...)) for v in s.input_vals],
        ) for s in base.segments
    ]
    pwl = PLO.PWLFunction(segments, base.structure)
    model = Model(optimizer)
    @variable(model, x[1:2])
    z = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    for seg in segments
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        for k in 1:2, sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, z[k])
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(z[k]) ≈ sum(v[k] for v in seg.output_vals) / 3 atol =
                1e-5
        end
    end
end

@testset "Trivariate: $method" for method in [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
]
    # Kuhn triangulation of the unit cube into 6 tetrahedra
    g = v -> v[1] + 2v[2]^2 - 3v[1] * v[3]
    segments = PLO.SegmentPointRep{3,1}[]
    for σ in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
        vertices = [(0.0, 0.0, 0.0)]
        for j in σ
            push!(vertices, ntuple(i -> i == j ? 1.0 : vertices[end][i], 3))
        end
        push!(
            segments,
            PLO.SegmentPointRep{3,1}(vertices, [(g(v),) for v in vertices]),
        )
    end
    pwl = PLO.PWLFunction(segments, PLO.Unstructured{3}())
    model = Model(optimizer)
    @variable(model, x[1:3])
    y = piecewiselinear(model, (x[1], x[2], x[3]), pwl; method = method)[1]
    for seg in segments
        for i in 1:3
            fix(x[i], sum(v[i] for v in seg.input_vals) / 4)
        end
        y0 = sum(v[1] for v in seg.output_vals) / 4
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y)
            optimize!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ y0 atol = 1e-5
        end
    end
end

@testset "Output variable" begin
    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    ret = piecewiselinear(model, x, [0, 1, 2], [0, 2, 1]; output_var = y)
    @test ret === y
    fix(x, 1.5)
    @objective(model, Max, y)
    optimize!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y) ≈ 1.5 atol = 1e-5
end

@testset "Errors" begin
    model = Model()
    @variable(model, x[1:3])
    # No segments
    pwl = PLO.PWLFunction(PLO.SegmentPointRep{1,1}[], PLO.Intervals())
    @test_throws ErrorException piecewiselinear(model, (x[1],), pwl)
    # Discontinuous functions require MultipleChoice
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    for method in (ConvexCombination(), Logarithmic(), Incremental())
        @test_throws ErrorException piecewiselinear(
            model,
            (x[1],),
            pwl;
            method = method,
        )
    end
    # Grid based methods require a connected domain
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    @test_throws ErrorException piecewiselinear(
        model,
        (x[1],),
        pwl;
        method = Logarithmic(),
    )
    # Unsupported combinations of dimension and method
    pwl = BivariatePWLFunction(0:1, 0:1, (a, b) -> a + b)
    @test_throws ErrorException piecewiselinear(
        model,
        (x[1], x[2]),
        pwl;
        method = Incremental(),
    )
    s = PLO.SegmentPointRep{3,1}(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (3.0,)],
    )
    pwl = PLO.PWLFunction([s], PLO.Unstructured{3}())
    @test_throws ErrorException piecewiselinear(
        model,
        (x[1], x[2], x[3]),
        pwl;
        method = SixStencil(),
    )
end

# Several functions in one model: all methods must agree on the optimum
@testset "Methods agree (univariate)" begin
    supply, demand = [3.0, 4.0], [2.0, 2.5, 2.5]
    xs = 0:4
    cost = [[0.0, 2.0, 3.0, 3.5, 3.8], [0.0, 1.0, 2.5, 3.0, 3.2]]
    objectives = map(methods_1D) do method
        model = Model(optimizer)
        @variable(model, x[1:2, 1:3] ≥ 0)
        @constraint(model, [i in 1:2], sum(x[i, :]) == supply[i])
        @constraint(model, [j in 1:3], sum(x[:, j]) == demand[j])
        @objective(
            model,
            Min,
            sum(
                piecewiselinear(model, x[i, j], xs, cost[i]; method = method) for i in 1:2, j in 1:3
            )
        )
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
end

@testset "Methods agree (bivariate)" begin
    d = range(0; stop = 1, length = 3)
    f = (a, b) -> (a - 0.4)^2 - sin(3b) + a * b
    methods = [
        methods_2D_selection
        OptimalIndependentBranching(optimizer)
        K1()
    ]
    objectives = map(methods) do method
        model = Model(optimizer)
        @variable(model, 0 ≤ a[1:2] ≤ 1)
        @variable(model, 0 ≤ b[1:2] ≤ 1)
        @constraint(model, sum(a) == 1.1)
        @constraint(model, sum(b) == 0.7)
        @objective(
            model,
            Min,
            sum(
                piecewiselinear(model, a[i], b[i], d, d, f; method = method) for i in 1:2
            )
        )
        optimize!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
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
