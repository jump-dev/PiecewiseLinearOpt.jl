# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Formulations of bivariate functions on grid triangulations

@testset "Simple bivariate: $method" for method in [
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
    pwl = PLO.PWLFunction([s1, s2], PLO.UnstructuredTriangulation())

    y, _ = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    @objective(model, Min, y[1])

    solve!(model)

    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.0 rtol = 1e-4
    @test value(x[2]) ≈ 0.0 rtol = 1e-4
    @test value(y[1]) ≈ 0.0 rtol = 1e-4
end

@testset "Bivariate convenience methods" begin
    d = 0:0.05:1
    f = (xi, yi) -> xi^2 + yi^2
    pwl = BivariatePWLFunction(d, d, f)

    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    z1, _ = piecewiselinear(model, (x, y), pwl)
    @constraint(model, x ≤ 0.75)
    @constraint(model, y ≤ 0.75)
    @objective(model, Max, z1[1])
    solve!(model)
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
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 0.75 rtol = 1e-4
    @test value(y) ≈ 0.75 rtol = 1e-4
    @test value(z2) ≈ 1.125 rtol = 1e-4
end

@testset "2D: $method, $triangulation" for (method, triangulation) in [
    vec(collect(Iterators.product(methods_2D_gen, triangulations)))
    [(K1(m), K1Triangulation()) for m in axis_methods]
    [(UnionJack(m), UnionJackTriangulation()) for m in axis_methods]
]
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
        triangulation = triangulation,
    )
    @objective(model, Min, z)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.285714 rtol = 1e-4
    @test value(x[2]) ≈ 0.571429 rtol = 1e-4
    @test value(z) ≈ 0.004535 rtol = 1e-3
    @test objective_value(model) ≈ 0.004535 rtol = 1e-3
    @test objective_value(model) ≈ value(z) rtol = 1e-3

    @constraint(model, x[1] ≥ 0.6)
    solve!(model)

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
function test_triangle_selection(method, triangulation, n)
    d = range(0; stop = 1, length = n)
    f = (a, b) -> a^2 + 3b^2 - a * b + sin(3a)
    pwl = BivariatePWLFunction(d, d, f; triangulation = triangulation)
    model = Model(optimizer)
    @variable(model, x[1:2])
    (z,), _ = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    for seg in pwl.segments
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        z0 = sum(v[1] for v in seg.output_vals) / 3
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, z)
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(z) ≈ z0 atol = 1e-5
        end
    end
    return
end

@testset "Triangle selection: $method, $triangulation" for method in
                                                           methods_2D_selection,
    triangulation in triangulations

    test_triangle_selection(method, triangulation, 4)
end

@testset "Triangle selection: $method, $n points" for (method, triangulation) in
                                                      [
        (K1(), K1Triangulation()),
        (UnionJack(), UnionJackTriangulation()),
    ],
    n in (4, 5)

    test_triangle_selection(method, triangulation, n)
end

# The cover computed for OptimalIndependentBranching gets expensive quickly
@testset "Triangle selection: OptimalIndependentBranching, $triangulation" for triangulation in
                                                                               triangulations
    test_triangle_selection(
        OptimalIndependentBranching(optimizer),
        triangulation,
        3,
    )
end

@testset "Triangulation checks" begin
    d = range(0; stop = 1, length = 4)
    # K1 and UnionJack check the triangulation itself, not how it was built
    for (method, triangulation) in (
        (K1(), UnionJackTriangulation()),
        (K1(), UpperTriangulation()),
        (UnionJack(), K1Triangulation()),
    )
        pwl = BivariatePWLFunction(
            d,
            d,
            (a, b) -> a * b;
            triangulation = triangulation,
        )
        model = Model()
        @variable(model, x[1:2])
        test_rejected(model, (x[1], x[2]), pwl; method = method)
    end
    # For a * b, the lower triangulation splits every cell like K1
    pwl = BivariatePWLFunction(
        d,
        d,
        (a, b) -> a * b;
        triangulation = LowerTriangulation(),
    )
    model = Model()
    @variable(model, x[1:2])
    _, formulation = piecewiselinear(model, x, pwl; method = K1())
    @test formulation.method == K1()
end

@testset "Grid geometry checks: $method" for method in [
    SixStencil(),
    OptimalTriangleSelection(optimizer),
    OptimalIndependentBranching(optimizer),
]
    f = (a, b) -> a * b
    function grid_function(triangles)
        segments = [
            SegmentPointRep{2,1}(t, [(f(v...),) for v in t]) for t in triangles
        ]
        return PWLFunction(segments, UnstructuredTriangulation())
    end
    # A valid triangulation of a 3×2 grid, with triangles spanning two cells
    spanning = grid_function([
        [(0.0, 0.0), (1.0, 0.0), (2.0, 1.0)],
        [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0, 0.0), (1.0, 1.0), (2.0, 1.0)],
    ])
    # The first cell is split along both diagonals
    crossing = grid_function([
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)],
        [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
        [(1.0, 0.0), (1.0, 1.0), (2.0, 1.0)],
    ])
    # The first cell has three triangles and the second only one
    unbalanced = grid_function([
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
    ])
    for pwl in (spanning, crossing, unbalanced)
        model = Model()
        @variable(model, x[1:2])
        test_rejected(model, x, pwl; method = method)
    end
    # The general methods formulate the spanning triangulation correctly
    for general in (ConvexCombination(), MultipleChoice())
        model = Model(optimizer)
        @variable(model, x[1:2])
        (z,), _ = piecewiselinear(model, x, spanning; method = general)
        fix(x[1], 1.0)
        fix(x[2], 0.5)
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, z)
            solve!(model)
            @test value(z) ≈ spanning(1.0, 0.5) atol = 1e-5
        end
    end
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
    pwl = PLO.PWLFunction([s1, s2], PLO.UnstructuredTriangulation())
    model = Model(optimizer)
    @variable(model, x[1:2])
    y, _ = piecewiselinear(model, (x[1], x[2]), pwl; method = MultipleChoice())
    @objective(model, Max, y[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 1.0 atol = 1e-6
    @test value(x[2]) ≈ 1.0 atol = 1e-6
    @test value(y[1]) ≈ 7.0 atol = 1e-6
    # On the shared edge the lower triangle (s2) gives the larger value
    @constraint(model, x[1] == 0.5)
    @constraint(model, x[2] == 0.5)
    solve!(model)
    @test value(y[1]) ≈ 6.0 atol = 1e-6
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
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(z) ≈ f(1, 1) atol = 1e-5
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
        ((Epigraph(), MIN_SENSE), (Hypograph(), MAX_SENSE))
        model = Model(optimizer)
        @variable(model, x[1:2])
        (z,), _ = piecewiselinear(
            model,
            (x[1], x[2]),
            pwl;
            method = method,
            direction = direction,
        )
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        @objective(model, sense, z)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(z) ≈ z0 atol = 1e-5
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
    z, _ = piecewiselinear(model, (x[1], x[2]), pwl; method = method)
    for seg in segments
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        for k in 1:2, sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, z[k])
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(z[k]) ≈ sum(v[k] for v in seg.output_vals) / 3 atol =
                1e-5
        end
    end
end

# Several functions in one model: all methods must agree on the optimum
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
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
end
