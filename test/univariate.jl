# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Formulations of univariate functions

@testset "Simple univariate: $method" for method in methods_1D
    model = Model(optimizer)
    @variable(model, x)

    s1 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.5,), (1.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())

    y, _ = piecewiselinear(model, (x,), pwl; method = method)
    @objective(model, Min, y[1])

    solve!(model)
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
                solve!(model)
                @test termination_status(model) == MOI.OPTIMAL
                @test value(y) ≈ y0 atol = 1e-5
            end
        end
    end
end

@testset "Discontinuous univariate: $method" for method in methods_discontinuous
    model = Model(optimizer)
    @variable(model, x)
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    y, _ = piecewiselinear(model, (x,), pwl; method = method)
    @objective(model, Max, y[1])
    solve!(model)
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
    (y,), _ = piecewiselinear(model, (x,), pwl; method = method)
    for ((a, b), (c, e)) in zip(xs, ys)
        fix(x, 0.3 * a + 0.7 * b)
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y)
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ 0.3 * c + 0.7 * e atol = 1e-5
        end
    end
    # Points in the gaps of the domain are infeasible
    fix(x, 3.2)
    solve!(model)
    @test termination_status(model) == MOI.INFEASIBLE
end

@testset "Univariate convenience methods" begin
    d = 0:0.01:1
    f = (xi -> xi^2)
    fd = [f(xi) for xi in d]
    pwl = UnivariatePWLFunction(d, f)

    model = Model(optimizer)
    @variable(model, x)
    y1, _ = piecewiselinear(model, (x,), pwl)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y1[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y1[1]) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(optimizer)
    @variable(model, x)
    y2 = piecewiselinear(model, x, d, f)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y2)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y2) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(optimizer)
    @variable(model, x)
    y3 = piecewiselinear(model, x, d, fd)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y3)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y3) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4
end

@testset "1D: $method" for method in methods_1D
    model = Model(optimizer)
    @variable(model, x)
    d = 7
    xs = collect(range(1; stop = 2π, length = (d + 1)))
    zs = sin.(xs)
    y = piecewiselinear(model, x, xs, zs; method = method)
    @objective(model, Max, y)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.75474 rtol = 1e-4
    @test value(y) ≈ 0.98313 rtol = 1e-4
    @test objective_value(model) ≈ 0.98313 rtol = 1e-4
    @test objective_value(model) ≈ value(y) rtol = 1e-4
    @constraint(model, x ≤ 1.5y)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 1.36495 rtol = 1e-4
    @test value(y) ≈ 0.90997 rtol = 1e-4
    @test objective_value(model) ≈ 0.90997 rtol = 1e-4
    @test objective_value(model) ≈ value(y) rtol = 1e-4
end

@testset "Incremental size" begin
    for d in 1:5
        model = Model()
        @variable(model, x)
        pwl = UnivariatePWLFunction(0:d, t -> t^2)
        _, formulation =
            piecewiselinear(model, [x], pwl; method = Incremental())
        # d fill variables, d - 1 binary variables and the output
        @test count(is_binary, formulation.variables) == d - 1
        @test length(formulation.variables) == 2d
    end
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
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ y0 atol = 1e-5
        end
    end
end

@testset "Direction: $method" for method in methods_1D
    xs = [0.0, 1.0, 2.0, 3.0]
    zs = [0.0, 2.0, 1.0, 3.0]
    for (direction, sense, beyond) in (
        (PLO.Epigraph, MIN_SENSE, maximum(zs) + 10),
        (PLO.Hypograph, MAX_SENSE, minimum(zs) - 10),
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
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ 1.5 atol = 1e-5
        # In the other direction, y is unbounded, also beyond the range of
        # the function
        @objective(model, sense == MIN_SENSE ? MAX_SENSE : MIN_SENSE, y)
        solve!(model)
        @test termination_status(model) in
              (MOI.DUAL_INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        @constraint(model, y == beyond)
        @objective(model, sense, y)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ beyond atol = 1e-5
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
    y, _ = piecewiselinear(model, (x,), pwl; method = method)
    @test length(y) == 2
    for x0 in (0.5, 1.25, 2.9)
        fix(x, x0)
        i = floor(Int, x0) + 1
        t = x0 - xs[i]
        for k in 1:2, sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y[k])
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y[k]) ≈ (1 - t) * fs[k][i] + t * fs[k][i+1] atol = 1e-5
        end
    end
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
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
end
