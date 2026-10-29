# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Helpers shared by the test files. Each test file includes this file unless it
# is already loaded, so that the file can also be run on its own.

module PLOTestUtilities

import HiGHS

using JuMP
using PiecewiseLinearOpt
using Test

export FULL,
    OPTIMIZER,
    METHODS_1D,
    AXIS_METHODS,
    TRIANGULATIONS,
    solve!,
    nonbound_constraints,
    check_rejected,
    check_formulation,
    check_univariate,
    check_discontinuous,
    check_relaxation,
    check_bivariate_simple,
    check_bivariate_grid,
    check_triangle_selection,
    check_bivariate_formulation,
    check_bivariate_affine,
    check_bivariate_direction,
    check_bivariate_multiple_outputs,
    grid_function,
    check_grid_geometry,
    check_spanning_triangulation,
    check_trivariate

# The standard run combines the bivariate grid methods only with a few
# representative axis methods. The full run uses every combination, and is set
# by runtests.jl for the --full argument.
const FULL = get(ENV, "PLO_TEST_FULL", "false") == "true"

const OPTIMIZER =
    optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

# HiGHS occasionally reports a solve error on these small models, when the
# postsolved solution violates a constraint by exactly the MIP feasibility
# tolerance, whatever its value. Solving again without presolve avoids it.
# Presolve stays off for the model, as changing it back would discard the
# solution.
function solve!(model)
    optimize!(model)
    if termination_status(model) == MOI.OTHER_ERROR
        set_attribute(model, "presolve", "off")
        optimize!(model)
    end
    return
end

# All constraints except bounds and integrality restrictions of variables
function nonbound_constraints(model)
    bounds = (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.ZeroOne,
        MOI.Integer,
    )
    constraints = ConstraintRef[]
    for (F, S) in list_of_constraint_types(model)
        if !(F == VariableRef && S in bounds)
            append!(constraints, all_constraints(model, F, S))
        end
    end
    return constraints
end

# Checks that adding pwl of inputs, formulated as given by kwargs, throws an
# error and leaves the model unchanged
function check_rejected(model, inputs, pwl; output_vars = nothing, kwargs...)
    variables = all_variables(model)
    constraints = nonbound_constraints(model)
    @test_throws ErrorException add_formulation(
        model,
        PiecewiseLinear(pwl; kwargs...),
        inputs...;
        output_vars = output_vars,
    )
    @test all_variables(model) == variables
    @test nonbound_constraints(model) == constraints
    return
end

const METHODS_1D = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    Incremental(),
    IndependentBranching(:biclique),
    IndependentBranching(:balanced),
    IndependentBranching(:brgc),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    MultipleChoice(),
    NativeSOS2(),
    ZigZagBinary(),
    ZigZagInteger(),
]

# The univariate methods used to select a cell in the bivariate grid methods.
# The representative ones cover an embedding formulation, an independent
# branching formulation, a native SOS2 constraint and the generic SOS2
# formulation through a dummy function.
const AXIS_METHODS =
    FULL ? METHODS_1D :
    [Logarithmic(), IndependentBranching(), NativeSOS2(), Incremental()]

const TRIANGULATIONS = [
    UpperTriangulation(),
    LowerTriangulation(),
    BestFitTriangulation(),
    K1Triangulation(),
    UnionJackTriangulation(),
    RandomTriangulation(),
]

# The formulation must record exactly what was added to a model that only
# contained the inputs, and deleting it must remove all of it.
function check_formulation(model, inputs, pwl, method)
    y, formulation =
        add_formulation(model, PiecewiseLinear(pwl; method = method), inputs...)
    @test formulation isa PWLFormulation
    @test formulation.method === method
    @test all(yi -> yi in formulation.variables, y)
    @test Set(formulation.variables) ==
          setdiff(Set(all_variables(model)), inputs)
    @test length(unique(formulation.variables)) == length(formulation.variables)
    @test Set(formulation.constraints) == Set(nonbound_constraints(model))
    @test all(
        v -> startswith(name(v), formulation.name * "_"),
        formulation.variables,
    )
    @test allunique(name.(formulation.variables))
    @test occursin(formulation.name, sprint(show, formulation))
    delete(model, formulation)
    @test Set(all_variables(model)) == Set(inputs)
    @test isempty(nonbound_constraints(model))
    return
end

# Continuous univariate functions

function check_univariate(method)
    @testset "Simple" begin
        _check_univariate_simple(method)
    end
    @testset "Sizes" begin
        _check_univariate_sizes(method)
    end
    @testset "Sine" begin
        _check_univariate_sine(method)
    end
    @testset "Affine inputs" begin
        _check_univariate_affine(method)
    end
    @testset "Direction" begin
        _check_univariate_direction(method)
    end
    @testset "Multiple outputs" begin
        _check_univariate_multiple_outputs(method)
    end
    @testset "Formulation" begin
        pwl = UnivariatePWLFunction([0, 1, 2, 4, 5, 7], [0, 2, 1, 3, -1, 0])
        model = Model()
        @variable(model, x)
        check_formulation(model, [x], pwl, method)
    end
    return
end

function _check_univariate_simple(method)
    model = Model(OPTIMIZER)
    @variable(model, x)
    s1 = SegmentPointRep{1,1}([(1.0,), (2.0,)], [(2.5,), (3.5,)])
    s2 = SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.5,), (1.0,)])
    pwl = PWLFunction([s1, s2], Intervals())
    y, _ = add_formulation(model, PiecewiseLinear(pwl; method = method), x)
    @objective(model, Min, y[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 3.0 rtol = 1e-4
    @test value(y[1]) ≈ 1.0 rtol = 1e-4
    return
end

# The encodings of the SOS2 formulations depend on the number of segments.
# Fixing x inside a segment, y must be pinned to the interpolated value (a weak
# formulation would let y range over the convex hull).
function _check_univariate_sizes(method)
    for n in 2:10
        xs = collect(range(0; stop = 2π, length = n))
        zs = sin.(xs)
        model = Model(OPTIMIZER)
        @variable(model, x)
        y = PiecewiseLinear(xs, zs; method = method)(x)
        for i in 1:(n-1)
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
    return
end

function _check_univariate_sine(method)
    model = Model(OPTIMIZER)
    @variable(model, x)
    d = 7
    xs = collect(range(1; stop = 2π, length = (d + 1)))
    zs = sin.(xs)
    y = PiecewiseLinear(xs, zs; method = method)(x)
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
    return
end

function _check_univariate_affine(method)
    # y = f(2x + 1), with breakpoints of f at 1:5
    xs = 1:5
    zs = [0.0, 2.0, 1.0, 3.0, -1.0]
    model = Model(OPTIMIZER)
    @variable(model, x)
    y = PiecewiseLinear(xs, zs; method = method)(2x + 1)
    for (x0, y0) in ((0.25, 0.5 * 0.0 + 0.5 * 2.0), (1.75, 0.5 * 3.0 - 0.5))
        fix(x, x0)
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y)
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ y0 atol = 1e-5
        end
    end
    return
end

function _check_univariate_direction(method)
    xs = [0.0, 1.0, 2.0, 3.0]
    zs = [0.0, 2.0, 1.0, 3.0]
    for (direction, sense, beyond) in (
        (Epigraph(), MIN_SENSE, maximum(zs) + 10),
        (Hypograph(), MAX_SENSE, minimum(zs) - 10),
    )
        model = Model(OPTIMIZER)
        @variable(model, x)
        y = PiecewiseLinear(xs, zs; method = method, direction = direction)(x)
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
    return
end

function _check_univariate_multiple_outputs(method)
    xs = [0.0, 1.0, 2.0, 3.0]
    fs = ([0.0, 2.0, 1.0, 3.0], [5.0, -1.0, 4.0, 0.0])
    segments = [
        SegmentPointRep{1,2}(
            [(xs[i],), (xs[i+1],)],
            [(fs[1][i], fs[2][i]), (fs[1][i+1], fs[2][i+1])],
        ) for i in 1:3
    ]
    pwl = PWLFunction(segments, Intervals())
    model = Model(OPTIMIZER)
    @variable(model, x)
    y, _ = add_formulation(model, PiecewiseLinear(pwl; method = method), x)
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
    return
end

# Discontinuous univariate functions

# A relaxation of t^2 with a band of ±0.1 around its interpolant: inside each
# interval, y ranges between the bounds
function check_relaxation(method)
    f = t -> t^2
    xs = 0:0.5:2
    pwl = UnivariatePWLRelaxation(xs, t -> f(t) - 0.1, t -> f(t) + 0.1)
    model = Model(OPTIMIZER)
    @variable(model, x)
    (y,), _ = add_formulation(model, PiecewiseLinear(pwl; method = method), x)
    for t in (0.0, 0.3, 0.5, 1.25, 2.0)
        fix(x, t)
        i = min(searchsortedlast(xs, t), length(xs) - 1)
        a, b = xs[i], xs[i+1]
        interpolant = f(a) + (f(b) - f(a)) * (t - a) / (b - a)
        @objective(model, Min, y)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ interpolant - 0.1 atol = 1e-5
        @objective(model, Max, y)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ interpolant + 0.1 atol = 1e-5
    end
    # Points outside the domain are infeasible
    fix(x, 2.5)
    solve!(model)
    @test termination_status(model) == MOI.INFEASIBLE
    # Where the bounds coincide, the relaxation has a single value
    pwl = UnivariatePWLRelaxation([0, 1, 2], [0, 1, 0], [1, 1, 2])
    model = Model(OPTIMIZER)
    @variable(model, x)
    fix(x, 1)
    (y,), _ = add_formulation(model, PiecewiseLinear(pwl; method = method), x)
    for sense in (MIN_SENSE, MAX_SENSE)
        @objective(model, sense, y)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(y) ≈ 1 atol = 1e-5
    end
    return
end

function check_discontinuous(method)
    @testset "Jump" begin
        model = Model(OPTIMIZER)
        @variable(model, x)
        s1 = SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
        s2 = SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
        pwl = PWLFunction([s1, s2], Intervals())
        y, _ = add_formulation(model, PiecewiseLinear(pwl; method = method), x)
        @objective(model, Max, y[1])
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(x) ≈ 1.0 atol = 1e-6
        @test value(y[1]) ≈ 3.0 atol = 1e-6
    end
    @testset "Jumps and gaps" begin
        # A function with jumps, gaps in the domain and continuous stretches.
        # Inside each segment, y must equal the interpolated value.
        xs =
            [(0, 1), (1, 2), (2, 3), (3.5, 4), (4, 5), (5, 5.5), (6, 7), (7, 8)]
        ys = [(0, 2), (1, 0), (0, 3), (3, -1), (2, 2), (2, 4), (4, 1), (1, 0)]
        segments = [
            SegmentPointRep{1,1}(
                [(float(a),), (float(b),)],
                [(float(c),), (float(e),)],
            ) for ((a, b), (c, e)) in zip(xs, ys)
        ]
        pwl = PWLFunction(segments, Intervals())
        model = Model(OPTIMIZER)
        @variable(model, x)
        (y,), _ =
            add_formulation(model, PiecewiseLinear(pwl; method = method), x)
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
    @testset "Formulation" begin
        s1 = SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
        s2 = SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
        s3 = SegmentPointRep{1,1}([(2.0,), (3.0,)], [(2.0,), (5.0,)])
        model = Model()
        @variable(model, x)
        pwl = PWLFunction([s1, s2, s3], Intervals())
        check_formulation(model, [x], pwl, method)
    end
    return
end

# Bivariate functions

function check_bivariate_simple(method)
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    s1 = SegmentPointRep{2,1}(
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0,), (1.0,), (2.0,)],
    )
    s2 = SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(0.0,), (3.0,), (2.0,)],
    )
    pwl = PWLFunction([s1, s2], UnstructuredTriangulation())
    y, _ = add_formulation(
        model,
        PiecewiseLinear(pwl; method = method),
        x[1],
        x[2],
    )
    @objective(model, Min, y[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 0.0 rtol = 1e-4
    @test value(x[2]) ≈ 0.0 rtol = 1e-4
    @test value(y[1]) ≈ 0.0 rtol = 1e-4
    return
end

function check_bivariate_grid(method, triangulation)
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    d = range(0; stop = 1, length = 8)
    f = (x1, x2) -> 2 * (x1 - 1 / 3)^2 + 3 * (x2 - 4 / 7)^4
    z =
        PiecewiseLinear(f, d, d; method = method, triangulation = triangulation)(
            x[1],
            x[2],
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
    return
end

# Fixing the input at the centroid of each triangle, the output must be pinned
# to the interpolated value. A formulation that selects the wrong triangle in
# some cell lets the output deviate from it.
function check_triangle_selection(method, triangulation, n)
    d = range(0; stop = 1, length = n)
    f = (a, b) -> a^2 + 3b^2 - a * b + sin(3a)
    pwl = BivariatePWLFunction(d, d, f; triangulation = triangulation)
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    (z,), _ = add_formulation(
        model,
        PiecewiseLinear(pwl; method = method),
        x[1],
        x[2],
    )
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

function check_bivariate_formulation(method, triangulation)
    pwl = BivariatePWLFunction(
        0:2,
        0:2,
        (a, b) -> a^2 - b;
        triangulation = triangulation,
    )
    model = Model()
    @variable(model, x[1:2])
    check_formulation(model, x, pwl, method)
    return
end

function check_bivariate_affine(method)
    d = 0:2
    f = (a, b) -> a^2 - b
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    z = PiecewiseLinear(f, d, d; method = method)(2x[1], x[2] + 1)
    fix(x[1], 0.5)
    fix(x[2], 0.0)
    for sense in (MAX_SENSE, MIN_SENSE)
        @objective(model, sense, z)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(z) ≈ f(1, 1) atol = 1e-5
    end
    return
end

function check_bivariate_direction(method)
    d = 0:2
    f = (a, b) -> a^2 - b
    pwl = BivariatePWLFunction(d, d, f)
    seg = pwl.segments[3]
    z0 = sum(v[1] for v in seg.output_vals) / 3
    for (direction, sense) in
        ((Epigraph(), MIN_SENSE), (Hypograph(), MAX_SENSE))
        model = Model(OPTIMIZER)
        @variable(model, x[1:2])
        (z,), _ = add_formulation(
            model,
            PiecewiseLinear(pwl; method = method, direction = direction),
            x[1],
            x[2],
        )
        fix(x[1], sum(v[1] for v in seg.input_vals) / 3)
        fix(x[2], sum(v[2] for v in seg.input_vals) / 3)
        @objective(model, sense, z)
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        @test value(z) ≈ z0 atol = 1e-5
    end
    return
end

function check_bivariate_multiple_outputs(method)
    d = range(0; stop = 1, length = 3)
    f = (a, b) -> a^2 + 3b^2 - a * b
    g = (a, b) -> sin(3a) - b
    base = BivariatePWLFunction(d, d, f)
    segments = [
        SegmentPointRep{2,2}(
            s.input_vals,
            [(f(v...), g(v...)) for v in s.input_vals],
        ) for s in base.segments
    ]
    pwl = PWLFunction(segments, base.structure)
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    z, _ = add_formulation(
        model,
        PiecewiseLinear(pwl; method = method),
        x[1],
        x[2],
    )
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
    return
end

# A function of a * b on the given triangles
function grid_function(triangles)
    f = (a, b) -> a * b
    segments =
        [SegmentPointRep{2,1}(t, [(f(v...),) for v in t]) for t in triangles]
    return PWLFunction(segments, UnstructuredTriangulation())
end

# A valid triangulation of a 3×2 grid, with triangles spanning two cells
const _SPANNING = [
    [(0.0, 0.0), (1.0, 0.0), (2.0, 1.0)],
    [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
    [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
    [(0.0, 0.0), (1.0, 1.0), (2.0, 1.0)],
]

# The grid methods reject triangulations that do not split each cell of the
# grid into two triangles
function check_grid_geometry(method)
    # The first cell is split along both diagonals
    crossing = [
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0)],
        [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
        [(1.0, 0.0), (1.0, 1.0), (2.0, 1.0)],
    ]
    # The first cell has three triangles and the second only one
    unbalanced = [
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(1.0, 0.0), (2.0, 0.0), (2.0, 1.0)],
    ]
    for triangles in (_SPANNING, crossing, unbalanced)
        model = Model()
        @variable(model, x[1:2])
        check_rejected(model, x, grid_function(triangles); method = method)
    end
    return
end

# The general methods formulate triangulations with triangles spanning several
# cells of the grid
function check_spanning_triangulation(method)
    spanning = grid_function(_SPANNING)
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    (z,), _ =
        add_formulation(model, PiecewiseLinear(spanning; method = method), x...)
    fix(x[1], 1.0)
    fix(x[2], 0.5)
    for sense in (MAX_SENSE, MIN_SENSE)
        @objective(model, sense, z)
        solve!(model)
        @test value(z) ≈ spanning(1.0, 0.5) atol = 1e-5
    end
    return
end

# Functions of three variables

function check_trivariate(method)
    @testset "Kuhn triangulation" begin
        # Kuhn triangulation of the unit cube into 6 tetrahedra
        g = v -> v[1] + 2v[2]^2 - 3v[1] * v[3]
        segments = SegmentPointRep{3,1}[]
        for σ in
            ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
            vertices = [(0.0, 0.0, 0.0)]
            for j in σ
                push!(vertices, ntuple(i -> i == j ? 1.0 : vertices[end][i], 3))
            end
            push!(
                segments,
                SegmentPointRep{3,1}(vertices, [(g(v),) for v in vertices]),
            )
        end
        pwl = PWLFunction(segments, Unstructured{3}())
        model = Model(OPTIMIZER)
        @variable(model, x[1:3])
        (y,), _ = add_formulation(
            model,
            PiecewiseLinear(pwl; method = method),
            x[1],
            x[2],
            x[3],
        )
        for seg in segments
            for i in 1:3
                fix(x[i], sum(v[i] for v in seg.input_vals) / 4)
            end
            y0 = sum(v[1] for v in seg.output_vals) / 4
            for sense in (MAX_SENSE, MIN_SENSE)
                @objective(model, sense, y)
                solve!(model)
                @test termination_status(model) == MOI.OPTIMAL
                @test value(y) ≈ y0 atol = 1e-5
            end
        end
    end
    @testset "Formulation" begin
        s = SegmentPointRep{3,1}(
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (0.0, 1.0, 0.0),
                (0.0, 0.0, 1.0),
            ],
            [(0.0,), (1.0,), (2.0,), (3.0,)],
        )
        model = Model()
        @variable(model, x[1:3])
        pwl = PWLFunction([s], Unstructured{3}())
        check_formulation(model, x, pwl, method)
    end
    return
end

end # module
