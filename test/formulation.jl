# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The interface of piecewiselinear and the PWLFormulation it returns

# The formulation must record exactly what was added to a model that only
# contained the inputs, and deleting it must remove all of it.
function test_formulation(model, inputs, pwl, method)
    y, formulation = piecewiselinear(model, inputs, pwl; method = method)
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

@testset "Formulation: univariate $method" for method in methods_1D
    pwl = UnivariatePWLFunction([0, 1, 2, 4, 5, 7], [0, 2, 1, 3, -1, 0])
    model = Model()
    @variable(model, x)
    test_formulation(model, [x], pwl, method)
end

@testset "Formulation: discontinuous $method" for method in
                                                  methods_discontinuous
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    s3 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(2.0,), (5.0,)])
    model = Model()
    @variable(model, x)
    test_formulation(
        model,
        [x],
        PLO.PWLFunction([s1, s2, s3], PLO.Intervals()),
        method,
    )
end

@testset "Formulation: bivariate $method" for (method, triangulation) in [
    [(m, UpperTriangulation()) for m in methods_2D_selection]
    [(SixStencil(m), UpperTriangulation()) for m in axis_methods]
    (OptimalIndependentBranching(optimizer), UpperTriangulation())
    (K1(), K1Triangulation())
    (UnionJack(), UnionJackTriangulation())
]
    pwl = BivariatePWLFunction(
        0:2,
        0:2,
        (a, b) -> a^2 - b;
        triangulation = triangulation,
    )
    model = Model()
    @variable(model, x[1:2])
    test_formulation(model, x, pwl, method)
end

@testset "Formulation: trivariate $method" for method in [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
]
    s = PLO.SegmentPointRep{3,1}(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (3.0,)],
    )
    model = Model()
    @variable(model, x[1:3])
    test_formulation(
        model,
        x,
        PLO.PWLFunction([s], PLO.Unstructured{3}()),
        method,
    )
end

@testset "Formulation: names and inputs" begin
    pwl = UnivariatePWLFunction([0, 1, 2], [0, 2, 1])
    model = Model(optimizer)
    @variable(model, x)
    # Default names are distinct, and can be chosen
    _, f1 = piecewiselinear(model, [x], pwl)
    _, f2 = piecewiselinear(model, [x], pwl)
    _, f3 = piecewiselinear(model, [x], pwl; name = "cost")
    @test f1.name != f2.name
    @test f3.name == "cost"
    @test all(v -> startswith(name(v), "cost_"), f3.variables)
    # Inputs and outputs may be vectors or tuples of variables, expressions
    # and constants
    @variable(model, z)
    y, f = piecewiselinear(model, [2x - 1], pwl; output_vars = [z])
    @test y == [z]
    @test !(z in f.variables)
    y, _ = piecewiselinear(model, (1.5,), pwl)
    fix(x, 0.75)
    @objective(model, Max, z + y[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(z) ≈ 1.0 atol = 1e-5  # f(2 * 0.75 - 1)
    @test value(y[1]) ≈ 1.5 atol = 1e-5  # f(1.5)
    # The number of inputs and outputs must match the function
    test_rejected(model, [x, x], pwl)
    test_rejected(model, [x], pwl; output_vars = [z, z])
    # Deleting from another model is an error
    @test_throws ErrorException delete(Model(), f)
end

@testset "Output variable" begin
    model = Model(optimizer)
    @variable(model, x)
    @variable(model, y)
    ret = piecewiselinear(model, x, [0, 1, 2], [0, 2, 1]; output_var = y)
    @test ret === y
    fix(x, 1.5)
    @objective(model, Max, y)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y) ≈ 1.5 atol = 1e-5
end

@testset "Errors" begin
    model = Model()
    @variable(model, x[1:3])
    # No segments
    pwl = PLO.PWLFunction(PLO.SegmentPointRep{1,1}[], PLO.Intervals())
    test_rejected(model, (x[1],), pwl)
    # Discontinuous functions require MultipleChoice
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    for method in (ConvexCombination(), Logarithmic(), Incremental())
        test_rejected(model, (x[1],), pwl; method = method)
    end
    # Grid based methods require a connected domain
    s2 = PLO.SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2], PLO.Intervals())
    test_rejected(model, (x[1],), pwl; method = Logarithmic())
    # Grid based methods require segments between neighboring breakpoints,
    # also when the number of segments is right
    overlapping = PLO.PWLFunction(
        [
            PLO.SegmentPointRep{1,1}([(0.0,), (2.0,)], [(0.0,), (2.0,)]),
            PLO.SegmentPointRep{1,1}([(1.0,), (2.0,)], [(1.0,), (2.0,)]),
        ],
        PLO.Intervals(),
    )
    for method in (Logarithmic(), Incremental(), ZigZagBinary(), NativeSOS2())
        test_rejected(model, (x[1],), overlapping; method = method)
    end
    # Unsupported combinations of dimension and method
    pwl = BivariatePWLFunction(0:1, 0:1, (a, b) -> a + b)
    test_rejected(model, (x[1], x[2]), pwl; method = Incremental())
    s = PLO.SegmentPointRep{3,1}(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (3.0,)],
    )
    pwl = PLO.PWLFunction([s], PLO.Unstructured{3}())
    test_rejected(model, (x[1], x[2], x[3]), pwl; method = SixStencil())
end
