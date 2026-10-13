# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The interface of piecewiselinear and the PWLFormulation it returns

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "utilities.jl"))

module TestFormulation

using JuMP
using PiecewiseLinearOpt
using Test

using ..PLOTestUtilities

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_names_and_inputs()
    pwl = UnivariatePWLFunction([0, 1, 2], [0, 2, 1])
    model = Model(OPTIMIZER)
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
    check_rejected(model, [x, x], pwl)
    check_rejected(model, [x], pwl; output_vars = [z, z])
    # Deleting from another model is an error
    @test_throws ErrorException delete(Model(), f)
    return
end

function test_output_variable()
    model = Model(OPTIMIZER)
    @variable(model, x)
    @variable(model, y)
    ret = piecewiselinear(model, x, [0, 1, 2], [0, 2, 1]; output_var = y)
    @test ret === y
    fix(x, 1.5)
    @objective(model, Max, y)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y) ≈ 1.5 atol = 1e-5
    return
end

function test_errors()
    model = Model()
    @variable(model, x[1:3])
    # No segments
    pwl = PWLFunction(SegmentPointRep{1,1}[], Intervals())
    check_rejected(model, (x[1],), pwl)
    # Discontinuous functions require MultipleChoice
    s1 = SegmentPointRep{1,1}([(0.0,), (1.0,)], [(0.0,), (1.0,)])
    s2 = SegmentPointRep{1,1}([(1.0,), (2.0,)], [(3.0,), (2.0,)])
    pwl = PWLFunction([s1, s2], Intervals())
    for method in (ConvexCombination(), Logarithmic(), Incremental())
        check_rejected(model, (x[1],), pwl; method = method)
    end
    # Grid based methods require a connected domain
    s2 = SegmentPointRep{1,1}([(2.0,), (3.0,)], [(3.0,), (2.0,)])
    pwl = PWLFunction([s1, s2], Intervals())
    check_rejected(model, (x[1],), pwl; method = Logarithmic())
    # Grid based methods require segments between neighboring breakpoints,
    # also when the number of segments is right
    overlapping = PWLFunction(
        [
            SegmentPointRep{1,1}([(0.0,), (2.0,)], [(0.0,), (2.0,)]),
            SegmentPointRep{1,1}([(1.0,), (2.0,)], [(1.0,), (2.0,)]),
        ],
        Intervals(),
    )
    for method in (Logarithmic(), Incremental(), ZigZagBinary(), NativeSOS2())
        check_rejected(model, (x[1],), overlapping; method = method)
    end
    # Unsupported combinations of dimension and method
    pwl = BivariatePWLFunction(0:1, 0:1, (a, b) -> a + b)
    check_rejected(model, (x[1], x[2]), pwl; method = Incremental())
    s = SegmentPointRep{3,1}(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (3.0,)],
    )
    pwl = PWLFunction([s], Unstructured{3}())
    check_rejected(model, (x[1], x[2], x[3]), pwl; method = SixStencil())
    return
end

end # module

TestFormulation.runtests()
