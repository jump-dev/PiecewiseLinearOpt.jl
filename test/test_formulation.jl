# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The interface of PiecewiseLinear and add_formulation, and the PWLFormulation

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
    _, f1 = add_formulation(model, PiecewiseLinear(pwl), x)
    _, f2 = add_formulation(model, PiecewiseLinear(pwl), x)
    _, f3 = add_formulation(model, PiecewiseLinear(pwl), x; name = "cost")
    @test f1.name != f2.name
    @test f3.name == "cost"
    @test all(v -> startswith(name(v), "cost_"), f3.variables)
    # Inputs and outputs may be variables, expressions and constants
    @variable(model, z)
    y, f =
        add_formulation(model, PiecewiseLinear(pwl), 2x - 1; output_vars = [z])
    @test y == [z]
    @test !(z in f.variables)
    y, _ = add_formulation(model, PiecewiseLinear(pwl), 1.5)
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
    @variable(model, z)
    f = PiecewiseLinear([0, 1, 2], [0, 2, 1])
    ret = f(x; output_var = y)
    @test ret === y
    # A constant input is formulated when the output is given
    @test f(0.5; output_var = z) === z
    fix(x, 1.5)
    @objective(model, Max, y + z)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y) ≈ 1.5 atol = 1e-5
    @test value(z) ≈ 1.0 atol = 1e-5
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

function test_piecewise_linear()
    f = PiecewiseLinear(t -> abs(t - 1), [0, 1, 2]; direction = Epigraph())
    g = PiecewiseLinear(0:1, 0:1; method = SixStencil(ZigZagBinary())) do a, b
        return a + 2b
    end
    @test f.pwl isa UnivariatePWLFunction
    @test g.method == SixStencil(ZigZagBinary())
    @test sprint(show, g) ==
          "PiecewiseLinear R^2 -> R^1 with 2 segments, SixStencil(ZigZagBinary()) and Graph"
    @test PiecewiseLinear([0, 1, 2], [1, 0, 1]).pwl(0.5) ≈ 0.5
    @test PiecewiseLinear(0:1, 0:1, [0 2; 1 3]).pwl(1.0, 1.0) ≈ 3
    # Constant inputs are evaluated
    @test f(0.5) ≈ 0.5
    @test g(1, 0.5) ≈ 2
    @test f(AffExpr(1.5)) ≈ 0.5
    # Each call adds a formulation to the model of the inputs
    model = Model(OPTIMIZER)
    @variable(model, 0 ≤ x[1:2] ≤ 2)
    @expression(model, fx, f(x[1]; name = "fx"))
    @test fx isa VariableRef
    @test startswith(name(fx), "fx_")
    @test !has_upper_bound(fx)
    n = num_variables(model)
    gx = g(x[1] / 2, 0.5)
    @test num_variables(model) > n
    @objective(model, Min, fx - gx)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    # f(x) - g(x / 2, 0.5) = |x - 1| - x / 2 - 1 is smallest at x = 1
    @test objective_value(model) ≈ -1.5 atol = 1e-5
    @test value(x[1]) ≈ 1 atol = 1e-5
    # Several outputs are returned as a tuple
    segments =
        [SegmentPointRep{1,2}([(0.0,), (1.0,)], [(0.0, 1.0), (2.0, -1.0)])]
    h = PiecewiseLinear(PWLFunction(segments, Intervals()))
    y = h(x[2])
    @test y isa Tuple{VariableRef,VariableRef}
    @test h(0.25) == (0.5, 0.5)
    @test_throws ErrorException h(x[2]; output_var = x[1])
    # The number of inputs must match the function
    @test_throws ErrorException f(x[1], x[2])
    @test_throws ErrorException f(0.5, 0.5)
    # Inputs from different models
    @variable(Model(), z)
    @test_throws ErrorException g(x[1], z)
    return
end

end # module

TestFormulation.runtests()
