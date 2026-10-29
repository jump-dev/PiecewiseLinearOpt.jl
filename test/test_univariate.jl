# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Univariate functions across the formulations

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "utilities.jl"))

module TestUnivariate

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

function test_constructors()
    d = 0:0.01:1
    f = (xi -> xi^2)
    fd = [f(xi) for xi in d]
    pwl = UnivariatePWLFunction(d, f)

    model = Model(OPTIMIZER)
    @variable(model, x)
    y1, _ = add_formulation(model, PiecewiseLinear(pwl), x)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y1[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y1[1]) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(OPTIMIZER)
    @variable(model, x)
    y2 = PiecewiseLinear(f, d)(x)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y2)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y2) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(OPTIMIZER)
    @variable(model, x)
    y3 = PiecewiseLinear(d, fd)(x)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y3)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y3) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4
    return
end

# Several functions in one model: all methods must agree on the optimum
function test_methods_agree()
    supply, demand = [3.0, 4.0], [2.0, 2.5, 2.5]
    xs = 0:4
    cost = [[0.0, 2.0, 3.0, 3.5, 3.8], [0.0, 1.0, 2.5, 3.0, 3.2]]
    objectives = map(METHODS_1D) do method
        model = Model(OPTIMIZER)
        @variable(model, x[1:2, 1:3] ≥ 0)
        @constraint(model, [i in 1:2], sum(x[i, :]) == supply[i])
        @constraint(model, [j in 1:3], sum(x[:, j]) == demand[j])
        @objective(
            model,
            Min,
            sum(
                PiecewiseLinear(xs, cost[i]; method = method)(x[i, j]) for
                i in 1:2, j in 1:3
            )
        )
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
    return
end

# A function given in pieces, with a jump and a gap between them
function test_pieces()
    f = PiecewiseLinear([t -> -t, t -> 1 + t^2], [[-2, 0], 0:0.5:2])
    @test length(f.pwl.segments) == 5
    @test f.pwl(-1.0) ≈ 1
    @test f.pwl(1.5) ≈ 3.25
    # At the jump, evaluation gives the left piece and the model either piece
    @test f.pwl(0.0) ≈ 0
    for (x0, values) in
        ((-1.0, (1.0, 1.0)), (0.0, (0.0, 1.0)), (1.5, (3.25, 3.25)))
        for (sense, value0) in zip((MIN_SENSE, MAX_SENSE), values)
            model = Model(OPTIMIZER)
            @variable(model, x == x0)
            @objective(model, sense, f(x))
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test objective_value(model) ≈ value0 atol = 1e-6
        end
    end
    # A gap between the pieces is left out of the domain
    g = PiecewiseLinear([[0, 1], [2, 3]], [[0, 1], [5, 3]])
    model = Model(OPTIMIZER)
    @variable(model, 0 ≤ x ≤ 3)
    @constraint(model, x ≥ 1.5)
    @objective(model, Min, x + g(x))
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 3 atol = 1e-6
    # Invalid pieces
    @test_throws ErrorException PiecewiseLinear(
        [[0, 2], [1, 3]],
        [[0, 1], [1, 3]],
    )
    @test_throws ErrorException PiecewiseLinear(
        [[0, 1], [1, 0.5]],
        [[0, 1], [1, 3]],
    )
    @test_throws ErrorException PiecewiseLinear([[0]], [[1]])
    @test_throws ErrorException PiecewiseLinear([[0, 1]], [[1, 2, 3]])
    @test_throws ErrorException PiecewiseLinear([[0, 1]], [[1, 2], [3, 4]])
    @test_throws ErrorException PiecewiseLinear([sin], [0:1, 1:2])
    @test_throws ErrorException PiecewiseLinear([log], [[0, 1]])
    return
end

# The default method depends on whether the function is continuous
function test_default_method()
    continuous = PiecewiseLinear([[0, 1], [1, 2]], [[0, 1], [1, 3]])
    @test PiecewiseLinearOpt._method(continuous) == Logarithmic()
    # The order of the segments does not matter
    segments = reverse(continuous.pwl.segments)
    reversed = PiecewiseLinear(PWLFunction(segments, Intervals()))
    @test PiecewiseLinearOpt._method(reversed) == Logarithmic()
    jump = PiecewiseLinear([[0, 1], [1, 2]], [[0, 1], [2, 3]])
    @test PiecewiseLinearOpt._method(jump) == MultipleChoice()
    gap = PiecewiseLinear([[0, 1], [2, 3]], [[0, 1], [1, 3]])
    @test PiecewiseLinearOpt._method(gap) == MultipleChoice()
    relaxation = UnivariatePWLRelaxation(0:2, [0, 1, 0], [1, 1, 2])
    @test PiecewiseLinearOpt._method(PiecewiseLinear(relaxation)) ==
          DisaggregatedLogarithmic()
    # The default methods formulate the functions
    for f in (continuous, jump, gap, PiecewiseLinear(relaxation))
        model = Model(OPTIMIZER)
        @variable(model, x == 0.5)
        @objective(model, Max, f(x))
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
    end
    return
end

end # module

TestUnivariate.runtests()
