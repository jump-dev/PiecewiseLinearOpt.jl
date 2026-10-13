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

function test_convenience_methods()
    d = 0:0.01:1
    f = (xi -> xi^2)
    fd = [f(xi) for xi in d]
    pwl = UnivariatePWLFunction(d, f)

    model = Model(OPTIMIZER)
    @variable(model, x)
    y1, _ = piecewiselinear(model, (x,), pwl)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y1[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y1[1]) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(OPTIMIZER)
    @variable(model, x)
    y2 = piecewiselinear(model, x, d, f)
    @constraint(model, x ≤ 0.75)
    @objective(model, Max, y2)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(y2) ≈ f(value(x)) rtol = 1e-4
    @test objective_value(model) ≈ 0.5625 rtol = 1e-4

    model = Model(OPTIMIZER)
    @variable(model, x)
    y3 = piecewiselinear(model, x, d, fd)
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
                piecewiselinear(model, x[i, j], xs, cost[i]; method = method) for i in 1:2, j in 1:3
            )
        )
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
    return
end

end # module

TestUnivariate.runtests()
