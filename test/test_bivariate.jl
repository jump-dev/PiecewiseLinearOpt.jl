# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Bivariate functions across the formulations

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "utilities.jl"))

module TestBivariate

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
    d = 0:0.05:1
    f = (xi, yi) -> xi^2 + yi^2
    pwl = BivariatePWLFunction(d, d, f)

    model = Model(OPTIMIZER)
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

    model = Model(OPTIMIZER)
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
    return
end

# Several functions in one model: all methods must agree on the optimum
function test_methods_agree()
    d = range(0; stop = 1, length = 3)
    f = (a, b) -> (a - 0.4)^2 - sin(3b) + a * b
    methods = [
        ConvexCombination(),
        DisaggregatedLogarithmic(),
        MultipleChoice(),
        GreedyTriangleSelection(),
        NineStencil(),
        OptimalTriangleSelection(OPTIMIZER),
        SixStencil(),
        OptimalIndependentBranching(OPTIMIZER),
        K1(),
    ]
    objectives = map(methods) do method
        model = Model(OPTIMIZER)
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
    return
end

end # module

TestBivariate.runtests()
