# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the ConvexPieces formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestConvexPieces

using JuMP
using PiecewiseLinearOpt
using Test

import PiecewiseLinearOpt as PLO

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

const METHODS = (ConvexPieces(:mc), ConvexPieces(:log))

function test_discontinuous()
    for method in METHODS
        @testset "$method" begin
            check_discontinuous(method)
        end
    end
    return
end

function test_pieces()
    sigmoid(t) = 1 / (1 + exp(-t))
    xs = -4:0.5:4
    for method in METHODS
        # The number of binaries is that of the pieces that are convex in the
        # direction, encoded as given
        function binaries(p)
            return p == 1 ? 0 :
                   (method.encoding == :mc ? p : ceil(Int, log2(p)))
        end
        for (f, direction, pieces) in (
            (t -> t^2, Epigraph(), 1),
            (t -> t^2, Hypograph(), 16),
            (t -> -t^2, Hypograph(), 1),
            (t -> t^2, Graph(), 16),
            # Convex for t < 0 and concave for t > 0, where the segments next
            # to 0 have the same slope: one piece for one half and the segment
            # after 0, and one for each of the 7 other segments
            (sigmoid, Epigraph(), 8),
            (sigmoid, Hypograph(), 8),
            (t -> abs(t) + (t > 1) * (t - 1)^2, Epigraph(), 1),
        )
            pwl = UnivariatePWLFunction(xs, f)
            @test length(PLO._convex_pieces(pwl, direction)) == pieces
            model = Model(OPTIMIZER)
            @variable(model, x)
            (y,), formulation = add_formulation(
                model,
                PiecewiseLinear(pwl; method = method, direction = direction),
                x,
            )
            @test count(is_binary, formulation.variables) == binaries(pieces)
            # The output is pinned to the interpolated value in the direction
            senses = Dict(
                Epigraph() => (MIN_SENSE,),
                Hypograph() => (MAX_SENSE,),
                Graph() => (MIN_SENSE, MAX_SENSE),
            )[direction]
            for x0 in (-4.0, -3.3, -0.2, 0.0, 1.7, 4.0), sense in senses
                fix(x, x0)
                @objective(model, sense, y)
                solve!(model)
                @test termination_status(model) == MOI.OPTIMAL
                @test value(y) ≈ pwl(x0) atol = 1e-6
            end
            # The domain is enforced
            fix(x, 4.5)
            solve!(model)
            @test termination_status(model) in
                  (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
        end
    end
    return
end

# Discontinuous functions in the relaxed directions, where the pieces end at
# the jumps and gaps
function test_relaxed_discontinuous()
    xs = [(0, 1), (1, 2), (2, 3), (3.5, 4), (4, 5), (5, 5.5), (6, 7), (7, 8)]
    ys = [(0, 2), (1, 0), (0, 3), (3, -1), (2, 2), (2, 4), (4, 1), (1, 0)]
    segments = [
        SegmentPointRep{1,1}(
            [(float(a),), (float(b),)],
            [(float(c),), (float(e),)],
        ) for ((a, b), (c, e)) in zip(xs, ys)
    ]
    pwl = PWLFunction(segments, Intervals())
    for method in METHODS,
        (direction, sense, best) in
        ((Epigraph(), MIN_SENSE, min), (Hypograph(), MAX_SENSE, max))

        model = Model(OPTIMIZER)
        @variable(model, x)
        (y,), _ = add_formulation(
            model,
            PiecewiseLinear(pwl; method = method, direction = direction),
            x,
        )
        @objective(model, sense, y)
        for ((a, b), (c, e)) in zip(xs, ys)
            fix(x, 0.3 * a + 0.7 * b)
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ 0.3 * c + 0.7 * e atol = 1e-5
        end
        # At a jump, the best of the two values
        fix(x, 1.0)
        solve!(model)
        @test value(y) ≈ best(2.0, 1.0) atol = 1e-5
        fix(x, 3.2)
        solve!(model)
        @test termination_status(model) in
              (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
    end
    return
end

# Same optimum as the other methods on concave costs
function test_concave_costs()
    supply, demand = [3.0, 4.0], [2.0, 2.5, 2.5]
    cost = [[0.0, 2.0, 3.0, 3.5, 3.8], [0.0, 1.0, 2.5, 3.0, 3.2]]
    objectives = map((METHODS..., Logarithmic())) do m
        model = Model(OPTIMIZER)
        @variable(model, x[1:2, 1:3] ≥ 0)
        @constraint(model, [i in 1:2], sum(x[i, :]) == supply[i])
        @constraint(model, [j in 1:3], sum(x[:, j]) == demand[j])
        @objective(
            model,
            Min,
            sum(
                PiecewiseLinear(
                    0:4,
                    cost[i];
                    method = m,
                    direction = Epigraph(),
                )(
                    x[i, j],
                ) for i in 1:2, j in 1:3
            )
        )
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[end]; atol = 1e-5))
    return
end

function test_unsupported_functions()
    model = Model()
    @variable(model, x)
    vertical = PWLFunction(
        [SegmentPointRep{1,1}([(1.0,), (1.0,)], [(0.0,), (1.0,)])],
        Intervals(),
    )
    two_outputs = PWLFunction(
        [SegmentPointRep{1,2}([(0.0,), (1.0,)], [(0.0, 1.0), (1.0, 0.0)])],
        Intervals(),
    )
    for method in METHODS, pwl in (vertical, two_outputs)
        check_rejected(model, [x], pwl; method = method)
    end
    return
end

function test_options()
    @test ConvexPieces() == ConvexPieces(:mc)
    @test_throws ErrorException ConvexPieces(:unknown)
    @test sprint(show, ConvexPieces(:log)) == "ConvexPieces(:log)"
    return
end

end # module

TestConvexPieces.runtests()
