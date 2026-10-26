# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the OptimalTriangleSelection formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestOptimalTriangleSelection

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

method(axis_method) = OptimalTriangleSelection(OPTIMIZER, axis_method)

function test_bivariate_simple()
    for axis_method in AXIS_METHODS
        @testset "$axis_method" begin
            check_bivariate_simple(method(axis_method))
        end
    end
    return
end

function test_bivariate_grid()
    for axis_method in AXIS_METHODS, triangulation in TRIANGULATIONS
        @testset "$axis_method, $triangulation" begin
            check_bivariate_grid(method(axis_method), triangulation)
        end
    end
    return
end

function test_triangle_selection()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_triangle_selection(
                OptimalTriangleSelection(OPTIMIZER),
                triangulation,
                4,
            )
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(
        OptimalTriangleSelection(OPTIMIZER),
        UpperTriangulation(),
    )
    return
end

function test_grid_geometry()
    check_grid_geometry(OptimalTriangleSelection(OPTIMIZER))
    return
end

function test_direction()
    check_bivariate_direction(OptimalTriangleSelection(OPTIMIZER))
    return
end

function test_multiple_outputs()
    check_bivariate_multiple_outputs(OptimalTriangleSelection(OPTIMIZER))
    return
end

function test_time_limit()
    empty!(PLO._OPTIMAL_COVER_CACHE)
    pwl = BivariatePWLFunction(
        0:5,
        0:5,
        (a, b) -> a * b;
        triangulation = RandomTriangulation(12345),
    )
    nw_se =
        PLO._triangle_directions(pwl, PLO._continuous_gridpoints_or_die(pwl))
    greedy = PLO._greedy_biclique_cover(nw_se, 200)
    @test length(greedy) > 1
    # Falls back to the greedy cover, in addition to 3 binaries for each axis
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    method = OptimalTriangleSelection(OPTIMIZER; time_limit = 0)
    (z,), formulation = @test_logs (:warn,) add_formulation(
        model,
        PiecewiseLinear(pwl; method = method),
        x...,
    )
    @test count(is_binary, formulation.variables) == 6 + length(greedy)
    fix(x[1], 2.25)
    fix(x[2], 3.5)
    @objective(model, Max, z)
    solve!(model)
    @test value(z) ≈ pwl(2.25, 3.5) atol = 1e-5
    # The fallback is cached, and only computed again with a longer time limit
    @test_logs add_formulation(
        model,
        PiecewiseLinear(pwl; method = method),
        x...,
    )
    @test !PLO._OPTIMAL_COVER_CACHE[(nw_se, false)].optimal
    add_formulation(
        model,
        PiecewiseLinear(
            pwl;
            method = OptimalTriangleSelection(OPTIMIZER; time_limit = 60),
        ),
        x...,
    )
    @test PLO._OPTIMAL_COVER_CACHE[(nw_se, false)].optimal
    return
end

end # module

TestOptimalTriangleSelection.runtests()
