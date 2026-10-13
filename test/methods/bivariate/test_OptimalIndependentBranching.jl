# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the OptimalIndependentBranching formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestOptimalIndependentBranching

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

function test_bivariate_simple()
    check_bivariate_simple(OptimalIndependentBranching(OPTIMIZER))
    return
end

# The cover computed for OptimalIndependentBranching gets expensive quickly
function test_triangle_selection()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_triangle_selection(
                OptimalIndependentBranching(OPTIMIZER),
                triangulation,
                3,
            )
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(
        OptimalIndependentBranching(OPTIMIZER),
        UpperTriangulation(),
    )
    return
end

function test_grid_geometry()
    check_grid_geometry(OptimalIndependentBranching(OPTIMIZER))
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
    # Throws an error, also when the result is cached
    model = Model()
    @variable(model, x[1:2])
    for _ in 1:2
        check_rejected(
            model,
            x,
            pwl;
            method = OptimalIndependentBranching(OPTIMIZER; time_limit = 0),
        )
    end
    @test PLO._OPTIMAL_COVER_CACHE[(nw_se, true)].cover === nothing
    return
end

end # module

TestOptimalIndependentBranching.runtests()
