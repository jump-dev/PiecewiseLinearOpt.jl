# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the ConvexCombination formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestConvexCombination

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

function test_univariate()
    check_univariate(ConvexCombination())
    return
end

function test_bivariate_simple()
    check_bivariate_simple(ConvexCombination())
    return
end

function test_bivariate_grid()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_bivariate_grid(ConvexCombination(), triangulation)
        end
    end
    return
end

function test_triangle_selection()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_triangle_selection(ConvexCombination(), triangulation, 4)
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(ConvexCombination(), UpperTriangulation())
    return
end

function test_trivariate()
    check_trivariate(ConvexCombination())
    return
end

function test_direction()
    check_bivariate_direction(ConvexCombination())
    return
end

function test_multiple_outputs()
    check_bivariate_multiple_outputs(ConvexCombination())
    return
end

function test_spanning_triangulation()
    check_spanning_triangulation(ConvexCombination())
    return
end

end # module

TestConvexCombination.runtests()
