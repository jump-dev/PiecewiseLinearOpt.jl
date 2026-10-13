# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the SixStencil formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestSixStencil

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

method(axis_method) = SixStencil(axis_method)

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
            check_triangle_selection(SixStencil(), triangulation, 4)
        end
    end
    return
end

function test_formulation()
    for axis_method in AXIS_METHODS
        @testset "$axis_method" begin
            check_bivariate_formulation(
                method(axis_method),
                UpperTriangulation(),
            )
        end
    end
    return
end

function test_grid_geometry()
    check_grid_geometry(SixStencil())
    return
end

function test_affine_inputs()
    check_bivariate_affine(SixStencil())
    return
end

function test_direction()
    check_bivariate_direction(SixStencil())
    return
end

function test_multiple_outputs()
    check_bivariate_multiple_outputs(SixStencil())
    return
end

end # module

TestSixStencil.runtests()
