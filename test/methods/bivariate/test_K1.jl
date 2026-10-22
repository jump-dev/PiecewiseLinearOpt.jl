# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the K1 formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestK1

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

function test_bivariate_grid()
    for axis_method in AXIS_METHODS
        @testset "$axis_method" begin
            check_bivariate_grid(K1(axis_method), K1Triangulation())
        end
    end
    return
end

function test_triangle_selection()
    for n in (4, 5)
        @testset "$n points" begin
            check_triangle_selection(K1(), K1Triangulation(), n)
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(K1(), K1Triangulation())
    return
end

function test_triangulation_checks()
    d = range(0; stop = 1, length = 4)
    # The triangulation itself is checked, not how it was built
    for triangulation in (UnionJackTriangulation(), UpperTriangulation())
        pwl = BivariatePWLFunction(
            d,
            d,
            (a, b) -> a * b;
            triangulation = triangulation,
        )
        model = Model()
        @variable(model, x[1:2])
        check_rejected(model, (x[1], x[2]), pwl; method = K1())
    end
    # For a * b, the lower triangulation splits every cell like K1
    pwl = BivariatePWLFunction(
        d,
        d,
        (a, b) -> a * b;
        triangulation = LowerTriangulation(),
    )
    model = Model()
    @variable(model, x[1:2])
    _, formulation = piecewiselinear(model, x, pwl; method = K1())
    @test formulation.method == K1()
    return
end

function test_affine_inputs()
    check_bivariate_affine(K1())
    return
end

function test_multiple_outputs()
    check_bivariate_multiple_outputs(K1())
    return
end

end # module

TestK1.runtests()
