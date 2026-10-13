# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the UnionJack formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestUnionJack

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
            check_bivariate_grid(
                UnionJack(axis_method),
                UnionJackTriangulation(),
            )
        end
    end
    return
end

function test_triangle_selection()
    for n in (4, 5)
        @testset "$n points" begin
            check_triangle_selection(UnionJack(), UnionJackTriangulation(), n)
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(UnionJack(), UnionJackTriangulation())
    return
end

function test_triangulation_checks()
    # The triangulation itself is checked, not how it was built
    d = range(0; stop = 1, length = 4)
    pwl = BivariatePWLFunction(
        d,
        d,
        (a, b) -> a * b;
        triangulation = K1Triangulation(),
    )
    model = Model()
    @variable(model, x[1:2])
    check_rejected(model, (x[1], x[2]), pwl; method = UnionJack())
    return
end

end # module

TestUnionJack.runtests()
