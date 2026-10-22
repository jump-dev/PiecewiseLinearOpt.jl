# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the DisaggregatedLogarithmic formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestDisaggregatedLogarithmic

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
    check_univariate(DisaggregatedLogarithmic())
    return
end

function test_relaxation()
    check_relaxation(DisaggregatedLogarithmic())
    return
end

function test_bivariate_simple()
    check_bivariate_simple(DisaggregatedLogarithmic())
    return
end

function test_bivariate_grid()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_bivariate_grid(DisaggregatedLogarithmic(), triangulation)
        end
    end
    return
end

function test_triangle_selection()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_triangle_selection(
                DisaggregatedLogarithmic(),
                triangulation,
                4,
            )
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(
        DisaggregatedLogarithmic(),
        UpperTriangulation(),
    )
    return
end

function test_trivariate()
    check_trivariate(DisaggregatedLogarithmic())
    return
end

function test_discontinuous()
    check_discontinuous(DisaggregatedLogarithmic())
    return
end

function test_multiple_outputs()
    check_bivariate_multiple_outputs(DisaggregatedLogarithmic())
    return
end

# The default method for functions of three or more variables
function test_default()
    s = SegmentPointRep{3,1}(
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (3.0,)],
    )
    model = Model()
    @variable(model, x[1:3])
    _, formulation =
        piecewiselinear(model, x, PWLFunction([s], Unstructured{3}()))
    @test formulation.method == DisaggregatedLogarithmic()
    return
end

end # module

TestDisaggregatedLogarithmic.runtests()
