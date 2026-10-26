# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the Incremental formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestIncremental

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
    check_univariate(Incremental())
    return
end

function test_size()
    for d in 1:5
        model = Model()
        @variable(model, x)
        pwl = UnivariatePWLFunction(0:d, t -> t^2)
        _, formulation = add_formulation(
            model,
            PiecewiseLinear(pwl; method = Incremental()),
            x,
        )
        # d fill variables, d - 1 binary variables and the output
        @test count(is_binary, formulation.variables) == d - 1
        @test length(formulation.variables) == 2d
    end
    return
end

end # module

TestIncremental.runtests()
