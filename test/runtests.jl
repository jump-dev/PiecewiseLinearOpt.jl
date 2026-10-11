# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The standard run tests every method, but combines the bivariate grid methods
# only with a few representative axis methods. The full run uses every
# combination; start it with
#
#     Pkg.test("PiecewiseLinearOpt"; test_args = ["full"])
#
# or by setting the environment variable PLO_TEST_FULL=true.

using PiecewiseLinearOpt
using HiGHS
using JuMP
using Test

import Aqua

const PLO = PiecewiseLinearOpt

const FULL = "full" in ARGS || get(ENV, "PLO_TEST_FULL", "false") == "true"

optimizer = optimizer_with_attributes(HiGHS.Optimizer, MOI.Silent() => true)

# HiGHS occasionally reports a solve error on these small models, when the
# postsolved solution violates a constraint by exactly the MIP feasibility
# tolerance, whatever its value. Solving again without presolve avoids it.
# Presolve stays off for the model, as changing it back would discard the
# solution.
function solve!(model)
    optimize!(model)
    if termination_status(model) == MOI.OTHER_ERROR
        set_attribute(model, "presolve", "off")
        optimize!(model)
    end
    return
end

# All constraints except bounds and integrality restrictions of variables
function nonbound_constraints(model)
    bounds = (
        MOI.GreaterThan{Float64},
        MOI.LessThan{Float64},
        MOI.EqualTo{Float64},
        MOI.ZeroOne,
        MOI.Integer,
    )
    constraints = ConstraintRef[]
    for (F, S) in list_of_constraint_types(model)
        if !(F == VariableRef && S in bounds)
            append!(constraints, all_constraints(model, F, S))
        end
    end
    return constraints
end

# Tests that piecewiselinear(model, args...; kwargs...) throws an error and
# leaves the model unchanged
function test_rejected(model, args...; kwargs...)
    variables = all_variables(model)
    constraints = nonbound_constraints(model)
    @test_throws ErrorException piecewiselinear(model, args...; kwargs...)
    @test all_variables(model) == variables
    @test nonbound_constraints(model) == constraints
    return
end

const methods_1D = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    Incremental(),
    IndependentBranching(:biclique),
    IndependentBranching(:balanced),
    IndependentBranching(:brgc),
    LogarithmicEmbedding(),
    LogarithmicIndependentBranching(),
    MultipleChoice(),
    NativeSOS2(),
    ZigZagBinary(),
    ZigZagInteger(),
]

const methods_discontinuous = [
    DisaggregatedLogarithmic(),
    IndependentBranching(:biclique),
    IndependentBranching(:balanced),
    IndependentBranching(:brgc),
    MultipleChoice(),
]

# The univariate methods used to select a cell in the bivariate grid methods.
# The representative ones cover an embedding formulation, an independent
# branching formulation, a native SOS2 constraint and the generic SOS2
# formulation through a dummy function.
const axis_methods =
    FULL ? methods_1D :
    [Logarithmic(), IndependentBranching(), NativeSOS2(), Incremental()]

# The bivariate methods that work with any grid triangulation
const methods_2D_gen = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
    [NineStencil(m) for m in axis_methods]...,
    [OptimalTriangleSelection(optimizer, m) for m in axis_methods]...,
    [SixStencil(m) for m in axis_methods]...,
]

# One bivariate method of each kind, for tests of the triangle selection
const methods_2D_selection = [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
    NineStencil(),
    OptimalTriangleSelection(optimizer),
    SixStencil(),
]

const triangulations = [
    UpperTriangulation(),
    LowerTriangulation(),
    BestFitTriangulation(),
    K1Triangulation(),
    UnionJackTriangulation(),
    RandomTriangulation(),
]

@testset verbose = true "PiecewiseLinearOpt" begin
    @testset "Aqua" begin
        Aqua.test_all(PiecewiseLinearOpt)
    end
    include("functions.jl")
    include("covers.jl")
    include("univariate.jl")
    include("bivariate.jl")
    include("multivariate.jl")
    include("formulation.jl")
    include("extending.jl")
end
