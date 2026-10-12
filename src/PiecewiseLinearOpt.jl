# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module PiecewiseLinearOpt

using JuMP

import LinearAlgebra
import Random

export PWLFunction,
    UnivariatePWLFunction, BivariatePWLFunction, PWLFormulation, piecewiselinear

export SegmentPointRep, SegmentStructure, Intervals, Unstructured

# Graph is the default direction, and would collide with Graphs.jl
export Epigraph, Hypograph

export GridTriangulation,
    BestFitTriangulation,
    K1Triangulation,
    LowerTriangulation,
    RandomTriangulation,
    UnionJackTriangulation,
    UnstructuredTriangulation,
    UpperTriangulation

include("types.jl")

# Counts the piecewise linear functions added to a model, to give their
# variables and constraints distinct default names
mutable struct PWLData
    counter::Int
    PWLData(counter::Int = 0) = new(counter)
end

function JuMP.copy_extension_data(
    data::PWLData,
    ::JuMP.AbstractModel,
    ::JuMP.AbstractModel,
)
    return PWLData(data.counter)
end

function _next_default_name(model::JuMP.GenericModel)
    data = get!(model.ext, :PWL, PWLData())::PWLData
    data.counter += 1
    return "pwl$(data.counter)"
end

"""
    PiecewiseLinearOpt.Scalar

The type of the inputs and outputs of a function in a formulation: variables,
expressions or constants.
"""
const Scalar = Union{Real,JuMP.AbstractJuMPScalar}

include("formulation.jl")
include("methods/util.jl")

export Incremental,
    IndependentBranching,
    LogarithmicEmbedding,
    LogarithmicIndependentBranching,
    NativeSOS2,
    ZigZagBinary,
    ZigZagInteger
include("methods/univariate/incremental.jl")
include("methods/univariate/independent_branching.jl")
include("methods/univariate/logarithmic_embedding.jl")
include("methods/univariate/logarithmic_independent_branching.jl")
include("methods/univariate/native_sos2.jl")
include("methods/univariate/zig_zag_binary.jl")
include("methods/univariate/zig_zag_integer.jl")
# ConvexCombination has an SOS2 formulation, so defer this until after the
# multivariate formulations are defined
include("methods/univariate/sos2_formulation_base.jl")

# Consider the colloqial "log" to refer to the embedding formulation
const Logarithmic = LogarithmicEmbedding
export Logarithmic

export GreedyTriangleSelection,
    K1,
    NineStencil,
    OptimalIndependentBranching,
    OptimalTriangleSelection,
    SixStencil,
    UnionJack
include("methods/bivariate/greedy_triangle_selection.jl")
include("methods/bivariate/k1.jl")
include("methods/bivariate/nine_stencil.jl")
include("methods/bivariate/optimal_independent_branching.jl")
include("methods/bivariate/optimal_triangle_selection.jl")
include("methods/bivariate/six_stencil.jl")
include("methods/bivariate/union_jack.jl")
include("methods/bivariate/common.jl")

export ConvexCombination, DisaggregatedLogarithmic, MultipleChoice
include("methods/multivariate/convex_combination.jl")
include("methods/multivariate/disaggregated_logarithmic.jl")
include("methods/multivariate/multiple_choice.jl")

include("pwlinear.jl")

# The extension interface is public but not exported
@static if VERSION >= v"1.11.0-DEV.469"
    eval(
        Meta.parse(
            """
            public Method,
                SOS2Method,
                TriangleSelectionMethod,
                Direction,
                Graph,
                Scalar,
                formulate_pwl!,
                formulate_sos2!,
                formulate_triangle_selection!,
                axis_method,
                record!,
                variable_name,
                constrain_output!,
                add_biclique_cover!,
                cell_diagonals,
                segments,
                input_values,
                output_values
            """,
        ),
    )
end

end # module
