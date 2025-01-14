# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module PiecewiseLinearOpt

using JuMP

import LinearAlgebra
import Random

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear

include("types.jl")

mutable struct PWLData
    counter::Int
    PWLData() = new(0)
end

function initPWL!(m::JuMP.Model)
    if !haskey(m.ext, :PWL)
        m.ext[:PWL] = PWLData()
    end
    return nothing
end

const VarOrAff = Union{JuMP.VariableRef,JuMP.AffExpr}

include("methods/util.jl")

export Incremental,
    LogarithmicEmbedding,
    LogarithmicIndependentBranching,
    NativeSOS2,
    ZigZagBinary,
    ZigZagInteger
include("methods/univariate/incremental.jl")

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

export K1,
    NineStencil,
    OptimalIndependentBranching,
    OptimalTriangleSelection,
    SixStencil,
    UnionJack
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

end # module
