# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module PiecewiseLinearOpt

using JuMP

import LinearAlgebra
import Random

include("types.jl")

# Counts the piecewise linear functions added to a model, to give their
# variables and constraints distinct default names
mutable struct _PWLData
    counter::Int
    _PWLData(counter::Int = 0) = new(counter)
end

function JuMP.copy_extension_data(
    data::_PWLData,
    ::JuMP.AbstractModel,
    ::JuMP.AbstractModel,
)
    return _PWLData(data.counter)
end

function _next_default_name(model::JuMP.GenericModel)
    data = get!(model.ext, :PWL, _PWLData())::_PWLData
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

include("methods/univariate/ConvexPieces.jl")
include("methods/univariate/Incremental.jl")
include("methods/univariate/IndependentBranching.jl")
include("methods/univariate/LogarithmicEmbedding.jl")
include("methods/univariate/LogarithmicIndependentBranching.jl")
include("methods/univariate/NativeSOS2.jl")
include("methods/univariate/ZigZagBinary.jl")
include("methods/univariate/ZigZagInteger.jl")
# ConvexCombination has an SOS2 formulation, so defer this until after the
# multivariate formulations are defined
include("methods/univariate/sos2_formulation_base.jl")

# Consider the colloqial "log" to refer to the embedding formulation
const Logarithmic = LogarithmicEmbedding

include("methods/bivariate/GreedyTriangleSelection.jl")
include("methods/bivariate/K1.jl")
include("methods/bivariate/NineStencil.jl")
include("methods/bivariate/OptimalIndependentBranching.jl")
include("methods/bivariate/OptimalTriangleSelection.jl")
include("methods/bivariate/SixStencil.jl")
include("methods/bivariate/UnionJack.jl")
include("methods/bivariate/common.jl")

include("methods/multivariate/ConvexCombination.jl")
include("methods/multivariate/DisaggregatedLogarithmic.jl")
include("methods/multivariate/MultipleChoice.jl")

include("pwlinear.jl")

# Export every symbol that does not begin with an underscore, following JuMP
const _EXCLUDE_SYMBOLS = [Symbol(@__MODULE__), :eval, :include]

for sym in names(@__MODULE__; all = true)
    sym_string = string(sym)
    if sym in _EXCLUDE_SYMBOLS ||
       startswith(sym_string, "_") ||
       startswith(sym_string, "@_")
        continue
    end
    if !(
        Base.isidentifier(sym) ||
        (startswith(sym_string, "@") && Base.isidentifier(sym_string[2:end]))
    )
        continue
    end
    @eval export $sym
end

end # module
