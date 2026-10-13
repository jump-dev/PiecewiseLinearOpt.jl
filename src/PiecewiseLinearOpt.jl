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

include("methods/univariate/convex_pieces.jl")
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

include("methods/bivariate/greedy_triangle_selection.jl")
include("methods/bivariate/k1.jl")
include("methods/bivariate/nine_stencil.jl")
include("methods/bivariate/optimal_independent_branching.jl")
include("methods/bivariate/optimal_triangle_selection.jl")
include("methods/bivariate/six_stencil.jl")
include("methods/bivariate/union_jack.jl")
include("methods/bivariate/common.jl")

include("methods/multivariate/convex_combination.jl")
include("methods/multivariate/disaggregated_logarithmic.jl")
include("methods/multivariate/multiple_choice.jl")

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
