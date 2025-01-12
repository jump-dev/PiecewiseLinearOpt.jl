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

function formulate_pwl!(
    model::JuMP.Model,
    input_vals::Vector{NTuple{D,VarOrAff}},
    output_vals::Vector{NTuple{F,VarOrAff}},
    pwl::PWLFunction,
    method::Method,
    direction::DIRECTION,
) where {D,F}
    return error(
        "No support for a R^$D -> R^$F piecewise linear function using the $method method.",
    )
end

_default_method(::Val{1}) = Logarithmic()
_default_method(::Val{2}) = SixStencil()
# _default_method(::Val) = MultipleChoice()

function piecewiselinear(
    model::JuMP.Model,
    input_vars::NTuple{D,VarOrAff},
    pwl::PWLFunction{D,F,SegmentPointRep{D,F}};
    method::Method = _default_method(Val(D)),
    direction::DIRECTION = Graph,
    output_vars::Union{Nothing,NTuple{F,VarOrAff}} = nothing,
) where {D,F}
    initPWL!(model)
    counter = model.ext[:PWL].counter
    counter += 1
    model.ext[:PWL].counter = counter

    if isempty(pwl.segments)
        error(
            "I don't know how to handle a piecewise linear function with no breakpoints.",
        )
    end

    output_lb =
        minimum(minimum(segment.output_vals) for segment in pwl.segments)
    output_ub =
        maximum(maximum(segment.output_vals) for segment in pwl.segments)

    if output_vars === nothing
        output_vars = tuple(
            JuMP.@variable(
                model,
                [i in 1:F],
                lower_bound = output_lb[i],
                upper_bound = output_ub[i],
                base_name = "y_$counter"
            )...,
        )
    end

    formulate_pwl!(model, input_vars, output_vars, pwl, method, direction)
    return output_vars
end

function piecewiselinear(
    model::JuMP.Model,
    input_var::VarOrAff,
    pwl::PWLFunction{1,1,SegmentPointRep{1,1}};
    method::Method = _default_method(Val(1)),
    direction::DIRECTION = Graph,
    output_var::Union{Nothing,VarOrAff} = nothing,
)
    return piecewiselinear(
        model,
        (input_var,),
        pwl;
        method = method,
        direction = direction,
        output_vars = isnothing(output_var) ? nothing : (output_var,),
    )[1]
end

function piecewiselinear(
    model::JuMP.Model,
    input_var_x::VarOrAff,
    input_var_y::VarOrAff,
    pwl::PWLFunction{2,1,SegmentPointRep{2,1}};
    method::Method = _default_method(Val(2)),
    direction::DIRECTION = Graph,
    output_var::Union{Nothing,VarOrAff} = nothing,
)
    return piecewiselinear(
        model,
        (input_var_x, input_var_y),
        pwl;
        method = method,
        direction = direction,
        output_vars = isnothing(output_var) ? nothing : (output_var,),
    )[1]
end

end # module
