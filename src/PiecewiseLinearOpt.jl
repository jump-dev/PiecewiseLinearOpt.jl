__precompile__()

module PiecewiseLinearOpt

import JuMP
import MathOptInterface
const MOI = MathOptInterface
using LinearAlgebra
using Random

export piecewiselinear

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

include(joinpath("methods", "util.jl"))

export Incremental, LogarithmicEmbedding, LogarithmicIndependentBranching, NativeSOS2, ZigZagBinary, ZigZagInteger
include(joinpath("methods", "univariate", "incremental.jl"))
include(joinpath("methods", "univariate", "logarithmic_embedding.jl"))
include(joinpath("methods", "univariate", "logarithmic_independent_branching.jl"))
include(joinpath("methods", "univariate", "native_sos2.jl"))
include(joinpath("methods", "univariate", "zig_zag_binary.jl"))
include(joinpath("methods", "univariate", "zig_zag_integer.jl"))

# Consider the colloqial "log" to refer to the embedding formulation
const Logarithmic = LogarithmicEmbedding
const SOS2Method = Union{LogarithmicEmbedding, LogarithmicIndependentBranching, NativeSOS2, ZigZagBinary, ZigZagInteger}
export Logarithmic

export K1, NineStencil, SixStencil, UnionJack
include(joinpath("methods", "bivariate", "k1.jl"))
include(joinpath("methods", "bivariate", "nine_stencil.jl"))
include(joinpath("methods", "bivariate", "six_stencil.jl"))
include(joinpath("methods", "bivariate", "union_jack.jl"))

export ConvexCombination, DisaggregatedLogarithmic, MultipleChoice, OptimalIndependentBranching, OptimalTriangleSelection
include(joinpath("methods", "multivariate", "convex_combination.jl"))
include(joinpath("methods", "multivariate", "disaggregated_logarithmic.jl"))
include(joinpath("methods", "multivariate", "multiple_choice.jl"))
include(joinpath("methods", "multivariate", "optimal_independent_branching.jl"))
include(joinpath("methods", "multivariate", "optimal_triangle_selection.jl"))

function formulate_pwl!(model::JuMP.Model, input_vals::Vector{NTuple{D,VarOrAff}}, output_vals::Vector{NTuple{F,VarOrAff}}, pwl::PWLFunction, method::Method, direction::DIRECTION) where {D,F}
    error("No support for a R^$D -> R^$F piecewise linear function using the $method method.")
end

_default_method(::Val{1}) = Logarithmic()
_default_method(::Val{2}) = SixStencil()
# _default_method(::Val) = MultipleChoice()

function piecewiselinear(model::JuMP.Model,
                         input_vars::NTuple{D,VarOrAff},
                         pwl::PWLFunction{D,F,SegmentPointRep{D,F}};
                         method::Method = _default_method(Val(D)),
                         direction::DIRECTION = Graph,
                         output_vars::Union{Nothing,NTuple{F,VarOrAff}} = nothing) where {D,F}
    initPWL!(model)
    counter = model.ext[:PWL].counter
    counter += 1
    model.ext[:PWL].counter = counter

    if isempty(pwl.segments)
        error(
            "I don't know how to handle a piecewise linear function with no breakpoints."
        )
    end

    output_lb = minimum(minimum(segment.output_vals) for segment in pwl.segments)
    output_ub = maximum(maximum(segment.output_vals) for segment in pwl.segments)

    if output_vars === nothing
        output_vars = tuple(JuMP.@variable(model, [i in 1:F], lower_bound=output_lb[i], upper_bound=output_ub[i], base_name="y_$counter")...)
    end

    formulate_pwl!(model, input_vars, output_vars, pwl, method, direction)
    return output_vars
end

end # module
