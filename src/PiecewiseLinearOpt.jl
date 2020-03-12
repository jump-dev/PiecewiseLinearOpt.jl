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

export Incremental, Logarithmic
include(joinpath("methods", "incremental.jl"))
include(joinpath("methods", "logarithmic.jl"))

function _formulate!(model::JuMP.Model, input_vals::Vector{NTuple{D,VarOrAff}}, output_vals::Vector{NTuple{F,VarOrAff}}, pwl::PWLFunction, method::Method, direction::DIRECTION) where {D,F}
    error("No support for a R^$D -> R^$F piecewise linear function using the $method method.")
end

function piecewiselinear(model::JuMP.Model,
                         input_var::VarOrAff,
                         pwl::UnivariatePWLFunction;
                         method = Logarithmic(),
                         direction::DIRECTION = Graph,
                         output_var::Union{VarOrAff,Nothing} = nothing)
    initPWL!(model)
    counter = model.ext[:PWL].counter
    counter += 1
    model.ext[:PWL].counter = counter

    if isempty(pwl.segments)
        error(
            "I don't know how to handle a piecewise linear function with no breakpoints."
        )
    end

    output_lb = minimum(minimum(segment.output_vals) for segment in pwl.segments)[1]
    output_ub = maximum(maximum(segment.output_vals) for segment in pwl.segments)[1]

    if output_var === nothing
        output_var = JuMP.@variable(model, lower_bound=output_lb, upper_bound=output_ub, base_name="y_$counter")
    end

    _formulate!(model, (input_var,), (output_var,), pwl, method, direction)
    return output_var
end

end # module
