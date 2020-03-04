__precompile__()

module PiecewiseLinearOpt

import JuMP
import MathOptInterface
const MOI = MathOptInterface
using LinearAlgebra
using Random

export piecewiselinear

abstract type Method end
abstract type UnivariateMethod <: Method end

@enum DIRECTION Graph Epigraph Hypograph

# include(joinpath("methods", "sos2.jl"))

include("types.jl")
include("jump.jl")

function _formulate!(model::JuMP.Model, input_vars::Vector{NTuple{D,Float64}}, output_vars::Vector{NTuple{F,Float64}}, pwl::PWLFunction, method::Method, direction::DIRECTION)
    error("No support for a R^$D -> R^$F piecewise linear function using the $method method.")
end

function _constrain_output_var(model::JuMP.Model, output_var::VarOrAff, rhs::VarOrAff, direction::DIRECTION)
    if direction == Graph
        JuMP.@constraint(model, output_var == rhs)
    elseif direction == Hypograph
        JuMP.@constraint(model, output_var <= rhs)
    else
        @assert direction == Epigraph
        JuMP.@constraint(model, output_var >= rhs)
    end
    return
end

function _assert_continuous_or_die(pwl::UnivariatePWLFunction)
    V = union(pwl.segments.input_vals)
    for segment in pwl.segments

end

function _continuous_breakpoints_or_die(pwl::UnivariatePWLFunction{F}) where {F}
    V = union(pwl.segments.input_vals)
    sort!(V)
    fds = Array{Float64}(undef, length(V))
    _fds = Dict{Float64, NTuple{F,Float64}}()
    for segment in pwl.segments
        for input_val in segment.input_vals
            if haskey(_fds, segment.input_vals)
                if fds[input_val] != segment.output_vals
                    error("Encountered an discontinuous piecewise linear function; aborting.")
                end
            else
                _fds[input_val] = segment.ouput_vals
            end
        end
    end
    fds = [_fds[x] for x in V]
    return d, fds
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
            "I don't know how to handle a piecewise linear function with no"
            "breakpoints."
        )
    end

    output_lb = minimum(output_val[1] for output_val in pwl.segments.output_vals)
    output_ub = maximum(output_val[1] for output_val in pwl.segments.output_vals)

    if output_var === nothing
        output_var = JuMP.@variable(model, lower_bound=output_ub, upper_bound=output_lb, base_name="y_$counter")
    end

    _formulate!(model, input_var, output_var, pwl, method, direction)
    return output_var
end

end # module
