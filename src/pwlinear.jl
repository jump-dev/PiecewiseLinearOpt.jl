# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

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

"""
    piecewiselinear(model, input_vars, pwl::PWLFunction; method, direction, output_vars)

    piecewiselinear(model, input_var, x, f::Function; method, direction, output_var)
    piecewiselinear(model, input_var, x, fd; method, direction, output_var)

    piecewiselinear(model, input_var_x, input_var_y, x, y, f::Function; method, direction, output_var, pattern)


Formulates a piecewise linear function in the given JuMP model.

# Arguments
- `model::JuMP.Model`: The JuMP model to which the piecewise linear function will be added.
- `input_vars::NTuple{D,VarOrAff}`: A tuple of input variables or affine expressions.
- `pwl::PWLFunction{D,F,SegmentPointRep{D,F}}`: The piecewise linear function to be added.
- `method::Method`: The method to be used for formulating the piecewise linear function.
   Defaults to `Logarithmic()` for 1D and `SixStencil()` for 2D.
- `direction::DIRECTION`: The direction of the piecewise linear function. Defaults to `Graph`.
- `output_vars::Union{Nothing,NTuple{F,VarOrAff}}`: A tuple of output variables or affine expressions.
   If `nothing`, new variables will be created.

# Returns
- `output_vars::NTuple{F,VarOrAff}`: The output variables of the piecewise linear function.

In addition to the general constructor, there are specialized constructors for univariate
and bivariate problems.

# Example
```julia
using JuMP
model = Model()
@variable(model, x >= 0)
@variable(model, y >= 0)

pwl = UnivariatePWLFunction(0:0.1:1, xi -> xi^2)

# General constructor
output_vars = piecewiselinear(model, (x,), pwl)

# Specialized constructor for univariate problems
output_var = piecewiselinear(model, x, 0:0.1:1, xi -> xi^2)

# Specialized constructor for bivariate problems
output_var = piecewiselinear(model, x, y, 0:0.1:1, 0:0.1:1, (xi, yi) -> xi^2 + yi^2; pattern = :BestFit)
```
"""
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

# Specialization for univariate problems
function piecewiselinear(
    model::JuMP.Model,
    input_var::VarOrAff,
    d,
    f::Function;
    method::Method = _default_method(Val(1)),
    direction::DIRECTION = Graph,
    output_var::Union{Nothing,VarOrAff} = nothing,
)
    return piecewiselinear(
        model,
        (input_var,),
        UnivariatePWLFunction(d, f);
        method = method,
        direction = direction,
        output_vars = isnothing(output_var) ? nothing : (output_var,),
    )[1]
end

function piecewiselinear(
    model::JuMP.Model,
    input_var::VarOrAff,
    d,
    fd;
    method::Method = _default_method(Val(1)),
    direction::DIRECTION = Graph,
    output_var::Union{Nothing,VarOrAff} = nothing,
)
    return piecewiselinear(
        model,
        (input_var,) ,
        UnivariatePWLFunction(d, fd);
        method = method,
        direction = direction,
        output_vars = isnothing(output_var) ? nothing : (output_var,),
    )[1]
end

# Specialization for bivariate problems
function piecewiselinear(
    model::JuMP.Model,
    input_var_x::VarOrAff,
    input_var_y::VarOrAff,
    x,
    y,
    f::Function;
    method::Method = _default_method(Val(2)),
    direction::DIRECTION = Graph,
    output_var::Union{Nothing,VarOrAff} = nothing,
    pattern = :K1,
)
    return piecewiselinear(
        model,
        (input_var_x, input_var_y),
        BivariatePWLFunction(x, y, f; pattern = UnstructuredTriangulation());
        method = method,
        direction = direction,
        output_vars = isnothing(output_var) ? nothing : (output_var,),
    )[1]
end
