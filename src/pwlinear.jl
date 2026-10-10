# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{D,_Scalar},
    output_vars::NTuple{F,_Scalar},
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

"""
    piecewiselinear(
        model,
        input_vars,
        pwl::PWLFunction;
        method,
        direction = Graph,
        output_vars = nothing,
        name = nothing,
    )

Adds the piecewise linear function `pwl` of `input_vars` to `model`, and
returns the tuple `(output_vars, formulation)`.

## Arguments

- `model`: the JuMP model.
- `input_vars`: a vector or tuple with an input variable, expression or
  constant for each input of `pwl`.
- `pwl`: the piecewise linear function.

## Keyword arguments

- `method`: the formulation method, defaults to `Logarithmic()` for univariate
  and `SixStencil()` for bivariate functions.
- `direction`: whether the outputs equal (`Graph`), are at least (`Epigraph`)
  or are at most (`Hypograph`) the function value.
- `output_vars`: a vector or tuple with an output variable or expression for
  each output of `pwl`. If `nothing`, new variables are created.
- `name`: the prefix of the names of the variables created. If `nothing`, a
  name of the form `pwl1`, `pwl2`, ... is used.

## Returns

- `output_vars::Vector`: the outputs of the function.
- `formulation::PWLFormulation`: the variables and constraints that were added.

There are also convenience methods for univariate and bivariate functions,
which return the output only:

    piecewiselinear(model, input_var, x, f::Function; kwargs...)
    piecewiselinear(model, input_var, x, fx; kwargs...)
    piecewiselinear(model, input_var_x, input_var_y, x, y, f::Function; triangulation, kwargs...)
    piecewiselinear(model, input_var_x, input_var_y, x, y, z::AbstractMatrix; triangulation, kwargs...)

## Example

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x[1:2])

pwl = UnivariatePWLFunction(0:0.1:1, xi -> xi^2)
y, formulation = piecewiselinear(model, x[1:1], pwl; name = "square")

z = piecewiselinear(model, x[1], 0:0.1:1, xi -> xi^2)
w = piecewiselinear(model, x[1], x[2], 0:0.1:1, 0:0.1:1, (a, b) -> a * b)
```
"""
function piecewiselinear(
    model::JuMP.GenericModel,
    input_vars::Union{Tuple,AbstractVector},
    pwl::PWLFunction{D,F};
    method::Method = _default_method(Val(D)),
    direction::DIRECTION = Graph,
    output_vars::Union{Nothing,Tuple,AbstractVector} = nothing,
    name::Union{Nothing,AbstractString} = nothing,
) where {D,F}
    if length(input_vars) != D
        error(
            "Expected $D input variables for a R^$D -> R^$F piecewise linear function, got $(length(input_vars)).",
        )
    end
    if output_vars !== nothing && length(output_vars) != F
        error(
            "Expected $F output variables for a R^$D -> R^$F piecewise linear function, got $(length(output_vars)).",
        )
    end
    if isempty(pwl.segments)
        error(
            "I don't know how to handle a piecewise linear function with no breakpoints.",
        )
    end
    if name === nothing
        name = _next_default_name(model)
    end
    formulation = PWLFormulation(model, String(name), method)

    if output_vars === nothing
        output_lb = ntuple(
            i -> minimum(v[i] for seg in pwl.segments for v in seg.output_vals),
            F,
        )
        output_ub = ntuple(
            i -> maximum(v[i] for seg in pwl.segments for v in seg.output_vals),
            F,
        )
        output_vars = _record!(
            formulation,
            JuMP.@variable(
                model,
                [i in 1:F],
                lower_bound = output_lb[i],
                upper_bound = output_ub[i],
                base_name = _name(formulation, "y")
            ),
        )
    end

    formulate_pwl!(
        formulation,
        Tuple(input_vars),
        Tuple(output_vars),
        pwl,
        method,
        direction,
    )
    return collect(output_vars), formulation
end

# Specialization for univariate problems
function piecewiselinear(
    model::JuMP.GenericModel,
    input_var::_Scalar,
    d,
    f::Function;
    output_var::Union{Nothing,_Scalar} = nothing,
    kwargs...,
)
    return piecewiselinear(
        model,
        input_var,
        d,
        [f(xi) for xi in d];
        output_var = output_var,
        kwargs...,
    )
end

function piecewiselinear(
    model::JuMP.GenericModel,
    input_var::_Scalar,
    d,
    fd;
    output_var::Union{Nothing,_Scalar} = nothing,
    kwargs...,
)
    y, _ = piecewiselinear(
        model,
        (input_var,),
        UnivariatePWLFunction(d, fd);
        output_vars = isnothing(output_var) ? nothing : (output_var,),
        kwargs...,
    )
    return only(y)
end

# Specialization for bivariate problems
function piecewiselinear(
    model::JuMP.GenericModel,
    input_var_x::_Scalar,
    input_var_y::_Scalar,
    x,
    y,
    f::Union{Function,AbstractMatrix};
    output_var::Union{Nothing,_Scalar} = nothing,
    triangulation::GridTriangulation = K1Triangulation(),
    kwargs...,
)
    z, _ = piecewiselinear(
        model,
        (input_var_x, input_var_y),
        BivariatePWLFunction(x, y, f; triangulation = triangulation);
        output_vars = isnothing(output_var) ? nothing : (output_var,),
        kwargs...,
    )
    return only(z)
end
