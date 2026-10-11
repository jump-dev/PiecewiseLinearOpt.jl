# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The fallback for any method and function it does not support. Its arguments
# other than the method are untyped, so that any implementation for a specific
# method is more specific than it.
function formulate_pwl!(
    formulation,
    input_vars,
    output_vars,
    pwl,
    method::Method,
    direction,
)
    return error(
        "No support for a R^$(length(input_vars)) -> R^$(length(output_vars)) piecewise linear function using the $method method.",
    )
end

_default_method(::Val{1}) = Logarithmic()
_default_method(::Val{2}) = SixStencil()
_default_method(::Val) = DisaggregatedLogarithmic()

"""
    piecewiselinear(
        model,
        input_vars,
        pwl::PWLFunction;
        method,
        direction = Graph(),
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
  functions, `SixStencil()` for bivariate functions and
  `DisaggregatedLogarithmic()` for functions of more variables.
- `direction`: whether the outputs equal (`PiecewiseLinearOpt.Graph()`), are
  at least (`Epigraph()`) or are at most (`Hypograph()`) the function value.
- `output_vars`: a vector or tuple with an output variable or expression for
  each output of `pwl`. If `nothing`, new variables are created, bounded by
  the smallest and largest function value (only from below for `Epigraph()`
  and only from above for `Hypograph()`).
- `name`: the prefix of the names of the variables created. If `nothing`, a
  name of the form `pwl1`, `pwl2`, ... is used.

If the function cannot be formulated, an error is thrown and the model is left
unchanged.

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
    direction::Direction = Graph(),
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

    try
        if output_vars === nothing
            output_vars = _create_output_vars(formulation, pwl, direction)
        end
        formulate_pwl!(
            formulation,
            Tuple(input_vars),
            Tuple(output_vars),
            pwl,
            method,
            direction,
        )
    catch
        # Leave the model as it was
        JuMP.delete(model, formulation)
        rethrow()
    end
    return collect(output_vars), formulation
end

# The outputs are bounded by the function values, from the sides that the
# direction implies.
function _create_output_vars(
    formulation::PWLFormulation,
    pwl::PWLFunction{D,F},
    direction::Direction,
) where {D,F}
    y = record!(
        formulation,
        JuMP.@variable(
            formulation.model,
            [1:F],
            base_name = variable_name(formulation, "y")
        ),
    )
    for i in 1:F
        if !(direction isa Hypograph)
            JuMP.set_lower_bound(
                y[i],
                minimum(v[i] for seg in pwl.segments for v in seg.output_vals),
            )
        end
        if !(direction isa Epigraph)
            JuMP.set_upper_bound(
                y[i],
                maximum(v[i] for seg in pwl.segments for v in seg.output_vals),
            )
        end
    end
    return y
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
