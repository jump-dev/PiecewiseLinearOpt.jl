# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    PiecewiseLinear(f::Function, x; method = nothing, direction = Graph())
    PiecewiseLinear(x, z::AbstractVector; method = nothing, direction = Graph())
    PiecewiseLinear(f::Function, x, y; triangulation = K1Triangulation(), method = nothing, direction = Graph())
    PiecewiseLinear(x, y, z::AbstractMatrix; triangulation = K1Triangulation(), method = nothing, direction = Graph())
    PiecewiseLinear(f::Vector{<:Function}, x::Vector{<:AbstractVector}; method = nothing, direction = Graph())
    PiecewiseLinear(x::Vector{<:AbstractVector}, z::Vector{<:AbstractVector}; method = nothing, direction = Graph())
    PiecewiseLinear(pwl::PWLFunction; method = nothing, direction = Graph())

A piecewise linear function together with the method and direction to
formulate it with, which is called like a function. The univariate and
bivariate forms build the function as [`UnivariatePWLFunction`](@ref) and
[`BivariatePWLFunction`](@ref) do, but take the function first, so that it can
be given with a `do` block.

A univariate function can also be given in pieces, e.g. with a discontinuity
between them: piece `k` interpolates the function `f[k]`, or the values
`z[k]`, at the strictly increasing breakpoints `x[k]`. Each piece must start
where the previous piece ends, or later, which leaves a gap in the domain.
Where two pieces meet with different values, the model may take either value,
as a strict inequality such as `x < 0` cannot be modeled, while evaluating the
function gives the value of the left piece.

## Keyword arguments

- `method`: the formulation method. If `nothing`, the default is
  `Logarithmic()` for continuous univariate functions, `MultipleChoice()` for
  discontinuous ones, `DisaggregatedLogarithmic()` for relaxations, see
  [`UnivariatePWLRelaxation`](@ref), and `DisaggregatedLogarithmic()` for
  functions of three or more variables. For bivariate functions it depends
  on the triangulation: `K1()` for a `K1Triangulation`, `UnionJack()` for a
  `UnionJackTriangulation`, `GreedyTriangleSelection()` for any other grid
  triangulation and `DisaggregatedLogarithmic()` for unstructured segments.
- `direction`: whether the outputs equal (`Graph()`), are at least
  (`Epigraph()`) or are at most (`Hypograph()`) the function value.
- `triangulation`: how each cell of the grid of a bivariate function is split
  into two triangles, see [`GridTriangulation`](@ref).

## Calling

The object belongs to no model, and can be used in several. Calling it as

    f(x...; name = nothing, output_var = nothing)

with an input variable, affine expression or constant for each input of the
function adds the function to the model of the inputs, as
[`add_formulation`](@ref) does, and returns the output, or a tuple of outputs
for a function with several outputs. Each call adds a new formulation, so a
function of the same inputs that is used several times should be added once,
e.g. with `@expression(model, fx, f(x))`. If all inputs are constants and no
`output_var` is given, the function is evaluated instead.

The `name` is the prefix of the names of the variables created, and
`output_var` an existing output variable or expression of a function with one
output, see [`add_formulation`](@ref).

## Example

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x[1:2])
f = PiecewiseLinear(sin, 0:0.5:2π)
g = PiecewiseLinear(0:0.25:1, 0:0.25:1; method = SixStencil()) do a, b
    return (a - 1)^2 + (b - 0.5)^2
end
@expression(model, fx, f(x[1]))
@objective(model, Min, fx + g(x[1], x[2]))
@constraint(model, fx ≥ -0.5)
```
"""
struct PiecewiseLinear{
    D,
    F,
    P<:PWLFunction{D,F},
    M<:Union{Nothing,PWLMethod},
    R<:Direction,
}
    pwl::P
    method::M
    direction::R
end

function PiecewiseLinear(
    pwl::PWLFunction;
    method::Union{Nothing,PWLMethod} = nothing,
    direction::Direction = Graph(),
)
    return PiecewiseLinear(pwl, method, direction)
end

function PiecewiseLinear(f::Function, x::AbstractVector; kwargs...)
    return PiecewiseLinear(UnivariatePWLFunction(x, f); kwargs...)
end

function PiecewiseLinear(x::AbstractVector, z::AbstractVector; kwargs...)
    return PiecewiseLinear(UnivariatePWLFunction(x, z); kwargs...)
end

function PiecewiseLinear(
    f::AbstractVector{<:Function},
    x::AbstractVector{<:AbstractVector};
    kwargs...,
)
    if length(f) != length(x)
        error("Mismatch in the number of pieces of functions and breakpoints.")
    end
    z = [[fk(xi) for xi in xk] for (fk, xk) in zip(f, x)]
    return PiecewiseLinear(x, z; kwargs...)
end

function PiecewiseLinear(
    x::AbstractVector{<:AbstractVector},
    z::AbstractVector{<:AbstractVector};
    kwargs...,
)
    return PiecewiseLinear(_pieces_pwl_function(x, z); kwargs...)
end

function PiecewiseLinear(
    f::Function,
    x::AbstractVector,
    y::AbstractVector;
    triangulation::GridTriangulation = K1Triangulation(),
    kwargs...,
)
    pwl = BivariatePWLFunction(x, y, f; triangulation = triangulation)
    return PiecewiseLinear(pwl; kwargs...)
end

function PiecewiseLinear(
    x::AbstractVector,
    y::AbstractVector,
    z::AbstractMatrix;
    triangulation::GridTriangulation = K1Triangulation(),
    kwargs...,
)
    pwl = BivariatePWLFunction(x, y, z; triangulation = triangulation)
    return PiecewiseLinear(pwl; kwargs...)
end

# The default method depends on the number of inputs and, for univariate
# functions, on whether they are continuous, discontinuous or relaxations. For
# bivariate functions, it depends on the structure of the triangulation: the
# methods for the K1 and Union Jack triangulations need fewer binary variables
# than the general ones.
function _default_method(pwl::PWLFunction{1})
    if any(s -> length(s.input_vals) > 2, pwl.segments)
        return DisaggregatedLogarithmic()
    elseif !_is_continuous(pwl)
        return MultipleChoice()
    end
    return Logarithmic()
end
_default_method(::PWLFunction{2,F,K1Triangulation}) where {F} = K1()
function _default_method(::PWLFunction{2,F,UnionJackTriangulation}) where {F}
    return UnionJack()
end
function _default_method(::PWLFunction{2,F,<:GridTriangulation}) where {F}
    return GreedyTriangleSelection()
end
_default_method(::PWLFunction) = DisaggregatedLogarithmic()

# Whether the segments of a univariate function, each with two points, join
# into a continuous function on an interval
function _is_continuous(pwl::PWLFunction{1})
    left(s) = argmin(first.(s.input_vals))
    right(s) = argmax(first.(s.input_vals))
    sorted = sort(pwl.segments; by = s -> s.input_vals[left(s)][1])
    for (a, b) in zip(sorted[1:(end-1)], sorted[2:end])
        if a.input_vals[right(a)] != b.input_vals[left(b)] ||
           a.output_vals[right(a)] != b.output_vals[left(b)]
            return false
        end
    end
    return true
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

_method(f::PiecewiseLinear) = something(f.method, _default_method(f.pwl))

function Base.show(io::IO, f::PiecewiseLinear{D,F}) where {D,F}
    print(io, "PiecewiseLinear R^$D -> R^$F with $(length(f.pwl.segments))")
    print(io, " segments, ", _method(f))
    print(io, " and ", nameof(typeof(f.direction)))
    return
end

"""
    add_formulation(
        model,
        f::PiecewiseLinear,
        x...;
        name = nothing,
        output_vars = nothing,
    )

Adds the piecewise linear function `f` of the inputs `x` to `model`, and
returns the tuple `(output_vars, formulation)`. Calling `f(x...)` does the
same, but returns the output only.

## Arguments

- `model`: the JuMP model.
- `f`: the piecewise linear function, with the method and direction to
  formulate it with.
- `x`: an input variable, expression or constant for each input of `f`.

## Keyword arguments

- `name`: the prefix of the names of the variables created. If `nothing`, a
  name of the form `pwl1`, `pwl2`, ... is used.
- `output_vars`: a vector or tuple with an output variable or expression for
  each output of `f`. If `nothing`, new variables are created, bounded by
  the smallest and largest function value (only from below for `Epigraph()`
  and only from above for `Hypograph()`).

If the function cannot be formulated, an error is thrown and the model is left
unchanged.

## Returns

- `output_vars::Vector`: the outputs of the function.
- `formulation::PWLFormulation`: the variables and constraints that were added.

## Example

```julia
using JuMP, PiecewiseLinearOpt
model = Model()
@variable(model, x)
f = PiecewiseLinear(xi -> xi^2, 0:0.1:1)
(y,), formulation = add_formulation(model, f, x; name = "square")
```
"""
function add_formulation(
    model::JuMP.GenericModel,
    f::PiecewiseLinear{D,F},
    x::Scalar...;
    name::Union{Nothing,AbstractString} = nothing,
    output_vars::Union{Nothing,Tuple,AbstractVector} = nothing,
) where {D,F}
    if length(x) != D
        error(
            "Expected $D input variables for a R^$D -> R^$F piecewise linear function, got $(length(x)).",
        )
    end
    if output_vars !== nothing && length(output_vars) != F
        error(
            "Expected $F output variables for a R^$D -> R^$F piecewise linear function, got $(length(output_vars)).",
        )
    end
    pwl = f.pwl
    if isempty(pwl.segments)
        error(
            "I don't know how to handle a piecewise linear function with no breakpoints.",
        )
    end
    if name === nothing
        name = _next_default_name(model)
    end
    method = _method(f)
    formulation = PWLFormulation(model, String(name), method)

    try
        if output_vars === nothing
            output_vars = _create_output_vars(formulation, pwl, f.direction)
        end
        formulate_pwl!(
            formulation,
            x,
            Tuple(output_vars),
            pwl,
            method,
            f.direction,
        )
    catch
        # Leave the model as it was
        JuMP.delete(model, formulation)
        rethrow()
    end
    return collect(output_vars), formulation
end

function (f::PiecewiseLinear{D,F})(
    x::Scalar...;
    name::Union{Nothing,AbstractString} = nothing,
    output_var::Union{Nothing,Scalar} = nothing,
) where {D,F}
    if output_var !== nothing && F != 1
        error(
            "An output_var can only be given for a function with one output, use add_formulation with output_vars instead.",
        )
    end
    model = _owner_model((x..., something(output_var, 0)))
    if model === nothing
        if length(x) != D
            error(
                "Expected $D inputs for a R^$D -> R^$F piecewise linear function, got $(length(x)).",
            )
        end
        return f.pwl(map(_constant_value, x)...)
    end
    y, _ = add_formulation(
        model,
        f,
        x...;
        name = name,
        output_vars = output_var === nothing ? nothing : (output_var,),
    )
    return F == 1 ? only(y) : Tuple(y)
end

# The model of the inputs, or nothing if they are all constants
function _owner_model(x::Tuple)
    model = nothing
    for xi in x
        m = xi isa Real ? nothing : JuMP.owner_model(xi)
        if m === nothing
            continue
        elseif model === nothing
            model = m
        elseif m !== model
            error(
                "The inputs of a piecewise linear function belong to different models.",
            )
        end
    end
    return model
end

_constant_value(x::Real) = x

# An expression without variables
_constant_value(x::JuMP.AbstractJuMPScalar) = JuMP.value(_ -> NaN, x)
