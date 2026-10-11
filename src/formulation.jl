# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    PWLFormulation

The variables and constraints that [`piecewiselinear`](@ref) added to a model
to formulate a piecewise linear function.

## Fields

- `model`: the model
- `name`: the prefix of the names of the variables
- `method`: the formulation method
- `variables`: the variables added, including the output variables unless
  they were given by the user
- `constraints`: the constraints added, excluding bounds and integrality
  restrictions of the variables

Use `JuMP.delete(model, formulation)` to remove the function from the model.
"""
struct PWLFormulation{
    M<:JuMP.GenericModel,
    V<:JuMP.AbstractVariableRef,
    T<:Method,
}
    model::M
    name::String
    method::T
    variables::Vector{V}
    constraints::Vector{JuMP.ConstraintRef}
end

function PWLFormulation(model::JuMP.GenericModel, name::String, method::Method)
    return PWLFormulation(
        model,
        name,
        method,
        JuMP.variable_ref_type(model)[],
        JuMP.ConstraintRef[],
    )
end

function Base.show(io::IO, formulation::PWLFormulation)
    print(
        io,
        "PWLFormulation \"$(formulation.name)\" using $(formulation.method) ",
        "with $(length(formulation.variables)) variables and ",
        "$(length(formulation.constraints)) constraints",
    )
    return
end

function JuMP.delete(model::JuMP.GenericModel, formulation::PWLFormulation)
    if model !== formulation.model
        error("The formulation does not belong to the model.")
    end
    # The constraints are of different types, so delete them one at a time
    for c in formulation.constraints
        JuMP.delete(model, c)
    end
    JuMP.delete(model, formulation.variables)
    empty!(formulation.constraints)
    empty!(formulation.variables)
    return
end

# Every variable and constraint the methods create is registered in the
# formulation by passing it through _record!.
function _record!(
    formulation::PWLFormulation,
    x::Union{JuMP.AbstractVariableRef,JuMP.ConstraintRef},
)
    _push!(formulation, x)
    return x
end

function _record!(formulation::PWLFormulation, x::AbstractArray)
    for xi in x
        _push!(formulation, xi)
    end
    return x
end

function _record!(
    formulation::PWLFormulation,
    x::JuMP.Containers.SparseAxisArray,
)
    for xi in values(x.data)
        _push!(formulation, xi)
    end
    return x
end

function _push!(formulation::PWLFormulation, x::JuMP.AbstractVariableRef)
    push!(formulation.variables, x)
    return
end

function _push!(formulation::PWLFormulation, c::JuMP.ConstraintRef)
    push!(formulation.constraints, c)
    return
end

# The base name of a variable created by the formulation
_name(formulation::PWLFormulation, name::String) = "$(formulation.name)_$name"

# A formulation for a part of the model, e.g. one axis of a bivariate function,
# whose variables have names with the additional prefix `name`. It records into
# the same lists as `formulation`.
function _subformulation(formulation::PWLFormulation, name::String)
    return PWLFormulation(
        formulation.model,
        _name(formulation, name),
        formulation.method,
        formulation.variables,
        formulation.constraints,
    )
end
