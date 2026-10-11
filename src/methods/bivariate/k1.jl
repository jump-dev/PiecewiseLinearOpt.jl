# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    K1(axis_method = Logarithmic())

Formulation of bivariate continuous piecewise linear functions on a grid where
every cell is split along the diagonal from its upper left to its lower right
corner, e.g. [`K1Triangulation`](@ref).

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with 2 binary variables to select the triangle in the chosen cell.
"""
struct K1{M<:Method} <: TriangleSelectionMethod
    axis_method::M
end
K1() = K1(Logarithmic())

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::K1,
)
    model = formulation.model
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)

    # The formulation assumes every cell is split along its NW-SE diagonal.
    if !all(triangle_direction)
        error(
            "The K1 method requires a K1 triangulation, where every cell is split along the diagonal from its upper left to its lower right corner.",
        )
    end

    w = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:2],
            Bin,
            base_name = variable_name(formulation, "w")
        ),
    )
    # Points on the diagonals i + j = const of the same parity, alternating
    # between the two sides of each level
    for (k, parity, offset) in ((1, true, 2), (2, false, 3))
        A = [
            λ[i, j] for
            i in 1:n_1, j in 1:n_2 if (mod(i, 2) == mod(j, 2)) == parity &&
            (i + j) in offset:4:(n_1+n_2)
        ]
        B = [
            λ[i, j] for
            i in 1:n_1, j in 1:n_2 if (mod(i, 2) == mod(j, 2)) == parity &&
            (i + j) in offset+2:4:(n_1+n_2)
        ]
        if !isempty(A)
            record!(formulation, JuMP.@constraint(model, sum(A) ≤ w[k]))
        end
        if !isempty(B)
            record!(formulation, JuMP.@constraint(model, sum(B) ≤ 1 - w[k]))
        end
    end
    return nothing
end
