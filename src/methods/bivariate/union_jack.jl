# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    UnionJack()
    UnionJack(axis_method::Method)

Bivariate formulation for Union Jack-triangulated grids.

Uses only 1 binary variable for triangle selection plus the binary variables
from `axis_method` for the axis-aligned SOS2 constraints. Requires the
piecewise linear function to use a Union Jack triangulation
(`pattern = :UnionJack`).

The default `axis_method` is [`Logarithmic`](@ref).
"""
struct UnionJack <: Method
    axis_method::Method
end
UnionJack() = UnionJack(Logarithmic())

axis_method(method::UnionJack) = method.axis_method

function formulate_triangle_selection!(
    model::JuMP.Model,
    λ::Matrix{JuMP.VariableRef},
    triangle_direction::Matrix{Bool},
    method::UnionJack,
    _::UnionJackTriangulation,
)
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)
    # TODO: Certify triangulation is Union Jack

    # j_start = number of triangle segments incident to (1, 1), the lower-left
    # most point in the grid.
    j_start = triangle_direction[1, 1] ? 1 : 2

    # diagonal lines running from SW to NE. Grouped with an offset of 3.
    z = JuMP.@variable(model, binary = true, base_name = _pwl_name(model, "w"))

    JuMP.@constraints(
        model,
        begin
            sum(λ[i, j] for i in 1:2:n_1, j in j_start:2:n_2) ≤ z
            sum(λ[i, j] for i in 2:2:n_1, j in (3-j_start):2:n_2) ≤ 1 - z
        end
    )

    return nothing
end
