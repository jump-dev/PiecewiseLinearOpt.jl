# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct UnionJack <: Method
    axis_method::Method
end
UnionJack() = UnionJack(Logarithmic())

axis_method(method::UnionJack) = method.axis_method

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::UnionJack,
    _::UnionJackTriangulation,
)
    model = formulation.model
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)

    # In a Union Jack triangulation the diagonals alternate between neighboring
    # cells.
    for i in 1:n_1-1, j in 1:n_2-1
        if triangle_direction[i, j] !=
           xor(triangle_direction[1, 1], isodd(i + j))
            error(
                "The UnionJack method requires a Union Jack triangulation, but the diagonals of cell ($i, $j) and cell (1, 1) do not match.",
            )
        end
    end

    # j_start = number of triangle segments incident to (1, 1), the lower-left
    # most point in the grid.
    j_start = triangle_direction[1, 1] ? 1 : 2

    # diagonal lines running from SW to NE. Grouped with an offset of 3.
    z = _record!(
        formulation,
        JuMP.@variable(
            model,
            binary = true,
            base_name = _name(formulation, "w")
        ),
    )
    _record!(
        formulation,
        JuMP.@constraint(
            model,
            sum(λ[i, j] for i in 1:2:n_1, j in j_start:2:n_2) ≤ z
        ),
    )
    _record!(
        formulation,
        JuMP.@constraint(
            model,
            sum(λ[i, j] for i in 2:2:n_1, j in (3-j_start):2:n_2) ≤ 1 - z
        ),
    )

    return nothing
end
