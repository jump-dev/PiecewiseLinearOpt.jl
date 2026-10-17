# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    UnionJack(axis_method = Logarithmic())

Formulation of bivariate continuous piecewise linear functions on a
[`UnionJackTriangulation`](@ref), the default method for it (Vielma and
Nemhauser, 2011). The diagonals of such a grid alternate between neighboring
cells.

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with 1 binary variable to select the triangle in the chosen cell.
"""
struct UnionJack{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
end
UnionJack() = UnionJack(Logarithmic())

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::UnionJack,
)
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)

    # In a Union Jack triangulation the diagonals alternate between neighboring
    # cells.
    for i in 1:(n_1-1), j in 1:(n_2-1)
        if triangle_direction[i, j] !=
           xor(triangle_direction[1, 1], isodd(i + j))
            error(
                "The UnionJack method requires a Union Jack triangulation, but the diagonals of cell ($i, $j) and cell (1, 1) do not match.",
            )
        end
    end

    # One level with the points of odd i on one side and those of even i on
    # the other, each taking every other j, so that the ends of the diagonal
    # that is not an edge are on opposite sides in every cell. j_start is the
    # number of triangle segments incident to (1, 1), the lower left point of
    # the grid.
    j_start = triangle_direction[1, 1] ? 1 : 2
    L = vec([(i, j) for i in 1:2:n_1, j in j_start:2:n_2])
    R = vec([(i, j) for i in 2:2:n_1, j in (3-j_start):2:n_2])
    add_biclique_cover!(formulation, λ, [(L, R)])
    return nothing
end
