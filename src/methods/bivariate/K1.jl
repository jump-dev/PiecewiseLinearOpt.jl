# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    K1(axis_method = Logarithmic())

Formulation of bivariate continuous piecewise linear functions on a
[`K1Triangulation`](@ref), the default method for it. Every cell of such a
grid is split along the diagonal from its upper left to its lower right
corner.

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with 2 binary variables to select the triangle in the chosen cell.
"""
struct K1{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
end
K1() = K1(Logarithmic())

function _triangle_selection_cover(::K1, triangle_direction::Matrix{Bool})
    n_1, n_2 = size(triangle_direction) .+ 1

    # The formulation assumes every cell is split along its NW-SE diagonal.
    if !all(triangle_direction)
        error(
            "The K1 method requires a K1 triangulation, where every cell is split along the diagonal from its upper left to its lower right corner.",
        )
    end

    # Each level takes the points on the diagonals i + j = const of one
    # parity, alternating between the two sides of the level
    cover = _BicliqueCover()
    for (parity, offset) in ((true, 2), (false, 3))
        A = [
            (i, j) for
            i in 1:n_1, j in 1:n_2 if (mod(i, 2) == mod(j, 2)) == parity &&
                (i + j) in offset:4:(n_1+n_2)
        ]
        B = [
            (i, j) for
            i in 1:n_1, j in 1:n_2 if (mod(i, 2) == mod(j, 2)) == parity &&
                (i + j) in (offset+2):4:(n_1+n_2)
        ]
        push!(cover, (A, B))
    end
    return cover
end
