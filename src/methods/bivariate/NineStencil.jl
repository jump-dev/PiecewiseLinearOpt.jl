# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    NineStencil(axis_method = Logarithmic())

Formulation of bivariate continuous piecewise linear functions on any grid
triangulation (Huchette and Vielma, 2023).

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with 9 binary variables to select the triangle in the chosen cell. Fewer are
needed on small grids, as levels that separate nothing are left out.
"""
struct NineStencil{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
end
NineStencil() = NineStencil(Logarithmic())

function _triangle_selection_cover(
    ::NineStencil,
    triangle_direction::Matrix{Bool},
)
    n_1, n_2 = size(triangle_direction) .+ 1
    cover = _BicliqueCover()

    # A level for each of the nine residue classes (o_1, o_2) of the grid
    # points modulo 3. It separates the points of the class from their
    # diagonal neighbors that are not joined to them by a triangle edge.
    for o_1 in 1:3, o_2 in 1:3
        stencil = vec([(i, j) for i in o_1:3:n_1, j in o_2:3:n_2])
        has_edge_across_stencil = Set{Tuple{Int,Int}}()
        for (i, j) in stencil
            if i < n_1 && j < n_2 && triangle_direction[i, j]
                push!(has_edge_across_stencil, (i + 1, j + 1))
            end
            if i < n_1 && j > 1 && !triangle_direction[i, j-1]
                push!(has_edge_across_stencil, (i + 1, j - 1))
            end
            if i > 1 && j > 1 && triangle_direction[i-1, j-1]
                push!(has_edge_across_stencil, (i - 1, j - 1))
            end
            if i > 1 && j < n_2 && !triangle_direction[i-1, j]
                push!(has_edge_across_stencil, (i - 1, j + 1))
            end
        end
        push!(cover, (sort!(collect(has_edge_across_stencil)), stencil))
    end
    return cover
end
