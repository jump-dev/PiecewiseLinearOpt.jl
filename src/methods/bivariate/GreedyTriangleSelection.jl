# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    GreedyTriangleSelection(axis_method = Logarithmic(); tries = 200)

Formulation of bivariate continuous piecewise linear functions on any grid
triangulation, which selects the triangle with a few binary variables found by
a greedy heuristic.

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with a triangle selection like that of [`OptimalTriangleSelection`](@ref), but
the levels that select the triangle are found greedily instead of by a
mixed-integer program. The cells are inserted into the first level where they
fit, and the best of `tries` random orders of the cells is kept. This takes
milliseconds, needs no solver, and typically gives the optimal number of binary
variables: 2 for [`K1Triangulation`](@ref), 1 for
[`UnionJackTriangulation`](@ref), and a few more for irregular triangulations.
The result is cached for each triangulation.
"""
struct GreedyTriangleSelection{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
    tries::Int
end

function GreedyTriangleSelection(
    axis_method::PWLMethod = Logarithmic();
    tries::Integer = 200,
)
    return GreedyTriangleSelection(axis_method, Int(tries))
end

function _triangle_selection_cover(
    method::GreedyTriangleSelection,
    triangle_direction::Matrix{Bool},
)
    return _greedy_biclique_cover(triangle_direction, method.tries)
end
