# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
"""
    OptimalTriangleSelection(sub_solver, axis_method = Logarithmic())

Formulation of bivariate continuous piecewise linear functions on any grid
triangulation (Huchette and Vielma, 2023).

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with the smallest number of binary variables that select the triangle in the
chosen cell, typically 2 or 3. That number is found by solving a mixed-integer
program with `sub_solver`, e.g. `HiGHS.Optimizer`, which takes about a second
for an 8×8 grid. The result is cached for each triangulation.
"""
struct OptimalTriangleSelection{S,M<:Method} <: Method
    sub_solver::S
    axis_method::M
end
function OptimalTriangleSelection(sub_solver)
    return OptimalTriangleSelection(sub_solver, Logarithmic())
end

axis_method(method::OptimalTriangleSelection) = method.axis_method

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::OptimalTriangleSelection,
)
    cover =
        _optimal_biclique_cover(method.sub_solver, triangle_direction, false)
    _add_biclique_cover!(formulation, λ, cover)
    return nothing
end
