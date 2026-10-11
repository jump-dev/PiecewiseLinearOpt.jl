# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
"""
    OptimalTriangleSelection(
        sub_solver,
        axis_method = Logarithmic();
        time_limit = 60.0,
    )

Formulation of bivariate continuous piecewise linear functions on any grid
triangulation (Huchette and Vielma, 2023).

It combines an SOS2 constraint along each axis, formulated with `axis_method`,
with the smallest number of binary variables that select the triangle in the
chosen cell, typically 2 to 4. That number is found by solving a mixed-integer
program with `sub_solver`, e.g. `HiGHS.Optimizer`, and the result is cached
for each triangulation.

This takes about a second for an 8×8 grid with a regular triangulation, but can
take much longer for irregular ones. If no solution is found within
`time_limit` seconds, a warning is given and the triangle is selected as in
[`SixStencil`](@ref) instead, with 6 binary variables. As the optimal number of
binary variables is typically only slightly smaller, `SixStencil` is usually
the better choice for large or irregular grids.
"""
struct OptimalTriangleSelection{S,M<:Method} <: Method
    sub_solver::S
    axis_method::M
    time_limit::Float64
end

function OptimalTriangleSelection(
    sub_solver,
    axis_method::Method = Logarithmic();
    time_limit::Real = 60.0,
)
    return OptimalTriangleSelection(
        sub_solver,
        axis_method,
        Float64(time_limit),
    )
end

axis_method(method::OptimalTriangleSelection) = method.axis_method

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::OptimalTriangleSelection,
)
    cover = _optimal_biclique_cover(
        method.sub_solver,
        triangle_direction,
        false,
        method.time_limit,
    )
    if cover === nothing
        @warn(
            "OptimalTriangleSelection found no triangle selection within the time limit of $(method.time_limit) seconds, using SixStencil instead.",
        )
        return formulate_triangle_selection!(
            formulation,
            λ,
            triangle_direction,
            SixStencil(method.axis_method),
        )
    end
    _add_biclique_cover!(formulation, λ, cover)
    return nothing
end
