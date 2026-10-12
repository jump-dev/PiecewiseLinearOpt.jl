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

The program only searches for fewer binary variables than the greedy selection
of [`GreedyTriangleSelection`](@ref) uses, which is often optimal. Proving this
takes about a second for regular triangulations of grids up to 8×8, but can
take much longer for irregular ones. If it is not done within `time_limit`
seconds, a warning is given and the greedy selection is used. This result is
also cached, so the time limit is only spent again with a longer time limit.
`GreedyTriangleSelection` is therefore usually the better choice for large or
irregular grids.

HiGHS sets up its threads at the first solve in a Julia session. If HiGHS is
both the sub-solver and the solver of a main model with a single thread, also
give the sub-solver a single thread, i.e.
`optimizer_with_attributes(HiGHS.Optimizer, "threads" => 1)`, or the main
model fails to solve.
"""
struct OptimalTriangleSelection{S,M<:Method} <: TriangleSelectionMethod
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

function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::OptimalTriangleSelection,
)
    greedy = _greedy_biclique_cover(triangle_direction, 200)
    cover, optimal, computed = _optimal_biclique_cover(
        method.sub_solver,
        triangle_direction,
        false,
        method.time_limit;
        upper_bound = greedy,
    )
    if !optimal && computed
        @warn(
            "OptimalTriangleSelection could not find the optimal triangle selection within the time limit of $(method.time_limit) seconds, using a greedy one with $(length(cover)) binary variables instead.",
        )
    end
    add_biclique_cover!(formulation, λ, cover)
    return nothing
end
