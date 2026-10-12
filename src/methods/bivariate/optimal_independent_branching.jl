# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
"""
    OptimalIndependentBranching(sub_solver; time_limit = 60.0)

Independent branching formulation with the fewest binary variables.

It formulates bivariate continuous piecewise linear functions on any grid
triangulation with the smallest number of binary variables (Huchette and
Vielma, 2019), at least `ceil(log2(t))` for `t` triangles.

That number is found by solving a mixed-integer program with `sub_solver`,
e.g. `HiGHS.Optimizer`, and the result is cached for each triangulation.
Computing it is only practical for small grids, up to about 4×4 points with a
regular triangulation, and an error is thrown if no solution is found within
`time_limit` seconds. This failure is also cached, so it is only tried again
with a longer time limit. For larger grids, use e.g.
[`GreedyTriangleSelection`](@ref) or [`OptimalTriangleSelection`](@ref).

See [`OptimalTriangleSelection`](@ref) on the number of threads when the
sub-solver is HiGHS.
"""
struct OptimalIndependentBranching{S} <: Method
    sub_solver::S
    time_limit::Float64
end

function OptimalIndependentBranching(sub_solver; time_limit::Real = 60.0)
    return OptimalIndependentBranching(sub_solver, Float64(time_limit))
end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{2,Scalar},
    output_vars::NTuple{F,Scalar},
    pwl::PWLFunction{2,F},
    method::OptimalIndependentBranching,
    direction::Direction,
) where {F}
    λ, triangle_direction = _grid_convex_multipliers!(
        formulation,
        input_vars,
        output_vars,
        pwl,
        direction,
    )
    cover, _, _ = _optimal_biclique_cover(
        method.sub_solver,
        triangle_direction,
        true,
        method.time_limit,
    )
    if cover === nothing
        n_1, n_2 = size(λ)
        error(
            "OptimalIndependentBranching found no formulation for the $(n_1)×$(n_2) grid within the time limit of $(method.time_limit) seconds. It is only practical for small grids; use e.g. SixStencil() or OptimalTriangleSelection(sub_solver) instead.",
        )
    end
    add_biclique_cover!(formulation, λ, cover)
    return nothing
end
