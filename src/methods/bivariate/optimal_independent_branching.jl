# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
"""
    OptimalIndependentBranching(sub_solver)

Independent branching formulation with the fewest binary variables.

It formulates bivariate continuous piecewise linear functions on any grid
triangulation with the smallest number of binary variables (Huchette and
Vielma, 2019), at least `ceil(log2(t))` for `t` triangles.

That number is found by solving a mixed-integer program with `sub_solver`,
e.g. `HiGHS.Optimizer`. The result is cached for each triangulation, but
computing it is only practical for small grids, up to about 4×4 points.
"""
struct OptimalIndependentBranching{S} <: Method
    sub_solver::S
end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{2,_Scalar},
    output_vars::NTuple{F,_Scalar},
    pwl::PWLFunction{2,F},
    method::OptimalIndependentBranching,
    direction::DIRECTION,
) where {F}
    λ, triangle_direction = _grid_convex_multipliers!(
        formulation,
        input_vars,
        output_vars,
        pwl,
        direction,
    )
    cover = _optimal_biclique_cover(method.sub_solver, triangle_direction, true)
    _add_biclique_cover!(formulation, λ, cover)
    return nothing
end
