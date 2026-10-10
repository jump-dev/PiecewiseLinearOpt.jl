# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
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
