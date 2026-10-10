# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
struct OptimalIndependentBranching <: Method
    sub_solver::Any
end

function formulate_pwl!(
    model::JuMP.Model,
    input_vars::NTuple{2,VarOrAff},
    output_vars::NTuple{F,VarOrAff},
    pwl::PWLFunctionPointRep{2,F},
    method::OptimalIndependentBranching,
    direction::DIRECTION,
) where {F}
    λ, triangle_direction = _grid_convex_multipliers!(
        model,
        input_vars,
        output_vars,
        pwl,
        direction,
    )
    cover = _optimal_biclique_cover(method.sub_solver, triangle_direction, true)
    _add_biclique_cover!(model, λ, cover)
    return nothing
end
