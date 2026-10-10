# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Generalize to multivariate case.
struct OptimalTriangleSelection <: Method
    sub_solver::Any
    axis_method::Method
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
    _::GridTriangulation,
)
    cover =
        _optimal_biclique_cover(method.sub_solver, triangle_direction, false)
    _add_biclique_cover!(formulation, λ, cover)
    return nothing
end
