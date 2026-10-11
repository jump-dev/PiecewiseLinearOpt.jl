# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    NativeSOS2()

Formulates univariate continuous piecewise linear functions with the convex
multipliers in a special ordered set of type 2 (`MOI.SOS2`), leaving the
branching to the solver. It adds no integer variables.

For solvers without support for SOS2 constraints, such as HiGHS, MathOptInterface
reformulates the constraint with binary variables. For solvers with support,
e.g. Gurobi, it can also be a good axis method for the bivariate formulations.
"""
struct NativeSOS2 <: Method end

function formulate_sos2!(
    formulation::PWLFormulation,
    λ::Vector{T},
    method::NativeSOS2,
) where {T<:_Scalar}
    _record!(
        formulation,
        JuMP.@constraint(
            formulation.model,
            λ in JuMP.SOS2([k for k in 1:length(λ)])
        ),
    )
    return nothing
end
