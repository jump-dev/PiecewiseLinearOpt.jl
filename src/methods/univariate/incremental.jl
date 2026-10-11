# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Implement a bivariate version of the incremental formulation, see
# https://github.com/jump-dev/PiecewiseLinearOpt.jl/issues/23
"""
    Incremental()

The incremental formulation of univariate continuous piecewise linear
functions (see Vielma, Ahmed and Nemhauser, 2010).

For `d` segments, it uses `d - 1` binary variables and `d` continuous variables
for the fill of each segment.
"""
struct Incremental <: Method end

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::Tuple{_Scalar},
    output_vars::NTuple{F,_Scalar},
    pwl::PWLFunction{1,F},
    method::Incremental,
    direction::DIRECTION,
) where {F}
    model = formulation.model
    grid = _continuous_gridpoints_or_die(pwl)
    xs, ys = grid.input_vals, grid.output_vals

    d = length(pwl.segments)
    @assert length(xs) == length(ys) == d + 1

    # δ[i] is the fill of segment i, and z[i] = 1 if segment i is full
    δ = _record!(
        formulation,
        JuMP.@variable(
            model,
            [1:d],
            lower_bound = 0,
            upper_bound = 1,
            base_name = _name(formulation, "δ")
        ),
    )
    z = _record!(
        formulation,
        JuMP.@variable(
            model,
            [1:d-1],
            Bin,
            base_name = _name(formulation, "z")
        ),
    )
    _record!(
        formulation,
        JuMP.@constraint(
            model,
            input_vars[1] ==
            xs[1][1] + sum(δ[i] * (xs[i+1][1] - xs[i][1]) for i in 1:d)
        ),
    )
    for j in 1:F
        rhs = ys[1][j] + sum(δ[i] * (ys[i+1][j] - ys[i][j]) for i in 1:d)
        _constrain_output_var(formulation, output_vars[j], rhs, direction)
    end
    for i in 1:d-1
        _record!(formulation, JuMP.@constraint(model, δ[i+1] ≤ z[i]))
        _record!(formulation, JuMP.@constraint(model, z[i] ≤ δ[i]))
    end
    return
end
