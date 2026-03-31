# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# TODO: Implement bivariate version of the incremental formulation
"""
    Incremental()

The incremental MIP formulation for univariate piecewise linear functions.

Uses `n - 1` binary variables for a function with `n` breakpoints. The formulation
is based on an incremental representation where each binary variable indicates
whether the function has "passed" a given breakpoint.
"""
struct Incremental <: Method end

function formulate_pwl!(
    model::JuMP.Model,
    input_vars::Tuple{VarOrAff},
    output_vars::NTuple{F,VarOrAff},
    pwl::PWLFunctionPointRep{1,F},
    method::Incremental,
    direction::DIRECTION,
) where {F}
    grid = _continuous_gridpoints_or_die(pwl)
    xs, ys = grid.input_vals, grid.output_vals

    counter = model.ext[:PWL].counter
    n = length(pwl.segments) + 1
    @assert length(xs) == length(ys) == n

    δ = JuMP.@variable(
        model,
        [1:n],
        lower_bound = 0,
        upper_bound = 1,
        base_name = "δ_$counter"
    )
    z = JuMP.@variable(model, [1:n-1], Bin, base_name = "z_$counter")
    JuMP.@constraint(
        model,
        input_vars[1] ==
        xs[1][1] + sum(δ[i] * (xs[i+1][1] - xs[i][1]) for i in 1:n-1)
    )
    for j in 1:F
        rhs = ys[1][j] + sum(δ[i] * (ys[i+1][j] - ys[i][j]) for i in 1:n-1)
        _constrain_output_var(model, output_vars[j], rhs, direction)
    end
    for i in 1:n-1
        JuMP.@constraint(model, δ[i+1] ≤ z[i])
        JuMP.@constraint(model, z[i] ≤ δ[i])
    end
    return
end
