# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

struct NineStencil <: Method
    axis_method::Method
end
NineStencil() = NineStencil(Logarithmic())

axis_method(method::NineStencil) = method.axis_method

# TODO: Unit tests for biclique cover
function formulate_triangle_selection!(
    formulation::PWLFormulation,
    λ::Matrix{<:JuMP.AbstractVariableRef},
    triangle_direction::Matrix{Bool},
    method::NineStencil,
)
    model = formulation.model
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)

    w = _record!(
        formulation,
        JuMP.@variable(
            model,
            [1:3, 1:3],
            Bin,
            base_name = _name(formulation, "w")
        ),
    )
    for o_1 in 1:3, o_2 in 1:3
        has_edge_across_stencil = Set{Tuple{Int,Int}}()
        for i in o_1:3:n_1, j in o_2:3:n_2
            let k = i + 1, l = j + 1
                if (1 ≤ k ≤ n_1) && (1 ≤ l ≤ n_2)
                    if triangle_direction[i, j]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
            let k = i + 1, l = j - 1
                if (1 ≤ k ≤ n_1) && (1 ≤ l ≤ n_2)
                    if !triangle_direction[i, j-1]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
            let k = i - 1, l = j - 1
                if (1 ≤ k ≤ n_1) && (1 ≤ l ≤ n_2)
                    if triangle_direction[i-1, j-1]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
            let k = i - 1, l = j + 1
                if (1 ≤ k ≤ n_1) && (1 ≤ l ≤ n_2)
                    if !triangle_direction[i-1, j]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
        end
        _record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(λ[i, j] for i in o_1:3:n_1, j in o_2:3:n_2) ≤
                1 - w[o_1, o_2]
            ),
        )
        _record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(λ[i, j] for (i, j) in has_edge_across_stencil) ≤
                w[o_1, o_2]
            ),
        )
    end
    return nothing
end
