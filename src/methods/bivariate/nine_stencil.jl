struct NineStencil <: Method
    axis_method::Method
end
NineStencil() = NineStencil(Logarithmic())

axis_method(method::NineStencil) = method.axis_method

# TODO: Unit tests for biclique cover
function formulate_triangle_selection!(model::JuMP.Model, λ::Matrix{JuMP.VariableRef}, triangle_direction::Matrix{Bool}, method::NineStencil)
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)
    counter = model.ext[:PWL].counter

    w = JuMP.@variable(model, [1:3, 1:3], Bin, base_name = "w_$counter")
    for o_1 in 1:3, o_2 in 1:3
        has_edge_across_stencil = Set{Tuple{Int, Int}}()
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
                    if !triangle_direction[i, j - 1]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
            let k = i - 1, l = j - 1
                if (1 ≤ k ≤ n_1) && (1 ≤ l ≤ n_2)
                    if triangle_direction[i - 1, j - 1]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
            let k = i - 1, l = j + 1
                if (1 ≤ k ≤ n_1) && (1 ≤ l ≤ n_2)
                    if !triangle_direction[i - 1, j]
                        push!(has_edge_across_stencil, (k, l))
                    end
                end
            end
        end
        JuMP.@constraints(model, begin
            sum(λ[i, j] for i in o_1:3:n_1, j in o_2:3:n_2) ≤  1 - w[o_1, o_2]
            sum(λ[i, j] for (i, j) in has_edge_across_stencil) ≤ w[o_1, o_2]
        end)
    end
    return nothing
end
