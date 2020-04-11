# TODO: Generalize to multivariate case.
struct OptimalTriangleSelection <: Method
    sub_solver
    axis_method::Method
end
OptimalTriangleSelection(sub_solver) = OptimalTriangleSelection(sub_solver, Logarithmic())

axis_method(method::OptimalTriangleSelection) = method.axis_method

function formulate_triangle_selection!(model::JuMP.Model, λ::Matrix{JuMP.VariableRef}, triangle_direction::Matrix{Bool}, method::OptimalTriangleSelection)
    counter = model.ext[:PWL].counter
    n_1, n_2 = size(λ)
    J = Set((i, j) for i in 1:n_1, j in 1:n_2)
    x_val = Matrix{Float64}(undef, n_1, n_2)
    y_val = Matrix{Float64}(undef, n_1, n_2)

    t = 1
    while true
        @show method.sub_solver
        sub_model = JuMP.Model(method.sub_solver)
        JuMP.@variable(sub_model, x[1:t, J], Bin)
        JuMP.@variable(sub_model, y[1:t, J], Bin)
        JuMP.@variable(sub_model, z[1:t, J, J], Bin)
        for j in 1:t
            for r in J, s in J
                # lexicographic ordering on points on grid
                if r[1] > s[1] || (r[1] == s[1] && r[2] ≥ s[2])
                    continue
                end
                JuMP.@constraints(sub_model, begin
                    z[j, r, s] ≤ x[j, r] + x[j, s]
                    z[j, r, s] ≤ x[j, r] + y[j, r]
                    z[j, r, s] ≤ x[j, s] + y[j, s]
                    z[j, r, s] ≤ y[j, r] + y[j, s]
                    z[j, r, s] ≥ x[j, r] + y[j, s] - 1
                    z[j, r, s] ≥ x[j, s] + y[j, r] - 1
                end)
            end
            for r in J
                JuMP.@constraint(sub_model, x[j, r] + y[j, r] ≤ 1)
            end
        end

        for r in J, s in J
            # lexicographic ordering on points on grid
            if r[1] > s[1] || (r[1] == s[1] && r[2] ≥ s[2])
                continue
            end
            edge_is_far_away = norm(r .- s, Inf) > 1
            edge_is_diagonal = !edge_is_far_away && (abs(r[1] - s[1]) == abs(r[2] - s[2]) == 1)
            if edge_is_diagonal
                edge_is_sw_to_ne = edge_is_diagonal && sign(r[1] - s[1]) == sign(r[1] - s[2])
                subrect_sw_corner = min.(r, s)
                triangle_cuts_se_to_nw = triangle_direction[subrect_sw_corner...]
                if (edge_is_sw_to_ne && triangle_cuts_se_to_nw) ||
                   (edge_is_diagonal && !edge_is_sw_to_ne && !triangle_cuts_se_to_nw)
                   JuMP.@constraint(sub_model, sum(z[j, r, s] for j in 1:t) >= 1)
                else
                    JuMP.@constraint(sub_model, sum(z[j, r, s] for j in 1:t) == 0)
                end
            elseif !edge_is_far_away
                @assert r[1] == s[1] || r[2] == s[2]
                JuMP.@constraint(sub_model, sum(z[j, r, s] for j in 1:t) == 0)
            end
        end

        JuMP.@objective(sub_model, Min, sum(x) + sum(y))
        JuMP.optimize!(sub_model)
        if JuMP.primal_status(sub_model) == MOI.FEASIBLE_POINT
            x_val = JuMP.value.(x)
            y_val = JuMP.value.(y)
            break
        else
            t += 1
        end
    end
    z = JuMP.@variable(model, [1:t], Bin, base_name = "z_$counter")

    for k in 1:t
        JuMP.@constraints(model, begin
            sum(λ[i, j] for (i, j) in J if x_val[k, (i, j)] ≈ 1) ≤ z[k]
            sum(λ[i, j] for (i, j) in J if y_val[k, (i, j)] ≈ 1) ≤ 1 - z[k]
        end)
    end
    return nothing
end
