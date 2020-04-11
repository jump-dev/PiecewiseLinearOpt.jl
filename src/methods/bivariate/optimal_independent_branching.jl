struct OptimalIndendentBranching <: Method
    sub_solver
    axis_method::Method
end

OptimalIndendentBranching(sub_solver) = OptimalIndendentBranching(, sub_solver, Logarithmic())

axis_method(method::OptimalIndendentBranching) = method.axis_method

function formulate_triangle_selection!(model::JuMP.Model, λ::Matrix{JuMP.VariableRef}, triangle_direction::Matrix{Bool}, method::OptimalIndendentBranching)
    n_1, n_2 = size(λ)
    # Number of segments
    t = ceil(Int, log2(2 * (n_1 - 1) * (n_2 - 1)))
    # Number of breakpoints
    m = n_1 * n_2
    P = Set((r, s) for r in 1:m, s in 1:m if r < s)
    x_val = Matrix{Float64}()
    y_val = Matrix{Float64}()
    while true
        model = JuMP.Model(solver = method.sub_solver)
        JuMP.@variable(model, x[1:t, P], Bin)
        JuMP.@variable(model, y[1:t, P], Bin)
        JuMP.@variable(model, z[1:t, P, P], Bin)

        for j in 1:t
            for (i, j)(r, s) in P
                JuMP.@constraints(model, begin
                    z[j, r, s] ≤ x[j, r] + x[j, s]
                    z[j, r, s] ≤ x[j, r] + y[j, r]
                    z[j, r, s] ≤ x[j, s] + y[j, s]
                    z[j, r, s] ≤ y[j, r] + y[j, s]
                    z[j, r, s] ≥ x[j, r] + y[j, s] - 1
                    z[j, r, s] ≥ x[j, s] + y[j, r] - 1
                end)
            end
            for r in 1:m
                JuMP.@constraint(model, x[j, r] + y[j, r] ≤ 1)
            end
        end

        for (r, s) in P
            if E[r, s]
                JuMP.@constraint(model, sum(z[j, r, s] for j in 1:t) == 0)
            else
                JuMP.@constraint(model, sum(z[j, r, s] for j in 1:t) ≥ 1)
            end
        end

        JuMP.@objective(model, Min, sum(x) + sum(y))
        JuMP.optimize!(model)
        if JuMP.solution_status == MOI.FEASIBLE_POINT
            x_val = JuMP.value.(x)
            y_val = JuMP.value.(y)
            break
        else
            t += 1
        end
    end

    w = JuMP.@variable(m, [1:t], Bin)
    for k in 1:t
        JuMP.@constraints(m, begin
            sum(λ[i, j] for (i, j) in E if x_val[k, i, j] == 1) ≤ w[k]
            sum(λ[i, j] for (i, j) in E if y_val[k, i, j] == 1) ≤ 1 - w[k]
        end)
    end
    return nothing
end
