# TODO: Generalize to multivariate case.
struct OptimalIndendentBranching <: Method
    sub_solver
end

function formulate_pwl!(model::JuMP.Model, input_vars::NTuple{2, VarOrAff}, output_vars::NTuple{F, VarOrAff}, pwl::BivariatePWLFunction{F}, method::OptimalIndendentBranching, direction::DIRECTION) where {F}
    initPWL!(model)
    counter = model.ext[:PWL].counter
    counter += 1
    model.ext[:PWL].counter = counter

    grid = _continuous_gridpoints_or_die(pwl)

    xs, ys = grid.input_vals, grid.output_vals
    u_1 = unique(x[1] for x in xs)
    u_2 = unique(x[2] for x in xs)
    @assert issorted(u_1)
    @assert issorted(u_2)

    n_1, n_2 = size(xs)

    x_1_to_i = Dict(u_1[i] => i for i in 1:n_1)
    x_2_to_j = Dict(u_2[j] => j for j in 1:n_2)

    λ = JuMP.@variable(model, [1:n_1,1:n_2], lower_bound=0, upper_bound=1, base_name="λ_$counter")
    JuMP.@constraint(model, sum(λ) == 1)
    JuMP.@constraint(model, sum(λ[i,j]*u_1[i] for i in 1:n_1, j in 1:n_2) == input_vars[1])
    JuMP.@constraint(model, sum(λ[i,j]*u_2[j] for i in 1:n_1, j in 1:n_2) == input_vars[2])
    for k in 1:F
        rhs = sum(λ[i,j]*ys[i,j][k] for i in 1:n_1, j in 1:n_2)
        _constrain_output_var(model, output_vars[k], rhs, direction)
    end

    canonical_input_segments = _canonicalize_triangulation(pwl, grid)
    # triangle_direction[i,j] = true means that triangle is of the form
    # -----
    # | \ |
    # -----
    # If false, it is of the form
    # -----
    # | / |
    # -----
    triangle_direction = fill(false, (n_1 - 1, n_2 - 1))
    for input_seg in canonical_input_segments
        # TODO: Remove; assertion redundant with _check_triangulation above.
        @assert length(input_seg) == 3
        i_min, i_max = extrema([input_seg[1][1], input_seg[2][1], input_seg[3][1]])
        j_min, j_max = extrema([input_seg[1][2], input_seg[2][2], input_seg[3][2]])
        if ((i_min, j_max) in input_seg) && ((i_max, j_min) in input_seg)
            triangle_direction[i_min, j_min] = true
        else
            @assert (i_min, j_min) in input_seg
            @assert (i_max, j_max) in input_seg
        end
    end

    n_1, n_2 = size(λ)
    J = Set((i, j) for i in 1:n_1, j in 1:n_2)
    x_val = Matrix{Float64}(undef, n_1, n_2)
    y_val = Matrix{Float64}(undef, n_1, n_2)

    t = ceil(Int, log2(2 * (n_1 - 1) * (n_2 - 1)))
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

        @show J
        for r in J, s in J
            # lexicographic ordering on points on grid
            if r[1] > s[1] || (r[1] == s[1] && r[2] ≥ s[2])
                continue
            end
            edge_is_far_away = norm(r .- s, Inf) > 1
            edge_is_diagonal = !edge_is_far_away && (abs(r[1] - s[1]) == abs(r[2] - s[2]) == 1)
            if edge_is_far_away
                JuMP.@constraint(sub_model, sum(z[j, r, s] for j in 1:t) >= 1)
            elseif edge_is_diagonal
                edge_is_sw_to_ne = edge_is_diagonal && sign(r[1] - s[1]) == sign(r[1] - s[2])
                subrect_sw_corner = min.(r, s)
                triangle_cuts_se_to_nw = triangle_direction[subrect_sw_corner...]
                if (edge_is_sw_to_ne && triangle_cuts_se_to_nw) ||
                   (edge_is_diagonal && !edge_is_sw_to_ne && !triangle_cuts_se_to_nw)
                   JuMP.@constraint(sub_model, sum(z[j, r, s] for j in 1:t) >= 1)
                else
                    JuMP.@constraint(sub_model, sum(z[j, r, s] for j in 1:t) == 0)
                end
            else
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
