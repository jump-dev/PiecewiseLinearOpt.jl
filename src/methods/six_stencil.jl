struct SixStencil <: Method
    axis_method::Method
end

# Logarithmic is default method for x_1 and x_2 axis sos2 constraints.
SixStencil() = SixStencil(Logarithmic())

# TODO: Unit tests for biclique cover
function formulate_pwl!(model::JuMP.Model, input_vars::NTuple{2, VarOrAff}, output_vars::NTuple{F,VarOrAff}, pwl::BivariatePWLFunction{F}, method::SixStencil, direction::DIRECTION) where {F}
    # Verify we have a triangulation
    _check_triangulation(pwl)

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

    # formulations with SOS2 along each dimension
    T_x = [sum(λ[i,j] for i in 1:n_1) for j in 1:n_2]
    T_y = [sum(λ[i,j] for j in 1:n_2) for i in 1:n_1]

    formulate_sos2!(model, T_x, method.axis_method)
    formulate_sos2!(model, T_y, method.axis_method)

    # Eⁿᵉ[i,j] = true means that we must cover the edge {(i,j),(i+1,j+1)}.
    # Otherwise, we must cover the edge {(i+1,j),(i,j+1)}.
    E_sw_ne = fill(false, n_1 - 1, n_2 - 1)
    for segment in pwl.segments
        # TODO: Remove; assertion redundant with _check_triangulation above.
        @assert length(segment.input_vals) == length(segment.output_vals) == 3
        i_a = x_1_to_i[segment.input_vals[1][1]]
        i_b = x_1_to_i[segment.input_vals[2][1]]
        i_c = x_1_to_i[segment.input_vals[3][1]]
        j_a = x_2_to_j[segment.input_vals[1][2]]
        j_b = x_2_to_j[segment.input_vals[2][2]]
        j_c = x_2_to_j[segment.input_vals[3][2]]
        IJ = Set([(i_a, j_a), (i_a, j_a), (i_c, j_c)])
        i_min, i_max = extrema([i_a, i_b, i_c])
        j_min, j_max = extrema([j_a, j_b, j_c])
        if ((i_min, j_max) in IJ) && ((i_max, j_min) in IJ)
            E_sw_ne[i_min, j_min] = true
        else
            @assert (i_min, j_min) in IJ
            @assert (i_max, j_max) in IJ
        end
    end

    # diagonal lines running from SW to NE. Grouped with an offset of 3.
    w_sw_ne = JuMP.@variable(model, [0:2], Bin, base_name="w_sw_ne_$counter")
    for o in 0:2
        A_o = Set{Tuple{Int,Int}}()
        B_o = Set{Tuple{Int,Int}}()
        for off_1 in o:3:(n_1 - 2)
            sw_in_A = true # whether we put the SW corner of the next triangle to cover in set A
            for i in (1 + off_1):(n_1 - 1)
                j = i - off_1
                if !(1 ≤ i ≤ n_1 - 1)
                    continue
                end
                if !(1 ≤ j ≤ n_2-1)
                    continue # should never happen
                end
                if E_sw_ne[i,j] # if we need to cover the edge...
                    if sw_in_A # figure out which set we need to put it in; this depends on previous triangle in our current line
                        push!(A_o, (i  ,j  ))
                        push!(B_o, (i+1,j+1))
                    else
                        push!(A_o, (i+1,j+1))
                        push!(B_o, (i  ,j  ))
                    end
                    sw_in_A = !sw_in_A
                end
            end
        end
        for off_2 in (3-o):3:(n_2-1)
            sw_in_A = true
            for j in (off_2+1):(n_2-1)
                i = j - off_2
                if !(1 ≤ i ≤ n_1-1)
                    continue
                end
                if E_sw_ne[i,j]
                    if sw_in_A
                        push!(A_o, (i  ,j  ))
                        push!(B_o, (i+1,j+1))
                    else
                        push!(A_o, (i+1,j+1))
                        push!(B_o, (i  ,j  ))
                    end
                    sw_in_A = !sw_in_A
                end
            end
        end
        JuMP.@constraints(model, begin
            sum(λ[i,j] for (i,j) in A_o) ≤     w_sw_ne[o]
            sum(λ[i,j] for (i,j) in B_o) ≤ 1 - w_sw_ne[o]
        end)
    end

    w_se_nw = JuMP.@variable(model, [0:2], Bin, base_name="w_se_nw_$counter")
    for o in 0:2
        A_o = Set{Tuple{Int,Int}}()
        B_o = Set{Tuple{Int,Int}}()
        for off_1 in o:3:(n_1 - 2)
            se_in_A = true
            for j in 1:(n_2-1)
                i = n_1 - j - off_1
                if !(1 ≤ i ≤ n_1-1)
                    continue
                end
                if !E_sw_ne[i,j]
                    if se_in_A
                        push!(A_o, (i+1,j  ))
                        push!(B_o, (i  ,j+1))
                    else
                        push!(A_o, (i  ,j+1))
                        push!(B_o, (i+1,j  ))
                    end
                    se_in_A = !se_in_A
                end
            end
        end
        for off_2 in (3-o):3:(n_2-1)
            se_in_A = true
            for j in (off_2+1):(n_2-1)
                i = n_1 - j + off_2
                if !(1 ≤ i ≤ n_1-1)
                    continue
                end
                if !E_sw_ne[i,j]
                    if se_in_A
                        push!(A_o, (i+1,j  ))
                        push!(B_o, (i  ,j+1))
                    else
                        push!(A_o, (i  ,j+1))
                        push!(B_o, (i+1,j  ))
                    end
                    se_in_A = !se_in_A
                end
            end
        end
        JuMP.@constraints(model, begin
            sum(λ[i,j] for (i,j) in A_o) ≤     w_se_nw[o]
            sum(λ[i,j] for (i,j) in B_o) ≤ 1 - w_se_nw[o]
        end)
    end
    return
end
