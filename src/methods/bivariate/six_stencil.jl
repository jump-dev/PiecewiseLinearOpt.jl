struct SixStencil <: Method
    axis_method::Method
end
SixStencil() = SixStencil(Logarithmic())

axis_method(method::SixStencil) = method.axis_method

# TODO: Unit tests for biclique cover
function formulate_triangle_selection!(model::JuMP.Model, λ::Matrix{JuMP.VariableRef}, triangle_direction::Matrix{Bool}, method::SixStencil)
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)
    counter = model.ext[:PWL].counter

    # diagonal lines running from SW to NE. Grouped with an offset of 3.
    w_sw_ne = JuMP.@variable(model, [0:2], Bin, base_name="w_sw_ne_$counter")
    for o in 0:2
        A_o = Set{Tuple{Int, Int}}()
        B_o = Set{Tuple{Int, Int}}()
        for off_1 in o:3:(n_1 - 2)
            sw_in_A = true # whether we put the SW corner of the next triangle to cover in set A
            for i in (1 + off_1):(n_1 - 1)
                j = i - off_1
                if !(1 ≤ i ≤ n_1 - 1)
                    continue
                end
                if !(1 ≤ j ≤ n_2 - 1)
                    continue # should never happen
                end
                if triangle_direction[i, j] # if we need to cover the edge...
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
                if triangle_direction[i,j]
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
                if !triangle_direction[i,j]
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
                if !triangle_direction[i,j]
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
    return nothing
end
