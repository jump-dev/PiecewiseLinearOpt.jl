const BivariateSOS2Method = Union{K1, OptimalTriangleSelection, NineStencil, SixStencil, UnionJack}

function formulate_pwl!(model::JuMP.Model, input_vars::NTuple{2, VarOrAff}, output_vars::NTuple{F, VarOrAff}, pwl::BivariatePWLFunction{F}, method::BivariateSOS2Method, direction::DIRECTION) where {F}
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

    formulate_sos2!(model, T_x, axis_method(method))
    formulate_sos2!(model, T_y, axis_method(method))

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

    formulate_triangle_selection!(model, λ, triangle_direction, method)
end
