struct K1 <: Method
    axis_method::Method
end
K1() = K1(Logarithmic())

axis_method(method::K1) = method.axis_method

function formulate_triangle_selection!(model::JuMP.Model, λ::Matrix{JuMP.VariableRef}, triangle_direction::Matrix{Bool}, method::K1)
    n_1, n_2 = size(λ)
    @assert size(triangle_direction) == (n_1 - 1, n_2 - 1)
    counter = model.ext[:PWL].counter

    # TODO: Certify triangulation is "uniform"

    w = JuMP.@variable(model, [1:2], Bin, base_name="w_$counter")
    JuMP.@constraints(model, begin
        sum(λ[i, j] for i in 1:n_1, j in 1:n_2 if mod(i, 2) == mod(j, 2) && (i + j) in 2:4:(n_1 + n_2)) ≤ w[1]
        sum(λ[i, j] for i in 1:n_1, j in 1:n_2 if mod(i, 2) == mod(j, 2) && (i + j) in 4:4:(n_1 + n_2)) ≤ 1 - w[1]
        sum(λ[i, j] for i in 1:n_1, j in 1:n_2 if mod(i, 2) != mod(j, 2) && (i + j) in 3:4:(n_1 + n_2)) ≤ w[2]
        sum(λ[i, j] for i in 1:n_1, j in 1:n_2 if mod(i, 2) != mod(j, 2) && (i + j) in 5:4:(n_1 + n_2)) ≤ 1 - w[2]
    end)

    return nothing
end
