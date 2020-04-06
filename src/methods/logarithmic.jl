struct Logarithmic <: Method end

function _reflected_gray_codes(k::Int)
    if k == 0
        return Vector{Float64}[]
    elseif k == 1
        return [[0.0], [1.0]]
    else
        codes′ = _reflected_gray_codes(k-1)
        return vcat([vcat(code, 0.0) for code in codes′],
                    [vcat(code, 1.0) for code in reverse(codes′)])
    end
end

function formulate_pwl!(model::JuMP.Model, input_vars::Tuple{VarOrAff}, output_vars::NTuple{F,VarOrAff}, pwl::UnivariatePWLFunction{F}, method::Logarithmic, direction::DIRECTION) where {F}
    grid = _continuous_gridpoints_or_die(pwl)

    λ = _create_convex_multiplier_vars(model, grid, input_vars, output_vars, direction)
    formulate_sos2!(model, λ, method)
    return
end

function formulate_sos2!(model::JuMP.Model, λ::Vector{T}, method::Logarithmic) where {T <: VarOrAff}
    counter = model.ext[:PWL].counter
    n = length(λ)
    d = n - 1
    k = ceil(Int, log2(d))
    y = JuMP.@variable(model, [1:k], Bin, base_name="y_$counter")
    _sos2_encoding_constraints!(model, λ, y, _reflected_gray_codes(k), _unit_vector_hyperplanes(k))
    return nothing
end
