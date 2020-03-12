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

function _formulate!(model::JuMP.Model, input_vars::Tuple{VarOrAff}, output_vars::NTuple{F,VarOrAff}, pwl::UnivariatePWLFunction, method::Logarithmic, direction::DIRECTION) where {F}
    xs, ys = _continuous_breakpoints_or_die(pwl)

    n = length(xs)
    d = n - 1
    k = ceil(Int, log2(n))
    counter = model.ext[:PWL].counter
    y = JuMP.@variable(model, [1:k], Bin, base_name="y_$counter")

    λ = _create_convex_multiplier_vars(model, xs, ys, input_vars, output_vars, direction)

    _sos2_encoding_constraints!(model, λ, y, _reflected_gray_codes(k), _unit_vector_hyperplanes(k))
    return
end
