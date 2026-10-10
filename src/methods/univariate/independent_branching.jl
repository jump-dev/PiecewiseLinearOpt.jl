# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Formulations for generalized 1D-ordered combinatorial disjunctive constraints
# (CDCs) from
#
#   B. Lyu, I. V. Hicks and J. Huchette. Building formulations for piecewise
#   linear relaxations of nonlinear functions. arXiv:2304.14542, 2023.
#
# A CDC is given by sets S^1, ..., S^d of indices into a vector of convex
# multipliers λ, and requires the support of λ to lie within one of them. It is
# generalized 1D-ordered if S^i and S^j are disjoint whenever |i - j| ≥ 2. The
# SOS2 constraint is one example, and so is the graph of a univariate piecewise
# linear function, which may be discontinuous, when each segment is given by
# its two end points.

"""
    IndependentBranching(scheme::Symbol = :biclique)

Logarithmically sized ideal formulation of univariate piecewise linear
functions, based on independent branching over the ordered segments (Lyu, Hicks
and Huchette, 2023). It supports discontinuous functions, and can be used as the
axis method of the bivariate grid formulations.

The `scheme` gives the construction of the branching levels:

- `:biclique`: biclique cover from a balanced reversed edge ranking
  (Theorem 4).
- `:balanced`: balanced Gray code (Theorem 1 and Remark 4).
- `:brgc`: binary reflected Gray code (Theorem 1), equivalent to
  [`LogarithmicIndependentBranching`](@ref) for continuous functions.
"""
struct IndependentBranching <: Method
    scheme::Symbol

    function IndependentBranching(scheme::Symbol = :biclique)
        if !(scheme in (:biclique, :balanced, :brgc))
            error(
                "Unknown independent branching scheme $(repr(scheme)). Supported schemes are :biclique, :balanced and :brgc.",
            )
        end
        return new(scheme)
    end
end

function formulate_sos2!(
    model::JuMP.Model,
    λ::Vector{T},
    method::IndependentBranching,
) where {T<:VarOrAff}
    S = [[i, i + 1] for i in 1:length(λ)-1]
    _add_biclique_cover!(model, λ, _ordered_cdc_cover(S, method.scheme))
    return nothing
end

function formulate_pwl!(
    model::JuMP.Model,
    input_vars::Tuple{VarOrAff},
    output_vars::NTuple{F,VarOrAff},
    pwl::PWLFunctionPointRep{1,F},
    method::IndependentBranching,
    direction::DIRECTION,
) where {F}
    counter = model.ext[:PWL].counter

    # The multipliers are indexed by the distinct (input, output) points, so
    # a jump in the function gives two points at the same input value.
    segments = sort(pwl.segments; by = seg -> minimum(seg.input_vals))
    point_index = Dict{Tuple{NTuple{1,Float64},NTuple{F,Float64}},Int}()
    points = Tuple{NTuple{1,Float64},NTuple{F,Float64}}[]
    S = Vector{Int}[]
    for seg in segments
        push!(
            S,
            [
                get!(point_index, p) do
                    push!(points, p)
                    return length(points)
                end for p in zip(seg.input_vals, seg.output_vals)
            ],
        )
    end
    _check_generalized_1d_ordered(S)

    n = length(points)
    λ = JuMP.@variable(
        model,
        [1:n],
        lower_bound = 0,
        upper_bound = 1,
        base_name = "λ_$counter"
    )
    JuMP.@constraint(model, sum(λ) == 1)
    JuMP.@constraint(
        model,
        input_vars[1] == sum(λ[v] * points[v][1][1] for v in 1:n)
    )
    for k in 1:F
        rhs = sum(λ[v] * points[v][2][k] for v in 1:n)
        _constrain_output_var(model, output_vars[k], rhs, direction)
    end
    _add_biclique_cover!(model, λ, _ordered_cdc_cover(S, method.scheme))
    return nothing
end

function _check_generalized_1d_ordered(S::Vector{Vector{Int}})
    for i in 1:length(S), j in i+2:length(S)
        if !isdisjoint(S[i], S[j])
            error(
                "The segments must be ordered so that only neighboring segments share points, but segments $i and $j share a point.",
            )
        end
    end
    return
end

# Returns the levels (L^j, R^j) of an ideal independent branching formulation
# of the generalized 1D-ordered CDC given by S.
function _ordered_cdc_cover(S::Vector{Vector{Int}}, scheme::Symbol)
    d = length(S)
    if d ≤ 1
        return Tuple{Vector{Int},Vector{Int}}[]
    end
    if scheme == :brgc
        t = ceil(Int, log2(d))
        return _gray_code_cover(S, _reflected_gray_codes(t)[1:d])
    end
    φ = _balanced_reversed_edge_ranking(d)
    if scheme == :balanced
        return _gray_code_cover(S, _gray_codes_from_ranking(φ))
    end
    @assert scheme == :biclique
    return _biclique_cover_from_ranking(S, φ)
end

# Theorem 1: given a Gray code h^1, ..., h^d, level j separates the segments
# whose code has a zero in position j from those with a one.
function _gray_code_cover(
    S::Vector{Vector{Int}},
    codes::Vector{<:AbstractVector},
)
    cover = Tuple{Vector{Int},Vector{Int}}[]
    for j in 1:length(codes[1])
        L = reduce(
            union,
            (S[i] for i in 1:length(S) if codes[i][j] == 0);
            init = Int[],
        )
        R = reduce(
            union,
            (S[i] for i in 1:length(S) if codes[i][j] == 1);
            init = Int[],
        )
        push!(cover, (setdiff(L, R), setdiff(R, L)))
    end
    return cover
end

# Algorithm 1 with balanced cuts (Remark 4): a reversed edge ranking φ of the
# path S^1 - S^2 - ... - S^d, where φ[i] is the rank of the edge between S^i
# and S^(i+1). It uses ceil(log2(d)) ranks.
function _balanced_reversed_edge_ranking(d::Int)
    φ = zeros(Int, d - 1)
    function label!(first, last, level)
        n = last - first + 1
        if n ≤ 1
            return
        end
        cut = first + div(n, 2) - 1
        φ[cut] = level
        label!(first, cut, level + 1)
        label!(cut + 1, last, level + 1)
        return
    end
    label!(1, d, 1)
    return φ
end

# Theorem 2: walking along the path, the code flips bit φ(e) when crossing the
# edge e.
function _gray_codes_from_ranking(φ::Vector{Int})
    r = maximum(φ)
    codes = [zeros(Int, r)]
    for rank in φ
        code = copy(codes[end])
        code[rank] = 1 - code[rank]
        push!(codes, code)
    end
    return codes
end

# Algorithms 2 and 3 with Theorem 4: cut the path recursively at the edge of
# lowest rank, and combine the cuts of each rank into one level, alternating
# the sides so that neighboring cuts of the same rank do not overlap.
function _biclique_cover_from_ranking(S::Vector{Vector{Int}}, φ::Vector{Int})
    # cuts[e] = (segments before edge e, segments after) within the subpath
    # that edge e was used to separate.
    cuts = Dict{Int,Tuple{UnitRange{Int},UnitRange{Int}}}()
    function separate!(first, last)
        if last ≤ first
            return
        end
        e = first - 1 + argmin(view(φ, first:last-1))
        cuts[e] = (first:e, e+1:last)
        separate!(first, e)
        separate!(e + 1, last)
        return
    end
    separate!(1, length(S))

    cover = Tuple{Vector{Int},Vector{Int}}[]
    for rank in 1:maximum(φ)
        A, B = Int[], Int[]
        for (count, e) in enumerate(e for e in 1:length(φ) if φ[e] == rank)
            I, J = cuts[e]
            if isodd(count)
                append!(A, I)
                append!(B, J)
            else
                append!(A, J)
                append!(B, I)
            end
        end
        L = reduce(union, (S[i] for i in A); init = Int[])
        R = reduce(union, (S[i] for i in B); init = Int[])
        push!(cover, (setdiff(L, R), setdiff(R, L)))
    end
    return cover
end
