# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

"""
    constrain_output!(formulation::PWLFormulation, y, value, direction::Direction)

Adds the constraint that the output `y` equals, is at least or is at most the
expression `value` of the function value, as given by `direction`.
"""
function constrain_output!(
    formulation::PWLFormulation,
    output_var::Scalar,
    func_val::Scalar,
    direction::Direction,
)
    record!(
        formulation,
        _output_constraint(formulation.model, output_var, func_val, direction),
    )
    return
end

function _output_constraint(model, y, value, ::Graph)
    return JuMP.@constraint(model, y == value)
end

function _output_constraint(model, y, value, ::Epigraph)
    return JuMP.@constraint(model, y >= value)
end

function _output_constraint(model, y, value, ::Hypograph)
    return JuMP.@constraint(model, y <= value)
end

struct _Grid{D,F}
    input_vals::Array{NTuple{D,Float64},D}
    output_vals::Array{NTuple{F,Float64},D}

    function _Grid(
        input_vals::Array{NTuple{D,Float64},D},
        output_vals::Array{NTuple{F,Float64},D},
    ) where {D,F}
        if size(input_vals) != size(output_vals)
            error(
                "Incompatible input/output sizes $(size(input_vals)) and $(size(output_vals)). Must be the same.",
            )
        end
        return new{D,F}(input_vals, output_vals)
    end
end

function _continuous_gridpoints_or_die(pwl::PWLFunction{D,F}) where {D,F}
    # Step 1: Collect all input points
    xs = Set{NTuple{D,Float64}}()
    for segment in pwl.segments
        for val in segment.input_vals
            push!(xs, val)
        end
    end

    # Step 2: Verify that the function is single-valued and continuous
    ys = Dict{NTuple{D,Float64},NTuple{F,Float64}}()
    for segment in pwl.segments
        _check_single_valued(segment)
        for (x, y) in zip(segment.input_vals, segment.output_vals)
            if haskey(ys, x)
                if y != ys[x]
                    error(
                        "Expected piecewise linear function to be continuous.",
                    )
                end
            else
                ys[x] = y
            end
        end
    end

    @assert length(xs) == length(ys)

    # Step 3: Build grid, and verify every point is included in PWL function
    axes = [sort!(unique(x[i] for x in xs)) for i in 1:D]

    x_grid = collect(Base.Iterators.product(axes...))
    y_grid = map(x -> get(ys, x, nothing), x_grid)

    if length(x_grid) != length(xs)
        error("Decomposition is not aligned along a grid.")
    end

    # Step 4: Verify that the segments cover the grid. Only univariate and
    # bivariate functions are formulated on grids.
    if D == 1
        # Every pair of neighboring breakpoints must form exactly one segment
        index = Dict(x[1] => i for (i, x) in enumerate(x_grid))
        covered = fill(false, length(x_grid) - 1)
        for segment in pwl.segments
            if length(segment.input_vals) != 2
                error(
                    "Every segment of a univariate piecewise linear function must have two points.",
                )
            end
            i, j = minmax((index[v[1]] for v in segment.input_vals)...)
            if j != i + 1 || covered[i]
                error(
                    "The segments of a univariate piecewise linear function must join neighboring breakpoints without overlapping, but a segment has the points $(segment.input_vals).",
                )
            end
            covered[i] = true
        end
        if !all(covered)
            error(
                "Univariate piecewise linear function does not have a connected domain.",
            )
        end
    else
        @assert D == 2
        # Each cell is checked to be split into two triangles by
        # _triangle_directions
        if length(pwl.segments) !=
           2 * (size(x_grid, 1) - 1) * (size(x_grid, 2) - 1)
            error(
                "Bivariate piecewise linear function domain is not rectangular.",
            )
        end
    end

    return _Grid(x_grid, y_grid)
end

function _canonicalize_triangulation(
    pwl::PWLFunction{D,F},
    grid::_Grid{D,F},
) where {D,F}
    _check_triangulation(pwl)

    xs = grid.input_vals
    U = [unique(x[j] for x in xs) for j in 1:D]
    for j in 1:D
        @assert issorted(U[j])
    end
    x_to_i = [Dict(U[j][i] => i for i in 1:size(xs, j)) for j in 1:D]

    canonical_input_segments = Vector{NTuple{D,Int64}}[]
    for segment in pwl.segments
        push!(
            canonical_input_segments,
            [ntuple(j -> x_to_i[j][v[j]], D) for v in segment.input_vals],
        )
    end
    return canonical_input_segments
end

function _create_convex_multiplier_vars(
    formulation::PWLFormulation,
    grid::_Grid{D,F},
    input_vars::NTuple{D,Scalar},
    output_vars::NTuple{F,Scalar},
    direction::Direction,
) where {D,F}
    model = formulation.model
    λ = similar(Array{JuMP.variable_ref_type(model)}, axes(grid.input_vals))
    for I in CartesianIndices(λ)
        λ[I] = JuMP.@variable(
            model,
            lower_bound = 0,
            upper_bound = 1,
            base_name =
                variable_name(formulation, "λ[$(join(Tuple(I), ","))]")
        )
    end
    record!(formulation, λ)

    record!(formulation, JuMP.@constraint(model, sum(λ) == 1))
    for j in 1:D
        record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(λ[I] * grid.input_vals[I][j] for I in eachindex(λ)) ==
                input_vars[j]
            ),
        )
    end
    for j in 1:F
        constrain_output!(
            formulation,
            output_vars[j],
            sum(λ[I] * grid.output_vals[I][j] for I in eachindex(λ)),
            direction,
        )
    end
    return λ
end

_multiplier(λ::AbstractArray, k::Int) = λ[k]
_multiplier(λ::AbstractArray, k::Tuple) = λ[k...]

"""
    add_biclique_cover!(formulation::PWLFormulation, λ::AbstractArray, cover)

Adds the independent branching formulation given by a biclique cover.

For each level `(L, R)` in `cover`, it adds a binary variable `z` and the
constraints `sum(λ[L]) ≤ z` and `sum(λ[R]) ≤ 1 - z`, so that the multipliers
in `L` and those in `R` cannot be nonzero at the same time (Huchette and
Vielma, 2019). The levels are pairs of vectors of indices into `λ`, either
integers or tuples for a matrix. A level with an empty side separates nothing
and is left out.

The formulation is valid for a disjunctive constraint on the support of `λ` if
no level separates two multipliers that may be nonzero together, and every
pair that may not is separated by some level.
"""
function add_biclique_cover!(
    formulation::PWLFormulation,
    λ::AbstractArray,
    cover::Vector{<:Tuple{AbstractVector,AbstractVector}},
)
    levels = [(L, R) for (L, R) in cover if !isempty(L) && !isempty(R)]
    if isempty(levels)
        return nothing
    end
    model = formulation.model
    z = record!(
        formulation,
        JuMP.@variable(
            model,
            [1:length(levels)],
            Bin,
            base_name = variable_name(formulation, "z")
        ),
    )
    for (j, (L, R)) in enumerate(levels)
        record!(
            formulation,
            JuMP.@constraint(model, sum(_multiplier(λ, k) for k in L) ≤ z[j]),
        )
        record!(
            formulation,
            JuMP.@constraint(
                model,
                sum(_multiplier(λ, k) for k in R) ≤ 1 - z[j]
            ),
        )
    end
    return nothing
end

# Errors if a segment has several values at the same point, i.e. is a piece of
# a relaxation rather than of a function
function _check_single_valued(segment::SegmentPointRep)
    for (i, x) in enumerate(segment.input_vals), j in 1:(i-1)
        if x == segment.input_vals[j] &&
           segment.output_vals[i] != segment.output_vals[j]
            error(
                "A segment has several values at the point $x, so it is a piece of a relaxation rather than of a function. Relaxations are formulated by IndependentBranching and DisaggregatedLogarithmic.",
            )
        end
    end
    return
end

function _check_triangle(segment::SegmentPointRep{D}) where {D}
    if !(length(segment.input_vals) == length(segment.output_vals) == D + 1)
        error(
            "Encountered something that was not a simplex while expecting a trianglulation.",
        )
    end
    return
end

function _check_triangulation(pwl::PWLFunction)
    for segment in pwl.segments
        _check_triangle(segment)
    end
    return
end
