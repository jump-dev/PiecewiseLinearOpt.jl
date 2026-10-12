# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

function formulate_pwl!(
    formulation::PWLFormulation,
    input_vars::NTuple{2,Scalar},
    output_vars::NTuple{F,Scalar},
    pwl::PWLFunction{2,F},
    method::TriangleSelectionMethod,
    direction::Direction,
) where {F}
    λ, triangle_direction = _grid_convex_multipliers!(
        formulation,
        input_vars,
        output_vars,
        pwl,
        direction,
    )

    # SOS2 constraints along each axis select a cell
    n_1, n_2 = size(λ)
    μ_1 = [sum(λ[i, j] for j in 1:n_2) for i in 1:n_1]
    μ_2 = [sum(λ[i, j] for i in 1:n_1) for j in 1:n_2]
    formulate_sos2!(
        _subformulation(formulation, "axis1"),
        μ_1,
        axis_method(method),
    )
    formulate_sos2!(
        _subformulation(formulation, "axis2"),
        μ_2,
        axis_method(method),
    )

    formulate_triangle_selection!(formulation, λ, triangle_direction, method)
    return nothing
end

# Adds a convex multiplier for each point of the grid underlying `pwl`, and
# returns them together with the direction of the diagonal in each cell.
function _grid_convex_multipliers!(
    formulation::PWLFormulation,
    input_vars::NTuple{2,Scalar},
    output_vars::NTuple{F,Scalar},
    pwl::PWLFunction{2,F},
    direction::Direction,
) where {F}
    grid = _continuous_gridpoints_or_die(pwl)
    triangle_direction = _triangle_directions(pwl, grid)
    λ = _create_convex_multiplier_vars(
        formulation,
        grid,
        input_vars,
        output_vars,
        direction,
    )
    return λ, triangle_direction
end

# triangle_direction[i,j] = true means that triangle is of the form
# -----
# | \ |
# -----
# If false, it is of the form
# -----
# | / |
# -----
function _triangle_directions(pwl::PWLFunction{2}, grid::Grid{2})
    n_1, n_2 = size(grid.input_vals)
    triangle_direction = fill(false, (n_1 - 1, n_2 - 1))
    # The first triangle found in each cell, and the number of triangles in it
    first_triangle = Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}()
    num_triangles = zeros(Int, n_1 - 1, n_2 - 1)
    for input_seg in _canonicalize_triangulation(pwl, grid)
        i_min, i_max = extrema(v[1] for v in input_seg)
        j_min, j_max = extrema(v[2] for v in input_seg)
        corners = Set(input_seg)
        if i_max - i_min != 1 || j_max - j_min != 1 || length(corners) != 3
            error(
                "The grid formulations require every triangle to have three corners of a single cell of the grid, but a triangle has the grid points $input_seg.",
            )
        end
        # Three corners of a cell contain exactly one of its diagonals
        nw_se = (i_min, j_max) in corners && (i_max, j_min) in corners
        cell = (i_min, j_min)
        num_triangles[cell...] += 1
        if num_triangles[cell...] == 1
            first_triangle[cell] = corners
            triangle_direction[cell...] = nw_se
        elseif num_triangles[cell...] > 2 ||
               triangle_direction[cell...] != nw_se ||
               length(union(first_triangle[cell], corners)) != 4
            _not_split_along_diagonal(cell)
        end
    end
    for cell in CartesianIndices(num_triangles)
        if num_triangles[cell] != 2
            _not_split_along_diagonal(Tuple(cell))
        end
    end
    return triangle_direction
end

function _not_split_along_diagonal(cell::Tuple{Int,Int})
    return error(
        "The grid formulations require every cell to be split into two triangles along one of its diagonals, but cell $cell is not.",
    )
end

"""
    cell_diagonals(nw_se::Matrix{Bool}, i::Int, j::Int)

Returns the diagonals of cell `(i, j)` of a grid triangulation as
`(edge, non_edge)`, where `edge` is the diagonal that splits the cell into its
two triangles and `non_edge` is the other. Each is a pair of grid indices,
e.g. `((i, j + 1), (i + 1, j))`, into the matrix of multipliers passed to
[`formulate_triangle_selection!`](@ref) together with `nw_se`.

A triangle selection must prevent the two ends of `non_edge` from being
nonzero together.
"""
function cell_diagonals(triangle_direction::Matrix{Bool}, i::Int, j::Int)
    sw, se, nw, ne = (i, j), (i + 1, j), (i, j + 1), (i + 1, j + 1)
    return triangle_direction[i, j] ? ((nw, se), (sw, ne)) :
           ((sw, ne), (nw, se))
end

const _BicliqueCover =
    Vector{Tuple{Vector{Tuple{Int,Int}},Vector{Tuple{Int,Int}}}}

# The optimal covers only depend on the triangulation, and are expensive to
# compute, so they are shared between models.
const _OPTIMAL_COVER_CACHE = Dict{Tuple{Matrix{Bool},Bool},_BicliqueCover}()
const _OPTIMAL_COVER_LOCK = ReentrantLock()

# Finds the smallest number of levels (A_k, B_k) such that two grid points
# that share a triangle are never on opposite sides of a level, while every
# other pair of points is on opposite sides of at least one level. If
# `separate_all` is false, only the pairs that are endpoints of a cell diagonal
# need to be separated, as the rest is handled by SOS2 constraints on each axis.
#
# Returns `nothing` if no cover is found within `time_limit` seconds.
function _optimal_biclique_cover(
    sub_solver,
    triangle_direction::Matrix{Bool},
    separate_all::Bool,
    time_limit::Real,
)
    key = (copy(triangle_direction), separate_all)
    return lock(_OPTIMAL_COVER_LOCK) do
        if haskey(_OPTIMAL_COVER_CACHE, key)
            return _OPTIMAL_COVER_CACHE[key]
        end
        cover = _compute_optimal_biclique_cover(
            sub_solver,
            triangle_direction,
            separate_all,
            time_limit,
        )
        if cover !== nothing
            _OPTIMAL_COVER_CACHE[key] = cover
        end
        return cover
    end
end

function _compute_optimal_biclique_cover(
    sub_solver,
    triangle_direction::Matrix{Bool},
    separate_all::Bool,
    time_limit::Real,
)
    deadline = time() + time_limit
    n_1, n_2 = size(triangle_direction) .+ 1
    J = vec([(i, j) for i in 1:n_1, j in 1:n_2])

    ordered((r, s)) = r < s ? (r, s) : (s, r)
    edges = Set{NTuple{2,Tuple{Int,Int}}}()
    for i in 1:n_1, j in 1:n_2
        if i < n_1
            push!(edges, ((i, j), (i + 1, j)))
        end
        if j < n_2
            push!(edges, ((i, j), (i, j + 1)))
        end
    end
    to_separate = NTuple{2,Tuple{Int,Int}}[]
    for i in 1:n_1-1, j in 1:n_2-1
        edge, non_edge = cell_diagonals(triangle_direction, i, j)
        push!(edges, ordered(edge))
        push!(to_separate, ordered(non_edge))
    end
    if separate_all
        to_separate =
            [(r, s) for r in J for s in J if r < s && !((r, s) in edges)]
    end

    t = separate_all ? max(1, ceil(Int, log2(2 * (n_1 - 1) * (n_2 - 1)))) : 1
    while true
        sub_model = JuMP.Model(sub_solver)
        JuMP.set_silent(sub_model)
        JuMP.@variable(sub_model, x[1:t, J], Bin)
        JuMP.@variable(sub_model, y[1:t, J], Bin)
        JuMP.@variable(sub_model, z[1:t, to_separate], Bin)
        for k in 1:t
            for r in J
                JuMP.@constraint(sub_model, x[k, r] + y[k, r] ≤ 1)
            end
            for (r, s) in edges
                JuMP.@constraint(sub_model, x[k, r] + y[k, s] ≤ 1)
                JuMP.@constraint(sub_model, x[k, s] + y[k, r] ≤ 1)
            end
            # z[k, (r, s)] = 1 implies r and s are separated by level k
            for (r, s) in to_separate
                JuMP.@constraints(
                    sub_model,
                    begin
                        z[k, (r, s)] ≤ x[k, r] + x[k, s]
                        z[k, (r, s)] ≤ y[k, r] + y[k, s]
                        z[k, (r, s)] ≤ x[k, r] + y[k, r]
                        z[k, (r, s)] ≤ x[k, s] + y[k, s]
                    end
                )
            end
        end
        for p in to_separate
            JuMP.@constraint(sub_model, sum(z[k, p] for k in 1:t) ≥ 1)
        end
        remaining = deadline - time()
        if remaining ≤ 0
            return nothing
        end
        JuMP.set_time_limit_sec(sub_model, remaining)
        JuMP.optimize!(sub_model)

        # Any feasible cover gives a valid formulation, and it has the smallest
        # number of levels, as smaller numbers have been proven infeasible.
        status = JuMP.termination_status(sub_model)
        if JuMP.primal_status(sub_model) == JuMP.FEASIBLE_POINT
            cover = [
                (
                    [r for r in J if JuMP.value(x[k, r]) > 0.5],
                    [r for r in J if JuMP.value(y[k, r]) > 0.5],
                ) for k in 1:t
            ]
            return _sparsify_biclique_cover!(cover, to_separate)
        elseif status in (JuMP.INFEASIBLE, JuMP.INFEASIBLE_OR_UNBOUNDED)
            t += 1
        elseif status == JuMP.TIME_LIMIT
            return nothing
        else
            error("Unable to compute an optimal biclique cover: $status")
        end
    end
end

# Removes points from the levels of `cover` as long as every pair in
# `to_separate` is still separated by some level. Removing a point never
# separates a pair that shares a triangle, so the cover stays valid.
function _sparsify_biclique_cover!(
    cover::_BicliqueCover,
    to_separate::Vector{NTuple{2,Tuple{Int,Int}}},
)
    separated(A, B, (r, s)) = (r in A && s in B) || (s in A && r in B)
    count = Dict(
        p => sum(separated(A, B, p) for (A, B) in cover) for p in to_separate
    )
    for (A, B) in cover, side in (A, B)
        for r in reverse(side)
            needed = false
            for p in to_separate
                if r in p && count[p] == 1 && separated(A, B, p)
                    needed = true
                    break
                end
            end
            if !needed
                for p in to_separate
                    if r in p && separated(A, B, p)
                        count[p] -= 1
                    end
                end
                deleteat!(side, findfirst(==(r), side))
            end
        end
    end
    return cover
end
