# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the GreedyTriangleSelection formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestGreedyTriangleSelection

using JuMP
using PiecewiseLinearOpt
using Test

import PiecewiseLinearOpt as PLO
import Random

using ..PLOTestUtilities

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

method(axis_method) = GreedyTriangleSelection(axis_method)

function test_bivariate_simple()
    for axis_method in AXIS_METHODS
        @testset "$axis_method" begin
            check_bivariate_simple(method(axis_method))
        end
    end
    return
end

function test_bivariate_grid()
    for axis_method in AXIS_METHODS, triangulation in TRIANGULATIONS
        @testset "$axis_method, $triangulation" begin
            check_bivariate_grid(method(axis_method), triangulation)
        end
    end
    return
end

function test_triangle_selection()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_triangle_selection(
                GreedyTriangleSelection(),
                triangulation,
                4,
            )
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(GreedyTriangleSelection(), UpperTriangulation())
    return
end

# Checks that the levels (A, B) of a triangle selection on the grid with the
# cell diagonals nw_se never separate two points joined by a grid edge or a
# triangle diagonal, and separate the endpoints of the other diagonal of every
# cell.
function is_valid_triangle_selection(nw_se, cover)
    n_1, n_2 = size(nw_se) .+ 1
    joined = Set{NTuple{2,Tuple{Int,Int}}}()
    for i in 1:n_1, j in 1:n_2
        i < n_1 && push!(joined, ((i, j), (i + 1, j)))
        j < n_2 && push!(joined, ((i, j), (i, j + 1)))
    end
    separated(A, B, (r, s)) = (r in A && s in B) || (s in A && r in B)
    for i in 1:(n_1-1), j in 1:(n_2-1)
        edge, non_edge = cell_diagonals(nw_se, i, j)
        push!(joined, edge)
        if !any(((A, B),) -> separated(A, B, non_edge), cover)
            return false
        end
    end
    for (A, B) in cover
        if !isdisjoint(A, B) || any(p -> separated(A, B, p), joined)
            return false
        end
    end
    return true
end

function test_covers()
    for n_1 in 2:9, n_2 in 2:9
        # One level for Union Jack and two for K1 when there is more than one
        # cell, as for the optimal covers
        uj = [isodd(i + j) for i in 1:(n_1-1), j in 1:(n_2-1)]
        cover = PLO._greedy_biclique_cover(uj, 200)
        @test is_valid_triangle_selection(uj, cover)
        @test length(cover) == 1
        k1 = fill(true, n_1 - 1, n_2 - 1)
        cover = PLO._greedy_biclique_cover(k1, 200)
        @test is_valid_triangle_selection(k1, cover)
        @test length(cover) == (n_1 == n_2 == 2 ? 1 : 2)
        for seed in 1:3
            nw_se = rand(Random.Xoshiro(seed), Bool, n_1 - 1, n_2 - 1)
            for tries in (0, 20)
                cover = PLO._greedy_biclique_cover(nw_se, tries)
                @test is_valid_triangle_selection(nw_se, cover)
            end
        end
    end
    # More tries never give more levels
    nw_se = rand(Random.Xoshiro(7), Bool, 16, 16)
    @test length(PLO._greedy_biclique_cover(nw_se, 200)) <=
          length(PLO._greedy_biclique_cover(nw_se, 0))
    return
end

end # module

TestGreedyTriangleSelection.runtests()
