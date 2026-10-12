# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Independent branching formulations and their biclique covers

# Checks that the levels (L, R) form a biclique cover of the conflict graph of
# the CDC given by S: no level separates two points of the same set, and any
# two points that are not in a common set are separated by some level.
function is_valid_cover(S, cover)
    J = reduce(union, S)
    for (L, R) in cover
        if !isdisjoint(L, R) || !issubset(union(L, R), J)
            return false
        end
        for s in S
            if !isdisjoint(s, L) && !isdisjoint(s, R)
                return false
            end
        end
    end
    for u in J, v in J
        if u < v && !any(s -> u in s && v in s, S)
            if !any(
                ((L, R),) -> (u in L && v in R) || (v in L && u in R),
                cover,
            )
                return false
            end
        end
    end
    return true
end

@testset "Ordered CDC covers: $scheme" for scheme in
                                           (:biclique, :balanced, :brgc)
    for d in 1:20
        # SOS2
        S = [[i, i + 1] for i in 1:d]
        cover = PLO._ordered_cdc_cover(S, scheme)
        @test is_valid_cover(S, cover)
        @test length(cover) == (d == 1 ? 0 : ceil(Int, log2(d)))
        # Polytopes sharing none, one or two points with their neighbors
        for shared in (0, 1, 2)
            S, next = Vector{Int}[], 1
            for i in 1:d
                first = i == 1 ? next : next - shared
                push!(S, collect(first:first+2*shared+mod(i, 3)))
                next = S[end][end] + 1
            end
            cover = PLO._ordered_cdc_cover(S, scheme)
            @test is_valid_cover(S, cover)
            @test length(cover) == (d == 1 ? 0 : ceil(Int, log2(d)))
        end
    end
    # Remark 4: the first part of each cut gets floor(n / 2) of the n sets
    @test PLO._balanced_reversed_edge_ranking(6) == [2, 3, 1, 2, 3]
end

@testset "Independent branching: $scheme" for scheme in
                                              (:biclique, :balanced, :brgc)
    # ceil(log2(d)) binaries for d segments, also for discontinuous functions
    for d in (2, 3, 5, 8, 13), jump in (0.0, 1.0)
        segments = [
            PLO.SegmentPointRep{1,1}(
                [(float(i),), (i + 1.0,)],
                [(float(i),), (i + 1.0 + jump * isodd(i),)],
            ) for i in 1:d
        ]
        pwl = PLO.PWLFunction(segments, PLO.Intervals())
        model = Model()
        @variable(model, x)
        piecewiselinear(model, (x,), pwl; method = IndependentBranching(scheme))
        @test num_constraints(model, VariableRef, MOI.ZeroOne) ==
              ceil(Int, log2(d))
    end
    # Segments must be ordered so that only neighbors share points
    s1 = PLO.SegmentPointRep{1,1}([(0.0,), (2.0,)], [(0.0,), (1.0,)])
    s2 = PLO.SegmentPointRep{1,1}([(1.0,), (3.0,)], [(1.0,), (0.0,)])
    s3 = PLO.SegmentPointRep{1,1}([(2.0,), (4.0,)], [(1.0,), (2.0,)])
    pwl = PLO.PWLFunction([s1, s2, s3], PLO.Intervals())
    model = Model()
    @variable(model, x)
    test_rejected(model, (x,), pwl; method = IndependentBranching(scheme))
    @test_throws ErrorException IndependentBranching(:unknown)
end

# Checks that the levels (A, B) of a triangle selection on the grid with the
# cell diagonals nw_se never separate two points joined by a grid edge or a
# triangle diagonal, and separate the endpoints of the other diagonal of
# every cell.
function is_valid_triangle_selection(nw_se, cover)
    n_1, n_2 = size(nw_se) .+ 1
    joined = Set{NTuple{2,Tuple{Int,Int}}}()
    for i in 1:n_1, j in 1:n_2
        i < n_1 && push!(joined, ((i, j), (i + 1, j)))
        j < n_2 && push!(joined, ((i, j), (i, j + 1)))
    end
    separated(A, B, (r, s)) = (r in A && s in B) || (s in A && r in B)
    for i in 1:n_1-1, j in 1:n_2-1
        edge, non_edge = PLO.cell_diagonals(nw_se, i, j)
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

@testset "Greedy triangle selection covers" begin
    for n_1 in 2:9, n_2 in 2:9
        # One level for Union Jack and two for K1 when there is more than one
        # cell, as for the optimal covers
        uj = [isodd(i + j) for i in 1:n_1-1, j in 1:n_2-1]
        cover = PLO._greedy_biclique_cover(uj, 200)
        @test is_valid_triangle_selection(uj, cover)
        @test length(cover) == 1
        k1 = fill(true, n_1 - 1, n_2 - 1)
        cover = PLO._greedy_biclique_cover(k1, 200)
        @test is_valid_triangle_selection(k1, cover)
        @test length(cover) == (n_1 == n_2 == 2 ? 1 : 2)
        for seed in 1:3
            nw_se = rand(PLO.Random.Xoshiro(seed), Bool, n_1 - 1, n_2 - 1)
            for tries in (0, 20)
                cover = PLO._greedy_biclique_cover(nw_se, tries)
                @test is_valid_triangle_selection(nw_se, cover)
            end
        end
    end
    # More tries never give more levels
    nw_se = rand(PLO.Random.Xoshiro(7), Bool, 16, 16)
    @test length(PLO._greedy_biclique_cover(nw_se, 200)) <=
          length(PLO._greedy_biclique_cover(nw_se, 0))
end

@testset "Optimal covers: time limit" begin
    empty!(PLO._OPTIMAL_COVER_CACHE)
    pwl = BivariatePWLFunction(
        0:5,
        0:5,
        (a, b) -> a * b;
        triangulation = RandomTriangulation(12345),
    )
    nw_se = PLO._triangle_directions(pwl, PLO._continuous_gridpoints_or_die(pwl))
    greedy = PLO._greedy_biclique_cover(nw_se, 200)
    @test length(greedy) > 1
    # OptimalTriangleSelection falls back to the greedy cover, in addition to
    # 3 binaries for each axis
    model = Model(optimizer)
    @variable(model, x[1:2])
    method = OptimalTriangleSelection(optimizer; time_limit = 0)
    (z,), formulation =
        @test_logs (:warn,) piecewiselinear(model, x, pwl; method = method)
    @test count(is_binary, formulation.variables) == 6 + length(greedy)
    fix(x[1], 2.25)
    fix(x[2], 3.5)
    @objective(model, Max, z)
    solve!(model)
    @test value(z) ≈ pwl(2.25, 3.5) atol = 1e-5
    # The fallback is cached, and only computed again with a longer time limit
    @test_logs piecewiselinear(model, x, pwl; method = method)
    @test !PLO._OPTIMAL_COVER_CACHE[(nw_se, false)].optimal
    piecewiselinear(
        model,
        x,
        pwl;
        method = OptimalTriangleSelection(optimizer; time_limit = 60),
    )
    @test PLO._OPTIMAL_COVER_CACHE[(nw_se, false)].optimal
    # OptimalIndependentBranching throws an error, also when cached
    model = Model()
    @variable(model, x[1:2])
    for _ in 1:2
        test_rejected(
            model,
            x,
            pwl;
            method = OptimalIndependentBranching(optimizer; time_limit = 0),
        )
    end
    @test PLO._OPTIMAL_COVER_CACHE[(nw_se, true)].cover === nothing
end
