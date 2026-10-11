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

@testset "Optimal covers: time limit" begin
    # A triangulation that no other test uses, so that its cover is not cached
    pwl = BivariatePWLFunction(
        0:5,
        0:5,
        (a, b) -> a * b;
        triangulation = RandomTriangulation(12345),
    )
    # OptimalTriangleSelection falls back to SixStencil
    model = Model(optimizer)
    @variable(model, x[1:2])
    method = OptimalTriangleSelection(optimizer; time_limit = 0)
    (z,), formulation =
        @test_logs (:warn,) piecewiselinear(model, x, pwl; method = method)
    @test count(v -> occursin("w_sw_ne", name(v)), formulation.variables) == 3
    fix(x[1], 2.25)
    fix(x[2], 3.5)
    @objective(model, Max, z)
    solve!(model)
    @test value(z) ≈ pwl(2.25, 3.5) atol = 1e-5
    # OptimalIndependentBranching throws an error
    model = Model()
    @variable(model, x[1:2])
    test_rejected(
        model,
        x,
        pwl;
        method = OptimalIndependentBranching(optimizer; time_limit = 0),
    )
end
