# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the IndependentBranching formulation and its biclique covers

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestIndependentBranching

using JuMP
using PiecewiseLinearOpt
using Test

import PiecewiseLinearOpt as PLO

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

const SCHEMES = (:biclique, :balanced, :brgc)

function test_univariate()
    for scheme in SCHEMES
        @testset "$scheme" begin
            check_univariate(IndependentBranching(scheme))
        end
    end
    return
end

function test_discontinuous()
    for scheme in SCHEMES
        @testset "$scheme" begin
            check_discontinuous(IndependentBranching(scheme))
        end
    end
    return
end

function test_relaxation()
    for scheme in SCHEMES
        @testset "$scheme" begin
            check_relaxation(IndependentBranching(scheme))
        end
    end
    # Both bounds are modeled with the binary variables of a single function
    # of 4 segments
    pwl = UnivariatePWLRelaxation(0:0.5:2, t -> t^2 - 0.1, t -> t^2 + 0.1)
    model = Model()
    @variable(model, x)
    _, formulation = add_formulation(
        model,
        PiecewiseLinear(pwl; method = IndependentBranching()),
        x,
    )
    @test count(is_binary, formulation.variables) == 2
    return
end

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

function test_ordered_cdc_covers()
    for scheme in SCHEMES, d in 1:20
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
                push!(S, collect(first:(first+2*shared+mod(i, 3))))
                next = S[end][end] + 1
            end
            cover = PLO._ordered_cdc_cover(S, scheme)
            @test is_valid_cover(S, cover)
            @test length(cover) == (d == 1 ? 0 : ceil(Int, log2(d)))
        end
    end
    # Remark 4: the first part of each cut gets floor(n / 2) of the n sets
    @test PLO._balanced_reversed_edge_ranking(6) == [2, 3, 1, 2, 3]
    return
end

function test_binaries()
    # ceil(log2(d)) binaries for d segments, also for discontinuous functions
    for scheme in SCHEMES, d in (2, 3, 5, 8, 13), jump in (0.0, 1.0)
        segments = [
            SegmentPointRep{1,1}(
                [(float(i),), (i + 1.0,)],
                [(float(i),), (i + 1.0 + jump * isodd(i),)],
            ) for i in 1:d
        ]
        pwl = PWLFunction(segments, Intervals())
        model = Model()
        @variable(model, x)
        add_formulation(
            model,
            PiecewiseLinear(pwl; method = IndependentBranching(scheme)),
            x,
        )
        @test num_constraints(model, VariableRef, MOI.ZeroOne) ==
              ceil(Int, log2(d))
    end
    return
end

function test_segment_order()
    # Segments must be ordered so that only neighbors share points
    s1 = SegmentPointRep{1,1}([(0.0,), (2.0,)], [(0.0,), (1.0,)])
    s2 = SegmentPointRep{1,1}([(1.0,), (3.0,)], [(1.0,), (0.0,)])
    s3 = SegmentPointRep{1,1}([(2.0,), (4.0,)], [(1.0,), (2.0,)])
    pwl = PWLFunction([s1, s2, s3], Intervals())
    for scheme in SCHEMES
        model = Model()
        @variable(model, x)
        check_rejected(model, (x,), pwl; method = IndependentBranching(scheme))
    end
    return
end

function test_options()
    @test_throws ErrorException IndependentBranching(:unknown)
    return
end

end # module

TestIndependentBranching.runtests()
