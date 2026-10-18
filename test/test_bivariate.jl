# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Bivariate functions across the formulations

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "utilities.jl"))

module TestBivariate

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

function test_convenience_methods()
    d = 0:0.05:1
    f = (xi, yi) -> xi^2 + yi^2
    pwl = BivariatePWLFunction(d, d, f)

    model = Model(OPTIMIZER)
    @variable(model, x)
    @variable(model, y)
    z1, _ = piecewiselinear(model, (x, y), pwl)
    @constraint(model, x ≤ 0.75)
    @constraint(model, y ≤ 0.75)
    @objective(model, Max, z1[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 0.75 rtol = 1e-4
    @test value(y) ≈ 0.75 rtol = 1e-4
    @test value(z1[1]) ≈ 1.125 rtol = 1e-4

    model = Model(OPTIMIZER)
    @variable(model, x)
    @variable(model, y)
    z2 = piecewiselinear(model, x, y, d, d, f)
    @constraint(model, x ≤ 0.75)
    @constraint(model, y ≤ 0.75)
    @objective(model, Max, z2)
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x) ≈ 0.75 rtol = 1e-4
    @test value(y) ≈ 0.75 rtol = 1e-4
    @test value(z2) ≈ 1.125 rtol = 1e-4
    return
end

# Several functions in one model: all methods must agree on the optimum
function test_methods_agree()
    d = range(0; stop = 1, length = 3)
    f = (a, b) -> (a - 0.4)^2 - sin(3b) + a * b
    methods = [
        ConvexCombination(),
        DisaggregatedLogarithmic(),
        MultipleChoice(),
        GreedyTriangleSelection(),
        NineStencil(),
        OptimalTriangleSelection(OPTIMIZER),
        SixStencil(),
        OptimalIndependentBranching(OPTIMIZER),
        K1(),
    ]
    objectives = map(methods) do method
        model = Model(OPTIMIZER)
        @variable(model, 0 ≤ a[1:2] ≤ 1)
        @variable(model, 0 ≤ b[1:2] ≤ 1)
        @constraint(model, sum(a) == 1.1)
        @constraint(model, sum(b) == 0.7)
        @objective(
            model,
            Min,
            sum(
                piecewiselinear(model, a[i], b[i], d, d, f; method = method) for i in 1:2
            )
        )
        solve!(model)
        @test termination_status(model) == MOI.OPTIMAL
        return objective_value(model)
    end
    @test all(isapprox.(objectives, objectives[1]; atol = 1e-5))
    return
end

# The default method depends on the triangulation
function test_default_method()
    f = (a, b) -> a^2 - b
    for (triangulation, method) in (
        (K1Triangulation(), K1()),
        (UnionJackTriangulation(), UnionJack()),
        (UpperTriangulation(), GreedyTriangleSelection()),
        (BestFitTriangulation(), GreedyTriangleSelection()),
    )
        pwl = BivariatePWLFunction(0:3, 0:3, f; triangulation = triangulation)
        model = Model(OPTIMIZER)
        @variable(model, x[1:2])
        (z,), formulation = piecewiselinear(model, x, pwl)
        @test formulation.method == method
        fix(x[1], 1.25)
        fix(x[2], 2.5)
        @objective(model, Min, z)
        solve!(model)
        @test value(z) ≈ pwl(1.25, 2.5) atol = 1e-5
    end
    # The convenience method uses the K1 triangulation by default
    model = Model()
    @variable(model, x)
    @variable(model, y)
    piecewiselinear(model, x, y, 0:3, 0:3, f)
    @test count(is_binary, all_variables(model)) == 2 + 2 + 2
    # Segments without a known structure get a method for any triangulation
    pwl = BivariatePWLFunction(0:3, 0:3, f)
    unstructured = PWLFunction(segments(pwl), UnstructuredTriangulation())
    model = Model()
    @variable(model, x[1:2])
    _, formulation = piecewiselinear(model, x, unstructured)
    @test formulation.method == GreedyTriangleSelection()
    unstructured = PWLFunction(segments(pwl), Unstructured{2}())
    _, formulation = piecewiselinear(model, x, unstructured)
    @test formulation.method == DisaggregatedLogarithmic()
    return
end

# Levels of the stencils that separate nothing get no binary variable
function test_stencil_levels()
    f = (a, b) -> a * b
    # One cell: one level of the six-stencil separates its two corners, and
    # the axes need no binary variables
    pwl = BivariatePWLFunction(0:1, 0:1, f)
    model = Model()
    @variable(model, x[1:2])
    _, formulation = piecewiselinear(model, x, pwl; method = SixStencil())
    @test count(is_binary, formulation.variables) == 1
    # With every diagonal in the same direction, the three levels for the
    # other direction are left out
    pwl = BivariatePWLFunction(0:4, 0:4, f)
    model = Model()
    @variable(model, x[1:2])
    _, formulation = piecewiselinear(model, x, pwl; method = SixStencil())
    @test count(is_binary, formulation.variables) == 2 + 2 + 3
    _, formulation = piecewiselinear(model, x, pwl; method = NineStencil())
    @test count(is_binary, formulation.variables) == 2 + 2 + 9
    # Every constraint of a level has multipliers on both sides
    for method in (SixStencil(), NineStencil(), K1(), UnionJack())
        for n in (2, 3, 4, 5)
            pwl = BivariatePWLFunction(
                0:(n-1),
                0:(n-1),
                f;
                triangulation = method isa UnionJack ?
                                UnionJackTriangulation() :
                                K1Triangulation(),
            )
            model = Model()
            @variable(model, x[1:2])
            _, formulation = piecewiselinear(model, x, pwl; method = method)
            for c in formulation.constraints
                func = constraint_object(c).func
                @test !(func isa AffExpr) || length(func.terms) ≥ 2
            end
        end
    end
    return
end

# Checks that a triangle selection is a valid biclique cover: no level
# separates two points that share a triangle, and the ends of the diagonal
# that is not an edge are separated in every cell.
function is_valid_triangle_selection(nw_se, cover)
    n_1, n_2 = size(nw_se) .+ 1
    for (L, R) in cover
        if !isdisjoint(L, R)
            return false
        end
        side = Dict((p => 1 for p in L)..., (p => 2 for p in R)...)
        separated((p, q)) = get(side, p, 0) + get(side, q, 0) == 3
        for i in 1:n_1, j in 1:n_2
            if (i < n_1 && separated(((i, j), (i + 1, j)))) ||
               (j < n_2 && separated(((i, j), (i, j + 1))))
                return false
            end
        end
        for i in 1:(n_1-1), j in 1:(n_2-1)
            if separated(first(PLO.cell_diagonals(nw_se, i, j)))
                return false
            end
        end
    end
    for i in 1:(n_1-1), j in 1:(n_2-1)
        p, q = last(PLO.cell_diagonals(nw_se, i, j))
        if !any(((L, R),) -> (p in L && q in R) || (q in L && p in R), cover)
            return false
        end
    end
    return true
end

function test_triangle_selection_covers()
    rng = Random.MersenneTwister(5)
    for n_1 in 2:6, n_2 in 2:6
        k1 = fill(true, n_1 - 1, n_2 - 1)
        union_jack = [isodd(i + j) for i in 1:(n_1-1), j in 1:(n_2-1)]
        random = rand(rng, Bool, n_1 - 1, n_2 - 1)
        flipped = [!x for x in union_jack]
        @test is_valid_triangle_selection(
            k1,
            PLO._triangle_selection_cover(K1(), k1),
        )
        @test is_valid_triangle_selection(
            union_jack,
            PLO._triangle_selection_cover(UnionJack(), union_jack),
        )
        @test is_valid_triangle_selection(
            flipped,
            PLO._triangle_selection_cover(UnionJack(), flipped),
        )
        for method in (SixStencil(), NineStencil(), GreedyTriangleSelection())
            for nw_se in (k1, union_jack, random)
                cover = PLO._triangle_selection_cover(method, nw_se)
                @test is_valid_triangle_selection(nw_se, cover)
            end
        end
    end
    @test_throws ErrorException PLO._triangle_selection_cover(
        K1(),
        [true false; true true],
    )
    @test_throws ErrorException PLO._triangle_selection_cover(
        UnionJack(),
        [true true; true true],
    )
    return
end

end # module

TestBivariate.runtests()
