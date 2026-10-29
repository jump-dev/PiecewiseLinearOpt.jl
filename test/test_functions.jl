# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Construction and evaluation of piecewise linear functions

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "utilities.jl"))

module TestFunctions

using JuMP
using PiecewiseLinearOpt
using Test

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

function test_triangulations()
    # A single cell with f = (a + b)^2: the corner values are SW = 0,
    # SE = NW = 1 and NE = 4, so the value at the center is 2 when the cell is
    # split along SW-NE, and 1 along NW-SE. The true value there is 1.
    f = (a, b) -> (a + b)^2
    function center_value(t)
        return BivariatePWLFunction(0:1, 0:1, f; triangulation = t)(0.5, 0.5)
    end
    @test center_value(UpperTriangulation()) ≈ 2
    @test center_value(LowerTriangulation()) ≈ 1
    @test center_value(BestFitTriangulation()) ≈ 1
    @test center_value(K1Triangulation()) ≈ 1
    # Union Jack alternates, starting with SW-NE in the first cell, while K1
    # always splits along NW-SE
    d = 0:3
    pwl =
        BivariatePWLFunction(d, d, f; triangulation = UnionJackTriangulation())
    @test pwl(0.5, 0.5) ≈ 2  # (f(0, 0) + f(1, 1)) / 2
    @test pwl(1.5, 0.5) ≈ 4  # (f(1, 1) + f(2, 0)) / 2
    @test pwl(1.5, 1.5) ≈ 10  # (f(1, 1) + f(2, 2)) / 2
    pwl = BivariatePWLFunction(d, d, f; triangulation = K1Triangulation())
    @test pwl(0.5, 0.5) ≈ 1  # (f(0, 1) + f(1, 0)) / 2
    @test pwl(1.5, 1.5) ≈ 9  # (f(1, 2) + f(2, 1)) / 2
    # Random triangulations are reproducible
    function vertices(t)
        return [
            s.input_vals for
            s in BivariatePWLFunction(d, d, f; triangulation = t).segments
        ]
    end
    @test vertices(RandomTriangulation(1)) == vertices(RandomTriangulation(1))
    # The triangulation is stored as the structure of the function
    @test pwl.structure == K1Triangulation()
    return
end

function test_matrix_values()
    # Function values can be given as a matrix
    f = (a, b) -> (a + b)^2
    d = 0:3
    z = [f(a, b) for a in d, b in d]
    for t in (K1Triangulation(), UpperTriangulation(), UnionJackTriangulation())
        pwl_f = BivariatePWLFunction(d, d, f; triangulation = t)
        pwl_z = BivariatePWLFunction(d, d, z; triangulation = t)
        for a in 0:0.25:3, b in 0:0.25:3
            @test pwl_z(a, b) ≈ pwl_f(a, b)
        end
    end
    @test_throws ErrorException BivariatePWLFunction(d, d, z[1:3, :])
    @test_throws ErrorException BivariatePWLFunction(
        d,
        d,
        z;
        triangulation = BestFitTriangulation(),
    )
    # Matrix input in PiecewiseLinear
    pwl = BivariatePWLFunction(d, d, f; triangulation = K1Triangulation())
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    w = PiecewiseLinear(d, d, z)(x[1], x[2])
    fix(x[1], 1.5)
    fix(x[2], 0.5)
    @objective(model, Max, w)
    solve!(model)
    @test value(w) ≈ pwl(1.5, 0.5) atol = 1e-5
    return
end

function test_evaluation()
    pwl = UnivariatePWLFunction([0, 1, 2], [0, 2, 1])
    @test pwl(0.5) ≈ 1
    @test pwl(1.5) ≈ 1.5
    @test_throws ErrorException pwl(2.5)
    # Several outputs are returned as a tuple
    segments =
        [SegmentPointRep{1,2}([(0.0,), (1.0,)], [(0.0, 1.0), (2.0, -1.0)])]
    @test PWLFunction(segments, Intervals())(0.25) == (0.5, 0.5)
    # Outside the domain
    f = (a, b) -> (a + b)^2
    @test_throws ErrorException BivariatePWLFunction(0:3, 0:3, f)(4.0, 0.0)
    # Evaluation is fast after the first call, which prepares the segments
    big = BivariatePWLFunction(0:0.02:1, 0:0.02:1, f)
    @test big(0.123, 0.456) ≈ big(0.123, 0.456)
    @test (@elapsed for _ in 1:100
        big(0.99, 0.99)
    end) < 1.0
    return
end

function test_input_validation()
    @test_throws ErrorException UnivariatePWLFunction([0, 2, 1], [0, 1, 2])
    @test_throws ErrorException UnivariatePWLFunction([0, 1, 1], [0, 1, 2])
    @test_throws ErrorException BivariatePWLFunction(
        [0, 1],
        [1, 0],
        (a, b) -> a + b,
    )
    @test_throws ErrorException BivariatePWLFunction(
        [0, 1],
        [0, 1],
        (a, b) -> a + b;
        triangulation = UnstructuredTriangulation(),
    )
    return
end

function test_relaxation_construction()
    f = t -> t^2
    xs = 0:0.5:2
    pwl = UnivariatePWLRelaxation(xs, t -> f(t) - 0.1, t -> f(t) + 0.1)
    @test length(segments(pwl)) == 4
    @test input_values(segments(pwl)[1]) == [(0.0,), (0.0,), (0.5,), (0.5,)]
    @test output_values(segments(pwl)[1]) == [(-0.1,), (0.1,), (0.15,), (0.35,)]
    @test pwl.structure == Intervals()
    # Bounds as vectors give the same segments
    lower = [f(t) - 0.1 for t in xs]
    upper = [f(t) + 0.1 for t in xs]
    pwl_z = UnivariatePWLRelaxation(xs, lower, upper)
    @test [s.input_vals for s in segments(pwl_z)] == [s.input_vals for s in segments(pwl)]
    @test [s.output_vals for s in segments(pwl_z)] == [s.output_vals for s in segments(pwl)]
    # Coinciding bounds give a single point
    pwl = UnivariatePWLRelaxation([0, 1, 2], [0, 1, 0], [1, 1, 2])
    @test input_values(segments(pwl)[1]) == [(0.0,), (0.0,), (1.0,)]
    @test output_values(segments(pwl)[2]) == [(1.0,), (0.0,), (2.0,)]
    # A relaxation cannot be evaluated
    @test_throws ErrorException pwl(0.5)
    # Invalid input
    @test_throws r"above the upper bound" UnivariatePWLRelaxation(
        [0, 1],
        [0, 2],
        [1, 1],
    )
    @test_throws r"number of breakpoints" UnivariatePWLRelaxation(
        [0, 1],
        [0, 1, 2],
        [1, 2],
    )
    @test_throws r"number of breakpoints" UnivariatePWLRelaxation(
        [0, 1],
        [0, 1],
        [1],
    )
    @test_throws r"strictly increasing" UnivariatePWLRelaxation(
        [1, 0],
        [0, 1],
        [1, 2],
    )
    @test_throws r"not finite" UnivariatePWLRelaxation([0, 1], [0, 1], [1, Inf])
    # The methods for functions reject relaxations, and say so
    pwl = UnivariatePWLRelaxation(xs, lower, upper)
    model = Model()
    @variable(model, x)
    for method in (
        Logarithmic(),
        ZigZagBinary(),
        NativeSOS2(),
        Incremental(),
        ConvexCombination(),
    )
        @test_throws r"relaxation" add_formulation(
            model,
            PiecewiseLinear(pwl; method = method),
            x,
        )
        check_rejected(model, (x,), pwl; method = method)
    end
    for method in (MultipleChoice(), ConvexPieces())
        check_rejected(model, (x,), pwl; method = method)
    end
    return
end

function test_non_finite_values()
    # Breakpoints and function values must be finite
    @test_throws r"function value -Inf at the breakpoint 0.0 is not finite" UnivariatePWLFunction(
        0:2,
        log,
    )
    @test_throws ErrorException UnivariatePWLFunction([0, 1, Inf], [0, 1, 2])
    @test_throws ErrorException UnivariatePWLFunction([0, 1, 2], [0, NaN, 2])
    @test_throws ErrorException BivariatePWLFunction(0:1, 0:1, (a, b) -> a / b)
    @test_throws ErrorException SegmentPointRep{1,2}(
        [(0.0,), (1.0,)],
        [(0.0, 1.0), (1.0, Inf)],
    )
    @test_throws ErrorException PiecewiseLinear(log, 0:2)
    return
end

end # module

TestFunctions.runtests()
