# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Construction and evaluation of piecewise linear functions

@testset "Triangulations" begin
    # A single cell with f = (a + b)^2: the corner values are SW = 0,
    # SE = NW = 1 and NE = 4, so the value at the center is 2 when the cell is
    # split along SW-NE, and 1 along NW-SE. The true value there is 1.
    f = (a, b) -> (a + b)^2
    center_value(t) =
        BivariatePWLFunction(0:1, 0:1, f; triangulation = t)(0.5, 0.5)
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
    vertices(t) = [
        s.input_vals for
        s in BivariatePWLFunction(d, d, f; triangulation = t).segments
    ]
    @test vertices(RandomTriangulation(1)) == vertices(RandomTriangulation(1))
    # The triangulation is stored as the structure of the function
    @test pwl.structure == K1Triangulation()
    # Function values can be given as a matrix
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
    # Matrix input in the convenience method
    model = Model(optimizer)
    @variable(model, x[1:2])
    w = piecewiselinear(model, x[1], x[2], d, d, z)
    fix(x[1], 1.5)
    fix(x[2], 0.5)
    @objective(model, Max, w)
    solve!(model)
    @test value(w) ≈ pwl(1.5, 0.5) atol = 1e-5
end

@testset "Evaluation" begin
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
end

@testset "Input validation" begin
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
    # Segments that are not simplices cannot be used with MultipleChoice
    s = PLO.SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (1.0,)],
    )
    pwl = PLO.PWLFunction([s], PLO.UnstructuredTriangulation())
    model = Model()
    @variable(model, x[1:2])
    test_rejected(model, (x[1], x[2]), pwl; method = MultipleChoice())
end

@testset "MultipleChoice coefficients" begin
    # Breakpoints like 0.1 are not exact in floating point, which must not
    # leave round-off errors in coefficients that should be zero
    for pwl in (
        UnivariatePWLFunction(0:0.1:1, t -> t^2),
        BivariatePWLFunction(
            0:0.1:0.5,
            0:0.3:0.9,
            (a, b) -> a * b;
            triangulation = UnionJackTriangulation(),
        ),
    )
        model = Model()
        @variable(model, x[1:length(pwl.segments[1].input_vals)-1])
        _, formulation =
            piecewiselinear(model, x, pwl; method = MultipleChoice())
        coefficients = Float64[]
        for c in formulation.constraints
            object = constraint_object(c)
            append!(coefficients, values(object.func.terms))
            push!(coefficients, MOI.constant(object.set))
        end
        @test !isempty(coefficients)
        @test all(a -> iszero(a) || abs(a) > 1e-9, coefficients)
    end
end
