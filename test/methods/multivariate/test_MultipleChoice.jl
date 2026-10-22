# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Tests of the MultipleChoice formulation

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "..", "..", "utilities.jl"))

module TestMultipleChoice

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

function test_univariate()
    check_univariate(MultipleChoice())
    return
end

function test_bivariate_simple()
    check_bivariate_simple(MultipleChoice())
    return
end

function test_bivariate_grid()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_bivariate_grid(MultipleChoice(), triangulation)
        end
    end
    return
end

function test_triangle_selection()
    for triangulation in TRIANGULATIONS
        @testset "$triangulation" begin
            check_triangle_selection(MultipleChoice(), triangulation, 4)
        end
    end
    return
end

function test_formulation()
    check_bivariate_formulation(MultipleChoice(), UpperTriangulation())
    return
end

function test_trivariate()
    check_trivariate(MultipleChoice())
    return
end

function test_discontinuous()
    check_discontinuous(MultipleChoice())
    return
end

function test_discontinuous_bivariate()
    # Two triangles sharing the diagonal from (0, 0) to (1, 1), with a jump
    # across it.
    s1 = SegmentPointRep{2,1}(
        [(0.0, 0.0), (0.0, 1.0), (1.0, 1.0)],
        [(0.0,), (1.0,), (2.0,)],
    )
    s2 = SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0)],
        [(5.0,), (3.0,), (7.0,)],
    )
    pwl = PWLFunction([s1, s2], UnstructuredTriangulation())
    model = Model(OPTIMIZER)
    @variable(model, x[1:2])
    y, _ = piecewiselinear(model, (x[1], x[2]), pwl; method = MultipleChoice())
    @objective(model, Max, y[1])
    solve!(model)
    @test termination_status(model) == MOI.OPTIMAL
    @test value(x[1]) ≈ 1.0 atol = 1e-6
    @test value(x[2]) ≈ 1.0 atol = 1e-6
    @test value(y[1]) ≈ 7.0 atol = 1e-6
    # On the shared edge the lower triangle (s2) gives the larger value
    @constraint(model, x[1] == 0.5)
    @constraint(model, x[2] == 0.5)
    solve!(model)
    @test value(y[1]) ≈ 6.0 atol = 1e-6
    return
end

function test_affine_inputs()
    check_bivariate_affine(MultipleChoice())
    return
end

function test_direction()
    check_bivariate_direction(MultipleChoice())
    return
end

function test_multiple_outputs()
    check_bivariate_multiple_outputs(MultipleChoice())
    return
end

function test_spanning_triangulation()
    check_spanning_triangulation(MultipleChoice())
    return
end

function test_non_simplices()
    # Segments that are not simplices cannot be used
    s = SegmentPointRep{2,1}(
        [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)],
        [(0.0,), (1.0,), (2.0,), (1.0,)],
    )
    pwl = PWLFunction([s], UnstructuredTriangulation())
    model = Model()
    @variable(model, x[1:2])
    check_rejected(model, (x[1], x[2]), pwl; method = MultipleChoice())
    return
end

function test_coefficients()
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
        @variable(model, x[1:(length(pwl.segments[1].input_vals)-1)])
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
    return
end

end # module

TestMultipleChoice.runtests()
