# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The extension interface, with the methods from the manual

isdefined(@__MODULE__, :PLOTestUtilities) ||
    include(joinpath(@__DIR__, "utilities.jl"))

module TestExtending

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

# An SOS2 formulation with a binary for each k that separates λ[1:k] from
# λ[k+2:n]
struct LinearBranching <: SOS2Method end

function PiecewiseLinearOpt.formulate_sos2!(
    formulation,
    λ::Vector,
    ::LinearBranching,
)
    n = length(λ)
    cover = [(collect(1:k), collect((k+2):n)) for k in 1:(n-2)]
    add_biclique_cover!(formulation, λ, cover)
    return
end

# A triangle selection with a binary for each cell
struct CellBinaries{M<:PWLMethod} <: TriangleSelectionMethod
    axis_method::M
end
CellBinaries() = CellBinaries(Logarithmic())

function PiecewiseLinearOpt.formulate_triangle_selection!(
    formulation,
    λ::Matrix,
    nw_se::Matrix{Bool},
    ::CellBinaries,
)
    cover = Tuple{Vector{Tuple{Int,Int}},Vector{Tuple{Int,Int}}}[]
    for i in 1:size(nw_se, 1), j in 1:size(nw_se, 2)
        _, (a, b) = cell_diagonals(nw_se, i, j)
        push!(cover, ([a], [b]))
    end
    add_biclique_cover!(formulation, λ, cover)
    return
end

# The epigraph of a convex univariate function, without binaries
struct ConvexEpigraph <: PWLMethod end

function PiecewiseLinearOpt.formulate_pwl!(
    formulation,
    (x,)::Tuple{Scalar},
    (y,)::Tuple{Scalar},
    pwl::UnivariatePWLFunction,
    ::ConvexEpigraph,
    direction::Direction,
)
    if !(direction isa Epigraph)
        error("ConvexEpigraph only formulates the epigraph of a function.")
    end
    slope = -Inf
    for segment in segments(pwl)
        (a,), (b,) = input_values(segment)
        (fa,), (fb,) = output_values(segment)
        if (fb - fa) / (b - a) < slope
            error("ConvexEpigraph requires a convex function.")
        end
        slope = (fb - fa) / (b - a)
        constrain_output!(formulation, y, fa + slope * (x - a), direction)
    end
    breakpoints = [only(p) for s in segments(pwl) for p in input_values(s)]
    model = formulation.model
    record!(formulation, @constraint(model, x >= minimum(breakpoints)))
    record!(formulation, @constraint(model, x <= maximum(breakpoints)))
    return
end

function test_sos2_method()
    # Pinned to the interpolated value inside every segment, with n - 2
    # binaries for n breakpoints
    for n in 2:8
        xs = collect(range(0; stop = 2π, length = n))
        zs = sin.(xs)
        model = Model(OPTIMIZER)
        @variable(model, x)
        (y,), formulation = piecewiselinear(
            model,
            [x],
            UnivariatePWLFunction(xs, zs);
            method = LinearBranching(),
        )
        @test count(is_binary, formulation.variables) == max(n - 2, 0)
        for i in 1:(n-1)
            fix(x, 0.3 * xs[i] + 0.7 * xs[i+1])
            for sense in (MAX_SENSE, MIN_SENSE)
                @objective(model, sense, y)
                solve!(model)
                @test termination_status(model) == MOI.OPTIMAL
                @test value(y) ≈ 0.3 * zs[i] + 0.7 * zs[i+1] atol = 1e-5
            end
        end
    end
    # As the axis method of a grid formulation
    check_triangle_selection(
        SixStencil(LinearBranching()),
        K1Triangulation(),
        4,
    )
    # Recorded and deleted like the built-in methods
    model = Model()
    @variable(model, x)
    check_formulation(
        model,
        [x],
        UnivariatePWLFunction(0:5, t -> t^2),
        LinearBranching(),
    )
    return
end

function test_triangle_selection_method()
    for triangulation in TRIANGULATIONS
        check_triangle_selection(CellBinaries(), triangulation, 4)
    end
    # One binary per cell, in addition to those of the axis method
    model = Model()
    @variable(model, x[1:2])
    pwl = BivariatePWLFunction(0:3, 0:3, (a, b) -> a * b)
    _, formulation = piecewiselinear(model, x, pwl; method = CellBinaries())
    _, reference = piecewiselinear(model, x, pwl; method = K1())
    @test count(is_binary, formulation.variables) ==
          count(is_binary, reference.variables) - 2 + 9
    # The axis method is taken from the field by default
    @test axis_method(CellBinaries(NativeSOS2())) == NativeSOS2()
    model = Model()
    @variable(model, x[1:2])
    check_formulation(model, x, pwl, CellBinaries(IndependentBranching()))
    return
end

function test_general_method()
    pwl = UnivariatePWLFunction(-2:2, t -> t^2)
    model = Model(OPTIMIZER)
    @variable(model, x)
    (y,), formulation = piecewiselinear(
        model,
        [x],
        pwl;
        method = ConvexEpigraph(),
        direction = Epigraph(),
    )
    @test !any(is_binary, all_variables(model))
    for t in (-2.0, -0.5, 1.0, 1.75)
        fix(x, t)
        @objective(model, Min, y)
        solve!(model)
        @test value(y) ≈ pwl(t) atol = 1e-6
    end
    # The domain is enforced
    fix(x, 2.5)
    solve!(model)
    @test termination_status(model) in
          (MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED)
    # Unsupported directions and functions are rejected, leaving the model as
    # it was
    model = Model()
    @variable(model, x)
    for (f, direction) in
        ((t -> t^2, Graph()), (t -> t^2, Hypograph()), (t -> -t^2, Epigraph()))
        @test_throws ErrorException piecewiselinear(
            model,
            [x],
            UnivariatePWLFunction(-2:2, f);
            method = ConvexEpigraph(),
            direction = direction,
        )
        @test num_variables(model) == 1
        @test num_constraints(
            model;
            count_variable_in_set_constraints = true,
        ) == 0
    end
    return
end

function test_public_interface()
    for name in (
        :PWLMethod,
        :SOS2Method,
        :TriangleSelectionMethod,
        :Direction,
        :Graph,
        :Scalar,
        :formulate_pwl!,
        :formulate_sos2!,
        :formulate_triangle_selection!,
        :axis_method,
        :record!,
        :variable_name,
        :constrain_output!,
        :add_biclique_cover!,
        :cell_diagonals,
        :segments,
        :input_values,
        :output_values,
    )
        @test Base.isexported(PiecewiseLinearOpt, name)
    end
    # The diagonals of a cell, as an edge and a non-edge of the triangulation
    nw_se = [true false]
    @test cell_diagonals(nw_se, 1, 1) == (((1, 2), (2, 1)), ((1, 1), (2, 2)))
    @test cell_diagonals(nw_se, 1, 2) == (((1, 2), (2, 3)), ((1, 3), (2, 2)))
    # Accessors of the segments
    pwl = UnivariatePWLFunction([0, 1, 3], [1, 2, 0])
    @test length(segments(pwl)) == 2
    @test input_values(segments(pwl)[2]) == [(1.0,), (3.0,)]
    @test output_values(segments(pwl)[2]) == [(2.0,), (0.0,)]
    return
end

end # module

TestExtending.runtests()
