# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# The extension interface, with the methods from the manual

import PiecewiseLinearOpt:
    SOS2Method,
    TriangleSelectionMethod,
    formulate_sos2!,
    formulate_triangle_selection!,
    add_biclique_cover!

# An SOS2 formulation with a binary for each k that separates λ[1:k] from
# λ[k+2:n]
struct LinearBranching <: SOS2Method end

function formulate_sos2!(formulation, λ::Vector, ::LinearBranching)
    n = length(λ)
    cover = [(collect(1:k), collect(k+2:n)) for k in 1:n-2]
    add_biclique_cover!(formulation, λ, cover)
    return
end

# A triangle selection with a binary for each cell
struct CellBinaries{M<:PiecewiseLinearOpt.Method} <: TriangleSelectionMethod
    axis_method::M
end
CellBinaries() = CellBinaries(Logarithmic())

function formulate_triangle_selection!(
    formulation,
    λ::Matrix,
    nw_se::Matrix{Bool},
    ::CellBinaries,
)
    cover = [
        nw_se[i, j] ? ([(i, j)], [(i + 1, j + 1)]) :
        ([(i, j + 1)], [(i + 1, j)]) for i in 1:size(nw_se, 1),
        j in 1:size(nw_se, 2)
    ]
    add_biclique_cover!(formulation, λ, vec(cover))
    return
end

@testset "Extension: SOS2 method" begin
    # Pinned to the interpolated value inside every segment, with n - 2
    # binaries for n breakpoints
    for n in 2:8
        xs = collect(range(0; stop = 2π, length = n))
        zs = sin.(xs)
        model = Model(optimizer)
        @variable(model, x)
        (y,), formulation = piecewiselinear(
            model,
            [x],
            UnivariatePWLFunction(xs, zs);
            method = LinearBranching(),
        )
        @test count(is_binary, formulation.variables) == max(n - 2, 0)
        for i in 1:n-1
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
    test_triangle_selection(SixStencil(LinearBranching()), K1Triangulation(), 4)
    # Recorded and deleted like the built-in methods
    model = Model()
    @variable(model, x)
    test_formulation(
        model,
        [x],
        UnivariatePWLFunction(0:5, t -> t^2),
        LinearBranching(),
    )
end

@testset "Extension: triangle selection method" begin
    for triangulation in triangulations
        test_triangle_selection(CellBinaries(), triangulation, 4)
    end
    # One binary per cell, in addition to those of the axis method
    model = Model()
    @variable(model, x[1:2])
    pwl = BivariatePWLFunction(0:3, 0:3, (a, b) -> a * b)
    _, formulation = piecewiselinear(model, x, pwl; method = CellBinaries())
    _, reference = piecewiselinear(model, x, pwl; method = SixStencil())
    @test count(is_binary, formulation.variables) ==
          count(is_binary, reference.variables) - 6 + 9
    # The axis method is taken from the field by default
    @test PiecewiseLinearOpt.axis_method(CellBinaries(NativeSOS2())) ==
          NativeSOS2()
    model = Model()
    @variable(model, x[1:2])
    test_formulation(model, x, pwl, CellBinaries(IndependentBranching()))
end
