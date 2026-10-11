# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

# Formulations of functions of three or more variables

@testset "Trivariate: $method" for method in [
    ConvexCombination(),
    DisaggregatedLogarithmic(),
    MultipleChoice(),
]
    # Kuhn triangulation of the unit cube into 6 tetrahedra
    g = v -> v[1] + 2v[2]^2 - 3v[1] * v[3]
    segments = PLO.SegmentPointRep{3,1}[]
    for σ in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (2, 3, 1), (3, 1, 2), (3, 2, 1))
        vertices = [(0.0, 0.0, 0.0)]
        for j in σ
            push!(vertices, ntuple(i -> i == j ? 1.0 : vertices[end][i], 3))
        end
        push!(
            segments,
            PLO.SegmentPointRep{3,1}(vertices, [(g(v),) for v in vertices]),
        )
    end
    pwl = PLO.PWLFunction(segments, PLO.Unstructured{3}())
    model = Model(optimizer)
    @variable(model, x[1:3])
    (y,), _ = piecewiselinear(model, (x[1], x[2], x[3]), pwl; method = method)
    for seg in segments
        for i in 1:3
            fix(x[i], sum(v[i] for v in seg.input_vals) / 4)
        end
        y0 = sum(v[1] for v in seg.output_vals) / 4
        for sense in (MAX_SENSE, MIN_SENSE)
            @objective(model, sense, y)
            solve!(model)
            @test termination_status(model) == MOI.OPTIMAL
            @test value(y) ≈ y0 atol = 1e-5
        end
    end
    # Functions of three or more variables have a default method
    model = Model()
    @variable(model, x[1:3])
    _, formulation = piecewiselinear(model, x, pwl)
    @test formulation.method == DisaggregatedLogarithmic()
end
