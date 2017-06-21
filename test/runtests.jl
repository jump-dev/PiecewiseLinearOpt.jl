using PiecewiseLinearOpt
using Base.Test

using JuMP, Cbc

const solver = CbcSolver()
methods_1D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:Incremental,:DisaggLogarithmic)
methods_2D = (:CC,:MC,:Logarithmic,:LogarithmicIB,:ZigZag,:ZigZagInteger,:SOS2,:GeneralizedCelaya,:SymmetricCelaya,:DisaggLogarithmic)

let d = linspace(1,2π,8), f = sin
    for method in methods_1D
        println("Method: $method")
        model = Model(solver=solver)
        @variable(model, x)
        z = piecewiselinear(model, x, d, sin, method=method)
        @objective(model, Max, z)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(x), 1.75474, rtol=1e-4)
        @test isapprox(getvalue(z), 0.98313, rtol=1e-4)

        @constraint(model, x ≤ 1.5z)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(x), 1.36495, rtol=1e-4)
        @test isapprox(getvalue(z), 0.90997, rtol=1e-4)
    end
end

let dˣ = linspace(0,1,8), dʸ = linspace(0,1,8), f = (x,y) -> 2*(x-1/3)^2 + 3*(y-4/7)^4
    for method in methods_2D
        println("Method: $method")
        model = Model(solver=solver)
        @variable(model, x)
        @variable(model, y)
        # z = piecewiselinear(model, x, y, dˣ, dʸ, f, method=method)
        z = piecewiselinear(model, x, y, BivariatePWLFunction(dˣ, dʸ, f, pattern=:UnionJack), method=method)
        @objective(model, Min, z)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(x), 0.285714, rtol=1e-4)
        @test isapprox(getvalue(y), 0.571429, rtol=1e-4)
        @test isapprox(getvalue(z), 0.004535, atol=1e-4)

        @constraint(model, x ≥ 0.6)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(x), 0.6, rtol=1e-4)
        @test isapprox(getvalue(y), 0.571428, rtol=1e-4)
        @test isapprox(getvalue(z), 0.148753, rtol=1e-4)
    end
end
