using PiecewiseLinearOpt
using Base.Test

using JuMP, Cbc

const solver = CbcSolver()
methods_1D = (:CC,:MC,:Logarithmic,:ZigZag,:ZigZagInteger,:SOS2,:Incremental)
methods_2D = (:CC,:MC,:Logarithmic,:ZigZag,:ZigZagInteger,:SOS2)

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

let dˣ = linspace(1,π,8), dʸ = linspace(1,π,8), f = (x,y) -> (π-x^2)*(π-y)*sin(x+2y)
    for method in methods_2D
        println("Method: $method")
        model = Model(solver=solver)
        @variable(model, x)
        @variable(model, y)
        z = piecewiselinear(model, x, y, dˣ, dʸ, f, method=method)
        @objective(model, Max, z)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(x),   3.14159, rtol=1e-3)
        @test isapprox(getvalue(y),   1.00000, rtol=1e-3)
        @test isapprox(getvalue(z), 13.101758, rtol=1e-3)

        @constraint(model, x+y ≤ 3)
        @test solve(model) == :Optimal
        @test isapprox(getvalue(x), 2.00000, rtol=1e-3)
        @test isapprox(getvalue(y), 1.00000, rtol=1e-3)
        @test isapprox(getvalue(z), 1.50475, rtol=1e-3)
    end
end
