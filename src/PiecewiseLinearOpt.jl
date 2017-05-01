module PiecewiseLinearOpt

import JuMP, Gurobi, CPLEX

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear

include("types.jl")
include("jump.jl")

end # module
