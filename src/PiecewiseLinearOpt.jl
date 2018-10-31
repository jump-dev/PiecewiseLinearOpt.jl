__precompile__()

module PiecewiseLinearOpt

using JuMP

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear

include("types.jl")
include("jump.jl")

end # module
