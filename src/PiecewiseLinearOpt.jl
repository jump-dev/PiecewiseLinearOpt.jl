__precompile__()

module PiecewiseLinearOpt

import JuMP

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear

include("types.jl")
include("jump.jl")

end # module
