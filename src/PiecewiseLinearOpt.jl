__precompile__()

module PiecewiseLinearOpt

using JuMP
using LinearAlgebra
using Random

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear

include("types.jl")
include("jump.jl")

end # module
