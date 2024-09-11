# Copyright (c) 2016: Joey Huchette and contributors
#
# Use of this source code is governed by an MIT-style license that can be found
# in the LICENSE.md file or at https://opensource.org/licenses/MIT.

module PiecewiseLinearOpt

using JuMP
import MathOptInterface
const MOI = MathOptInterface
using LinearAlgebra
using Random

export PWLFunction, UnivariatePWLFunction, BivariatePWLFunction, piecewiselinear

include("types.jl")
include("jump.jl")

end # module
