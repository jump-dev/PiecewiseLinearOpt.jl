abstract type Method end
abstract type UnivariateMethod <: Method end

@enum DIRECTION Graph Epigraph Hypograph

# TODO: Make eltypes of input_vals and output_vals a type parameter
abstract type Segment{D, F} end

struct SegmentPointRep{D, F} <: Segment{D, F}
    input_vals::Vector{NTuple{D, Float64}}
    output_vals::Vector{NTuple{F, Float64}}

    function SegmentPointRep{D,F}(input_vals::Vector{NTuple{D,Float64}}, output_vals::Vector{NTuple{F,Float64}}) where {D,F}
        if length(input_vals) != length(output_vals)
            error("Must specify the same number of input and output values.")
        end
        # TODO: Run verifier to ensure this is actually a PWL function
        return new{D,F}(input_vals, output_vals)
    end
end

struct AffineFunction{D}
    coeffs::NTuple{D, Float64}
    offset::Float64
end

struct SegmentHyperplaneRep{D,F} <: Segment{D,F}
    # Domain given by f_i(x) >= 0 where f_i is i-th constraint in constraints
    constraints::Vector{AffineFunction{D}}
    funcs::NTuple{F, AffineFunction{D}}
end

struct PWLFunction{D, F, T <: Segment{D, F}}
    segments::Vector{T}
    meta::Dict

    function PWLFunction{D, F, T}(segments::Vector{T}) where {D, F, T <: Segment}
        return new(segments, Dict())
    end
end

const PWLFunctionPointRep{D, F} = PWLFunction{D, F, SegmentPointRep{D, F}}
const PWLFunctionHyperplaneRep{D, F} = PWLFunction{D, F, SegmentHyperplaneRep{D, F}}

const UnivariatePWLFunction{F} = PWLFunctionPointRep{1, F}
const BivariatePWLFunction{F} = PWLFunctionPointRep{2, F}
