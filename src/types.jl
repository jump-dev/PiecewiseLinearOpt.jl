abstract type Method end
abstract type UnivariateMethod <: Method end

@enum DIRECTION Graph Epigraph Hypograph

abstract type Segment{D,F} end

struct SegmentPointRep{D,F} <: Segment{D,F}
    input_vals::Vector{NTuple{D,Float64}}
    output_vals::Vector{NTuple{F,Float64}}

    function SegmentPointRep{D,F}(input_vals::Vector{NTuple{D,Float64}}, output_vals::Vector{NTuple{F,Float64}}) where {D,F}
        if length(input_vals) != length(output_vals)
            error("Must specify the same number of input and output values.")
        end
        # TODO: Run verifier to ensure this is actually a PWL function
        return new{D,F}(input_vals, output_vals)
    end
end

struct SegmentHyperplaneRep{D,F} <: Segment{D,F}
    breakpoints::Vector{NTuple{D,Float64}}
    coeffs::NTuple{F,NTuple{D,Float64}}
    offsets::NTuple{F,Float64}
end

struct PWLFunction{D, F, T <: Segment{D, F}}
    segments::Vector{T}
    meta::Dict

    function PWLFunction{D, F, T}(segments::Vector{T}) where {D, F, T <: Segment}
        return new(segments, Dict())
    end
end

const PWLFunctionPointRep{D, F} = PWLFunction{D, F, SegmentPointRep{D, F}}
const PWLFunctionSegmentRep{D, F} = PWLFunction{D, F, SegmentHyperplaneRep{D, F}}

const UnivariatePWLFunction{F} = PWLFunctionPointRep{1, F}
const BivariatePWLFunction{F} = PWLFunctionPointRep{2, F}
