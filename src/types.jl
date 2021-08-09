struct PWLFunction{D}
    x::Vector{NTuple{D,Float64}}
    z::Vector{Float64}
    T::Vector{Vector{Int}}
    meta::Dict
end

function PWLFunction{D}(x::Vector{NTuple{D}}, z::Vector, T::Vector{Vector}, meta::Dict) where {D}
    @assert length(x) == length(z)
    for t in T
        @assert minimum(t) > 0 && maximum(t) <= length(x)
    end
    return PWLFunction{D}(x, z, T, meta)
end

PWLFunction(x, z, T) = PWLFunction(x, z, T, Dict())

const UnivariatePWLFunction = PWLFunction{1}

function UnivariatePWLFunction(x, z)
    @assert issorted(x)
    return PWLFunction(Tuple{Float64}[(xx,) for xx in x], convert(Vector{Float64}, z), [[i,i+1] for i in 1:length(x)-1])
end

function UnivariatePWLFunction(x, fz::Function)
    @assert issorted(x)
    return PWLFunction(Tuple{Float64}[(xx,) for xx in x], map(t->convert(Float64,fz(t)), x), [[i,i+1] for i in 1:length(x)-1])
end

const BivariatePWLFunction = PWLFunction{2}

function BivariatePWLFunction(x, y, fz::Function; pattern=:BestFit, seed=hash((length(x),length(y))))
    @assert issorted(x)
    @assert issorted(y)
    X = vec(collect(Base.product(x,y)))
    # X = vec(Tuple{Float64,Float64}[(_x,_y) for _x in x, _y in y])
    Z = map(t -> convert(Float64,fz(t...)), X)
    T = Vector{Vector{Int}}()
    m = length(x)
    n = length(y)

    mt = MersenneTwister(seed)
    # run for each square on [x[i],x[i+1]] × [y[i],y[i+1]]
    for i in 1:length(x)-1, j in 1:length(y)-1
        SWt, NWt, NEt, SEt = LinearIndices((m,n))[i,j], LinearIndices((m,n))[i,j+1], LinearIndices((m,n))[i+1,j+1], LinearIndices((m,n))[i+1,j]
        xL, xU, yL, yU = x[i], x[i+1], y[j], y[j+1]
        @assert xL == X[SWt][1] == X[NWt][1]
        @assert xU == X[SEt][1] == X[NEt][1]
        @assert yL == X[SWt][2] == X[SEt][2]
        @assert yU == X[NWt][2] == X[NEt][2]
        SW, NW, NE, SE = Z[SWt], Z[NWt], Z[NEt], Z[SEt]
        mid1 = 0.5*(SW+NE)
        mid2 = 0.5*(NW+SE)

        if pattern == :Upper
            if mid1 > mid2
                t1 = [SWt,NWt,NEt]
                t2 = [SWt,NEt,SEt]
            else
                t1 = [SWt,NWt,SEt]
                t2 = [SEt,NWt,NEt]
            end
        elseif pattern == :Lower
            if mid1 > mid2
                t1 = [SWt,NWt,SEt]
                t2 = [SEt,NWt,NEt]
            else
                t1 = [SWt,NWt,NEt]
                t2 = [SWt,NEt,SEt]
            end
        elseif pattern == :BestFit
            mid3 = fz(0.5*(xL+xU), 0.5*(yL+yU))
            if abs(mid1-mid3) < abs(mid2-mid3)
                t1 = [SWt,NWt,NEt]
                t2 = [SWt,NEt,SEt]
            else
                t1 = [SWt,NWt,SEt]
                t2 = [SEt,NWt,NEt]
            end
        elseif pattern == :UnionJack
            t1 = [SWt,SEt]
            t2 = [NWt,NEt]
            if iseven(i+j)
                push!(t1, NWt)
                push!(t2, SEt)
            else
                push!(t1, NEt)
                push!(t2, SWt)
            end
        elseif pattern == :K1
            t1 = [SEt,SWt,NWt]
            t2 = [NWt,NEt,SEt]
        elseif pattern == :Random
            if rand(mt, Bool)
                t1 = [NWt,NEt,SEt]
                t2 = [SEt,SWt,NWt]
            else
                t1 = [SWt,NWt,NEt]
                t2 = [NEt,SEt,SWt]
            end
        else
            error("pattern $pattern not currently supported")
        end

        push!(T, t1)
        push!(T, t2)
    end

    return PWLFunction{2}(X, Z, T, Dict(:structure=>pattern))
end
