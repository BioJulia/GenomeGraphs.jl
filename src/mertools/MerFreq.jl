
"""
    Kmer frequency type
"""
struct MerFreq{M<:AbstractMer}
    mer::M
    count::UInt8
    function MerFreq{M}(mer::M, count::Integer) where {M<:AbstractMer}
        return new(mer, convert(UInt8, min(typemax(UInt8), count)))
    end
end

@inline mer(x::MerFreq{M}) where {M<:AbstractMer} = x.mer
@inline freq(x::MerFreq{M}) where {M<:AbstractMer} = x.count
@inline freq(::Type{R}, x::MerFreq{M}) where {R<:Real,M<:AbstractMer} = convert(R, freq(x))

const DNAMerFreq{K} = MerFreq{DNAMer{K}}
const RNAMerFreq{K} = MerFreq{RNAMer{K}}

Base.isless(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer} = mer(x) < mer(y)
Base.:(>)(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer} = mer(x) > mer(y)
#Base.:(==)(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer} = mer(x) == mer(y)

function merge(x::MerFreq{M}, y::MerFreq{M}) where {M<:AbstractMer}
    return MerFreq{M}(x.mer, UInt16(x.count) + UInt16(y.count))
end

function Base.show(io::IO, mfreq::MerFreq{<:AbstractMer})
    print(io, mer(mfreq), " occurs ", freq(mfreq), " times")
end

"""
    merge_into_sorted!(a::Vector{MerFreq{M}}, b::Vector{MerFreq{M}}) where {M<:AbstractMer}

Collapse the 
"""
function merge_into_sorted!(a::Vector{MerFreq{M}}, b::Vector{MerFreq{M}}) where {M<:AbstractMer}
    a_i = firstindex(a)
    a_end = lastindex(a) + 1
    b_i = b_i2 = firstindex(b)
    b_end = lastindex(b) + 1
    
    # Merge, accumulating counts on `a`.
    @inbounds while b_i < b_end
        while a_i < a_end && a[a_i] < b[b_i]
            a_i = a_i + 1
        end
        if a_i < a_end && mer(a[a_i]) == mer(b[b_i])
            # Combine entries
            a[a_i] = merge(a[a_i], b[b_i])
            a_i = a_i + 1
            b_i = b_i + 1
        end
        while b_i < b_end && (a_i == a_end || b[b_i] < a[a_i])
            b[b_i2] = b[b_i]
            b_i2 = b_i2 + 1
            b_i = b_i + 1
        end
        
    end
    # Shrink `b` to the size of the remaining contents.
    resize!(b, b_i2 - 1)
    
    # Expand `a` to allow the insertion of unique values in `b`.
    oldsize = length(a)
    resize!(a, oldsize + length(b))
    r_a = oldsize
    
    # Merge-sort from the bottom into `a`.
    wr_a = lastindex(a)
    rend_a = firstindex(a)
    r_b = lastindex(b)
    r_end_b = firstindex(b)
    @inbounds while wr_a >= rend_a
        if r_b >= r_end_b && (r_a < rend_a || b[r_b] > a[r_a])
            a[wr_a] = b[r_b]
            r_b = r_b - 1
        else
            a[wr_a] = a[r_a]
            r_a = r_a - 1
        end
        wr_a = wr_a - 1
    end
    empty!(b)
    return a
end

function merge_into!(a::Vector{MerFreq{M}}, b::Vector{MerFreq{M}}) where {M<:AbstractMer}
    sort!(a)
    sort!(b)
    return merge_into_sorted!(a, b)
end

function collapse_sorted!(freqs::Vector{MerFreq{M}}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    pi = 1
    stop = lastindex(freqs) + 1
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mer(freqs[ri]) == mer(freqs[pi])
            ci = ci + one(UInt16)
        end
        freqs[wi] = MerFreq{M}(mer(freqs[pi]), ci)
        pi = ri
        wi = wi + 1
    end
    resize!(freqs, wi - 1)
    return freqs
end

collapse!(freqs::Vector{MerFreq{M}}) where {M<:AbstractMer} = collapse_sorted!(sort!(freqs))

function collapse_into_freqs_sorted!(mers::Vector{M}, result::Vector{MerFreq{M}}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    stop = lastindex(mers) + 1
    empty!(result)
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mers[ri] == mers[wi]
            ci = ci + one(UInt16)
        end
        push!(result, MerFreq{M}(mers[wi], ci))
        wi = ri
    end
    return result
end

function collapse_into_freqs!(mers::Vector{M}, result::Vector{MerFreq{M}}) where {M<:AbstractMer}
    sort!(mers)
    return collapse_into_freqs_sorted!(mers, result)
end

function collapse_into_freqs_sorted(mers::Vector{M}) where {M<:AbstractMer}
    return collapse_into_freqs_sorted!(mers, Vector{MerFreq{M}}())
end

function collapse_into_freqs(mers::Vector{M}) where {M<:AbstractMer}
    return collapse_into_freqs!(mers, Vector{MerFreq{M}}())
end

struct MerFreqHist
    data::Vector{UInt64}
    min::UInt8
end

function Base.summary(io::IO, hist::MerFreqHist)
    print(io, "Frequency histogram of motifs appearing more than ", hist.min, "times")
end

function hist(freqs::Vector{MerFreq{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    hist = zeros(UInt64, 256)
    for x in freqs
        f = freq(x)
        if f ≥ min_count
            old = hist[f]
            hist[f] = old + 1
        end
    end
    return MerFreqHist(hist, convert(UInt8, min_count))
end

function hist!(freqs::Vector{MerFreq{M}}, min_count::Integer = 0) where {M<:AbstractMer}
    hist = zeros(UInt64, 256)
    wi = firstindex(freqs)
    used = 0
    for x in freqs
        f = freq(x)
        if f ≥ min_count
            old = hist[f]
            hist[f] = old + 1
            freqs[wi] = x
            wi = wi + 1
            used = used + 1
        end
        resize!(freqs, used)
    end
    return MerFreqHist(hist, convert(UInt8, min_count))
end
