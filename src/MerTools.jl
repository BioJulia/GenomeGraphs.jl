module MerTools

export
    # The MerCount type.
    MerCount,
    DNAMerCount,
    RNAMerCount,
    MerCountHist,
    mer,
    freq,
    hist,
    hist!,
    collapse_into_counts,
    collapse_into_counts!,
    merge_into!


using BioSequences
using ReadDatastores

import BioSequences: ksize

###
### Counting Mers
###

"""
A simple mer count struct.

MerCount is a simple struct that binds a mer value to a count of the number
of times it has been observed. This type, (sorted) vectors of them, and some
additional utility methods, form the basic building blocks of the higher-level
mer counting functionality of the MerTools sub-module. 

!!! note
    The count is stored as an UInt8 because often once the count is more than
    255 we hardly care anymore.
"""
struct MerCount{M<:AbstractMer}
    mer::M
    count::UInt8
    function MerCount{M}(mer::M, count::Integer) where {M<:AbstractMer}
        return new(mer, convert(UInt8, min(typemax(UInt8), count)))
    end
end

"Get the mer from a `MerCount`."
@inline mer(x::MerCount{M}) where {M<:AbstractMer} = x.mer

"Get the count from a `MerCount`."
@inline freq(x::MerCount{M}) where {M<:AbstractMer} = x.count
"Get the count from a `MerCount`, and convert it to type R."
@inline freq(::Type{R}, x::MerCount{M}) where {R<:Real,M<:AbstractMer} = convert(R, freq(x))

"Shorthand for `MerCount{DNAMer{K}}`"
const DNAMerCount{K} = MerCount{DNAMer{K}}

"Shorthand for `MerCount{RNAMer{K}}`"
const RNAMerCount{K} = MerCount{RNAMer{K}}

Base.isless(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer} = mer(x) < mer(y)
Base.:(>)(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer} = mer(x) > mer(y)
#Base.:(==)(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer} = mer(x) == mer(y)

function merge(x::MerCount{M}, y::MerCount{M}) where {M<:AbstractMer}
    return MerCount{M}(x.mer, UInt16(x.count) + UInt16(y.count))
end

function Base.show(io::IO, mfreq::MerCount{<:AbstractMer})
    print(io, mer(mfreq), " occurs ", freq(mfreq), " times")
end

"""
    unsafe_collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}

!!! warning
    This method is marked as unsafe because it assumes that the `mers` input
    vector is already sorted.
"""
function unsafe_collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    stop = lastindex(mers) + 1
    empty!(result)
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mers[ri] == mers[wi]
            ci = ci + one(UInt16)
        end
        push!(result, MerCount{M}(mers[wi], ci))
        wi = ri
    end
    return result
end

"""
    collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}

Build a vector of sorted `MerCount`s from a Vector of a mer type.

This is a basic kernel function used for any higher level and more complex
kmer counting procedures.

This is like `collapse_into_counts`, except it's first argument is a `result`
vector that is cleared and filled with the result.

!!! note
    The input vector `mers` will be sorted by this method.
"""
function collapse_into_counts!(result::Vector{MerCount{M}}, mers::Vector{M}) where {M<:AbstractMer}
    sort!(mers)
    return unsafe_collapse_into_counts!(result, mers)
end

"""
    collapse_into_counts(mers::Vector{M}) where {M<:AbstractMer}

Build a vector of sorted `MerCount`s from a Vector of a mer type.

This is a basic kernel function used for any higher level and more complex
kmer counting procedures.
"""
function collapse_into_counts(mers::Vector{M}) where {M<:AbstractMer}
    return collapse_into_counts!(Vector{MerCount{M}}(), mers)
end

"""
    unsafe_merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}

Merge the `MerCount`s from vector `b` into the vector `a`.

!!! warning
    This method is marked as unsafe as it assumes both of the input vectors `a`
    and `b` are already sorted.
"""
function unsafe_merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}
    a_i = firstindex(a)
    a_end = lastindex(a) + 1
    b_i = b_i2 = firstindex(b)
    b_end = lastindex(b) + 1
    
    # Merge, accumulating counts on `a`.
    @inbounds while b_i < b_end
        # Move cursor a_i forward so long as mers in a are less than the mer
        # in b[b_i]
        while a_i < a_end && a[a_i] < b[b_i]
            a_i = a_i + 1
        end
        # Check if the mers at the a_i and b_i cursors are the same. If so,
        # the count in b needs to be merged into the count in a.
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

"""
    merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}

Merge the `MerCount`s from vector `b` into the vector `a`.

!!! note
    This will sort the input vectors `a` and `b`.
"""
function merge_into!(a::Vector{MerCount{M}}, b::Vector{MerCount{M}}) where {M<:AbstractMer}
    sort!(a)
    sort!(b)
    return unsafe_merge_into!(a, b)
end

function unsafe_collapse!(freqs::Vector{MerCount{M}}) where {M<:AbstractMer}
    wi = 1
    ri = 1
    pi = 1
    stop = lastindex(freqs) + 1
    @inbounds while ri < stop
        ci = one(UInt16)
        while (ri += 1) < stop && mer(freqs[ri]) == mer(freqs[pi])
            ci = ci + one(UInt16)
        end
        freqs[wi] = MerCount{M}(mer(freqs[pi]), ci)
        pi = ri
        wi = wi + 1
    end
    resize!(freqs, wi - 1)
    return freqs
end

collapse!(freqs::Vector{MerCount{M}}) where {M<:AbstractMer} = collapse_sorted!(sort!(freqs))

###
### Mer count histogram (kmer spectra) type.
###

"""
A type for storing a frequency histogram of MerCounts, also sometimes referred
to as a Kmer spectra. 
"""
struct MerCountHist
    data::Vector{UInt64}
    min::UInt8
end

function Base.summary(io::IO, hist::MerCountHist)
    print(io, "Count histogram of motifs appearing more than ", hist.min, "times")
end

"""
    hist!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

Build a histogram of kmer frequencies, excluding any kmer counts that don't meet
`min_count`.
"""
function hist(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
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

"""
    hist!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}

A version of `hist` that builds a MerCountHist, AND mutates the input `freqs` by
filtering it to remove MerCounts that did not meet `min_count` and so were not
included in the constructed MerCountHist.
"""
function hist!(freqs::Vector{MerCount{M}}, min_count::Integer = 0) where {M<:AbstractMer}
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
    return MerCountHist(hist, convert(UInt8, min_count))
end


"""
    build_freq_list(::Type{M}, sbuf::DatastoreBuffer{PairedReads}, range::UnitRange{Int}) where {M<:AbstractMer}

Build a sorted list (vector) of kmer counts (MerFreq), serially and in memory.

This function is a serial and in memory `MerFreq` list builder that can build a
kmer count from a PairedReads datastore on its own (if you have memory and time),
but it is also intended to be composed into other multi-process or multi-threaded
kmer counting strategies.

This method estimates roughly how many kmers will be generated by the reads
specified by `range` in the dataset. It then pre-allocates an array to contain
them. It then collects the kmers, sorts, them, and then collapses them into a
list of counts sorted by the kmer.
"""
function build_freq_list(::Type{M}, sbuf::DatastoreBuffer{<:PairedReads}, range::UnitRange{Int}) where {M<:AbstractMer}
    max_read_size = max_read_length(ReadDatastores.datastore(sbuf))
    chunk_mers = Vector{M}(undef, length(range) * (max_read_size - ksize(M) + 1))
    wi = firstindex(chunk_mers)
    read_sequence = eltype(sbuf)()
    @inbounds for i in range
        for mer in each(M, load_sequence!(sbuf, i, read_sequence))
            chunk_mers[wi] = canonical(mer)
            wi = wi + 1
        end
    end
    resize!(chunk_mers, wi - 1)
    return collapse_into_counts(chunk_mers)
end

"""
    build_freq_list(::Type{M}, sbuf::SequenceBuffer{PairedReads}, range::UnitRange{Int}, chunk_size::Int) where {M<:AbstractMer}

Build a sorted list (vector) of kmer counts (MerFreq), serially and in memory.

This function is a serial and in memory `MerFreq` list builder that can build a
kmer count from a PairedReads datastore on it own (if you have memory and time),
but it is also intended to be composed into other multi-process or multi-threaded
kmer counting strategies.

This method pre-allocates space for `chunk_size` kmers, and iterates over kmers
in the reads in the dataset specified by `range` until the buffer is filled.
The mers are then collapsed into a list of counts, sorted by the kmer.
This list is then merged into another output list.
This process repeats for many chunks of kmers, building up the output list.

This method is useful for situations where you don't want (or have the space) to
allocate a buffer to collect all the kmers in the dataset all in one go.
"""
function build_freq_list(::Type{M}, sbuf::DatastoreBuffer{PairedReads}, range::UnitRange{Int}, chunk_size::Int) where {M<:AbstractMer}
    chunk_mers = Vector{M}(undef, chunk_size)
    chunk_freqs = Vector{MerCount{M}}()
    sizehint!(chunk_freqs, chunk_size)
    output = Vector{MerCount{M}}()
    
    read = first(range)
    lastread = last(range)
    read_sequence = eltype(sbuf)()
    
    mergen = each(M, load_sequence!(sbuf, read, read_sequence))
    mernext = iterate(mergen)
    chunkfill = 1
    
    @inbounds while read <= lastread
        # Fill the `chunk_mers` buffer with a load of kmers.
        while mernext !== nothing && chunkfill <= chunk_size
            chunk_mers[chunkfill] = canonical(mernext[1])
            chunkfill = chunkfill + 1
            mernext = iterate(mergen, mernext[2])
        end
        if isnothing(mernext) # Current read is consumed.
            read = read + 1
            if read <= lastread # There are more reads to come.
                mergen = each(M, load_sequence!(sbuf, read, read_sequence))
                mernext = iterate(mergen)
            else # There are no more reads to come.
                resize!(chunk_mers, chunkfill - 1)
                collapse_into_freqs!(chunk_mers, chunk_freqs)
                merge_into_sorted!(output, chunk_freqs)
                return output
            end
        end
        if chunkfill > chunk_size # The buffer is full, time to merge into output.
            collapse_into_freqs!(chunk_mers, chunk_freqs)
            merge_into_sorted!(output, chunk_freqs)
            chunkfill = 1
        end
    end
end



end # module




