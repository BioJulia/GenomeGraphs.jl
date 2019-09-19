module TestGenomeGraphs

using GenomeGraphs, BioSequences, Test

import GenomeGraphs.MerFreq

# write your own tests here
@testset "MerFreq" begin
    v = DNAMer{4}[
        mer"AAAA", mer"AAAA", mer"AAAA", mer"ATAG",
        mer"GGGG", mer"GGGG", mer"GGGT", mer"AGGT",
        mer"AGGT", mer"AGGT", mer"AGGT"
    ]
    
    v2 = DNAMer{4}[
        mer"AAAA", mer"AAAA", mer"AAAA", mer"ATAG",
        mer"GGGG", mer"GGGG", mer"GGGT", mer"AGGT",
        mer"ACGT", mer"AAGT", mer"AAGT"
    ]
    
    function test_merge_into(a, b, a′)
        merge_into!(a, b)
        return isempty(b) && a == a′
    end
    
    @test collapse_into_freqs(v) == [
        MerFreq{DNAMer{4}}(mer"AAAA", 3),
        MerFreq{DNAMer{4}}(mer"AGGT", 4),
        MerFreq{DNAMer{4}}(mer"ATAG", 1),
        MerFreq{DNAMer{4}}(mer"GGGG", 2),
        MerFreq{DNAMer{4}}(mer"GGGT", 1)
    ]
    
    @test collapse_into_freqs(v2) == [
        MerFreq{DNAMer{4}}(mer"AAAA", 3),
        MerFreq{DNAMer{4}}(mer"AAGT", 2),
        MerFreq{DNAMer{4}}(mer"ACGT", 1),
        MerFreq{DNAMer{4}}(mer"AGGT", 1),
        MerFreq{DNAMer{4}}(mer"ATAG", 1),
        MerFreq{DNAMer{4}}(mer"GGGG", 2),
        MerFreq{DNAMer{4}}(mer"GGGT", 1)
    ]
    
    @test test_merge_into(collapse_into_freqs(v), collapse_into_freqs(v2), [
        MerFreq{DNAMer{4}}(mer"AAAA", 6),
        MerFreq{DNAMer{4}}(mer"AAGT", 2),
        MerFreq{DNAMer{4}}(mer"ACGT", 1),
        MerFreq{DNAMer{4}}(mer"AGGT", 5),
        MerFreq{DNAMer{4}}(mer"ATAG", 2),
        MerFreq{DNAMer{4}}(mer"GGGG", 4),
        MerFreq{DNAMer{4}}(mer"GGGT", 2)
    ])
end

end # module