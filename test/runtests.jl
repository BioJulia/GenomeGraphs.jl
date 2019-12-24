module TestGenomeGraphs

using GenomeGraphs, BioSequences, Test

using GenomeGraphs.MerTools

# write your own tests here
@testset "MerCount" begin
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
    
    @test collapse_into_counts(v) == [
        MerCount{DNAMer{4}}(mer"AAAA", 3),
        MerCount{DNAMer{4}}(mer"AGGT", 4),
        MerCount{DNAMer{4}}(mer"ATAG", 1),
        MerCount{DNAMer{4}}(mer"GGGG", 2),
        MerCount{DNAMer{4}}(mer"GGGT", 1)
    ]
    
    @test collapse_into_counts(v2) == [
        MerCount{DNAMer{4}}(mer"AAAA", 3),
        MerCount{DNAMer{4}}(mer"AAGT", 2),
        MerCount{DNAMer{4}}(mer"ACGT", 1),
        MerCount{DNAMer{4}}(mer"AGGT", 1),
        MerCount{DNAMer{4}}(mer"ATAG", 1),
        MerCount{DNAMer{4}}(mer"GGGG", 2),
        MerCount{DNAMer{4}}(mer"GGGT", 1)
    ]
    
    @test test_merge_into(collapse_into_counts(v), collapse_into_counts(v2), [
        MerCount{DNAMer{4}}(mer"AAAA", 6),
        MerCount{DNAMer{4}}(mer"AAGT", 2),
        MerCount{DNAMer{4}}(mer"ACGT", 1),
        MerCount{DNAMer{4}}(mer"AGGT", 5),
        MerCount{DNAMer{4}}(mer"ATAG", 2),
        MerCount{DNAMer{4}}(mer"GGGG", 4),
        MerCount{DNAMer{4}}(mer"GGGT", 2)
    ])
end

end # module