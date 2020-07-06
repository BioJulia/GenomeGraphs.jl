
function remove_tips!(sdg::Graphs.SDG, min_size::Integer)
    @info "Beginning tip removal process"
    #sdg = graph(ws)
    pass = 1
    tips = Graphs.find_tip_nodes(sdg, min_size)
    ntips = length(tips)
    utgs = Vector{Graphs.SequenceDistanceGraphPath{typeof(sdg)}}()
    newnodes = Vector{Graphs.NodeID}()
    while true
        @info string("Pass number: ", pass)
        if isempty(tips)
            @info "No more tips in graph"
            break
        end
        @info string("Found ", ntips, " tips")
        for tip in tips
            Graphs.remove_node!(sdg, tip)
        end
        @info string("Removed ", ntips, " tips")
        @info "Collapsing any resulting transient paths"
        Graphs.collapse_all_unitigs!(utgs, newnodes, sdg, 2, true)
        pass = pass + 1
        tips = Graphs.find_tip_nodes(sdg, min_size)
        ntips = length(tips)
    end
    @info string("Finished tip removal process in ", pass, " passes")
    #return ws
    return sdg
end