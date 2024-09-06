


function Base.append!(h::AbstractHistogram{T,1}, 
    si::Union{AbstractSegmentalsIterator, Vector{Segmentals.Segmental{S}}}, 
    func = segment_length) where {T, S}

    # collect the counts
    for s in si
        val = func(s)
        push!(h, val)
    end
    h    
end

