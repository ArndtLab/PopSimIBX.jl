function Base.append!(
        h::AbstractHistogram{T, 1},
        si::AbstractSegmentalsIterator,
        func = segment_length
    ) where {T}

    # collect the counts
    for s in si
        val = func(s)
        push!(h, val)
    end
    h
end

