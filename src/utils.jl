


function Base.append!(h::AbstractHistogram{T,1}, si::Union{AbstractSegmentalsIterator, Vector{Segmentals.Segmental{S}}}, 
    func = segment_length) where {T, S}

    @assert h.closed == :left
    edges = h.edges[1]

    lo = edges[1]
    hi = edges[end] - 1
    

    # collect the counts
    countv = zeros(Int, hi)
    for s in si
        l = func(s)
        if 1 <= l <= hi
            countv[l] += 1
        end
    end

    # add the counts to the histogram
    k = 1
    for i in eachindex(countv)
        iszero(countv[i]) && continue
        i < edges[k] && continue
        while (k + 1 <= length(edges)) && !(edges[k] <= i < edges[k+1])
            k += 1
        end
        if k <= length(edges) - 1
            h.weights[k] += countv[i]
        end
    end
    h    
end

