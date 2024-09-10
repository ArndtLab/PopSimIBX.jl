function Base.append!(
        h::AbstractHistogram{T, 1},
        si::Union{AbstractSegmentalsIterator, Vector{Segmentals.Segmental{S}}},
        func = segment_length
    ) where {T, S}

    # collect the counts
    for s in si
        val = func(s)
        push!(h, val)
    end
    h
end


function multi_threaded_append!(
        h::AbstractHistogram{T, 1},
        si::Union{AbstractSegmentalsIterator, Vector{Segmentals.Segmental{S}}},
        len::Int,
        func = segment_length
    ) where {T, S}

    hs = map(1:Threads.nthreads()) do i
        Histogram(deepcopy(h.edges), Int)
    end

    chunk_len = cld(len, Threads.nthreads())
    sis = map(1:Threads.nthreads()) do i
        max_length = i == Threads.nthreads() ? len - chunk_len * (Threads.nthreads() - 1) : chunk_len
        PopSimIBX.MaxLenIterator(deepcopy(si), max_length)
    end

    tasks = map(1:Threads.nthreads()) do i
        Threads.@spawn begin
            for s in sis[i]
                push!(hs[i], func(s))
            end
        end
    end
    fetch.(tasks)

    merge!(h, hs...)
end
