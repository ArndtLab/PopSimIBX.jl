module Hudson 

using PopSimBase
using Distributions

function nextinterval(vis::Iterators.Stateful, posmax::T) where {T <: Integer}
    if isempty(vis)
        return posmax, posmax
    end
    i = popfirst!(vis)
    return start(i), stop(i)
end



function distribute(vi::Vector{Segment{T}}, bps::Vector{T}) where {T <: Integer}
    v1 = Vector{Segment{T}}()
    v2 = Vector{Segment{T}}()
    length(vi) == 0 && return (v1, v2)

    vis = Iterators.Stateful(vi)

    length(bps) == 0 && return (vi, v2)
    nextpushv1 = true

    bpss = Iterators.Stateful(sort(bps))
    nextbppos = popfirst!(bpss)

    i = popfirst!(vis)
    nextintervalstart = i.start
    nextintervalstop = i.stop

    k = 1
    pos = min(nextintervalstop, nextbppos)
    posmax = vi[end].stop + 1
    while pos < posmax
        # @show (k, pos) (nextintervalstart, nextintervalstop, nextbppos)

        if (pos == nextintervalstop || pos == nextbppos) && nextintervalstart <= pos
            if nextpushv1
                push!(v1, Segment(nextintervalstart, pos))
            else
                push!(v2, Segment(nextintervalstart, pos))
            end
            nextintervalstart = pos + 1
        end
        if pos == nextintervalstop
            nextintervalstart, nextintervalstop = nextinterval(vis, posmax)
        end
        if pos == nextbppos
            nextbppos = isempty(bpss) ? posmax : popfirst!(bpss)
            nextpushv1 = !nextpushv1
        end
        # @show (nextintervalstart, nextintervalstop, nextbppos)
            
        # @show v1 v2

        pos = min(nextintervalstop, nextbppos)
        k += 1
        # k > 20 && break
    end

    v1, v2
end

function coalesce(
        v1::Vector{Segment{T}}, 
        v2::Vector{Segment{T}}, 
        vc::Vector{SegItem{T, CoalescentTrees.SimpleCoalescentTree}},
        tau::F) where {T <: Integer, F <: Real}


    length(v1) == 0 && return (v2, vc)
    length(v2) == 0 && return (v1, vc)
    
    v = Vector{Segment{T}}()

    posmax = max(v1[end].stop, v2[end].stop) + 1
    
    v1s = Iterators.Stateful(v1)
    v2s = Iterators.Stateful(v2)
    next1start, next1stop = nextinterval(v1s, posmax)
    next2start, next2stop = nextinterval(v2s, posmax)

    pos = min(next1start, next2start)

    while pos < posmax
     
        if pos == next1start && next2start == pos
            nextpos = min(next1stop, next2stop)
        elseif pos == next1start && next2start > pos
            nextpos = min(next1stop, next2start)
        elseif pos == next2start && next1start > pos
            nextpos = min(next2stop, next1start)
        elseif pos == next1start && next2start < pos
            push!(v, Segment(next2start, pos - 1))
            next2start = pos
            nextpos = min(next1stop, next2stop)
        elseif pos == next2start && next1start < pos
            push!(v, Segment(next1start, pos - 1))
            next1start = pos
            nextpos = min(next1stop, next2stop)
        end

        if pos == next1stop && pos == next2stop
            @assert next1start == next2start
            tree = CoalescentTrees.SimpleCoalescentTree(tau)
            push!(vc, SegItem(Segment(next1start, pos), tree))
            next1start, next1stop = nextinterval(v1s, posmax)
            next2start, next2stop = nextinterval(v2s, posmax)
            nextpos = min(next1start, next2start)
        elseif pos == next1stop && next2start > pos
            push!(v, Segment(next1start, pos))
            next1start, next1stop = nextinterval(v1s, posmax)
            nextpos = min(next1start, next2start)
        elseif pos == next1stop && next2stop > pos
            @assert next1start == next2start
            tree = CoalescentTrees.SimpleCoalescentTree(tau)
            push!(vc, SegItem(Segment(next1start, pos), tree))
            next1start, next1stop = nextinterval(v1s, posmax)
            next2start = pos + 1
            nextpos = min(next1start, next2stop)
        elseif pos == next2stop && next1start > pos
            push!(v, Segment(next2start, pos))
            next2start, next2stop = nextinterval(v2s, posmax)
            nextpos = min(next1start, next2start)
        elseif pos == next2stop && next1stop > pos
            @assert next1start == next2start
            tree = CoalescentTrees.SimpleCoalescentTree(tau)
            push!(vc, SegItem(Segment(next1start, pos), tree))
            next2start, next2stop = nextinterval(v2s, posmax)
            next1start = pos + 1
            nextpos = min(next1stop, next2start)
        end

        pos = nextpos
    end
   
    return v
end



function IBDIterator(anc::StationaryPopulation)

    L = genome_length(anc)
    N = anc.population_size

    T = Int 
    v1 = [ [Segment{T}(1, L)], [Segment{T}(1, L)] ]
    v2 = similar(v1, 0)
    vc = Vector{SegItem{T, CoalescentTrees.SimpleCoalescentTree}}()

    t = 1
    while true
        empty!(v2)

        for vi in v1
            n_breaks = rand(Poisson(L * anc.recombination_rate))
            bps = sort!(rand(1:L, n_breaks))
            vi1, vi2 = distribute(vi, bps)

            if !isempty(vi1)
                k = rand(1: 2 * N)
                if k > length(v2)
                    push!(v2, vi1)
                else
                    v2[k] = coalesce(v2[k], vi1, vc, t)
                end
            end

            if !isempty(vi2)
                k = rand(1: 2 * N)
                if k > length(v2)
                    push!(v2, vi2)
                else
                    v2[k] = coalesce(v2[k], vi2, vc, t)
                end
            end
        end

        v1 = filter(x -> length(x) > 0, v2)
        # lc = sum(vci -> segment_length(vci), vc, init = 0)
        # n = length(v1)
        # l = sum(v -> sum(segment_length, v), v1, init = 0)
        # @assert 2 * lc + l == 2* L
        if isempty(v1)
            break
        end
        
        t += 1

    end
    sort!(vc, by = start)
    @assert sum(segment_length, vc) == L
    return vc
end



end # module Hudson