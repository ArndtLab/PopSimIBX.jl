module Iterators



export IBAIterator, IBMIterator, IBSIterator, IBSIteratorMutated

using PopSimBase
using Distributions





mutable struct IBAIterator{T} <: AbstractSegmentalsIterator
    ibds::T
end


Base.IteratorSize(::Type{IBAIterator{T}}) where {T} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBAIterator{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBAIterator{T}}) where {T} =
    Segmentals.Segmental{Nothing,CoalescentTrees.SimpleCoalescentTree{Int64}}

function Base.iterate(ai::IBAIterator{T}) where {T}
    f = iterate(ai.ibds)
    iterate(ai, f)
end

function Base.iterate(ai::IBAIterator{T}, state) where {T}
    isnothing(state) && return nothing

    seg = state[1]
    aiibdsstate = state[2]

    mystart = start(seg)
    mystop = stop(seg)

    # lens = [segment_length(seg)]

    local next
    while true
        # println("state: ", state)
        next = iterate(ai.ibds, aiibdsstate)

        isnothing(next) && break
        aiibdsstate = next[2]
        cond = data(seg) == data(next[1])
        # cond = all(f -> f(seg.data) ≈ f(next[1].data), [time_span])

        # @show seg.data, next[1].data, seg.data == next[1].data, cond
        if cond
            # push!(lens, segment_length(next[1]))
            mystop = stop(next[1])
        else
            break
        end
    end

    Segmental(mystart, mystop, data(seg)), next
end

# -----------------------------------------------------------------------------


function break_segment(L::Int64, prob::Float64, allow_multiple_hits = true)
    if prob == Inf
        return collect(1:L)
    end
    n_breaks = rand(Poisson(L * prob))
    n_breaks == 0 && return Int64[]
    if allow_multiple_hits
        r = rand(1:L, n_breaks)
        return unique!(sort!(r))
    else
        n_breaks >= L && return collect(1:L)
        r = sample(1:L, n_breaks, replace = false)
        return sort!(r)
    end
end




mutable struct IBSIterator{T} <: AbstractSegmentalsIterator
    ibxs::T
    breaks::Vector{Int64}
    lastibxstop::Int64
    mutation_rate::Float64
end

IBSIterator(ibx, mutation_rate) = IBSIterator(ibx, Int[], 0, mutation_rate)


Base.IteratorSize(::Type{IBSIterator{T}}) where {T} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIterator{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBSIterator{T}}) where {T} = Segmentals.Segmental{Nothing,Int64}


mutable struct IBSIteratorState{T}
    ibxstate::Union{Nothing,T}
    laststop::Int64
    b::Int64
end



function Base.iterate(si::IBSIterator)
    ibx = iterate(si.ibxs)
    isnothing(ibx) && return nothing # no ibx to iterate

    seg = ibx[1]

    dt = CoalescentTrees.time_span(data(seg))
    si.breaks =
        (start(seg) - 1) .+ break_segment(segment_length(seg), 2 * si.mutation_rate * dt)
    si.lastibxstop = stop(seg)
    state = IBSIteratorState(ibx[2], 0, 1)
    iterate(si, state)
end


function Base.iterate(si::IBSIterator, state)
    isnothing(state) && return nothing

    mystart = state.laststop + 1
    # @show mystart, state.b, si.breaks

    if state.b <= length(si.breaks)
        mystop = si.breaks[state.b]
        state.b += 1
        state.laststop = mystop
        return Segmental(mystart, mystop, 0), state
    else
        r = 0
        while true
            ibx = iterate(si.ibxs, state.ibxstate)
            if isnothing(ibx) # last ibd reached, emit last interval
                return Segmental(mystart, si.lastibxstop, r), nothing
            end

            state.ibxstate = ibx[2]
            seg = ibx[1]
            r += 1
            # println("seg: ", seg)

            dt = CoalescentTrees.time_span(data(seg))
            si.breaks =
                (start(seg) - 1) .+
                break_segment(segment_length(seg), 2 * si.mutation_rate * dt)
            si.lastibxstop = stop(seg)

            if isempty(si.breaks) # no breaks in ibx
                continue
            end

            # @show si.breaks
            mystop = si.breaks[1]
            state.b = 2
            state.laststop = mystop

            return Segmental(mystart, mystop, r), state
        end
    end
end

# -----------------------------------------------------------------------------


function sprinckle_mutations(s::Segmental{}, mutation_rate::Float64)

    slen = segment_length(s)
    branches = data(s).tree

    breaks = map(branches) do branch
        r = Int64[]
        if branch.vid_anc > 0
            dtime = branch.time - branches[branch.vid_anc].time 
            r = break_segment(slen, mutation_rate * dtime)
        end
        r
    end

    d = CoalescentTrees.MutatedCoalescentTree(data(s).ids, data(s).first_id, 
        data(s).first_time, data(s).last_time, 
        branches, slen, breaks)
    Segmental(s.start, s.stop, d)
end

function IBMIterator(iter, mutation_rate)
    Iterators.map(s->sprinckle_mutations(s, mutation_rate),iter)
end

# -----------------------------------------------------------------------------


mutable struct IBSIteratorMutated{T}
    ibxs::T
    vids::Vector{Int64}
    breaks::Vector{Int64}
    lastibxstop::Int64
end

IBSIteratorMutated(ibx, vids::Vector{Int64}) = IBSIteratorMutated(ibx, vids, Int64[], 0)
IBSIteratorMutated(ibx, vid1, vid2) = IBSIteratorMutated(ibx, [vid1, vid2], Int64[], 0)


function get_breaks(ct::CoalescentTrees.MutatedCoalescentTree, vids::Vector{Int64})
    breaks = Int64[]
    vs = sort(vids)
    while true
        if (length(vs) >= 2) && (vs[1] == vs[2])
            vs = vs[2:end]
        end
        if length(vs) == 1
            break
        end
        vid_anc = ct.tree[vs[1]].vid_anc
        if vid_anc > 0
            breaks = vcat(breaks, ct.mutations[vs[1]])
            vs[1] = vid_anc
        else
            vs = vs[2:end]
        end
        sort!(vs)
    end
    unique!(sort!(breaks))
end

function get_breaks(ct::CoalescentTrees.MutatedSimpleCoalescentTree, vids::Vector{Int64})
    ct.breaks
end





Base.IteratorSize(::Type{IBSIteratorMutated{T}}) where {T} = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIteratorMutated{T}}) where {T} = Base.HasEltype()
Base.eltype(::Type{IBSIteratorMutated{T}}) where {T} = Segmentals.Segmental{Nothing,Int64}


mutable struct IBSIteratorMutatedState{T}
    ibxstate::Union{Nothing,T}
    laststop::Int64
    b::Int64
end



function Base.iterate(si::IBSIteratorMutated)
    ibx = iterate(si.ibxs)
    isnothing(ibx) && return nothing # no ibx to iterate

    seg = ibx[1]
    all(0 .<= si.vids .<= length(data(seg).ids)) || throw(ArgumentError("The vids must be a vector of integers between 0 and $(length(data(seg).ids))"))


    si.breaks = get_breaks(data(seg), si.vids)
    si.lastibxstop = stop(seg)
    state = IBSIteratorMutatedState(ibx[2], 0, 1)
    iterate(si, state)
end


function Base.iterate(si::IBSIteratorMutated, state)
    isnothing(state) && return nothing

    mystart = state.laststop + 1
    # @show mystart, state.b, si.breaks

    if state.b <= length(si.breaks)
        mystop = si.breaks[state.b]
        state.b += 1
        state.laststop = mystop
        return Segmental(mystart, mystop, 0), state
    else
        r = 0
        while true
            ibx = iterate(si.ibxs, state.ibxstate)
            if isnothing(ibx) # last ibd reached, emit last interval
                return Segmental(mystart, si.lastibxstop, r), nothing
            end

            state.ibxstate = ibx[2]
            seg = ibx[1]
            r += 1
            # println("seg: ", seg)

            si.breaks =
                (start(seg) - 1) .+ get_breaks(data(seg), si.vids)
            si.lastibxstop = stop(seg)

            if isempty(si.breaks) # no breaks in ibx
                continue
            end

            # @show si.breaks
            mystop = si.breaks[1]
            state.b = 2
            state.laststop = mystop

            return Segmental(mystart, mystop, r), state
        end
    end
end


allIBS(collection, args...) = _f(collection, Base.IteratorEltype(collection), args...)
_f(collection, ::Base.EltypeUnknown, args...) = throw(ArgumentError("The eltype of the collection is unknown"))
_f(collection, ::Base.HasEltype, args...) = 
    (println("eltype = $(eltype(collection))"); _f(collection, eltype(collection), args...))



_f(collection, ::Type{T}, args...) where {T<:Segmental{CoalescentTrees.MutatedSimpleCoalescentTree}} = 0
_f(collection, ::Type{T}, args...) where {T<:Segmental{CoalescentTrees.MutatedCoalescentTree}} = 0


end
