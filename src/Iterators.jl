module Iterators



export IBAIterator, IBSIterator

using ..Segmentals
using ..CoalescentTrees
using Distributions





mutable struct IBAIterator{T} <: AbstractSegmentalsIterator
    ibds::T
end


Base.IteratorSize(::Type{IBAIterator{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBAIterator{T}}) where T = Base.HasEltype()
Base.eltype(::Type{IBAIterator{T}}) where {T} = Segmentals.Segmental{CoalescentTrees.SimpleCoalescentTree{Int64}}

function Base.iterate(ai::IBAIterator{T}) where T
    f = iterate(ai.ibds)
    iterate(ai, f)
end

function Base.iterate(ai::IBAIterator{T}, state) where T
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
        return unique(sort(rand(1:L, n_breaks)))
    else
        n_breaks >=L && return collect(1:L)
        return sort(sample(1:L, n_breaks, replace = false))
    end
end




mutable struct IBSIterator{T} <: AbstractSegmentalsIterator
    ibxs::T
    breaks::Vector{Int64}
    lastibxstop::Int64
    mutation_rate::Float64
end

IBSIterator(ibx, mutation_rate) = IBSIterator(ibx, Int[], 0, mutation_rate)


Base.IteratorSize(::Type{IBSIterator{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIterator{T}}) where T = Base.HasEltype()
Base.eltype(::Type{IBSIterator{T}}) where {T} = Segmentals.Segmental{Int64}


mutable struct IBSIteratorState{T}
    ibxstate::Union{Nothing, T}
    laststop::Int64
    b::Int64
end



function Base.iterate(si::IBSIterator)
    ibx = iterate(si.ibxs)
    isnothing(ibx) && return nothing # no ibx to iterate

    seg = ibx[1]

    dt = time_span(data(seg))
    si.breaks = (start(seg) - 1) .+ break_segment(segment_length(seg), 2 * si.mutation_rate * dt)
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

            dt = time_span(data(seg))
            si.breaks = (start(seg) - 1) .+ break_segment(segment_length(seg), 2 * si.mutation_rate * dt)
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






end
