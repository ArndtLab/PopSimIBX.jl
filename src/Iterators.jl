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




mutable struct IBSIterator{T} <: AbstractSegmentalsIterator
    ibxs::T
    mutation_rate::Float64
end

IBSIterator(ibx, mutation_rate) = IBSIterator(ibx, Int[], 0, mutation_rate)


Base.IteratorSize(::Type{IBSIterator{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBSIterator{T}}) where T = Base.HasEltype()
Base.eltype(::Type{IBSIterator{T}}) where {T} = Segmentals.Segmental{Int64}


mutable struct IBSIteratorState{T}
    ibxstate::Union{Nothing, T}
    lastibxstop::Int64
    lastibxdt::Real
    laststop::Int64
end



function Base.iterate(si::IBSIterator)
    ibx = iterate(si.ibxs)
    isnothing(ibx) && return nothing # no ibx to iterate

    seg = ibx[1]

    dt = time_span(data(seg))
    state = IBSIteratorState(ibx[2], stop(seg), dt,  0)
    iterate(si, state)
end


function randnextstop(start, p)
    if p < 1
        return start + rand(Geometric(p)) 
    else 
        return start 
    end
end



function Base.iterate(si::IBSIterator, state)
    isnothing(state) && return nothing

    mystart = state.laststop + 1
    nextstop = randnextstop(mystart, 2 * si.mutation_rate * state.lastibxdt) # first break


    if nextstop <= state.lastibxstop
        mystop = nextstop
        state.laststop = mystop
        return Segmental(mystart, mystop, 0), state  # 0 recombination events in between mutations
    else
        r = 0
        while true
            ibx = iterate(si.ibxs, state.ibxstate)
            if isnothing(ibx) # last ibd reached, emit last interval
                return Segmental(mystart, state.lastibxstop, r), nothing
            end
            r += 1

            seg = ibx[1]
            state.ibxstate = ibx[2]
            state.lastibxstop = stop(seg)
            state.lastibxdt = time_span(data(seg))
            nextstop = randnextstop(start(seg) - 1, 2 * si.mutation_rate * state.lastibxdt) 
            
            if nextstop > state.lastibxstop
                continue
            end

            mystop = nextstop
            state.laststop = mystop 

            return Segmental(mystart, mystop, r), state
        end
    end
end



end
