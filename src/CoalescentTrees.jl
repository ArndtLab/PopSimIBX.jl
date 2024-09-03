module CoalescentTrees


export AbstractCoalescentTree, 
    SimpleCoalescentTree, CoalescentTree, 
    first_time, first_id, last_time, iscoalescent, time_span

abstract type AbstractCoalescentTree end


struct SimpleCoalescentTree{T} <: AbstractCoalescentTree
    timespan::T
end

time_span(ct::SimpleCoalescentTree{T}) where {T} = ct.timespan
iscoalescent(ct::SimpleCoalescentTree{T}) where {T} = true


"""
    CoalescentTree{T}

A coalescent tree. Object of this type are returned by the IBDTreeIterator.
"""
struct CoalescentTree{T} <: AbstractCoalescentTree
    ids::Vector{Int64}
    first_id::Int64
    first_time::Float64
    last_time::Float64
    tree::T

    function CoalescentTree{T}(ids::Vector{Int64}, firstid::Int64, firsttime::Float64, lasttime::Float64, tree::T) where T
        @assert length(ids) == 2 "The ids must be a vector of length 2 - more needs to be implemented"
        new(ids, firstid, firsttime, lasttime, tree)
    end
end

CoalescentTree(ids::Vector{Int64}, first_id::Int64, first_time::Float64, last_time::Float64) = CoalescentTree{Nothing}(ids, first_id, first_time, last_time, nothing)

first_time(ct::CoalescentTree{T}) where {T} = ct.first_time
last_time(ct::CoalescentTree{T}) where {T} = ct.last_time
time_span(ct::CoalescentTree{T}) where {T} = last_time(ct) - first_time(ct)
first_id(ct::CoalescentTree{T}) where {T}   = ct.first_id
iscoalescent(ct::CoalescentTree{T}) where {T} = first_id(ct) > 0

Base.show(io::IO, ct::CoalescentTree{T}) where {T} = print(io, "CoalescentTree starting at $(ct.first_time) in $(ct.first_id) for $(length(ct.ids)) individuals")



end    # module CoalescentTrees

