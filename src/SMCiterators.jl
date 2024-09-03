module SMC


using Distributions

using ..CoalescentTrees
using ..Segmentals
using ..Populations





mutable struct IBDIterrator <: AbstractSegmentalsIterator
    anc::StationaryPopulation
    tau_recombination::Int
end

IBDIterrator(anc::StationaryPopulation) = IBDIterrator(anc, 1)


Base.IteratorSize(::Type{IBDIterrator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBDIterrator}) = Base.HasEltype()
Base.eltype(::Type{IBDIterrator}) = Segmentals.Segmental{CoalescentTrees.SimpleCoalescentTree{Int64}}


function Base.iterate(ti::IBDIterrator, pos = 1)
    if pos > ti.anc.genome_length
        return nothing
    end
    
    tau = rand(Geometric(1/(2 * ti.anc.population_size))) + ti.tau_recombination
    ti.tau_recombination = rand(DiscreteUniform(1, tau))
    len = rand(Geometric(2 * ti.anc.recombination_rate * tau))

    tree = SimpleCoalescentTree(tau) # for two individuals no actual tree is returned

    stop = min(pos + len - 1, ti.anc.genome_length)    
    return Segmental(pos, stop, tree), stop + 1
end



end # module




module SMCprime


using Distributions

using ..CoalescentTrees
using ..Segmentals
using ..Populations


mutable struct IBDIterrator{T} <: AbstractSegmentalsIterator
    anc::T
    tau_recombination::Int
    tau_previous::Int
end

IBDIterrator(anc::StationaryPopulation) = IBDIterrator(anc, 1, 1)
IBDIterrator(anc::VaryingPopulation) = IBDIterrator(anc, 1, 1)


Base.IteratorSize(::Type{IBDIterrator{T}}) where T = Base.SizeUnknown()
Base.IteratorEltype(::Type{IBDIterrator{T}}) where T = Base.HasEltype()
Base.eltype(::Type{IBDIterrator{T}}) where {T} = Segmentals.Segmental{CoalescentTrees.SimpleCoalescentTree{Int64}}



function Base.iterate(ti::IBDIterrator{StationaryPopulation}, pos = 1)
    if pos > ti.anc.genome_length
        return nothing
    end

    tau_back = rand(Geometric(1/ti.anc.population_size))
    if (ti.tau_recombination + tau_back) > ti.tau_previous
        tau = rand(Geometric(1/(2 * ti.anc.population_size))) + ti.tau_previous
    else
        equal = rand(Bool)
        if equal
            tau = ti.tau_previous
        else
            tau = ti.tau_recombination + tau_back
        end
    end
    ti.tau_previous = tau
    ti.tau_recombination = rand(DiscreteUniform(1, tau))
    len = rand(Geometric(2 * ti.anc.recombination_rate * tau))

    tree = SimpleCoalescentTree(tau) # for two individuals no actual tree is returned

    stop = min(pos + len - 1, ti.anc.genome_length)    
    return Segmental(pos, stop, tree), stop + 1
end



function Base.iterate(ti::IBDIterrator{VaryingPopulation}, pos = 1)
    if pos > ti.anc.genome_length
        return nothing
    end

    (; tau_recombination, tau_previous) = ti
    (; population_sizes, times) = ti.anc
    
    epoch = findlast(tau_recombination .>= times)
    tau_back = tau_recombination + rand(Geometric(1/population_sizes[epoch]))
    while (epoch < length(times)) && (tau_back >= times[epoch+1]) && (tau_previous > times[epoch+1])
        epoch += 1
        tau_back = times[epoch] + rand(Geometric(1/population_sizes[epoch]))
    end
    if tau_back > tau_previous
        epoch_pr = findlast(tau_previous .>= times)
        tau = tau_previous + rand(Geometric(1/(2 * population_sizes[epoch_pr])))
        while (epoch_pr<length(times)) && (tau >= times[epoch_pr+1])
            epoch_pr += 1
            tau = times[epoch_pr] + rand(Geometric(1/(2 * population_sizes[epoch_pr])))
        end
    else
        equal = rand(Bool)
        if equal
            tau = tau_previous
        else
            tau = tau_back
        end
    end
    ti.tau_previous = ceil(Int, tau)
    ti.tau_recombination = rand(DiscreteUniform(1, ti.tau_previous))
    len = rand(Geometric(2 * ti.anc.recombination_rate * tau))


    tree = SimpleCoalescentTree(ti.tau_previous) # for two individuals no actual tree is returned

    stop = min(pos + len - 1, ti.anc.genome_length)    
    return Segmental(pos, stop, tree), stop + 1
end



end # module



