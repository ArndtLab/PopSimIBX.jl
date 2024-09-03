
module Populations

export StationaryPopulation, VaryingPopulation

struct StationaryPopulation
    ploidy::Int
    genome_length::Int
    recombination_rate::Float64
    mutation_rate::Float64
    population_size::Int
end

function StationaryPopulation(; 
        ploidy = 2, 
        population_size = 1_000, 
        genome_length = 1_000_000, 
        recombination_rate = 1e-8, 
        mutation_rate = 1e-8)
    StationaryPopulation(ploidy, genome_length, recombination_rate, mutation_rate, population_size)    
end


struct VaryingPopulation
    ploidy::Int
    genome_length::Int
    recombination_rate::Float64
    mutation_rate::Float64

    population_sizes::Vector{Int}
    times::Vector{Float64}
end


function VaryingPopulation(; 
        ploidy = 2, 
        population_sizes = [1_000], 
        times = [0.0], 
        genome_length = 1_000_000, 
        recombination_rate = 1e-8, 
        mutation_rate = 1e-8)
    if length(population_sizes) != length(times)
        throw(ArgumentError("population_sizes and times must have the same length"))
    end
    if !issorted(times)
        throw(ArgumentError("times must be sorted"))
    end
    if times[1] != 0.0
        throw(ArgumentError("times must start at 0.0"))
    end
    VaryingPopulation(ploidy, genome_length, recombination_rate, mutation_rate, population_sizes, times)
end

end    # module Populations
