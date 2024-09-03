# PopSimIBX

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ArndtLab.github.io/PopSimIBX.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ArndtLab.github.io/PopSimIBX.jl/dev/)
[![Build Status](https://github.com/ArndtLab/PopSimIBX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArndtLab/PopSimIBX.jl/actions/workflows/CI.yml?query=branch%3Amain)


Example

```julia
using PopSimIBX
using StatsBase


pop = StationaryPopulation(;
    popoulation_size = 1000,
    genome_length = 100_000_000, 
    mutation_rate = 1e-7, recombination_rate = 1e-7)


hist = Histogram(1:1000:1_000_001)
append!(hist, IBSIterator(SMCprime.IBDIterator(pop), pop.mutation_rate))
```




```julia
using PopSimIBX
using StatsBase

Ts = [0, 200, 400]
Ns = [1_000, 100, 1_000]

pop = VaryingPopulation(;
    genome_length = 100_000_000, 
    mutation_rate = 1e-7, recombination_rate = 1e-7,
    population_sizes = Ns,
    times = Ts
    )

hist = Histogram(1:1000:1_000_001)
append!(hist, IBSIterator(SMCprime.IBDIterator(pop), pop.mutation_rate))
```