# PopSimIBX

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ArndtLab.github.io/PopSimIBX.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ArndtLab.github.io/PopSimIBX.jl/dev/)
[![Build Status](https://github.com/ArndtLab/PopSimIBX.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/ArndtLab/PopSimIBX.jl/actions/workflows/CI.yml?query=branch%3Amain)


Example

```julia
using PopSimIBX
using StatsBase


pop = StationaryPopulation(;genome_length = 100_000_000, mutation_rate = 1e-7, recombination_rate = 1e-7)


hist = Histogram(1:1000:1_000_001)
append!(hist, IBSIterator(SMCprime.IBDIterator(pop), pop.mutation_rate))


```