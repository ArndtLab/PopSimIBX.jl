using TestItems
using TestItemRunner
@run_package_tests  verbose = true




@testitem "SMC" begin
    using .PopSimBase

    for genome_length in [10, 1000, 1_000_000],
            mutation_rate in [1.0e-3, 1.0e-9, 1.0e-10]

        pop = StationaryPopulation(; genome_length, mutation_rate)

        h = sum(segment_length, SMC.IBDIterator(pop))
        @test h == pop.genome_length

        h = collect(segment_length.(SMC.IBDIterator(pop)))
        @test length(h) > 0
        @test sum(h) == pop.genome_length

        ibds = collect(SMC.IBDIterator(pop))
        ibas = collect(IBAIterator(ibds))
        ibss = collect(IBSIterator(ibds, pop.mutation_rate))
        ibms = collect(IBMIterator(ibds, pop.mutation_rate))
        ibmss = collect(IBSIterator(ibms))
        ibmss2= collect(IBSIterator(ibms))
        @test ibmss == ibmss2

        # @show length(ibds), length(ibas), length(ibss)

        @test genome_length == sum(segment_length, ibds)
        @test genome_length == sum(segment_length, ibas)
        @test genome_length == sum(segment_length, ibss)
        @test genome_length == sum(segment_length, ibms)
        @test genome_length == sum(segment_length, ibmss)

        @test genome_length == mapreduce(segment_length, +, SMC.IBDIterator(pop))
        @test genome_length == mapreduce(segment_length, +, IBAIterator(ibds))
        @test genome_length == mapreduce(segment_length, +, IBAIterator(SMC.IBDIterator(pop)))
        @test genome_length == mapreduce(segment_length, +, IBSIterator(SMC.IBDIterator(pop), pop.mutation_rate))
        @test genome_length == mapreduce(segment_length, +, IBMIterator(SMC.IBDIterator(pop), pop.mutation_rate))
        @test genome_length == mapreduce(segment_length, +, IBSIterator(IBMIterator(SMC.IBDIterator(pop), pop.mutation_rate)))
    end

end

@testitem "SMCprime" begin
    using PopSimIBX.CoalescentTrees
    # for genome_length  in [1_000, 1_000_000, 100_000],
    for genome_length in [1_000],
            mutation_rate in [1.0e-4, 1.0e-5, 1.0e-9],
            population_size in [1_000, 10_000, 100_000],
            rep in 1:3

        pop = StationaryPopulation(; genome_length, mutation_rate, population_size)


        ibds = collect(SMCprime.IBDIterator(pop))
        ibas = collect(IBAIterator(ibds))
        ibss = collect(IBSIterator(ibds, pop.mutation_rate))

        # @show length(ibds), length(ibas), length(ibss)

        meantau = mapreduce(+, ibds) do x
            time_span(PopSimIBX.Segments.data(x))
        end / length(ibds)
        # @show meantau

        @test genome_length == sum(segment_length, ibds)
        @test genome_length == sum(segment_length, ibas)
        @test genome_length == sum(segment_length, ibss)

        @test genome_length == mapreduce(segment_length, +, SMCprime.IBDIterator(pop))
        @test genome_length == mapreduce(segment_length, +, IBAIterator(ibds))
        @test genome_length == mapreduce(segment_length, +, IBAIterator(SMCprime.IBDIterator(pop)))
        @test genome_length == mapreduce(segment_length, +, IBSIterator(SMCprime.IBDIterator(pop), pop.mutation_rate))

        rtotal = mapreduce(+, ibss) do x
            PopSimIBX.Segments.data(x)
        end
        @test rtotal + 1 == length(ibds)

    end
end



@testitem "SMCprime const VaryingPopulation" begin
    using PopSimIBX.CoalescentTrees
    for genome_length in [10_000_000_000],
            population_size in [100_000],
            rep in 1:10

        pop = StationaryPopulation(; genome_length)
        ibds = collect(SMCprime.IBDIterator(pop))

        meantau1 = mapreduce(+, ibds) do x
            time_span(PopSimIBX.Segments.data(x))
        end / length(ibds)

        nepochs = 5
        Ns = fill(pop.population_size, nepochs)
        Ts = [i * meantau1 for i in 0:(nepochs - 1)]
        pop = VaryingPopulation(; genome_length, population_sizes = Ns, times = Ts)
        ibds = collect(SMCprime.IBDIterator(pop))
        @test length(ibds) > 1000

        meantau2 = mapreduce(+, ibds) do x
            time_span(PopSimIBX.Segments.data(x))
        end / length(ibds)
        @test abs(meantau1 - meantau2) / meantau1 < 0.05
    end
end

@testitem "Histogram append" begin
    using StatsBase
    using Random
    pop = StationaryPopulation(; genome_length = 1_000_000_000)

    h = Histogram(1:1000:1_000_001)
    si = SMCprime.IBDIterator(pop)
    append!(h, si)
    @test sum(h.weights) > 0

    h = Histogram(1:1000:1_000_001)
    si = IBAIterator(SMCprime.IBDIterator(pop))
    append!(h, si)
    @test sum(h.weights) > 0

    h = Histogram(1:1000:1_000_001)
    si = IBSIterator(SMCprime.IBDIterator(pop), pop.mutation_rate)
    append!(h, si)
    @test sum(h.weights) > 0

    # for i in 1:100
    #     local h = Histogram(1:1000:1_000_001)
    #     local si = collect(SMCprime.IBDIterator(pop))
    #     append!(h, si)
    #     @test sum(x -> segment_length(x), si) == pop.genome_length
    #     toolo = sum(x -> segment_length(x) == 0, si)
    #     toohi = sum(x -> segment_length(x) > 1_000_000, si)

    #     @test sum(h.weights) + toolo + toohi == length(si)
    # end

    for i in 1:10
        local pop = StationaryPopulation(; genome_length = 1_000_000, recombination_rate = 1.0e-6)
        local h = Histogram(1:1_000)
        Random.seed!(i)
        append!(h, SMCprime.IBDIterator(pop))

        Random.seed!(i)
        local si = collect(SMCprime.IBDIterator(pop))

        c1 = sum(x -> segment_length(x) == 1, si)
        c2 = sum(x -> segment_length(x) == 2, si)
        c3 = sum(x -> segment_length(x) == 3, si)
        c50 = sum(x -> segment_length(x) == 50, si)
        @test c1 == h.weights[1]
        @test c2 == h.weights[2]
        @test c3 == h.weights[3]
        @test c50 == h.weights[50]
    end
end
