using TestItems
using TestItemRunner
@run_package_tests  verbose=true



@testitem "StationaryPopulation" begin


    genome_length = 1_000_000

    pop = StationaryPopulation(;genome_length)
    # pop = Populations.StationaryPopulation(;genome_length)
    
    @test pop.genome_length == genome_length
end

@testitem "VaryingPopulation" begin
    genome_length = 1_000_000
    Ts = [0.0, 1.0, 2.0]
    Ns = [1_000, 1_000, 1_000]
    pop = VaryingPopulation(;genome_length, population_sizes = Ns, times = Ts)

    @test pop.genome_length == genome_length
    @test pop.population_sizes == Ns
    @test pop.times == Ts
end



@testitem "SMC" begin
    for genome_length  in [10, 1000, 1_000_000],
        mutation_rate in [1e-3, 1e-9, 1e-10]

        pop = StationaryPopulation(;genome_length, mutation_rate)

        ibds = collect(SMC.IBDIterator(pop))
        ibas = collect(IBAIterator(ibds))
        ibss = collect(IBSIterator(ibds, pop.mutation_rate))

        # @show length(ibds), length(ibas), length(ibss)

        @test genome_length == sum(segment_length, ibds)
        @test genome_length == sum(segment_length, ibas)
        @test genome_length == sum(segment_length, ibss)

        @test genome_length == mapreduce(segment_length, +, SMC.IBDIterator(pop))
        @test genome_length == mapreduce(segment_length, +, IBAIterator(ibds))
        @test genome_length == mapreduce(segment_length, +, IBAIterator(SMC.IBDIterator(pop)))
        @test genome_length == mapreduce(segment_length, +, IBSIterator(SMC.IBDIterator(pop), pop.mutation_rate))
    end

end

@testitem "SMCprime" begin
    using PopSimIBX.CoalescentTrees
    # for genome_length  in [1_000, 1_000_000, 100_000],
    for genome_length  in [1_000],
        mutation_rate in [1e-4, 1e-5, 1e-9],
        population_size in [1_000, 10_000, 100_000],
        rep = 1:3

        pop = StationaryPopulation(;genome_length, mutation_rate, population_size)

        
        ibds = collect(SMCprime.IBDIterator(pop))
        ibas = collect(IBAIterator(ibds))
        ibss = collect(IBSIterator(ibds, pop.mutation_rate))

        # @show length(ibds), length(ibas), length(ibss)
        
        meantau = mapreduce(+, ibds) do x
            time_span(PopSimIBX.Segmentals.data(x))
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
            PopSimIBX.Segmentals.data(x)
        end
        @test rtotal + 1 == length(ibds)

    end
end


@testitem "SMCprime const VaryingPopulation" begin
    using PopSimIBX.CoalescentTrees
    for genome_length  in [10_000_000_000],
        population_size in [100_000],
        rep = 1:10

        pop = StationaryPopulation(;genome_length)
        ibds = collect(SMCprime.IBDIterator(pop))
        
        meantau1 = mapreduce(+, ibds) do x
            time_span(PopSimIBX.Segmentals.data(x))
        end / length(ibds)

        nepochs = 5
        Ns = fill(pop.population_size, nepochs)
        Ts = [i * meantau1 for i in 0:nepochs-1]
        pop = VaryingPopulation(;genome_length, population_sizes = Ns, times = Ts)
        ibds = collect(SMCprime.IBDIterator(pop))
        @test length(ibds) > 1000

        meantau2 = mapreduce(+, ibds) do x
            time_span(PopSimIBX.Segmentals.data(x))
        end / length(ibds)
        @test abs(meantau1 - meantau2)/meantau1 < 0.05
    end
end

@testitem "Histogram" begin
    using StatsBase
    pop = StationaryPopulation(;genome_length = 1_000_000_000)
    
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

    for i in 1:100
        local h = Histogram(1:1000:1_000_001)
        local si = collect(SMCprime.IBDIterator(pop))
        append!(h, si)
        @test sum(x->segment_length(x), si) == pop.genome_length
        toolo = sum(x -> segment_length(x)==0, si)
        toohi = sum(x -> segment_length(x)>1_000_000, si)

        @test sum(h.weights) + toolo + toohi == length(si)
    end

    for i in 1:10
        local pop = StationaryPopulation(;genome_length = 1_000_000, recombination_rate = 1e-6)
        local h = Histogram(1:1_000)
        local si = collect(SMCprime.IBDIterator(pop))
        append!(h, si)
        c1 = sum(x->segment_length(x)==1, si)
        c2 = sum(x->segment_length(x)==2, si)
        c3 = sum(x->segment_length(x)==3, si)
        c50 = sum(x->segment_length(x)==50, si)
        @test c1 == h.weights[1]
        @test c2 == h.weights[2]
        @test c3 == h.weights[3]
        @test c50 == h.weights[50]
    end
end