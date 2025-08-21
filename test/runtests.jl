using TestItems
using TestItemRunner
@run_package_tests  verbose = true


@testitem "IBS kwargs" begin
    tau = 100.0
    tree = CoalescentTrees.SimpleCoalescentTree(tau)

    ibds = [SegItem(Segment(1, 1000), tree), SegItem(Segment(1001, 2000), tree)]
    ibs = collect(IBSIterator(ibds, 0.1))
    @test length(ibs) > 3000
    @test sum(segment_length, ibs) === 2000

    ibs = collect(IBSIterator(ibds, 0.1, multiple_hits = :as_one))
    @test length(ibs) in [1999,2000,2001]
    @test sum(segment_length, ibs) === 2000

    ibs = collect(IBSIterator(ibds, 0.1, multiple_hits = :JCcorrected))
    @test length(ibs) < 1750
    @test sum(segment_length, ibs) === 2000

end

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




@testitem "Hudson distribute" begin
    using PopSimIBX.Hudson
    bps = [10, 20, 30, 100, 110, 120, 130, 140, 150, 160]

    vi = [Segment(2, 4)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0

    vi = [Segment(1, 1)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0


    vi = [Segment(172, 175)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) + length(v2) == 1


    vi = [Segment(12, 14)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 1

    vi = [Segment(2, 4), Segment(12, 14)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1


    vi = [Segment(2, 10)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0

    vi = [Segment(11, 14)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 1

    vi = [Segment(2, 10), Segment(11, 14)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1

    vi = [Segment(2, 14)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1


    vi = [Segment(2, 24)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 2 && length(v2) == 1

    vi = [Segment(2, 34)]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 2 && length(v2) == 2

    vi = Segment{Int}[]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 0

    vi = [Segment(2, 34)]
    v1, v2 = Hudson.distribute(vi, Int[])
    @test length(v1) == 1 && length(v2) == 0

end


@testitem "Hudson distribute rate" begin
    using PopSimIBX.Hudson

    vi = [Segment(1, 1000000)]
    v1, v2 = Hudson.distribute(vi, 1e-3)

    @test length(v1) > 0
    @test length(v2) > 0
    @test sum(segment_length, v1) + sum(segment_length, v2) == 1000000
end


@testitem "Hudson coalesce" begin
    using PopSimIBX.Hudson


    t = -1    
    vc = Vector{SegItem{Int, CoalescentTrees.SimpleCoalescentTree}}()
    v1 = [Segment(1, 2), Segment(3, 4)]
    v2 = [Segment(5, 6), Segment(7, 8)]

    v = Hudson.coalesce(v1, v2, vc, t)
    @test length(v) == 4 && length(vc) == 0

    v1 = [Segment(1, 2), Segment(3, 4)]
    v2 = [Segment(15, 16), Segment(17, 18)]

    v = Hudson.coalesce(v1, v2, empty!(vc), t)
    @test length(v) == 4 && length(vc) == 0
    @test vcat(v1, v2) == v
    
    v1 = [Segment(1, 20)]
    v2 = [Segment(4,6)] 
    v = Hudson.coalesce(v1, v2, empty!(vc), t)
    @test length(v) == 2 && length(vc) == 1
    @test start(vc[1]) == 4
    @test stop(vc[1]) == 6
    @test data(vc[1]).timespan ==  1
    
    v = Hudson.coalesce(v2, v1, empty!(vc), t)
    @test length(v) == 2 && length(vc) == 1


    v1 = [Segment(1, 20)]
    v2 = [Segment(4,4)] 
    v = Hudson.coalesce(v1, v2, empty!(vc), t)
    @test length(v) == 2 && length(vc) == 1

    v = Hudson.coalesce(v2, v1, empty!(vc), t)
    @test length(v) == 2 && length(vc) == 1

    v1 = [Segment(1, 2), Segment(3, 4)]
    v2 = [Segment(2, 3), Segment(4, 5)]
    v = Hudson.coalesce(v1, v2, empty!(vc), t)
    @test length(v) == 2 && length(vc) == 3

    v1 = [Segment(1, 2), Segment(3, 4)]
    v = Hudson.coalesce(v1, v1, empty!(vc), t)
    @test length(v) == 0 && length(vc) == 2

    v1 = [Segment(1,1), Segment(3, 3)]
    v = Hudson.coalesce(v1, v1, empty!(vc), t)
    @test length(v) == 0 && length(vc) == 2
end




@testitem "coalesce loop" begin
    using PopSimIBX.Hudson
    t = 1
    
    a1 = 100
    for l1 in 1:5
        for l2 in 1:5
            for a2 in a1-l2-1 : a1 + l1 + 1
                v1 = [Segment(a1, a1 + l1)]
                v2 = [Segment(a2, a2 + l2)]
                vc = Vector{SegItem{Int, CoalescentTrees.SimpleCoalescentTree}}()
                v = Hudson.coalesce(v1, v2, vc, t)
                @test length(vc) <= 1
                @test length(v) <= 2
                @test length(v) + length(vc) >= 1
            end
        end
    end
end



@testitem "Hudson distribute & coalesce" begin
    using PopSimIBX.Hudson

    v1 = [Segment(1, 100)]
    v2 = [Segment(1, 100)]
    vc = Vector{SegItem{Int, CoalescentTrees.SimpleCoalescentTree}}()


    bps = [30,50,70,90]
    v11, v12 = Hudson.distribute(v1, bps)

    bps = [10,50,60]
    v21, v22 = Hudson.distribute(v2, bps)

    t = 1
    v3 = Hudson.coalesce(v11, v21, vc, t)
    v4 = Hudson.coalesce(v12, v22, vc, t)

    t = 2
    v5 = Hudson.coalesce(v3, v4, vc, t)

    sort!(vc, by = start)
    @test length(v5) == 0
    @test sum(segment_length, vc) == 100

end


@testitem "Hudson StationaryPopulation" begin
    using PopSimIBX.Hudson


    for genome_length in [1000, 10000, 1000000],
            population_size in [100, 1000],
            recombination_rate in [1.0e-8, 1.0e-9],
            mutation_rate in [1.0e-9, 1.0e-10]
            
        pop = StationaryPopulation(; genome_length, recombination_rate, population_size)

        ibds = Hudson.IBDIterator(pop)
        @test sum(segment_length, ibds) == genome_length

        @test sum(segment_length, IBAIterator(ibds)) == genome_length
        @test sum(segment_length, IBSIterator(ibds, pop.mutation_rate)) == genome_length

    end
end



@testitem "Hudson VaryingPopulation" begin
    using PopSimIBX.Hudson


    for genome_length in [1000, 10000, 1000000],
            population_size in [100, 1000],
            recombination_rate in [1.0e-8, 1.0e-9],
            mutation_rate in [1.0e-9, 1.0e-10]
            

        ps = [population_size*2, population_size ÷ 2, population_size]
        ts = [0.0, 100, 200]
        pop = VaryingPopulation(; genome_length, recombination_rate, population_sizes = ps, times = ts)
        @test length(pop.times) == length(pop.population_sizes)

        ibds = Hudson.IBDIterator(pop)
        @test sum(segment_length, ibds) == genome_length

        @test sum(segment_length, IBAIterator(ibds)) == genome_length
        @test sum(segment_length, IBSIterator(ibds, pop.mutation_rate)) == genome_length

    end
end

@testitem "HudsonMulti" begin
    h1 = Hudson.HudsonARG{Int}(1, 0)
    h2 = Hudson.HudsonARG{Int}(1, 0)
    @show h1 h2

    h = Hudson.HudsonARG{Int}(1, 3, h1, h2)

    @show h


    genome_length = 1000
    population_size = 100
    recombination_rate = 1.0e-8
    mutation_rate = 1.0e-9
    pop = StationaryPopulation(; genome_length, recombination_rate, population_size)

    ibds = Hudson.IBDIteratorMulti(pop, 4)
    @show ibds
end




@testitem "Hudson Multi distribute" begin
    using PopSimIBX.Hudson
    bps = [10, 20, 30, 100, 110, 120, 130, 140, 150, 160]

    vi = [SegItem(Segment(2, 4), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0


    vi = [SegItem(Segment(2, 4), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0

    vi = [SegItem(Segment(1, 1), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0


    vi = [SegItem(Segment(172, 175), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) + length(v2) == 1


    vi = [SegItem(Segment(12, 14), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 1

    vi = [SegItem(Segment(2, 4), Hudson.HudsonARG{Int}(1, 0)), SegItem(Segment(12, 14), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1


    vi = [SegItem(Segment(2, 10), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 0

    vi = [SegItem(Segment(11, 14), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 1

    vi = [SegItem(Segment(2, 10), Hudson.HudsonARG{Int}(1, 0)), SegItem(Segment(11, 14), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1

    vi = [SegItem(Segment(2, 14), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 1 && length(v2) == 1


    vi = [SegItem(Segment(2, 24), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 2 && length(v2) == 1

    vi = [SegItem(Segment(2, 34), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 2 && length(v2) == 2

    vi = similar([SegItem(Segment(2, 34), Hudson.HudsonARG{Int}(1, 0))], 0)
    v1, v2 = Hudson.distribute(vi, bps)
    @test length(v1) == 0 && length(v2) == 0

    vi = [SegItem(Segment(2, 34), Hudson.HudsonARG{Int}(1, 0))]
    v1, v2 = Hudson.distribute(vi, Int[])
    @test length(v1) == 1 && length(v2) == 0

end


@testitem "Hudson Multi coalesce" begin
    using PopSimIBX.Hudson


    t = 1 
    
    h1 = Hudson.HudsonARG{Int}(1, 0)
    h2 = Hudson.HudsonARG{Int}(2, 0)
    h3 = Hudson.HudsonARG{Int}(3, 0)
   
    v1 = [SegItem(Segment(1, 2), h1), SegItem(Segment(3, 4), h2)]
    v2 = [SegItem(Segment(5, 6), h3), SegItem(Segment(7, 8), h1)]


    vc = similar(v1, 0)
    v = Hudson.coalesce(v1, v2, vc, t, 3)
    @test length(v) == 4 && length(vc) == 0
   
    
    v1 = [SegItem(Segment(1, 2), h1)]
    v2 = [SegItem(Segment(1, 2), h2)]
    v3 = [SegItem(Segment(1, 2), h3)]
    vc = similar(v1, 0)
    v4 = Hudson.coalesce(v1, v2, vc, t, 3)
    @test length(v4) == 1 && length(vc) == 0

    t += 1
    v5 = Hudson.coalesce(v4, v3, vc, t, 3)
    @test length(v5) == 0 && length(vc) == 1

end
   


@testitem "Hudson Multi StationaryPopulation" begin
    using PopSimIBX.Hudson


    for genome_length in [1000, 10000, 1000000],
            population_size in [100, 1000],
            recombination_rate in [1.0e-8, 1.0e-9],
            mutation_rate in [1.0e-4, 1.0e-10]
            
        pop = StationaryPopulation(; genome_length, recombination_rate, population_size)

        ibds = Hudson.IBDIteratorMulti(pop, 4)
        @test sum(segment_length, ibds) == genome_length

        ibdm = IBMIterator(ibds, pop.mutation_rate)
        @test sum(segment_length, ibdm) == genome_length

        ibs = IBSIterator(ibdm, 1, 2)
        @test sum(segment_length, ibs) == genome_length

        ibs = IBSIterator(ibdm, 1, 4)
        @test sum(segment_length, ibs) == genome_length

    end

end

@testitem "Hudson Multi createbranches" begin
    using PopSimIBX.Hudson

    h1 = Hudson.HudsonARG{Int}(1, 0)
    h2 = Hudson.HudsonARG{Int}(2, 0)
    h3 = Hudson.HudsonARG{Int}(3, 0)
    h4 = Hudson.HudsonARG{Int}(4, 0)
    h5 = Hudson.HudsonARG{Int}(5, 0)
    h6 = Hudson.HudsonARG{Int}(6, 6, h2, h3)
    h7 = Hudson.HudsonARG{Int}(7, 7, h4, h5)
    h8 = Hudson.HudsonARG{Int}(8, 8, h6, h1)
    h9 = Hudson.HudsonARG{Int}(9, 9, h7, h8)



    n = 5
    ids = collect(1:n)
    first_id = 2 * n - 1 
    idc = first_id
    first_time = float(h9.time)
    last_time = 0.0 
    branches = fill((pid = -1, time = 0.0, vid_anc = -1), idc)

    nextinternal = idc - 1
    parent_idc = -1
    Hudson.createbranches!(h9, branches, nextinternal, idc, parent_idc)

    tree = CoalescentTrees.CoalescentTree(ids, first_id, first_time, last_time, branches)
    @show tree
    @show branches
    @test 1 == 1 
end

