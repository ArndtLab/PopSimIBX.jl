module PopSimIBX

import StatsBase: AbstractHistogram

export Populations
export StationaryPopulation, VaryingPopulation

export SMC, SMCprime


export IBAIterator, IBSIterator, AbstractSegmentalsIterator
# export Segmentals
export segment_length, time_span
export append!


include("Populations.jl")
using .Populations   # This is necessary to get StationaryPopulation into the namespace


include("CoalescentTrees.jl")
include("Segmentals.jl")
using .Segmentals   # This is necessary to get Segmental into the namespace

include("SMCiterators.jl")


include("Iterators.jl")
using .Iterators   # This is necessary to get IBAIterator into the namespace

include("utils.jl")

end    # module PopSimIBX
