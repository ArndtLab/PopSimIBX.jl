module PopSimIBX

using Reexport

@reexport using PopSimBase

import StatsBase: AbstractHistogram, Histogram

export Populations
export StationaryPopulation, VaryingPopulation

export SMC, SMCprime, Hudson


export IBAIterator, IBMIterator, IBSIterator, AbstractSegmentsIterator
export segment_length, time_span
export append!


# using PopSimBase.Populations   # This is necessary to get StationaryPopulation into the namespace


# using PopSimBase.Segmentals   # This is necessary to get Segmental into the namespace




include("SMCiterators.jl")
include("Hudson.jl")

include("Iterators.jl")
using .Iterators   # This is necessary to get IBAIterator into the namespace

include("utils.jl")

end    # module PopSimIBX
