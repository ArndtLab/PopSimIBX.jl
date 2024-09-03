module Segmentals

export Segmental, start, stop, data, segment_length,
    AbstractSegmentalsIterator

struct Segmental{T}
    start::Int64
    stop::Int64
    data::T
end


start(s::Segmental) = s.start
stop(s::Segmental) = s.stop
data(s::Segmental) = s.data
# segmentlength(s::Segmental) = s.stop - s.start + 1
segment_length(s::Segmental) = s.stop - s.start + 1


abstract type AbstractSegmentalsIterator end

end    # module Segmentals