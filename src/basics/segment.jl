export Segment

struct Segment{N, T}
    p₁::Point{N, T}
    p₂::Point{N, T}
end

Segment(ds::Segment{N, T}) where {N, T} = Segment{N, T}(ds.p₁, ds.p₂)

Base.eltype(::Segment{N, T}) where {N, T} = T

function Base.getindex(s::Segment, i::Integer)
    return if i == 1
        s.p₁
    elseif i == 2
        s.p₂
    else
        throw(BoundsError(s, i))
    end
end

Base.reverse(s::Segment) = Segment(s[1], s[2])
