export DirSegment

struct DirSegment{N,T}
    p₁::Point{N,T}
    p₂::Point{N,T}
end

DirSegment(s::DirSegment{N,T}) where {N,T} = DirSegment{N,T}(s.p₁, s.p₂)

Base.eltype(::DirSegment{N,T}) where {N,T} = T

function Base.getindex(s::DirSegment, i::Integer)
    if i == 1
        s.p₁
    elseif i == 2
        s.p₂
    else
        throw(BoundsError(s, i))
    end
end

Base.reverse(ds::DirSegment) = DirSegment(ds[2], ds[1])
