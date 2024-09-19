export Ray, RayF64

struct Ray{N,T}
    point::Point{N,T}
    vector::Vec{N,T}
    Vec{N,T}(p, v) where {N,T} = new(p, normalize!(v))
end
const RayF64 = Ray{N,Float64} where {N}

struct Line{N,T}
    point::Point{N,T}
    vector::Vec{N,T}
    Vec{N,T}(p, v) where {N,T} = new(p, normalize!(v))
end
const LineF64 = Line{N,Float64} where {N}


Line(p₁::Point{N}, p₂::Point{N}) where {N} = Line{N}(p₁, p₂-p₁)
Line(s::Union{Segment{N},DirSegment{N}}) where {N} = Line{N}(s[1], s[2]-s[1])
