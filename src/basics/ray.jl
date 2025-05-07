export Ray, RayF64,
    Line, LineF64

struct Ray{N, T}
    point::Point{N, T}
    vector::Vec{N, T}
    Ray(p::Point{N, T}, v::Vec{N, T}) where {N, T} = new{N, T}(p, normalize(v))
end
const RayF64 = Ray{N, Float64} where {N}

Ray(p₁::Point, p₂::Point) = Ray(p₁, p₂ - p₁)
Ray(s::Union{Segment, DirSegment}) = Ray(s[1], s[2] - s[1])

point(r::Ray, x) = if x ≥ 0
    r.point + r.vector * x
else
    throw(ArgumentError(lazy"`r` must be non-negative"))
end
point(r::Ray, x, dim) = if !iszero(r.vector[dim])
    point(r, (x - r.point[dim]) / r.vector[dim])
else
    point(r, 0)
end

struct Line{N, T}
    point::Point{N, T}
    vector::Vec{N, T}
    Line(p::Point{N, T}, v::Vec{N, T}) where {N, T} = new{N, T}(p, normalize(v))
end
const LineF64 = Line{N, Float64} where {N}

Line(p₁::Point, p₂::Point) = Line(p₁, p₂ - p₁)
Line(s::Union{Segment, DirSegment}) = Line(s[1], s[2] - s[1])

point(l::Line, x) = l.point + l.vector * x
point(l::Line, x, dim) = if !iszero(l.vector[dim])
    point(l, (x - l.point[dim]) / l.vector[dim])
else
    point(l, 0)
end
