include("basics/point.jl")
include("basics/vec.jl")

Vec(p::Point{N,T}) where {N,T} = Vec(p.coords)
Point(v::Vec{N,T}) where {N,T} = Point(v.coords)

include("basics/segment.jl")
include("basics/dirsegment.jl")

DirSegment(s::Segment{N,T}) where {N,T} = DirSegment{N,T}(s.p₁, s.p₂)
Segment(ds::DirSegment{N,T}) where {N,T} = Segment{N,T}(ds.p₁, ds.p₂)

include("basics/polygon.jl")

include("basics/dcel_table.jl")

export signapprox, orientation, isintersect, slope


signapprox(x::Number; tol=0) = ifelse(x < -tol, -1, ifelse(x > tol, 1, 0))


# Math and LinearAlgebra

Base.:+(p::Point{N}, v::Vec{N}) where {N} = Point(p.coords + v.coords)
Base.:+(v₁::Vec{N}, v₂::Vec{N}) where {N} = Vec(v₁.coords + v₂.coords)
Base.:-(p₁::Point{N}, p₂::Point{N}) where {N} = Vec(p₁.coords - p₂.coords)
Base.:-(p::Point{N}) where {N} = Vec(-p.coords)
Base.:-(v₁::Vec{N}, v₂::Vec{N}) where {N} = Vec(v₁.coords - v₂.coords)

LinearAlgebra.norm(p::Union{Vec,Point}) = norm(p.coords)
LinearAlgebra.norm(s::Union{DirSegment,Segment}) = norm(s.p₁ - s.p₂)

LinearAlgebra.cross(v₁::Vec{2}, v₂::Vec{2}) = v₁.x * v₂.y - v₁.y * v₂.x
LinearAlgebra.cross(v₁::Vec{3}, v₂::Vec{3}) = v₁.coords × v₂.coords


slope(s::Union{Segment,DirSegment}) = (s[1].y-s[2].y)/(s[1].x-s[2].x)


# O(1) operations

"""
    orientation(p₀::Point{2}, p₁::Point{2}, p₂::Point{2}; tol=0)

Return the orientation of turning from `p₁` to `p₂`, with `p₀` as pivot.

# Examples

```
"
       · p₁
      /
     /
    /
p₀ ·--------· p₂
"
```

The example bellow is illustrated by the ascii draw above.

```jldoctest
julia> p₀ = Point(0, 0)
Point{2, Int64}([0, 0])

julia> p₁ = Point(1, 2)
Point{2, Int64}([1, 2])

julia> p₂ = Point(2, 1)
Point{2, Int64}([2, 1])

julia> orientation(p₀, p₁, p₂)
ClockWise::Orientation = 1

julia> orientation(p₀, p₂, p₁)
CounterClockWise::Orientation = -1
```
"""
function orientation(p₀::Point{2}, p₁::Point{2}, p₂::Point{2}; tol=0)
    Orientation(Int(signapprox((p₂ - p₀) × (p₁ - p₀); tol)))
end

"""
    orientation(s::DirSegment{2}, p::Point{2}; [tol=0])

Return the orientation of the point `p` respect to the directed segment `s`.

# Examples

```
"
       ^
      /
   s /     · p
    /
   ·
"
```

The example bellow is illustrated by the ascii draw above.

```jldoctest
julia> s = DirSegment(Point(0, 0), Point(1, 2))
DirSegment{2, Int64}(Point{2, Int64}([0, 0]), Point{2, Int64}([1, 2]))

julia> p = Point(2, 1)
Point{2, Int64}([2, 1])

julia> orientation(s, p)
ClockWise::Orientation = 1
```
"""
orientation(s::DirSegment{2}, p::Point{2}; tol=0) = orientation(s[1], s[2], p; tol)

"""
    orientation(s₁::Union{Segment{2},DirSegment{2}}, s₂::Union{Segment{2},DirSegment{2}}; tol=0)

Return the orientation of `s₂` respect to the segment `s₁`, as if both segments
shared origin, that is `s₁[1] == s₂[1]`.

# Examples

```
"
       ^
      /
  s₁ /
    /
   ·--------> s₂
"
```

The example bellow is illustrated by the ascii draw above.

```jldoctest
julia> s₁ = DirSegment(Point(0, 0), Point(1, 2))
DirSegment{2, Int64}(Point{2, Int64}([0, 0]), Point{2, Int64}([1, 2]))

julia> s₂ = DirSegment(Point(0, 0), Point(2, 1))
DirSegment{2, Int64}(Point{2, Int64}([0, 0]), Point{2, Int64}([1, 2]))

julia> orientation(s₁, s₂)
ClockWise::Orientation = 1

julia> orientation(s₂, s₁)
CounterClockWise::Orientation = -1
```
"""
orientation(s₁::Union{Segment{2},DirSegment{2}}, s₂::Union{Segment{2},DirSegment{2}}; tol=0) = orientation(Point(0,0), Point(s₁[2] - s₁[1]), Point(s₂[2] - s₂[1]))

"""
    isintersect(s₁::Union{Segment{2},DirSegment{2}}, s₂::Union{Segment{2},DirSegment{2}}; tol=0)

Return wheter `s₁` intersects with `s₂` or not.

# Examples

```jldoctest
julia> s₁ = Segment(Point(0, 0), Point(2, 2))
Segment{2, Int64}(Point{2, Int64}([0, 0]), Point{2, Int64}([2, 2]))

julia> s₂ = Segment(Point(0, 2), Point(2, 0))
Segment{2, Int64}(Point{2, Int64}([0, 2]), Point{2, Int64}([2, 0]))

julia> isintersect(s₁, s₂)
true

julia> s₁ = Segment(Point(0, 0), Point(0, 1))
Segment{2, Int64}(Point{2, Int64}([0, 0]), Point{2, Int64}([0, 1]))

julia> s₂ = Segment(Point(1, 1), Point(1, 0))
Segment{2, Int64}(Point{2, Int64}([1, 1]), Point{2, Int64}([1, 0]))

julia> isintersect(s₁, s₂)
false
```
"""
function isintersect(s₁::DirSegment{2}, s₂::DirSegment{2}; tol=0)
    d₁ = orientation(s₂, s₁[1]; tol)
    d₂ = orientation(s₂, s₁[2]; tol)
    d₃ = orientation(s₁, s₂[1]; tol)
    d₄ = orientation(s₁, s₂[2]; tol)
    if Int(d₁)*Int(d₂) == -1 && Int(d₃)*Int(d₄) == -1
        return true
    elseif d₁ == Colinear && _onsegment(s₂, s₁[1])
        return true
    elseif d₂ == Colinear && _onsegment(s₂, s₁[2])
        return true
    elseif d₃ == Colinear && _onsegment(s₁, s₂[1])
        return true
    elseif d₄ == Colinear && _onsegment(s₁, s₂[2])
        return true
    end
    return false
end
isintersect(s₁::Union{Segment{2},DirSegment{2}}, s₂::Union{Segment{2},DirSegment{2}}; tol=0) =
    isintersect(DirSegment(s₁), DirSegment(s₂); tol)

# ONLY CALL IF YOU UNDERSTAND
#
# Determines wheter a point is inside a region delimited by two extrema points.
# It is helpful if you know `(s[1], p, s[2])` are colineal the three.
#
#       ·- - - - - -· s[2]
#       |       __/ |
#       |    __· p  |
#       | __/       |
#       |/          |
#  s[1] ·- - - - - -·
#
_onsegment(s::Union{Segment{2},DirSegment{2}}, p::Point{2}) =
    (min(s[1].x, s[2].x) ≤ p.x ≤ max(s[1].y, s[2].y)) &&
    (min(s[1].x, s[2].x) ≤ p.x ≤ max(s[1].y, s[2].y))


function isintersect(s₁::Union{Segment{2},DirSegment{2}}, r::Ray{2}; tol=0)
    p₁ = r.point
    xₘᵢₙ, xₘₐₓ, yₘᵢₙ, yₘₐₓ = extrema((p₁.x, s₁[1].x, s₁[2].x)), extrema((p₁.y, s₁[1].y, s₁[2].y))
    l = √((xₘₐₓ - xₘᵢₙ)^2 + (yₘₐₓ - yₘᵢₙ)^2)
    p₂ = r.point + l*r.vector
    s₂ = DirSegment(p₁, p₂)
    isintersect(s₁, s₂)
end

function isintersect(s₁::Union{Segment{2},DirSegment{2}}, r::Line{2}; tol=0)
    p₁, p₂ = r.point, r.point+r.vector
    s₂ = DirSegment(p₁, p₂)
    d₁ = orientation(s₂, s₁[1]; tol)
    d₂ = orientation(s₂, s₁[2]; tol)
    return (Int(d₁)*Int(d₂) == -1)
end
