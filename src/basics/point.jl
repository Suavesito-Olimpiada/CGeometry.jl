export Point, PointF64

struct Point{N, T}
    coords::SVector{N, T}
    Point{N, T}(x) where {N, T} = new(x)
end
const PointF64 = Point{N, Float64} where {N}

Point(x...) = Point{length(x), Base.promote_typeof(x...)}(x)
Point(x::AbstractArray) = Point{length(x), Base.promote_typeof(x...)}(x)

(::Type{PointF64})(x...) = Point{length(x), Float64}(x)
(::Type{PointF64})(p::Point) = Point{length(p), Float64}(p.coords)

Base.zero(::Type{Point{N, T}}) where {N, T} = Point(ntuple(_ -> zero(T), N)...)
Base.zero(::Point{N, T}) where {N, T} = Point(ntuple(_ -> zero(T), N)...)

Base.eltype(::Point{N, T}) where {N, T} = T

Base.length(::Point{N}) where {N} = N

Base.getindex(x::Point, i::Integer) = x.coords[i]

function Base.getproperty(p::Point, sym::Symbol)
    return if sym == :x
        p.coords.x
    elseif sym == :y
        p.coords.y
    elseif sym == :z
        p.coords.z
    else
        getfield(p, sym)
    end
end

Base.isless(p::Point{2}, q::Point{2}) = (p.y, -p.x) < (q.y, -q.x)
