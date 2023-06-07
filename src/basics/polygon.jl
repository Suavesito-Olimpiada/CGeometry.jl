export Polygon

# always assumed to be counter-clockwise
struct Polygon{T}
    points::Vector{Point{2,T}}
end

Polygon(points::AbstractArray{Point{N,T}}) where {N,T} = Polygon{T}(points)

Base.eltype(::Polygon{T}) where {T} = T

Base.length(pol::Polygon) = length(pol.points)

Base.getindex(p::Polygon, i::Integer) = p.points[i]
function Base.getindex(p::Polygon, ::Type{Segment}, i::Integer)
    if i < length(p)
        Segment(p.points[i], p.points[i+1])
    else
        Segment(p.points[i], p.points[1])
    end
end
Base.getindex(p::Polygon, ::Type{DirSegment}, i::Integer) = DirSegment(p[Segment, i])

