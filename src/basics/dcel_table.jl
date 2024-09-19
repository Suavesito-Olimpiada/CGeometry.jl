module DCELs

import StaticArrays: MVector
import ..CGeometry: DirSegment, Polygon

export DCEL,
    Vertex, VertexHandle,
    Halfedge, HalfedgeHandle,
    Face, FaceHandle,
    nvertices, vertices, points,
    nsegments, segments,
    nedges, nhalfedges, halfedges,
    nfaces, faces,
    point, halfedge, halfedges, faces,
    segment, origin, vertex, face, prev, next, twin, facet,
    vertices, halfedges, segments,
    insertdiagonal!, updatefaces!, explode

# - `segments` is the list of half edges represented by a directed segment
# - `vertices` contains the indices for outgoing `halfedges` and the index** in `segments`
#   of the point they correspond to (the first in the `DirSegment`)
# - `halfedges` contains a vector of tuples which indices represent the following
#    [1] self: index on `segments` to the half-edge**
#    [2] vertex: index on `vertices` to the starting vertex
#    [3] face: index on `faces` to the incident face
#    [4] prev: index on `halfedges` to the previous half-edge
#    [5] next: index on `halfedges` to the next half-edge
#    [6] twin: index on `halfedges` to the twin half-edge !!! must be in the (2i - i&2)
# - `faces` contains the indices for `halfedges` representing half-edges of the incident
#   face, one half-edge for every independent facet. A last face is always kept as the
#   infinite (or exterior) face.
#
# ** NOTE: every negative index to `segments` means that the directed segment needs to be
# `reverse`d.
struct DCEL{T}
    segments::Vector{DirSegment{2,T}}
    vertices::Vector{Tuple{Int,Vector{Int}}}
    halfedges::NTuple{6,Vector{Int}}
    faces::Vector{Vector{Int}}
end

struct Vertex end
struct VertexHandle{T}
    dcel::DCEL{T}
    i::Int
end

struct Halfedge end
struct HalfedgeHandle{T}
    dcel::DCEL{T}
    i::Int
end

struct Face end
struct FaceHandle{T}
    dcel::DCEL{T}
    i::Int
end


# Printing

function _show(io, vertex::VertexHandle{T}, fmt=false) where {T}
    summary(io, vertex)
    fmt && println(io)
    p = point(vertex)
    fmt || print(io, "(")
    fmt && print(io, " ")
    print(io, "[", p.x, ",", p.y, "]")
    fmt || print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", vertex::VertexHandle{T}) where {T} = _show(io, vertex, true)
Base.show(io::IO, vertex::VertexHandle{T}) where {T} = _show(io, vertex)

function _show(io, halfedge::HalfedgeHandle{T}, fmt=false) where {T}
    summary(io, halfedge)
    b = origin(halfedge)
    e = origin(next(halfedge))
    fmt || print(io, "(")
    fmt && print(io, "\n ")
    print(io, "[", b.x, ",", b.y, "]")
    fmt || print(io, ",")
    fmt && print(io, "\n ")
    print(io, "[", e.x, ",", e.y, "]")
    fmt || print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", halfedge::HalfedgeHandle{T}) where {T} = _show(io, halfedge, true)
Base.show(io::IO, halfedge::HalfedgeHandle{T}) where {T} = _show(io, halfedge)

function _show(io, face::FaceHandle{T}, fmt=false) where {T}
    summary(io, face)
    fmt || print(io, "(")
    fmt && println(io)
    for facet in vertices(face)
        print(io, "[")
        fmt && println(io)
        for vertex in facet
            fmt && print(io, " ")
            p = point(vertex)
            print(io, "[", p.x, ",", p.y, "],")
            fmt && println(io)
        end
        print(io, "],")
    end
    fmt || print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", face::FaceHandle{T}) where {T} = _show(io, face, true)
Base.show(io::IO, face::FaceHandle{T}) where {T} = _show(io, face)


# basics

Base.eltype(::DCEL{T}) where T = T

Base.getindex(dcel::DCEL, i) = i < 0 ? reverse(dcel.segments[-i]) : dcel.segments[i]
Base.getindex(dcel::DCEL, ::Type{Vertex}, i) = VertexHandle(dcel, i)
Base.getindex(dcel::DCEL, ::Type{Halfedge}, i) = HalfedgeHandle(dcel, i)
Base.getindex(dcel::DCEL, ::Type{Face}, i) = FaceHandle(dcel, i)

Base.getindex(h::VertexHandle, j) = h.dcel.vertices[h.i][j]
Base.getindex(h::VertexHandle) = h.dcel.vertices[h.i]
Base.getindex(h::HalfedgeHandle, j) = h.dcel.halfedges[j][h.i]
Base.getindex(h::HalfedgeHandle) = ntuple(i -> h[i], 6)
Base.getindex(h::FaceHandle) = h.i == 0 ? last(h.dcel.faces) : h.dcel.faces[h.i]

Base.setindex!(h::HalfedgeHandle, v, j) = h.dcel.halfedges[j][h.i] = v

function addhalfedge!(dcel::DCEL, h::NTuple{6, Int})
    for i in 1:6
        push!(dcel.halfedges[i], h[i])
    end
end

# dcel
nvertices(dcel::DCEL) = length(dcel.vertices)
vertices(dcel::DCEL) = map(i -> dcel[Vertex, i], eachindex(dcel.vertices))
points(dcel::DCEL) = map(i -> point(dcel[Vertex, i]), eachindex(dcel.vertices))

nsegments(dcel::DCEL) = length(dcel.segments)
segments(dcel::DCEL) = copy(dcel.segments)

nedges(dcel::DCEL) = nhalfedges(dcel) >> 1
nhalfedges(dcel::DCEL) = length(dcel.halfedges[1])
halfedges(dcel::DCEL) = map(i -> dcel[Halfedge, i], eachindex(dcel.halfedges[1]))

nfaces(dcel::DCEL) = length(dcel.faces)
faces(dcel::DCEL) = map(i -> dcel[Face, i], eachindex(dcel.faces))

# vertex
point(h::VertexHandle) = h.dcel[h[1]][1]
halfedge(v::VertexHandle, f::FaceHandle) = v.dcel[Halfedge, findfirst(==(f.i), v[2])]
halfedges(h::VertexHandle) = map(i -> h.dcel[Halfedge, i], h[2])
faces(h::VertexHandle) = map(i -> face(h.dcel[Halfedge, i]), h[2])

# halfedge
segment(h::HalfedgeHandle) = h.dcel[h[1]]
origin(h::HalfedgeHandle) = segment(h)[1]
vertex(h::HalfedgeHandle) = h.dcel[Vertex, h[2]]
face(h::HalfedgeHandle) = h.dcel[Face, h[3]]
prev(h::HalfedgeHandle) = h.dcel[Halfedge, h[4]]
next(h::HalfedgeHandle) = h.dcel[Halfedge, h[5]]
twin(h::HalfedgeHandle) = h.dcel[Halfedge, h[6]]

struct FacetIterator{T}
    start::HalfedgeHandle{T}
end

function Base.length(b::FacetIterator)
    len = 1
    tmp = next(b.start)
    while tmp != b.start
        tmp = next(tmp)
        len += 1
    end
    return len
end
# Base.IteratorSize(::FacetIterator) = Base.SizeUnknown()
Base.iterate(b::FacetIterator) = (b.start, next(b.start))
Base.iterate(b::FacetIterator, s::HalfedgeHandle) =
    if b.start == s
        return nothing
    else
        return (s, next(s))
    end

function facet(::Type{Vertex}, h::HalfedgeHandle)
    start = h
    border = [vertex(h)]
    tmp = next(h)
    while start != tmp
        push!(border, vertex(tmp))
        tmp = next(tmp)
    end
    return border
end
function facet(::Type{Halfedge}, h::HalfedgeHandle)
    start = h
    border = [h]
    tmp = next(h)
    while start != tmp
        push!(border, tmp)
        tmp = next(tmp)
    end
    return border
end
facet(h::HalfedgeHandle) = facet(Halfedge, h)

# face
vertices(h::FaceHandle) = map(i -> facet(Vertex, h.dcel[Halfedge, i]), h[])
halfedges(h::FaceHandle) = map(i -> facet(Halfedge, h.dcel[Halfedge, i]), h[])
segments(h::FaceHandle) = map(i -> segment.(Iterators.facet(Halfedge, h.dcel[Halfedge, i])), h[])


# receives a polygon edges in order
function DCEL(segments::AbstractVector{DirSegment{2,T}}) where {T}
    n = length(segments)
    vertices = [(i, [2i - 1, 2mod(i - 1, 1:n)]) for i in 1:n]
    d(i) = (i+1)>>1
    halfedges = (
        [isodd(i) ? d(i) : -d(i) for i in 1:2n], # self
        [isodd(i) ? d(i) : mod(d(i)+1, 1:n) for i in 1:2n], # vertex
        [isodd(i) ? 1 : 0 for i in 1:2n], # face
        [isodd(i) ? 2mod(d(i)-1, 1:n) - 1 : 2mod(d(i)+1, 1:n) for i in 1:2n], # prev
        [isodd(i) ? 2mod(d(i)+1, 1:n) - 1 : 2mod(d(i)-1, 1:n) for i in 1:2n], # next
        [isodd(i) ? 2d(i) : 2d(i)-1 for i in 1:2n], #twin
    )
    faces = [[1], [2]] # outside face is last
    return DCEL(copy(segments), vertices, halfedges, faces)
end

DCEL(p::Polygon{T}) where {T} = DCEL(map(i -> p[DirSegment, i], 1:length(p)))


"""
Lazy operations, similar to `Base.Iterators`.
"""
module Iterators

using ..DCELs: DCEL, Vertex, VertexHandle, Halfedge, HalfedgeHandle, Face, FaceHandle, FacetIterator, face, vertex, segment, point

vertices(dcel::DCEL) = Base.Iterators.map(i -> dcel[Vertex, i], eachindex(dcel.vertices))
points(dcel::DCEL) = Base.Iterators.map(i -> point(dcel[Vertex, i]), eachindex(dcel.vertices))
segments(dcel::DCEL) = Base.Iterators.map(identity, dcel.segments)
halfedges(dcel::DCEL) = Base.Iterators.map(i -> dcel[Halfedge, i], eachindex(dcel.halfedges[1]))
faces(dcel::DCEL) = Base.Iterators.map(i -> dcel[Face, i], eachindex(dcel.faces))

halfedges(h::VertexHandle) = Base.Iterators.map(i -> h.dcel[Halfedge, i], h[2])
faces(h::VertexHandle) = Base.Iterators.map(i -> face(h.dcel[Halfedge, i]), h[2])

facet(::Type{Vertex}, h::HalfedgeHandle) = Base.Iterators.map(vertex, FacetIterator(h))
facet(::Type{Halfedge}, h::HalfedgeHandle) = FacetIterator(h)

vertices(h::FaceHandle) = Base.Iterators.map(i -> facet(Vertex, h.dcel[Halfedge, i]), h[])
halfedges(h::FaceHandle) = Base.Iterators.map(i -> facet(Halfedge, h.dcel[Halfedge, i]), h[])
segments(h::FaceHandle) = Base.Iterators.map(i -> Base.Iterators.map(segment, facet(Halfedge, h.dcel[Halfedge, i])), h[])

end

function insertdiagonal!(dcel::DCEL{T}, v1::VertexHandle, v2::VertexHandle, f::FaceHandle{T}) where {T}
    if f ∉ Iterators.faces(v1) || f ∉ Iterators.faces(v2)
        throw(ArgumentError("Vertex $(v1) or $(v2) does not belong to $(f)"))
    end
    e1 = filter(h -> face(h) == f, Iterators.halfedges(v1))
    e2 = filter(h -> face(h) == f, Iterators.halfedges(v2))
    insertdiagonal!(e1, e2)
end

insertdiagonal!(dcel::DCEL{T}, e1::Integer, e2::Integer) where {T} = insertdiagonal!(dcel[Halfedge, e1], dcel[Halfedge, e2])

function insertdiagonal!(e1::HalfedgeHandle{T}, e2::HalfedgeHandle{T}) where {T}
    if iszero(face(e1).i) || iszero(face(e2).i) || face(e1) != face(e2)
        println(e1, ": ", e1.i, " ", face(e1).i)
        println(e2, ": ", e2.i, " ", face(e2).i)
        throw(ArgumentError("Halfedges that are not from the same face cannot form a diagonal"))
    end
    dcel = e1.dcel
    ns = nsegments(dcel)
    nv = nvertices(dcel)
    nh = nhalfedges(dcel)
    pe1 = prev(e1)
    pe2 = prev(e2)
    # insert new segments, halfedges, and add the later to vertices
    # segment
    segment = DirSegment(origin(e1), origin(e2))
    push!(dcel.segments, segment)
    # vertices
    v1 = vertex(e1)
    v2 = vertex(e2)
    push!(v1[2], nh + 1)
    push!(v2[2], nh + 2)
    # halfedges
    addhalfedge!(dcel, (
        ns + 1,
        v1.i,
        face(e2).i,
        prev(e1).i,
        e2.i,
        nh + 2,
    ))
    addhalfedge!(dcel, (
        -(ns + 1),
        v2.i,
        face(e1).i,
        prev(e2).i,
        e1.i,
        nh + 1,
    ))
    # modify existing ones
    # next
    pe1[5] = nh + 1
    pe2[5] = nh + 2
    # prev
    e1[4] = nh + 2
    e2[4] = nh + 1
    # return
    dcel[Halfedge, nh+1]
end

function Base.in(p::Point{2,T}, face::FaceHandle{T}) where {T}
    cruses = 0
    for facet in face
        for hedge in facet
            seg = segment(hedge)
            seg = seg[1].y < seg[2].y ? seg : reverse(seg)
            if seg[1].y < p.y < seg[2].y
                cruses += clamp(Int(orientation(segment(hedge), p)), 0, 1)
            end
        end
    end
    isodd(cruses)
end

function updatefaces!(dcel::DCEL{T}) where {T}
    nf = nfaces(dcel)
    deleteat!(dcel.faces, 1:nf-1) # the infinite side is left
    tag = falses(nhalfedges(dcel))
    # first do the outside face
    n = dcel[Halfedge, dcel.faces[1][1]]
    while !tag[n.i]
        n[3] = 0
        tag[n.i] = true
        n = next(n)
    end
    # then the rest
    f = 1
    for h in Iterators.halfedges(dcel)
        if tag[h.i]
            continue
        end
        n = h
        while !tag[n.i]
            n[3] = f
            tag[n.i] = true
            n = next(n)
        end
        splice!(dcel.faces, f:f-1, [[h.i]])
        f += 1
    end
end

function explode(dcel::DCEL{T}) where {T}
    nfaces = length(Iterators.faces(dcel))
    polys = Vector{DCEL{T}}(undef, nfaces)
    head = nfaces
    for face in Iterators.faces(dcel)
        segments = collect(first(Iterators.segments(face)))
        if head != 0
            polys[nfaces-head+1] = DCEL(segments)
        end
        head -= 1
    end
    return polys
end

end
