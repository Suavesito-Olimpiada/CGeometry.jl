module DCELs

const Its = Base.Iterators

import StaticArrays: MVector
import ..CGeometry: DirSegment, Segment, Polygon, Point, point, BrokenInvariant

export DCEL,
    Vertex, Halfedge, Face,
    FaceIterator, VertexIterator,
    nvertices, vertices, points,
    nedges, nhalfedges, edges, halfedges,
    nfaces, faces, exterior, vertices,
    point, halfedge, halfedges,
    segment, tail, head, vertex, face, prev, next, twin,
    degree, isexterior, isinterior, isboundary,
    isboundingbox

# - `points` is the list of points of half edges
# - `vertices` contains the index for one outgoing `halfedge` of the point they correspond
#   to (the first in the `halfedge`, i.e. the `halfedge` pointing outward of the vertex)
# - `halfedges` contains a vector of tuples which indices represent the following
#    [1] vertex: index on `vertices` (and `points`) to the starting vertex
#    [2] face: index on `faces` to the incident face
#    [3] prev: index on `halfedges` to the previous half-edge
#    [4] next: index on `halfedges` to the next half-edge
#   *[5] twin: index on `halfedges` to the twin half-edge # it is not stored, it
#              is the 2n-1 for 2n and 2n for 2n-1.
# - `faces` contains the index for the `halfedge` representing the incident face. A last
#   face is always kept as the infinite (or exterior) face.
#
struct DCEL{T}
    points::Vector{Point{2, T}}
    vertices::Vector{Int}
    halfedges::NTuple{4, Vector{Int}}
    faces::Vector{Int}
end

struct Vertex{T}
    dcel::DCEL{T}
    i::Int
end

struct Halfedge{T}
    dcel::DCEL{T}
    i::Int
end

struct Face{T}
    dcel::DCEL{T}
    i::Int
end


#  ____       _       _   _
# |  _ \ _ __(_)_ __ | |_(_)_ __   __ _
# | |_) | '__| | '_ \| __| | '_ \ / _` |
# |  __/| |  | | | | | |_| | | | | (_| |
# |_|   |_|  |_|_| |_|\__|_|_| |_|\__, |
#                                 |___/
# (printing)

function _show(io, v::Vertex{T}, fmt = false) where {T}
    summary(io, v)
    if v.i == -1
        return nothing
    end
    fmt && println(io)
    fmt && print(io, " ")
    p = v.point
    print(io, "[", p.x, ",", p.y, "]")
    return nothing
end

Base.show(io::IO, ::MIME"text/plain", v::Vertex{T}) where {T} = _show(io, v, true)
Base.show(io::IO, v::Vertex{T}) where {T} = _show(io, v)

function _show(io, h::Halfedge{T}, fmt = false) where {T}
    summary(io, h)
    if h.i == -1
        return nothing
    end
    t = h.twin
    fmt || print(io, "(")
    fmt && print(io, "\n ")
    if h[1] > 0
        b = tail(h)
        print(io, "[", b.x, ",", b.y, "]")
    else
        print(io, "[,]")
    end
    fmt || print(io, ",")
    fmt && print(io, "\n ")
    if t[1] > 0
        e = tail(t)
        print(io, "[", e.x, ",", e.y, "]")
    else
        print(io, "[,]")
    end
    return fmt || print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", h::Halfedge{T}) where {T} = _show(io, h, true)
Base.show(io::IO, h::Halfedge{T}) where {T} = _show(io, h)

function _show(io, face::Face{T}, fmt = false) where {T}
    summary(io, face)
    if face.i == -1 || face.halfedge.i == -1
        return nothing
    end
    fmt || print(io, "(")
    fmt && println(io)
    for vertex in Iterators.vertices(face)
        fmt && print(io, " ")
        p = vertex.point
        print(io, "[", p.x, ",", p.y, "],")
        fmt && println(io)
    end
    return fmt || print(io, ")")
end

Base.show(io::IO, ::MIME"text/plain", face::Face{T}) where {T} = _show(io, face, true)
Base.show(io::IO, face::Face{T}) where {T} = _show(io, face)


#  ____       _        __   ____      _
# / ___|  ___| |_     / /  / ___| ___| |_
# \___ \ / _ \ __|   / /  | |  _ / _ \ __|
#  ___) |  __/ |_   / /   | |_| |  __/ |_
# |____/ \___|\__| /_/     \____|\___|\__|
# (set/get)

# basics

Base.eltype(::DCEL{T}) where {T} = T

Base.getindex(dcel::DCEL, i) = dcel.points[i]
Base.getindex(dcel::DCEL, ::Type{Vertex}, i) = Vertex(dcel, i)
Base.getindex(dcel::DCEL, ::Type{Halfedge}, i) = Halfedge(dcel, i)
Base.getindex(dcel::DCEL, ::Type{Face}, i) = Face(dcel, i)

Base.getindex(v::Vertex) = v.dcel.vertices[v.i]
Base.getindex(v::Vertex, ::Colon) = v.dcel[v.i]
Base.getproperty(v::Vertex, sym::Symbol) =
if sym == :halfedge
    v.dcel[Halfedge, v[]]
elseif sym == :point
    v[:]
else
    getfield(v, sym)
end
# set the field of the Halfedge {vertex,face,prev,next}
_twin(h::Halfedge) = h.i + ifelse(isodd(h.i), 1, -1)
Base.getindex(h::Halfedge, j) = j == 5 ? _twin(h) : h.dcel.halfedges[j][h.i]
Base.getproperty(h::Halfedge, sym::Symbol) =
if sym == :point
    h.dcel[h[1]]
elseif sym == :vertex
    h.dcel[Vertex, h[1]]
elseif sym == :face
    h.dcel[Face, h[2]]
elseif sym == :prev
    h.dcel[Halfedge, h[3]]
elseif sym == :next
    h.dcel[Halfedge, h[4]]
elseif sym == :twin
    h.dcel[Halfedge, _twin(h)]
else
    getfield(h, sym)
end
Base.getindex(f::Face) = f.i == 0 ? f.dcel.faces[end] : f.dcel.faces[f.i]
Base.getproperty(f::Face, sym::Symbol) =
if sym == :halfedge
    f.dcel[Halfedge, f[]]
else
    getfield(f, sym)
end

# set the index of one halfedge representing the vertex
Base.setindex!(v::Vertex, i) = v.dcel.vertices[v.i] = i
Base.setindex!(v::Vertex, p::Point, ::Colon) = v.dcel.points[v.i] = p
_setproperty!(v::Vertex, ::Val{:halfedge}, h::Halfedge) =
let i = h.i
    v[] = i
    h
end
_setproperty!(v::Vertex, ::Val{:point}, p::Point) = v[:] = p
Base.setproperty!(v::Vertex, sym::Symbol, x) =
if sym == :halfedge
    _setproperty!(v, Val(sym), x::Halfedge)::Halfedge
elseif sym == :point
    _setproperty!(v, Val(sym), x::Point)::Point
else
    setfield!(v, sym, x)
end
# set the field of the Halfedge {vertex,face,prev,next}
Base.setindex!(h::Halfedge, i, j) = h.dcel.halfedges[j][h.i] = i
_setproperty!(h::Halfedge, ::Val{:vertex}, v::Vertex) =
let i = v.i
    h[1] = i
    h.dcel[Vertex, i]
end
_setproperty!(h::Halfedge, ::Val{:face}, f::Face) =
let i = f.i
    h[2] = i
    h.dcel[Face, i]
end
_setproperty!(h::Halfedge, ::Val{:prev}, p::Halfedge) =
let i = p.i
    h[3] = i
    h.dcel[Halfedge, i]
end
_setproperty!(h::Halfedge, ::Val{:next}, n::Halfedge) =
let i = n.i
    h[4] = i
    h.dcel[Halfedge, i]
end
Base.setproperty!(h::Halfedge, sym::Symbol, x) =
if sym == :vertex
    _setproperty!(h, Val(sym), x::Vertex)::Vertex
elseif sym == :face
    _setproperty!(h, Val(sym), x::Face)::Face
elseif sym ∈ (:prev, :next)
    _setproperty!(h, Val(sym), x::Halfedge)::Halfedge
else
    setfield!(h, sym, x)
end
# set the Halfedge that identifies the face
Base.setindex!(f::Face, i) =
if f.i == 0
    f.dcel.faces[end] = i
else
    f.dcel.faces[f.i] = i
end
Base.setproperty!(f::Face, sym::Symbol, h::Halfedge) =
if sym == :halfedge
    f[] = h.i
    h
else
    setfield!(f, sym, h)
end


#     _    ____ ___
#    / \  |  _ \_ _|
#   / _ \ | |_) | |
#  / ___ \|  __/| |
# /_/   \_\_|  |___|
# (api)

# dcel

nvertices(dcel::DCEL) = count(!=(-1), dcel.vertices)
vertices(dcel::DCEL) = map(i -> dcel[Vertex, i], findall(!=(-1), dcel.vertices))
points(dcel::DCEL) = dcel.points[findall(!=(-1), dcel.vertices)]

nedges(dcel::DCEL) = nhalfedges(dcel) >> 1
nhalfedges(dcel::DCEL) = count(!=(-1), dcel.halfedges[1])
edges(dcel::DCEL) = map(i -> segment(dcel[Halfedge, i]), filter!(isodd, findall(i -> (i != -1), dcel.halfedges[1])))
halfedges(dcel::DCEL) = map(i -> dcel[Halfedge, i], findall(!=(-1), dcel.halfedges[1]))

nfaces(dcel::DCEL) = count(!=(-1), @view(dcel.faces[begin:(end - 1)]))
faces(dcel::DCEL) = map(i -> dcel[Face, i], findall(!=(-1), @view(dcel.faces[begin:(end - 1)])))
exterior(dcel::DCEL) = dcel[Face, 0]

# vertices

point(v::Vertex) = v.dcel[v.i]
halfedge(v::Vertex) = v.dcel[Halfedge, v[]]
degree(v::Vertex) = length(v)
function halfedge(v::Vertex, f::Face)
    h = o = v.halfedge
    while h.face != f
        h = h.twin.next
        if h == o
            throw(ArgumentError(lazy"vertex $(v.i) is not in face $(f.i)"))
        end
    end
    return h
end
function halfedges(v::Vertex)
    hs = [v.halfedge]
    h = hs[1].twin.next
    while h != hs[1]
        push!(hs, h)
        h = h.twin.next
    end
    return hs
end
function faces(v::Vertex)
    o = v.halfedge
    fs = [o.face]
    h = o.twin.next
    while h != o
        push!(fs, h.face)
        h = h.twin.next
    end
    return fs
end
isboundingbox(v::Vertex) = 1 <= v.i <= 4
isinterior(v::Vertex) = all(isinterior, v)
isboundary(v::Vertex) = any(isboundary, v)

# halfedge

segment(h::Halfedge) = DirSegment(h.dcel[h[1]], h.dcel[h.twin[1]])
tail(h::Halfedge) = h.dcel[h[1]]
head(h::Halfedge) = h.dcel[h.next[1]]
vertex(h::Halfedge) = h.dcel[Vertex, h[1]]
face(h::Halfedge) = h.dcel[Face, h[2]]
prev(h::Halfedge) = h.dcel[Halfedge, h[3]]
next(h::Halfedge) = h.dcel[Halfedge, h[4]]
twin(h::Halfedge) = h.dcel[Halfedge, h.i + ifelse(isodd(h.i), 1, -1)]
isinterior(h::Halfedge) = !isexterior(h.face) && !isexterior(h.twin.face)
isboundary(h::Halfedge) = isexterior(h.face) || isexterior(h.twin.face)
isexterior(h::Halfedge) = isexterior(h.face)

# faces

face(::Type{Vertex}, h::Halfedge) = vertex.(FaceIterator(h))
face(::Type{Halfedge}, h::Halfedge) = collect(FaceIterator(h))

points(f::Face) = tail.(f)
segments(f::Face) = segment.(f)
vertices(f::Face) = vertex.(f)
halfedge(f::Face) = f.dcel[halfedge, f[]]
halfedges(f::Face) = collect(f)
isexterior(f::Face) = iszero(f.i)


#  ___ _                 _
# |_ _| |_ ___ _ __ __ _| |_ ___  _ __ ___
#  | || __/ _ \ '__/ _` | __/ _ \| '__/ __|
#  | || ||  __/ | | (_| | || (_) | |  \__ \
# |___|\__\___|_|  \__,_|\__\___/|_|  |___/
# (iterators)

# vertices

struct VertexIterator{T}
    start::Halfedge{T}
end

VertexIterator(v::Vertex) = VertexIterator(v.halfedge)

Base.eltype(::VertexIterator{T}) where {T} = Halfedge{T}
function Base.length(v::VertexIterator)
    len = 0
    for _ in v
        len += 1
    end
    return len
end
# Base.IteratorSize(::VertexIterator) = Base.SizeUnknown()
Base.iterate(v::VertexIterator) = (v.start, v.start.twin.next)
Base.iterate(v::VertexIterator, s::Halfedge) =
if v.start == s
    return nothing
else
    return (s, s.twin.next)
end

Base.eltype(::Vertex{T}) where {T} = Halfedge{T}
Base.length(v::Vertex) = length(VertexIterator(v))
Base.iterate(v::Vertex) = iterate(VertexIterator(v))
Base.iterate(v::Vertex, s::Halfedge) = iterate(VertexIterator(v), s)

# faces

struct FaceIterator{T}
    start::Halfedge{T}
end
FaceIterator(f::Face) = FaceIterator(f.halfedge)

Base.eltype(::FaceIterator{T}) where {T} = Halfedge{T}
function Base.length(b::FaceIterator)
    len = 0
    for _ in b
        len += 1
    end
    return len
end
# Base.IteratorSize(::FaceIterator) = Base.SizeUnknown()
Base.iterate(b::FaceIterator) = (b.start, b.start.next)
Base.iterate(b::FaceIterator, s::Halfedge) =
if b.start == s
    return nothing
else
    return (s, s.next)
end

Base.eltype(::Face{T}) where {T} = Halfedge{T}
Base.length(f::Face) = length(FaceIterator(f))
Base.iterate(f::Face) = iterate(FaceIterator(f))
Base.iterate(f::Face, s::Halfedge) = iterate(FaceIterator(f), s)

# base iterators

"""
Lazy operations, similar to `Base.Iterators`.
"""
module Iterators

    using ..DCELs: DCEL, Vertex, Halfedge, Halfedge, Face, FaceIterator, vertex, segment, point, nvertices, nhalfedges, nfaces

    const Its = Base.Iterators

    vertices(dcel::DCEL) = Its.map(i -> dcel[Vertex, i], findall(!=(-1), dcel.vertices))
    points(dcel::DCEL) = Its.map(i -> dcel.points[i], findall(!=(-1), dcel.vertices))
    edges(dcel::DCEL) = Its.map(i -> segment(dcel[Halfedge, i]), filter!(isodd, findall(!=(-1), dcel.halfedge[1])))
    halfedges(dcel::DCEL) = Its.map(i -> dcel[Halfedge, i], findall(!=(-1), dcel.halfedges[1]))
    faces(dcel::DCEL) = Its.map(i -> dcel[Face, i], findall(!=(-1), @view(dcel.faces[begin:(end - 1)])))

    halfedges(v::Vertex) = Its.map(identity, v)
    faces(v::Vertex) = Its.map(face, v)

    face(::Type{Vertex}, h::Halfedge) = Its.map(vertex, FaceIterator(h))
    face(::Type{Halfedge}, h::Halfedge) = FaceIterator(h)

    points(h::Face) = Its.map(tail, h)
    segments(h::Face) = Its.map(segment, h)
    vertices(h::Face) = Its.map(vertex, h)
    halfedges(h::Face) = Its.map(identity, h)

end

# construction

# receives a polygon points in counter-clockwise order
function DCEL(points::AbstractVector{Point{2, T}}) where {T}
    n = length(points)
    # references to the halfedge that have the vertex as its tail
    # as the twin will be next to the halfedge, we chose the odd ones
    # as the interior face, being the even ones the twins
    vertices = collect(1:2:(2n - 1))
    # tranforms [2n-1, 2n] -> n
    d(i) = i >> 1 + 1
    halfedges = (
        [isodd(i) ? d(i) : mod1(d(i), n) for i in 1:2n], # vertex
        [isodd(i) ? 1 : 0 for i in 1:2n], # face
        [isodd(i) ? mod(i - 2, 1:2n) : mod(i + 2, 1:2n) for i in 1:2n], # prev
        [isodd(i) ? mod(i + 2, 1:2n) : mod(i - 2, 1:2n) for i in 1:2n], # next
    )
    faces = [1, 2] # outside face is last
    return DCEL(copy(points), vertices, halfedges, faces)
end
DCEL(segments::AbstractVector{Union{Segment{2, T}, DirSegment{2, T}}}) where {T} =
    DCEL(map(x -> x[1], segments))
DCEL(p::Polygon{T}) where {T} = DCEL(p[:])

# geometric operations

# point/vertex: add, delete, modify, connect 2
# edge: create, delete, modify
# face: update

function new!(dcel::DCEL{T}, ::Type{Vertex}) where {T}
    k = length(dcel.vertices)
    z = zero(T)
    push!(dcel.points, Point(z, z))
    push!(dcel.vertices, -1)  # means no halfedge assigned
    # return
    return dcel[Vertex, k + 1]
end

function new!(dcel::DCEL, p::Point)
    v = new!(dcel, Vertex)
    v.point = p
    return v
end

function unlink!(v::Vertex)
    dcel = v.dcel
    if v.i == -1 || v.halfedge.i == -1
        return nothing
    end
    hv = v.halfedge
    # de-reference the halfedge
    v.halfedge = dcel[Halfedge, -1]
    for h in VertexIterator(hv)
        unlink!(h)
    end
    return nothing
end

function new!(dcel::DCEL, ::Type{Halfedge})
    k = length(dcel.halfedges[1])
    # halfedge
    push!(dcel.halfedges[1], -1)   # vertex
    push!(dcel.halfedges[2], -1)   # face
    push!(dcel.halfedges[3], k + 2)  # prev
    push!(dcel.halfedges[4], k + 2)  # next
    # twin
    push!(dcel.halfedges[1], -1)
    push!(dcel.halfedges[2], -1)
    push!(dcel.halfedges[3], k + 1)
    push!(dcel.halfedges[4], k + 1)
    # return
    return dcel[Halfedge, k + 1]
end

# unlink on the vertex origin of h
# p->h / t->nt  =>  p->nt, h<->t
# pt->t / h->n  =>  pt->n, h<->t
function unlink!(h::Halfedge)
    p = h.prev
    t = h.twin
    nt = t.next
    p.next = nt
    nt.prev = p
    h.prev = t
    t.next = h
    return nothing
end

function delete!(h::Halfedge)
    dcel = h.dcel
    # possible cases
    #
    #  \ p     n_/    |\ p     n_/            n_/    |\ p     n_/             n_/
    #   \|  h̸ \ /    nt \|  h̸ \ /       nt=h̸ \ /    nt \|  h̸ \ /_pt     nt=h̸ \ /_pt  nt=h̸=pt \
    #  v ·-----· w     v ·-----· w   v ·------· w     v ·-----· w     v ·-----· w     v ·-----· w
    #   / \ t̸  |\         \ t̸  |\       \ t̸=p |\         \ t̸             \ t̸=p           \ t̸=p=n
    #  /_nt    pt\             pt\           pt \
    #
    # de-reference the vertices
    t = h.twin
    v, w = h.vertex, t.vertex
    if v.i != -1
        # if h is the only halfedge of v, unlink it
        tp = h.prev.twin
        if tp == h
            v.halfedge = dcel[Halfedge, -1]
        else
            v.halfedge = tp
        end
    end
    if w.i != -1
        # if t is the only halfedge of w, unlink it
        tpt = t.prev.twin
        if tpt == t
            w.halfedge = dcel[Halfedge, -1]
        else
            v.halfedge = tpt
        end
    end
    h.vertex = dcel[Vertex, -1]
    t.vertex = dcel[Vertex, -1]
    # unlink v
    unlink!(h)
    # unlink w
    unlink!(t)
    # # delete faces
    h.face = dcel[Face, -1]
    t.face = dcel[Face, -1]
    return nothing
end

#  \         n_/
#   \  +t \   /
#    ·-----· · v
#   / \ +h   |\
#  /         p \
function connect_to!(h::Halfedge, p::Halfedge, n::Halfedge)
    dcel = h.dcel
    if dcel != n.dcel || dcel != p.dcel
        throw(ArgumentError("Not all arguments are from the same DCEL"))
    end
    t = h.twin
    v = n.vertex
    h.vertex = v
    h.prev = p
    p.next = h
    t.next = n
    n.prev = t
    return h
end

#  \ p         n_/
#   \|   +h \   /
#  v · ·-----· · w
#   /   \ +t   |\
#  /_nt        pt\
function connect!(h::Halfedge, p::Halfedge, nt::Halfedge, pt::Halfedge, n::Halfedge)
    t = h.twin
    connect_to!(h, p, nt)
    connect_to!(t, pt, n)
    return h
end
connect(p::Halfedge, nt::Halfedge, pt::Halfedge, n::Halfedge) =
    connect!(new!(p.dcel, Halfedge), p, nt, pt, n)

#  \ p         n_/
#   \|   +h \   /
#  v · ·-----· · w
#   /   \ +t   |\
#  /_nt        pt\
connect_diagonal!(h::Halfedge, p::Halfedge, n::Halfedge) =
    connect!(h, p, p.next, n.prev, n)
connect_diagonal(p::Halfedge, n::Halfedge) =
    connect_diagonal!(new!(p.dcel, Halfedge), p, n)

#  p                n
#   \      +h \      \
# ---· v ·-----· w ·---
#  \      \ +t      \
#   nt               pt
connect_edges!(h::Halfedge, p::Halfedge, n::Halfedge) =
    connect!(h, p, p.twin, n.twin, n)
connect_edges(p::Halfedge, n::Halfedge) =
    connect_edges!(new!(p.dcel, Halfedge), p, n)

#          n_/
#    +t \   /
#  ·-----· · v
#   \ +h   |\
#          p \
connect_on_vertex!(h::Halfedge, p::Halfedge) =
    connect_to!(h, p, p.next)
connect_on_vertex(p::Halfedge) =
    connect_on_vertex!(new!(p.dcel, Halfedge), p)

#
#    +t \
#  ·-----· · v
#   \ +h
#
function connect_vertex!(h::Halfedge, v::Vertex)
    dcel = h.dcel
    if dcel != v.dcel
        throw(ArgumentError("Not all arguments are from the same DCEL"))
    end
    h.vertex = v
    v.halfedge = h
    return h
end
connect_vertex(v::Vertex) =
    connect_vertex!(new!(v.dcel, Halfedge), v)

delete!(f::Face) = f.halfedge = f.dcel[Halfedge, -1]

function update!(dcel::DCEL)
    nf = nfaces(dcel)
    # delete the internal ones
    deleteat!(dcel.faces, 1:nf)
    tag = falses(length(dcel.halfedges[1]))
    # first we mark the outside face
    for h in FaceIterator(exterior(dcel))
        tag[h.i] = true
    end
    # then the inside ones
    f = 1
    for h in Iterators.halfedges(dcel)
        if tag[h.i]
            continue
        end
        for n in FaceIterator(h)
            tag[n.i] = true
            n[2] = f
        end
        insert!(dcel.faces, f, h.i)
        f += 1
    end
    for h in FaceIterator(dcel[Face, f])
        h[2] = 0
    end
    return
end

function split!(h::Halfedge, v::Vertex)
    dcel = h.dcel
    # relatives
    t = h.twin
    w = t.vertex
    n = h.next
    # add new halfedges
    hp = new!(dcel, Halfedge)
    ht = hp.twin
    # unlink and set new vertex
    unlink!(t)
    connect_vertex!(hp, v)
    connect_vertex!(t, v)
    connect_vertex!(ht, w)
    if n != t
        connect_diagonal!(hp, h, n)
    else
        connect_on_vertex!(hp, h)
    end
    # set the faces
    hp.face = h.face
    hp.twin.face = t.face
    # return
    return hp
end
split!(h::Halfedge, p::Point{2}) = split!(h, new!(h.dcel, p))
function split!(h::Halfedge{T}, x::T) where {T}
    s = segment(h)
    # calculate the point in the segment
    p = s[1] + x * (s[2] - s[1])
    return split!(h, p)
end

function splice!(g1::Halfedge, g2::Halfedge, h1::Halfedge, h2::Halfedge)
    u1, u2 = g1.twin, g2.twin
    v1, v2 = h1.twin, h2.twin
    # |\ g1    /_v2
    # u1\| h2_/
    #    \   /
    #      ∘
    # h1_/   \ g2
    #   /_v1 |\|
    #  /    u2 \
    #
    # set g1 and h2
    g1.next = h2
    h2.prev = g1
    # set v2 and g2
    v2.next = g2
    g2.prev = v2
    # set u2 and v1
    u2.next = v1
    v1.prev = u2
    # set h1 and u1
    h1.next = u1
    u1.prev = h1
    return nothing
end

function new!(dcel::DCEL, ::Type{Face})
    push!(dcel.faces, dcel.faces[end])
    dcel.faces[end - 1] = -1
    return dcel[Face, length(dcel.faces) - 1]
end

function explode(dcel::DCEL)
    DCELs.update!(dcel)
    nf = length(dcel.faces)
    polys = Vector{typeof(dcel)}(undef, nf)
    head = nf
    for f in 1:nf
        face = dcel[Face, f]
        ps = points(face)
        if head != 0
            polys[nf - head + 1] = DCEL(ps)
        end
        head -= 1
    end
    return polys
end

end
