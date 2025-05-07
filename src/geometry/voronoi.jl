using Random: shuffle!
using Base: Iterators

using .DCELs


export Voronoi, incrementalvoronoi


struct Voronoi{T, V <: AbstractVector{Point{2, T}}}
    points::V
    bbox::@NamedTuple{x::Tuple{T, T}, y::Tuple{T, T}}
    geometry::DCEL{T}
    cells::Dict{Face{T}, Point{2, T}}
end

function findface(p, geometry::DCEL)
    for cell in DCELs.Iterators.faces(geometry)
        if isinside(p, cell)
            return cell
        end
    end
    return nothing
end

function _split_intersection!(h, v)
    face = h.face
    # c is the nearest side to site, and C the furthest
    c, C = if isexterior(face)
        DCELs.split!(h.twin, v).twin, h
    else
        h, DCELs.split!(h, v)
    end
    face.halfedge = C
    return c, C
end

function _close_behind!(hp, h, prev, nf)
    hp = if isexterior(prev)
        DCELs.connect_diagonal!(hp, prev.prev, h)  # connects inside
    else
        DCELs.connect_diagonal!(hp.twin, h.twin, prev).twin  # connects outside
    end
    hp.face = nf
    hp.twin.face = prev.face
    return hp
end

function gencell(voronoi::Voronoi, site::Point{2}; dvts = Set{Int}(), dedges = Set{Int}(), cb = nothing)
    geometry, cells = voronoi.geometry, voronoi.cells
    cvts = Set{Int}()
    oface = CGeometry.findface(site, geometry)
    center = cells[oface]
    ray = Ray(mediatrix(center, site))
    h = oface.halfedge
    while !isintersect(ray, segment(h))
        h = h.next
    end
    q = intersection(ray, segment(h))
    v = if isempty(dvts)
        DCELs.new!(geometry, q)
    else
        dv = geometry[Vertex, pop!(dvts)]
        dv.point = q
        dv
    end
    vts = [v]
    bedges = Dict(v => h)
    start = h
    prev = h
    isnothing(cb) || cb(:gen, 1, voronoi; vts, start, v, center, ray)
    while true
        fc = h.face
        center = if isexterior(fc)
            cells[h.twin.face]
        else
            cells[h.face]
        end
        mtx = if isexterior(fc)
            mediatrix(center, site)
        else
            mediatrix(center, site; turn = ClockWise)
        end
        isnothing(cb) || cb(:gen, 2, voronoi; vts, h, v, center, mtx)
        while true
            v = h.vertex
            h = h.prev
            if isexterior(fc)
                mtx = mediatrix(cells[h.twin.face], site)
            end
            if v.i <= 4 && isexterior(fc)
                push!(vts, v)
            end
            isnothing(cb) || cb(:gen, 2, voronoi; vts, start, h, v, center, mtx)
            if v.i > 4
                push!(cvts, v.i)
            end
            if isintersect(mtx, segment(h)) || h.twin == start
                break
            end
            push!(dedges, h.i - iseven(h.i))
        end
        fc.halfedge = h
        if h.twin == start
            break
        end
        h = h.twin
        s = segment(h)
        q = intersection(mtx, s)
        v = if isempty(dvts)
            DCELs.new!(geometry, q)
        else
            dv = geometry[Vertex, pop!(dvts)]
            dv.point = q
            dv
        end
        bedges[v] = h
        prev = h
        isnothing(cb) || cb(:gen, 3, voronoi; vts, start, h, v, center, mtx)
        push!(vts, v)
    end
    return vts, bedges, union(dvts, cvts), dedges
end

function add_site!(voronoi::Voronoi, site::Point{2}; dvts = Set{Int}(), dedges = Set{Int}(), cb = nothing)
    geometry, cells = voronoi.geometry, voronoi.cells
    vts, bedges, dvts, dedges = CGeometry.gencell(voronoi, site; dvts, dedges, cb)
    isnothing(cb) || cb(:add, 1, voronoi; vts)
    for (_, h) in bedges
        DCELs.unlink!(h)
    end
    for i in dedges
        DCELs.delete!(geometry[Halfedge, i])
    end
    for i in dvts
        DCELs.unlink!(geometry[Vertex, i])
    end
    push!(vts, first(vts))
    nf = DCELs.new!(geometry, Face)
    prev = bedges[vts[begin]]
    hp = prev
    isnothing(cb) || cb(:add, 2, voronoi; vts, prev)
    for (i, v) in enumerate(vts)
        if !isboundingbox(v)
            h = bedges[v]
            DCELs.connect_vertex!(h, v)
            isnothing(cb) || cb(:add, 4, voronoi; vts, prev, h, v)
        else
            isnothing(cb) || cb(:add, 3, voronoi; vts, prev, v)
        end
        if i == 1
            continue
        end
        hp = if isempty(dedges)
            DCELs.new!(geometry, Halfedge)
        else
            geometry[Halfedge, pop!(dedges)]
        end
        if isboundingbox(v)
            DCELs.connect_vertex!(hp.twin, v)
            DCELs.connect_on_vertex!(hp, prev.prev)
            hp.face = nf
            hp.twin.face = prev.face
            prev = hp.twin
            isnothing(cb) || cb(:add, 5, voronoi; vts, prev, v, hp)
        else
            if i == length(vts) && isexterior(h.twin)
                CGeometry._close_behind!(hp, h.twin.next, prev, nf)  # prev outside, h is start
                isnothing(cb) || cb(:add, 6, voronoi; vts, prev, v, h, hp)
                break
            end
            CGeometry._close_behind!(hp, h, prev, nf)
            prev = h
            isnothing(cb) || cb(:add, 3, voronoi; vts, prev, v, hp)
        end
    end
    nf.halfedge = hp
    for h in nf
        h.face = nf
    end
    isnothing(cb) || cb(:add, 7, voronoi; vts, nf, hp)
    isnothing(cb) || cb(:add, 8, voronoi)
    cells[nf] = site
    return voronoi
end

function incrementalvoronoi(P::AbstractVector{Point{2, T}}; cb = nothing) where {T}
    pts = unique!(copy(P))
    # shuffle the points in the voronoi diagram
    # shuffle!(pts)
    # set the borders
    x₋, x₊ = extrema(p -> p.x, pts)
    y₋, y₊ = extrema(p -> p.y, pts)
    Δx, Δy = (x₊ - x₋), (y₊ - y₋)
    x = (x₋ - 0.1Δx, x₊ + 0.1Δx)
    y = (y₋ - 0.1Δy, y₊ + 0.1Δy)
    # construct the bounding box, counter-clockwise
    geometry = DCEL(
        [
            Point(x[1], y[2]),  # -, +
            Point(x[1], y[1]),  # -, -
            Point(x[2], y[1]),  # +, -
            Point(x[2], y[2]),  # +, +
        ]
    )
    # the map face => point, as faces will only be added, and not removed
    # the first points is assigned to the box face
    cells = Dict(geometry[Face, 1] => first(pts))
    voronoi = Voronoi(pts, (; x, y), geometry, cells)
    dvts = Set{Int}()
    dedges = Set{Int}()
    for i in 2:length(pts)
        site = pts[i]
        add_site!(voronoi, site; dvts, dedges, cb = (args...; kws...) -> cb(i, args...; kws...))
    end
    return voronoi
end
