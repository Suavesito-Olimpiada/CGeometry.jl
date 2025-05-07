using Plots, PlotThemes, LinearAlgebra, StaticArrays
using CGeometry
using CGeometry.DCELs


theme(:dark)

function autolims(points)
    x₋, x₊ = extrema(p -> p.x, points)
    y₋, y₊ = extrema(p -> p.y, points)
    Δx, Δy = (x₊ - x₋), (y₊ - y₋)
    x₋, x₊ = x₋ - 0.2Δx, x₊ + 0.2Δx
    y₋, y₊ = y₋ - 0.2Δy, y₊ + 0.2Δy
    return (x₋, x₊), (y₋, y₊)
end

function newradialpoly(n)
    points = [Point(rand(2)) for _ in 1:n]
    pivot = mean(points)
    poly = Polygon(radialsort!(points; pivot))
    dcel = DCEL(poly)
    return poly, dcel
end

function Plots.plot!(p::Point{2}; ann = nothing)
    scatter!(@SArray[p.x], @SArray[p.y]; legend = :none, ms = 2, msw = 0.5)
    if !isnothing(ann)
        annotate!([p.x + 0.025], [p.y + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot!(l::Line{2}; ann = nothing)
    xs = xlims()
    o = l.point
    p, q = point(l, xs[1], 1), point(l, xs[2], 1)
    x = @SArray[p.x, q.x]
    y = @SArray[p.y, q.y]
    plot!(x, y; leg = :none, lw = 1)
    if !isnothing(ann)
        annotate!([o.x + 0.025], [o.y + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot!(r::Ray{2}; ann = nothing)
    xs, ys = xlims(), ylims()
    o = r.point
    v = r.vector
    c, dim = if sign(v.x) == 1 && xs[1] != xs[2]
        (max(o.x, xs[1]), xs[2]), 1
    elseif sign(v.x) == -1 && xs[1] != xs[2]
        (xs[1], min(o.x, xs[2])), 1
    else
        if sign(v.y) == 1 && ys[1] != ys[2]
            (max(o.y, ys[1]), ys[2]), 2
        elseif sign(v.y) == -1 && ys[1] != ys[2]
            (ys[1], min(o.y, ys[2])), 2
        else
            (ys[1], ys[2]), 0
        end
    end
    p, q = if !iszero(dim)
        (
            point(r, c[1], dim),
            point(r, c[2], dim),
        )
    else
        (
            Point(xs[1], ys[1]),
            Point(xs[2], ys[2]),
        )
    end
    # p, q = point(r), r.vector
    x = @SArray[p.x, q.x]
    y = @SArray[p.y, q.y]
    plot!(x, y; leg = :none, lw = 1)
    plot!(o)
    if !isnothing(ann)
        annotate!([o.x + 0.025], [o.y + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot!(v::Vertex; ann = nothing, ms = 3)
    if v.i == -1
        return current()
    end
    p = v.point
    scatter!(@SArray[p.x], @SArray[p.y]; legend = :none, ms = ms, msw = 0.5)
    if !isnothing(ann)
        annotate!([p.x + 0.025], [p.y + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot!(h::Halfedge; heads = false, ann = nothing, lw = 1)
    if h.i == -1 || h.vertex.i == -1
        return current()
    end
    s = segment(h)
    x = @SArray[s[1].x, s[2].x]
    y = @SArray[s[1].y, s[2].y]
    plot!(x, y; leg = :none, lw = lw, arrow = heads ? (:closed, :head) : nothing)
    if !isnothing(ann)
        istwin = isodd(h.i) ? -1 : 1
        annotate!(@SArray[mean(x) + istwin * 0.025], @SArray[mean(y) + istwin * 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot!(f::Face{T}; heads = false, ann = nothing, lw = 1) where {T}
    if f.halfedge.i == -1
        return nothing
    end
    P = points(f)
    x = (p -> p.x).(P)
    y = (p -> p.y).(P)
    annx = mean(x)
    anny = mean(y)
    push!(x, x[begin])
    push!(y, y[begin])
    plot!(x, y; leg = :none, lw, arrow = heads ? (:closed, :head) : nothing)
    if !isnothing(ann)
        annotate!(@SArray[annx + 0.025], @SArray[anny + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot(f::Face; lims = nothing, heads = false, ann = nothing, lw = 1)
    plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
    pts = points(f)
    xs, ys = if isnothing(lims)
        autolims(pts)
    else
        lims.x, lims.y
    end
    xlims!(xs)
    ylims!(ys)
    plot!(f; heads, ann, lw)
    return nothing
end

function Plots.plot!(poly::Polygon; ann = nothing)
    n = 1
    _x = @SArray zeros(eltype(f), 2)
    _y = @SArray zeros(eltype(f), 2)
    for i in eachindex(poly)
        n += 1
        ds = poly[Segment, i]
        _x += x = @SArray[ds[1][1], ds[2][1]]
        _y += y = @SArray[ds[1][2], ds[2][2]]
        plot!(x, y; legend = :none)
    end
    if !isnothing(ann)
        annotate!(@SArray[mean(_x) / n + 0.025], @SArray[mean(_y) / n + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot(poly::Polygon; lims = nothing, ann = nothing)
    plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
    xs, ys = if isnothing(lims)
        autolims(poly.points)
    else
        lims.x, lims.y
    end
    xlims!(xs)
    ylims!(ys)
    plot!(poly; ann)
    return nothing
end

function Plots.plot!(dcel::DCEL; heads = false, ann = nothing, color = :auto)
    for e in edges(dcel)
        x = @SArray[e[1].x, e[2].x]
        y = @SArray[e[1].y, e[2].y]
        plot!(x, y; leg = :none, lw = 1, arrow = heads ? (:closed, :head) : nothing, color)
    end
    x = map(p -> p.point.x, DCELs.Iterators.vertices(dcel))
    y = map(p -> p.point.y, DCELs.Iterators.vertices(dcel))
    scatter!(x, y; legend = :none, ms = 2, msw = 0.5, color)
    if !isnothing(ann)
        annotate!(@SArray[maximum(x) + 0.025], @SArray[maximum(y) + 0.025], string(ann); annotationfontsize = 8)
    end
    return nothing
end

function Plots.plot(dcel::DCEL; lims = nothing, heads = false, ann = nothing)
    plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
    xs, ys = if isnothing(lims)
        autolims(dcel.points)
    else
        lims.x, lims.y
    end
    xlims!(xs)
    ylims!(ys)
    plot!(dcel; heads, ann)
    return nothing
end

function animdcel(dcel::DCEL, t = 0.0; lims = nothing, heads = false, ann = nothing)
    xs, ys = if isnothing(lims)
        autolims(dcel.points)
    else
        lims.x, lims.y
    end
    s2a(p1, p2) = @SArray[p1.x, p2.x], @SArray[p1.y, p2.y]
    for i in dcel.faces
        if i < 1
            break
        end
        h = dcel[Halfedge, i]
        plot(dcel; lims)
        xlims!(xs)
        ylims!(ys)
        for h in CGeometry.FaceIterator(h)
            plot!(h; heads, ann = isnothing(ann) ? nothing : "h ($(h.i))")
            gui()
        end
        t > 0 && sleep(t)
    end
    return nothing
end

function plotvor(dcel, pts, vts, m; heads = false, lw = 2, ann = nothing, lims = nothing, kw...)
    plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
    xs, ys = if isnothing(lims)
        autolims(pts)
    else
        lims.x, lims.y
    end
    xlims!(xs)  # xlims!((-0.2, 1.2))
    ylims!(ys)  # ylims!((-0.2, 1.2))
    plot!(dcel, color = :gray)
    scatter!((x -> x.x).(pts), (y -> y.y).(pts); ms = 2, msw = 0)
    scatter!((x -> x.x).(@view(pts[1:m])), (y -> y.y).(@view(pts[1:m])); ms = 2, msw = 0)
    scatter!((x -> x.point.x).(vts), (y -> y.point.y).(vts); ms = 3)
    isnothing(ann) || annotate!((x -> x.x + 0.025).(@view(pts[1:m])), (y -> y.y + 0.025).(@view(pts[1:m])), string.(1:m); annotationfontsize = 8)
    kw = values(kw)
    for (sym, val) in zip(keys(kw), kw)
        if val isa Halfedge
            plot!(val; lw, heads, ann = isnothing(ann) ? nothing : string(sym) * " ($(val.i))")
        elseif val isa Vertex
            plot!(val; ann = isnothing(ann) ? nothing : string(sym) * " ($(val.i))")
        elseif val isa Face
            plot!(val; lw, ann = isnothing(ann) ? nothing : string(sym) * " ($(val.i))")
        else
            plot!(val; ann = isnothing(ann) ? nothing : string(sym))
        end
    end
    return nothing
end
