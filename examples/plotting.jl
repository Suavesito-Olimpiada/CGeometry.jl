using Plots, PlotThemes, LinearAlgebra, StaticArrays, DataStructures
using CGeometry
using CGeometry.DCELs

theme(:dark)

function testradialpoly(n)
    points = [Point(rand(2)) for _ ∈ 1:n]
    pol = Polygon(radialsort!(points))
    dcel = DCEL(pol)
end

function plotpoly!(poly; size=(x=0:1,y=0:1))
    for i in 1:length(poly)
        ds = poly[Segment, i]
        x = [ds[1][1], ds[2][1]]
        y = [ds[1][2], ds[2][2]]
        plot!(x, y; legend=:none)
    end
end

function plotpoly(poly; size=(x=0:1,y=0:1))
    plot(; aspect_ratio=:equal, showaxis=false)
    xlims!(1.05 .* tuple(size.x...)); ylims!(1.05 .* tuple(size.y...))
    plotpoly!(poly; size)
    gui()
end

function plotdcel!(dcel; size=(x=0:1,y=0:1), heads=true)
    hs = halfedges(dcel)
    for i in 1:length(hs)
        hs[i][][1] < 0 && continue
        ds = segment(hs[i])
        x = [ds[1].x, ds[2].x]
        y = [ds[1].y, ds[2].y]
        plot!(x, y; leg=:none, lw=1, arrow= !heads ? nothing : (:closed, :head))
    end
    verts = point.(vertices(dcel))
    x = map(p -> p.x, verts)
    y = map(p -> p.y, verts)
    scatter!(x, y; legend=:none, ms=1, msw=0.5)
end

function plotdcel(dcel; size=(x=0:1,y=0:1), heads=false)
    plot(; aspect_ratio=:equal, showaxis=false)
    xlims!(1.05 .* tuple(size.x...)); ylims!(1.05 .* tuple(size.y...))
    plotdcel!(dcel; size, heads)
    gui()
end

function animdcel(dcel, t=0.1; size=(x=0:1,y=0:1), heads=false)
    plot(; aspect_ratio=:equal, showaxis=false, leg=:none)
    xlims!(1.05 .* tuple(size.x...)); ylims!(1.05 .* tuple(size.y...))
    hs = halfedges(dcel)
    s2a(s) = [s[1].x, s[2].x], [s[1].y, s[2].y]
    for h in hs
        h[][1] < 0 && continue
        x, y = s2a(segment(h))
        px, py = s2a(segment(prev(h)))
        nx, ny = s2a(segment(next(h)))
        plot!(px, py; lw=0.6, arrow= !heads ? nothing : (:closed, :head), c=:red)
        plot!(x, y; lw=0.6, arrow= !heads ? nothing : (:closed, :head), c=:green)
        plot!(nx, ny; lw=0.6, arrow= !heads ? nothing : (:closed, :head), c=:blue)
        gui()
        readline()
        plot!(px, py; lw=0.5, arrow= !heads ? nothing : (:closed, :head), c=:gray)
        plot!(x, y; lw=0.5, arrow= !heads ? nothing : (:closed, :head), c=:gray)
        plot!(nx, ny; lw=0.5, arrow= !heads ? nothing : (:closed, :head), c=:gray)
    end
end

function animmakemonotone(dcel, t=0.1; size=(x=0:1,y=0:1), heads=false)
    plotdcel(dcel; size, heads)
    buffer = 1 + 0.05 # 5% buffer
    x = buffer * size.x
    y = buffer * size.y
    xlims!((x[1], x[2])); ylims!((y[1], y[2]))
    flag = true
    makemonotone(dcel; cb=(dcel, tree, helper, v, pk) -> begin
        if v isa VertexHandle
            p = point(v)
            plot!([-size.x, size.x], [p.y, p.y]; c=pk == -3 ? :red : :white, ls=:dot, lw=0.5, la=1) # linestyle, linewidth, linealpha
            scatter!([p.x], [p.y]; c=:red, ms=2, msw=0.0)
        else
            p = segment(v)
            if pk
                plot!([p[1].x, p[2].x], [p[1].y, p[2].y]; lw=1, arrow=heads ? (:closed, :head, 1, 1) : nothing)
            else
                plot!([p[1].x, p[2].x], [p[1].y, p[2].y]; c=:yellow, lw=1, arrow= heads ? (:closed, :head, 1, 1) : nothing)
            end
        end
        gui()
        sleep(t)
    end)
end

