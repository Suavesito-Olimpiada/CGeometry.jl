include("./plotting.jl")

using CGeometry, CGeometry.DCELs


function animmakemonotone(dcel::DCEL, t = 0.0; lims = nothing, heads = false, anim = nothing)
    plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
    xs, ys = if isnothing(lims)
        autolims(dcel.points)
    else
        lims.x, lims.y
    end
    xlims!(xs)
    ylims!(ys)
    plot!(dcel; heads)
    flag = true
    ip = 0
    np = nvertices(dcel)
    return makemonotone(
        dcel; cb = (i, dcel, tree, helper, v, pk) -> begin
            ip += 1
            @info "$ip of $np"
            if v isa Vertex
                p = point(v)
                plot!([xs[1], xs[2]], [p.y, p.y]; c = pk == -3 ? :red : :white, ls = :dot, lw = 0.5, la = 1) # linestyle, linewidth, linealpha
                scatter!([p.x], [p.y]; c = :red, ms = 2, msw = 0.0)
            else
                p = segment(v)
                if pk
                    plot!([p[1].x, p[2].x], [p[1].y, p[2].y]; lw = 1, arrow = heads ? (:closed, :head, 1, 1) : nothing)
                else
                    plot!([p[1].x, p[2].x], [p[1].y, p[2].y]; c = :yellow, lw = 1, arrow = heads ? (:closed, :head, 1, 1) : nothing)
                end
            end
            isnothing(anim) || frame(anim)
            isnothing(anim) && gui()
            isnothing(anim) && sleep(t)
        end
    )
end

function testtriangulate!(dcel::DCEL, t = 0.0; heads = false, anim = nothing)
    nf = nfaces(dcel)
    ip = 0
    for i in 1:nf
        triangulatemonotone!(
            dcel, i;
            cb = (h, i, s, d) -> begin
                if ip != i
                    ip = i
                    @info "$i of $nf"
                end
                plot(d; heads)
                x, y = map(x -> tail(x[1])[1], s), map(x -> tail(x[1])[2], s)
                plot!(x, y; c = :white, lw = 1, arrow = !heads ? nothing : (:closed, :head))
                scatter!(x, y; c = :yellow, ms = 2)
                x, y = map(x -> tail(x[1])[1], h), map(y -> tail(y[1])[2], h)
                plot!(x, y; c = :red, lw = 0.5)
                u, _ = h[i]
                uo = tail(u)
                scatter!([uo[1]], [uo[2]]; c = :white, ms = 2)
                isnothing(anim) || frame(anim)
                isnothing(anim) && gui()
                isnothing(anim) && sleep(t)
            end
        )
    end
    DCELs.update!(dcel)
    plot(dcel; heads)
    println("dcel3 eulers characteristic: ", eulerschar(dcel))
    vc = tricolorvertices(dcel)
    nc, minc = findmin(count(==(i), vc) for i in 1:3)
    icolor = (findall(==(i), vc) for i in 1:3)
    color = (map(j -> dcel[Vertex, j], ic) for ic in icolor)
    xycolor = [[[point(v).coords[i] for v in c] for i in 1:2] for c in color]
    colors = (:blue, :red, :green)
    for (i, c) in enumerate(colors)
        scatter!(xycolor[i][1], xycolor[i][2]; c, ms = 2)
    end
    isnothing(anim) || frame(anim)
    isnothing(anim) && gui()
    isnothing(anim) && sleep(t)
    scatter!(xycolor[minc][1], xycolor[minc][2]; c = colors[minc], ms = 4)
    println("The ", nc, " ", colors[minc], " points are the guards!")
    isnothing(anim) || for _ in 1:30
        frame(anim)
    end
    isnothing(anim) && gui()
    return dcel, vc, anim
end

function testartgallery(n, t = 0.0; heads = false, anim = nothing)
    pol, dcel = newradialpoly(n)
    plot(dcel; heads)
    dcel2 = animmakemonotone(dcel, t; heads, anim)
    DCELs.update!(dcel2)
    plot(dcel2; heads)
    println("dcel eulers characteristic: ", eulerschar(dcel))
    println("dcel2 eulers characteristic: ", eulerschar(dcel2))
    dcel3, vc, anim = testtriangulate!(deepcopy(dcel2), t; heads, anim)
    return dcel, dcel2, dcel3, vc, anim
end


# dcel, dcel2, dcel3, colors, ann = testartgallery(40; anim = Animation())
# gif(ann, "triangulation.gif"; fps = 30)
