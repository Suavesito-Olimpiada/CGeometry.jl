include("./plotting.jl")

using CGeometry
using CGeometry.DCELs

function testartgallery(n, t=0.0; heads=false)
    dcel = testradialpoly(n)
    plotdcel(dcel; heads)
    dcel2 = animmakemonotone(dcel, t; heads)
    updatefaces!(dcel2)
    plotdcel(dcel2; heads)
    dcel3 = deepcopy(dcel2)
    nf = nfaces(dcel3)
    for i in 1:nf-1
        triangulatemonotone!(dcel3, i;
            cb=(h, i, s, d) -> begin
                plotdcel(d; heads)
                x, y = map(x -> origin(x[1])[1], s), map(x -> origin(x[1])[2], s)
                plot!(x, y; c=:white, lw=1, arrow= !heads ? nothing : (:closed, :head))
                scatter!(x, y; c=:yellow)
                x, y = map(x -> origin(x[1])[1], h), map(y -> origin(y[1])[2], h)
                plot!(x, y; c=:red, lw=0.5)
                u, _ = h[i]
                uo = origin(u)
                scatter!([uo[1]], [uo[2]]; c=:white)
                gui()
                sleep(t)
            end
        )
    end
    updatefaces!(dcel3)
    plotdcel(dcel3; heads)
    println("dcel eulers characteristic: ", eulerschar(dcel))
    println("dcel2 eulers characteristic: ", eulerschar(dcel2))
    println("dcel3 eulers characteristic: ", eulerschar(dcel3))
    vc = tricolorvertices(dcel3)
    nc, minc = findmin(count(==(i), vc) for i = 1:3)
    icolor = (findall(==(i), vc) for i = 1:3)
    color = (map(j -> dcel3[Vertex, j], ic) for ic in icolor)
    xycolor = [[[point(v).coords[i] for v in c] for i in 1:2] for c in color]
    colors = (:blue, :red, :green)
    for (i,c) in enumerate(colors)
        scatter!(xycolor[i][1], xycolor[i][2]; c, ms=2)
    end
    gui()
    sleep(t)
    scatter!(xycolor[minc][1], xycolor[minc][2]; c=colors[minc], ms=3)
    println("The ", nc, " ", colors[minc], " points are the guards!")
    gui()
    dcel, dcel2, dcel3, vc
end
