using Plots, LinearAlgebra, StaticArrays
using CGeometry, CGeometry.DCELs
using Random: shuffle!

include("./plotting.jl")

theme(:dark)

function animvoronoi(pts, t = 0.0; check = false, lims = nothing, anim = nothing)
    n = length(pts)
    X = (x -> x.x).(pts)
    Y = (y -> y.y).(pts)
    ann = nothing
    ip = 0
    lims = if isnothing(lims)
        xs, ys = autolims(pts)
        (; x = xs, y = ys)
    else
        lims
    end
    voronoi = CGeometry.incrementalvoronoi(
        pts;
        cb = (i, stage, step, voronoi; kw...) -> begin
            if ip != i
                ip = i
                @info "$i of $n"
            end
            site = pts[i]
            geometry, cells = voronoi.geometry, voronoi.cells
            if stage == :add
                if step == 1
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, ann)
                elseif step == 2
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, prev = kw[:prev], ann)
                elseif step == 3
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, prev = kw[:prev], v = kw[:v], ann)
                elseif step == 4
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, prev = kw[:prev], h = kw[:h], v = kw[:v], ann)
                elseif step == 5
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, prev = kw[:prev], v = kw[:v], hp = kw[:hp], ann)
                elseif step == 6
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, prev = kw[:prev], h = kw[:h], v = kw[:v], hp = kw[:hp], ann)
                elseif step == 7
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, nf = kw[:nf], lw = 3, hp = kw[:hp], ann)
                elseif step == 8 && check
                    animdcel(geometry; heads = true)
                end
            elseif stage == :gen
                if step == 1
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, var"\n\ncenter" = kw[:center], start = kw[:start], ray = kw[:ray], v = kw[:v], ann)
                elseif step == 2
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, var"\n\ncenter" = kw[:center], h = kw[:h], mtx = kw[:mtx], v = kw[:v], ann)
                elseif step == 3
                    plotvor(geometry, pts, kw[:vts], i; lims, heads = true, var"\n\nsite" = site, var"\n\ncenter" = kw[:center], h = kw[:h], start = kw[:start], mtx = kw[:mtx], v = kw[:v], ann)
                end
            end
            # readline() == "q" && throw(ArgumentError("Damn"))
            isnothing(anim) || frame(anim)
            isnothing(anim) && gui()
            isnothing(anim) && sleep(t)
        end
    )
    vts = typeof(voronoi.geometry[Vertex, -1])[]
    plotvor(voronoi.geometry, pts, vts, length(pts); lims, ann = true)
    isnothing(anim) && gui()
    for _ in 1:30
        isnothing(anim) || frame(anim)
    end
    return voronoi, anim
end

# pts = [Point(rand(2)) for _ in 1:2^4]
# pts = [Point(randn(2)) for _ in 1:2^5]

# pts = [Point(normalize!(randn(2))) for _ in 1:2^4]
# pts = vec([Point(x, y) for x in 0.125:0.25:1.0, y in 0.125:0.25:1.0])
# pts = [Point(cos(θ) / 3 + 0.5, sin(θ) / 3 + 0.5) for θ in LinRange(0, 2π*(15/16), 16)]

# pts = [Point(normalize!(randn(2)) + 0.001rand(2)) for _ in 1:2^4]
# pts = vec([Point(x + 0.001rand(), y + 0.001rand()) for x in 0.125:0.25:1.0, y in 0.125:0.25:1.0])
# pts = [Point(cos(θ) / 3 + 0.5 + 0.001rand(), sin(θ) / 3 + 0.5 + 0.001rand()) for θ in LinRange(0, 2π*(15/16), 16)]

# shuffle!(pts)

# vor, ann = animvoronoi(pts; anim = Animation())
# gif(ann, "voronoi.gif"; fps = 30)
