using Plots, LinearAlgebra, StaticArrays
using CGeometry

include("./plotting.jl")

theme(:dark)

testjarvismarch(points, delay = 0.0; anim = nothing) = jarvismarch(
    points; cb = (c, ps, ip, in, it) -> begin
        t = ps[it]
        p = ps[ip]
        n = ps[in]
        _x = p -> p.x
        _y = p -> p.y
        Xc = _x.(c)
        Yc = _y.(c)
        plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
        xs, ys = autolims(points)
        xlims!(xs)
        ylims!(ys)
        # all the points
        scatter!(_x.(points), _y.(points); ms = 2)
        # plot!(_x.(@view(points[i:end])), _y.(@view(points[i:end])))
        scatter!(_x.(@view(ps[it:end])), _y.(@view(ps[it:end])); ms = 2)
        plot!(Xc, Yc; leg = :none)
        scatter!(Xc, Yc; leg = :none, ms = 2)
        X = _x.([t, p, n])
        Y = _y.([t, p, n])
        plot!(X, Y; leg = :none)
        scatter!(X, Y; leg = :none, ms = 2)
        isnothing(anim) || frame(anim)
        isnothing(anim) && gui()
        isnothing(anim) && sleep(delay)
    end
)

testquickhull(points, delay = 0.0; anim = nothing) = quickhull(
    points; cb = (f, l, r) -> begin
        _x = p -> p.x
        _y = p -> p.y
        X = _x.(f)
        Y = _y.(f)
        if length(f) == 2
            plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
            xs, ys = autolims(points)
            xlims!(xs)
            ylims!(ys)
            scatter!(_x.(points), _y.(points); ms = 2)
            scatter!(X, Y; c = :red, leg = :none, ms = 2)
            plot!(X, Y)
        else
            scatter!(X, Y; c = :blue, leg = :none, ms = 2)
            plot!(X, Y)
        end
        isnothing(anim) || frame(anim)
        isnothing(anim) && gui()
        isnothing(anim) && sleep(delay)
    end
)

testgrahamscan(points, delay = 0.0; anim = nothing) = grahamscan(
    points; cb = (c, p, i) -> begin
        t = p[i]
        _x = p -> p.x
        _y = p -> p.y
        Xc = _x.(c)
        Yc = _y.(c)
        plot(; aspect_ratio = :equal, showaxis = false, leg = :none)
        xs, ys = autolims(points)
        xlims!(xs)
        ylims!(ys)
        scatter!(_x.(points), _y.(points); ms = 2)
        plot!(_x.(p[i:end] ∪ [c[begin]]), _y.(p[i:end] ∪ [c[begin]]))
        scatter!(_x.(@view(p[begin:i])), _y.(@view(p[begin:i])); ms = 2)
        plot!(Xc, Yc; leg = :none)
        scatter!(Xc, Yc; leg = :none, ms = 2)
        X = _x.([c[end], t])
        Y = _y.([c[end], t])
        plot!(X, Y; leg = :none)
        scatter!(X, Y; leg = :none, ms = 2)
        isnothing(anim) || frame(anim)
        isnothing(anim) && gui()
        isnothing(anim) && sleep(delay)
    end
)

# pts = [Point(rand(2)) for _ in 1:2^5]
# pts = [Point(normalize!(randn(2))) for _ in 1:2^6]

# ann = Animation()
# poljm = testjarvismarch(pts; anim = ann)
# gif(ann, "jarvismarch.gif"; fps = 30,)

# ann = Animation()
# polqh = testquickhull(pts; anim = ann)
# gif(ann, "quickhull.gif"; fps = 30)

# ann = Animation()
# polgs = testgrahamscan(pts; anim = ann)
# plot(DCEL(polgs))
# X = (x->x.x).(pts)
# Y = (y->y.y).(pts)
# scatter!(X, Y; ms = 2)
# gif(ann, "grahamscan.gif"; fps = 30)
