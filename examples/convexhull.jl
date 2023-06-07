using Plots, LinearAlgebra, StaticArrays
using CGeometry

points = [Point(rand(2)) for _ in 1:10^2]

jarvismarch(points; cb=(c,p,n,t)->begin
    _x = p->p.x
    _y = p->p.y
    Xc = _x.(c)
    Yc = _y.(c)
    plot(Xc, Yc; leg=:none)
    scatter!(Xc, Yc; ms=2, leg=:none)
    X = _x.([p,n,t])
    Y = _y.([p,n,t])
    plot!(X, Y; leg=:none)
    scatter!(X, Y; ms=2, leg=:none)
    xlims!((-0.1, 1.1))
    ylims!(-0.1, 1.1)
    gui()
end)

quickhull(points; cb=(f,l,r)->begin
    X = (p->p.x).(f)
    Y = (p->p.y).(f)
    if length(f) == 2
        scatter(X, Y; ms=2, c=:red, leg=:none)
        xlims!((-0.1, 1.1))
        ylims!(-0.1, 1.1)
        plot!(X, Y)
        gui()
    else
        scatter!(X, Y; ms=2, c=:blue, leg=:none)
        plot!(X, Y)
        gui()
    end
end)

