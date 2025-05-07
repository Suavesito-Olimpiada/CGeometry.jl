export radialsort, radialsort!, eulerschar, mediatrix

# sorting points

using Statistics #, DataStructures

Statistics.mean(points::AbstractArray{<:Point{2}}) = Point(Statistics.mean(p -> p.coords, points))
Statistics.mean(points::AbstractArray{<:Vec{2}}) = Vec(Statistics.mean(p -> p.coords, points))

# sort counter-clockwise from x+, with origin in pivot
function radialsort!(points::AbstractArray{<:Point{2, T}}; pivot = nothing, direction = 0.0, turn = CounterClockWise, revnorm = false, rev = false) where {T}
    pt = if isnothing(pivot)
        zero(eltype(points))
    else
        pivot
    end
    function _less(p, q)
        o = Point(rotate(Vec(one(T), zero(T)), direction))
        z = zero(p)
        p = Point(p - pt)
        q = Point(q - pt)
        po = orientation(z, o, p)
        qo = orientation(z, o, q)
        if po == CounterClockWise && qo == ClockWise
            return turn == CounterClockWise
        elseif po == ClockWise && qo == CounterClockWise
            return turn == ClockWise
        elseif po == Colinear && qo == Colinear
            if sign(p.x * o.x) == -1 && sign(q.x * o.x) == 1
                return false
            elseif sign(p.x * o.x) == 1 && sign(q.x * o.x) == -1
                return true
            end
        end
        ori = orientation(z, p, q)
        if ori == Colinear
            np = norm(p)
            nq = norm(q)
            # sort by norm, nearest are 'less'
            return revnorm ? np > nq : np < nq
        end
        return ori == turn
    end
    return sort!(points; lt = _less, rev)
end
radialsort(points; pivot = nothing, direction = 0.0, turn = CounterClockWise, revnorm = false, rev = false) = radialsort!(copy(points); pivot, direction, turn, revnorm, rev)

eulerschar(dcel) = nvertices(dcel) - nedges(dcel) + nfaces(dcel) + 1

function mediatrix(p1::Point, p2::Point; turn = CounterClockWise)
    if p1 == p2
        throw(ArgumentError("The points must be different!"))
    end
    p = (p1 + Vec(p2)) / 2
    O = if turn == CounterClockWise
        Vec(0, 0, 1)
    else
        Vec(0, 0, -1)
    end
    q = p2 - p1
    Q = Vec(q.x, q.y, 0)
    V = O × Q  # CounterClockWise
    v = Vec(V.x, V.y)
    return Line(p, v)
end


include("geometry/convex_hull.jl")
include("geometry/triangulation.jl")
include("geometry/location.jl")
include("geometry/voronoi.jl")
