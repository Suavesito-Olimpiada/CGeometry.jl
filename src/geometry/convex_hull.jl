using .DCELs

export jarvismarch, quickhull, grahamscan, boundingbox

function boundingbox(pts)
    X = extrema(p -> p.x, pts)
    Y = extrema(p -> p.y, pts)
    return X, Y
end

boundingbox(dcel::DCEL) = boundingbox(DCELs.Iterators.points(dcel))

# receives points and extrema
function _half_jarvismarch!(chull, points, ip, iq; cb = nothing, tol = 0)
    prev = ip
    next = iq
    while prev != iq
        push!(chull, points[prev])
        for i in eachindex(points)
            # the three are different
            if i == prev || i == next || prev == next
                continue
            end
            if !isnothing(cb)
                cb(chull, points, prev, next, i)
            end
            ori = orientation(points[prev], points[next], points[i]; tol)
            if ori == CounterClockWise
                next = i
            elseif ori == Colinear
                # keep the furthest colinear point
                n = norm(points[prev] - points[i])
                m = norm(points[prev] - points[next])
                if n > m
                    next = i
                end
            end
        end
        prev = next
        next = iq
    end
    return
end

@doc raw"""
    jarvismarch(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)

Compute the convex hull of `points` given the tolerance `tol` using the Gift Wrapping
algorithm with complexity $O(nh)$.

Returns a vector of `Points` that are in clockwise order and correspond to
$ConvexHull(points)$.

If the callback `cb` is different from `nothing` then it will be called on every
iteration as follow

    cb(chull, prev, next, testing)

Where

 1. `chull` is the collection of points already in the convex hull
 2. `prev` is the previous vertex in the convex hull
 3. `next` is the current best guess for vertex in the convex hull
 4. `index` is the current vertex testing to be in the convex hull

See also [`quickhull`](@ref), [`grahamscan`](@ref).

# Examples

```julia-repl
julia> points = [Point(rand(2)) for _ in 1:10^4];

julia> chull = jarvismarch(points; cb=nothing, tol=sqrt(eps()));
```
"""
function jarvismarch(points::AbstractArray{<:Union{Point{2}, Vec{2}}}; cb = nothing, tol = 0)
    if length(points) ≤ 3
        return radialsort(points)
    end
    _, ip = findmin(p -> (p.x, -p.y), points)
    _, iq = findmax(p -> (p.x, -p.y), points)
    chull = eltype(points)[]
    _half_jarvismarch!(chull, points, ip, iq; cb, tol)
    _half_jarvismarch!(chull, points, iq, ip; cb, tol)
    return Polygon(reverse!(chull))
end


_leftof(points, ds; tol) = filter(p -> orientation(ds, p; tol) == CounterClockWise, points)

_furthesttoline(points, ds) = argmin(p -> (p - ds[1]) × (ds[2] - ds[1]), points)

function _quickhalfhull(points, a, b; cb, tol)
    if length(points) == 0
        return points
    end
    c = _furthesttoline(points, DirSegment(a, b))
    A = _leftof(points, DirSegment(a, c); tol)
    B = _leftof(points, DirSegment(c, b); tol)
    if !isnothing(cb)
        cb([a, c, b], A, B)
    end
    Qa = _quickhalfhull(A, a, c; cb, tol)
    Qb = _quickhalfhull(B, c, b; cb, tol)
    return Qa ∪ [c] ∪ Qb
end

@doc raw"""
    quickhull(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)

Compute the convex hull of `points` given the tolerance `tol` using the Quick Hull algorithm
with mean complexity of $O(n\log n)$ and worst case complexity of $O(n^2)$.

Returns a vector of `Points` that are in clockwise order and correspond to
$ConvexHull(points)$

If the callback `cb` is different from `nothing` then it will be called on every
iteration as follow

    cb(furthest, left, right)

Where

 1. `furthest` is the collection of points `[a,c,b]` where $a̅b̅$ is the segment that divides
    the plane in two and `c` is the furthest point to the left of $a̅b̅$. In
    the first call it will only be `[a,b]`, the extrema of `points` on the $x$ axis
 2. `left` is the collection of points to the left of the triangle
 3. `right` is the collection of points to the right of the triangle

See also [`jarvismarch`](@ref), [`grahamscan`](@ref).

# Examples

```julia-repl
julia> points = [Point(rand(2)) for _ in 1:10^4];

julia> chull = quickhull(points; cb=(f,l,r)->nothing, tol=sqrt(eps()));
```
"""
function quickhull(points::AbstractArray{<:Union{Point{2}, Vec{2}}}; cb = nothing, tol = 0)
    if length(points) ≤ 3
        return radialsort(points)
    end
    a, b = argmin(p -> (p.x, -p.y), points), argmax(p -> (p.x, -p.y), points)
    A = _leftof(points, DirSegment(a, b); tol)
    B = _leftof(points, DirSegment(b, a); tol)
    if !isnothing(cb)
        cb([a, b], A, B)
    end
    Qa = _quickhalfhull(A, a, b; cb, tol)
    Qb = _quickhalfhull(B, b, a; cb, tol)
    return Polygon(reverse!([a] ∪ Qa ∪ [b] ∪ Qb))
end


@doc raw"""
    grahamscan(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)

Compute the convex hull of `points` given the tolerance `tol` using the Graham Scan
algorithm with complexity of $O(n\log n)$.

Returns a vector of `Points` that are in clockwise order and correspond to
$ConvexHull(points)$

If the callback `cb` is different from `nothing` then it will be called on every
iteration as follow

    cb(chull, pivot, points, index)

Where

 1. `chull` is the collection of points already in the convex hull
 2. `pivot` is the y⁻-most, x⁻-most point of the original set
 3. `points` is the vertices considered to be in the convex hull
 4. `index` is the current veretex in points to be tested

See also [`jarvismarch`](@ref), [`quickhull`](@ref).

# Examples

```julia-repl
julia> points = [Point(rand(2)) for _ in 1:10^4];

julia> chull = grahamscan(points);
```
"""
grahamscan(points::AbstractArray{<:Union{Point{2}, Vec{2}}}; cb = nothing, tol = 0) = grahamscan!(copy(points); cb, tol)
function grahamscan!(points::AbstractArray{<:Union{Point{2}, Vec{2}}}; cb = nothing, tol = 0)
    npoints = length(points)
    if npoints ≤ 3
        # sort in counter-clockwise order
        pivot = mean(points)
        return radialsort!(points; pivot)
    end
    # find the y⁻-most, x⁻-most point (bottom-left)
    pivot = argmin(p -> (p.x, -p.y), points)
    ip = findfirst(==(pivot), points)
    deleteat!(points, ip)
    # sort by angle first and by reverse-norm second
    radialsort!(points; pivot, direction = π, turn = ClockWise, revnorm = true)
    # only kept the furthest (first) points if several have the same angle
    unique!(p -> atan(-(p - pivot).y, -(p - pivot).x), points)
    # the stack
    chull = [pivot]
    for i in eachindex(points)
        p = points[i]
        while length(chull) > 2 && orientation(p, chull[end - 1], chull[end]) != ClockWise
            pop!(chull)
        end
        push!(chull, p)
        if !isnothing(cb)
            cb(chull, points, i)
        end
    end
    return Polygon(reverse!(chull))
end
