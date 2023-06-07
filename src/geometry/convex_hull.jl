export jarvismarch, quickhull, grahamscan

# receives sorted points
function _half_jarvismarch!(chull, spoints; cb=nothing, tol=0)
    npoints = length(spoints)
    prev = 1
    next = npoints
    while prev != npoints
        push!(chull, spoints[prev])
        # we only need to check for those points _after_ the point where we are
        for i ∈ prev+1:npoints
            if !isnothing(cb)
                cb(chull, spoints[prev], spoints[next], spoints[i])
            end
            ori = orientation(spoints[prev], spoints[next], spoints[i]; tol)
            if ori == CounterClockWise || ori == Colinear
                next = i
            end
        end
        prev = next
        next = npoints
    end
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
 4. `testing` is the current guess for the vertex in the convex hull

See also [`quickhull`](@ref), [`grahamscan`](@ref).

# Examples

```julia-repl
julia> points = [Point(rand(2)) for _ in 1:10^4];

julia> chull = jarvismarch(points; cb=(c,p,n,t)->nothing, tol=sqrt(eps()));
```
"""
function jarvismarch(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)
    if length(points) ≤ 3
        return radialsort(points)
    end
    spoints = sort(points; by=p -> (p.x, -p.y))
    chull = eltype(points)[]
    _half_jarvismarch!(chull, spoints; cb, tol)
    reverse!(spoints)
    _half_jarvismarch!(chull, spoints; cb, tol)
    return chull
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
function quickhull(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)
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
    return [a] ∪ Qa ∪ [b] ∪ Qb
end


@doc raw"""
    grahamscan(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)

Compute the convex hull of `points` given the tolerance `tol` using the Graham Scan
algorithm with complexity of $O(n\log n)$.

Returns a vector of `Points` that are in clockwise order and correspond to
$ConvexHull(points)$

If the callback `cb` is different from `nothing` then it will be called on every
iteration as follow

    cb(chull, prev, next, testing)

Where

 1. `chull` is the collection of points already in the convex hull
 2. `prev` is the previous vertex in the convex hull
 3. `next` is the current best guess for vertex in the convex hull
 4. `testing` is the current guess for the vertex in the convex hull
    cb(chull, prev, next, testing)

See also [`jarvismarch`](@ref), [`quickhull`](@ref).

# Examples

```julia-repl
julia> points = [Point(rand(2)) for _ in 1:10^4];

julia> chull = grahamscan(points; cb=(c,p,n,t)->nothing, tol=sqrt(eps()));
```
"""
function grahamscan(points::AbstractArray{<:Union{Point{2},Vec{2}}}; cb=nothing, tol=0)
    if length(points) ≤ 3
        return radialsort(points)
    end
    p = argmin(p -> (p.y, p.x), points)
    spoints = radialsort(points, p)
    spoints = spoints[filter(
        i -> orientation(spoints[i-1], p, spoints[i]) != Colinear,
        3:length(spoints)
    )]
    chull = [p, spoints[2]]
    for i ∈ 3:length(spoints)
        if orientation(chull[end], chull[end-1], spoints[i]) == ClockWise
            push!(chull, spoints[i])
        else
            pop!(chull)
        end
    end
    return chull
end
