export radialsort, radialsort!, eulerschar

# sorting points

using Statistics #, DataStructures

Statistics.mean(points::AbstractArray{<:Point{2}}) = Point(Statistics.mean(p -> p.coords, points))
Statistics.mean(points::AbstractArray{<:Vec{2}}) = Vec(Statistics.mean(p -> p.coords, points))

function radialsort!(points::AbstractArray{<:Union{Point{2},Vec{2}}}, pivot=nothing; rev=false)
    if isnothing(pivot)
        pivot = mean(points)
    end
    # sort by angle first and by length second
    sort!(points; by=p -> (atan(-(p - pivot).y, -(p - pivot).x), norm(p - pivot)), rev)
end

radialsort(points, pivot=nothing) = radialsort!(copy(points), pivot)

eulerschar(dcel) = nvertices(dcel) - nedges(dcel) + nfaces(dcel)

include("geometry/convex_hull.jl")
include("geometry/triangulation.jl")
# include("geometry/location.jl")
include("geometry/voronoi.jl")
