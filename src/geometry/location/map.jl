using StaticArrays
using StructArrays

struct DataTrapezoid
    node::Int   # index to dag[node]
    leftp::Int  # index to dcel[Vertex, leftp], or 0 for BB
    rightp::Int # index to dcel[Vertex, rightp], or 0 for BB
    bottom::Int # index to dcel[Halfedge, bottom], or 0 for BB
    top::Int    # index to dcel[Halfedge, top], or 0 for BB
    # order is left-bottom, left-top, right-bottom, right-top
    neightbours::SVector{4, Int} # -1 for NA, 0 for BB
end

DataTrapezoid(node) = Trapezoid(node, 0, 0, 0, 0, 0, @SVector[-1, -1, -1, -1])

struct Map <: AbstractVector{DataTrapezoid}
    data::StructVector{DataTrapezoid, @NamedTuple{node::Vector{Int}, leftp::Vector{Int}, rightp::Vector{Int}, bottom::Vector{Int}, top::Vector{Int}, neightbours::Vector{SVector{4, Int}}}, Int}
end

struct Trapezoid
    map::Map
    i::Int
end

Map() = Map(StructVector{DataTrapezoid}((node = Int[], leftp = Int[], rightp = Int[], bottom = Int[], top = Int[], neightbours = SVector{4, Int}[])))

Base.size(map::Map) = size(map.data)

Base.@propagate_inbounds function Base.setindex!(map::Map, i::Int, v::Union{DataTrapezoid, Trapezoid})
    @boundscheck checkbounds(map.data, i)
    map.data.node[i] = v.node
    map.data.leftp[i] = v.leftp
    map.data.rightp[i] = v.rightp
    map.data.bottom[i] = v.bottom
    map.data.top[i] = v.top
    return map.data.neightbours[i] = v.neightbours
end

Base.@propagate_inbounds function Base.getindex(map::Map, i::Int)
    @boundscheck checkbounds(map.data, i)
    return Trapezoid(map, i)
end

function Base.setproperty!(trapezoid::Trapezoid, sym::Symbol)
    return if sym == :node
        trapezoid.map.data[trapezoid.i].node
    elseif sym == :leftp
        trapezoid.map.data[trapezoid.i].leftp
    elseif sym == :rightp
        trapezoid.map.data[trapezoid.i].rightp
    elseif sym == :bottom
        trapezoid.map.data[trapezoid.i].bottom
    elseif sym == :top
        trapezoid.map.data[trapezoid.i].top
    elseif sym == :neightbours
        trapezoid.map.data[trapezoid.i].neightbours
    else
        getfield(trapezoid.map.data[trapezoid.i], sym)
    end
end

function Base.setproperty!(trapezoid::Trapezoid, sym::Symbol, v)
    return if sym == :node
        trapezoid.map.data.node[trapezoid.i] = v
    elseif sym == :leftp
        trapezoid.map.data.leftp[trapezoid.i] = v
    elseif sym == :rightp
        trapezoid.map.data.rightp[trapezoid.i] = v
    elseif sym == :bottom
        trapezoid.map.data.bottom[trapezoid.i] = v
    elseif sym == :top
        trapezoid.map.data.top[trapezoid.i] = v
    elseif sym == :neightbours
        trapezoid.map.data.neightbours[trapezoid.i] = v
    else
        setfield!(trapezoid.map.data, sym, v)
    end
end

function Base.push!(map::Map, trapezoid::Trapezoid)
    push!(map.data, trapezoid)
    return map.data[end]
end

function Base.replace!(t1::Trapezoid, t2::Trapezoid)
    t1.data[t1.i] = trapezoid
    return t1
end
