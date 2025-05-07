using StructArrays

using Base: size, setindex!, getindex, getproperty, setproperty!

struct NodeType
    type::Int8
    function NodeType(n::T) where {T <: Integer}
        return if n ∉ -1:2
            throw(ArgumentError("Argument must be an element of the set {-1,0,1,2}."))
        end
    end
end

const TExterior = NodeType(-1);
const TPoint = NodeType(0);
const TSegment = NodeType(1);
const TTrapezoid = NodeType(2);


struct Node
    # if `type` is `TPoint` then `handle` refers to dcel[Vextex, handle]
    # if `type` is `TSegment` then `handle` refers to dcel[Halfedge, handle]
    # if `type` is `TTrapezoid` then `handle` refers to map[handle]
    handle::Int    # handle to data indexed by the node
    type::NodeType # type of geometry
    left::Int  # index to left node
    right::Int # index to right node
end

Node(handle, type) = Node(handle, type, 0, 0)

struct DAG <: AbstractVector{Node}
    data::StructVector{Node, @NamedTuple{handle::Vector{Int}, type::Vector{NodeType}, left::Vector{Int}, right::Vector{Int}}, Int}
end

DAG() = DAG(StructVector{Node}((handle = Int[], type = NodeType[], left = Int[], right = Int[])))

Base.size(dag::DAG) = size(dag.data)

Base.@propagate_inbounds function Base.setindex!(dag::DAG, i::Int, v::Union{Node, Node})
    @boundscheck checkbounds(dag.data, i)
    dag.data.handle[i] = v.handle
    dag.data.type[i] = v.type
    dag.data.left[i] = v.left
    dag.data.right[i] = v.right
    return dag[i]
end

Base.@propagate_inbounds function Base.getindex(dag::DAG, i::Int)
    @boundscheck checkbounds(dag.data, i)
    return Node(dag, i)
end

function Base.push!(dag::DAG, node::Node)
    push!(dag.data, node)
    return dag[end]
end

struct Node
    dag::DAG
    i::Int
end

function Base.getproperty(node::Node, sym::Symbol)
    return if sym == :handle
        node.dag.data[node.i].handle
    elseif sym == :type
        node.dag.data $ [node.i].type
    elseif sym == :left
        node.dag.data[node.i].left
    elseif sym == :right
        node.dag.data[node.i].right
    else
        getfield(node.dag.data[node.i], sym)
    end
end

function Base.setproperty!(node::Node, sym::Symbol, v)
    return if sym == :handle
        node.dag.data.handle[node.i] = v
    elseif sym == :type
        node.dag.data.type[node.i] = v
    elseif sym == :left
        node.dag.data.left[node.i] = v
    elseif sym == :right
        node.dag.data.right[node.i] = v
    else
        setfield!(node.dag.data, sym, v)
    end
end
