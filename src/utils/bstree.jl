module BSTrees

using Base.Iterators
using Base.Order

export BSTree, find, value, nextnode, prevnode

@enum Traverse begin
    left = -1
    self = 0
    right = 1
end

@inline Base.:(~)(traverse::Traverse) = ifelse(traverse == left, right, ifelse(traverse == right, left, traverse))

# height > 0 are elements of the tree that can be deleted an inserted to,
# height == 0 means back references to the node which would have it as
# ancestor or successor, or to itself in case of extrema
mutable struct BSTreeNode{T}
    height::Int
    value::T
    left::BSTreeNode{T}
    right::BSTreeNode{T}
    BSTreeNode{T}(value::T) where {T} = new{T}(0, value)
end

height(node::BSTreeNode) = node.height
value(node::BSTreeNode) = node.value

isleaf(node::BSTreeNode) = node.height == 1
nchildren(node::BSTreeNode) = !iszero(height(node.left)) + !iszero(height(node.right))
balance(node::BSTreeNode) = height(node.right) - height(node.left)

updateheight!(node::BSTreeNode) = (node.height = max(height(node.left), height(node.right)) + 1)

@inline function Base.getindex(node::BSTreeNode, traverse::Traverse)
    if traverse == self
        node
    elseif traverse == left
        node.left
    elseif traverse == right
        node.right
    end
end

@inline function Base.setindex!(node::BSTreeNode, child::BSTreeNode, traverse::Traverse)
    if traverse == left
        node.left = child
    elseif traverse == right
        node.right = child
    else
        throw(ArgumentError("$node[$traverse] = $child does not makes sense."))
    end
end

mutable struct BSTree{T,Ord<:Ordering}
    count::Int
    root::Union{BSTreeNode{T}, Nothing}
    o::Ord
    BSTree{T}(;ord::Ord=Forward) where {T,Ord} = new{T,Ord}(0,nothing,ord)
end

height(tree::BSTree) = isempty(tree) ? 0 : height(tree.root)
Base.isempty(tree::BSTree) = tree.count == 0

# left := minimum, right := maximum
function extrema(node::BSTreeNode, side::Traverse)
    height(node) == 0 && return node
    while height(node[side]) != 0
        node = node[side]
    end
    return node
end

# left := prev, right := next
function nearest(node::BSTreeNode, side::Traverse)
    height(node) == 0 && return node
    tr = side
    while height(node[tr]) != 0
        node = node[tr]
        tr = ~side
    end
    return node
end

# left := prev, right := next
function pathtonearest(node::BSTreeNode{T}, side::Traverse) where {T}
    path = Vector{Tuple{BSTreeNode{T},Traverse}}(undef, height(node))
    head = 0
    tr = side
    while height(node) != 0
        head += 1
        path[head] = (node, tr)
        node = node[tr]
        tr = ~side
    end
    resize!(path, head)
    return path
end

function prevnode(tree::BSTree{T}, val::T) where {T}
    node = find(tree, val)
    isnothing(node) && return nothing
    nnode = nearest(node, left)
    if nnode == node
        if nnode.left.left == nnode.left
            return nothing
        else
            return nnode.left.left
        end
    end
    return nnode
end

function nextnode(tree::BSTree{T}, val::T) where {T}
    node = find(tree, val)
    isnothing(node) && return nothing
    nnode = nearest(node, right)
    if nnode == node
        if nnode.right.right == nnode.right
            return nothing
        else
            return nnode.right.right
        end
    end
    return nnode
end

find(tree::BSTree{T}, val::T) where {T} = findnode(tree.root, BSTreeNode{T}(val), tree.o)

Base.in(val::T, tree::BSTree{T}) where {T} = !isnothing(findnode(tree.root, BSTreeNode{T}(val), tree.o))

findnode(tree::BSTree{T,Ord}, node::BSTreeNode{T}) where {T,Ord} = findnodepath(tree.root, node, tree.o)
function findnode(root::BSTreeNode{T}, node::BSTreeNode{T}, ord) where {T}
    next = root
    local traverse::Traverse
    while height(next) != 0
        traverse = if value(node) == value(next)
            self
        else
            lt(ord, value(node), value(next)) ? left : right
        end
        next = if traverse ∈ (left, right) && height(next[traverse]) != 0
            next[traverse]
        else
            break
        end
    end
    return traverse == self ? next : nothing
end

findnodepath(tree::BSTree{T,Ord}, node::BSTreeNode{T}) where {T,Ord} = findnodepath(tree.root, node, tree.o)
function findnodepath(root::BSTreeNode{T}, node::BSTreeNode{T}, ord) where {T}
    path = Vector{Tuple{BSTreeNode{T},Traverse}}(undef, height(root))
    head = 0
    next = root
    local traverse::Traverse
    while height(next) != 0
        head += 1
        traverse = if value(node) == value(next)
            self
        else
            lt(ord, value(node), value(next)) ? left : right
        end
        path[head] = (next, traverse)
        next = if traverse ∈ (left, right) && height(next[traverse]) != 0
            next[traverse]
        else
            break
        end
    end
    resize!(path, head)
    return path, traverse == self
end

function rotate!(node::BSTreeNode, traverse::Traverse)
    child = node[~traverse]
    node[~traverse] = child[traverse]
    child[traverse] = node
    updateheight!(node)
    updateheight!(child)
    return child
end

function fixrotate!(root, trc, trgc)
    if trc != trgc
        root[trc] = rotate!(root[trc], ~trgc)
    end
    return rotate!(root, ~trc)
end

function initroot!(node::BSTreeNode{T}) where {T}
    nodemin = BSTreeNode{T}(value(node))
    nodemin.left = nodemin
    nodemin.right = node
    nodemax = BSTreeNode{T}(value(node))
    nodemax.left = node
    nodemax.right = nodemax
    node.left = nodemin
    node.right = nodemax
    updateheight!(node)
    return node
end

@inline function insertleaf!(root::BSTreeNode{T}, leaf::BSTreeNode{T}, traverse::Traverse) where {T}
    nleaf = BSTreeNode{T}(value(leaf))
    nleaf[traverse] = leaf
    nleaf[~traverse] = root
    leaf[traverse] = root[traverse]
    leaf[traverse][~traverse] = leaf
    leaf[~traverse] = nleaf
    root[traverse] = leaf
    updateheight!(leaf)
    return leaf
end

function insertnode!(tree::BSTree{T,Ord}, node::BSTreeNode{T}) where {T,Ord}
    if isempty(tree)
        tree.root = initroot!(node)
        tree.count += 1
        return tree.root
    end
    path, eq = findnodepath(tree.root, node, tree.o)
    eq && return path[end][1]
    parent, traverse = path[end]
    node = insertleaf!(parent, node, traverse)
    tree.count += 1
    child = node
    reverse!(path)
    # this loop goes from the parent of the new node to root
    for i ∈ eachindex(path)
        parent, trc = path[i] # _tr_averse to _c_hild
        parent[trc] = child
        updateheight!(parent)
        if abs(balance(parent)) > 1
            trgc = path[i-1][2] # _tr_averse to _g_rand_c_hild
            child = fixrotate!(parent, trc, trgc)
        else
            child = parent
        end
    end
    tree.root = child
    return node
end

function deletenode!(tree::BSTree{T,Ord}, node::BSTreeNode{T}) where {T,Ord}
    isempty(tree) && return
    path, eq = findnodepath(tree.root, node, tree.o)
    eq || throw(KeyError(value(node)))
    node, _ = path[end]
    isroot = length(path) == 1
    nch = nchildren(node)
    child = if nch == 2 # inside tree
        npath = pathtonearest(node, right)
        pnext, trgc = npath[end-1]
        next, trc = npath[end]
        prev = nearest(node, left)
        prev.right.right = next
        pnext[trgc] = next.right
        next.right = node.right
        next.left = node.left
        path[end] = (next, right)
        append!(path, @view(npath[begin+1:end-1]))
        updateheight!(pnext)
        updateheight!(next)
        pnext[trgc]
    elseif nch == 1 # almost leaf
        side = height(node.right) == 0 ? left : right # side with node
        node[~side][side] = node[side]
        node[side][~side] = node[~side]
        pop!(path)
        node[side]
    else # leaf
        if isroot
            pop!(path)
            nothing
        else
            parent, ptr = path[end-1]
            node[ptr][~ptr] = parent
            pop!(path)
            node[ptr]
        end
    end
    reverse!(path)
    for i ∈ eachindex(path)
        parent, trc = path[i] # _tr_averse to _c_hild
        parent[trc] = child
        updateheight!(parent)
        if abs(balance(parent)) > 1
            trgc = Traverse(balance(parent[~trc])) in (self, ~trc) ? ~trc : trc
            child = fixrotate!(parent, ~trc, trgc)::BSTreeNode
        else
            child = parent
        end
    end
    tree.root = child
    tree.count -= 1
    return node
end

Base.insert!(tree::BSTree{T,Ord}, val::T) where {T,Ord} = value(insertnode!(tree, BSTreeNode{T}(val)))

Base.delete!(tree::BSTree{T,Ord}, val::T) where {T,Ord} = value(deletenode!(tree, BSTreeNode{T}(val)))

function _print(io, node::BSTreeNode, height, deep, active, leaves)
    deep > height && return
    node.height == 0 || _print(io, node.right, height, deep+1, active, leaves)
    if deep == 1
        println(io, "   ", node.value, "(", node.height, ")")
    else
        active[deep] = ~active[deep]
        str = reduce((str, b) -> str*(b ? "|   " : "    "), active[begin:deep-1]; init="")
        if node.height != 0
            println(io, str*"+-- ", node.value, "(", node.height, ")")
        elseif leaves
            println(io, str*"+-- ", node.left.value, ":", node.value, ":", node.right.value, "(", node.height, ")")
        end
    end
    node.height == 0 || _print(io, node.left, height, deep+1, active, leaves)
    nothing
end

printfull(node::BSTreeNode) = _print(stdout, node, node.height+1, 1, falses(node.height+1), true)

printfull(tree::BSTree) = _print(stdout, tree.root, tree.root.height+1, 1, falses(tree.root.height+1), true)

function Base.show(io::IO, ::MIME"text/plain", tree::BSTree{T,Ord}) where {T,Ord}
    summary(io, tree)
    println(io, " with ", tree.count, " elements")
    isempty(tree) || _print(io, tree.root, tree.root.height, 1, falses(tree.root.height), false)
    nothing
end

function Base.show(io::IO, ::MIME"text/plain", node::BSTreeNode{T}) where {T}
    summary(io, node)
    println(io, " with height ", node.height)
    _print(io, node, node.height, 1, falses(node.height), false)
    nothing
end

end
