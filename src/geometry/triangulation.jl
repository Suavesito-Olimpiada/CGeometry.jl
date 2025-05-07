using Base.Order

using .BSTrees
using .DCELs

export makemonotone, makemonotone!, triangulatemonotone, triangulatemonotone!, tricolorvertices

@enum PointType VSplit = -2 VStart = -1 VRegular = 0 VEnd = 1 VMerge = 2

# the return is an `Int` with the following codes
#
#  split   := -2
#  start   := -1
#  regular := 0
#  end     := 1
#  merge   := 2
#
function pointkind(he; tol = 0)
    p = tail(prev(he))
    s = tail(he)
    n = tail(next(he))
    up = if n < s && p < s # all bellow
        -1
    elseif n > s && p > s # all above
        1
    else
        0 # regular
    end
    ori = orientation(s, p, n; tol)
    rot = if ori == ClockWise # convex
        1
    elseif ori == CounterClockWise # reflex
        2
    else # colinear
        0
    end
    return PointType(rot * up)
end

function makemonotone!(dcel::DCEL{T}; tol = 0, cb = nothing) where {T}
    # sweep line priority queue, it will not change
    Q = sort!(vertices(dcel); by = point)
    n = length(Q)
    # sorting of half-edges
    # TODO: explain the sorting and draw
    order = Lt() do h1, h2
        # if both half-edges are the same the first is not less than
        # the second
        h1 == h2 && return false
        e = segment(h2)
        # flat line
        e[1].y == e[2].y && return false
        o = let s = segment(h1)
            b, t = minmax(e[1].y, e[2].y)
            if b <= s[1].y <= t # choose the one that fits
                s[1]
            else
                s[2]
            end
        end
        c = (o.y - e[2].y) / (e[1].y - e[2].y)
        p = e[2] + c * (e[1] - e[2])
        isapprox(o.x, p.x; rtol = tol) && return false
        return o.x < p.x
    end
    tree = BSTree{Halfedge{T}}(; ord = order)
    helper = Dict{Halfedge{T}, Tuple{Halfedge{T}, PointType}}()
    face = first(faces(dcel))
    while !isempty(Q)
        n -= 1
        v = pop!(Q)
        e = halfedge(v, face)
        pk = pointkind(e; tol)
        isnothing(cb) || cb(n, dcel, tree, helper, v, Int(pk))
        if pk == VStart # start
            insert!(tree, e)
            helper[e] = (e, pk)
        elseif pk == VEnd # end
            pe = e.prev
            if helper[pe][2] == VMerge
                hpe = helper[pe][1]
                d = DCELs.connect_diagonal(e.prev, hpe)
                isnothing(cb) || cb(n, dcel, tree, helper, d, true)
            end
            delete!(tree, pe)
            delete!(helper, pe)
        elseif pk == VSplit # split
            insert!(tree, e)
            ej = value(prevnode(tree, e)) # always exist
            isnothing(cb) || cb(n, dcel, tree, helper, ej, false)
            delete!(tree, e)
            hej = helper[ej][1]
            d = DCELs.connect_diagonal(e.prev, hej)
            isnothing(cb) || cb(n, dcel, tree, helper, d, true)
            helper[ej] = (d, pk)
            insert!(tree, e)
            helper[e] = (e, pk)
        elseif pk == VMerge # merge
            pe = e.prev
            if helper[pe][2] == VMerge
                hpe = helper[pe][1]
                d = DCELs.connect_diagonal(e.prev, hpe)
                isnothing(cb) || cb(n, dcel, tree, helper, d, true)
            end
            delete!(tree, pe)
            delete!(helper, pe)
            insert!(tree, e)
            ej = value(prevnode(tree, e))
            isnothing(cb) || cb(n, dcel, tree, helper, ej, false)
            delete!(tree, e)
            if helper[ej][2] == VMerge
                hej = helper[ej][1]
                d = DCELs.connect_diagonal(e.prev, hej)
                isnothing(cb) || cb(n, dcel, tree, helper, d, true)
                helper[ej] = (d, pk)
            else
                helper[ej] = (e, pk)
            end
        elseif pk == VRegular # regular
            s = segment(e)
            if s[1] > s[2] # inside is to the right
                pe = e.prev
                if helper[pe][2] == VMerge
                    hpe = helper[pe][1]
                    d = DCELs.connect_diagonal(e.prev, hpe)
                    isnothing(cb) || cb(n, dcel, tree, helper, d, true)
                    helper[e] = (e.prev, pk)
                else
                    helper[e] = (e, pk)
                end
                delete!(tree, pe)
                delete!(helper, pe)
                insert!(tree, e)
            else
                insert!(tree, e)
                ej = value(prevnode(tree, e))
                isnothing(cb) || cb(n, dcel, tree, helper, ej, false)
                delete!(tree, e)
                if helper[ej][2] == VMerge
                    hej = helper[ej][1]
                    d = DCELs.connect_diagonal(e.prev, hej)
                    isnothing(cb) || cb(n, dcel, tree, helper, d, true)
                    helper[ej] = (d, pk)
                else
                    helper[ej] = (e, pk)
                end
            end
        end
        isnothing(cb) || cb(n, dcel, tree, helper, v, -3)
    end
    return dcel
end

makemonotone(dcel::DCEL{T}; tol = 0, cb = nothing) where {T} = makemonotone!(deepcopy(dcel); tol, cb)

makemonotone(poly::Polygon{T}; tol = 0, cb = nothing) where {T} = makemonotone!(DCEL(poly); tol, cb)


function triangulatemonotone!(dcel::DCEL{T}, f = 1; tol = 0, cb = nothing) where {T}
    fs = DCELs.Iterators.halfedges(dcel[Face, f])
    bot, top = argmin(tail, fs), argmax(tail, fs)
    hedges = Vector{Tuple{Halfedge{T}, Bool}}(undef, length(fs))
    p, n = top.prev, top.next
    hedges[begin] = top, true # left
    hedges[end] = bot, false # right
    # make a list of halfedges by y-coordinate
    for i in 2:length(hedges)
        if argmax(tail, (p, n)) == n
            hedges[i] = n, true
            n = n.next
        else
            hedges[i] = p, false
            p = p.prev
        end
    end
    stack = Vector{Tuple{Halfedge{T}, Bool}}()
    push!(stack, hedges[begin])
    push!(stack, hedges[begin + 1])
    for i in 3:(length(hedges) - 1)
        ui, uit = hedges[i]
        s, st = stack[end]
        if !isnothing(cb)
            cb(hedges, i, stack, dcel)
        end
        if st != uit
            if st # left
                for i in 2:length(stack)
                    v, vt = stack[i]
                    ui = DCELs.connect_diagonal(ui.prev, v)
                end
            else # right
                st = uit
                for i in 2:length(stack)
                    v, vt = stack[i]
                    s = DCELs.connect_diagonal(v.prev, ui)
                end
            end
            splice!(stack, 1:length(stack))
            push!(stack, (s, st))
            push!(stack, (ui, uit))
        else
            h, ht = pop!(stack)
            while true
                v, vt = pop!(stack)
                ori = orientation(tail(h), tail(ui), tail(v); tol)
                if (st && ori == ClockWise) || (!st && ori == CounterClockWise)
                    push!(stack, (v, vt))
                    push!(stack, (h, ht))
                    push!(stack, (ui, uit))
                    break
                end
                if st # left
                    v = DCELs.connect_diagonal(v.prev, ui)
                else # right
                    ui = DCELs.connect_diagonal(ui.prev, v)
                end
                h, ht = v, vt
                if isempty(stack)
                    push!(stack, (h, ht))
                    push!(stack, (ui, uit))
                    break
                end
            end
        end
    end
    un, unt = hedges[end]
    for i in 2:(length(stack) - 1)
        v, vt = stack[i]
        if vt # left
            un = DCELs.connect_diagonal(un.prev, v)
        else # right
            v = DCELs.connect_diagonal(v.prev, un)
        end
    end
    return dcel
end

triangulatemonotone(dcel::DCEL{T}, f = 1; tol = 0, cb = nothing) where {T} = triangulatemonotone!(deepcopy(dcel), f; tol, cb)

function tricolorvertices(dcel::DCEL{T}) where {T}
    # check that every face, except the last one (infinite or external face) has three halfedges
    if !all(f -> length(f) == 3, DCELs.Iterators.faces(dcel))
        throw(ArgumentError("There are non triangular faces!"))
    end
    nf = nfaces(dcel)
    fs = falses(nf)
    # colors of the vertices
    vs = zeros(Int, length(dcel.vertices))
    # which face is next
    backtrack = Vector{Face{T}}()
    # first halfedge
    h = dcel[Halfedge, 1]
    # draw the firsts two vertices
    vs[h.vertex.i] = 1
    vs[h.next.vertex.i] = 2
    # put the first face
    push!(backtrack, h.face)
    # the invariants in the code are that the current face will have at least two vertices
    # set to some color (so it will be solvable), and every face already visited will be
    # fully colored
    while !isempty(backtrack)
        f = pop!(backtrack)
        fs[f.i] = true
        h = first(DCELs.Iterators.halfedges(f))
        colors = MVector(false, false, false)
        for _ in 1:3
            colors .= false
            vi = h.vertex.i
            if vs[vi] == 0
                colors[vs[h.prev.vertex.i]] = true
                colors[vs[h.next.vertex.i]] = true
                vs[vi] = findfirst(!, colors)
            end
            ft = h.twin.face
            if ft.i != 0 && !fs[ft.i]
                push!(backtrack, ft)
            end
            h = h.next
        end
    end
    return vs
end
