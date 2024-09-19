using Base.Order

using .BSTrees
using .DCELs

export makemonotone, makemonotone!, triangulatemonotone, triangulatemonotone!, tricolorvertices

# the return is an `Int` with the following codes
#
#  split   := -2
#  start   := -1
#  regular := 0
#  end     := 1
#  merge   := 2
#
function pointkind(v; tol=0)
    he = v |> halfedges |> first
    p = origin(prev(he))
    s = origin(he)
    n = origin(next(he))
    up = if n < s && p < s # all bellow
        -1
    elseif n > s && p > s # all above
        1
    else
        0 # regular
    end
    rot = if orientation(s, p, n; tol) == ClockWise # convex
        1
    elseif orientation(s, p, n; tol) == CounterClockWise # reflex
        2
    else # colinear
        0
    end
    rot * up
end

function makemonotone!(dcel::DCEL{T}; tol=0, cb=nothing) where {T}
    Q = sort!(vertices(dcel); by=point)
    o = Lt() do h1, h2
        h1 == h2 && return false
        e = segment(h2)
        e[1].y == e[2].y && return false # flat line
        o = let s = segment(h1)
            b, t = minmax(e[1].y, e[2].y)
            if t >= s[1].y >= b # choose the one that fits
                s[1]
            else
                s[2]
            end
        end
        c = (o.y - e[2].y) / (e[1].y - e[2].y)
        p = e[2] + c * (e[1] - e[2])
        isapprox(o.x, p.x; rtol=tol) && return false
        return o.x < p.x
    end
    tree = BSTree{HalfedgeHandle{T}}(; ord=o)
    helper = Dict{HalfedgeHandle{T},Tuple{HalfedgeHandle{T},Int}}()
    while !isempty(Q)
        v = pop!(Q)
        pk = pointkind(v; tol)
        e = halfedges(v)[1]
        isnothing(cb) || cb(dcel, tree, helper, v, pk)
        if pk == -1 # start
            insert!(tree, e)
            helper[e] = (e, pk)
        elseif pk == 1 # end
            pe = prev(e)
            if helper[pe][2] == 2
                hpe = helper[pe][1]
                d = insertdiagonal!(e, hpe)
                isnothing(cb) || cb(dcel, tree, helper, d, true)
            end
            delete!(tree, pe)
            delete!(helper, pe)
        elseif pk == -2 # split
            insert!(tree, e)
            ej = value(prevnode(tree, e)) # always exist
            isnothing(cb) || cb(dcel, tree, helper, ej, false)
            delete!(tree, e)
            hej = helper[ej][1]
            d = insertdiagonal!(e, hej)
            isnothing(cb) || cb(dcel, tree, helper, d, true)
            helper[ej] = (d, pk)
            insert!(tree, e)
            helper[e] = (e, pk)
        elseif pk == 2 # merge
            pe = prev(e)
            if helper[pe][2] == 2
                hpe = helper[pe][1]
                d = insertdiagonal!(e, hpe)
                isnothing(cb) || cb(dcel, tree, helper, d, true)
            end
            delete!(tree, pe)
            delete!(helper, pe)
            insert!(tree, e)
            ej = value(prevnode(tree, e))
            isnothing(cb) || cb(dcel, tree, helper, ej, false)
            delete!(tree, e)
            if helper[ej][2] == 2
                hej = helper[ej][1]
                d = insertdiagonal!(e, hej)
                isnothing(cb) || cb(dcel, tree, helper, d, true)
                helper[ej] = (d, pk)
            else
                helper[ej] = (e, pk)
            end
        elseif pk == 0 # regular
            s = segment(e)
            if s[1] > s[2] # inside to the right
                pe = prev(e)
                if helper[pe][2] == 2
                    hpe = helper[pe][1]
                    d = insertdiagonal!(e, hpe)
                    isnothing(cb) || cb(dcel, tree, helper, d, true)
                    helper[e] = (prev(e), pk)
                else
                    helper[e] = (e, pk)
                end
                delete!(tree, pe)
                delete!(helper, pe)
                insert!(tree, e)
            else
                insert!(tree, e)
                ej = value(prevnode(tree, e))
                isnothing(cb) || cb(dcel, tree, helper, ej, false)
                delete!(tree, e)
                if helper[ej][2] == 2
                    hej = helper[ej][1]
                    d = insertdiagonal!(e, hej)
                    isnothing(cb) || cb(dcel, tree, helper, d, true)
                    helper[ej] = (d, pk)
                else
                    helper[ej] = (e, pk)
                end
            end
        end
        isnothing(cb) || cb(dcel, tree, helper, v, -3)
    end
    dcel
end

makemonotone(dcel::DCEL{T}; tol=0, cb=nothing) where {T} = makemonotone!(deepcopy(dcel); tol, cb)

makemonotone(poly::Polygon{T}; tol=0, cb=nothing) where {T} = makemonotone!(DCEL(poly); tol, cb)


function triangulatemonotone!(dcel::DCEL{T}, f=1; tol=0, cb=nothing) where {T}
    fs = first(halfedges(dcel[Face, f]))
    bot, top = argmin(origin, fs), argmax(origin, fs)
    hedges = Vector{Tuple{HalfedgeHandle{T},Bool}}(undef, length(fs))
    p, n = prev(top), next(top)
    hedges[begin] = top, true # left
    hedges[end] = bot, false # right
    # make a list of halfedges by y-coordinate
    for i in 2:length(hedges)
        if argmax(origin, (p, n)) == n
            hedges[i] = n, true
            n = next(n)
        else
            hedges[i] = p, false
            p = prev(p)
        end
    end
    stack = Vector{Tuple{HalfedgeHandle{T},Bool}}()
    push!(stack, hedges[begin])
    push!(stack, hedges[begin+1])
    for i in 3:length(hedges)-1
        ui, uit = hedges[i]
        s, st = stack[end]
        if !isnothing(cb)
            cb(hedges, i, stack, dcel)
        end
        if st != uit
            if st # left
                for i in 2:length(stack)
                    v, vt = stack[i]
                    ui = insertdiagonal!(ui, v)
                end
            else # right
                st = uit
                for i in 2:length(stack)
                    v, vt = stack[i]
                    s = insertdiagonal!(v, ui)
                end
            end
            splice!(stack, 1:length(stack))
            push!(stack, (s, st))
            push!(stack, (ui, uit))
        else
            h, ht = pop!(stack)
            while true
                v, vt = pop!(stack)
                ori = orientation(origin(h), origin(ui), origin(v); tol)
                if (st && ori == ClockWise) || (!st && ori == CounterClockWise)
                    push!(stack, (v, vt))
                    push!(stack, (h, ht))
                    push!(stack, (ui, uit))
                    break
                end
                if st # left
                    v = insertdiagonal!(v, ui)
                else # right
                    ui = insertdiagonal!(ui, v)
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
    for i in 2:length(stack)-1
        v, vt = stack[i]
        if vt # left
            un = insertdiagonal!(un, v)
        else # right
            v = insertdiagonal!(v, un)
        end
    end
    return dcel
end

triangulatemonotone(dcel::DCEL{T}, f=1; tol=0, cb=nothing) where {T} = triangulatemonotone!(deepcopy(dcel), f; tol, cb)

function tricolorvertices(dcel::DCEL{T}) where {T}
    # check that every face, except the last one (infinite or external face) has three halfedges
    @assert all(
        f -> length(first(DCELs.Iterators.halfedges(f))) == 3,
        @view(faces(dcel)[begin:end-1])
    ) "Face with more or less than 3 edges!"
    fs = falses(nfaces(dcel)-1)
    # colors of the vertices
    vs = zeros(Int, nvertices(dcel))
    # which face is next
    backtrack = Vector{FaceHandle{T}}()
    # first halfedge
    h = halfedges(faces(dcel)[1])[1][1]
    # draw the firsts two vertices
    vs[vertex(h).i] = 1
    vs[vertex(next(h)).i] = 2
    # put the first face
    push!(backtrack, face(h))
    # the invariants in the code are that the current face will have at least two vertices
    # set to some color (so it will be solvable), and every face already visited will be
    # fully colored
    while !isempty(backtrack)
        f = pop!(backtrack)
        fs[f.i] = true
        h = first(first(DCELs.Iterators.halfedges(f)))
        colors = MVector(false, false, false)
        for _ in 1:3
            colors .= false
            vi = vertex(h).i
            if vs[vi] == 0
                colors[vs[vertex(prev(h)).i]] = true
                colors[vs[vertex(next(h)).i]] = true
                vs[vi] = findfirst(!, colors)
            end
            ft = face(twin(h))
            if ft.i != 0 && !fs[ft.i]
                push!(backtrack, ft)
            end
            h = next(h)
        end
    end
    vs
end
