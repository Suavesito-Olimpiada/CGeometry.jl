using Base.Order

using .BSTrees
using .DCELs

export makemonotone, makemonotone!

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
    # pks = pointkind.(Q; tol)
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
        c = (o.y - e[2].y )/(e[1].y - e[2].y)
        p = e[2] + c*(e[1] - e[2])
        isapprox(o.x, p.x; rtol=tol) && return false
        return o.x < p.x
    end
    tree = BSTree{HalfedgeHandle{T}}(;ord=o)
    helper = Dict{HalfedgeHandle{T}, Tuple{HalfedgeHandle{T},Int}}()
    while !isempty(Q)
        v = pop!(Q)
        # pk = pop!(pks)
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
                    helper[ej]= (d, pk)
                else
                    helper[ej]= (e, pk)
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
    bot, top = extrema(origin, halfedges(dcel[Face, f])[1])
    hedges = Vector{Tuple{HalfedgeHandle{T},Bool}}(undef, nhalfedges(dcel))
    p, n = prev(top), next(top)
    hedges[begin] = top, true # left
    hedges[end] = bot, false # right
    for i in 2:length(hedges)
        if argmax(origin, (p, n)) == n
            hedges[i] = n, true
            n = next(n)
        else
            hedges[i] = n, false
            p = prev(p)
        end
    end
    stack = Vector{Tuple{HalfedgeHandle{T},Bool}}()
    push!(stack, hedges[begin])
    push!(stack, hedges[begin+1])
    for j in 3:length(hedges)-1
        uj, ujt = hedges[j]
        s, st = stack[end]
        if st != ujt
            while length(stack) != 1
                v, vt = pop!(stack)
                uj = vt ? insertdiagonal!(v, uj) : insertdiagonal!(uj, v)
            end
            pop!(stack)
            push!(stack, hedges[j-1])
            push!(stack, uj)
        else
            pop!(stack)
            while true
                v, vt = pop!(stack)
                helper = vt ? prev(v) : next(v)
                ori = orientation(origin(uj), orientation(v), orientation(helper))
                if vt && ori == ClockWise || !vt && ori == CounterClockWise
                    break
                end
                insertdiagonal!()
            end
        end
    end
end
