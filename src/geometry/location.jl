include("../basics.jl")
include("../geometry.jl")

using Base: @NamedTuple
using Random: shuffle!

using StructArrays
using StaticArrays

using .DCELs

export TrapezoidalMap



find(tm::TrapezoidalMap, p::Point{2}; tol=0) = find(tm, DirSegment(p, p); tol, _build=false)

# _build flag is internally used during building time
function find(tm::TrapezoidalMap, ds::DirSegment{2}; tol=0, _build=false)
    p = ds[1]
    bb = tm.bb
    dcel = tm.dcel
    # check outside of bounding box
    if !(bb.X[1] ≤ p.x ≤ bb.X[2]) || !(bb.Y[1] ≤ p.y ≤ bb.Y[2])
        return 0, TExterior
    end
    dag = tm.dag.data
    root = first(dag)
    node = root
    # main loop
    while true
        # x-node, check left-right
        if node.type == TPoint
            v = point(dcel[Vertex, node.index])
            # if p is right or above v, then (1) true,  (2) false =  1
            #    p is the same as    v, then (1) false, (2) false =  0
            #    p is left or bellow v, then (1) false, (2) true  = -1
            # +-------(1)---------+   +---------(2)---------+
            ord = (p.x, p.y) > (v.x, v.y) - (p.x, p.y) < (v.x, v.y)
            if ord == 1 || _build
                node = dag[node.right]
            elseif ord == -1
                node = dag[node.left]
            elseif ord == 0
                break
            end
        elseif node.type == TSegment
            s = segment(dcel[Halfedge, node.index])
            # check if segment is pointing right or up,
            # otherwise reverse its direction
            s = (s[1].x, s[1].y) < (s[2].x, s[2].y) ? s : reverse(s)
            # then the point _must_ be at its right (clockwise)
            # to be under the segment
            ori = orientation(s, p; tol)
            if ori == ClockWise            # under
                node = dag[node.right]
            elseif ori == CounterClockWise # over
                node = dag[node.left]
            elseif !_build                 # on the line
                break
            else                           # on the line and build time
                # when building the trapezoidal map a point can be on `s` if
                # it is the extreme left of it, we _have_ to choose
                # a trapezoid, so we check whether the other extreme of the segment
                # being inserted is above or bellow `s`
                if orientation(s, ds[2]; tol) == CounterClockWise
                    node = dag[node.left]  # over
                else
                    node = dag[node.right] # under
                end
            end
        elseif node.type == TTrapezoid
            break
        end
    end
    return node.index, node.type
end

function followsegment(tm, seg; tol=0)
    bb = tm.bb
    dcel = tm.dcel
    map = tm.map.data
    # find segment's left point in the trapezoidal map
    # we need to use `_build = true` so that it returns
    # a trapezoid
    trapezoids = Int[]
    index, _ = find(tm, seg[1]; tol, _build=true)
    push!(trapezoids, nindex)
    # continue until `seg` ends
    while true
        # get right point of the current trapezoid
        rightp = if index != 0
            point(dcel[Vertex, map[index].rightp])
        else
            Point(bb.X[2], bb.Y[2])
        end
        # if the right point of the new segment is
        # to the left of the rightmost point of the
        # current trapezoid then the line has finished
        if seg[2].x < rightp.x
            break
        end
        # we know that `seg` is looking to the right (or up),
        # therefore, for rightp to be above `seg` it must be
        # counterclockwise from it
        ori = orientation(seg, rightp; tol)
        if ori == CounterClockWise  # right-bottom
            index = map[index].neightbours[3]
            nindex = map[index].node
        elseif ori == ClockWise     # right-top
            index = map[index].neightbours[4]
            nindex = map[index].node
        else
            error("Unexpected data.")
        end
        push!(trapezoids, nindex)
    end
    return trapezoids
end
