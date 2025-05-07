include("./dag.jl")
include("./map.jl")

const BI = Base.Iterators
const DI = DCELs.Iterators

struct TrapezoidalMap{T}
    dcel::DCEL{T}
    # low-high
    bb::@NamedTuple{X::Tuple{T, T}, Y::Tuple{T, T}}
    dag::DAG
    map::Map
end

function TrapezoidalMap(dcel::DCEL{T}) where {T}
    # empty search and map structures
    dag = DAG() # index, type, left, right
    map = Map() # face, leftp, rightp, bottom, top, neightbours
    bb = boundingbox(dcel)
    return TrapezoidalMap{T}(deepcopy(dcel), (X = bb[1], Y = bb[2]), dag, map)
end

function solve!(tm::TrapezoidalMap; tol)
    # add the bounding box to both structures to create
    # the trapezoidal map
    push!(tm.dag, Node(1, TTrapezoid))
    push!(tm.map, Trapezoid(1))
    # randomize selection of segments
    halfs = tm.dcel |>
        DI.halfedges |>  # all halfedges, lazy
        BI.enumerate |>  # enumerate them, lazy
        v -> BI.filter(iv -> isodd(iv[1]), v) |>  # filter only odd ones
        v -> BI.map(iv -> iv[2], v) |>  # take only the values
        shuffle! ∘ collect  # finally, collect and shuffle
    # enters main loop
    for half in halfs
        # check whether the segment is looking right or up,
        # otherwise reverse its direction so that `seg[1]` is the
        # left-bottom point of `seg`
        #
        #  • <──╮        · ·      #
        #   ╲   │       ╱  │      #    >  ╭──────
        #    ╲ this    ╱   │      #       • <─ p
        #     ╲  │    ╱    │      #  ─────╯  <
        #      · ├─> • ╭─> •      #
        #        ╰─────╯          # Order defined for the plane
        #
        half = (half[1].x, half[1].y) < (half[2].x, half[2].y) ? half : twin(half)
        step!(tm, half; tol)
    end
    # return the finished trapezoidal map
    return tm
end

function step!(tm::TrapezoidalMap, half; tol)
    dag = tm.dag
    map = tm.map
    # compute the trapezoids that the segment goes through
    # it is a vector of `Node.index` of `Node.type=TTrapezoid`
    return traps = followsegment(tm, seg; tol)
end

find(tm::TrapezoidalMap, p::Point{2}; tol = 0) = find(tm, DirSegment(p, p); tol, _build = false)

# _build flag is internally used during building time
function find(tm::TrapezoidalMap, ds::DirSegment{2}; tol = 0, _build = false)
    p = ds[1]
    bb = tm.bb
    dcel = tm.dcel
    # check outside of bounding box
    if !(bb.X[1] ≤ p.x ≤ bb.X[2]) || !(bb.Y[1] ≤ p.y ≤ bb.Y[2])
        return 0, TExterior
    end
    dag = tm.dag
    root = first(dag)
    node = root
    # main loop
    while true
        # x-node, check left-right
        if node.type == TPoint
            v = point(dcel[Vertex, node.index])
            # if p is right or above v, then (*1) - (*2) = true  - false =  1
            #    p is the same as    v, then (*1) - (*2) = false - false =  0
            #    p is left or bellow v, then (*1) - (*2) = false - true  = -1
            #      +-------(1)---------+     +-------(2)---------+
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

function followsegment(tm::TrapezoidalMap, seg; tol = 0)
    bb = tm.bb
    dcel = tm.dcel
    map = tm.map
    # find segment's left point in the trapezoidal map
    # we need to use `_build = true` so that it returns
    # a trapezoid
    trapezoids = Int[]
    index, _ = find(tm, seg[1]; tol, _build = true)
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

function TrapezoidalMap(dcel::DCEL{T}; tol) where {T}
    # enters main loop
    for half in halfs
        n = length(traps)

        if n == 0
            error("Not trapezoids found: broken invariant during execution.")
        end

        ndag = generatenewdag(tm, traps)

        nn = length(dag)
        nt = length(map)

        if n == 1 # ─│┌┐└┘├┤┬┴┼╱╲╭╮╰╯
            # Δ   D(Sᵢ₋₁)   D(Sᵢ₋₁)                     ─┬────────┬─
            # ⇓      │    ⇒    │                       ──┤ oldt Δ │
            # s ┐    Δ         ╰────────────────────╮   ─┴────────┴─
            # C │ 🠅               ╭─────l─────╮  ╭l╮│         ⇓
            # D │ │              [s, C, D, B,╭q, A,╰p]  ─┬──┬───┬──┬─
            # B?├── dag           ├l─╯  │  ╰r╯╰──r──╯    │A │  C│  │
            # q?│ │               ╰──r──╯              ──┤  p─s─q  │
            # A?┘ │                                      │  │D  │ B│
            # p?<─┘                                     ─┴──┴───┴──┴─
            i = first(traps)
            # we get the first trapezoid Δ
            oldt = map[i]
            # we create the new segment that will divide Δ
            s = Node(half.i, TSegment, nn + 2, nn + 3)
            push!(dag, s); nn += 1
            # then we create the new nodes for the DAG and the Map (first C and D)
            C = Node(nt + 1, TTrapezoid, 0, 0)
            push!(dag, C); nn += 1
            Ct = Trapezoid(
                nn, face(dcel[Halfedge, oldt.top]).i, vertex(half).i, vertex(twin(half)).i, twin(half).i, oldt.top,
                @SVector[-1, oldt.neightbours[2] == -1 ? -1 : n + 1 #=A=#, -1, oldt.neightbours[2] == -1 ? -1 : n + 4 #=B=#]
            )
            push!(map, Ct); nt += 1
            D = Node(nt + 2, TTrapezoid, 0, 0)
            push!(dag, D); nn += 1
            Dt = Trapezoid(
                nn, face(dcel[Halfedge, oldt.bottom]).i, vertex(half).i, vertex(twin(half)).i, oldt.bottom, half.i,
                @SVector[oldt.neightbours[1] == -1 ? -1 : n + 1 #=A=#, -1, oldt.neightbours[3] == -1 ? -1 : n + 4 #=B=#, -1]
            )
            push!(map, Ct); nt += 1
            # if one of right neightbours exist (top-right, bottom-right)
            # the DAG node and the trapezoid for B
            if oldt.neightbours[3] != -1 || oldt.neightbours[4] != -1
                B = Node(nt + 4, TTrapezoid, 0, 0)
                push!(dag, B); nn += 1
                Bt = Trapezoid(
                    n + 4, face(dcel[Halfedge, oldt.top]).i, oldt.leftp, oldt.rightp,
                    oldt.bottom, oldt.top, @SVector[n + 3 #=D=#, n + 2 #=C=#, oldt.neightbours[3], oldt.neightbours[4]]
                )
                push!(map, Bt); nt += 1
                q = Node(vertex(twin(half)).i, TPoint, nn + 5, nn + 4)
                push!(dag, q); nn += 1
            end
            # if one of left neightbours exist (top-left, bottom-right) then we create
            # the DAG node and the trapezoid for A
            if oldt.neightbours[1] != -1 || oldt.neightbours[2] != -1
                A = Node(nt + 1, TTrapezoid, 0, 0)
                push!(dag, A); n += 1
                At = Trapezoid(
                    n + 1, face(dcel[Halfedge, oldt.bottom]).i, oldt.leftp, vertex(half).i,
                    oldt.bottom, oldt.top, @SVector[oldt.neightbours[1], oldt.neightbours[2], n + 2 #=C=#, n + 3 #=D=#]
                )
                push!(map, At); nt += 1
                p = Node(vertex(half).i, TPoint, nn + 1, nn + 6)
                push!(dag, p); n += 1
            end
            @setnode!(dag[i] = p)
        else
            i = first(traps)
            # check if the leftmost trapezoid left-point already existed
            for i in @view(traps[(begin + 1):(end - 1)])
            end
        end
    end
    # return the finished trapezoidal map
    return tm
end
