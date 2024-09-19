include("./dag.jl")
include("./map.jl")

struct TrapezoidalMap{T}
    dcel::DCEL{T}
    # low-high
    bb::@NamedTuple{X::Tuple{T,T}, Y::Tuple{T,T}}
    dag::DAG
    map::Map
end

function TrapezoidalMap(dcel::DCEL{T}; tol) where {T}
    # empty search and map structures
    _dag = DAG() # index, type, left, right
    dag = _dag.data
    _map = Map() # face, leftp, rightp, bottom, top, neightbours
    map = _map.data
    bb = boundingbox(dcel)
    # add the bounding box to both structures to create
    # the trapezoidal map
    push!(dag, Node(1, TTrapezoid, 0, 0))
    push!(map, Trapezoid(1, 0, 0, 0, 0, 0, @SVector[-1, -1, -1, -1]))
    tm = TrapezoidalMap{T}(deepcopy(dcel), (X=bb[1], Y=bb[2]), _dag, _map)
    const BI = Base.Iterators
    const DI = DCELs.Iterators
    # randomize selection of segments
    halfs = dcel |> DI.halfedges |> BI.enumerate |> v->BI.filter(iv->isodd(iv[1]),v) |> v->BI.map(iv->iv[2],v) |> shuffle!∘collect
    # enters main loop
    for half in halfs
        seg = segment(half)
        # check whether the segment is looking right or up,
        # otherwise reverse its direction so that `seg[1]` is the
        # bottom-left point of `seg`
        seg, half = (seg[1].x, seg[1].y) < (seg[2].x, seg[2].y) ? (seg, half) : (reverse(seg), twin(half))
        # compute the trapezoids that the segment goes through
        # it is a vector of `Node.index` of `Node.type=TTrapezoid`
        traps = followsegment(tm, seg; tol)
        n = length(traps)
        nn = length(dag)
        nt = length(map)
        if n == 0
            # TODO: this is an error
            continue
        elseif n == 1 # ─│┌┐└┘├┤┬┴┼╱╲╭╮╰╯
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
            push!(dag, s); nn+=1
            # then we create the new nodes for the DAG and the Map (first C and D)
            C = Node(nt + 1, TTrapezoid, 0, 0)
            push!(dag, C); nn+=1
            Ct = Trapezoid(nn, face(dcel[Halfedge, oldt.top]).i, vertex(half).i, vertex(twin(half)).i, twin(half).i, oldt.top,
                @SVector[-1, oldt.neightbours[2] == -1 ? -1 : n + 1 #=A=#, -1, oldt.neightbours[2] == -1 ? -1 : n + 4 #=B=#])
            push!(map, Ct); nt+=1
            D = Node(nt + 2, TTrapezoid, 0, 0)
            push!(dag, D); nn+=1
            Dt = Trapezoid(nn, face(dcel[Halfedge, oldt.bottom]).i, vertex(half).i, vertex(twin(half)).i, oldt.bottom, half.i,
                @SVector[oldt.neightbours[1] == -1 ? -1 : n + 1 #=A=#, -1, oldt.neightbours[3] == -1 ? -1 : n + 4 #=B=#, -1])
            push!(map, Ct); nt+=1
            # if one of right neightbours exist (top-right, bottom-right)
            # the DAG node and the trapezoid for B
            if oldt.neightbours[3] != -1 || oldt.neightbours[4] != -1
                B = Node(nt + 4, TTrapezoid, 0, 0)
                push!(dag, B); nn+=1
                Bt = Trapezoid(n + 4, face(dcel[Halfedge, oldt.top]).i, oldt.leftp, oldt.rightp,
                    oldt.bottom, oldt.top, @SVector[n + 3 #=D=#, n + 2 #=C=#, oldt.neightbours[3], oldt.neightbours[4]])
                push!(map, Bt); nt+=1
                q = Node(vertex(twin(half)).i, TPoint, nn + 5, nn + 4)
                push!(dag, q); nn+=1
            end
            # if one of left neightbours exist (top-left, bottom-right) then we create
            # the DAG node and the trapezoid for A
            if oldt.neightbours[1] != -1 || oldt.neightbours[2] != -1
                A = Node(nt + 1, TTrapezoid, 0, 0)
                push!(dag, A); n+=1
                At = Trapezoid(n + 1, face(dcel[Halfedge, oldt.bottom]).i, oldt.leftp, vertex(half).i,
                    oldt.bottom, oldt.top, @SVector[oldt.neightbours[1], oldt.neightbours[2], n + 2 #=C=#, n + 3 #=D=#])
                push!(map, At); nt+=1
                p = Node(vertex(half).i, TPoint, nn + 1, nn + 6)
                push!(dag, p); n+=1
            end
            @setnode!(dag[i] = p)
        else
            i = first(traps)
            # check if the leftmost trapezoid left-point already existed
            for i in @view(traps[begin+1:end-1])
            end
        end
    end
    # return the finished trapezoidal map
    return tm
end
