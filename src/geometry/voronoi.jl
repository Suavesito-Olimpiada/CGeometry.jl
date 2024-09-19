using Random

export incrementalvoronoi


function findface(geometry::DCEL{T}, p::Point{2,T}) where {T}
    for cell in Iterators.faces(geometry)
        if p in cell
            return cell.i
        end
    end
    return nothing
end

function splitcell()
end

function incrementalvoronoi(P::AbstractVector{Point{2,T}}) where {T}
    xₘᵢₙ, xₘₐₓ = 1.05 .* extrema(p -> p.x, P)
    yₘᵢₙ, yₘₐₓ = 1.05 .* extrema(p -> p.y, P)
    shuffle!(P)
    v = (Point(xₘᵢₙ, yₘₐₓ), Point(xₘᵢₙ, yₘᵢₙ), Point(xₘₐₓ, yₘᵢₙ), Point(xₘₐₓ, yₘₐₓ))
    geometry = DCEL{T}([DirSegment(v[1], v[2]), DirSegment(v[2], v[3]), DirSegment(v[3], v[4]), DirSegment(v[4], v[1])])
    updatefaces!(geometry)
    voronoi = Dict{HalfedgeHandle{T}, Point{2,T}}(geometry[1] => first(P))
    for point in @view(P[begin+1:end])
        cell = findface(geometry, point)
        vorpoint = voronoi[cell]
        mdtx = mediatrix(point, vorpoint)
        for facet in Iterators.halfedges(cell)
            for hedge in facet
                s = segment(hedge)
                isintersect(mdtx, s)
                ip = intersection(mdtx, s)
                if ip ∉ (s[1], s[2])
                    cut(s, ip)
                end
            end
        end
    end
end
