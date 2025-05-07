using Base: @NamedTuple
using Random: shuffle!

using StructArrays
using StaticArrays

using .DCELs

# export TrapezoidalMap
export isinside

function isinside(p::Point{2, T}, face::Face{T}; tol = 0) where {T}
    cruses = 0
    for hedge in face
        d = segment(hedge)
        # looking upwards
        s = d[1].y ≤ d[2].y ? d : reverse(d)
        if s[1].y ≤ p.y ≤ s[2].y
            if orientation(s, p; tol) == CounterClockWise
                cruses += -Int(orientation(d, p; tol))
            end
        end
    end
    return cruses > 0
end

# include("../geometry/location/trapezoidalmap.jl")
