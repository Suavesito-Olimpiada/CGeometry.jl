module CGeometry

import Base.Iterators
using StaticArrays
using LinearAlgebra

export Orientation, CounterClockWise, ClockWise, NoOrientation, Colinear

struct Orientation
    type::Int8
    function Orientation(n::T ) where {T<:Integer}
        if n ∉ -1:2
            throw(ArgumentError("Argument must be an element of the set {-1,0,1,2}."))
        end
    end
end
const CounterClockWise = Orientation(-1);
const Colinear = Orientation(0);
const ClockWise = Orientation(1);
const NoOrientation = Orientation(2);

include("basics.jl")
include("utils.jl")
include("geometry.jl")

end # module CGeometry
