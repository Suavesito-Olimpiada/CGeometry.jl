module CGeometry

import Base.Iterators
using StaticArrays
using LinearAlgebra

export Orientation, CounterClockWise, ClockWise, NoOrientation, Colinear

struct BrokenInvariant <: Exception
    msg::String
end

@enum Orientation::Int8 CounterClockWise = -1 Colinear = 0 ClockWise = 1 NoOrientation = 2

include("basics.jl")
include("utils.jl")
include("geometry.jl")

end # module CGeometry
