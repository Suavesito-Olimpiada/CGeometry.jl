module CGeometry

import Base.Iterators
using StaticArrays
using LinearAlgebra

export Orientation, CounterClockWise, ClockWise, NO, Colinear

# ClockWise  := Clockwise
# CounterClockWise := Counter-clockwise
# Colinear  := Colinear
# NO  := No Orientation
@enum Orientation begin
    CounterClockWise = -1
    Colinear = 0
    ClockWise = 1
    NO = 2
end

include("basics.jl")
include("utils.jl")
include("geometry.jl")

end # module CGeometry
