using Test
using CGeometry
using LinearAlgebra

@testset "Basics" verbose=true begin
    include("./basics.jl")
end

@testset "Geometry" verbose=true begin
    include("./geometry.jl")
end
