@testset verbose=true "Convex Hull" begin
    points = Point.([[0.0, 0.0], [0.5, 0.0], [1.0, 0.0], [1.0, 0.5], [1.0, 1.0], [0.5, 1.0], [0.0, 1.0], [0.0, 0.5]])
    chull = Point.([[0.0, 1.0], [1.0, 1.0], [1.0, 0.0], [0.0, 0.0]])

    @testset verbose=true "Jarvis March" begin
        @test jarvismarch(points) == chull
    end

    @testset verbose=true "Quick Hull" begin
        @test quickhull(points) == chull
    end
end


