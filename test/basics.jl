@testset verbose=true "$T" for (T, TF64) in ((Point, PointF64), (Vec, VecF64))
    p = T(1, 2, 3)
    pf = TF64(p)

    @test length(p) == 3
    @test p[1] == 1
    @test_throws BoundsError p[4]
    @test p.x == 1
    @test p.y == 2
    @test p.z == 3
end

@testset verbose=true "Point ↔ Vec" begin
    p₁ = Point(1, 0)
    p₂ = Point(2, 0)
    v = Vec(1, 0)

    @test v == Vec(p₁)
    @test p₁ == Point(v)

    @test_throws MethodError p₁ + p₂
    @test p₁ + v == p₂
    @test v + v == Vec(2, 0)
    @test p₂ - p₁ == v
    @test v - v == Vec(0, 0)
end

@testset verbose=true "$T" for T in (Segment, DirSegment)
    p₁ = Point(1, 0)
    p₂ = Point(0, 1)
    s = T(p₁, p₂)

    @test s.p₁ == p₁
    @test s.p₂ == p₂
    @test s[1] == s.p₁
    @test s[2] == s.p₂
end

@testset verbose=true "Segment ↔ DirSegment" begin
    p₁ = Point(1, 0)
    p₂ = Point(0, 1)
    s = Segment(p₁, p₂)
    ds = DirSegment(p₁, p₂)

    @test DirSegment(s) == ds
    @test Segment(ds) == s

    @test reverse(ds) == DirSegment(p₂, p₁)
    @test reverse(s) == s
end

@testset verbose=true "Cross product" begin
    v₁ = Vec(1, 0)
    v₂ = Vec(0, 1)

    @test v₁ × v₂ == 1
    @test v₂ × v₁ == -1
    @test v₁ × v₁ == 0
end

@testset verbose=true "O(1) operations" begin
    @testset verbose=true "Orientation angle" begin
        p₀ = Point(0, 0)
        p₁ = Point(1, 2)
        p₂ = Point(2, 1)

        @test orientation(p₀, p₁, p₂) == ClockWise
        @test orientation(p₀, p₂, p₁) == CounterClockWise
    end

    @testset verbose=true "Orientation segment-point" begin
        s = DirSegment(Point(0, 0), Point(1, 2))
        p = Point(2, 1)

        @test orientation(s, p) == ClockWise
    end

    @testset verbose=true "Orientation adjacent segments" begin
        ds₁ = DirSegment(Point(0, 0), Point(1, 2))
        ds₂ = DirSegment(Point(0, 0), Point(2, 1))
        s₁ = Segment(Point(0, 0), Point(1, 2))
        s₂ = Segment(Point(0, 0), Point(2, 1))

        @test orientation(ds₁, ds₂) == ClockWise
        @test orientation(ds₂, ds₁) == CounterClockWise
        @test orientation(s₁, s₂) == ClockWise
        @test orientation(s₂, s₁) == CounterClockWise
    end

    @testset verbose=true "Intersection of segments" begin
        ds₁ = DirSegment(Point(0, 0), Point(2, 2))
        ds₂ = DirSegment(Point(0, 2), Point(2, 0))
        s₁ = Segment(Point(0, 0), Point(0, 1))
        s₂ = Segment(Point(1, 1), Point(1, 0))


        @test isintersect(ds₁, ds₂) == true
        @test isintersect(s₁, s₂) == false
        @test isintersect(ds₁, s₂) == true
    end
end
