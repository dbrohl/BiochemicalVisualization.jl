@testitem "frame_construction" begin
    using BiochemicalVisualization: rmf, frames_from_two_splines
    using LinearAlgebra

    points = [0.0 1 2 3 4 5
    0 0 0 0 0 0
    0 0 0 0 0 0]

    tangents = Matrix{Float64}(undef, 3, 6)
    for i=2:6
        tangents[:, i-1] = points[:, i]-points[:, i-1]
    end
    tangents[:, end] = [1, 0, 0]

    t, n, b = rmf(points, tangents)
    for i=axes(t, 2)
        # test for unit vectors
        @test norm(t[:, i])≈1
        @test norm(n[:, i])≈1
        @test norm(b[:, i])≈1

        # orthogonal frame
        @test dot(t[:, i], n[:, i]) ≈ 0
        @test dot(t[:, i], b[:, i]) ≈ 0
        @test dot(n[:, i], b[:, i]) ≈ 0
    end

    minor_points = [0.0 1 2 3 4 5
    1 1 1 1 1 1
    0 0 0 0 0 0]

    t, n, b = frames_from_two_splines(points, tangents, minor_points)
    println("n$n, b$b, t$t")
    for i=axes(t, 2)
        # test for unit vectors
        @test norm(t[:, i])≈1
        @test norm(n[:, i])≈1
        @test norm(b[:, i])≈1

        # orthogonal frame
        @test dot(t[:, i], n[:, i]) ≈ 0
        @test dot(t[:, i], b[:, i]) ≈ 0
        @test dot(n[:, i], b[:, i]) ≈ 0
    end
end