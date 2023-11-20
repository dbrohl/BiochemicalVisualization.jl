@testitem "filter_points_threshold" begin
    using BiochemicalVisualization: filter_points_threshold
    using LinearAlgebra
    
    constants = repeat([1.0,0,0], 1, 6)
    constants2 = repeat([0.0,1,0], 1, 6)


    # no change -> only one index remains
    @test length(filter_points_threshold(constants, constants2, Vector{Int}()))==1

    # fixed_indices are contained in the result
    fixed = [1, 4, 5]
    res = filter_points_threshold(constants, constants2, fixed)
    for a in fixed
        @test a ∈ res
    end
        
    # changing color results in more remaining points
    @test length(filter_points_threshold(constants, constants2, Vector{Int}(), [(255, 0, 0), (255, 255, 0), (0, 255, 0), (0, 255, 255), (0, 0, 255), (255, 0, 255)]))>1

    # changing tangents or normals lead to more remaining points
    up = Matrix{Float64}(undef, 3, 6)
    for (i,a) in enumerate(range(0, 1, 6))
        up[:, i] = [1-a, 0, a]./norm([1-a, 0, a])
    end
    @test length(filter_points_threshold(up, constants2, Vector{Int}()))>1
    @test length(filter_points_threshold(constants2, up, Vector{Int}()))>1



end