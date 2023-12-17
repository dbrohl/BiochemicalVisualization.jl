
@testitem "filter_points_threshold" begin
    using BiochemicalVisualization: filter_points_threshold
    using LinearAlgebra

    function mapping_to_containing_list(mapping)
        res = []
        for i=eachindex(mapping)
            if(mapping[i]!=-1)
                push!(res, i)
            end
        end
        return res
    end
    
    constants = repeat([1.0,0,0], 1, 6)
    constants2 = repeat([0.0,1,0], 1, 6)
    fixed = [1, 4, 5]

    # second element really shows the number of indices!=-1
    tup = filter_points_threshold(constants, constants2, Vector{Int}())
    @test length(filter(x -> x!=-1, tup[1]))==tup[2]

    tup = filter_points_threshold(constants, constants2, fixed)
    @test length(filter(x -> x!=-1, tup[1]))==tup[2]

    tup = filter_points_threshold(constants, constants2, Vector{Int}(), with_color=true)
    @test length(filter(x -> x!=-1, tup[1]))==tup[2]

    # no change -> only one index remains
    @test filter_points_threshold(constants, constants2, Vector{Int}())[2]==1

    # fixed_indices are contained in the result
    mapped_indices = filter_points_threshold(constants, constants2, fixed)[1]
    res = mapping_to_containing_list(mapped_indices)

    for a in fixed
        @test a ∈ res
    end
        
    # changing color results in more remaining points
    @test filter_points_threshold(repeat([1.0,0,0], 1, 12), repeat([0.0,1,0], 1, 12), Vector{Int}(), with_color=true)[2]>1

    # changing tangents or normals lead to more remaining points
    up = Matrix{Float64}(undef, 3, 6)
    for (i,a) in enumerate(range(0, 1, 6))
        up[:, i] = [1-a, 0, a]./norm([1-a, 0, a])
    end
    @test filter_points_threshold(up, constants2, Vector{Int}())[2]>1
    @test filter_points_threshold(constants2, up, Vector{Int}())[2]>1

end