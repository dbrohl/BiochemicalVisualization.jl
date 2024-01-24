@testitem "HashGrid" begin
    using BiochemicalVisualization: HashGrid, push!, world_to_index, each_neighbor, iterate, IteratorSize, eltype, isdone
    using LinearAlgebra
    # test the constructor (has the grid the right dimensions?)
    bounding_box = [-2.0 2.0; -4 4; -6 6]
    grid = HashGrid{Float64, Vector{Float64}}(bounding_box, [1.0, 0.5, 2.0], identity, 20)
    @test size(grid.grid)==(5, 17, 7)
    @test grid.data==[]
    @test grid.origin == [-2.0, -4, -6]
    @test grid.box_size == [1.0, 0.5, 2.0]

    # test world_to_index
    @test world_to_index(grid, [-2.0, -4, -6]) == [1, 1, 1]
    @test world_to_index(grid, [2.0, -4, -6]) == [5, 1, 1]
    @test world_to_index(grid, [0, 0, 0]) == [3, 9, 4]
    @test_throws DomainError world_to_index(grid, [-2.1, 0, 0])
    @test_throws DomainError world_to_index(grid, [0, -5, 0])


    # construct new grid with cube boxes for easier testing
    grid = HashGrid{Float64, Vector{Float64}}(bounding_box, [1.0, 1.0, 1.0], identity, 20)

    # test data insertion
    push!(grid, [-2.0, -4, -6])
    @test length(grid.data)==1
    @test grid.grid[1, 1, 1] == [1]

    push!(grid, [-1.0, -4, -6])
    @test length(grid.data)==2
    @test grid.grid[2, 1, 1] == [2]

    push!(grid, [-1.5, -4, -6])
    @test length(grid.data)==3
    @test grid.grid[1, 1, 1] == [1, 3]

    @test_throws DomainError push!(grid, [-2.5, -4, -6])

    # more data for later tests
    push!(grid, [-0.25, -4, -6])
    push!(grid, [0.0, -4, -6])
    

    # test data access
    #result = Set()
    result = collect(each_neighbor(grid, [-2.0, -4, -6], 1.0))
    @test result == [[-2.0, -4, -6], [-1.5, -4, -6],[-1.0, -4, -6]]

    # more test data (compare results of grid with result of brute force)
    points = rand(Float64, 3, 400)
    points[1, :] .= points[1, :] .* 4 .- 2
    points[2, :] .= points[2, :] .* 8 .- 4
    points[3, :] .= points[3, :] .* 12 .- 6

    reference_point = rand(Float64, 3)
    dist_one = Set()
    dist_zero_five = Set()

    for i in axes(points, 2)
        push!(grid, points[:, i])
        dist = norm(points[:, i] .- reference_point)
        if dist < 1
            push!(dist_one, points[:, i])
        end
        if dist < 0.5 
            push!(dist_zero_five, points[:, i])
        end
    end

    result_one = Set()
    temp = collect(each_neighbor(grid, reference_point, 1.0))
    if(length(temp)>0)
        push!(result_one, temp...)
    end
    result_zero_five = Set()
    temp = collect(each_neighbor(grid, reference_point, 0.5))
    if(length(temp)>0)
        push!(result_zero_five, temp...)
    end

    @test issetequal(result_one, dist_one)
    @test issetequal(result_zero_five, dist_zero_five)


    # test position_access_function for complex data
    struct Container
        data
    end
    grid = HashGrid{Float64, Container}(bounding_box, [1.0, 1.0, 1.0], item -> item.data)
    c1 = Container([-2.0, -4, -6])
    c2 = Container([1.5, 0.1, 1.4])
    push!(grid, c1, c2)
    @test [c1]==collect(each_neighbor(grid, [-2.0, -4, -6], 1.0))
    @test [c2]==collect(each_neighbor(grid, [1.5, 0.1, 1.4], 1.0))

    # test edge cases
    # one dimension has 0 width and box_size larger than width (should still work)
    bounding_box = [-2.0 -2.0; -4 4; -6 6]
    grid = HashGrid{Float64, Vector{Float64}}(bounding_box, [1.0, 1.0, 1.0], identity, 20)
    push!(grid, [-2.0, 0, 0])
    @test 1==1


    # box_size 0 or negative
    @test_throws ArgumentError HashGrid{Float64, Vector{Float64}}(bounding_box, [0.0, 1.0, 1.0], identity, 20)
    @test_throws ArgumentError HashGrid{Float64, Vector{Float64}}(bounding_box, [1.0, -1.0, 1.0], identity, 20)
end