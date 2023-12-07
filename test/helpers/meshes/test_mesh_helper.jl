@testitem "color!" begin
    using BiochemicalVisualization: PlainMesh, PlainNonStdMesh, color!

    points = [0 0 1
    0 1 0
    0 0 0]

    normals = [0 0 0
    0 0 0
    1 1 1]

    colors = [(0, 0, 0), (255, 0, 0), (0, 255, 0)]
    white = [(255, 255, 255), (255, 255, 255), (255, 255, 255)]
    
    plain_mesh = PlainMesh(points, normals, reshape([1, 2, 3], (3, 1)), copy(colors))
    @test plain_mesh.colors!=white
    color!(plain_mesh, (255, 255, 255))
    @test plain_mesh.colors==white

    println("A", colors)

    nonstd_plain_mesh = PlainNonStdMesh(points, normals, [[1,2,3]], copy(colors))
    @test nonstd_plain_mesh.colors!=white
    color!(nonstd_plain_mesh, (255, 255, 255))
    @test nonstd_plain_mesh.colors==white

    @test_throws MethodError color!(plain_mesh, 255)

end

@testitem "create_in_local_frame" begin
    using BiochemicalVisualization: create_circle_in_local_frame, create_ellipse_in_local_frame, create_rectangle_in_local_frame
    using LinearAlgebra


    function allapproxequal(array)
        for item in array
            if(!(item ≈ array[1]))
                return false
            end
        end
        return true
    end

    points = [0.0 1 0
    0 0 1
    0 0 0]

    points2 = [1.0 -2 0
    1 0 0 
    1 0 2]

    # all points of a circle have the same distance to the origin
    c0,nc0 = create_circle_in_local_frame(points2[:, 1], points2[:, 2], points2[:, 3], 20, 1)
    distances = norm.(map(p->p.-points2[:,1], eachcol(c0)))
    @test allapproxequal(distances)

    c1,nc1 = create_circle_in_local_frame(points[:, 1], points[:, 2], points[:, 3], 20, 1)
    distances = norm.(map(p->p.-points[:,1], eachcol(c1)))
    @test allapproxequal(distances)

    # circles and ellipses can be equivalent when applying the right scales
    e1,ne1 = create_ellipse_in_local_frame(points[:, 1], points[:, 2], points[:, 3], 20, 1, 1)
    @test c1==e1

    c2,nc2 = create_circle_in_local_frame(points[:, 1], points[:, 2], 2*points[:, 3], 20, 1)
    e2,ne2 = create_ellipse_in_local_frame(points[:, 1], points[:, 2], points[:, 3], 20, 1, 2)
    @test c2==e2

    # rectangle
    r1, nr1 = create_rectangle_in_local_frame(points[:, 1], points[:, 2], 2*points[:, 3], 20, 1, 2)
    min_max_x = extrema(r1[1, :])
    min_max_y = extrema(r1[2, :])
    for point in eachcol(r1)
        @test point[1] ∈ min_max_x || point[2] ∈ min_max_y
    end

    # when using a frame without z component, all points lie in the xy-plain_mesh
    @test all(c1[3, :].==0)
    @test all(e1[3, :].==0)
    @test all(r1[3, :].==0)

    # resolution parameter is respected
    @test size(c1)==(3, 20)
    @test size(c1)==size(e1)
    @test size(c1)==size(r1)

    # zero radius => all points are the center point
    c3,nc3 = create_circle_in_local_frame(points[:, 1], points[:, 2], points[:, 3], 20, 0)
    e3,ne3 = create_ellipse_in_local_frame(points[:, 1], points[:, 2], points[:, 3], 20, 0, 0)
    r3,nr3 = create_rectangle_in_local_frame(points[:, 1], points[:, 2], points[:, 3], 20, 0, 0)
    for obj in [c3, e3, r3]
        for point in eachcol(obj)
            @test point ≈ points[:, 1]
        end
    end
end

@testitem "merge_multiple_meshes" begin
    using BiochemicalVisualization: PlainMesh, merge_multiple_meshes

    triangle = PlainMesh([0 0 1
    0 1 0
    0 0 0], [0 0 0
    0 0 0
    1 1 1],reshape([1, 2, 3], (3, 1)), [(0, 0, 0), (255, 0, 0), (0, 255, 0)])

    square = PlainMesh([0 0 1 1
    0 1 0 1
    0 0 0 0], [0 0 0 0
    0 0 0 0
    1 1 1 1],[1 2
    2 3
    3 4], [(0, 0, 0), (255, 0, 0), (0, 255, 0), (0, 0, 255)])

    empty = PlainMesh(Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Vector{NTuple{3, Int}}())

    ts = merge_multiple_meshes([triangle, square])
    @test size(ts.vertices) == (3, 7)
    @test size(ts.connections) == (3, 3)
    @test size(ts.colors) == (7, )

    st = merge_multiple_meshes([square, triangle])
    @test Set(eachcol(ts.vertices)) == Set(eachcol(st.vertices))

    @test merge_multiple_meshes([triangle, empty]) == triangle
    @test merge_multiple_meshes([triangle]) == triangle

    @test merge_multiple_meshes(Vector{PlainMesh{Int64}}()) == empty
end

@testitem "connect_circles_to_tube" begin
    using BiochemicalVisualization: PlainMesh, PlainNonStdMesh, connect_circles_to_tube, merge_multiple_meshes

    points = [0 0 1
    0 1 0
    0 0 0]
    normals = [0 0 0
    0 0 0
    1 1 1]
    colors = [(0, 0, 0), (255, 0, 0), (0, 255, 0)]
    a = PlainNonStdMesh(points, normals, [[1,2,3]], colors)

    a_std = PlainMesh(points, normals, reshape([1,2,3], (3, 1)), colors)

    empty = PlainMesh(Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Vector{NTuple{3, Int}}())

    @test connect_circles_to_tube(Vector{PlainNonStdMesh{Int}}())==empty
    res_single = connect_circles_to_tube([a])
    @test res_single.vertices==a.vertices
    @test res_single.colors==a.colors
    @test size(res_single.connections)==(3, 0)
    
    res_connect = connect_circles_to_tube([a, a])
    @test size(res_connect.vertices, 2) == 2*size(a.vertices, 2)
    @test length(res_connect.colors) == 2*length(a.colors)
    @test size(res_connect.connections, 2) == 6

    res_endpoints = connect_circles_to_tube([a, a], (([0, 0, 0], [0,0,1]), ([1,1,1], [0,0,1])))
    @test size(res_endpoints.vertices, 2)==size(res_connect.vertices, 2)+2
    @test length(res_endpoints.colors) == length(res_connect.colors)+2
    @test size(res_endpoints.connections, 2) > size(res_connect.connections, 2)
end

@testitem "DebugOutput" begin
    using BiochemicalVisualization: ColoredMesh, local_frame_mesh, local_arrow_mesh

    points = [0.0 1 0 0
    0 0 1 0
    0 0 0 1]
    res = local_frame_mesh(points[:, 1], points[:, 2], points[:, 3], points[:, 4])
    @test typeof(res)<:ColoredMesh
    @test size(res.vertices, 2)!=0

    res = local_arrow_mesh(points[:, 1], points[:, 2], (255, 0, 0))
    @test typeof(res)<:ColoredMesh
    @test size(res.vertices, 2)!=0
end
