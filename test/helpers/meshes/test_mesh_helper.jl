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
    using BiochemicalVisualization: create_circle_in_local_frame!, create_ellipse_in_local_frame!, create_rectangle_in_local_frame!
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
    c0 = Matrix{Float64}(undef, 3, 20)
    nc0 = Matrix{Float64}(undef, 3, 20)
    create_circle_in_local_frame!(c0,nc0, points2[:, 1], points2[:, 2], points2[:, 3], 20, 1)
    distances = norm.(map(p->p.-points2[:,1], eachcol(c0)))
    @test allapproxequal(distances)

    c1 = Matrix{Float64}(undef, 3, 20)
    nc1 = Matrix{Float64}(undef, 3, 20)
    create_circle_in_local_frame!(c1,nc1, points[:, 1], points[:, 2], points[:, 3], 20, 1)
    distances = norm.(map(p->p.-points[:,1], eachcol(c1)))
    @test allapproxequal(distances)

    # circles and ellipses can be equivalent when applying the right scales
    e1 = Matrix{Float64}(undef, 3, 20)
    ne1 = Matrix{Float64}(undef, 3, 20)
    create_ellipse_in_local_frame!(e1,ne1, points[:, 1], points[:, 2], points[:, 3], 20, 1, 1)
    @test c1==e1

    c2 = Matrix{Float64}(undef, 3, 20)
    nc2 = Matrix{Float64}(undef, 3, 20)
    e2 = Matrix{Float64}(undef, 3, 20)
    ne2 = Matrix{Float64}(undef, 3, 20)
    create_circle_in_local_frame!(c2,nc2, points[:, 1], points[:, 2], 2*points[:, 3], 20, 1)
    create_ellipse_in_local_frame!(e2,ne2, points[:, 1], points[:, 2], points[:, 3], 20, 1, 2)
    @test c2==e2

    # rectangle
    r1 = Matrix{Float64}(undef, 3, 20)
    nr1 = Matrix{Float64}(undef, 3, 20)
    create_rectangle_in_local_frame!(r1, nr1, points[:, 1], points[:, 2], 2*points[:, 3], 20, 1, 2)
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
    c3 = Matrix{Float64}(undef, 3, 20)
    nc3 = Matrix{Float64}(undef, 3, 20)
    e3 = Matrix{Float64}(undef, 3, 20)
    ne3 = Matrix{Float64}(undef, 3, 20)
    r3 = Matrix{Float64}(undef, 3, 20)
    nr3 = Matrix{Float64}(undef, 3, 20)
    create_circle_in_local_frame!(c3,nc3, points[:, 1], points[:, 2], points[:, 3], 20, 0)
    create_ellipse_in_local_frame!(e3,ne3, points[:, 1], points[:, 2], points[:, 3], 20, 0, 0)
    create_rectangle_in_local_frame!(r3,nr3, points[:, 1], points[:, 2], points[:, 3], 20, 0, 0)
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

@testitem "add_faces_to_tube_mesh!" begin
    using BiochemicalVisualization: PlainMesh, add_faces_to_tube_mesh!

    # nothing happens to an empty mesh
    empty = PlainMesh(Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Vector{NTuple{3, Int}}())
    m1 = PlainMesh(Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Matrix{Int}(undef, 3, 0), Vector{NTuple{3, Int}}())
    add_faces_to_tube_mesh!(m1, 3, 0)
    @test m1==empty

    # only one cross-section + endcaps
    small = PlainMesh([0 0 0 -1 1
                        1 -1 -1 0 0
                        0 -1 1 0 0], 
                    copy(repeat([1 0 0], 5)'), 
                    [1 2 3 1 2 3
                    3 1 2 3 1 2
                    4 4 4 5 5 5], 
                    repeat([(0, 0, 0)], 5))
    m2 = PlainMesh([0 0 0 -1 1
                    1 -1 -1 0 0
                    0 -1 1 0 0], 
                copy(repeat([1 0 0], 5)'), 
                Matrix{Int}(undef, 3, 0), 
                repeat([(0, 0, 0)], 5))
    add_faces_to_tube_mesh!(m2, 3, 1)
    println("Faces", m2.connections)
    @test m2==small
    
    # one tube segment
    medium = PlainMesh([0 0 0 1 1 1 -1 1
                        1 -1 -1 1 -1 -1 0 0
                        0 -1 1 0 -1 1 0 0], 
                    copy(repeat([1 0 0], 8)'), 
                    [4 5 6 4 5 6 1 2 3 4 5 6
                    1 2 3 3 1 2 3 1 2 6 4 5
                    3 1 2 6 4 5 7 7 7 8 8 8], 
                    repeat([(0, 0, 0)], 8))
    m3 = PlainMesh([0 0 0 1 1 1 -1 1
                    1 -1 -1 1 -1 -1 0 0
                    0 -1 1 0 -1 1 0 0], 
                copy(repeat([1 0 0], 8)'), 
                Matrix{Int}(undef, 3, 0), 
                repeat([(0, 0, 0)], 8))
    add_faces_to_tube_mesh!(m3, 3, 2)
    println("Faces", m3.connections)
    @test m3==medium
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
