@testitem "PlainNonStdMesh" begin
    using BiochemicalVisualization: PlainNonStdMesh
    import Meshes

    points = [0 1 0
    0 0 1
    0 0 0]

    simplemesh = Meshes.SimpleMesh(map(p->Tuple(p), eachcol(points)), [Meshes.connect((1, 2, 3))])
    plainmesh = PlainNonStdMesh(simplemesh, (0, 127, 0))

    @test plainmesh.vertices == points
    println(plainmesh.connections, typeof(plainmesh.connections))
    @test plainmesh.connections == [[1,2,3]]
    @test plainmesh.colors == [(0, 127, 0), (0, 127, 0), (0, 127, 0)]

    # faces that are not triangles are allowed
    points2 = [0 1 0 1
    0 0 1 1
    0 0 0 0]

    simplemesh = Meshes.SimpleMesh(map(p->Tuple(p), eachcol(points2)), [Meshes.connect((1, 2, 3, 4))])
    plainmesh = PlainNonStdMesh(simplemesh, (0, 127, 0))

    @test plainmesh.vertices == points2
    @test plainmesh.connections == [[1,2,3,4]]
    @test plainmesh.colors == [(0, 127, 0), (0, 127, 0), (0, 127, 0), (0, 127, 0)]
end