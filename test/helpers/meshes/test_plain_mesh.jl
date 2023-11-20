@testitem "PlainMesh" begin
    using BiochemicalVisualization: PlainMesh
    import Meshes

    points = [0 1 0
    0 0 1
    0 0 0]

    simplemesh = Meshes.SimpleMesh(map(p->Tuple(p), eachcol(points)), [Meshes.connect((1, 2, 3))])
    plainmesh = PlainMesh(simplemesh, (0, 127, 0))

    @test plainmesh.vertices == points
    @test plainmesh.connections == reshape([1 2 3], (3, 1))
    @test plainmesh.colors == [(0, 127, 0), (0, 127, 0), (0, 127, 0)]


    # faces that are not triangles are not allowed
    points2 = [0 1 0 1
    0 0 1 1
    0 0 0 0]

    simplemesh = Meshes.SimpleMesh(map(p->Tuple(p), eachcol(points)), [Meshes.connect((1, 2, 3, 4))])
    @test_throws ErrorException PlainMesh(simplemesh, (0, 127, 0))

end