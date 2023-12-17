@testitem "PlainNonStdMesh" begin
    using BiochemicalVisualization: PlainNonStdMesh

    meshes = []
    for i=1:5
        push!(meshes, PlainNonStdMesh([0 1 0;0 0 1; 0 0 0 ], [0 0 0; 0 0 0; 1 1 1], [[1,2,3]], repeat([(255, 0, 0)], 3)))
    end
    @test meshes[1]==meshes[2]

    meshes[2].vertices = [0 1 1;0 0 1; 0 0 0 ]
    @test meshes[1]!=meshes[2]

    meshes[3].normals = [0 0 0; 0 0 0; -1 -1 -1]
    @test meshes[1]!=meshes[3]

    meshes[4].connections = [[3,2,1]]
    @test meshes[1]!=meshes[4]

    meshes[5].colors = [(255, 0, 0), (0, 255, 0), (0, 255, 0)]
    @test meshes[1]!=meshes[5]
end