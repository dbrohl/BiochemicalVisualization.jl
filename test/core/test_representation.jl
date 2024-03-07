@testitem "ismesh + isprimitivecollection" begin
    using BiochemicalVisualization
    using BiochemicalVisualization: PlainMesh, ismesh, isprimitivecollection
    import GeometryBasics

    # recognises correct meshes
    triangle = PlainMesh([0.0 0 1
    0 1 0
    0 0 0], [0.0 0 0
    0 0 0
    1 1 1],reshape([1, 2, 3], (3, 1)), [(0, 0, 0), (255, 0, 0), (0, 255, 0)])

    tri_rep = Representation(triangle)
    @test ismesh(tri_rep) && !isprimitivecollection(tri_rep)

    # recognises correct primitive collections
    spheres = [GeometryBasics.Sphere(GeometryBasics.Point{3}([0.0, 0, 0]), 1.0)]
    sphere_colors = ["#FF00FF"]

    sph_rep = Representation{Float64}(
        primitives=Dict([("spheres", spheres)]), 
        colors=Dict([("spheres", sphere_colors)]))

    @test !ismesh(sph_rep) && isprimitivecollection(sph_rep)
    


    # empty representations are neither
    empty = Representation{Float64}()
    @test !ismesh(empty) && !isprimitivecollection(empty)

    # mixed representations are neither
    mixed = Representation{Float64}(
        vertices = vec([0.0 0 1
        0 1 0
        0 0 0]),
        primitives=Dict([("spheres", spheres)]), 
        colors=Dict([("spheres", sphere_colors)]))
    @test !ismesh(mixed) && !isprimitivecollection(mixed)
end