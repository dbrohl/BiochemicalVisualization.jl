@testitem "c_alphas_to_points" begin
    using BiochemicalAlgorithms
    using BiochemicalVisualization: c_alphas_to_points

    path = normpath(joinpath(@__DIR__, "..", "data", "2lnf.pdb"))
    pdb = load_pdb(path)
    positions, indices = c_alphas_to_points(first(eachchain(pdb)))

    @test length(positions)==13
    @test length(indices)==13



    path = normpath(joinpath(@__DIR__, "..", "data", "AlaAla.pdb"))
    pdb = load_pdb(path, Float64)
    positions, indices = c_alphas_to_points(first(eachchain(pdb)))

    @test positions == [[0, 0, 0], [0.605,0.980,3.639]]
    @test length(indices)==2



end