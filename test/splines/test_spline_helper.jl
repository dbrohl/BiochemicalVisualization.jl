@testitem "get_c_alpha_positions" begin
    using BiochemicalAlgorithms
    using BiochemicalVisualization: get_c_alpha_positions

    path = normpath(joinpath(@__DIR__, "..", "data", "2lnf.pdb"))
    pdb = load_pdb(path)
    positions, indices = get_c_alpha_positions(first(eachchain(pdb)))

    @test size(positions, 2)==13
    @test length(indices)==13

    path = normpath(joinpath(@__DIR__, "..", "data", "AlaAla.pdb"))
    pdb = load_pdb(path, Float64)
    positions, indices = get_c_alpha_positions(first(eachchain(pdb)))

    @test positions == [0 0 0
                        0.605 0.980 3.639]'
    @test length(indices)==2
end