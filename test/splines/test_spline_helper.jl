@testitem "get_c_alpha_positions" begin
    using BiochemicalAlgorithms
    using BiochemicalVisualization: get_c_alpha_positions

    function loadPDB(path)
        fdb = FragmentDB()
        pdb = load_pdb(path)
        normalize_names!(pdb, fdb);
        reconstruct_fragments!(pdb, fdb);
        add_secondary_structures!(pdb, path)
        return pdb
    end

    path = normpath(joinpath(@__DIR__, "..", "data", "2lnf.pdb"))
    pdb = loadPDB(path)
    positions, indices = get_c_alpha_positions(first(eachchain(pdb)))

    @test size(positions, 2)==13
    @test length(indices)==13

    path = normpath(joinpath(@__DIR__, "..", "data", "AlaAla.pdb"))
    pdb = loadPDB(path)
    positions, indices = get_c_alpha_positions(first(eachchain(pdb)))

    @test positions == [0 0 0
                        0.605f0 0.980f0 3.639f0]'
    @test length(indices)==2
end