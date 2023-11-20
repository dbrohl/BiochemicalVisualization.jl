@testitem "ColorMappings" begin
    using BiochemicalAlgorithms: SecondaryStructure
    using BiochemicalVisualization: get_structure_color_mapping, get_amino_acid_color_mapping

    structure_keys = keys(get_structure_color_mapping())
    @test SecondaryStructure.NONE ∈ structure_keys
    @test SecondaryStructure.HELIX ∈ structure_keys
    @test SecondaryStructure.SHEET ∈ structure_keys
    for structure in [SecondaryStructure.NONE, SecondaryStructure.HELIX, SecondaryStructure.SHEET]
        @test typeof(get_structure_color_mapping()[structure]) == NTuple{3, Int}
        for color_channel in get_structure_color_mapping()[structure]
            @test 0<=color_channel && color_channel<=255
        end
    end


    @test length(get_amino_acid_color_mapping())>=20
    aa_keys = keys(get_amino_acid_color_mapping())
    for k in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        @test k ∈ aa_keys
    end
end