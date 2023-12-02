@testitem "ColorMappings" begin
    using BiochemicalAlgorithms: SecondaryStructure
    using BiochemicalVisualization: SS_COLORS, AA_COLORS

    structure_keys = keys(SS_COLORS)
    @test SecondaryStructure.NONE ∈ structure_keys
    @test SecondaryStructure.HELIX ∈ structure_keys
    @test SecondaryStructure.SHEET ∈ structure_keys
    for structure in [SecondaryStructure.NONE, SecondaryStructure.HELIX, SecondaryStructure.SHEET]
        @test typeof(SS_COLORS[structure]) == NTuple{3, Int}
        for color_channel in SS_COLORS[structure]
            @test 0<=color_channel && color_channel<=255
        end
    end


    @test length(AA_COLORS)>=20
    aa_keys = keys(AA_COLORS)
    for k in ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        @test k ∈ aa_keys
    end
end