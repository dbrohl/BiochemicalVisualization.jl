function get_structure_color_mapping()
    return Dict(
        SecondaryStructure.NONE => (255, 255, 255), 
        SecondaryStructure.HELIX => (255, 75, 120), 
        SecondaryStructure.SHEET => (255, 150, 0))
end

function get_amino_acid_color_mapping()
    result = Dict{String, NTuple{3, Int}}()
    for ((aa, aa_three_letters, aa_one_letter), color) in zip(values(BiochemicalAlgorithms.AminoAcidProperties), distinguishable_colors(length(AminoAcidProperties)))
        result[aa_three_letters] = map(channel->Int(channel*255), (color.r, color.g, color.b))
    end
    return result
end